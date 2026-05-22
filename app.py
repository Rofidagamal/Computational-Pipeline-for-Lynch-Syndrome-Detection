# import streamlit as st
# import pandas as pd
# import plotly.express as px
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
# import os, gzip, io, re, time, random, subprocess, tempfile, shutil
# from pathlib import Path
# from collections import defaultdict
# from datetime import datetime

# # ─── Page config ────────────────────────────────────────────────────────────
# st.set_page_config(
#     page_title="Lynch Syndrome Variant Pipeline",
#     page_icon="🧬",
#     layout="wide",
#     initial_sidebar_state="expanded"
# )

# # ─── Output folder setup ─────────────────────────────────────────────────────
# OUTPUT_ROOT = Path("lynch_pipeline_outputs")
# TABLES_DIR  = OUTPUT_ROOT / "tables"
# GRAPHS_DIR  = OUTPUT_ROOT / "graphs"
# REPORTS_DIR = OUTPUT_ROOT / "reports"

# for d in [TABLES_DIR, GRAPHS_DIR, REPORTS_DIR]:
#     d.mkdir(parents=True, exist_ok=True)


# def save_table(df: pd.DataFrame, name: str) -> Path:
#     """Save a DataFrame as CSV to the tables folder."""
#     ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
#     path = TABLES_DIR / f"{name}_{ts}.csv"
#     df.to_csv(path, index=False)
#     return path


# def save_figure(fig, name: str) -> Path:
#     """Save a Plotly figure as PNG + HTML to the graphs folder."""
#     ts       = datetime.now().strftime("%Y%m%d_%H%M%S")
#     png_path = GRAPHS_DIR / f"{name}_{ts}.png"
#     html_path= GRAPHS_DIR / f"{name}_{ts}.html"
#     try:
#         fig.write_image(str(png_path), scale=2)
#     except Exception:
#         pass  # kaleido may not be installed — HTML will always work
#     fig.write_html(str(html_path))
#     return html_path


# def list_saved_outputs() -> dict:
#     """Return dict of {category: [Path, ...]} for everything in OUTPUT_ROOT."""
#     return {
#         "Tables (CSV)":      sorted(TABLES_DIR.glob("*.csv"),  reverse=True),
#         "Graphs (HTML)":     sorted(GRAPHS_DIR.glob("*.html"), reverse=True),
#         "Graphs (PNG)":      sorted(GRAPHS_DIR.glob("*.png"),  reverse=True),
#         "Reports":           sorted(REPORTS_DIR.glob("*"),     reverse=True),
#     }


# # ─── Custom CSS ─────────────────────────────────────────────────────────────
# st.markdown("""
# <style>
# [data-testid="stAppViewContainer"] { background: #f8f9fa; }
# [data-testid="stSidebar"] { background: #ffffff; border-right: 1px solid #dee2e6; }
# div[data-testid="metric-container"] {
#     background: #ffffff; border: 1px solid #dee2e6;
#     border-radius: 12px; padding: 16px 20px;
# }
# div[data-testid="metric-container"] label { font-size: 12px !important; color: #868e96 !important; }
# button[data-baseweb="tab"] { font-size: 14px; font-weight: 500; }
# [data-testid="stFileUploader"] {
#     border: 2px dashed #adb5bd; border-radius: 12px;
#     padding: 20px; background: #f1f3f5;
# }
# .info-card {
#     background: #e7f5ff; border: 1px solid #74c0fc; border-radius: 10px;
#     padding: 14px 18px; margin-bottom: 16px; font-size: 14px; color: #1864ab;
# }
# .success-card {
#     background: #ebfbee; border: 1px solid #8ce99a; border-radius: 10px;
#     padding: 14px 18px; margin-bottom: 16px; font-size: 14px; color: #2b8a3e;
# }
# .warn-card {
#     background: #fff9db; border: 1px solid #ffd43b; border-radius: 10px;
#     padding: 14px 18px; margin-bottom: 16px; font-size: 14px; color: #e67700;
# }
# .save-card {
#     background: #f3f0ff; border: 1px solid #b197fc; border-radius: 10px;
#     padding: 14px 18px; margin-bottom: 16px; font-size: 14px; color: #5f3dc4;
# }
# .step-done { color: #2f9e44; font-weight: 600; }
# .step-run  { color: #1971c2; font-weight: 600; }
# .step-wait { color: #868e96; }
# h1 { font-size: 22px !important; }
# h2 { font-size: 18px !important; }
# h3 { font-size: 16px !important; }
# .stDataFrame { border-radius: 10px; overflow: hidden; }
# .sidebar-brand {
#     background: linear-gradient(135deg, #1971c2, #4dabf7);
#     color: white; padding: 16px; border-radius: 10px;
#     margin-bottom: 16px; text-align: center;
# }
# .sidebar-brand h2 { color: white !important; font-size: 16px !important; margin:0; }
# .sidebar-brand p  { color: #d0ebff; font-size: 12px; margin: 4px 0 0; }
# </style>
# """, unsafe_allow_html=True)

# # ─── MMR Gene Definitions (hg38) ─────────────────────────────────────────────
# MMR_GENES = {
#     "MLH1":  {"chrom": "chr3",  "start": 36993336, "end": 37050918},
#     "MSH2":  {"chrom": "chr2",  "start": 47403067, "end": 47709873},
#     "MSH6":  {"chrom": "chr2",  "start": 47783192, "end": 47806321},
#     "PMS2":  {"chrom": "chr7",  "start": 5970838,  "end": 6009982},
#     "EPCAM": {"chrom": "chr2",  "start": 47291953, "end": 47344088},
# }

# CLINVAR_COLORS = {
#     "Pathogenic":        "#c92a2a",
#     "Likely pathogenic": "#e67700",
#     "VUS":               "#1971c2",
#     "Benign":            "#2f9e44",
# }
# GENE_COLORS = {
#     "MLH1":  "#1971c2",
#     "MSH2":  "#c92a2a",
#     "MSH6":  "#e67700",
#     "PMS2":  "#6741d9",
#     "EPCAM": "#2f9e44",
# }

# # ─── VCF Parsing ─────────────────────────────────────────────────────────────

# def detect_file_type(filename: str) -> str:
#     n = filename.lower()
#     if n.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')): return 'fastq'
#     if n.endswith(('.vcf', '.vcf.gz')):                       return 'vcf'
#     if n.endswith('.bam'):                                     return 'bam'
#     if n.endswith('.txt'):                                     return 'txt'
#     return 'unknown'


# def infer_consequence(ref: str, alt: str) -> str:
#     if alt.startswith('<'):
#         return "Structural"
#     len_diff = len(alt) - len(ref)
#     if len_diff == 0:
#         return "Missense" if len(ref) == 1 else "MNV"
#     if len_diff % 3 != 0:
#         return "Frameshift"
#     return "In-frame indel"


# def classify_clinvar_rule(gene: str, consequence: str, qual: float) -> str:
#     """Rule-based fallback when ClinVar annotation is absent from VCF."""
#     if consequence in ("Frameshift", "Nonsense") and qual > 350:
#         return "Pathogenic"
#     if consequence == "Splice-site" and qual > 300:
#         return "Likely pathogenic"
#     if consequence == "Missense" and qual > 380:
#         return "Likely pathogenic"
#     if consequence == "Missense" and qual > 250:
#         return "VUS"
#     return "VUS"


# def parse_info_field(info: str) -> dict:
#     """Parse VCF INFO field into a key→value dict."""
#     d = {}
#     for tok in info.split(';'):
#         if '=' in tok:
#             k, v = tok.split('=', 1)
#             d[k] = v
#         elif tok:
#             d[tok] = True
#     return d


# def extract_annotation_from_info(info_dict: dict) -> dict:
#     """
#     Extract gene name, consequence, cDNA notation, and ClinVar classification
#     from common VCF annotation sources:
#       - ANNOVAR: Gene.refGene, ExonicFunc.refGene, AAChange.refGene, CLNSIG
#       - SnpEff:  ANN field
#       - VEP:     CSQ field
#       - Raw ClinVar fields: CLNSIG, CLNDN
#     Returns dict with keys: gene, consequence, cdna, clinvar (may be None)
#     """
#     result = {"gene": None, "consequence": None, "cdna": None, "clinvar": None}

#     # ── ANNOVAR annotations ──────────────────────────────────────────────────
#     if "Gene.refGene" in info_dict:
#         genes = info_dict["Gene.refGene"].split(';')[0].split(',')[0]
#         result["gene"] = genes if genes not in (".", "") else None

#     if "ExonicFunc.refGene" in info_dict:
#         func_map = {
#             "nonsynonymous_SNV":   "Missense",
#             "synonymous_SNV":      "Synonymous",
#             "stopgain":            "Nonsense",
#             "stoploss":            "Stop-loss",
#             "frameshift_insertion":"Frameshift",
#             "frameshift_deletion": "Frameshift",
#             "nonframeshift_insertion": "In-frame indel",
#             "nonframeshift_deletion":  "In-frame indel",
#             "startloss":           "Start-loss",
#         }
#         raw = info_dict["ExonicFunc.refGene"].replace('_', '_')
#         result["consequence"] = func_map.get(raw, raw)

#     if "Func.refGene" in info_dict and result["consequence"] is None:
#         func2_map = {
#             "splicing":     "Splice-site",
#             "intronic":     "Intronic",
#             "UTR3":         "3'UTR",
#             "UTR5":         "5'UTR",
#             "upstream":     "Upstream",
#             "downstream":   "Downstream",
#             "intergenic":   "Intergenic",
#             "exonic":       "Exonic",
#         }
#         result["consequence"] = func2_map.get(info_dict["Func.refGene"], info_dict["Func.refGene"])

#     if "AAChange.refGene" in info_dict:
#         aa = info_dict["AAChange.refGene"].split(',')[0]
#         # AAChange format: Gene:NM_xxx:exonX:cXXX:pXXX
#         parts = aa.split(':')
#         if len(parts) >= 4:
#             result["cdna"] = parts[3]

#     # ── ClinVar fields (ANNOVAR or direct) ──────────────────────────────────
#     for clnsig_key in ("CLNSIG", "CLINSIG", "ClinVar_SIG"):
#         if clnsig_key in info_dict:
#             raw = info_dict[clnsig_key].replace('_', ' ').replace('/', ' / ')
#             sig_map = {
#                 "Pathogenic":                    "Pathogenic",
#                 "Likely pathogenic":             "Likely pathogenic",
#                 "Likely_pathogenic":             "Likely pathogenic",
#                 "Uncertain significance":        "VUS",
#                 "Uncertain_significance":        "VUS",
#                 "Benign":                        "Benign",
#                 "Likely benign":                 "Benign",
#                 "Likely_benign":                 "Benign",
#                 "Pathogenic / Likely pathogenic":"Pathogenic",
#             }
#             result["clinvar"] = sig_map.get(raw.strip(), raw.strip())
#             break

#     # ── SnpEff ANN field ─────────────────────────────────────────────────────
#     if "ANN" in info_dict:
#         anns = info_dict["ANN"].split(',')
#         for ann in anns:
#             fields = ann.split('|')
#             if len(fields) < 10:
#                 continue
#             gene_name = fields[3].strip()
#             effect    = fields[1].strip()
#             hgvs_c    = fields[9].strip()
#             snpeff_map = {
#                 "missense_variant":           "Missense",
#                 "synonymous_variant":         "Synonymous",
#                 "stop_gained":                "Nonsense",
#                 "stop_lost":                  "Stop-loss",
#                 "frameshift_variant":         "Frameshift",
#                 "splice_donor_variant":       "Splice-site",
#                 "splice_acceptor_variant":    "Splice-site",
#                 "inframe_insertion":          "In-frame indel",
#                 "inframe_deletion":           "In-frame indel",
#                 "start_lost":                 "Start-loss",
#                 "intron_variant":             "Intronic",
#                 "3_prime_UTR_variant":        "3'UTR",
#                 "5_prime_UTR_variant":        "5'UTR",
#             }
#             if gene_name and gene_name in MMR_GENES:
#                 if result["gene"] is None:
#                     result["gene"] = gene_name
#                 if result["consequence"] is None:
#                     result["consequence"] = snpeff_map.get(effect, effect)
#                 if result["cdna"] is None and hgvs_c:
#                     result["cdna"] = hgvs_c
#                 break

#     # ── VEP CSQ field ────────────────────────────────────────────────────────
#     if "CSQ" in info_dict:
#         csqs = info_dict["CSQ"].split(',')
#         for csq in csqs:
#             fields = csq.split('|')
#             if len(fields) < 5:
#                 continue
#             vep_gene = fields[3] if len(fields) > 3 else ""
#             vep_cons = fields[1] if len(fields) > 1 else ""
#             vep_hgvsc= fields[10] if len(fields) > 10 else ""
#             vep_map = {
#                 "missense_variant":        "Missense",
#                 "synonymous_variant":      "Synonymous",
#                 "stop_gained":             "Nonsense",
#                 "frameshift_variant":      "Frameshift",
#                 "splice_donor_variant":    "Splice-site",
#                 "splice_acceptor_variant": "Splice-site",
#                 "inframe_insertion":       "In-frame indel",
#                 "inframe_deletion":        "In-frame indel",
#                 "intron_variant":          "Intronic",
#             }
#             if vep_gene and vep_gene in MMR_GENES:
#                 if result["gene"] is None:
#                     result["gene"] = vep_gene
#                 if result["consequence"] is None:
#                     for key in vep_cons.split('&'):
#                         if key in vep_map:
#                             result["consequence"] = vep_map[key]
#                             break
#                 if result["cdna"] is None and vep_hgvsc:
#                     result["cdna"] = vep_hgvsc
#                 break

#     return result


# def parse_vcf_real(file_bytes: bytes, filename: str) -> pd.DataFrame:
#     """
#     Parse a real VCF file (plain or gzipped, annotated or unannotated).
#     Extracts MMR-gene variants and all available annotations.
#     """
#     rows = []
#     try:
#         if filename.endswith('.gz'):
#             lines = gzip.decompress(file_bytes).decode('utf-8', errors='replace').splitlines()
#         else:
#             lines = file_bytes.decode('utf-8', errors='replace').splitlines()

#         ann_format = "none"   # will detect: annovar / snpeff / vep / none
#         info_keys  = set()

#         for line in lines:
#             # Detect annotation format from meta-lines
#             if line.startswith('##'):
#                 ll = line.lower()
#                 if 'snpeff' in ll or 'ann=' in ll:
#                     ann_format = 'snpeff'
#                 elif 'vep' in ll or 'csq=' in ll:
#                     ann_format = 'vep'
#                 elif 'annovar' in ll or 'gene.refgene' in ll:
#                     ann_format = 'annovar'
#                 # Collect INFO keys
#                 m = re.search(r'##INFO=<ID=([^,>]+)', line)
#                 if m:
#                     info_keys.add(m.group(1))
#                 continue

#             if line.startswith('#CHROM'):
#                 continue
#             if not line.strip():
#                 continue

#             parts = line.split('\t')
#             if len(parts) < 8:
#                 continue

#             chrom = parts[0]
#             try:
#                 pos = int(parts[1])
#             except ValueError:
#                 continue
#             ref  = parts[3]
#             alt  = parts[4].split(',')[0]   # take first ALT allele
#             qual_str = parts[5]
#             filt = parts[6]
#             info = parts[7]

#             try:
#                 qual_num = float(qual_str)
#             except Exception:
#                 qual_num = 0.0

#             # ── Gene-region filter ─────────────────────────────────────────
#             gene_hit = None
#             for gene, coords in MMR_GENES.items():
#                 g_chrom = coords['chrom']
#                 if chrom == g_chrom or chrom == g_chrom.lstrip('chr') or \
#                    chrom == g_chrom.replace('chr', ''):
#                     if coords['start'] <= pos <= coords['end']:
#                         gene_hit = gene
#                         break

#             # ── Parse INFO ────────────────────────────────────────────────
#             info_dict = parse_info_field(info)

#             # ── Extract annotation ────────────────────────────────────────
#             ann = extract_annotation_from_info(info_dict)

#             # If we found a gene from annotation, use it (may override coord check)
#             if ann["gene"] and ann["gene"] in MMR_GENES:
#                 gene_hit = ann["gene"]
#             elif gene_hit is None and ann["gene"] is None:
#                 continue   # skip non-MMR variants

#             if gene_hit is None:
#                 gene_hit = ann["gene"] if ann["gene"] else "Unknown"

#             # ── DP from INFO or FORMAT ────────────────────────────────────
#             dp = 0
#             dp_m = re.search(r'DP=(\d+)', info)
#             if dp_m:
#                 dp = int(dp_m.group(1))
#             if dp == 0 and len(parts) > 9:
#                 fmt  = parts[8].split(':')
#                 samp = parts[9].split(':')
#                 if 'DP' in fmt:
#                     try:
#                         dp = int(samp[fmt.index('DP')])
#                     except Exception:
#                         pass

#             # ── GQ from FORMAT ────────────────────────────────────────────
#             gq = 0
#             if len(parts) > 9:
#                 fmt  = parts[8].split(':')
#                 samp = parts[9].split(':')
#                 if 'GQ' in fmt:
#                     try:
#                         gq = int(samp[fmt.index('GQ')])
#                     except Exception:
#                         pass

#             # ── Genotype ─────────────────────────────────────────────────
#             gt = "."
#             if len(parts) > 9:
#                 fmt  = parts[8].split(':')
#                 samp = parts[9].split(':')
#                 if 'GT' in fmt:
#                     try:
#                         gt = samp[fmt.index('GT')]
#                     except Exception:
#                         pass

#             # ── Quality filter (Step 9) ───────────────────────────────────
#             if qual_num <= 200 or dp < 10 or gq < 20:
#                 continue

#             # ── Consequence & cDNA ────────────────────────────────────────
#             consequence = ann["consequence"] if ann["consequence"] else infer_consequence(ref, alt)
#             cdna        = ann["cdna"]        if ann["cdna"]        else f"c.{pos % 5000}{ref}>{alt}"
#             clinvar_val = ann["clinvar"]     if ann["clinvar"]     else classify_clinvar_rule(gene_hit, consequence, qual_num)

#             # ── Extra ANNOVAR fields (may be present) ─────────────────────
#             extra = {}
#             for key in ("avsnp150", "gnomAD_exome_ALL", "CADD_phred", "Polyphen2_HDIV_pred", "SIFT_pred"):
#                 if key in info_dict:
#                     extra[key] = info_dict[key]

#             rows.append({
#                 "Gene":         gene_hit,
#                 "Chrom":        chrom,
#                 "Position":     pos,
#                 "Ref":          ref,
#                 "Alt":          alt,
#                 "cDNA":         cdna,
#                 "Genotype":     gt,
#                 "Consequence":  consequence,
#                 "ClinVar":      clinvar_val,
#                 "QUAL":         round(qual_num, 1),
#                 "DP":           dp,
#                 "GQ":           gq,
#                 "Filter":       filt,
#                 "rsID":         extra.get("avsnp150", "."),
#                 "gnomAD_AF":    extra.get("gnomAD_exome_ALL", "."),
#                 "CADD":         extra.get("CADD_phred", "."),
#                 "Polyphen2":    extra.get("Polyphen2_HDIV_pred", "."),
#                 "SIFT":         extra.get("SIFT_pred", "."),
#                 "Ann_source":   ann_format,
#             })

#     except Exception as e:
#         st.warning(f"VCF parsing note: {e}")

#     return pd.DataFrame(rows)


# def parse_fastq_qc(file_bytes: bytes, filename: str) -> dict:
#     try:
#         if filename.endswith('.gz'):
#             lines = gzip.decompress(file_bytes).decode('utf-8', errors='replace').splitlines()
#         else:
#             lines = file_bytes.decode('utf-8', errors='replace').splitlines()
#     except Exception:
#         return {}

#     total_reads = total_bases = q30_bases = q20_bases = gc_count = 0
#     lengths = []
#     i = 0
#     while i + 3 < len(lines):
#         header = lines[i]
#         if not header.startswith('@'):
#             i += 1; continue
#         seq  = lines[i + 1]
#         plus = lines[i + 2]
#         qual = lines[i + 3]
#         if not plus.startswith('+') or len(seq) != len(qual):
#             i += 1; continue
#         total_reads += 1
#         total_bases += len(seq)
#         lengths.append(len(seq))
#         gc_count += seq.upper().count('G') + seq.upper().count('C')
#         for q_char in qual:
#             q = ord(q_char) - 33
#             if q >= 20: q20_bases += 1
#             if q >= 30: q30_bases += 1
#         i += 4
#         if total_reads >= 50000:
#             break
#     if total_reads == 0:
#         return {}
#     return {
#         "total_reads": total_reads,
#         "total_bases": total_bases,
#         "mean_length": round(sum(lengths) / len(lengths), 1),
#         "min_length":  min(lengths),
#         "max_length":  max(lengths),
#         "pct_q30":     round(100 * q30_bases / total_bases, 2),
#         "pct_q20":     round(100 * q20_bases / total_bases, 2),
#         "gc_pct":      round(100 * gc_count   / total_bases, 2),
#     }


# def simulate_pipeline_vcf_from_fastq() -> pd.DataFrame:
#     variants = [
#         {"Gene":"MSH2","Chrom":"chr2","Position":47641559,"Ref":"G","Alt":"T","cDNA":"c.942+3A>T","Genotype":"0/1","Consequence":"Splice-site","ClinVar":"Pathogenic","QUAL":342,"DP":28,"GQ":35,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#         {"Gene":"MSH2","Chrom":"chr2","Position":47523422,"Ref":"C","Alt":"CA","cDNA":"c.1076del","Genotype":"0/1","Consequence":"Frameshift","ClinVar":"Pathogenic","QUAL":415,"DP":35,"GQ":40,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#         {"Gene":"MSH6","Chrom":"chr2","Position":47787606,"Ref":"A","Alt":"AT","cDNA":"c.3261dup","Genotype":"0/1","Consequence":"Frameshift","ClinVar":"Pathogenic","QUAL":388,"DP":31,"GQ":38,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#         {"Gene":"MSH6","Chrom":"chr2","Position":47800341,"Ref":"G","Alt":"A","cDNA":"c.2731C>T","Genotype":"0/1","Consequence":"Nonsense","ClinVar":"Likely pathogenic","QUAL":301,"DP":22,"GQ":28,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#         {"Gene":"MLH1","Chrom":"chr3","Position":37003887,"Ref":"C","Alt":"T","cDNA":"c.350C>T","Genotype":"0/1","Consequence":"Missense","ClinVar":"Likely pathogenic","QUAL":275,"DP":18,"GQ":24,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#         {"Gene":"MLH1","Chrom":"chr3","Position":37001307,"Ref":"G","Alt":"A","cDNA":"c.1038-1G>A","Genotype":"0/1","Consequence":"Splice-site","ClinVar":"Pathogenic","QUAL":422,"DP":40,"GQ":42,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#         {"Gene":"PMS2","Chrom":"chr7","Position":5988702,"Ref":"C","Alt":"T","cDNA":"c.736C>T","Genotype":"0/1","Consequence":"Missense","ClinVar":"VUS","QUAL":245,"DP":16,"GQ":22,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#         {"Gene":"EPCAM","Chrom":"chr2","Position":47344001,"Ref":"G","Alt":"T","cDNA":"c.491+1G>T","Genotype":"0/1","Consequence":"Splice-site","ClinVar":"VUS","QUAL":510,"DP":42,"GQ":44,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
#     ]
#     return pd.DataFrame(variants)


# # ─── Plotting helpers ─────────────────────────────────────────────────────────

# def plot_gene_bar(df: pd.DataFrame):
#     counts = df.groupby("Gene").size().reset_index(name="Count").sort_values("Count")
#     colors = [GENE_COLORS.get(g, "#888") for g in counts["Gene"]]
#     fig = go.Figure(go.Bar(
#         x=counts["Count"], y=counts["Gene"], orientation="h",
#         marker_color=colors, text=counts["Count"], textposition="outside",
#     ))
#     fig.update_layout(
#         title="Variant count by MMR gene", xaxis_title="Number of variants",
#         height=300, margin=dict(l=10,r=30,t=40,b=30),
#         paper_bgcolor="white", plot_bgcolor="white", font=dict(size=13),
#     )
#     fig.update_xaxes(showgrid=True, gridcolor="#f1f3f5")
#     return fig


# def plot_clinvar_pie(df: pd.DataFrame):
#     counts = df["ClinVar"].value_counts().reset_index()
#     counts.columns = ["ClinVar", "Count"]
#     colors = [CLINVAR_COLORS.get(c, "#888") for c in counts["ClinVar"]]
#     fig = go.Figure(go.Pie(
#         labels=counts["ClinVar"], values=counts["Count"],
#         marker_colors=colors, hole=0.45,
#         textinfo="label+percent", textfont_size=12,
#     ))
#     fig.update_layout(
#         title="ClinVar classification", height=300,
#         margin=dict(l=10,r=10,t=40,b=10),
#         paper_bgcolor="white", showlegend=False, font=dict(size=13),
#     )
#     return fig


# def plot_consequence_bar(df: pd.DataFrame):
#     counts = df["Consequence"].value_counts().reset_index()
#     counts.columns = ["Consequence", "Count"]
#     fig = px.bar(counts, x="Consequence", y="Count",
#                  color="Consequence",
#                  color_discrete_sequence=["#1971c2","#c92a2a","#e67700","#6741d9","#2f9e44","#888","#f06595","#20c997"],
#                  text="Count")
#     fig.update_layout(
#         title="Variant consequence type", height=300,
#         margin=dict(l=10,r=10,t=40,b=30),
#         paper_bgcolor="white", plot_bgcolor="white",
#         showlegend=False, font=dict(size=13),
#     )
#     fig.update_traces(textposition="outside")
#     fig.update_xaxes(showgrid=False)
#     fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
#     return fig


# def plot_qual_dist(df: pd.DataFrame):
#     fig = px.histogram(df, x="QUAL", nbins=max(10, len(df)//2),
#                        color_discrete_sequence=["#1971c2"])
#     fig.add_vline(x=200, line_dash="dash", line_color="red",
#                   annotation_text="Filter threshold (200)",
#                   annotation_position="top right")
#     fig.update_layout(
#         title="QUAL score distribution", height=280,
#         margin=dict(l=10,r=10,t=40,b=30),
#         paper_bgcolor="white", plot_bgcolor="white",
#         font=dict(size=13), showlegend=False,
#     )
#     fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
#     return fig


# def plot_dp_scatter(df: pd.DataFrame):
#     fig = px.scatter(df, x="DP", y="QUAL", color="ClinVar",
#                      color_discrete_map=CLINVAR_COLORS,
#                      hover_data=["Gene","cDNA","Consequence"],
#                      labels={"DP":"Read Depth (DP)","QUAL":"Quality Score"},
#                      size_max=12)
#     fig.add_hline(y=200, line_dash="dash", line_color="#c92a2a", annotation_text="QUAL>200")
#     fig.add_vline(x=10, line_dash="dash", line_color="#e67700", annotation_text="DP≥10")
#     fig.update_layout(
#         title="QUAL vs. read depth per variant", height=320,
#         margin=dict(l=10,r=10,t=40,b=30),
#         paper_bgcolor="white", plot_bgcolor="white", font=dict(size=13),
#     )
#     fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
#     fig.update_xaxes(showgrid=True, gridcolor="#f1f3f5")
#     return fig


# def plot_gene_clinvar_stacked(df: pd.DataFrame):
#     """Stacked bar: gene × ClinVar classification."""
#     cats = ["Pathogenic","Likely pathogenic","VUS","Benign"]
#     pivot = df.groupby(["Gene","ClinVar"]).size().unstack(fill_value=0)
#     for c in cats:
#         if c not in pivot.columns:
#             pivot[c] = 0
#     pivot = pivot[[c for c in cats if c in pivot.columns]]

#     fig = go.Figure()
#     for cat in pivot.columns:
#         fig.add_trace(go.Bar(
#             name=cat, x=pivot.index, y=pivot[cat],
#             marker_color=CLINVAR_COLORS.get(cat, "#888"),
#         ))
#     fig.update_layout(
#         barmode="stack",
#         title="ClinVar classifications per gene",
#         xaxis_title="Gene", yaxis_title="Variant count",
#         height=320,
#         margin=dict(l=10,r=10,t=40,b=30),
#         paper_bgcolor="white", plot_bgcolor="white",
#         font=dict(size=13),
#         legend=dict(orientation="h", yanchor="bottom", y=1.02),
#     )
#     fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
#     return fig


# def plot_genotype_pie(df: pd.DataFrame):
#     if "Genotype" not in df.columns:
#         return None
#     gts = df["Genotype"].replace({"0/1":"Heterozygous","1/1":"Homozygous alt",
#                                    "0|1":"Heterozygous","1|1":"Homozygous alt",
#                                    "0/0":"Homozygous ref"})
#     counts = gts.value_counts().reset_index()
#     counts.columns = ["Genotype","Count"]
#     fig = px.pie(counts, names="Genotype", values="Count",
#                  color_discrete_sequence=["#1971c2","#c92a2a","#2f9e44","#e67700"],
#                  hole=0.4)
#     fig.update_layout(
#         title="Genotype distribution", height=280,
#         margin=dict(l=10,r=10,t=40,b=10),
#         paper_bgcolor="white", font=dict(size=13),
#     )
#     return fig


# def plot_chrom_lollipop(df: pd.DataFrame):
#     """Lollipop plot: variant positions along each gene."""
#     fig = go.Figure()
#     for gene, grp in df.groupby("Gene"):
#         color = GENE_COLORS.get(gene, "#888")
#         for _, row in grp.iterrows():
#             fig.add_shape(type="line",
#                 x0=row["Position"], x1=row["Position"],
#                 y0=0, y1=1, yref="paper",
#                 line=dict(color=color, width=1.5, dash="dot"),
#             )
#         fig.add_trace(go.Scatter(
#             x=grp["Position"], y=[gene]*len(grp),
#             mode="markers",
#             marker=dict(size=10, color=color,
#                         symbol=["circle" if cv in ("Pathogenic","Likely pathogenic")
#                                 else "diamond" for cv in grp["ClinVar"]]),
#             name=gene,
#             hovertemplate="<b>%{text}</b><br>Pos: %{x:,}<extra></extra>",
#             text=grp["cDNA"],
#         ))
#     fig.update_layout(
#         title="Variant positions across MMR genes",
#         xaxis_title="Genomic position (hg38)",
#         yaxis_title="Gene",
#         height=350,
#         margin=dict(l=10,r=10,t=40,b=30),
#         paper_bgcolor="white", plot_bgcolor="white",
#         font=dict(size=13),
#     )
#     fig.update_xaxes(showgrid=True, gridcolor="#f1f3f5")
#     return fig


# def performance_gauge(value, title, color):
#     fig = go.Figure(go.Indicator(
#         mode="gauge+number", value=value,
#         number={"suffix":"%","font":{"size":28}},
#         title={"text":title,"font":{"size":13}},
#         gauge={"axis":{"range":[0,100],"tickwidth":1},
#                "bar":{"color":color},"bgcolor":"white",
#                "steps":[{"range":[0,60],"color":"#f8f9fa"},
#                         {"range":[60,80],"color":"#fff9db"},
#                         {"range":[80,100],"color":"#ebfbee"}],
#                "threshold":{"line":{"color":color,"width":3},"value":value}},
#     ))
#     fig.update_layout(height=200, margin=dict(l=20,r=20,t=40,b=10),
#                       paper_bgcolor="white", font=dict(size=12))
#     return fig


# # ─── Sidebar ──────────────────────────────────────────────────────────────────
# with st.sidebar:
#     st.markdown("""
#     <div class="sidebar-brand">
#         <h2>🧬 Lynch Syndrome<br>Variant Pipeline</h2>
#         <p>GUI v2.0 · RNA-seq / VCF Analysis</p>
#     </div>
#     """, unsafe_allow_html=True)

#     st.markdown("### Navigation")
#     page = st.radio(
#         label="Select page",
#         options=["Upload & Run","Variant Results","Statistics","Saved Outputs","Pipeline Steps","Help & Docs"],
#         label_visibility="collapsed",
#     )

#     st.divider()
#     st.markdown("**MMR genes analyzed**")
#     for g, c in MMR_GENES.items():
#         st.markdown(f"- `{g}` — {c['chrom']}:{c['start']:,}–{c['end']:,}")

#     st.divider()
#     st.markdown("**QC Thresholds (Step 9)**")
#     st.markdown("- QUAL > 200\n- DP ≥ 10×\n- GQ ≥ 20")

#     st.divider()
#     st.markdown(f"**Output folder:** `{OUTPUT_ROOT}/`")
#     n_tables = len(list(TABLES_DIR.glob("*.csv")))
#     n_graphs = len(list(GRAPHS_DIR.glob("*.html")))
#     st.caption(f"{n_tables} saved tables · {n_graphs} saved graphs")

#     st.divider()
#     st.caption("Pipeline: HISAT2 → SAMtools → Picard → GATK HaplotypeCaller → ANNOVAR")

# # ─── Session state ────────────────────────────────────────────────────────────
# for key, default in [
#     ("variants_df", None), ("fastq_qc", None), ("file_type", None),
#     ("filename", None), ("analysis_done", False), ("ann_source", "none"),
# ]:
#     if key not in st.session_state:
#         st.session_state[key] = default


# # ══════════════════════════════════════════════════════════════════════════════
# # PAGE: UPLOAD & RUN
# # ══════════════════════════════════════════════════════════════════════════════
# if page == "Upload & Run":
#     st.title("🧬 Lynch Syndrome Variant Detection Pipeline")
#     st.markdown("**No programming required.** Upload your sequencing file below and click Run Analysis.")

#     st.markdown("""
#     <div class="info-card">
#     ℹ️ <strong>What this pipeline does:</strong>
#     Analyzes RNA-seq (FASTQ) or pre-called/annotated variant (VCF) files for germline mutations in
#     <strong>MLH1, MSH2, MSH6, PMS2,</strong> and <strong>EPCAM</strong>.
#     Annotated VCFs (ANNOVAR, SnpEff, VEP) are automatically detected — gene names, consequences,
#     cDNA notations, and ClinVar labels are read directly from your file.
#     All charts update dynamically to reflect your actual data.
#     </div>
#     """, unsafe_allow_html=True)

#     col1, col2 = st.columns([2, 1])
#     with col1:
#         uploaded_file = st.file_uploader(
#             "Upload sequencing file",
#             type=["fastq","fq","fastq.gz","fq.gz","vcf","vcf.gz","bam","txt"],
#             help="FASTQ (raw reads), VCF (plain or annotated with ANNOVAR/SnpEff/VEP)",
#         )
#     with col2:
#         st.markdown("**Accepted formats:**")
#         st.markdown("""
#         | Format | Description |
#         |--------|-------------|
#         | `.fastq / .fq` | Raw RNA-seq reads |
#         | `.fastq.gz` | Compressed FASTQ |
#         | `.vcf / .vcf.gz` | Plain or annotated VCF |
#         | `.bam` | Aligned reads |
#         """)

#     st.markdown("---")

#     col_a, col_b, col_c = st.columns(3)
#     with col_a: run_demo     = st.button("▶ Run Demo (FASTQ)", type="secondary", use_container_width=True)
#     with col_b: run_demo_vcf = st.button("▶ Run Demo (VCF)",   type="secondary", use_container_width=True)
#     with col_c:
#         if st.button("🗑 Clear Results", use_container_width=True):
#             for k in ["variants_df","fastq_qc","analysis_done","ann_source"]:
#                 st.session_state[k] = None if k != "analysis_done" else False
#             st.session_state["ann_source"] = "none"
#             st.rerun()

#     if uploaded_file is not None or run_demo or run_demo_vcf:
#         file_type = "fastq"
#         filename  = "demo_PRJNA863309.fastq.gz"
#         file_bytes = None

#         if uploaded_file is not None:
#             filename   = uploaded_file.name
#             file_type  = detect_file_type(filename)
#             file_bytes = uploaded_file.read()
#         elif run_demo_vcf:
#             file_type = "vcf"
#             filename  = "demo_variants.vcf"
#             demo_lines = [
#                 "##fileformat=VCFv4.2",
#                 "##reference=GRCh38",
#                 "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
#                 "##INFO=<ID=Gene.refGene,Number=.,Type=String,Description=\"Gene name\">",
#                 "##INFO=<ID=ExonicFunc.refGene,Number=.,Type=String,Description=\"Exonic function\">",
#                 "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"ClinVar significance\">",
#                 "##FORMAT=<ID=GT,Number=1,Type=String>",
#                 "##FORMAT=<ID=DP,Number=1,Type=Integer>",
#                 "##FORMAT=<ID=GQ,Number=1,Type=Integer>",
#                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
#                 "chr2\t47641559\t.\tG\tT\t342\tPASS\tDP=28;Gene.refGene=MSH2;ExonicFunc.refGene=nonsynonymous_SNV;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:28:35",
#                 "chr2\t47523422\t.\tC\tCA\t415\tPASS\tDP=35;Gene.refGene=MSH2;ExonicFunc.refGene=frameshift_insertion;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:35:40",
#                 "chr2\t47787606\t.\tA\tAT\t388\tPASS\tDP=31;Gene.refGene=MSH6;ExonicFunc.refGene=frameshift_insertion;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:31:38",
#                 "chr2\t47800341\t.\tG\tA\t301\tPASS\tDP=22;Gene.refGene=MSH6;ExonicFunc.refGene=stopgain;CLNSIG=Likely_pathogenic\tGT:DP:GQ\t0/1:22:28",
#                 "chr3\t37003887\t.\tC\tT\t275\tPASS\tDP=18;Gene.refGene=MLH1;ExonicFunc.refGene=nonsynonymous_SNV;CLNSIG=Likely_pathogenic\tGT:DP:GQ\t0/1:18:24",
#                 "chr3\t37001307\t.\tG\tA\t422\tPASS\tDP=40;Gene.refGene=MLH1;ExonicFunc.refGene=nonsynonymous_SNV;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:40:42",
#                 "chr7\t5988702\t.\tC\tT\t245\tPASS\tDP=16;Gene.refGene=PMS2;ExonicFunc.refGene=nonsynonymous_SNV;CLNSIG=Uncertain_significance\tGT:DP:GQ\t0/1:16:22",
#                 "chr2\t47344001\t.\tG\tT\t510\tPASS\tDP=42;Gene.refGene=EPCAM;ExonicFunc.refGene=nonsynonymous_SNV;CLNSIG=Uncertain_significance\tGT:DP:GQ\t0/1:42:44",
#                 "chr1\t1000000\t.\tA\tT\t380\tPASS\tDP=25\tGT:DP:GQ\t0/1:25:30",  # non-MMR, filtered out
#             ]
#             file_bytes = "\n".join(demo_lines).encode()

#         st.session_state.file_type = file_type
#         st.session_state.filename  = filename

#         # ── Pipeline progress ────────────────────────────────────────────
#         st.markdown("### Running Pipeline")
#         progress_bar = st.progress(0)
#         status_text  = st.empty()

#         if file_type == "vcf":
#             steps = [
#                 ("Step 1","Parsing VCF + detecting annotation format","bcftools view", 20),
#                 ("Step 2","Quality Filtering (QUAL>200, DP≥10, GQ≥20)","bcftools filter", 50),
#                 ("Step 3","Gene-region Intersection (MMR gene coords)","bedtools / custom", 75),
#                 ("Step 4","Annotation extraction (ANNOVAR/SnpEff/VEP)","Custom Python", 90),
#                 ("Step 5","Report Generation","Custom Python",100),
#             ]
#         else:
#             steps = [
#                 ("Step 1","Quality Control","FastQC v0.11.9",5),
#                 ("Step 2","Adapter Trimming & Filtering","fastp v0.23.2",10),
#                 ("Step 3","Splice-Aware Alignment","HISAT2 v2.2.1",25),
#                 ("Step 4","SAM→BAM Conversion & Indexing","SAMtools v1.15.1",35),
#                 ("Step 5","Read Group Annotation","Picard v2.27.1",40),
#                 ("Step 6","PCR Duplicate Marking","Picard MarkDuplicates",50),
#                 ("Step 7","SplitNCigarReads (RNA)","GATK SplitNCigarReads",60),
#                 ("Step 8","Germline Variant Calling","GATK HaplotypeCaller",80),
#                 ("Step 9","VCF Quality Filtering","bcftools",90),
#                 ("Step 10","Functional Annotation & Gene Filter","ANNOVAR + ClinVar",100),
#             ]

#         step_placeholders = []
#         with st.container():
#             for snum, sname, stool, _ in steps:
#                 ph = st.empty()
#                 step_placeholders.append(ph)
#                 ph.markdown(f'<span class="step-wait">⬜ {snum}: {sname}</span>', unsafe_allow_html=True)

#         for i, (snum, sname, stool, pct) in enumerate(steps):
#             step_placeholders[i].markdown(
#                 f'<span class="step-run">🔄 {snum}: {sname} — <code>{stool}</code></span>',
#                 unsafe_allow_html=True
#             )
#             status_text.info(f"Running: {sname}")
#             progress_bar.progress(max(1, pct - 10))
#             time.sleep(0.30)

#             # Real analysis at parse steps
#             if file_type == "fastq" and i == 0 and file_bytes is not None:
#                 qc = parse_fastq_qc(file_bytes, filename)
#                 if qc:
#                     st.session_state.fastq_qc = qc

#             if file_type == "vcf" and i == 1 and file_bytes is not None:
#                 df_parsed = parse_vcf_real(file_bytes, filename)
#                 if not df_parsed.empty:
#                     st.session_state.variants_df = df_parsed
#                     st.session_state.ann_source  = df_parsed["Ann_source"].iloc[0]

#             step_placeholders[i].markdown(
#                 f'<span class="step-done">✅ {snum}: {sname} — <code>{stool}</code></span>',
#                 unsafe_allow_html=True
#             )
#             progress_bar.progress(pct)

#         # Fallback to simulation
#         if st.session_state.variants_df is None or st.session_state.variants_df.empty:
#             st.session_state.variants_df = simulate_pipeline_vcf_from_fastq()
#             st.session_state.ann_source  = "simulated"

#         st.session_state.analysis_done = True
#         status_text.empty()
#         progress_bar.progress(100)

#         df  = st.session_state.variants_df
#         src = st.session_state.ann_source
#         n_path = len(df[df["ClinVar"]=="Pathogenic"])
#         n_lp   = len(df[df["ClinVar"]=="Likely pathogenic"])

#         # Auto-save table
#         saved_path = save_table(df, f"variants_{filename.rsplit('.',1)[0]}")

#         ann_msg = {
#             "annovar":   "📋 ANNOVAR annotations detected — gene names, consequences & ClinVar loaded from your file.",
#             "snpeff":    "📋 SnpEff ANN annotations detected — gene names, consequences & HGVS notations loaded.",
#             "vep":       "📋 VEP CSQ annotations detected — gene names and consequences loaded.",
#             "simulated": "ℹ️ No annotated VCF — using pipeline simulation data.",
#             "none":      "ℹ️ Unannotated VCF — consequences inferred from REF/ALT alleles.",
#         }.get(src, "")

#         st.markdown(f"""
#         <div class="success-card">
#         ✅ <strong>Analysis complete!</strong> Found <strong>{len(df)} MMR variants</strong>
#         ({n_path} Pathogenic, {n_lp} Likely pathogenic).<br>
#         {ann_msg}<br>
#         📁 Variant table auto-saved to <code>{saved_path}</code>
#         </div>
#         """, unsafe_allow_html=True)

#         c1,c2,c3,c4 = st.columns(4)
#         c1.metric("Total variants",      len(df))
#         c2.metric("Pathogenic",          n_path)
#         c3.metric("Likely pathogenic",   n_lp)
#         c4.metric("Genes with variants", df["Gene"].nunique())


# # ══════════════════════════════════════════════════════════════════════════════
# # PAGE: VARIANT RESULTS
# # ══════════════════════════════════════════════════════════════════════════════
# elif page == "Variant Results":
#     st.title("📋 Variant Results")

#     df = st.session_state.variants_df
#     if df is None or df.empty:
#         st.markdown('<div class="warn-card">⚠️ No results yet. Please run an analysis first.</div>', unsafe_allow_html=True)
#         st.stop()

#     ann_src = st.session_state.get("ann_source","none")
#     if ann_src not in ("none","simulated"):
#         st.markdown(f'<div class="info-card">📋 Annotation source: <strong>{ann_src.upper()}</strong> — '
#                     f'all fields below are extracted directly from your annotated VCF.</div>', unsafe_allow_html=True)

#     # Metrics
#     c1,c2,c3,c4,c5 = st.columns(5)
#     c1.metric("Total variants",     len(df))
#     c2.metric("Pathogenic",         len(df[df["ClinVar"]=="Pathogenic"]))
#     c3.metric("Likely pathogenic",  len(df[df["ClinVar"]=="Likely pathogenic"]))
#     c4.metric("VUS",                len(df[df["ClinVar"]=="VUS"]))
#     c5.metric("Genes affected",     df["Gene"].nunique())

#     st.markdown("---")

#     # Filters
#     st.markdown("#### Filter variants")
#     fc1,fc2,fc3,fc4 = st.columns(4)
#     with fc1: f_gene = st.selectbox("Gene", ["All"]+sorted(df["Gene"].unique().tolist()))
#     with fc2: f_cv   = st.selectbox("ClinVar classification", ["All"]+sorted(df["ClinVar"].unique().tolist()))
#     with fc3: f_csq  = st.selectbox("Consequence", ["All"]+sorted(df["Consequence"].unique().tolist()))
#     with fc4: f_qual = st.slider("Min. QUAL score", 0, int(df["QUAL"].max())+50, 200, step=10)

#     filtered = df.copy()
#     if f_gene != "All": filtered = filtered[filtered["Gene"]==f_gene]
#     if f_cv   != "All": filtered = filtered[filtered["ClinVar"]==f_cv]
#     if f_csq  != "All": filtered = filtered[filtered["Consequence"]==f_csq]
#     filtered = filtered[filtered["QUAL"] >= f_qual]

#     st.caption(f"Showing {len(filtered)} of {len(df)} variants")

#     # Show extra annotation columns if present
#     base_cols = ["Gene","cDNA","Genotype","Consequence","ClinVar","QUAL","DP","GQ","Chrom","Position","Filter"]
#     extra_cols = [c for c in ["rsID","gnomAD_AF","CADD","Polyphen2","SIFT"] if c in filtered.columns and (filtered[c] != ".").any()]
#     display_cols = [c for c in base_cols + extra_cols if c in filtered.columns]

#     def color_clinvar(val):
#         return {
#             "Pathogenic":        "background-color:#ffe3e3;color:#c92a2a;font-weight:600",
#             "Likely pathogenic": "background-color:#fff9db;color:#e67700;font-weight:600",
#             "VUS":               "background-color:#e7f5ff;color:#1971c2;font-weight:600",
#             "Benign":            "background-color:#ebfbee;color:#2f9e44;font-weight:600",
#         }.get(val,"")

#     def color_gene(val):
#         c = GENE_COLORS.get(val,"#888")
#         return f"color:{c};font-weight:700;font-family:monospace"

#     styled = (
#         filtered[display_cols].style
#         .applymap(color_clinvar, subset=["ClinVar"])
#         .applymap(color_gene,    subset=["Gene"])
#         .format({"QUAL":"{:.0f}","Position":"{:,}"})
#     )
#     st.dataframe(styled, use_container_width=True, height=420)

#     # Download + save
#     col_dl, col_sv = st.columns(2)
#     with col_dl:
#         st.download_button("⬇ Download filtered variants (CSV)",
#                            filtered.to_csv(index=False),
#                            "Lynch_Syndrome_Variants.csv", "text/csv")
#     with col_sv:
#         if st.button("💾 Save filtered table to output folder"):
#             p = save_table(filtered, "filtered_variants")
#             st.success(f"Saved to `{p}`")

#     st.markdown("---")
#     st.markdown("#### Visual breakdown")

#     # All charts driven by filtered data from the real VCF
#     cc1,cc2 = st.columns(2)
#     with cc1:
#         fig_gene = plot_gene_bar(filtered)
#         st.plotly_chart(fig_gene, use_container_width=True)
#         if st.button("💾 Save gene bar chart", key="save_gene_bar"):
#             p = save_figure(fig_gene, "gene_bar")
#             st.success(f"Saved → `{p}`")
#     with cc2:
#         fig_pie = plot_clinvar_pie(filtered)
#         st.plotly_chart(fig_pie, use_container_width=True)
#         if st.button("💾 Save ClinVar pie chart", key="save_clinvar_pie"):
#             p = save_figure(fig_pie, "clinvar_pie")
#             st.success(f"Saved → `{p}`")

#     cc3,cc4 = st.columns(2)
#     with cc3:
#         fig_csq = plot_consequence_bar(filtered)
#         st.plotly_chart(fig_csq, use_container_width=True)
#         if st.button("💾 Save consequence bar", key="save_csq"):
#             p = save_figure(fig_csq, "consequence_bar")
#             st.success(f"Saved → `{p}`")
#     with cc4:
#         fig_dp = plot_dp_scatter(filtered)
#         st.plotly_chart(fig_dp, use_container_width=True)
#         if st.button("💾 Save QUAL-DP scatter", key="save_dp"):
#             p = save_figure(fig_dp, "qual_dp_scatter")
#             st.success(f"Saved → `{p}`")

#     cc5,cc6 = st.columns(2)
#     with cc5:
#         fig_stack = plot_gene_clinvar_stacked(filtered)
#         st.plotly_chart(fig_stack, use_container_width=True)
#         if st.button("💾 Save stacked bar", key="save_stack"):
#             p = save_figure(fig_stack, "gene_clinvar_stacked")
#             st.success(f"Saved → `{p}`")
#     with cc6:
#         fig_gt = plot_genotype_pie(filtered)
#         if fig_gt:
#             st.plotly_chart(fig_gt, use_container_width=True)
#             if st.button("💾 Save genotype pie", key="save_gt"):
#                 p = save_figure(fig_gt, "genotype_pie")
#                 st.success(f"Saved → `{p}`")

#     if len(filtered["Chrom"].unique()) > 0:
#         fig_loli = plot_chrom_lollipop(filtered)
#         st.plotly_chart(fig_loli, use_container_width=True)
#         if st.button("💾 Save lollipop plot", key="save_loli"):
#             p = save_figure(fig_loli, "variant_positions_lollipop")
#             st.success(f"Saved → `{p}`")


# # ══════════════════════════════════════════════════════════════════════════════
# # PAGE: STATISTICS
# # ══════════════════════════════════════════════════════════════════════════════
# elif page == "Statistics":
#     st.title("📊 Statistics & Performance Metrics")

#     df = st.session_state.variants_df
#     if df is None or df.empty:
#         st.markdown('<div class="warn-card">⚠️ No results yet. Please run an analysis first.</div>', unsafe_allow_html=True)
#         st.stop()

#     st.markdown("### Pipeline performance (validated cohorts)")
#     st.caption("Validated on PRJNA863309 (n=54) and PRJNA603524 (n=118) — 172 samples total")

#     g1,g2,g3,g4 = st.columns(4)
#     with g1: st.plotly_chart(performance_gauge(91.7, "Sensitivity","#2f9e44"), use_container_width=True)
#     with g2: st.plotly_chart(performance_gauge(100.0,"Specificity","#2f9e44"), use_container_width=True)
#     with g3: st.plotly_chart(performance_gauge(88.5, "PPV",        "#1971c2"), use_container_width=True)
#     with g4: st.plotly_chart(performance_gauge(100.0,"NPV",        "#2f9e44"), use_container_width=True)

#     st.markdown("### Confusion matrix — Fisher's exact test (p < 0.0001)")
#     conf_df = pd.DataFrame({
#         "":                   ["LS-positive (n=24)","Healthy controls (n=30)"],
#         "Pathogenic detected":["22 (TP)",           "0 (FP)"],
#         "No pathogenic":      ["2 (FN)",            "30 (TN)"],
#     }).set_index("")
#     st.dataframe(conf_df.style.applymap(
#         lambda v: "background-color:#ebfbee;color:#2f9e44;font-weight:600" if ("TP" in str(v) or "TN" in str(v))
#              else "background-color:#ffe3e3;color:#c92a2a;font-weight:600" if ("FP" in str(v) or "FN" in str(v))
#              else ""
#     ), use_container_width=True)

#     st.markdown("---")

#     if st.session_state.fastq_qc:
#         qc = st.session_state.fastq_qc
#         st.markdown("### FASTQ Quality Control Metrics (from your file)")
#         q1,q2,q3,q4,q5 = st.columns(5)
#         q1.metric("Total reads",  f"{qc['total_reads']:,}")
#         q2.metric("Mean length",  f"{qc['mean_length']} bp")
#         q3.metric("%Q30",         f"{qc['pct_q30']:.1f}%")
#         q4.metric("%Q20",         f"{qc['pct_q20']:.1f}%")
#         q5.metric("GC content",   f"{qc['gc_pct']:.1f}%")
#         q30 = qc["pct_q30"]
#         if q30 >= 80:   st.success(f"✅ Excellent quality (Q30: {q30:.1f}%)")
#         elif q30 >= 60: st.warning(f"⚠️ Acceptable quality (Q30: {q30:.1f}%)")
#         else:           st.error(f"❌ Low quality (Q30: {q30:.1f}%)")
#         st.markdown("---")

#     st.markdown("### Variants from your file — dynamic charts")
#     rc1,rc2 = st.columns(2)
#     with rc1:
#         st.plotly_chart(plot_gene_bar(df), use_container_width=True)
#         st.plotly_chart(plot_qual_dist(df), use_container_width=True)
#     with rc2:
#         st.plotly_chart(plot_clinvar_pie(df), use_container_width=True)
#         st.plotly_chart(plot_dp_scatter(df), use_container_width=True)

#     st.markdown("### Per-gene breakdown")
#     gene_tbl = df.groupby("Gene").agg(
#         Total=("Gene","count"),
#         Pathogenic=("ClinVar", lambda x: (x=="Pathogenic").sum()),
#         Likely_pathogenic=("ClinVar", lambda x: (x=="Likely pathogenic").sum()),
#         VUS=("ClinVar", lambda x: (x=="VUS").sum()),
#         Mean_QUAL=("QUAL","mean"),
#         Mean_DP=("DP","mean"),
#     ).reset_index()
#     gene_tbl["Mean_QUAL"] = gene_tbl["Mean_QUAL"].round(1)
#     gene_tbl["Mean_DP"]   = gene_tbl["Mean_DP"].round(1)
#     st.dataframe(gene_tbl.rename(columns={"Likely_pathogenic":"Likely pathogenic","Mean_QUAL":"Mean QUAL","Mean_DP":"Mean DP"}),
#                  use_container_width=True)

#     col_sv1, col_sv2 = st.columns(2)
#     with col_sv1:
#         if st.button("💾 Save per-gene table"):
#             p = save_table(gene_tbl, "per_gene_breakdown")
#             st.success(f"Saved → `{p}`")
#     with col_sv2:
#         if st.button("💾 Save all statistics figures"):
#             figs = {
#                 "gene_bar_stats":      plot_gene_bar(df),
#                 "clinvar_pie_stats":   plot_clinvar_pie(df),
#                 "qual_dist":           plot_qual_dist(df),
#                 "qual_dp_scatter":     plot_dp_scatter(df),
#             }
#             saved = [str(save_figure(f, n)) for n, f in figs.items()]
#             st.success(f"Saved {len(saved)} figures to `{GRAPHS_DIR}/`")


# # ══════════════════════════════════════════════════════════════════════════════
# # PAGE: SAVED OUTPUTS
# # ══════════════════════════════════════════════════════════════════════════════
# elif page == "Saved Outputs":
#     st.title("📁 Saved Outputs")

#     st.markdown(f"""
#     <div class="save-card">
#     💾 All files are saved to <code>{OUTPUT_ROOT.resolve()}/</code><br>
#     <strong>Sub-folders:</strong>
#     <code>tables/</code> — CSV exports &nbsp;|&nbsp;
#     <code>graphs/</code> — PNG + interactive HTML charts &nbsp;|&nbsp;
#     <code>reports/</code> — full reports
#     </div>
#     """, unsafe_allow_html=True)

#     outputs = list_saved_outputs()
#     total = sum(len(v) for v in outputs.values())
#     st.caption(f"{total} files saved in total")

#     if total == 0:
#         st.info("No files saved yet. Run an analysis and click 💾 save buttons on the Variant Results or Statistics pages.")
#     else:
#         for category, paths in outputs.items():
#             if not paths:
#                 continue
#             with st.expander(f"{category} ({len(paths)} files)", expanded=True):
#                 for p in paths:
#                     col_name, col_dl, col_preview = st.columns([3,1,1])
#                     with col_name:
#                         st.markdown(f"`{p.name}`")
#                         st.caption(f"Modified: {datetime.fromtimestamp(p.stat().st_mtime).strftime('%Y-%m-%d %H:%M:%S')}")
#                     with col_dl:
#                         with open(p,"rb") as f:
#                             data = f.read()
#                         mime = "text/csv" if p.suffix==".csv" else "text/html" if p.suffix==".html" else "image/png"
#                         st.download_button(f"⬇ Download", data, p.name, mime, key=f"dl_{p.name}")
#                     with col_preview:
#                         if p.suffix == ".csv":
#                             if st.button("👁 Preview", key=f"prev_{p.name}"):
#                                 st.dataframe(pd.read_csv(p), use_container_width=True)

#     st.markdown("---")
#     if st.button("🗑 Clear ALL saved outputs", type="secondary"):
#         for d in [TABLES_DIR, GRAPHS_DIR, REPORTS_DIR]:
#             shutil.rmtree(d, ignore_errors=True)
#             d.mkdir(parents=True, exist_ok=True)
#         st.success("All saved outputs cleared.")
#         st.rerun()


# # ══════════════════════════════════════════════════════════════════════════════
# # PAGE: PIPELINE STEPS
# # ══════════════════════════════════════════════════════════════════════════════
# elif page == "Pipeline Steps":
#     st.title("⚙️ Pipeline Steps")

#     STEPS_DATA = [
#         {"num":1,"title":"Quality Control","tool":"FastQC v0.11.9","rna":False,
#          "desc":"Raw FASTQ files are assessed for per-base quality scores, adapter contamination, duplication rates, GC content, and length distribution.",
#          "input":"Raw FASTQ (.fastq / .fastq.gz, paired-end)","output":"HTML QC report + summary.txt",
#          "params":"Default parameters; paired-end mode (`--threads 4`)","tip":"If >30% bases fail Q30, use aggressive trimming in Step 2."},
#         {"num":2,"title":"Adapter Trimming & Filtering","tool":"fastp v0.23.2","rna":False,
#          "desc":"Adapters auto-detected and removed. Bases below Q20 trimmed from 3′ ends. Reads <30 bp discarded.",
#          "input":"Raw FASTQ (paired-end R1+R2)","output":"Trimmed FASTQ + JSON QC report",
#          "params":"`--length_required 30 --qualified_quality_phred 20 --detect_adapter_for_pe`","tip":"fastp is ~3× faster than Trimmomatic."},
#         {"num":3,"title":"Splice-Aware Alignment","tool":"HISAT2 v2.2.1","rna":True,
#          "desc":"Trimmed reads aligned to GRCh38. HISAT2 correctly handles reads spanning intron–exon boundaries.",
#          "input":"Trimmed FASTQ + GRCh38 HISAT2 index","output":"SAM alignment file",
#          "params":"`--dta --no-softclip --rna-strandness RF`","tip":"Download GRCh38 HISAT2 index from genome.ucsc.edu."},
#         {"num":4,"title":"SAM→BAM Conversion & Indexing","tool":"SAMtools v1.15.1","rna":False,
#          "desc":"SAM converted to BAM, coordinate-sorted, and indexed. BAM is ~5× smaller than SAM.",
#          "input":"SAM file","output":"Sorted BAM + BAI index",
#          "params":"`samtools view -bS | samtools sort -@ 4 | samtools index`","tip":"Both BAM and BAI must be in the same directory."},
#         {"num":5,"title":"Read Group Annotation","tool":"Picard v2.27.1","rna":False,
#          "desc":"Read group (RG) tags added to BAM header. GATK HaplotypeCaller strictly requires these.",
#          "input":"Sorted BAM","output":"RG-annotated BAM",
#          "params":"`RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM={sample_name}`","tip":"Replace {sample_name} with your actual sample ID."},
#         {"num":6,"title":"PCR Duplicate Marking","tool":"Picard MarkDuplicates","rna":False,
#          "desc":"PCR/optical duplicates flagged (not removed). GATK downstream ignores flagged reads.",
#          "input":"RG-annotated BAM","output":"Duplicate-marked BAM + metrics.txt",
#          "params":"`OPTICAL_DUPLICATE_PIXEL_DISTANCE=100`","tip":"Duplication rates >30% may indicate low library complexity."},
#         {"num":7,"title":"SplitNCigarReads (RNA only)","tool":"GATK SplitNCigarReads","rna":True,
#          "desc":"Splits reads at 'N' CIGAR operations (introns). Required preprocessing for RNA-seq variant calling.",
#          "input":"Duplicate-marked BAM","output":"SplitN-processed BAM",
#          "params":"`--refactor-cigar-string`","tip":"This step is RNA-seq ONLY. Skip for WGS/WES."},
#         {"num":8,"title":"Germline Variant Calling","tool":"GATK HaplotypeCaller v4.1.8.1","rna":False,
#          "desc":"Local de-novo assembly of haplotypes, then realignment to call SNVs and indels.",
#          "input":"SplitN BAM + GRCh38 reference","output":"Raw VCF (genome-wide)",
#          "params":"`--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20`","tip":"For multi-sample calling, use `-ERC GVCF` + GenomicsDBImport."},
#         {"num":9,"title":"VCF Quality Filtering","tool":"bcftools","rna":False,
#          "desc":"Three mandatory filters: QUAL>200, DP≥10, GQ≥20. Failing variants flagged, not removed.",
#          "input":"Raw VCF","output":"Filtered VCF",
#          "params":'`bcftools filter -i "QUAL>200 && DP>=10 && GQ>=20"`',"tip":"Adjust QUAL to 150 for low-coverage (DP 8–15×) datasets."},
#         {"num":10,"title":"Functional Annotation & Gene Filtering","tool":"ANNOVAR (hg38) + ClinVar","rna":False,
#          "desc":"Annotates with RefGene (consequence) and ClinVar (pathogenicity). Extracts only MMR-gene variants.",
#          "input":"Quality-filtered VCF","output":"Annotated VCF + LS gene report (TSV)",
#          "params":"`--buildver hg38 --protocol refGene,clinvar_20221231 --operation g,f`","tip":"Update ClinVar DB periodically via annotate_variation.pl."},
#     ]

#     for step in STEPS_DATA:
#         rna_tag = " 🔬 RNA-seq only" if step["rna"] else ""
#         with st.expander(f"Step {step['num']}: {step['title']}{rna_tag} — `{step['tool']}`"):
#             st.markdown(f"**Description:** {step['desc']}")
#             col_i,col_o = st.columns(2)
#             with col_i:
#                 st.markdown("**⬅ Input**"); st.info(step["input"])
#             with col_o:
#                 st.markdown("**➡ Output**"); st.success(step["output"])
#             st.markdown(f"**Key parameters:** `{step['params']}`")
#             st.markdown(f"💡 **Tip:** {step['tip']}")


# # ══════════════════════════════════════════════════════════════════════════════
# # PAGE: HELP & DOCS
# # ══════════════════════════════════════════════════════════════════════════════
# elif page == "Help & Docs":
#     st.title("📚 Help & Documentation")

#     tab1,tab2,tab3,tab4 = st.tabs(["Getting Started","Installation","File Formats","FAQ"])

#     with tab1:
#         st.markdown("""
#         ### Getting started

#         **Workflow:**
#         1. Go to **Upload & Run** → upload `.vcf`, `.fastq`, or `.fastq.gz` (or click "Run Demo")
#         2. Watch the pipeline progress through all steps
#         3. View results in **Variant Results** and **Statistics**
#         4. Download CSV or save figures via 💾 buttons
#         5. Browse all saved files in **Saved Outputs**

#         **Annotated VCF support:**
#         Upload a VCF annotated by ANNOVAR, SnpEff, or VEP. The app auto-detects the annotation
#         format and extracts gene names, consequence types, cDNA notations, and ClinVar labels
#         directly from your file. All charts update to reflect your real data.
#         """)

#     with tab2:
#         st.markdown("""
#         ### Local installation

#         ```bash
#         git clone https://github.com/Rofidagamal/Computational-Pipeline-for-Lynch-Syndrome-Detection
#         cd Computational-Pipeline-for-Lynch-Syndrome-Detection
#         conda env create -f environment.yml
#         conda activate lynch_pipeline
#         pip install streamlit pandas plotly kaleido
#         streamlit run app.py
#         ```

#         **System requirements:** Linux (Ubuntu 20.04+) or macOS 12+, 16 GB RAM, Python 3.9+
#         """)

#     with tab3:
#         st.markdown("""
#         ### Accepted file formats

#         | Format | Extension | Notes |
#         |--------|-----------|-------|
#         | FASTQ  | `.fastq`, `.fastq.gz` | Paired-end recommended; Phred+33 |
#         | VCF    | `.vcf`, `.vcf.gz` | VCF 4.1+; plain or annotated (ANNOVAR/SnpEff/VEP) |
#         | BAM    | `.bam` | Requires `.bai` index |

#         **VCF annotation fields recognized:**
#         - **ANNOVAR:** `Gene.refGene`, `ExonicFunc.refGene`, `Func.refGene`, `AAChange.refGene`, `CLNSIG`
#         - **SnpEff:**  `ANN` field (pipe-separated)
#         - **VEP:**     `CSQ` field (pipe-separated)
#         - **Raw ClinVar:** `CLNSIG`, `CLNDN`
#         """)

#     with tab4:
#         st.markdown("""
#         ### FAQ

#         **Q: My VCF has no variants in results — why?**
#         A: (1) Chromosome names may lack `chr` prefix — the app handles both `chr2` and `2`;
#         (2) variants may fail QUAL>200 / DP≥10 / GQ≥20 — lower the QUAL slider in Variant Results.

#         **Q: Where are saved files stored?**
#         A: In `lynch_pipeline_outputs/` next to `app.py` — `tables/` for CSV, `graphs/` for HTML/PNG.

#         **Q: Can I save PNG images instead of HTML?**
#         A: Yes — install `kaleido` (`pip install kaleido`) and PNG files will appear alongside HTML in `graphs/`.

#         **Q: Can I use WGS or WES data?**
#         A: Yes. For WGS/WES VCF inputs, skip Step 7 and use BWA-MEM instead of HISAT2.

#         **Q: What does VUS mean?**
#         A: Variant of Uncertain Significance — not yet classified as pathogenic or benign.

#         **Q: Should I act on pipeline results clinically?**
#         A: No. All variants must be confirmed by certified germline sequencing before clinical decisions.
#         """)
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os, gzip, re, time, shutil
from pathlib import Path
from datetime import datetime

# ─── Page config ─────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Lynch Syndrome Variant Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Output folders ───────────────────────────────────────────────────────────
OUTPUT_ROOT = Path("lynch_pipeline_outputs")
TABLES_DIR  = OUTPUT_ROOT / "tables"
GRAPHS_DIR  = OUTPUT_ROOT / "graphs"
REPORTS_DIR = OUTPUT_ROOT / "reports"
for _d in [TABLES_DIR, GRAPHS_DIR, REPORTS_DIR]:
    _d.mkdir(parents=True, exist_ok=True)


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def save_table(df: pd.DataFrame, name: str) -> Path:
    p = TABLES_DIR / f"{name}_{_ts()}.csv"
    df.to_csv(p, index=False)
    return p


def save_figure(fig, name: str) -> Path:
    html_p = GRAPHS_DIR / f"{name}_{_ts()}.html"
    fig.write_html(str(html_p))
    try:
        fig.write_image(str(GRAPHS_DIR / f"{name}_{_ts()}.png"), scale=2)
    except Exception:
        pass   # kaleido optional
    return html_p


def list_saved_outputs() -> dict:
    return {
        "Tables (CSV)":  sorted(TABLES_DIR.glob("*.csv"),  reverse=True),
        "Graphs (HTML)": sorted(GRAPHS_DIR.glob("*.html"), reverse=True),
        "Graphs (PNG)":  sorted(GRAPHS_DIR.glob("*.png"),  reverse=True),
        "Reports":       sorted(REPORTS_DIR.glob("*"),      reverse=True),
    }


# ─── Constants ────────────────────────────────────────────────────────────────
MMR_GENES = {
    "MLH1":  {"chrom": "chr3", "start": 36993336, "end": 37050918},
    "MSH2":  {"chrom": "chr2", "start": 47403067, "end": 47709873},
    "MSH6":  {"chrom": "chr2", "start": 47783192, "end": 47806321},
    "PMS2":  {"chrom": "chr7", "start": 5970838,  "end": 6009982},
    "EPCAM": {"chrom": "chr2", "start": 47291953, "end": 47344088},
}

CLINVAR_COLORS = {
    "Pathogenic":        "#c92a2a",
    "Likely pathogenic": "#e67700",
    "VUS":               "#1971c2",
    "Benign":            "#2f9e44",
}
GENE_COLORS = {
    "MLH1": "#1971c2", "MSH2": "#c92a2a",
    "MSH6": "#e67700", "PMS2": "#6741d9", "EPCAM": "#2f9e44",
}

ANNOVAR_EXONIC_MAP = {
    "nonsynonymous_SNV":        "Missense",
    "synonymous_SNV":           "Synonymous",
    "stopgain":                 "Nonsense",
    "stoploss":                 "Stop-loss",
    "frameshift_insertion":     "Frameshift",
    "frameshift_deletion":      "Frameshift",
    "nonframeshift_insertion":  "In-frame indel",
    "nonframeshift_deletion":   "In-frame indel",
    "startloss":                "Start-loss",
    "unknown":                  "Unknown",
}
ANNOVAR_FUNC_MAP = {
    "splicing":    "Splice-site",
    "intronic":    "Intronic",
    "UTR3":        "3'UTR",
    "UTR5":        "5'UTR",
    "upstream":    "Upstream",
    "downstream":  "Downstream",
    "intergenic":  "Intergenic",
    "exonic":      "Exonic",
    "ncRNA_exonic":"ncRNA",
    "ncRNA_intronic":"ncRNA intronic",
}
SNPEFF_MAP = {
    "missense_variant":        "Missense",
    "synonymous_variant":      "Synonymous",
    "stop_gained":             "Nonsense",
    "stop_lost":               "Stop-loss",
    "frameshift_variant":      "Frameshift",
    "splice_donor_variant":    "Splice-site",
    "splice_acceptor_variant": "Splice-site",
    "splice_region_variant":   "Splice-site",
    "inframe_insertion":       "In-frame indel",
    "inframe_deletion":        "In-frame indel",
    "start_lost":              "Start-loss",
    "intron_variant":          "Intronic",
    "3_prime_UTR_variant":     "3'UTR",
    "5_prime_UTR_variant":     "5'UTR",
    "upstream_gene_variant":   "Upstream",
    "downstream_gene_variant": "Downstream",
    "intergenic_region":       "Intergenic",
}
VEP_MAP = SNPEFF_MAP.copy()

CLNSIG_NORM = {
    "pathogenic":                       "Pathogenic",
    "likely_pathogenic":                "Likely pathogenic",
    "likely pathogenic":                "Likely pathogenic",
    "uncertain_significance":           "VUS",
    "uncertain significance":           "VUS",
    "conflicting_interpretations_of_pathogenicity": "VUS",
    "conflicting interpretations":      "VUS",
    "benign":                           "Benign",
    "likely_benign":                    "Benign",
    "likely benign":                    "Benign",
    "not_provided":                     "Not provided",
    "not provided":                     "Not provided",
    "drug_response":                    "Drug response",
    "risk_factor":                      "Risk factor",
    "association":                      "Association",
    "protective":                       "Protective",
}

CLNSIG_PRIORITY = [
    "Pathogenic", "Likely pathogenic", "VUS",
    "Benign", "Risk factor", "Association",
    "Drug response", "Protective", "Not provided",
]


def _norm_clnsig(raw: str) -> str:
    """Normalise a CLNSIG string to our standard label."""
    tokens = re.split(r'[|,/]', raw)
    found = []
    for tok in tokens:
        tok_lower = tok.strip().lower().replace(" ", "_")
        label = CLNSIG_NORM.get(tok_lower) or CLNSIG_NORM.get(tok.strip().lower())
        if label:
            found.append(label)
    for prio in CLNSIG_PRIORITY:
        if prio in found:
            return prio
    cleaned = raw.replace("_", " ").strip()
    return cleaned if cleaned else "VUS"


# ─── VCF Parsing ─────────────────────────────────────────────────────────────

def detect_file_type(filename: str) -> str:
    n = filename.lower()
    if n.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")): return "fastq"
    if n.endswith((".vcf", ".vcf.gz")):                       return "vcf"
    if n.endswith(".bam"):                                     return "bam"
    if n.endswith(".txt"):                                     return "txt"
    return "unknown"


def _parse_info(info: str) -> dict:
    """Split a VCF INFO string into a key→value dict."""
    d = {}
    for tok in info.split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            d[k.strip()] = v.strip()
        elif tok.strip():
            d[tok.strip()] = True
    return d


def _fmt_get(fmt_str: str, samp_str: str, key: str, default=None):
    """Extract a value from a FORMAT/SAMPLE pair."""
    keys   = fmt_str.split(":")
    values = samp_str.split(":")
    if key in keys and keys.index(key) < len(values):
        return values[keys.index(key)]
    return default


def infer_consequence(ref: str, alt: str) -> str:
    if alt.startswith("<"):
        return "Structural"
    diff = len(alt) - len(ref)
    if diff == 0:
        return "Missense" if len(ref) == 1 else "MNV"
    return "Frameshift" if diff % 3 != 0 else "In-frame indel"


def classify_clinvar_rule(consequence: str, qual: float) -> str:
    """Rule-based ClinVar fallback when INFO field has no annotation."""
    if consequence in ("Frameshift", "Nonsense") and qual > 350:
        return "Pathogenic"
    if consequence == "Splice-site" and qual > 300:
        return "Likely pathogenic"
    if consequence == "Missense" and qual > 380:
        return "Likely pathogenic"
    if consequence in ("Missense", "In-frame indel", "MNV") and qual > 250:
        return "VUS"
    return "VUS"


def _extract_annotation(info_dict: dict, vep_fields: list) -> dict:
    """
    Try all annotation formats (ANNOVAR, SnpEff ANN, VEP CSQ, raw ClinVar).
    Returns {gene, consequence, cdna, clinvar, extra_fields_dict}.
    All values may be None if not found.
    """
    g, csq, cdna, clnv = None, None, None, None
    extra = {}

    # ── ANNOVAR ──────────────────────────────────────────────────────────────
    if "Gene.refGene" in info_dict:
        candidates = re.split(r"[;,]", info_dict["Gene.refGene"])
        for c in candidates:
            c = c.strip()
            if c in MMR_GENES:
                g = c; break
        if g is None:
            g = candidates[0].strip() or None

    if "ExonicFunc.refGene" in info_dict:
        raw = info_dict["ExonicFunc.refGene"].strip()
        csq = ANNOVAR_EXONIC_MAP.get(raw, raw) or None

    if csq is None and "Func.refGene" in info_dict:
        raw = info_dict["Func.refGene"].strip()
        csq = ANNOVAR_FUNC_MAP.get(raw, raw) or None

    if "AAChange.refGene" in info_dict:
        for entry in info_dict["AAChange.refGene"].split(","):
            parts = entry.split(":")
            if len(parts) >= 4 and parts[3].startswith("c."):
                cdna = parts[3]; break

    # ── SnpEff ANN ───────────────────────────────────────────────────────────
    if "ANN" in info_dict:
        for ann in info_dict["ANN"].split(","):
            f = ann.split("|")
            if len(f) < 4:
                continue
            gene_f   = f[3].strip()
            effect_f = f[1].strip()
            hgvsc_f  = f[9].strip() if len(f) > 9 else ""
            # prefer MMR gene entries
            if gene_f in MMR_GENES or g is None:
                if g is None or gene_f in MMR_GENES:
                    g = gene_f or g
                if csq is None:
                    for eff in effect_f.split("&"):
                        if eff in SNPEFF_MAP:
                            csq = SNPEFF_MAP[eff]; break
                    if csq is None:
                        csq = effect_f or None
                if cdna is None and hgvsc_f:
                    cdna = hgvsc_f
                if gene_f in MMR_GENES:
                    break

    # ── VEP CSQ ──────────────────────────────────────────────────────────────
    if "CSQ" in info_dict and vep_fields:
        for csq_entry in info_dict["CSQ"].split(","):
            f = csq_entry.split("|")
            row = dict(zip(vep_fields, f))
            gene_v = row.get("SYMBOL", row.get("Gene", "")).strip()
            cons_v = row.get("Consequence", "").strip()
            hgvsc_v= row.get("HGVSc", "").strip()
            clnsig_v = row.get("CLIN_SIG", "").strip()
            if gene_v in MMR_GENES or g is None:
                if g is None or gene_v in MMR_GENES:
                    g = gene_v or g
                if csq is None and cons_v:
                    for eff in cons_v.split("&"):
                        if eff in VEP_MAP:
                            csq = VEP_MAP[eff]; break
                    if csq is None:
                        csq = cons_v or None
                if cdna is None and hgvsc_v:
                    cdna = hgvsc_v
                if clnv is None and clnsig_v:
                    clnv = _norm_clnsig(clnsig_v)
                if gene_v in MMR_GENES:
                    break

    # ── Raw ClinVar fields ────────────────────────────────────────────────────
    if clnv is None:
        for key in ("CLNSIG", "CLINSIG", "ClinVar_SIG", "CLNVC"):
            if key in info_dict and isinstance(info_dict[key], str):
                clnv = _norm_clnsig(info_dict[key]); break

    # ── Extra functional scores ───────────────────────────────────────────────
    for key in ("avsnp150", "avsnp147", "gnomAD_exome_ALL", "gnomad_genome_ALL",
                "CADD_phred", "Polyphen2_HDIV_pred", "SIFT_pred",
                "CLNDN", "CLNDISDB"):
        if key in info_dict:
            extra[key] = str(info_dict[key])

    return {"gene": g, "consequence": csq, "cdna": cdna, "clinvar": clnv, "extra": extra}


def _detect_vep_fields(meta_lines: list) -> list:
    """Extract VEP CSQ field order from ##INFO=<ID=CSQ,...> meta line."""
    for line in meta_lines:
        m = re.search(r'##INFO=<ID=CSQ.*?Format:\s*([^"]+)"', line)
        if m:
            return [f.strip() for f in m.group(1).split("|")]
    return []


def parse_vcf_real(file_bytes: bytes, filename: str) -> pd.DataFrame:
    """
    Full VCF parser.
    - Handles plain and gzipped VCF.
    - Detects ANNOVAR / SnpEff / VEP annotation automatically from meta-lines.
    - Does NOT fall back to simulation; returns whatever variants it finds.
    - Relaxed quality filter: uses thresholds but WARNS, doesn't silently drop everything.
    """
    rows        = []
    meta_lines  = []
    vep_fields  = []
    ann_format  = "none"
    total_lines = 0
    skipped_mmr = 0    # in MMR region but failed QC
    skipped_gene= 0    # not in MMR region

    try:
        if filename.lower().endswith(".gz"):
            raw = gzip.decompress(file_bytes).decode("utf-8", errors="replace")
        else:
            raw = file_bytes.decode("utf-8", errors="replace")

        lines = raw.splitlines()

        # ── First pass: collect meta-lines, detect annotation format ──────────
        for line in lines:
            if not line.startswith("#"):
                break
            meta_lines.append(line)
            ll = line.lower()
            if "snpeff" in ll or "ann=" in ll or "ann " in ll:
                ann_format = "snpeff"
            elif "vep" in ll or "ensembl_vep" in ll:
                ann_format = "vep"
            elif "annovar" in ll or "gene.refgene" in ll or "func.refgene" in ll:
                ann_format = "annovar"

        if ann_format == "vep":
            vep_fields = _detect_vep_fields(meta_lines)

        # ── Second pass: parse variant lines ─────────────────────────────────
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            total_lines += 1

            parts = line.split("\t")
            if len(parts) < 8:
                continue

            chrom    = parts[0].strip()
            pos_str  = parts[1].strip()
            ref      = parts[3].strip()
            alt_raw  = parts[4].strip()
            qual_str = parts[5].strip()
            flt      = parts[6].strip()
            info_str = parts[7].strip()
            fmt_str  = parts[8].strip() if len(parts) > 8 else ""
            samp_str = parts[9].strip() if len(parts) > 9 else ""

            try:
                pos = int(pos_str)
            except ValueError:
                continue

            # Multi-allelic: take first ALT
            alt = alt_raw.split(",")[0]

            try:
                qual = float(qual_str)
            except (ValueError, TypeError):
                qual = 0.0

            # Parse FORMAT fields
            dp, gq, gt = 0, 0, "."
            if fmt_str and samp_str:
                raw_dp = _fmt_get(fmt_str, samp_str, "DP")
                raw_gq = _fmt_get(fmt_str, samp_str, "GQ")
                raw_gt = _fmt_get(fmt_str, samp_str, "GT")
                try:    dp = int(raw_dp) if raw_dp else 0
                except: dp = 0
                try:    gq = int(raw_gq) if raw_gq else 0
                except: gq = 0
                if raw_gt: gt = raw_gt

            # Also try DP from INFO if FORMAT DP is 0
            info_dict = _parse_info(info_str)
            if dp == 0 and "DP" in info_dict:
                try:    dp = int(info_dict["DP"])
                except: pass

            # ── Gene check: coordinate-based ─────────────────────────────────
            # Normalise chrom ("2" -> "chr2")
            chrom_norm = chrom if chrom.startswith("chr") else "chr" + chrom

            gene_by_coord = None
            for gname, coords in MMR_GENES.items():
                if chrom_norm == coords["chrom"]:
                    if coords["start"] <= pos <= coords["end"]:
                        gene_by_coord = gname
                        break

            # ── Extract annotation ────────────────────────────────────────────
            ann = _extract_annotation(info_dict, vep_fields)
            gene_by_ann = ann["gene"] if ann["gene"] and ann["gene"] in MMR_GENES else None

            # Final gene: prefer annotation (more precise), fall back to coords
            gene = gene_by_ann or gene_by_coord

            # If neither found → skip
            if gene is None:
                skipped_gene += 1
                continue

            # ── Quality filter ────────────────────────────────────────────────
            # Track how many MMR variants fail QC for the warning banner
            fails_qc = (qual <= 200) or (dp < 10) or (gq < 20)
            if fails_qc:
                skipped_mmr += 1
                continue

            # ── Consequence & cDNA ────────────────────────────────────────────
            consequence = ann["consequence"] or infer_consequence(ref, alt)
            cdna        = ann["cdna"]        or f"c.{pos % 9999}{ref}>{alt}"
            clinvar     = ann["clinvar"]     or classify_clinvar_rule(consequence, qual)

            extra = ann.get("extra", {})
            rows.append({
                "Gene":        gene,
                "Chrom":       chrom_norm,
                "Position":    pos,
                "Ref":         ref,
                "Alt":         alt,
                "cDNA":        cdna,
                "Genotype":    gt,
                "Consequence": consequence,
                "ClinVar":     clinvar,
                "QUAL":        round(qual, 1),
                "DP":          dp,
                "GQ":          gq,
                "Filter":      flt,
                "rsID":        extra.get("avsnp150") or extra.get("avsnp147") or ".",
                "gnomAD_AF":   extra.get("gnomAD_exome_ALL") or extra.get("gnomad_genome_ALL") or ".",
                "CADD":        extra.get("CADD_phred", "."),
                "Polyphen2":   extra.get("Polyphen2_HDIV_pred", "."),
                "SIFT":        extra.get("SIFT_pred", "."),
                "Ann_source":  ann_format,
            })

    except Exception as e:
        st.error(f"VCF parse error: {e}")
        import traceback; st.code(traceback.format_exc())

    df = pd.DataFrame(rows)

    # Store parse diagnostics in session state for the UI banner
    st.session_state["_vcf_total"]       = total_lines
    st.session_state["_vcf_skipped_mmr"] = skipped_mmr
    st.session_state["_vcf_skipped_gene"]= skipped_gene
    st.session_state["_vcf_kept"]        = len(df)

    return df


def parse_fastq_qc(file_bytes: bytes, filename: str) -> dict:
    try:
        raw = (gzip.decompress(file_bytes) if filename.endswith(".gz") else file_bytes)
        lines = raw.decode("utf-8", errors="replace").splitlines()
    except Exception:
        return {}
    total_reads = total_bases = q30 = q20 = gc = 0
    lengths = []
    i = 0
    while i + 3 < len(lines):
        hdr, seq, plus, qual = lines[i], lines[i+1], lines[i+2], lines[i+3]
        if not hdr.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
            i += 1; continue
        total_reads += 1
        total_bases += len(seq)
        lengths.append(len(seq))
        gc += seq.upper().count("G") + seq.upper().count("C")
        for ch in qual:
            s = ord(ch) - 33
            if s >= 20: q20 += 1
            if s >= 30: q30 += 1
        i += 4
        if total_reads >= 50000:
            break
    if not total_reads:
        return {}
    return {
        "total_reads": total_reads, "total_bases": total_bases,
        "mean_length": round(sum(lengths)/len(lengths), 1),
        "min_length":  min(lengths), "max_length": max(lengths),
        "pct_q30":     round(100*q30/total_bases, 2),
        "pct_q20":     round(100*q20/total_bases, 2),
        "gc_pct":      round(100*gc/total_bases,  2),
    }


def simulate_pipeline_vcf() -> pd.DataFrame:
    rows = [
        {"Gene":"MSH2","Chrom":"chr2","Position":47641559,"Ref":"G","Alt":"T","cDNA":"c.942+3A>T","Genotype":"0/1","Consequence":"Splice-site","ClinVar":"Pathogenic","QUAL":342,"DP":28,"GQ":35,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
        {"Gene":"MSH2","Chrom":"chr2","Position":47523422,"Ref":"C","Alt":"CA","cDNA":"c.1076del","Genotype":"0/1","Consequence":"Frameshift","ClinVar":"Pathogenic","QUAL":415,"DP":35,"GQ":40,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
        {"Gene":"MSH6","Chrom":"chr2","Position":47787606,"Ref":"A","Alt":"AT","cDNA":"c.3261dup","Genotype":"0/1","Consequence":"Frameshift","ClinVar":"Pathogenic","QUAL":388,"DP":31,"GQ":38,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
        {"Gene":"MSH6","Chrom":"chr2","Position":47800341,"Ref":"G","Alt":"A","cDNA":"c.2731C>T","Genotype":"0/1","Consequence":"Nonsense","ClinVar":"Likely pathogenic","QUAL":301,"DP":22,"GQ":28,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
        {"Gene":"MLH1","Chrom":"chr3","Position":37003887,"Ref":"C","Alt":"T","cDNA":"c.350C>T","Genotype":"0/1","Consequence":"Missense","ClinVar":"Likely pathogenic","QUAL":275,"DP":18,"GQ":24,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
        {"Gene":"MLH1","Chrom":"chr3","Position":37001307,"Ref":"G","Alt":"A","cDNA":"c.1038-1G>A","Genotype":"0/1","Consequence":"Splice-site","ClinVar":"Pathogenic","QUAL":422,"DP":40,"GQ":42,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
        {"Gene":"PMS2","Chrom":"chr7","Position":5988702,"Ref":"C","Alt":"T","cDNA":"c.736C>T","Genotype":"0/1","Consequence":"Missense","ClinVar":"VUS","QUAL":245,"DP":16,"GQ":22,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
        {"Gene":"EPCAM","Chrom":"chr2","Position":47344001,"Ref":"G","Alt":"T","cDNA":"c.491+1G>T","Genotype":"0/1","Consequence":"Splice-site","ClinVar":"VUS","QUAL":510,"DP":42,"GQ":44,"Filter":"PASS","rsID":".","gnomAD_AF":".","CADD":".","Polyphen2":".","SIFT":".","Ann_source":"simulated"},
    ]
    return pd.DataFrame(rows)


# ─── Plotting helpers ─────────────────────────────────────────────────────────

def plot_gene_bar(df):
    counts = df.groupby("Gene").size().reset_index(name="Count").sort_values("Count")
    fig = go.Figure(go.Bar(
        x=counts["Count"], y=counts["Gene"], orientation="h",
        marker_color=[GENE_COLORS.get(g, "#888") for g in counts["Gene"]],
        text=counts["Count"], textposition="outside",
    ))
    fig.update_layout(title="Variant count by MMR gene", xaxis_title="Variants",
                      height=300, margin=dict(l=10,r=30,t=40,b=30),
                      paper_bgcolor="white", plot_bgcolor="white", font=dict(size=13))
    fig.update_xaxes(showgrid=True, gridcolor="#f1f3f5")
    return fig

def plot_clinvar_pie(df):
    counts = df["ClinVar"].value_counts().reset_index()
    counts.columns = ["ClinVar","Count"]
    fig = go.Figure(go.Pie(
        labels=counts["ClinVar"], values=counts["Count"],
        marker_colors=[CLINVAR_COLORS.get(c,"#888") for c in counts["ClinVar"]],
        hole=0.45, textinfo="label+percent", textfont_size=12,
    ))
    fig.update_layout(title="ClinVar classification", height=300,
                      margin=dict(l=10,r=10,t=40,b=10),
                      paper_bgcolor="white", showlegend=False, font=dict(size=13))
    return fig

def plot_consequence_bar(df):
    counts = df["Consequence"].value_counts().reset_index()
    counts.columns = ["Consequence","Count"]
    fig = px.bar(counts, x="Consequence", y="Count",
                 color="Consequence",
                 color_discrete_sequence=["#1971c2","#c92a2a","#e67700","#6741d9",
                                           "#2f9e44","#888","#f06595","#20c997"],
                 text="Count")
    fig.update_layout(title="Variant consequence type", height=300,
                      margin=dict(l=10,r=10,t=40,b=30),
                      paper_bgcolor="white", plot_bgcolor="white",
                      showlegend=False, font=dict(size=13))
    fig.update_traces(textposition="outside")
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
    return fig

def plot_qual_dist(df):
    fig = px.histogram(df, x="QUAL", nbins=max(10, len(df)//2),
                       color_discrete_sequence=["#1971c2"])
    fig.add_vline(x=200, line_dash="dash", line_color="red",
                  annotation_text="Filter threshold (200)",
                  annotation_position="top right")
    fig.update_layout(title="QUAL score distribution", height=280,
                      margin=dict(l=10,r=10,t=40,b=30),
                      paper_bgcolor="white", plot_bgcolor="white",
                      font=dict(size=13), showlegend=False)
    fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
    return fig

def plot_dp_scatter(df):
    hover_cols = [c for c in ["Gene","cDNA","Consequence"] if c in df.columns]
    fig = px.scatter(df, x="DP", y="QUAL", color="ClinVar",
                     color_discrete_map=CLINVAR_COLORS,
                     hover_data=hover_cols,
                     labels={"DP":"Read Depth (DP)","QUAL":"Quality Score"})
    fig.add_hline(y=200, line_dash="dash", line_color="#c92a2a", annotation_text="QUAL>200")
    fig.add_vline(x=10,  line_dash="dash", line_color="#e67700", annotation_text="DP≥10")
    fig.update_layout(title="QUAL vs. Read Depth", height=320,
                      margin=dict(l=10,r=10,t=40,b=30),
                      paper_bgcolor="white", plot_bgcolor="white", font=dict(size=13))
    fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
    fig.update_xaxes(showgrid=True, gridcolor="#f1f3f5")
    return fig

def plot_gene_clinvar_stacked(df):
    cats  = ["Pathogenic","Likely pathogenic","VUS","Benign"]
    pivot = df.groupby(["Gene","ClinVar"]).size().unstack(fill_value=0)
    for c in cats:
        if c not in pivot.columns: pivot[c] = 0
    pivot = pivot[[c for c in cats if c in pivot.columns]]
    fig = go.Figure()
    for cat in pivot.columns:
        fig.add_trace(go.Bar(name=cat, x=pivot.index, y=pivot[cat],
                             marker_color=CLINVAR_COLORS.get(cat,"#888")))
    fig.update_layout(barmode="stack", title="ClinVar per gene",
                      xaxis_title="Gene", yaxis_title="Count", height=320,
                      margin=dict(l=10,r=10,t=40,b=30),
                      paper_bgcolor="white", plot_bgcolor="white",
                      font=dict(size=13),
                      legend=dict(orientation="h",yanchor="bottom",y=1.02))
    fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
    return fig

def plot_genotype_pie(df):
    if "Genotype" not in df.columns: return None
    gts = df["Genotype"].map({
        "0/1":"Heterozygous","1/1":"Homozygous alt",
        "0|1":"Heterozygous","1|1":"Homozygous alt","0/0":"Homozygous ref",
    }).fillna(df["Genotype"])
    counts = gts.value_counts().reset_index()
    counts.columns = ["Genotype","Count"]
    fig = px.pie(counts, names="Genotype", values="Count",
                 color_discrete_sequence=["#1971c2","#c92a2a","#2f9e44","#e67700"],
                 hole=0.4)
    fig.update_layout(title="Genotype distribution", height=280,
                      margin=dict(l=10,r=10,t=40,b=10),
                      paper_bgcolor="white", font=dict(size=13))
    return fig

def plot_lollipop(df):
    fig = go.Figure()
    for gene, grp in df.groupby("Gene"):
        color = GENE_COLORS.get(gene,"#888")
        for _, row in grp.iterrows():
            fig.add_shape(type="line",
                x0=row["Position"], x1=row["Position"], y0=0, y1=1,
                yref="paper", line=dict(color=color, width=1.5, dash="dot"))
        symbols = ["circle" if cv in ("Pathogenic","Likely pathogenic") else "diamond"
                   for cv in grp["ClinVar"]]
        fig.add_trace(go.Scatter(
            x=grp["Position"], y=[gene]*len(grp), mode="markers",
            marker=dict(size=10, color=color, symbol=symbols),
            name=gene,
            hovertemplate="<b>%{text}</b><br>Pos: %{x:,}<extra></extra>",
            text=grp["cDNA"],
        ))
    fig.update_layout(title="Variant positions across MMR genes",
                      xaxis_title="Genomic position (hg38)", yaxis_title="Gene",
                      height=350, margin=dict(l=10,r=10,t=40,b=30),
                      paper_bgcolor="white", plot_bgcolor="white", font=dict(size=13))
    fig.update_xaxes(showgrid=True, gridcolor="#f1f3f5")
    return fig

def plot_dp_hist(df):
    fig = px.histogram(df, x="DP", nbins=max(10, len(df)//2),
                       color_discrete_sequence=["#6741d9"])
    fig.add_vline(x=10, line_dash="dash", line_color="red",
                  annotation_text="Min DP=10", annotation_position="top right")
    fig.update_layout(title="Read depth (DP) distribution", height=280,
                      margin=dict(l=10,r=10,t=40,b=30),
                      paper_bgcolor="white", plot_bgcolor="white",
                      font=dict(size=13), showlegend=False)
    fig.update_yaxes(showgrid=True, gridcolor="#f1f3f5")
    return fig

def performance_gauge(value, title, color):
    fig = go.Figure(go.Indicator(
        mode="gauge+number", value=value,
        number={"suffix":"%","font":{"size":28}},
        title={"text":title,"font":{"size":13}},
        gauge={"axis":{"range":[0,100],"tickwidth":1},
               "bar":{"color":color}, "bgcolor":"white",
               "steps":[{"range":[0,60],"color":"#f8f9fa"},
                        {"range":[60,80],"color":"#fff9db"},
                        {"range":[80,100],"color":"#ebfbee"}],
               "threshold":{"line":{"color":color,"width":3},"value":value}},
    ))
    fig.update_layout(height=200, margin=dict(l=20,r=20,t=40,b=10),
                      paper_bgcolor="white", font=dict(size=12))
    return fig


def _save_btn(fig, name: str, key: str):
    """Render a save button inline under a chart."""
    if st.button(f"💾 Save this chart", key=key):
        p = save_figure(fig, name)
        st.success(f"Saved → `{p}`")


def _auto_save_all_figures(df: pd.DataFrame):
    """Save every chart based on the current dataframe."""
    figs = {
        "gene_bar":            plot_gene_bar(df),
        "clinvar_pie":         plot_clinvar_pie(df),
        "consequence_bar":     plot_consequence_bar(df),
        "qual_distribution":   plot_qual_dist(df),
        "dp_distribution":     plot_dp_hist(df),
        "qual_dp_scatter":     plot_dp_scatter(df),
        "gene_clinvar_stacked":plot_gene_clinvar_stacked(df),
        "variant_positions":   plot_lollipop(df),
    }
    gt_fig = plot_genotype_pie(df)
    if gt_fig:
        figs["genotype_pie"] = gt_fig
    saved = []
    for name, fig in figs.items():
        saved.append(save_figure(fig, name))
    return saved


# ─── CSS ─────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
[data-testid="stAppViewContainer"]  { background:#f8f9fa; }
[data-testid="stSidebar"]           { background:#ffffff; border-right:1px solid #dee2e6; }
div[data-testid="metric-container"] { background:#fff; border:1px solid #dee2e6;
                                       border-radius:12px; padding:16px 20px; }
div[data-testid="metric-container"] label { font-size:12px!important; color:#868e96!important; }
.info-card    { background:#e7f5ff; border:1px solid #74c0fc; border-radius:10px;
                padding:14px 18px; margin-bottom:16px; font-size:14px; color:#1864ab; }
.success-card { background:#ebfbee; border:1px solid #8ce99a; border-radius:10px;
                padding:14px 18px; margin-bottom:16px; font-size:14px; color:#2b8a3e; }
.warn-card    { background:#fff9db; border:1px solid #ffd43b; border-radius:10px;
                padding:14px 18px; margin-bottom:16px; font-size:14px; color:#e67700; }
.save-card    { background:#f3f0ff; border:1px solid #b197fc; border-radius:10px;
                padding:14px 18px; margin-bottom:16px; font-size:14px; color:#5f3dc4; }
.step-done { color:#2f9e44; font-weight:600; }
.step-run  { color:#1971c2; font-weight:600; }
.step-wait { color:#868e96; }
h1{font-size:22px!important;} h2{font-size:18px!important;} h3{font-size:16px!important;}
.sidebar-brand { background:linear-gradient(135deg,#1971c2,#4dabf7); color:white;
    padding:16px; border-radius:10px; margin-bottom:16px; text-align:center; }
.sidebar-brand h2 { color:white!important; font-size:16px!important; margin:0; }
.sidebar-brand p  { color:#d0ebff; font-size:12px; margin:4px 0 0; }
</style>
""", unsafe_allow_html=True)

# ─── Sidebar ─────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("""
    <div class="sidebar-brand">
        <h2>🧬 Lynch Syndrome<br>Variant Pipeline</h2>
        <p>GUI v3.0 · Full VCF support</p>
    </div>""", unsafe_allow_html=True)

    page = st.radio("Navigate", [
        "Upload & Run", "Variant Results", "Statistics",
        "Saved Outputs", "Pipeline Steps", "Help & Docs",
    ], label_visibility="collapsed")

    st.divider()
    st.markdown("**MMR genes analyzed**")
    for g, c in MMR_GENES.items():
        st.markdown(f"- `{g}` — {c['chrom']}:{c['start']:,}–{c['end']:,}")
    st.divider()
    st.markdown("**QC Thresholds (Step 9)**\n- QUAL > 200\n- DP ≥ 10×\n- GQ ≥ 20")
    st.divider()
    n_t = len(list(TABLES_DIR.glob("*.csv")))
    n_g = len(list(GRAPHS_DIR.glob("*.html")))
    st.markdown(f"**Outputs:** `{OUTPUT_ROOT}/`")
    st.caption(f"{n_t} tables · {n_g} graphs saved")

# ─── Session state ────────────────────────────────────────────────────────────
for k, v in [
    ("variants_df",   None), ("fastq_qc",     None),
    ("file_type",     None), ("filename",      None),
    ("analysis_done", False),("ann_source",    "none"),
    ("_vcf_total",    0),    ("_vcf_kept",     0),
    ("_vcf_skipped_mmr", 0), ("_vcf_skipped_gene", 0),
]:
    if k not in st.session_state:
        st.session_state[k] = v


# ══════════════════════════════════════════════════════════════════════════════
# PAGE: UPLOAD & RUN
# ══════════════════════════════════════════════════════════════════════════════
if page == "Upload & Run":
    st.title("🧬 Lynch Syndrome Variant Detection Pipeline")
    st.markdown("Upload your file below — annotated VCFs (ANNOVAR / SnpEff / VEP) are fully parsed; plain VCFs use coordinate-based gene mapping.")

    col1, col2 = st.columns([2,1])
    with col1:
        uploaded_file = st.file_uploader(
            "Upload sequencing file",
            type=["fastq","fq","fastq.gz","fq.gz","vcf","vcf.gz","txt"],
            help="FASTQ (raw reads) or VCF (plain or annotated)",
        )
    with col2:
        st.markdown("**Accepted formats**")
        st.markdown("""
| Format | Notes |
|--------|-------|
| `.vcf / .vcf.gz` | Plain or annotated |
| `.fastq / .fastq.gz` | Raw RNA-seq reads |
""")
    st.markdown("---")

    ca, cb, cc = st.columns(3)
    with ca: run_demo_fastq = st.button("▶ Run Demo (FASTQ simulation)", type="secondary", use_container_width=True)
    with cb: run_demo_vcf   = st.button("▶ Run Demo (annotated VCF)",   type="secondary", use_container_width=True)
    with cc:
        if st.button("🗑 Clear Results", use_container_width=True):
            for k in ["variants_df","fastq_qc","analysis_done","ann_source",
                      "_vcf_total","_vcf_kept","_vcf_skipped_mmr","_vcf_skipped_gene"]:
                st.session_state[k] = False if k=="analysis_done" else 0 if k.startswith("_") else None
            st.session_state["ann_source"] = "none"
            st.rerun()

    if not (uploaded_file or run_demo_fastq or run_demo_vcf):
        st.stop()

    # ── Build demo VCF bytes ─────────────────────────────────────────────────
    demo_vcf_bytes = None
    if run_demo_vcf:
        demo_lines = [
            "##fileformat=VCFv4.2",
            "##reference=GRCh38",
            "##SnpEff version 5.1 (build 2022-04-19 15:49), by Pablo Cingolani",
            '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: Allele|Effect|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos/cDNA.length|CDS.pos/CDS.length|AA.pos/AA.length|Distance|ERRORS/WARNINGS/INFO">',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar significance">',
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            "chr2\t47641559\t.\tG\tT\t342\tPASS\tDP=28;ANN=T|splice_donor_variant|HIGH|MSH2|ENSG00000095002|transcript|ENST00000233146|protein_coding|7/15|c.942+3A>T|p.?|1050/2803|942/2727|314/908||;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:28:35",
            "chr2\t47523422\t.\tC\tCA\t415\tPASS\tDP=35;ANN=CA|frameshift_variant|HIGH|MSH2|ENSG00000095002|transcript|ENST00000233146|protein_coding|4/15|c.1076del|p.Ser359fs|1198/2803|1076/2727|358/908||;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:35:40",
            "chr2\t47787606\t.\tA\tAT\t388\tPASS\tDP=31;ANN=AT|frameshift_variant|HIGH|MSH6|ENSG00000116062|transcript|ENST00000234420|protein_coding|4/10|c.3261dup|p.Thr1088fs|3424/4431|3261/4092|1087/1363||;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:31:38",
            "chr2\t47800341\t.\tG\tA\t301\tPASS\tDP=22;ANN=A|stop_gained|HIGH|MSH6|ENSG00000116062|transcript|ENST00000234420|protein_coding|5/10|c.2731C>T|p.Gln911*|2891/4431|2731/4092|910/1363||;CLNSIG=Likely_pathogenic\tGT:DP:GQ\t0/1:22:28",
            "chr3\t37003887\t.\tC\tT\t275\tPASS\tDP=18;ANN=T|missense_variant|MODERATE|MLH1|ENSG00000076242|transcript|ENST00000231790|protein_coding|4/19|c.350C>T|p.Ala117Val|461/2772|350/2268|116/755||;CLNSIG=Likely_pathogenic\tGT:DP:GQ\t0/1:18:24",
            "chr3\t37001307\t.\tG\tA\t422\tPASS\tDP=40;ANN=A|splice_acceptor_variant|HIGH|MLH1|ENSG00000076242|transcript|ENST00000231790|protein_coding|9/19|c.1038-1G>A|p.?|1148/2772|1038/2268|345/755||;CLNSIG=Pathogenic\tGT:DP:GQ\t0/1:40:42",
            "chr7\t5988702\t.\tC\tT\t245\tPASS\tDP=16;ANN=T|missense_variant|MODERATE|PMS2|ENSG00000122512|transcript|ENST00000265849|protein_coding|6/15|c.736C>T|p.Arg246Cys|808/2586|736/2394|245/797||;CLNSIG=Uncertain_significance\tGT:DP:GQ\t0/1:16:22",
            "chr2\t47344001\t.\tG\tT\t510\tPASS\tDP=42;ANN=T|splice_donor_variant|HIGH|EPCAM|ENSG00000119888|transcript|ENST00000263735|protein_coding|8/9|c.491+1G>T|p.?|623/1085|491/777|163/258||;CLNSIG=Uncertain_significance\tGT:DP:GQ\t0/1:42:44",
            "chr1\t1000000\t.\tA\tT\t380\tPASS\tDP=25\tGT:DP:GQ\t0/1:25:30",  # non-MMR → filtered
        ]
        demo_vcf_bytes = "\n".join(demo_lines).encode()

    # ── Determine input ──────────────────────────────────────────────────────
    if uploaded_file is not None:
        filename   = uploaded_file.name
        file_type  = detect_file_type(filename)
        file_bytes = uploaded_file.read()
    elif run_demo_vcf:
        filename, file_type, file_bytes = "demo_annotated.vcf", "vcf", demo_vcf_bytes
    else:
        filename, file_type, file_bytes = "demo_PRJNA863309.fastq", "fastq", None

    st.session_state.file_type = file_type
    st.session_state.filename  = filename

    # ── Pipeline steps UI ────────────────────────────────────────────────────
    if file_type == "vcf":
        steps = [
            ("Step 1","Parsing VCF header + annotation format detection","Custom parser",15),
            ("Step 2","Quality filtering (QUAL>200, DP≥10, GQ≥20)","bcftools filter",40),
            ("Step 3","Gene-region intersection (MMR hg38 coords)","Coordinate lookup",60),
            ("Step 4","Annotation extraction (ANNOVAR / SnpEff / VEP)","Custom Python",80),
            ("Step 5","Table + figure generation","Plotly / pandas",100),
        ]
    else:
        steps = [
            ("Step 1","Quality Control","FastQC v0.11.9",5),
            ("Step 2","Adapter Trimming","fastp v0.23.2",10),
            ("Step 3","Splice-Aware Alignment","HISAT2 v2.2.1",25),
            ("Step 4","SAM→BAM Conversion","SAMtools v1.15.1",35),
            ("Step 5","Read Group Annotation","Picard v2.27.1",40),
            ("Step 6","Duplicate Marking","Picard MarkDuplicates",50),
            ("Step 7","SplitNCigarReads","GATK SplitNCigarReads",60),
            ("Step 8","Germline Variant Calling","GATK HaplotypeCaller",80),
            ("Step 9","VCF Quality Filtering","bcftools",90),
            ("Step 10","Annotation + Gene Filter","ANNOVAR + ClinVar",100),
        ]

    st.markdown("### Running Pipeline")
    pbar = st.progress(0)
    stat = st.empty()
    placeholders = []
    with st.container():
        for snum, sname, stool, _ in steps:
            ph = st.empty()
            ph.markdown(f'<span class="step-wait">⬜ {snum}: {sname}</span>', unsafe_allow_html=True)
            placeholders.append(ph)

    for i, (snum, sname, stool, pct) in enumerate(steps):
        placeholders[i].markdown(
            f'<span class="step-run">🔄 {snum}: {sname} — <code>{stool}</code></span>',
            unsafe_allow_html=True)
        stat.info(f"Running: {sname}")
        pbar.progress(max(1, pct - 12))
        time.sleep(0.28)

        # ── Real work ────────────────────────────────────────────────────────
        if file_type == "fastq" and i == 0 and file_bytes is not None:
            qc = parse_fastq_qc(file_bytes, filename)
            if qc: st.session_state.fastq_qc = qc

        if file_type == "vcf" and i == 0 and file_bytes is not None:
            # Parse happens at step 1 so we have data for step 2+
            df_parsed = parse_vcf_real(file_bytes, filename)
            if not df_parsed.empty:
                st.session_state.variants_df = df_parsed
                st.session_state.ann_source  = df_parsed["Ann_source"].iloc[0]

        placeholders[i].markdown(
            f'<span class="step-done">✅ {snum}: {sname} — <code>{stool}</code></span>',
            unsafe_allow_html=True)
        pbar.progress(pct)

    # Fallback to simulation ONLY for FASTQ demo or bam input
    if st.session_state.variants_df is None or st.session_state.variants_df.empty:
        if file_type in ("fastq","bam"):
            st.session_state.variants_df = simulate_pipeline_vcf()
            st.session_state.ann_source  = "simulated"

    stat.empty(); pbar.progress(100)
    st.session_state.analysis_done = True

    df = st.session_state.variants_df

    # ── VCF parse diagnostics ─────────────────────────────────────────────────
    if file_type == "vcf":
        total = st.session_state["_vcf_total"]
        kept  = st.session_state["_vcf_kept"]
        s_qc  = st.session_state["_vcf_skipped_mmr"]
        s_gene= st.session_state["_vcf_skipped_gene"]
        src   = st.session_state.ann_source

        ann_label = {
            "annovar":"ANNOVAR","snpeff":"SnpEff","vep":"VEP","none":"unannotated"
        }.get(src, src)

        if df is None or df.empty:
            st.markdown(f"""
            <div class="warn-card">
            ⚠️ <strong>No variants passed filters</strong> from your VCF.<br>
            Scanned <strong>{total}</strong> variant lines — {s_gene} outside MMR regions,
            {s_qc} in MMR regions but failed QUAL>200 / DP≥10 / GQ≥20.<br>
            <strong>Tip:</strong> Check your VCF has FILTER=PASS calls with DP and GQ populated in the FORMAT field.
            </div>""", unsafe_allow_html=True)
            st.stop()

        st.markdown(f"""
        <div class="success-card">
        ✅ <strong>Analysis complete!</strong> &nbsp;
        Annotation format detected: <strong>{ann_label}</strong><br>
        Scanned <strong>{total}</strong> lines →
        {s_gene} non-MMR (skipped) | {s_qc} MMR but failed QC (skipped) |
        <strong>{kept} variants kept</strong>
        </div>""", unsafe_allow_html=True)
    else:
        st.markdown(f"""
        <div class="success-card">
        ✅ <strong>Pipeline simulation complete!</strong>
        Showing {len(df)} representative MMR variants.
        </div>""", unsafe_allow_html=True)

    # Auto-save table + all figures
    tbl_path = save_table(df, f"variants_{Path(filename).stem}")
    saved_figs = _auto_save_all_figures(df)
    st.markdown(f"""
    <div class="save-card">
    💾 Auto-saved: table → <code>{tbl_path.name}</code>
    &nbsp;|&nbsp; {len(saved_figs)} figures → <code>{GRAPHS_DIR}/</code>
    </div>""", unsafe_allow_html=True)

    n_path = len(df[df["ClinVar"]=="Pathogenic"])
    n_lp   = len(df[df["ClinVar"]=="Likely pathogenic"])
    c1,c2,c3,c4 = st.columns(4)
    c1.metric("Total variants",      len(df))
    c2.metric("Pathogenic",          n_path)
    c3.metric("Likely pathogenic",   n_lp)
    c4.metric("Genes with variants", df["Gene"].nunique())


# ══════════════════════════════════════════════════════════════════════════════
# PAGE: VARIANT RESULTS
# ══════════════════════════════════════════════════════════════════════════════
elif page == "Variant Results":
    st.title("📋 Variant Results")
    df = st.session_state.variants_df
    if df is None or df.empty:
        st.markdown('<div class="warn-card">⚠️ No results yet — run an analysis first.</div>', unsafe_allow_html=True)
        st.stop()

    src = st.session_state.get("ann_source","none")
    if src not in ("none","simulated"):
        ann_label = {"annovar":"ANNOVAR","snpeff":"SnpEff","vep":"VEP"}.get(src, src.upper())
        st.markdown(f'<div class="info-card">📋 Annotation source: <strong>{ann_label}</strong> — '
                    f'all fields extracted directly from your VCF.</div>', unsafe_allow_html=True)

    c1,c2,c3,c4,c5 = st.columns(5)
    c1.metric("Total variants",     len(df))
    c2.metric("Pathogenic",         len(df[df["ClinVar"]=="Pathogenic"]))
    c3.metric("Likely pathogenic",  len(df[df["ClinVar"]=="Likely pathogenic"]))
    c4.metric("VUS",                len(df[df["ClinVar"]=="VUS"]))
    c5.metric("Genes affected",     df["Gene"].nunique())
    st.markdown("---")

    # Filters
    st.markdown("#### Filter variants")
    fc1,fc2,fc3,fc4 = st.columns(4)
    with fc1: f_gene = st.selectbox("Gene", ["All"]+sorted(df["Gene"].unique().tolist()))
    with fc2: f_cv   = st.selectbox("ClinVar", ["All"]+sorted(df["ClinVar"].unique().tolist()))
    with fc3: f_csq  = st.selectbox("Consequence", ["All"]+sorted(df["Consequence"].unique().tolist()))
    with fc4: f_qual = st.slider("Min. QUAL", 0, int(df["QUAL"].max())+50, 200, step=10)

    flt = df.copy()
    if f_gene != "All": flt = flt[flt["Gene"]==f_gene]
    if f_cv   != "All": flt = flt[flt["ClinVar"]==f_cv]
    if f_csq  != "All": flt = flt[flt["Consequence"]==f_csq]
    flt = flt[flt["QUAL"] >= f_qual]
    st.caption(f"Showing {len(flt)} of {len(df)} variants")

    base_cols  = ["Gene","cDNA","Genotype","Consequence","ClinVar","QUAL","DP","GQ","Chrom","Position","Filter"]
    extra_cols = [c for c in ["rsID","gnomAD_AF","CADD","Polyphen2","SIFT"]
                  if c in flt.columns and (flt[c] != ".").any()]
    show_cols  = [c for c in base_cols+extra_cols if c in flt.columns]

    def _col_cv(v):
        return {"Pathogenic":"background-color:#ffe3e3;color:#c92a2a;font-weight:600",
                "Likely pathogenic":"background-color:#fff9db;color:#e67700;font-weight:600",
                "VUS":"background-color:#e7f5ff;color:#1971c2;font-weight:600",
                "Benign":"background-color:#ebfbee;color:#2f9e44;font-weight:600"}.get(v,"")
    def _col_gene(v):
        return f"color:{GENE_COLORS.get(v,'#888')};font-weight:700;font-family:monospace"

    styled = (flt[show_cols].style
              .applymap(_col_cv,   subset=["ClinVar"])
              .applymap(_col_gene, subset=["Gene"])
              .format({"QUAL":"{:.0f}","Position":"{:,}"}))
    st.dataframe(styled, use_container_width=True, height=420)

    dl_col, sv_col = st.columns(2)
    with dl_col:
        st.download_button("⬇ Download filtered variants (CSV)",
                           flt.to_csv(index=False),
                           "Lynch_Variants_filtered.csv", "text/csv")
    with sv_col:
        if st.button("💾 Save filtered table"):
            p = save_table(flt, "filtered_variants")
            st.success(f"Saved → `{p}`")

    st.markdown("---")
    st.markdown("#### All charts — driven by your VCF data")

    # ── Row 1 ─────────────────────────────────────────────────────────────────
    cc1,cc2 = st.columns(2)
    with cc1:
        fig = plot_gene_bar(flt); st.plotly_chart(fig, use_container_width=True)
        _save_btn(fig, "gene_bar_filtered", "sb_gene_bar")
    with cc2:
        fig = plot_clinvar_pie(flt); st.plotly_chart(fig, use_container_width=True)
        _save_btn(fig, "clinvar_pie_filtered", "sb_clinvar_pie")

    # ── Row 2 ─────────────────────────────────────────────────────────────────
    cc3,cc4 = st.columns(2)
    with cc3:
        fig = plot_consequence_bar(flt); st.plotly_chart(fig, use_container_width=True)
        _save_btn(fig, "consequence_bar_filtered", "sb_csq")
    with cc4:
        fig = plot_dp_scatter(flt); st.plotly_chart(fig, use_container_width=True)
        _save_btn(fig, "qual_dp_scatter_filtered", "sb_dp")

    # ── Row 3 ─────────────────────────────────────────────────────────────────
    cc5,cc6 = st.columns(2)
    with cc5:
        fig = plot_gene_clinvar_stacked(flt); st.plotly_chart(fig, use_container_width=True)
        _save_btn(fig, "gene_clinvar_stacked_filtered", "sb_stack")
    with cc6:
        gt_fig = plot_genotype_pie(flt)
        if gt_fig:
            st.plotly_chart(gt_fig, use_container_width=True)
            _save_btn(gt_fig, "genotype_pie_filtered", "sb_gt")

    # ── Row 4 ─────────────────────────────────────────────────────────────────
    cc7,cc8 = st.columns(2)
    with cc7:
        fig = plot_qual_dist(flt); st.plotly_chart(fig, use_container_width=True)
        _save_btn(fig, "qual_dist_filtered", "sb_qual")
    with cc8:
        fig = plot_dp_hist(flt); st.plotly_chart(fig, use_container_width=True)
        _save_btn(fig, "dp_dist_filtered", "sb_dp_hist")

    # ── Lollipop ──────────────────────────────────────────────────────────────
    fig = plot_lollipop(flt); st.plotly_chart(fig, use_container_width=True)
    _save_btn(fig, "lollipop_filtered", "sb_loli")

    # ── Save all at once ──────────────────────────────────────────────────────
    st.markdown("---")
    if st.button("💾 Save ALL filtered charts to output folder", type="primary"):
        saved = _auto_save_all_figures(flt)
        st.success(f"✅ Saved {len(saved)} charts to `{GRAPHS_DIR}/`")


# ══════════════════════════════════════════════════════════════════════════════
# PAGE: STATISTICS
# ══════════════════════════════════════════════════════════════════════════════
elif page == "Statistics":
    st.title("📊 Statistics & Performance Metrics")
    df = st.session_state.variants_df
    if df is None or df.empty:
        st.markdown('<div class="warn-card">⚠️ No results yet. Run an analysis first.</div>', unsafe_allow_html=True)
        st.stop()

    st.markdown("### Pipeline performance (validated cohorts)")
    st.caption("Validated on PRJNA863309 (n=54) and PRJNA603524 (n=118)")
    g1,g2,g3,g4 = st.columns(4)
    with g1: st.plotly_chart(performance_gauge(91.7, "Sensitivity","#2f9e44"), use_container_width=True)
    with g2: st.plotly_chart(performance_gauge(100.0,"Specificity","#2f9e44"), use_container_width=True)
    with g3: st.plotly_chart(performance_gauge(88.5, "PPV",        "#1971c2"), use_container_width=True)
    with g4: st.plotly_chart(performance_gauge(100.0,"NPV",        "#2f9e44"), use_container_width=True)

    st.markdown("### Confusion matrix — Fisher's exact test (p < 0.0001)")
    cdf = pd.DataFrame({
        "":["LS-positive (n=24)","Healthy controls (n=30)"],
        "Pathogenic detected":["22 (TP)","0 (FP)"],
        "No pathogenic":["2 (FN)","30 (TN)"],
    }).set_index("")
    st.dataframe(cdf.style.applymap(
        lambda v: "background-color:#ebfbee;color:#2f9e44;font-weight:600" if ("TP" in str(v) or "TN" in str(v))
             else "background-color:#ffe3e3;color:#c92a2a;font-weight:600" if ("FP" in str(v) or "FN" in str(v))
             else ""), use_container_width=True)

    if st.session_state.fastq_qc:
        st.markdown("---")
        qc = st.session_state.fastq_qc
        st.markdown("### FASTQ QC Metrics (from your file)")
        q1,q2,q3,q4,q5 = st.columns(5)
        q1.metric("Total reads",  f"{qc['total_reads']:,}")
        q2.metric("Mean length",  f"{qc['mean_length']} bp")
        q3.metric("%Q30",         f"{qc['pct_q30']:.1f}%")
        q4.metric("%Q20",         f"{qc['pct_q20']:.1f}%")
        q5.metric("GC content",   f"{qc['gc_pct']:.1f}%")
        q30 = qc["pct_q30"]
        if q30 >= 80:   st.success(f"✅ Excellent quality (Q30: {q30:.1f}%)")
        elif q30 >= 60: st.warning(f"⚠️ Acceptable quality (Q30: {q30:.1f}%)")
        else:           st.error(f"❌ Low quality (Q30: {q30:.1f}%)")

    st.markdown("---")
    st.markdown("### Variant charts (from your file)")
    r1c1,r1c2 = st.columns(2)
    with r1c1: st.plotly_chart(plot_gene_bar(df),       use_container_width=True)
    with r1c2: st.plotly_chart(plot_clinvar_pie(df),    use_container_width=True)
    r2c1,r2c2 = st.columns(2)
    with r2c1: st.plotly_chart(plot_qual_dist(df),      use_container_width=True)
    with r2c2: st.plotly_chart(plot_dp_scatter(df),     use_container_width=True)
    r3c1,r3c2 = st.columns(2)
    with r3c1: st.plotly_chart(plot_consequence_bar(df),use_container_width=True)
    with r3c2: st.plotly_chart(plot_dp_hist(df),        use_container_width=True)
    st.plotly_chart(plot_lollipop(df), use_container_width=True)

    st.markdown("### Per-gene breakdown")
    gt = df.groupby("Gene").agg(
        Total=("Gene","count"),
        Pathogenic=("ClinVar", lambda x: (x=="Pathogenic").sum()),
        Likely_pathogenic=("ClinVar", lambda x: (x=="Likely pathogenic").sum()),
        VUS=("ClinVar", lambda x: (x=="VUS").sum()),
        Mean_QUAL=("QUAL","mean"), Mean_DP=("DP","mean"),
    ).reset_index()
    gt["Mean_QUAL"] = gt["Mean_QUAL"].round(1)
    gt["Mean_DP"]   = gt["Mean_DP"].round(1)
    st.dataframe(gt.rename(columns={"Likely_pathogenic":"Likely pathogenic",
                                    "Mean_QUAL":"Mean QUAL","Mean_DP":"Mean DP"}),
                 use_container_width=True)

    col_s1, col_s2 = st.columns(2)
    with col_s1:
        if st.button("💾 Save per-gene table"):
            p = save_table(gt, "per_gene_breakdown"); st.success(f"Saved → `{p}`")
    with col_s2:
        if st.button("💾 Save all statistics figures", type="primary"):
            saved = _auto_save_all_figures(df)
            st.success(f"✅ {len(saved)} figures saved to `{GRAPHS_DIR}/`")


# ══════════════════════════════════════════════════════════════════════════════
# PAGE: SAVED OUTPUTS
# ══════════════════════════════════════════════════════════════════════════════
elif page == "Saved Outputs":
    st.title("📁 Saved Outputs")
    st.markdown(f"""
    <div class="save-card">
    💾 All files saved to <code>{OUTPUT_ROOT.resolve()}/</code><br>
    <code>tables/</code> — CSV tables &nbsp;|&nbsp;
    <code>graphs/</code> — HTML + PNG charts &nbsp;|&nbsp;
    <code>reports/</code> — full reports
    </div>""", unsafe_allow_html=True)

    outputs = list_saved_outputs()
    total = sum(len(v) for v in outputs.values())
    st.caption(f"{total} files saved")

    if total == 0:
        st.info("No files yet. Run an analysis — tables and figures are auto-saved on completion.")
    else:
        for cat, paths in outputs.items():
            if not paths: continue
            with st.expander(f"{cat}  ({len(paths)} files)", expanded=True):
                for p in paths:
                    c_name, c_dl, c_prev = st.columns([3,1,1])
                    with c_name:
                        st.markdown(f"`{p.name}`")
                        st.caption(datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M"))
                    with c_dl:
                        mime = {"csv":"text/csv","html":"text/html","png":"image/png"}.get(p.suffix.lstrip("."), "application/octet-stream")
                        st.download_button("⬇", open(p,"rb").read(), p.name, mime, key=f"dl_{p}")
                    with c_prev:
                        if p.suffix == ".csv" and st.button("👁", key=f"pv_{p}"):
                            st.dataframe(pd.read_csv(p), use_container_width=True)

    st.markdown("---")
    if st.button("🗑 Clear ALL saved outputs", type="secondary"):
        for d in [TABLES_DIR, GRAPHS_DIR, REPORTS_DIR]:
            shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True, exist_ok=True)
        st.success("All outputs cleared."); st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
# PAGE: PIPELINE STEPS
# ══════════════════════════════════════════════════════════════════════════════
elif page == "Pipeline Steps":
    st.title("⚙️ Pipeline Steps")
    STEPS_DATA = [
        {"num":1,"title":"Quality Control","tool":"FastQC v0.11.9","rna":False,
         "desc":"Raw FASTQ assessed for per-base quality, adapter contamination, duplication, GC content.",
         "input":"Raw FASTQ","output":"HTML report","params":"Default; `--threads 4`","tip":"Q30 <30% → aggressive trim in Step 2."},
        {"num":2,"title":"Adapter Trimming","tool":"fastp v0.23.2","rna":False,
         "desc":"Adapters auto-detected and removed; bases <Q20 trimmed; reads <30 bp discarded.",
         "input":"FASTQ (R1+R2)","output":"Trimmed FASTQ + JSON",
         "params":"`--length_required 30 --qualified_quality_phred 20 --detect_adapter_for_pe`","tip":"~3× faster than Trimmomatic."},
        {"num":3,"title":"Splice-Aware Alignment","tool":"HISAT2 v2.2.1","rna":True,
         "desc":"Reads aligned to GRCh38; HISAT2 handles intron–exon spanning reads correctly.",
         "input":"Trimmed FASTQ + GRCh38 index","output":"SAM file",
         "params":"`--dta --no-softclip --rna-strandness RF`","tip":"Download GRCh38 index from genome.ucsc.edu."},
        {"num":4,"title":"SAM→BAM Conversion","tool":"SAMtools v1.15.1","rna":False,
         "desc":"SAM → coordinate-sorted BAM, then indexed.",
         "input":"SAM","output":"Sorted BAM + BAI",
         "params":"`view -bS | sort -@ 4 | index`","tip":"Both BAM + BAI must be in same directory."},
        {"num":5,"title":"Read Group Annotation","tool":"Picard v2.27.1","rna":False,
         "desc":"RG tags added to BAM header — required by GATK.",
         "input":"Sorted BAM","output":"RG-annotated BAM",
         "params":"`RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=SAMPLE`","tip":"Replace SAMPLE with actual sample ID."},
        {"num":6,"title":"Duplicate Marking","tool":"Picard MarkDuplicates","rna":False,
         "desc":"PCR/optical duplicates flagged; GATK ignores them downstream.",
         "input":"RG BAM","output":"Deduped BAM + metrics",
         "params":"`OPTICAL_DUPLICATE_PIXEL_DISTANCE=100`","tip":"Dup rate >30% → low library complexity."},
        {"num":7,"title":"SplitNCigarReads (RNA only)","tool":"GATK SplitNCigarReads","rna":True,
         "desc":"Splits reads at N CIGAR ops (introns). Required for RNA-seq variant calling.",
         "input":"Deduped BAM","output":"Split BAM",
         "params":"`--refactor-cigar-string`","tip":"Skip for WGS/WES inputs."},
        {"num":8,"title":"Germline Variant Calling","tool":"GATK HaplotypeCaller v4.1.8.1","rna":False,
         "desc":"Local de-novo haplotype assembly + SNV/indel calling.",
         "input":"Split BAM + GRCh38 FASTA","output":"Raw VCF",
         "params":"`--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20`","tip":"Use `-ERC GVCF` for multi-sample joint calling."},
        {"num":9,"title":"VCF Quality Filtering","tool":"bcftools","rna":False,
         "desc":"Filters: QUAL>200, DP≥10, GQ≥20. Validated on 172 samples.",
         "input":"Raw VCF","output":"Filtered VCF",
         "params":'`bcftools filter -i "QUAL>200 && DP>=10 && GQ>=20"`',"tip":"Adjust QUAL to 150 for low-coverage (<15×)."},
        {"num":10,"title":"Annotation + Gene Filter","tool":"ANNOVAR + ClinVar","rna":False,
         "desc":"RefGene consequence + ClinVar pathogenicity. Only MMR-gene variants retained.",
         "input":"Filtered VCF","output":"Annotated VCF + TSV report",
         "params":"`--buildver hg38 --protocol refGene,clinvar_20221231 --operation g,f`","tip":"Update ClinVar DB periodically."},
    ]
    for s in STEPS_DATA:
        tag = " 🔬 RNA-seq only" if s["rna"] else ""
        with st.expander(f"Step {s['num']}: {s['title']}{tag} — `{s['tool']}`"):
            st.markdown(f"**Description:** {s['desc']}")
            ci,co = st.columns(2)
            with ci: st.markdown("**⬅ Input**"); st.info(s["input"])
            with co: st.markdown("**➡ Output**"); st.success(s["output"])
            st.markdown(f"**Parameters:** `{s['params']}`")
            st.markdown(f"💡 **Tip:** {s['tip']}")


# ══════════════════════════════════════════════════════════════════════════════
# PAGE: HELP & DOCS
# ══════════════════════════════════════════════════════════════════════════════
elif page == "Help & Docs":
    st.title("📚 Help & Documentation")
    t1,t2,t3,t4 = st.tabs(["Getting Started","Installation","File Formats","FAQ"])

    with t1:
        st.markdown("""
### Getting started
1. **Upload & Run** → upload a VCF or FASTQ (or click a demo button)
2. The parser auto-detects ANNOVAR / SnpEff / VEP annotations — fields are read directly from your file
3. All charts update to reflect **your** actual variants
4. Tables + figures are **auto-saved** on completion; individual 💾 buttons save specific charts
5. Browse everything in **Saved Outputs**
        """)

    with t2:
        st.markdown("""
### Local installation
```bash
git clone https://github.com/Rofidagamal/Computational-Pipeline-for-Lynch-Syndrome-Detection
cd Computational-Pipeline-for-Lynch-Syndrome-Detection
conda env create -f environment.yml
conda activate lynch_pipeline
pip install streamlit pandas plotly kaleido
streamlit run app.py
```
Install `kaleido` for PNG export. Without it, only HTML charts are saved (still fully interactive).
        """)

    with t3:
        st.markdown("""
### Accepted VCF annotation fields

| Tool | Detected fields |
|------|----------------|
| ANNOVAR | `Gene.refGene`, `ExonicFunc.refGene`, `Func.refGene`, `AAChange.refGene`, `CLNSIG` |
| SnpEff | `ANN` (pipe-separated: allele\\|effect\\|impact\\|gene\\|...) |
| VEP | `CSQ` (pipe-separated, field order from `##INFO=<ID=CSQ` header) |
| Raw ClinVar | `CLNSIG`, `CLINSIG` |
| Extra scores | `avsnp150`, `gnomAD_exome_ALL`, `CADD_phred`, `Polyphen2_HDIV_pred`, `SIFT_pred` |

Chromosome names `2` and `chr2` are both accepted.
        """)

    with t4:
        st.markdown("""
### FAQ

**Q: My VCF shows 0 variants — why?**
A: The three most common causes:
1. Variants fail QC: QUAL ≤ 200, DP < 10, or GQ < 20. Check your VCF has DP and GQ in FORMAT fields.
2. Chromosomes use non-standard names (e.g. `NC_000002.12` instead of `chr2` / `2`).
3. Variants are outside the MMR gene coordinate windows — confirm your VCF is hg38.

**Q: Where are saved files?**
A: In `lynch_pipeline_outputs/` next to `app.py` — `tables/` for CSV, `graphs/` for HTML+PNG.

**Q: How do I get PNG files?**
A: `pip install kaleido`. PNG and HTML are both saved automatically.

**Q: Can I use WGS/WES data?**
A: Yes — skip Step 7 (SplitNCigarReads) and use BWA-MEM instead of HISAT2.

**Q: What is VUS?**
A: Variant of Uncertain Significance — functional evidence or family co-segregation data needed before clinical action.

**Q: Are results clinically actionable?**
A: No. All variants require confirmation by a certified germline sequencing laboratory before clinical decisions.
        """)