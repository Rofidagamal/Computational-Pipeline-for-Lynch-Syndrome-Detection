#!/usr/bin/env bash

# =============================================================================
# RNA-Seq Variant Calling Pipeline — Lynch Syndrome
# Style: modular helper functions, USER CONFIGURATION block at top
# Tools: fastp → HISAT2 → Picard → GATK HaplotypeCaller → ANNOVAR
# =============================================================================

# ---------------------- USER CONFIGURATION ----------------------
PROJECT_HOME="${1:-$HOME/Lynch_RNA}"
FASTQ_DIR="${PROJECT_HOME}/Data_Folder"
REFERENCE_FILE="${PROJECT_HOME}/Alignment/Homo_sapiens_assembly38.fasta"
ANNOVAR_DIR="${PROJECT_HOME}/annovar.latest/annovar"
GATK_PATH="${PROJECT_HOME}/gatk-4.1.8.1/gatk"
PICARD_JAR="/usr/local/bin/picard.jar"
THREADS=20

RNA_FASTP="${PROJECT_HOME}/RNA_Fastp"
RNA_ALIGNED="${PROJECT_HOME}/RNA_Aligned"
RNA_PROCESSED="${PROJECT_HOME}/RNA_Processed"
RNA_VCF="${PROJECT_HOME}/RNA_VCF"
RNA_ANNOTATED="${PROJECT_HOME}/RNA_Annotated"
RNA_REPORT="${PROJECT_HOME}/RNA_Lynch_Report"

# Lynch syndrome gene list and hg38 coordinates
declare -A GENE_RANGES=(
    ["MLH1"]="chr3:3702015-3733317"
    ["MSH2"]="chr2:47650919-47674520"
    ["MSH6"]="chr2:48032114-48129989"
    ["PMS2"]="chr7:5985774-6004636"
    ["EPCAM"]="chr2:47472139-47490377"
)

ANNOTATION_LIST="MLH1|MSH2|MSH6|PMS2|EPCAM"
# ----------------------------------------------------------------


# =============================================================================
# HELPER: Install a tool if not already present (tries conda, then apt)
# Usage: install_if_missing <command_name> <apt_package_name>
# =============================================================================
install_if_missing() {
    local cmd="$1"
    local pkg="$2"
    if ! command -v "$cmd" &>/dev/null; then
        echo "Installing $cmd..."
        if command -v conda &>/dev/null; then
            conda install -y -c bioconda "$cmd" && return
        fi
        sudo apt-get update -qq
        sudo apt-get install -y "$pkg"
    else
        echo "$cmd is already installed."
    fi
}

# =============================================================================
# HELPER: Install Picard if not present
# =============================================================================
install_picard_if_missing() {
    if [ -f "$PICARD_JAR" ]; then
        echo "Picard is already installed."
    else
        echo "Installing Picard Tools..."
        wget -q https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar
        sudo mv picard.jar "$PICARD_JAR"
        sudo ln -sf "$PICARD_JAR" /usr/local/bin/picard
        echo "Picard installed."
    fi
}

# =============================================================================
# HELPER: Install GATK if not present
# =============================================================================
install_gatk_if_missing() {
    if [ -f "$GATK_PATH" ]; then
        echo "GATK4 is already installed."
    elif command -v gatk &>/dev/null; then
        echo "GATK4 is already installed (system)."
        GATK_PATH=$(command -v gatk)
    else
        echo "Installing GATK4..."
        wget -q https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip
        unzip -q gatk-4.1.8.1.zip -d "$PROJECT_HOME"
        rm gatk-4.1.8.1.zip
        echo "GATK4 installed at: $GATK_PATH"
    fi

    # Ensure Python symlink exists (required by GATK)
    if ! dpkg -s python3 &>/dev/null; then
        sudo apt-get install -y python3
    fi
    if [ ! -e /usr/bin/python ]; then
        sudo ln -s /usr/bin/python3 /usr/bin/python
    fi
}

# =============================================================================
# STEP 1: FASTP — quality control and adapter trimming (paired-end)
# =============================================================================
run_fastp() {
    local sample_id="$1"
    local r1="$2"
    local r2="$3"
    local trimmed_r1="${RNA_FASTP}/${sample_id}_R1.fil.fastq.gz"
    local trimmed_r2="${RNA_FASTP}/${sample_id}_R2.fil.fastq.gz"

    if [ -f "$trimmed_r1" ] && [ -f "$trimmed_r2" ]; then
        echo "Trimmed FASTQs already exist for $sample_id. Skipping fastp."
        return 0
    fi

    echo "Running fastp on: $sample_id"
    fastp \
        -i "$r1" -I "$r2" \
        -o "$trimmed_r1" -O "$trimmed_r2" \
        --length_required 30 \
        --json "${RNA_FASTP}/${sample_id}.json" \
        --html "${RNA_FASTP}/${sample_id}.html" \
        || { echo "fastp failed for $sample_id"; return 1; }

    echo "Trimming complete: $sample_id"
}

# =============================================================================
# STEP 2: ALIGNMENT — HISAT2 (splice-aware, required for RNA-Seq)
#   --rna-strandness RF : Illumina paired-end stranded RNA-Seq
#   --dta               : downstream transcript assembly compatible output
# Followed by SAM → BAM conversion via Picard SamFormatConverter
# =============================================================================
run_alignment() {
    local sample_id="$1"
    local sam_file="${RNA_ALIGNED}/${sample_id}.sam"
    local bam_file="${RNA_ALIGNED}/${sample_id}.bam"
    local ref_dir
    ref_dir="$(dirname "$REFERENCE_FILE")"

    # Build HISAT2 index if missing
    if [ -f "${ref_dir}/genome.1.ht2" ]; then
        echo "HISAT2 index found."
    else
        echo "Building HISAT2 index..."
        hisat2-build "$REFERENCE_FILE" "${ref_dir}/genome"
    fi

    # Align with HISAT2
    if [ -f "$sam_file" ]; then
        echo "SAM already exists for $sample_id. Skipping alignment."
    else
        echo "Aligning RNA-Seq reads with HISAT2: $sample_id"
        hisat2 \
            -x "${ref_dir}/genome" \
            -1 "${RNA_FASTP}/${sample_id}_R1.fil.fastq.gz" \
            -2 "${RNA_FASTP}/${sample_id}_R2.fil.fastq.gz" \
            -S "$sam_file" \
            --rna-strandness RF \
            --dta \
            2>> "${RNA_ALIGNED}/${sample_id}_hisat2.log" \
            || { echo "HISAT2 failed for $sample_id"; return 1; }
    fi

    # Convert SAM → BAM
    if [ -f "$bam_file" ]; then
        echo "BAM already exists for $sample_id. Skipping conversion."
    else
        echo "Converting SAM to BAM: $sample_id"
        java -Xmx8g -jar "$PICARD_JAR" SamFormatConverter \
            -I "$sam_file" \
            -O "$bam_file" \
            || { echo "SAM→BAM conversion failed for $sample_id"; return 1; }
    fi
}

# =============================================================================
# STEP 3: ADD READ GROUPS
# =============================================================================
run_add_read_groups() {
    local sample_id="$1"
    local rg_bam="${RNA_PROCESSED}/${sample_id}.rg.bam"

    if [ -f "$rg_bam" ]; then
        echo "Read groups BAM already exists for $sample_id. Skipping."
        return 0
    fi

    echo "Adding read groups: $sample_id"
    java -jar "$PICARD_JAR" AddOrReplaceReadGroups \
        I="${RNA_ALIGNED}/${sample_id}.bam" \
        O="$rg_bam" \
        RGID="ID_${sample_id}" \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM="$sample_id" \
        || { echo "AddReadGroups failed for $sample_id"; return 1; }
}

# =============================================================================
# STEP 4: SORT BAM BY COORDINATE
# =============================================================================
run_sort() {
    local sample_id="$1"
    local sorted_bam="${RNA_PROCESSED}/${sample_id}.sorted.bam"

    if [ -f "$sorted_bam" ]; then
        echo "Sorted BAM already exists for $sample_id. Skipping."
        return 0
    fi

    echo "Sorting BAM: $sample_id"
    java -Xmx8g -jar "$PICARD_JAR" SortSam \
        I="${RNA_PROCESSED}/${sample_id}.rg.bam" \
        O="$sorted_bam" \
        SORT_ORDER=coordinate \
        || { echo "SortSam failed for $sample_id"; return 1; }
}

# =============================================================================
# STEP 5: MARK PCR DUPLICATES
# =============================================================================
run_mark_duplicates() {
    local sample_id="$1"
    local marked_bam="${RNA_PROCESSED}/${sample_id}.sorted.marked.bam"

    if [ -f "$marked_bam" ]; then
        echo "Marked duplicates BAM already exists for $sample_id. Skipping."
        return 0
    fi

    echo "Marking duplicates: $sample_id"
    java -Xmx8g -jar "$PICARD_JAR" MarkDuplicates \
        I="${RNA_PROCESSED}/${sample_id}.sorted.bam" \
        O="$marked_bam" \
        M="${RNA_PROCESSED}/${sample_id}_dup_metrics.txt" \
        || { echo "MarkDuplicates failed for $sample_id"; return 1; }
}

# =============================================================================
# STEP 6: INDEX BAM
# =============================================================================
run_index() {
    local sample_id="$1"
    local marked_bam="${RNA_PROCESSED}/${sample_id}.sorted.marked.bam"

    if [ -f "${marked_bam}.bai" ]; then
        echo "BAM index already exists for $sample_id. Skipping."
        return 0
    fi

    echo "Indexing BAM: $sample_id"
    samtools index "$marked_bam" \
        || { echo "Indexing failed for $sample_id"; return 1; }
}

# =============================================================================
# STEP 7: SPLIT N CIGAR READS (RNA-Seq REQUIRED before HaplotypeCaller)
# RNA reads spanning splice junctions carry 'N' in the CIGAR string.
# GATK misinterprets these N-gaps as mismatches, producing false variants.
# SplitNCigarReads breaks spanning reads into continuous sub-reads at each N.
# =============================================================================
run_split_n_cigar() {
    local sample_id="$1"
    local split_bam="${RNA_PROCESSED}/${sample_id}.sorted.marked.split.bam"

    if [ -f "$split_bam" ]; then
        echo "SplitNCigarReads BAM already exists for $sample_id. Skipping."
        return 0
    fi

    echo "Running SplitNCigarReads: $sample_id"
    "$GATK_PATH" SplitNCigarReads \
        -R "$REFERENCE_FILE" \
        -I "${RNA_PROCESSED}/${sample_id}.sorted.marked.bam" \
        -O "$split_bam" \
        || { echo "SplitNCigarReads failed for $sample_id"; return 1; }
}

# =============================================================================
# STEP 8: HAPLOTYPECALLER (RNA mode)
# Input: split BAM from SplitNCigarReads
# RNA-specific flags:
#   --dont-use-soft-clipped-bases              : ignore soft-clipped bases common in RNA
#   --standard-min-confidence-threshold-for-calling 20 : lower threshold for RNA data
# =============================================================================
run_haplotypecaller() {
    local sample_id="$1"
    local vcf_file="${RNA_VCF}/${sample_id}.vcf"

    if [ -f "${vcf_file}.gz" ]; then
        echo "VCF already exists for $sample_id. Skipping HaplotypeCaller."
        return 0
    fi

    echo "Running HaplotypeCaller (RNA): $sample_id"
    "$GATK_PATH" --java-options "-Xmx16G" HaplotypeCaller \
        -R "$REFERENCE_FILE" \
        -I "${RNA_PROCESSED}/${sample_id}.sorted.marked.split.bam" \
        -O "$vcf_file" \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20 \
        || { echo "HaplotypeCaller failed for $sample_id"; return 1; }

    bgzip -f "$vcf_file"
    tabix -p vcf "${vcf_file}.gz"
}

# =============================================================================
# STEP 9: VCF FILTERING
# Thresholds for HaplotypeCaller RNA output:
#   QUAL > 200, DP >= 10, GQ >= 20
# =============================================================================
run_filter_vcf() {
    local sample_id="$1"
    local input_vcf="${RNA_VCF}/${sample_id}.vcf.gz"
    local filtered_vcf="${RNA_VCF}/${sample_id}.filtered.vcf.gz"

    if [ -f "$filtered_vcf" ]; then
        echo "Filtered VCF already exists for $sample_id. Skipping."
        return 0
    fi

    echo "Filtering VCF: $sample_id"
    bcftools filter \
        -O z \
        -o "$filtered_vcf" \
        -i 'QUAL > 200 && INFO/DP >= 10 && GQ >= 20' \
        "$input_vcf" \
        || { echo "bcftools filter failed for $sample_id"; return 1; }

    bcftools index "$filtered_vcf"
}

# =============================================================================
# STEP 10: ANNOVAR ANNOTATION + LYNCH GENE DETECTION
# =============================================================================
annotate_with_annovar_and_detect_lynch() {
    local vcf_file="$1"
    local base_name
    base_name=$(basename "${vcf_file%.vcf.gz}")
    local avinput_file="${RNA_ANNOTATED}/${base_name}.avinput"
    local geneanno_output="${RNA_ANNOTATED}/${base_name}.geneanno.txt"
    local lynch_output="${RNA_ANNOTATED}/${base_name}.lynch.txt"

    echo "Annotating with ANNOVAR: $base_name"
    perl "$ANNOVAR_DIR/convert2annovar.pl" \
        -format vcf4old "$vcf_file" \
        -outfile "$avinput_file" || return 1

    perl "$ANNOVAR_DIR/annotate_variation.pl" \
        -buildver hg38 \
        -geneanno \
        -dbtype refGene \
        "$avinput_file" \
        "$ANNOVAR_DIR/humandb/" > "$geneanno_output" || return 1

    grep -Ei "$ANNOTATION_LIST" "$geneanno_output" > "$lynch_output"
    echo "ANNOVAR annotation saved: $lynch_output"
}

# =============================================================================
# STEP 11: GENERATE LYNCH SYNDROME POSITIONAL REPORT
# Scans each filtered VCF against known Lynch gene hg38 coordinates.
# Reports which genes have overlapping variants per sample.
# =============================================================================
generate_lynch_position_report() {
    mkdir -p "$RNA_REPORT"

    for vcf_file in "${RNA_VCF}"/*.filtered.vcf.gz; do
        local base
        base=$(basename "$vcf_file" .filtered.vcf.gz)
        local report_file="${RNA_REPORT}/${base}_lynch_report.txt"
        local found_any=0

        echo "Checking Lynch syndrome regions in $base..." > "$report_file"

        for gene in "${!GENE_RANGES[@]}"; do
            local region="${GENE_RANGES[$gene]}"
            local chr start end
            chr=$(cut -d':' -f1 <<< "$region")
            start=$(cut -d':' -f2 <<< "$region" | cut -d'-' -f1)
            end=$(cut -d':' -f2 <<< "$region" | cut -d'-' -f2)

            local match
            match=$(zcat "$vcf_file" | awk \
                -v chr="$chr" -v start="$start" -v end="$end" -F'\t' \
                '$0 !~ /^#/ && $1 == chr && $2 >= start && $2 <= end { print }')

            if [[ -n "$match" ]]; then
                echo "  [+] Variant found in $gene ($region)" >> "$report_file"
                found_any=1
            fi
        done

        if [[ "$found_any" -eq 0 ]]; then
            echo "  [-] No Lynch syndrome variants found." >> "$report_file"
        fi

        echo "Report saved: $report_file"
    done
}

# =============================================================================
# MAIN
# =============================================================================
main() {
    echo "========================================================"
    echo " RNA-Seq Lynch Syndrome Pipeline"
    echo " Project: $PROJECT_HOME"
    echo "========================================================"

    # Install all required tools
    install_if_missing fastp fastp
    install_if_missing hisat2 hisat2
    install_if_missing samtools samtools
    install_if_missing bcftools bcftools
    install_if_missing tabix tabix
    install_if_missing bgzip tabix
    install_if_missing java openjdk-8-jdk
    install_picard_if_missing
    install_gatk_if_missing

    # Create output directories
    mkdir -p "$RNA_FASTP" "$RNA_ALIGNED" "$RNA_PROCESSED" "$RNA_VCF" "$RNA_ANNOTATED" "$RNA_REPORT"

    # Process each paired-end sample found in FASTQ_DIR
    for r1 in "$FASTQ_DIR"/*_1.fastq.gz; do
        [ -f "$r1" ] || { echo "No FASTQ files found in $FASTQ_DIR"; exit 1; }

        local sample_id
        sample_id=$(basename "$r1" _1.fastq.gz)
        local r2="${FASTQ_DIR}/${sample_id}_2.fastq.gz"

        if [ ! -f "$r2" ]; then
            echo "WARNING: Paired file not found for $sample_id — skipping."
            continue
        fi

        echo ""
        echo "--- Processing sample: $sample_id ---"

        run_fastp           "$sample_id" "$r1" "$r2"  || continue
        run_alignment       "$sample_id"              || continue
        run_add_read_groups "$sample_id"              || continue
        run_sort            "$sample_id"              || continue
        run_mark_duplicates "$sample_id"              || continue
        run_index           "$sample_id"              || continue
        run_split_n_cigar   "$sample_id"              || continue
        run_haplotypecaller "$sample_id"              || continue
        run_filter_vcf      "$sample_id"              || continue
        annotate_with_annovar_and_detect_lynch \
            "${RNA_VCF}/${sample_id}.filtered.vcf.gz" || continue

        echo "Finished: $sample_id"
    done

    echo ""
    echo "Generating Lynch syndrome positional reports..."
    generate_lynch_position_report

    echo ""
    echo "Pipeline complete. Results in:"
    echo "  Annotated variants : $RNA_ANNOTATED"
    echo "  Lynch reports      : $RNA_REPORT"
}

main

: <<'COMMENT'
Run with: bash RNA_Lynch_Pipeline.sh [/optional/project/path]

Default project home: $HOME/Lynch_RNA

Expected FASTQ layout:
  $PROJECT_HOME/Data_Folder/SAMPLE_1.fastq.gz
  $PROJECT_HOME/Data_Folder/SAMPLE_2.fastq.gz

Required external tools (auto-installed if missing):
  fastp, hisat2, samtools, bcftools, tabix/bgzip,
  openjdk-8, Picard 2.23.3, GATK 4.1.8.1, ANNOVAR

ANNOVAR must be pre-downloaded and placed at:
  $PROJECT_HOME/annovar.latest/annovar/

Reference genome (hg38) must be placed at:
  $PROJECT_HOME/Alignment/Homo_sapiens_assembly38.fasta

Lynch syndrome genes covered: MLH1, MSH2, MSH6, PMS2, EPCAM
COMMENT
