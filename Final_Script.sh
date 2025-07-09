#!/usr/bin/env bash

# ---------------------- USER CONFIGURATION ----------------------
PROJECT_HOME="${1:-$HOME/NILES3}"
FASTQ_DIR="${PROJECT_HOME}/Data_Folder"
REFERENCE_FILE="${PROJECT_HOME}/Alignment/Homo_sapiens_assembly38.fasta"
ANNOVAR_DIR="${PROJECT_HOME}/annovar.latest/annovar"
GATK_PATH="${PROJECT_HOME}/gatk-4.1.8.1/gatk"
THREADS=20
DNA_FASTP="${PROJECT_HOME}/DNA_Fastp"
DNA_ALIGNED="${PROJECT_HOME}/DNA_Aligned"
DNA_PROCESSED="${PROJECT_HOME}/DNA_Processed"
DNA_VCF="${PROJECT_HOME}/DNA_VCF"
DNA_ANNOTATED="${PROJECT_HOME}/DNA_Annotated"
DNA_REPORT="${PROJECT_HOME}/DNA_Annotated_Report"

# Lynch syndrome gene list and hg38 coordinates
declare -A GENE_RANGES=(
    ["MLH1"]="chr3:3702015-3733317"
    ["MSH2"]="chr2:47650919-47674520"
    ["MSH6"]="chr2:48032114-48129989"
    ["PMS2"]="chr7:5985774-6004636"
    ["EPCAM"]="chr2:47472139-47490377"
)

ANNOTATION_LIST="MLH1|MSH2|MSH6|PMS2|EPCAM"

install_if_missing() {
    if ! command -v "$1" &> /dev/null; then
        echo "Installing $1..."
        if command -v conda &> /dev/null; then
            conda install -y -c bioconda "$1" && return
        fi
        sudo apt-get update -qq
        sudo apt-get install -y "$2"
    else
        echo "$1 is already installed."
    fi
}

annotate_with_annovar_and_detect_lynch() {
    local vcf_file="$1"
    local base_name=$(basename "${vcf_file%.vcf.gz}")
    local avinput_file="${DNA_ANNOTATED}/${base_name}.avinput"
    local geneanno_output="${DNA_ANNOTATED}/${base_name}.geneanno.txt"
    local lynch_output="${DNA_ANNOTATED}/${base_name}.lynch.txt"

    echo "Annotating with ANNOVAR..."
    perl "$ANNOVAR_DIR/convert2annovar.pl" -format vcf4old "$vcf_file" -outfile "$avinput_file" || return 1
    perl "$ANNOVAR_DIR/annotate_variation.pl" -buildver hg38 -geneanno -dbtype refGene "$avinput_file" "$ANNOVAR_DIR/humandb/" > "$geneanno_output" || return 1
    grep -Ei "$ANNOTATION_LIST" "$geneanno_output" > "$lynch_output"
    echo "ANNOVAR results saved: $lynch_output"
}

generate_lynch_position_report() {
    mkdir -p "$DNA_REPORT"
    for file in "$DNA_VCF"/*.vcf.gz; do
        base=$(basename "$file" .vcf.gz)
        report_file="${DNA_REPORT}/${base}_lynch_report.txt"
        found_any=0
        echo "Checking Lynch syndrome regions in $base..." > "$report_file"

        for gene in "${!GENE_RANGES[@]}"; do
            region="${GENE_RANGES[$gene]}"
            chr=$(cut -d':' -f1 <<< "$region")
            start=$(cut -d':' -f2 <<< "$region" | cut -d'-' -f1)
            end=$(cut -d':' -f2 <<< "$region" | cut -d'-' -f2)

            match=$(zcat "$file" | awk -v chr="$chr" -v start="$start" -v end="$end" -F'\t' ' $0 !~ /^#/ && $1 == chr && $2 >= start && $2 <= end { print }')

            if [[ -n "$match" ]]; then
                echo "Found variant in $gene ($region)" >> "$report_file"
                found_any=1
            fi
        done

        if [[ "$found_any" -eq 0 ]]; then
            echo "No Lynch syndrome variants found in $base" >> "$report_file"
        fi

        echo "Report saved: $report_file"
    done
}

main() {
    echo "Starting pipeline in: $PROJECT_HOME"

    install_if_missing fastp fastp
    install_if_missing bwa bwa
    install_if_missing samtools samtools
    install_if_missing bcftools bcftools
    install_if_missing tabix tabix
    install_if_missing bgzip tabix
    install_if_missing picard picard

    mkdir -p "$DNA_FASTP" "$DNA_ALIGNED" "$DNA_PROCESSED" "$DNA_VCF" "$DNA_ANNOTATED"

    for r1 in "$FASTQ_DIR"/*_1.fastq.gz; do
        sample_id=$(basename "$r1" _1.fastq.gz)
        r2="${FASTQ_DIR}/${sample_id}_2.fastq.gz"

        if [[ -f "$r2" ]]; then
            echo "Processing Paired-End sample: $sample_id"
            trimmed_r1="${DNA_FASTP}/${sample_id}_R1_trimmed.fastq.gz"
            trimmed_r2="${DNA_FASTP}/${sample_id}_R2_trimmed.fastq.gz"
            sam_file="${DNA_ALIGNED}/${sample_id}.sam"
            bam_file="${DNA_PROCESSED}/${sample_id}.bam"
            sorted_bam="${DNA_PROCESSED}/${sample_id}.sorted.bam"
            rg_bam="${DNA_PROCESSED}/${sample_id}.rg.bam"
            vcf_file="${DNA_VCF}/${sample_id}.vcf"

            fastp -i "$r1" -I "$r2" -o "$trimmed_r1" -O "$trimmed_r2" --length_required 30 \
                --json "${DNA_FASTP}/${sample_id}.json" --html "${DNA_FASTP}/${sample_id}.html" || continue

            bwa mem -t "$THREADS" "$REFERENCE_FILE" "$trimmed_r1" "$trimmed_r2" > "$sam_file" || continue
        else
            echo "Processing Single-End sample: $sample_id"
            trimmed="${DNA_FASTP}/${sample_id}_trimmed.fastq.gz"
            sam_file="${DNA_ALIGNED}/${sample_id}.sam"
            bam_file="${DNA_PROCESSED}/${sample_id}.bam"
            sorted_bam="${DNA_PROCESSED}/${sample_id}.sorted.bam"
            rg_bam="${DNA_PROCESSED}/${sample_id}.rg.bam"
            vcf_file="${DNA_VCF}/${sample_id}.vcf"

            fastp -i "$r1" -o "$trimmed" --length_required 30 \
                --json "${trimmed}.json" --html "${trimmed}.html" || continue

            bwa mem -t "$THREADS" "$REFERENCE_FILE" "$trimmed" > "$sam_file" || continue
        fi

        samtools view -Sb "$sam_file" > "$bam_file" || continue
        samtools sort "$bam_file" -o "$sorted_bam" || continue
        samtools index "$sorted_bam" || continue

        picard AddOrReplaceReadGroups \
            I="$sorted_bam" O="$rg_bam" \
            RGID="$sample_id" RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$sample_id" || continue
        samtools index "$rg_bam" || continue

        "$GATK_PATH" --java-options "-Xmx16G" HaplotypeCaller \
            -R "$REFERENCE_FILE" -I "$rg_bam" -O "$vcf_file" || continue

        bgzip -f "$vcf_file"
        tabix -p vcf "${vcf_file}.gz"

        annotate_with_annovar_and_detect_lynch "${vcf_file}.gz"
        echo "Finished processing: $sample_id"
    done

    echo "Generating Lynch gene position reports..."
    generate_lynch_position_report
    echo "All done. Check results in:"
    echo "- $DNA_ANNOTATED"
    echo "- $DNA_REPORT"
}

main

