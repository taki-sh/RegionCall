#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N jobname
#$ -l d_rt=300:00:00
#$ -l s_rt=300:00:00
#$ -l s_vmem=12G
#$ -l mem_req=12G
#$ -pe def_slot 32

set -e
set -u
set -o pipefail

export MALLOC_ARENA_MAX=2

# sample list
sample_list=(
)


# info
target_name="" # target name
target_region="" # target region +-10 base
REF= # reference
BWA_REF= # index for bwa
threads="" # number of threads
fastp_threads="" # number of threads for fastp (up to 16)
add_threads="${threads}"-1
WORK_PATH=~/main/results # working directory


# Singularity
shopt -s expand_aliases
alias prefetch="singularity exec /usr/local/biotools/s/sra-tools:2.11.0--pl5321ha49a11a_3 prefetch"
alias fasterq-dump="singularity exec /usr/local/biotools/s/sra-tools:2.11.0--pl5321ha49a11a_3 fasterq-dump"
alias pigz="singularity exec /usr/local/biotools/p/pigz:2.3.4 pigz"
alias fastqc="singularity exec /usr/local/biotools/f/fastqc:0.11.9--hdfd78af_1 fastqc"
alias fastp="singularity exec /usr/local/biotools/f/fastp:0.23.2--hb7a2d85_2 fastp"
alias bwa-mem2="singularity exec /usr/local/biotools/b/bwa-mem2:2.2--he513fc3_0 bwa-mem2"
alias samtools="singularity exec /usr/local/biotools/s/samtools:1.15--h3843a85_0 samtools"
alias gatk="singularity exec /usr/local/biotools/g/gatk4:4.2.6.1--py36hdfd78af_1 gatk"
alias bgzip="singularity exec /usr/local/biotools/s/samtools:1.15--h3843a85_0 bgzip"
alias bcftools="singularity exec /usr/local/biotools/b/bcftools:1.15--haf5b3da_0 bcftools"


for fname in "${sample_list[@]}"
do
    # time
    start_time=$(date "+%Y-%m-%d %H:%M:%S")
    date=$(date "+%Y%m%d%H%M")
    
    # preparation
    mkdir -p "${WORK_PATH}"/"${fname}"/{"${target_name}",bam_ref,fastq,cleaned_fastq}
    mkdir "${WORK_PATH}"/"${fname}"/"${target_name}"/log
    PROJECT_PATH="${WORK_PATH}"/"${fname}"/"${target_name}"
    
    
    # log
    LOG_OUT=${PROJECT_PATH}/log/${date}_stdout.log
    LOG_ERR=${PROJECT_PATH}/log/${date}_stderr.log
    
    exec 1> >(tee -a "${LOG_OUT}")
    exec 2>>"${LOG_ERR}"
    
    echo "${fname}"
    
    
    # get SRA
    cd "${WORK_PATH}"/"${fname}"/fastq || exit
    
    echo "prefetch"
    prefetch \
    -O ./ \
    --max-size u \
    "${fname}"
    
    echo "fasterq-dump"
    fasterq-dump \
    ./"${fname}" \
    -e "${threads}" \
    --split-files
    
    echo "pigz"
    pigz \
    -p "${threads}" \
    *.fastq
    
    
    # fastqc
    echo "fastqc_1"
    fastqc -t "${threads}" ./*fastq.gz
    
    
    #fastp
    echo "fastp"
    fastp -w "${fastp_threads}" \
    -i "${fname}"_1.fastq.gz \
    -I "${fname}"_2.fastq.gz \
    -o "${WORK_PATH}"/"${fname}"/cleaned_fastq/"${fname}"_cleaned_1_"${date}".fastq.gz \
    -O "${WORK_PATH}"/"${fname}"/cleaned_fastq/"${fname}"_cleaned_2_"${date}".fastq.gz
    
    rm ./*.fastq.gz
    rm -rf "${fname}"
    
    
    # check trimming
    cd "${WORK_PATH}"/"${fname}"/cleaned_fastq || exit
    
    echo "fastqc_2"
    fastqc -t "${threads}" ./*_cleaned_*.gz
    
    
    # mapping
    echo "bwa-mem2"
    bamRG="@RG\tID:${fname}\tPL:ILLUMINA\tSM:"${fname}
    bwa-mem2 mem \
    -t "${threads}" \
    -R "${bamRG}" "${BWA_REF}" \
    "${fname}"_cleaned_1_"${date}".fastq.gz \
    "${fname}"_cleaned_2_"${date}".fastq.gz | \
    samtools sort \
    -@ "${add_threads}" \
    -O BAM \
    -o "${WORK_PATH}"/"${fname}"/bam_ref/"${fname}"_"${date}".bam
    
    rm ./*.fastq.gz
    
    (
        echo "samtools flagstat"
        samtools flagstat \
        "${WORK_PATH}"/"${fname}"/bam_ref/"${fname}"_"${date}".bam
    )&
    flagstat=$!
    
    
    # view target
    cd "${PROJECT_PATH}" || exit
    
    mkdir {bam,bqsr,fasta,vcf}
    bam=("${WORK_PATH}"/"${fname}"/bam_ref/*[0-9]_[0-9]*[0-9].bam)
    echo "${bam}"
    
    echo "samtools index"
    samtools index -@ "${threads}" "${bam}"
    
    cd "${PROJECT_PATH}"/bam || exit
    
    echo "samtools view"
    samtools view \
    -o "${fname}"_"${target_name}"_"${date}".bam \
    -h \
    -b \
    -@ "${add_threads}" \
    "${bam}" \
    "${target_region}"
    
    
    # duplication markup
    echo "gatk MarkDuplicates"
    gatk MarkDuplicates \
    -I "${fname}"_"${target_name}"_"${date}".bam \
    -O "${fname}"_"${target_name}"_"${date}"_markdup.bam \
    -M "${fname}"_"${target_name}"_"${date}"_markdup_metrics.txt
    
    
    echo "samtools index"
    samtools index \
    -@ "${threads}" \
    "${fname}"_"${target_name}"_"${date}"_markdup.bam
    
    
    # haplotype call
    echo "gatk HaplotypeCaller"
    gatk HaplotypeCaller \
    --native-pair-hmm-threads "${threads}" \
    -R "${REF}" \
    -I "${fname}"_"${target_name}"_"${date}"_markdup.bam \
    -O "${PROJECT_PATH}"/bqsr/"${fname}"_"${target_name}"_"${date}".vcf
    
    
    #Variant Filtration SNPs
    (
        cd "${PROJECT_PATH}"/bqsr || exit
        
        echo "Variant Filtration SNPs_1"
        gatk SelectVariants \
        -R "${REF}" \
        -V "${fname}"_"${target_name}"_"${date}".vcf \
        --select-type-to-include SNP \
        -O "${fname}"_"${target_name}"_"${date}"_snps.vcf
        
        
        gatk VariantFiltration \
        -R "${REF}" \
        -V "${fname}"_"${target_name}"_"${date}"_snps.vcf \
        -O "${fname}"_"${target_name}"_"${date}"_snps_filtered.vcf \
        -filter "QD < 2.0" --filter-name "QD2"       \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 4.0" --filter-name "SOR4"     \
        -filter "FS > 60.0" --filter-name "FS60"     \
        -filter "MQ < 40.0" --filter-name "MQ40"     \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
        
        
        gatk SelectVariants \
        --exclude-filtered \
        -V "${fname}"_"${target_name}"_"${date}"_snps_filtered.vcf \
        -O "${fname}"_"${target_name}"_"${date}"_snps_bqsr.vcf
        echo "Variant Filtration SNPs_1 done"
    )&
    slect_SNPs=$!
    
    
    #Variant Filtration indels
    cd "${PROJECT_PATH}"/bqsr || exit
    
    echo "Variant Filtration indels_1"
    gatk SelectVariants \
    -R "${REF}" \
    -V "${fname}"_"${target_name}"_"${date}".vcf \
    --select-type-to-include INDEL \
    -O "${fname}"_"${target_name}"_"${date}"_indels.vcf
    
    gatk VariantFiltration \
    -R "${REF}" \
    -V "${fname}"_"${target_name}"_"${date}"_indels.vcf \
    -O "${fname}"_"${target_name}"_"${date}"_indels_filtered.vcf \
    -filter "QD < 2.0" --filter-name "QD2"       \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200"   \
    -filter "SOR > 10.0" -filter-name "SOR10"    \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
    
    gatk SelectVariants \
    --exclude-filtered \
    -V "${fname}"_"${target_name}"_"${date}"_indels_filtered.vcf \
    -O "${fname}"_"${target_name}"_"${date}"_indels_bqsr.vcf
    echo "Variant Filtration indels_1 done"
    
    wait ${slect_SNPs}
    
    
    #BQSR
    cd "${PROJECT_PATH}"/bam || exit
    
    echo "gatk BaseRecalibrator"
    gatk BaseRecalibrator \
    -R "${REF}" \
    -I "${fname}"_"${target_name}"_"${date}"_markdup.bam \
    --known-sites "${PROJECT_PATH}"/bqsr/"${fname}"_"${target_name}"_"${date}"_snps_bqsr.vcf \
    --known-sites "${PROJECT_PATH}"/bqsr/"${fname}"_"${target_name}"_"${date}"_indels_bqsr.vcf \
    -O "${fname}"_"${target_name}"_"${date}"_recal_data.table
    
    echo "gatk ApplyBQSR"
    gatk ApplyBQSR \
    -R "${REF}" \
    -I "${fname}"_"${target_name}"_"${date}"_markdup.bam \
    -bqsr "${fname}"_"${target_name}"_"${date}"_recal_data.table \
    -O "${fname}"_"${target_name}"_"${date}"_bqsr.bam
    
    
    #Variant Call
    cd "${PROJECT_PATH}"/bam || exit
    
    echo "gatk HaplotypeCaller"
    gatk HaplotypeCaller \
    --native-pair-hmm-threads "${threads}" \
    -R "${REF}" \
    -I "${fname}"_"${target_name}"_"${date}"_bqsr.bam \
    -O "${PROJECT_PATH}"/vcf/"${fname}"_"${target_name}"_"${date}".vcf
    
    
    #check BQSR
    echo "gatk BaseRecalibrator"
    gatk BaseRecalibrator \
    -R "${REF}" \
    -I "${fname}"_"${target_name}"_"${date}"_bqsr.bam \
    --known-sites "${PROJECT_PATH}"/bqsr/"${fname}"_"${target_name}"_"${date}"_snps_bqsr.vcf \
    --known-sites "${PROJECT_PATH}"/bqsr/"${fname}"_"${target_name}"_"${date}"_indels_bqsr.vcf \
    -O "${fname}"_"${target_name}"_"${date}"_recal_data.table.2
    
    
    #Variant Filtration SNPs
    (
        cd "${PROJECT_PATH}"/vcf || exit
        
        echo "Variant Filtration SNPs_2"
        gatk SelectVariants \
        -R "${REF}" \
        -V "${fname}"_"${target_name}"_"${date}".vcf \
        --select-type-to-include SNP \
        -O "${fname}"_"${target_name}"_"${date}"_snps.vcf
        
        gatk VariantFiltration \
        -R "${REF}" \
        -V "${fname}"_"${target_name}"_"${date}"_snps.vcf \
        -O "${fname}"_"${target_name}"_"${date}"_snps_filtered.vcf \
        -filter "QD < 2.0" --filter-name "QD2"       \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 4.0" --filter-name "SOR4"     \
        -filter "FS > 60.0" --filter-name "FS60"     \
        -filter "MQ < 40.0" --filter-name "MQ40"     \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
        echo "Variant Filtration SNPs_2 done"
    )&
    slect_SNPs_2=$!
    
    
    #Variant Filtration indels
    cd "${PROJECT_PATH}"/vcf || exit
    
    echo "Variant Filtration indels_2"
    gatk SelectVariants \
    -R "${REF}" \
    -V "${fname}"_"${target_name}"_"${date}".vcf \
    --select-type-to-include INDEL \
    -O "${fname}"_"${target_name}"_"${date}"_indels.vcf
    
    gatk VariantFiltration \
    -R "${REF}" \
    -V "${fname}"_"${target_name}"_"${date}"_indels.vcf \
    -O "${fname}"_"${target_name}"_"${date}"_indels_filtered.vcf \
    -filter "QD < 2.0" --filter-name "QD2"       \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200"   \
    -filter "SOR > 10.0" -filter-name "SOR10"    \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
    echo "Variant Filtration indels_2 done"
    
    wait ${slect_SNPs_2}
    
    
    # merge SNPs and indels
    cd "${PROJECT_PATH}"/vcf || exit
    
    echo "merge SNPs and indels"
    bgzip -c \
    "${fname}"_"${target_name}"_"${date}"_snps_filtered.vcf \
    > "${fname}"_"${target_name}"_"${date}"_snps_filtered.vcf.gz
    
    bgzip -c \
    "${fname}"_"${target_name}"_"${date}"_indels_filtered.vcf \
    > "${fname}"_"${target_name}"_"${date}"_indels_filtered.vcf.gz
    
    
    bcftools index "${fname}"_"${target_name}"_"${date}"_snps_filtered.vcf.gz
    
    bcftools index "${fname}"_"${target_name}"_"${date}"_indels_filtered.vcf.gz
    
    
    bcftools concat -a -d all \
    -o "${fname}"_"${target_name}"_"${date}"_snps_indels_filtered.vcf \
    -O v \
    "${fname}"_"${target_name}"_"${date}"_snps_filtered.vcf.gz \
    "${fname}"_"${target_name}"_"${date}"_indels_filtered.vcf.gz
    echo "merge SNPs and indels done"
    
    # vcf to tsv
    echo "gatk VariantsToTable"
    gatk VariantsToTable \
    -V "${fname}"_"${target_name}"_"${date}"_snps_indels_filtered.vcf \
    -F CHROM \
    -F POS \
    -F REF \
    -F ALT \
    -F FILTER \
    -GF GT \
    -GF AD \
    -GF DP \
    -GF GQ \
    -GF PL \
    -O "${fname}"_"${target_name}"_"${date}"_snps_indels_filtered.tsv \
    --show-filtered
    
    
    
    # make fasta
    echo "bgzip"
    bgzip -c \
    "${fname}"_"${target_name}"_"${date}"_snps_indels_filtered.vcf \
    > "${fname}"_"${target_name}"_"${date}"_snps_indels_filtered.vcf.gz
    
    bcftools index "${fname}"_"${target_name}"_"${date}"_snps_indels_filtered.vcf.gz
    
    echo "bcftools consensus"
    bcftools consensus \
    --include 'FILTER="PASS"' \
    --fasta-ref "${REF}" \
    --output "${PROJECT_PATH}"/fasta/"${fname}"_"${target_name}"_"${date}"_m.fasta \
    -H I \
    "${fname}"_"${target_name}"_"${date}"_snps_indels_filtered.vcf.gz
    
    
    # view target
    echo "samtools faidx"
    cd "${PROJECT_PATH}"/fasta || exit
    samtools faidx "${fname}"_"${target_name}"_"${date}"_m.fasta \
    "${target_region}" \
    > "${fname}"_"${target_name}"_"${date}".fasta
    
    
    rm "${fname}"_"${target_name}"_"${date}"_m.fasta
    rm "${fname}"_"${target_name}"_"${date}"_m.fasta.fai
    
    wait ${flagstat}
    
    
    end_time=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$start_time"
    echo "$end_time"
done
