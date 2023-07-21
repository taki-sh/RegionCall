#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null
#$ -N RC_p_NW
#$ -l d_rt=72:00:00
#$ -l s_rt=72:00:00
#$ -pe def_slot 16
#$ -l s_vmem=24G
#$ -l mem_req=24G
#$ -t 1-96:1
#$ -tc 25

set -e
set -u
set -o pipefail

export MALLOC_ARENA_MAX=2

# List of SRA accession numbers to be analyzed
SAMPLE_LIST=(
    # SRR13192256
    # SRR13192257
    # SRR13192258
    # SRR16915804
    # SRR13192261
    # ERR583983
    # SRR16915805
    # ERR715541
    # ERR4189826
    # SRR16915803
    # SRR16915798
    # SRR16915796
    # SRR16915807
    # SRR16915806
    # SRR16915797
    # SRR16915802
    # SRR16915808
    # SRR16915800
    # SRR16915801
    # SRR16915799
    ERR10768208
    ERR10768192
    ERR10768247
    ERR10768202
    ERR10789849
    ERR10789854
    ERR10768198
    ERR10768215
    ERR10768334
    ERR10768188
    ERR10789894
    ERR10768288
    ERR10789884
    ERR10768336
    ERR10789835
    ERR10768242
    ERR10768190
    ERR10768312
    ERR10768331
    ERR10768271
    ERR10789853
    ERR10789850
    ERR10768311
    ERR10789848
    ERR10768241
    ERR10768328
    ERR10768304
    ERR10789874
    ERR10768194
    ERR10768339
    ERR10768193
    ERR10768352
    ERR10789911
    ERR10789886
    ERR10789872
    ERR10768324
    ERR10789846
    ERR10789865
    ERR10789901
    ERR10768296
    ERR10789871
    ERR10768218
    ERR10768349
    ERR10768199
    ERR10768261
    ERR10768294
    ERR10768277
    ERR10789909
    ERR10768232
    ERR10789845
    ERR10789917
    ERR580916
    ERR10768280
    ERR10768240
    ERR10768251
    ERR10789862
    ERR10789867
    ERR10768234
    ERR10768282
    ERR10768286
    ERR10768266
    ERR10789905
    ERR10768279
    ERR10789881
    ERR10768299
    ERR10768303
    ERR10789836
    ERR10768265
    ERR10789904
    ERR10768285
    ERR10789843
    ERR10789876
    ERR10768292
    ERR10768305
    ERR10789892
    ERR10768302
    ERR10768221
    ERR10768200
    ERR10789856
    ERR10768212
    ERR10789839
    ERR10768289
    ERR10789875
    ERR580921
    ERR10789903
    ERR10768315
    ERR10768358
    ERR580917
    ERR10768323
    ERR10789890
    ERR10768322
    ERR10768329
    ERR580918
    SRR9675353
    SRR9675354
    SRR9657532

)

# Select the sample for analysis based on the SGE_TASK_ID
SAMPLE="${SAMPLE_LIST[$((--SGE_TASK_ID))]}"

# Define lists of target names and regions
TARGET_NAME_LIST=(V1R1 V1R2 V1R3 V1R4 V1R5 V1R6)
TARGET_REGION_LIST=(
    CM024196.1:6251602-6252571
    CM024196.1:6246835-6247801
    CM024196.1:16709988-16712644
    CM024196.1:16708119-16709293
    JACBYM010000007.1:14908744-14909878
    CM024205.1:20660233-20661157
)

# If the length of TARGET_NAME and TARGET_REGION is not equal, output an error message and exit.
if [ ${#TARGET_NAME_LIST[@]} -ne ${#TARGET_REGION_LIST[@]} ]; then
    echo "Error: TARGET_NAME and TARGET_REGION must have the same length." | tee -a "${LOG_OUT}"
    exit 1
fi

# Other required input information
REF="$HOME/main/sequence/A_citrinellus/genome_masked/GCA_013435755.1_ASM1343575v1_genomic_masked.fasta" # # Path to reference
BWA_REF="$HOME/main/sequence/A_citrinellus/genome_masked/A_cit_m_bwa"                                   # Path to mapping index
WORK_DIR="$HOME/main/RegionCall_parallel/results/ref_new_world_2"                                       # Path to working directory

# Number of threads and memory for each tool
fasterq_threads="16"
pigz_threads="16"
fastp_threads="6" # ~6 threads
bwa_threads="16"
sort_mem="10G"    # ~10G
sort_threads="12" # ~12 threads
index_threads="16"
genozip_threads="16"

# For parallel processing
view_threads="4"
markdup_mem="2G" # ~60G
markdup_core="4" # ~16 cores
call_mem="4G"    # ~20G
call_threads="2" # ~2 threads
br_mem="4G"      # ~4G
br_threads="2"   # ~2 threads
bqsr_mem="2G"    # ~2G
bqsr_threads="2" # ~2 threads

# For using Singularity
shopt -s expand_aliases
alias prefetch="singularity exec /usr/local/biotools/s/sra-tools:2.11.0--pl5321ha49a11a_3 prefetch"
alias fasterq-dump="singularity exec /usr/local/biotools/s/sra-tools:2.11.0--pl5321ha49a11a_3 fasterq-dump"
alias pigz="singularity exec /usr/local/biotools/p/pigz:2.3.4 pigz"
alias fastp="singularity exec /usr/local/biotools/f/fastp:0.23.2--hb7a2d85_2 fastp"
alias bwa-mem2="singularity exec /usr/local/biotools/b/bwa-mem2:2.2--he513fc3_0 bwa-mem2"
alias samtools="singularity exec /usr/local/biotools/s/samtools:1.15--h3843a85_0 samtools"
alias bgzip="singularity exec /usr/local/biotools/s/samtools:1.15--h3843a85_0 bgzip"
alias bcftools="singularity exec /usr/local/biotools/b/bcftools:1.15--haf5b3da_0 bcftools"

################################################################################################################
################################################################################################################

# Get the path of the script
SCRIPT=$(basename "$0")
SCRIPT_PATH=$(
    cd "$(dirname "$0")"
    pwd
)/"${SCRIPT}"

# Pipeline
# Get the start time
START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
DATE=$(date "+%Y%m%d%H%M")

# Create necessary directories
SAMPLE_DIR="${WORK_DIR}"/"${SAMPLE}"
mkdir -p "${SAMPLE_DIR}"/{fastq,bam,log,lscratch}

# Logging
LOG_OUT="${SAMPLE_DIR}"/log/"${DATE}"_stdout.log
LOG_ERR="${SAMPLE_DIR}"/log/"${DATE}"_stderr.log
LOG_VER="${SAMPLE_DIR}"/log/"${DATE}"_version.log

#Common processing
{
    echo "${SAMPLE}" | tee -a "${LOG_OUT}"

    # Get the version of the tools used
    {
        prefetch --version
        fasterq-dump --version
        pigz --version
        fastp --version
        bwa-mem2 version
        samtools --version
        gatk --version
        bgzip --version
        bcftools --version
    } >"${LOG_VER}"

    # Copy the script
    cp "${SCRIPT_PATH}" "${SAMPLE_DIR}"/log/.

    # Get short read data
    cd "${SAMPLE_DIR}"/fastq || exit

    echo "prefetch" | tee -a "${LOG_OUT}"
    prefetch \
        -O ./ \
        --max-size u \
        "${SAMPLE}"

    echo "prefetch DONE!" | tee -a "${LOG_OUT}"

    echo "fasterq-dump" | tee -a "${LOG_OUT}"
    fasterq-dump \
        ./"${SAMPLE}" \
        -e "${fasterq_threads}" \
        --split-files

    echo "fasterq-dump DONE!" | tee -a "${LOG_OUT}"

    echo "pigz" | tee -a "${LOG_OUT}"
    pigz \
        -p "${pigz_threads}" \
        ./*.fastq

    echo "pigz DONE!" | tee -a "${LOG_OUT}"

    # Quality check and trimming
    cd "${SAMPLE_DIR}"/fastq || exit

    echo "fastp" | tee -a "${LOG_OUT}"
    fastp -w "${fastp_threads}" \
        -i "${SAMPLE}"_1.fastq.gz \
        -I "${SAMPLE}"_2.fastq.gz \
        -o "${SAMPLE}"_cleaned_1_"${DATE}".fastq.gz \
        -O "${SAMPLE}"_cleaned_2_"${DATE}".fastq.gz

    echo "fastp DONE!" | tee -a "${LOG_OUT}"

    # Mapping
    cd "${SAMPLE_DIR}"/fastq || exit

    echo "bwa-mem2 | samtools sort" | tee -a "${LOG_OUT}"
    bamRG="@RG\tID:${SAMPLE}\tPL:ILLUMINA\tSM:"${SAMPLE}
    bwa-mem2 mem \
        -v 2 \
        -t "${bwa_threads}" \
        -R "${bamRG}" "${BWA_REF}" \
        "${SAMPLE}"_cleaned_1_"${DATE}".fastq.gz \
        "${SAMPLE}"_cleaned_2_"${DATE}".fastq.gz |
        samtools sort \
            -T "${SAMPLE_DIR}"/lscratch \
            -m "${sort_mem}" \
            -@ "${sort_threads}" \
            -O BAM \
            -o "${SAMPLE_DIR}"/bam/"${SAMPLE}"_"${DATE}".bam

    echo "bwa-mem2 | samtools sort DONE!" | tee -a "${LOG_OUT}"

    # Get mapping results
    (
        cd "${SAMPLE_DIR}"/bam || exit

        echo "samtools flagstat" | tee -a "${LOG_OUT}"
        samtools flagstat \
            "${SAMPLE}"_"${DATE}".bam \
            >"${SAMPLE}"_"${DATE}"_flagstat.txt

        echo "samtools flagstat DONE!" | tee -a "${LOG_OUT}"
    ) &

    (
        cd "${SAMPLE_DIR}"/bam || exit

        echo "samtools coverage" | tee -a "${LOG_OUT}"
        samtools coverage \
            "${SAMPLE}"_"${DATE}".bam \
            -o "${SAMPLE}"_"${DATE}"_stats.tsv

        echo "samtools coverage DONE!" | tee -a "${LOG_OUT}"
    ) &

    echo "samtools index" | tee -a "${LOG_OUT}"
    samtools index \
        -@ "${index_threads}" \
        "${SAMPLE_DIR}"/bam/"${SAMPLE}"_"${DATE}".bam

    echo "samtools index DONE!" | tee -a "${LOG_OUT}"

    rm ./*.fastq.gz
    rm -rf "${SAMPLE}"
    rm -rf "${SAMPLE_DIR}"/lscratch
} 2>>"${LOG_ERR}"

# Define a function to be processed independently.
function process_target {
    local TARGET_NAME=$1
    local TARGET_REGION=$2

    echo "Processing ${TARGET_NAME} (${TARGET_REGION})" | tee -a "${LOG_OUT}"

    TARGET_DIR="${WORK_DIR}"/"${SAMPLE}"/"${TARGET_NAME}"
    mkdir -p "${TARGET_DIR}"/{lscratch,bqsr,vcf,results,log}
    TARGET_LOG_ERR="${TARGET_DIR}/log/${DATE}_${TARGET_NAME}_stderr.log"

    {
        # Extract the necessary regions from the bam file
        cd "${SAMPLE_DIR}"/bam || exit

        samtools view \
            -o "${TARGET_DIR}"/bqsr/"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".bam \
            -h \
            -b \
            -@ "${view_threads}" \
            "${SAMPLE}"_"${DATE}".bam \
            "${TARGET_REGION}"

        # Duplicate read detection
        cd "${TARGET_DIR}"/bqsr || exit

        gatk \
            --java-options \
            "-Djava.io.tmpdir=${TARGET_DIR}/lscratch \
        -Xms${markdup_mem} \
        -Xmx${markdup_mem}" \
            MarkDuplicatesSpark \
            -I "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".bam \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_markdup.bam \
            --spark-runner LOCAL \
            --conf "spark.executor.cores=${markdup_core}" \
            --verbosity ERROR

        samtools index \
            -@ "${index_threads}" \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_markdup.bam

        # First round of variant calling
        cd "${TARGET_DIR}"/bqsr || exit

        gatk \
            --java-options \
            "-Djava.io.tmpdir=${TARGET_DIR}/lscratch \
        -Xms${call_mem} \
        -Xmx${call_mem} \
        -XX:ParallelGCThreads=${call_threads}" \
            HaplotypeCaller \
            -R "${REF}" \
            -I "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_markdup.bam \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
            --verbosity ERROR

        # Filtering of detected SNPs
        (
            cd "${TARGET_DIR}"/bqsr || exit

            gatk SelectVariants \
                -R "${REF}" \
                -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
                --select-type-to-include SNP \
                -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs.vcf \
                --verbosity ERROR

            gatk VariantFiltration \
                -R "${REF}" \
                -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs.vcf \
                -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf \
                -filter "QD < 2.0" --filter-name "QD2" \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 4.0" --filter-name "SOR4" \
                -filter "FS > 60.0" --filter-name "FS60" \
                -filter "MQ < 40.0" --filter-name "MQ40" \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                --verbosity ERROR

            gatk SelectVariants \
                --exclude-filtered \
                -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf \
                -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_bqsr.vcf \
                --verbosity ERROR

        ) &
        slect_SNPs=$!

        # Filtering of detected INDELs
        cd "${TARGET_DIR}"/bqsr || exit

        gatk SelectVariants \
            -R "${REF}" \
            -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
            --select-type-to-include INDEL \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs.vcf \
            --verbosity ERROR

        gatk VariantFiltration \
            -R "${REF}" \
            -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs.vcf \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "SOR > 10.0" -filter-name "SOR10" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            --verbosity ERROR

        gatk SelectVariants \
            --exclude-filtered \
            -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_bqsr.vcf \
            --verbosity ERROR

        wait ${slect_SNPs}

        #BQSR
        cd "${TARGET_DIR}"/bqsr || exit

        gatk \
            --java-options \
            "-Djava.io.tmpdir=${TARGET_DIR}/lscratch \
        -Xms${br_mem} \
        -Xmx${br_mem} \
        -XX:ParallelGCThreads=${br_threads}" \
            BaseRecalibrator \
            -R "${REF}" \
            -I "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_markdup.bam \
            --known-sites "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_bqsr.vcf \
            --known-sites "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_bqsr.vcf \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_recal_data.table \
            --verbosity ERROR

        gatk \
            --java-options \
            "-Djava.io.tmpdir=${TARGET_DIR}/lscratch \
        -Xms${bqsr_mem} \
        -Xmx${bqsr_mem} \
        -XX:ParallelGCThreads=${bqsr_threads}" \
            ApplyBQSR \
            -R "${REF}" \
            -I "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_markdup.bam \
            -bqsr "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_recal_data.table \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_bqsr.bam \
            --verbosity ERROR

        # Second round of variant calling
        cd "${TARGET_DIR}"/bqsr || exit

        gatk \
            --java-options \
            "-Djava.io.tmpdir=${TARGET_DIR}/lscratch \
        -Xms${call_mem} \
        -Xmx${call_mem} \
        -XX:ParallelGCThreads=${call_threads}" \
            HaplotypeCaller \
            -R "${REF}" \
            -I "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_bqsr.bam \
            -O "${TARGET_DIR}"/vcf/"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
            --verbosity ERROR

        # Checking BQSR results
        (
            cd "${TARGET_DIR}"/bqsr || exit

            gatk \
                --java-options \
                "-Djava.io.tmpdir=${TARGET_DIR}/lscratch \
        -Xms${br_mem} \
        -Xmx${br_mem} \
            -XX:ParallelGCThreads=${br_threads}" \
                BaseRecalibrator \
                -R "${REF}" \
                -I "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_bqsr.bam \
                --known-sites "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_bqsr.vcf \
                --known-sites "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_bqsr.vcf \
                -O "${TARGET_DIR}"/results/"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_recal_data.table.2 \
                --verbosity ERROR

        ) &
        check_BQSR=$!

        # Filtering of detected SNPs
        (
            cd "${TARGET_DIR}"/vcf || exit

            gatk SelectVariants \
                -R "${REF}" \
                -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
                --select-type-to-include SNP \
                -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs.vcf \
                --verbosity ERROR

            gatk VariantFiltration \
                -R "${REF}" \
                -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs.vcf \
                -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf \
                -filter "QD < 2.0" --filter-name "QD2" \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 4.0" --filter-name "SOR4" \
                -filter "FS > 60.0" --filter-name "FS60" \
                -filter "MQ < 40.0" --filter-name "MQ40" \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                --verbosity ERROR

        ) &
        slect_SNPs_2=$!

        # Filtering of detected INDELs
        cd "${TARGET_DIR}"/vcf || exit

        gatk SelectVariants \
            -R "${REF}" \
            -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
            --select-type-to-include INDEL \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs.vcf \
            --verbosity ERROR

        gatk VariantFiltration \
            -R "${REF}" \
            -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs.vcf \
            -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "SOR > 10.0" -filter-name "SOR10" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            --verbosity ERROR

        wait ${slect_SNPs_2}

        # Merging final variant information
        cd "${TARGET_DIR}"/vcf || exit

        bgzip -c \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf \
            >"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf.gz

        bgzip -c \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf \
            >"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf.gz

        bcftools index "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf.gz
        bcftools index "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf.gz

        bcftools concat -a -d all \
            -o "${TARGET_DIR}"/results/"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_INDELs_filtered.vcf.gz \
            -O z \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf.gz \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf.gz

        # Creation of consensus sequence
        cd "${TARGET_DIR}"/results || exit

        bcftools index "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_INDELs_filtered.vcf.gz

        samtools faidx \
            "${REF}" \
            "${TARGET_REGION}" |
            bcftools consensus \
                --include 'FILTER="PASS"' \
                -H I \
                "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_INDELs_filtered.vcf.gz \
                >"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".fasta

        sed -i "s/^>${TARGET_REGION}/>${SAMPLE}_${TARGET_NAME}/" "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".fasta

        wait ${check_BQSR}

        rm -rf "${TARGET_DIR}"/lscratch

        echo "Finished processing ${TARGET_NAME} (${TARGET_REGION})" | tee -a "${LOG_OUT}"
    } 2>>"${TARGET_LOG_ERR}"
}

# Set TARGET_NAME and TARGET_REGION lists as a pair, and perform independent operations on each combination
for ((i = 0; i < ${#TARGET_NAME_LIST[@]}; i++)); do
    name=${TARGET_NAME_LIST[$i]}
    region=${TARGET_REGION_LIST[$i]}
    process_target "$name" "$region" &
done

# Waiting for processing
sleep 1
echo "Started all processes in the background." | tee -a "${LOG_OUT}"
echo "Waiting for all processes to finish..." | tee -a "${LOG_OUT}"
wait
echo "All processes have finished." | tee -a "${LOG_OUT}"

# Compressing BAM file
cd "${SAMPLE_DIR}"/bam || exit
echo "Compressing BAM file..." | tee -a "${LOG_OUT}"

genozip \
    --threads "${genozip_threads}" \
    "${SAMPLE}"_"${DATE}".bam

rm "${SAMPLE}"_"${DATE}".bam
rm "${SAMPLE}"_"${DATE}".bam.bai

echo "Compression of BAM files completed." | tee -a "${LOG_OUT}"

# Get the end time
END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
echo "$START_TIME" | tee -a "${LOG_OUT}"
echo "$END_TIME" | tee -a "${LOG_OUT}"
