#! /bin/bash


set -e
set -u
set -o pipefail

export MALLOC_ARENA_MAX=2


# Get the path of the text file from the argument
file_path="$1"

# Arrays for variable assignment
SAMPLE_LIST=()
TARGET_NAME_LIST=()
TARGET_REGION_LIST=()

# Other variables
REF=""
BWA_REF=""
WORK_DIR=""
fasterq_threads=""
pigz_threads=""
fastp_threads=""
bwa_threads=""
sort_mem=""
sort_threads=""
index_threads=""
genozip_threads=""
view_threads=""
markdup_mem=""
markdup_core=""
call_mem=""
call_threads=""
br_mem=""
br_threads=""
bqsr_mem=""
bqsr_threads=""

# Initial value for variable name
current_variable=""

# Read the text file line by line
while IFS= read -r line || [[ -n "$line" ]]
do
  # Extract lines for array elements
  if [[ "$line" =~ ^[[:space:]]*([[:alnum:]_]+)=\($ ]]; then
    # Set the current variable name
    current_variable="${BASH_REMATCH[1]}"
    continue
  fi

  # Detect the end line of array elements
  if [[ "$line" =~ ^[[:space:]]*\)$ ]]; then
    # Reset the current variable name
    current_variable=""
    continue
  fi

  # Extract lines for array elements
  if [[ "$current_variable" && "$line" =~ ^[[:space:]]*([[:alnum:]_.:-]+)$ ]]; then
    # Add the element to the current variable
    if [[ "$current_variable" == "SAMPLE_LIST" ]]; then
      SAMPLE_LIST+=("${BASH_REMATCH[1]}")
    elif [[ "$current_variable" == "TARGET_NAME_LIST" ]]; then
      TARGET_NAME_LIST+=("${BASH_REMATCH[1]}")
    elif [[ "$current_variable" == "TARGET_REGION_LIST" ]]; then
      TARGET_REGION_LIST+=("${BASH_REMATCH[1]}")
    fi
  fi

  # Extract lines for other variables
  if [[ "$line" =~ ^[[:space:]]*([[:alnum:]_]+)=\"(.*)\"[[:space:]]*#.*$ ]]; then
    # Assign values to variables
    variable_name="${BASH_REMATCH[1]}"
    variable_value="${BASH_REMATCH[2]}"
    declare "$variable_name=$variable_value"
  fi

done < <(printf '%s\n' "$(cat "$file_path")")



# If the length of TARGET_NAME and TARGET_REGION is not equal, output an error message and exit.
if [ ${#TARGET_NAME_LIST[@]} -ne ${#TARGET_REGION_LIST[@]} ]; then
    echo "Error: TARGET_NAME and TARGET_REGION must have the same length." | tee -a "${LOG_OUT}"
    exit 1
fi


################################################################################################################
################################################################################################################


for SAMPLE in "${SAMPLE_LIST[@]}";
do
    # Get the path of the script
    SCRIPT=$(basename "$0")
    SCRIPT_PATH=$(cd "$(dirname "$0")"; pwd)/"${SCRIPT}"


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
            prefetch --version;
            fasterq-dump --version;
            pigz --version;
            fastp --version;
            bwa-mem2 version;
            samtools --version;
            gatk --version;
            bgzip --version;
            bcftools --version;
        } > "${LOG_VER}"


        # Copy the script
        cp "${SCRIPT_PATH}" "${SAMPLE_DIR}"/log/.


        # Get short read data
        cd "${SAMPLE_DIR}"/fastq || exit

        echo "prefetch" | tee -a "${LOG_OUT}"
        prefetch \
        -O ./ \
        --max-size u \
        "${SAMPLE}" \

        echo "prefetch DONE!" | tee -a "${LOG_OUT}"

        echo "fasterq-dump" | tee -a "${LOG_OUT}"
        fasterq-dump \
        ./"${SAMPLE}" \
        -e "${fasterq_threads}" \
        --split-files \

        echo "fasterq-dump DONE!" | tee -a "${LOG_OUT}"

        echo "pigz" | tee -a "${LOG_OUT}"
        pigz \
        -p "${pigz_threads}" \
        ./*.fastq \

        echo "pigz DONE!" | tee -a "${LOG_OUT}"


        # Quality check and trimming
        cd "${SAMPLE_DIR}"/fastq || exit

        echo "fastp" | tee -a "${LOG_OUT}"
        fastp -w "${fastp_threads}" \
        -i "${SAMPLE}"_1.fastq.gz \
        -I "${SAMPLE}"_2.fastq.gz \
        -o "${SAMPLE}"_cleaned_1_"${DATE}".fastq.gz \
        -O "${SAMPLE}"_cleaned_2_"${DATE}".fastq.gz \

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
        "${SAMPLE}"_cleaned_2_"${DATE}".fastq.gz | \
        samtools sort \
        -T "${SAMPLE_DIR}"/lscratch \
        -m "${sort_mem}" \
        -@ "${sort_threads}" \
        -O BAM \
        -o "${SAMPLE_DIR}"/bam/"${SAMPLE}"_"${DATE}".bam \

        echo "bwa-mem2 | samtools sort DONE!" | tee -a "${LOG_OUT}"


        # Get mapping results
        (
            cd "${SAMPLE_DIR}"/bam || exit

            echo "samtools flagstat" | tee -a "${LOG_OUT}"
            samtools flagstat \
            "${SAMPLE}"_"${DATE}".bam \
            >"${SAMPLE}"_"${DATE}"_flagstat.txt \

            echo "samtools flagstat DONE!" | tee -a "${LOG_OUT}"
        )&

        (
            cd "${SAMPLE_DIR}"/bam || exit

            echo "samtools coverage" | tee -a "${LOG_OUT}"
            samtools coverage \
            "${SAMPLE}"_"${DATE}".bam \
            -o "${SAMPLE}"_"${DATE}"_stats.tsv \

            echo "samtools coverage DONE!" | tee -a "${LOG_OUT}"
        )&


        echo "samtools index" | tee -a "${LOG_OUT}"
        samtools index \
        -@ "${index_threads}" \
        "${SAMPLE_DIR}"/bam/"${SAMPLE}"_"${DATE}".bam \

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
            "${TARGET_REGION}" \


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
                -filter "QD < 2.0" --filter-name "QD2"       \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 4.0" --filter-name "SOR4"     \
                -filter "FS > 60.0" --filter-name "FS60"     \
                -filter "MQ < 40.0" --filter-name "MQ40"     \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                --verbosity ERROR

                gatk SelectVariants \
                --exclude-filtered \
                -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf \
                -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_bqsr.vcf \
                --verbosity ERROR

            )&
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
            -filter "QD < 2.0" --filter-name "QD2"       \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200"   \
            -filter "SOR > 10.0" -filter-name "SOR10"    \
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

            )&
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
                -filter "QD < 2.0" --filter-name "QD2"       \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 4.0" --filter-name "SOR4"     \
                -filter "FS > 60.0" --filter-name "FS60"     \
                -filter "MQ < 40.0" --filter-name "MQ40"     \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                --verbosity ERROR

            )&
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
            -filter "QD < 2.0" --filter-name "QD2"       \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200"   \
            -filter "SOR > 10.0" -filter-name "SOR10"    \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            --verbosity ERROR

            wait ${slect_SNPs_2}


            # Merging final variant information
            cd "${TARGET_DIR}"/vcf || exit

            bgzip -c \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf \
            > "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf.gz

            bgzip -c \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf \
            > "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf.gz


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
            "${TARGET_REGION}" | \
            bcftools consensus \
            --include 'FILTER="PASS"' \
            -H I \
            "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_INDELs_filtered.vcf.gz \
            > "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".fasta

            sed -i "s/^>${TARGET_REGION}/>${SAMPLE}_${TARGET_NAME}/" "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".fasta

            wait ${check_BQSR}

            rm -rf "${TARGET_DIR}"/lscratch

            echo "Finished processing ${TARGET_NAME} (${TARGET_REGION})" | tee -a "${LOG_OUT}"
        } 2>>"${TARGET_LOG_ERR}"
    }


    # Set TARGET_NAME and TARGET_REGION lists as a pair, and perform independent operations on each combination
    for (( i=0; i<${#TARGET_NAME_LIST[@]}; i++ )); do
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
    echo "Compressing BAM file..."| tee -a "${LOG_OUT}"

    genozip \
    --threads "${genozip_threads}" \
    "${SAMPLE}"_"${DATE}".bam

    rm "${SAMPLE}"_"${DATE}".bam
    rm "${SAMPLE}"_"${DATE}".bam.bai

    echo "Compression of BAM files completed."| tee -a "${LOG_OUT}"

    # Get the end time
    END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$START_TIME" | tee -a "${LOG_OUT}"
    echo "$END_TIME" | tee -a "${LOG_OUT}"



done
