#! /bin/bash
#$ -S /bin/bash
#$ -cwd

set -e
set -u
set -o pipefail

export MALLOC_ARENA_MAX=2


# 解析するSRA登録番号
SAMPLE="SRR7089217"


# 目的遺伝子情報
TARGET_NAME="V1R2" # target name
TARGET_REGION="NC_031970.2:10560082-10561023" # target region +-100bp base


# その他必要な入力情報
# リファレンスのパス
REF="/home/nakaharu/taki/ref/O_niloticus/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna"
# マッピング用インデックスのパス
BWA_REF="/home/nakaharu/taki/ref/O_niloticus/genome/Onilo_bwa"
# 作業ディレクトリパス
WORK_DIR="/home/nakaharu/taki/test/results"


# 各ツールスレッド数・メモリ数
fasterq_threads="16"
pigz_threads="16"
fastp_threads="6" # ~6 threads
bwa_threads="16"
sort_mem="10G" # ~10G
sort_threads="12" # ~12 threads
index_threads="16"
view_threads="16"
markdup_mem="10G" # ~60G
markdup_core="16" # ~16 cores
call_mem="20G" # ~20G
call_threads="2" # ~2 threads
br_mem="4G" # ~4G
br_threads="2" # ~2 threads
bqsr_mem="2G" # ~2G
bqsr_threads="2" # ~2 threads


# 実行スクリプトのパスを取得
SCRIPT=$(basename "$0")
SCRIPT_PATH=$(cd "$(dirname "$0")"; pwd)/"${SCRIPT}"


# パイプライン
# 開始時間を取得
START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
DATE=$(date "+%Y%m%d%H%M")


# 必要なディレクトリを作成
SAMPLE_DIR="${WORK_DIR}"/"${SAMPLE}"
TARGET_DIR="${WORK_DIR}"/"${SAMPLE}"/"${TARGET_NAME}"
mkdir -p "${SAMPLE_DIR}"/{fastq,bam}
mkdir -p "${TARGET_DIR}"/{lscratch,log,bqsr,vcf,results}



# ログ
LOG_OUT="${TARGET_DIR}"/log/"${DATE}"_stdout.log
LOG_ERR="${TARGET_DIR}"/log/"${DATE}"_stderr.log
LOG_VER="${TARGET_DIR}"/log/"${DATE}"_version.log

exec 2>>"${LOG_ERR}"

echo "${SAMPLE}" | tee -a "${LOG_OUT}"


# 使用するツールのバージョンを取得
prefetch --version > "${LOG_VER}"
{
    fasterq-dump --version
    pigz --version
    fastp --version
    bwa-mem2 version
    samtools --version
    gatk --version
    bgzip --version
    bcftools --version
} >> "${LOG_VER}"


#実行スクリプトを取得
cp "${SCRIPT_PATH}" "${TARGET_DIR}"/log/.


# ショートリードデータ取得
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


#クオリティチェック・トリミング
cd "${SAMPLE_DIR}"/fastq || exit

echo "fastp" | tee -a "${LOG_OUT}"
fastp -w "${fastp_threads}" \
-i "${SAMPLE}"_1.fastq.gz \
-I "${SAMPLE}"_2.fastq.gz \
-o "${SAMPLE}"_cleaned_1_"${DATE}".fastq.gz \
-O "${SAMPLE}"_cleaned_2_"${DATE}".fastq.gz
echo "fastp DONE!" | tee -a "${LOG_OUT}"


# マッピング
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
-o "${SAMPLE_DIR}"/bam/"${SAMPLE}"_"${DATE}".bam
echo "bwa-mem2 | samtools sort DONE!" | tee -a "${LOG_OUT}"


# マッピング結果を取得
(
    cd "${SAMPLE_DIR}"/bam || exit
    
    echo "samtools flagstat" | tee -a "${LOG_OUT}"
    samtools flagstat \
    "${SAMPLE}"_"${DATE}".bam \
    >"${SAMPLE}"_"${DATE}"_flagstat.txt
    echo "samtools flagstat DONE!" | tee -a "${LOG_OUT}"
)&
flagstat=$!

(
    cd "${SAMPLE_DIR}"/bam || exit
    
    echo "samtools coverage" | tee -a "${LOG_OUT}"
    samtools coverage \
    "${SAMPLE}"_"${DATE}".bam \
    -o "${SAMPLE}"_"${DATE}"_stats.tsv
    echo "samtools coverage DONE!" | tee -a "${LOG_OUT}"
)&
coverage=$!


echo "samtools index" | tee -a "${LOG_OUT}"
samtools index \
-@ "${index_threads}" \
"${SAMPLE_DIR}"/bam/"${SAMPLE}"_"${DATE}".bam
echo "samtools index DONE!" | tee -a "${LOG_OUT}"

rm ./*.fastq.gz
rm -rf "${SAMPLE}"


# view target
cd "${SAMPLE_DIR}"/bam || exit

echo "samtools view" | tee -a "${LOG_OUT}"
samtools view \
-o "${TARGET_DIR}"/bqsr/"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".bam \
-h \
-b \
-@ "${view_threads}" \
"${SAMPLE}"_"${DATE}".bam \
"${TARGET_REGION}"
echo "samtools view DONE!" | tee -a "${LOG_OUT}"


# 重複リード検出
cd "${TARGET_DIR}"/bqsr || exit

echo "gatk MarkDuplicatesSpark" | tee -a "${LOG_OUT}"
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
echo "gatk MarkDuplicatesSpark DONE!" | tee -a "${LOG_OUT}"

echo "samtools index" | tee -a "${LOG_OUT}"
samtools index \
-@ "${index_threads}" \
"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_markdup.bam
echo "samtools index DONE!" | tee -a "${LOG_OUT}"


# バリアントコール1回目
cd "${TARGET_DIR}"/bqsr || exit

echo "gatk HaplotypeCaller 1" | tee -a "${LOG_OUT}"
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
echo "gatk HaplotypeCaller 1 DONE!" | tee -a "${LOG_OUT}"


#検出されたSNPのフィルタリング
(
    cd "${TARGET_DIR}"/bqsr || exit
    
    echo "gatk SelectVariants SNPs 1" | tee -a "${LOG_OUT}"
    gatk SelectVariants \
    -R "${REF}" \
    -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
    --select-type-to-include SNP \
    -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs.vcf \
    --verbosity ERROR
    echo "gatk SelectVariants 1 DONE!" | tee -a "${LOG_OUT}"
    
    echo "gatk VariantFiltration SNPs 1" | tee -a "${LOG_OUT}"
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
    echo "gatk VariantFiltration SNPs 1 DONE!" | tee -a "${LOG_OUT}"
    
    echo "gatk SelectVariants filtered SNPs 1" | tee -a "${LOG_OUT}"
    gatk SelectVariants \
    --exclude-filtered \
    -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_filtered.vcf \
    -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_bqsr.vcf \
    --verbosity ERROR
    echo "gatk SelectVariants filtered SNPs 1 DONE!" | tee -a "${LOG_OUT}"
)&
slect_SNPs=$!


#検出されたINDELのフィルタリング
cd "${TARGET_DIR}"/bqsr || exit

echo "gatk SelectVariants INDELs 1" | tee -a "${LOG_OUT}"
gatk SelectVariants \
-R "${REF}" \
-V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
--select-type-to-include INDEL \
-O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs.vcf \
--verbosity ERROR
echo "gatk SelectVariants INDELs 1 DONE!" | tee -a "${LOG_OUT}"

echo "gatk VariantFiltration INDELs 1" | tee -a "${LOG_OUT}"
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
echo "gatk VariantFiltration INDELs 1 DONE!" | tee -a "${LOG_OUT}"

echo "gatk SelectVariants filtered INDELs 1" | tee -a "${LOG_OUT}"
gatk SelectVariants \
--exclude-filtered \
-V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_filtered.vcf \
-O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs_bqsr.vcf \
--verbosity ERROR
echo "gatk SelectVariants filtered INDELs 1 DONE!" | tee -a "${LOG_OUT}"

wait ${slect_SNPs}


#BQSR
cd "${TARGET_DIR}"/bqsr || exit

echo "gatk BaseRecalibrator 1" | tee -a "${LOG_OUT}"
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
echo "gatk BaseRecalibrator 1 DONE!" | tee -a "${LOG_OUT}"

echo "gatk ApplyBQSR" | tee -a "${LOG_OUT}"
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
echo "gatk ApplyBQSR DONE!" | tee -a "${LOG_OUT}"


#バリアントコール2回目
cd "${TARGET_DIR}"/bqsr || exit

echo "gatk HaplotypeCaller 2" | tee -a "${LOG_OUT}"
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
echo "gatk HaplotypeCaller 2 DONE!" | tee -a "${LOG_OUT}"


#BQSR結果確認
(
    cd "${TARGET_DIR}"/bqsr || exit
    
    echo "gatk BaseRecalibrator 2" | tee -a "${LOG_OUT}"
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
    echo "gatk BaseRecalibrator 2 DONE!" | tee -a "${LOG_OUT}"
)&
check_BQSR=$!


#検出されたSNPのフィルタリング
(
    cd "${TARGET_DIR}"/vcf || exit
    
    echo "gatk SelectVariants SNPs 2" | tee -a "${LOG_OUT}"
    gatk SelectVariants \
    -R "${REF}" \
    -V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
    --select-type-to-include SNP \
    -O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs.vcf \
    --verbosity ERROR
    echo "gatk SelectVariants SNPs 2 DONE!" | tee -a "${LOG_OUT}"
    
    echo "gatk VariantFiltration SNPs 2" | tee -a "${LOG_OUT}"
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
    echo "gatk VariantFiltration SNPs 2 DONE!" | tee -a "${LOG_OUT}"
)&
slect_SNPs_2=$!


#検出されたINDELのフィルタリング
cd "${TARGET_DIR}"/vcf || exit

echo "gatk SelectVariants INDELs 2" | tee -a "${LOG_OUT}"
gatk SelectVariants \
-R "${REF}" \
-V "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".vcf \
--select-type-to-include INDEL \
-O "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_INDELs.vcf \
--verbosity ERROR
echo "gatk SelectVariants INDELs 2 DONE!" | tee -a "${LOG_OUT}"

echo "gatk VariantFiltration INDELs 2" | tee -a "${LOG_OUT}"
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
echo "gatk VariantFiltration INDELs 2 DONE!" | tee -a "${LOG_OUT}"

wait ${slect_SNPs_2}


# 最終的なバリアント情報のマージ
cd "${TARGET_DIR}"/vcf || exit

echo "merge SNPs and INDELs" | tee -a "${LOG_OUT}"
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
echo "merge SNPs and INDELs DONE!" | tee -a "${LOG_OUT}"


# make fasta
cd "${TARGET_DIR}"/results || exit

echo "samtools faidx | bcftools consensus" | tee -a "${LOG_OUT}"
bcftools index "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_INDELs_filtered.vcf.gz

samtools faidx \
"${REF}" \
"${TARGET_REGION}" | \
bcftools consensus \
--include 'FILTER="PASS"' \
-H I \
"${SAMPLE}"_"${TARGET_NAME}"_"${DATE}"_SNPs_INDELs_filtered.vcf.gz \
> "${SAMPLE}"_"${TARGET_NAME}"_"${DATE}".fasta

echo "samtools faidx | bcftools consensus DONE!" | tee -a "${LOG_OUT}"


wait ${flagstat}
wait ${coverage}
wait ${check_BQSR}

rm -rf "${TARGET_DIR}"/lscratch


# 終了時間の取得
END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
echo "$START_TIME" | tee -a "${LOG_OUT}"
echo "$END_TIME" | tee -a "${LOG_OUT}"
