EN
# RegionCall
A pipeline that automates the process from downloading WGS data, mapping, variant calling, to extraction of target region sequences.

## Requirements
- SRA Toolkit (verified with 2.11.0)
- pigz (verified with 2.3.4)
- fastp (verified with 0.23.2)
- bwa-mem2 (verified with 2.2)
- samtools (verified with 1.15)
- bcftools (verified with 1.15)
- gatk (verified with 4.3.0.0)
- genozip (verified with 14.0.34)

## Usage
### Running locally
Create an input file:
```
# List of SRA accession numbers to be analyzed
SAMPLE_LIST=(
SRR0000000
SRR0000001
...
)

# Define lists of target names and regions
TARGET_NAME_LIST=(
hoge1
hoge2
...
)

TARGET_REGION_LIST=(
chr1:1-100
chr2:300-500
...
)

# Other required input information
REF="where/to/path" # Path to reference
BWA_REF="where/to/path" # Path to mapping index
WORK_DIR="where/to/path" # Path to working directory


# Number of threads and memory for each tool
fasterq_threads="16" # fasterq_dump
pigz_threads="16" # pigz
fastp_threads="6" # fastp ~6 threads
bwa_threads="16" # bwa-mem2
sort_mem="10G" # samtool sort ~10G
sort_threads="12" # samtools sort ~12 threads
index_threads="16" # samtools index
genozip_threads="16" # genozip

# For parallel processing
view_threads="4" # samtools view
markdup_mem="2G" # MarkDuplicatesSpark mark ~60G
markdup_core="4" # MarkDuplicatesSpark ~16 cores
call_mem="4G" # HaplotypeCaller ~20G
call_threads="2" # HaplotypeCaller ~2 threads
br_mem="4G" # BaseRecalibrator ~4G
br_threads="2" # BaseRecalibrator ~2 threads
bqsr_mem="2G" # ApplyBQSR ~2G
bqsr_threads="2" # ApplyBQSR ~2 threads
```

Usage:
```bash
$ bash RegionCall_local.sh RegionCall_input.txt
```

### Running on the DDBJ supercomputer
*Please write the input information directly into the script.
```bash
$ qsub_beta RegionCall_DDBJ.sh
```

## License
MIT


JN

WGSデータの取得からマッピング、バリアントコール、目的領域配列の抽出までを自動で行うパイプライン

## 必要なもの
- SRA Toolkit (2.11.0で動作確認)
- pigz (2.3.4で動作確認)
- fastp (0.23.2で動作確認)
- bwa-mem2 (2.2で動作確認)
- samtools (1.15で動作確認)
- bcftools (1.15で動作確認)
- gatk (4.3.0.0で動作確認)
- genozip (14.0.34で動作確認)

## 使い方
### ローカルで実行する場合
インプットファイルを作成する
```
# List of SRA accession numbers to be analyzed
SAMPLE_LIST=(
SRR0000000
SRR0000001
...
)

# Define lists of target names and regions
TARGET_NAME_LIST=(
hoge1
hoge2
...
)

TARGET_REGION_LIST=(
chr1:1-100
chr2:300-500
...
)

# Other required input information
REF="where/to/path" # Path to reference
BWA_REF="where/to/path" # Path to mapping index
WORK_DIR="where/to/path" # Path to working directory


# Number of threads and memory for each tool
fasterq_threads="16" # fasterq_dump
pigz_threads="16" # pigz
fastp_threads="6" # fastp ~6 threads
bwa_threads="16" # bwa-mem2
sort_mem="10G" # samtool sort ~10G
sort_threads="12" # samtools sort ~12 threads
index_threads="16" # samtools index
genozip_threads="16" # genozip

# For parallel processing
view_threads="4" # samtools view
markdup_mem="2G" # MarkDuplicatesSpark mark ~60G
markdup_core="4" # MarkDuplicatesSpark ~16 cores
call_mem="4G" # HaplotypeCaller ~20G
call_threads="2" # HaplotypeCaller ~2 threads
br_mem="4G" # BaseRecalibrator ~4G
br_threads="2" # BaseRecalibrator ~2 threads
bqsr_mem="2G" # ApplyBQSR ~2G
bqsr_threads="2" # ApplyBQSR ~2 threads
```

使用法
```bash
$ bash RegionCall_local.sh RegionCall_input.txt
```



### 遺伝研スパコンで実行する場合
*スクリプトに直接入力情報を書き込んでください。
```bash
$ qsub_beta RegionCall_DDBJ.sh
```

## ライセンス
MIT

