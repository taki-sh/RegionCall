# RegionCall

## はじめに
ショートリードデータの取得からマッピング、バリアントコール、目的遺伝子配列の単離までを行うパイプラインです。
- スクリプトはNikaidoLab_NAS/共有データ/Bioinformatics/Variant_callにあります。
- 遺伝研スパコンの利用を前提にしています。
- パイプラインの使用は**自己責任**でお願いします。
- マニュアルを読んで、それぞれの目的にあったパイプラインに修正してください。
  - [GATK4 ベストプラクティス](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
  - [NIH GATK4 チュートリアル](https://hpc.nih.gov/training/gatk_tutorial/)

## 必要なもの
- 解析対象のSRA登録番号
- リファレンス配列
- マッピング用インデックス（詳しくは[こちら](https://github.com/bwa-mem2/bwa-mem2)）  
例：`bwa-mem2 index -p Nbri_bwa ~/main/sequence/N_bri/GCF_000239395.1_NeoBri1.0_genomic.fna`
- GATK4用インデックスその１（詳しくは[こちら](https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-)）  
例：`gatk CreateSequenceDictionary -R ~/main/sequence/N_bri/GCF_000239395.1_NeoBri1.0_genomic.fna -O ~/main/sequence/N_bri/GCF_000239395.1_NeoBri1.0_genomic.dict`
- GATK4用インデックスその２（詳しくは[こちら](http://www.htslib.org/doc/samtools-faidx.html)）  
例：`samtools faidx ~/main/sequence/N_bri/GCF_000239395.1_NeoBri1.0_genomic.fna`
## ヘッダー

```bash
#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N test_1
#$ -l d_rt=300:00:00
#$ -l s_rt=300:00:00
#$ -l s_vmem=12G
#$ -l mem_req=12G
#$ -pe def_slot 32
```
遺伝研スパコンにジョブを投入する際に必要な情報です。ジョブに合わせて要求する計算リソースを変更してください。
- [遺伝研スパコン アレイジョブ](https://sc.ddbj.nig.ac.jp/software/grid_engine/batch_jobs)

`#$ -pe def_slot 32` は並列ジョブを実行する際に必要な情報です。
- [遺伝研スパコン パラレルジョブ](https://sc.ddbj.nig.ac.jp/software/grid_engine/parallel_jobs)

```bash
set -e
set -u
set -o pipefail

export MALLOC_ARENA_MAX=2
```
`set -e`  
>スクリプトの途中でエラーになったときに、次に進んでしまうことを防ぎます。  

`set -u`  
>変数の値が設定されていない場合にスクリプトが中止されます。  

`set -o pipefail`  
>パイプの中で実行されているプログラムの1つが失敗した場合にスクリプトが中止されます。
- [バイオインフォマティクスデータスキル](https://www.oreilly.co.jp/books/9784873118635/)(ラボにあります)  

`export MALLOC_ARENA_MAX=2`
>arena機能はJavaプログラムの実行時には不要ですのでJavaプログラムを使う際には環境変数 MALLOC_ARENA_MAX に小さな値を設定してください。  
- [遺伝研スパコン Javaの使い方](https://sc.ddbj.nig.ac.jp/software/java/)  

## SRA 登録番号入力
```bash
sample_list=(
    ERR4022800
)
```
解析したいショートリードデータのSRA登録番号を入れてください。  
複数サンプル入力可能です。

## その他必要な入力情報
```bash
# その他必要な入力情報
# 目的遺伝子名
target_name="V1R2"
# リファレンスにおける目的遺伝子の番地（余裕をもって+-100bpにしておくといいかも）
target_region="NW_006272005.1:5268493-5269634"
# リファレンスのパス
REF=~/main/sequence/N_bri/GCF_000239395.1_NeoBri1.0_genomic.fna
# マッピング用インデックスのパス
BWA_REF=~/main/sequence/N_bri/Nbri_bwa
# マッピング用スレッド数
threads="32"
# 作業ディレクトリパス
WORK_PATH=~/main/results
```
基本的にはコメントの通りです。  
- 目的遺伝子名はファイル名に使用するだけなので、何でも構いません。
- リファレンスにおける目的遺伝子の番地の書き方は[こちら](http://www.htslib.org/doc/samtools-view.html)を参照してください。
- マッピング用インデックスはあらかじめ作成しておく必要があります。
- 作業ディレクトリの中に結果が出力されます。

## Singularity
```bash
# Singularity
shopt -s expand_aliases
alias prefetch="singularity exec /usr/local/biotools/s/sra-tools:2.11.0--pl5321ha49a11a_3 prefetch"
alias fasterq-dump="singularity exec /usr/local/biotools/s/sra-tools:2.11.0--pl5321ha49a11a_3 fasterq-dump"
alias pigz="singularity exec /usr/local/biotools/p/pigz:2.3.4 pigz"
alias fastp="singularity exec /usr/local/biotools/f/fastp:0.23.2--hb7a2d85_2 fastp"
alias bwa-mem2="singularity exec /usr/local/biotools/b/bwa-mem2:2.2--he513fc3_0 bwa-mem2"
alias samtools="singularity exec /usr/local/biotools/s/samtools:1.15--h3843a85_0 samtools"
alias gatk="singularity exec /usr/local/biotools/g/gatk4:4.2.6.1--py36hdfd78af_1 gatk"
alias bgzip="singularity exec /usr/local/biotools/s/samtools:1.15--h3843a85_0 bgzip"
alias bcftools="singularity exec /usr/local/biotools/b/bcftools:1.15--haf5b3da_0 bcftools"
```
遺伝研スパコンでは、解析ソフトウェアのインストールの手間を軽減するために、BioContainers project が作成した Singularity コンテナイメージを、遺伝研スパコンの/usr/local/biotools/ディレクトリ以下に配置してあります。詳しくは[こちら](https://sc.ddbj.nig.ac.jp/software/BioContainers/)  

スクリプトの#Singularityではパイプラインで使用するツールのパスをあらかじめ指定してあります。ツールのバージョンを変更したい場合は、パスを変更してください。

## パイプライン
**必要な入力情報はありませんが、目を通しておいてください。**  
パイプラインの簡単な流れを説明します。  

1. **ショートリードデータ取得**
   - `prefetch` ：SRAファイルをダウンロード
   - `fasterq-dump`：SRAファイルをfastqファイルに変換
   - `pigz`：fastqファイルを圧縮

2. **クオリティチェック・トリミング**
   - `fastp`

3. **マッピング**
   - `bwa-mem2`：マッピング
     - シーケンサーがイルミナではない場合はbamRGのtPLを変更してください。詳しくは[こちら](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)
   - `samtools sort`：bamファイルを作成
   - `samtools flagstat`：マッピング結果を出力

4. **目的遺伝子領域を単離**
   - `samtools index`：bamファイルのインデックスを作成
   - `samtools view`：bamファイルから目的遺伝子領域を単離

5. **重複リード検出**
   - `gatk MarkDuplicatesSpark`

6. **バリアントコール1回目**
   - `gatk HaplotypeCaller`

7. **フィルタリング**
   - `gatk SelectVariants`：バリアントからSNPとINDELを分離
   - `gatk VariantFiltration`：バリアントのフィルタリング
     - フィルタリングの設定は緩めにしてあります。[こちら](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471)を参考に必要に応じて変更してください。

8. **BQSR（Base Quality Score Recalibration）** 詳しくは[こちら](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-)
   - `gatk BaseRecalibrator`：クオリティスコア補正用のモデル構築
   - `gatk ApplyBQSR`：クオリティスコアの補正
9.  **バリアントコール2回目**
10. **BQSR結果確認**
11. **フィルタリング**
12. **最終的なバリアント情報のマージ**
    - `bgzip`：vcfファイルを圧縮
    - `bcftools concat`：SNP情報とINDEL情報をマージ
      - 最終的なバリアント情報を見るときはtsvファイルがおすすめです。
13. **目的遺伝子配列の作成**
    - `bcftools consensus`：バリアント情報をもとにリファレンスを置換
    - `samtools faidx`：置換した配列から目的遺伝子領域を単離




## ディレクトリ構造
- sample
  - bam_ref
    >他の遺伝子を単離する場合は、このディレクトリにあるbamファイルを使用することができます。
  - fastq
    >ショートリードデータのクオリティチェック・トリミング結果があります。
  - target_name
    - bam
    - bqsr
    - fasta
      >目的遺伝子配列があります。
    - vcf
      >バリアント情報があります。


## 作成済みのbamファイルを使用する場合
同じサンプルの別の領域においてバリアントコールを行う場合、作成済みのbamファイルを使用することができます。  
- `# bamファイルパス`の変更  
- `### 作成済みbamファイルを使用する場合はコメントアウト`で挟まれた部分のコードのコメントアウト  

