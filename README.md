# RegionCall
WGSデータの取得からマッピング、バリアントコール、目的領域配列の抽出までを自動で行うパイプラインです。

## 必要なもの
- SRA Toolkit (2.11.0で動作確認)
- pigz (2.3.4で動作確認)
- fastp (0.23.2で動作確認)
- bwa-mem2 (2.2で動作確認)
- samtools (1.15で動作確認)
- bcftools (1.15で動作確認)
- gatk (4.3.0.0で動作確認)

## 使い方
### ローカルで実行する場合
インプットファイルを作成する
```
```

使用法
```bash
$ bash RegionCall_local.sh RegionCall_input.txt
```



### 遺伝研スパコンで実行する場合

```bash
$ qsub_beta RegionCall_DDBJ.sh
```

## License
MIT

