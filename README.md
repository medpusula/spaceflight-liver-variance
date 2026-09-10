# Uzay Uçuşunun Fare Karaciğerinde Gen İfadesi Varyansı Üzerindeki Etkisi

TÜBİTAK 2204-A Lise Öğrencileri Araştırma Projeleri Yarışması kapsamında hazırlanmıştır.

## Bu kod ne yapıyor?

NASA GeneLab'den alınan RNA-seq verisinde (OSD-379, RR-8/RRRM-1 görevi, fare karaciğeri), uzay uçuşunun (FLT) gen ifadesi varyansını yer kontrolüne (GC) kıyasla artırıp artırmadığını test eder.

## Yöntem

1. Düşük ifadeli genler filtrelenir (ortalama < 10 okuma)
2. Veri log2(x+1) dönüşümüne tabi tutulur
3. Her gen için Levene testi (medyan merkezli) ile varyans farkı test edilir
4. Genom çapı toplu etki, 10.000 tekrarlı bir permütasyon testiyle değerlendirilir
5. Benjamini-Hochberg (FDR) düzeltmesi uygulanır

## Nasıl çalıştırılır?

1. `GLDS-379_rna_seq_Normalized_Counts_rRNArm_GLbulkRNAseq.csv` dosyasını NASA OSDR'den indirin: https://osdr.nasa.gov/bio/repo/data/studies/OSD-379
2. Dosyayı bu script ile aynı klasöre koyun
3. `pip install pandas numpy scipy statsmodels`
4. `python variance_analysis.py`

## Veri kaynağı

NASA GeneLab Open Science Data Repository, OSD-379 (GLDS-379):
https://osdr.nasa.gov/bio/repo/data/studies/OSD-379


## GSEA (gen seti zenginleştirme analizi)

`gsea_analysis.py`, tekil gen anlamlılığı FDR düzeltmesiyle kaybolduğu için, tamamlayıcı bir gen seti düzeyinde analiz çalıştırır. Script, gen seti dosyasındaki Ensembl ID versiyon eklerini otomatik temizler; ayrı bir ön-işleme adımı gerekmez.

1. Gen seti dosyasını (GO Biyolojik Süreç, fare, Ensembl ID) indirin ve bu script ile aynı klasöre koyun (dosya adını değiştirmeyin):

https://raw.githubusercontent.com/ELTEbioinformatics/GMT_files_for_mulea/main/GMT_files/Mus_musculus_10090/GO_BP_Mus_musculus_EnsemblID.gmt

2. Çalıştırmadan önce: `pip install gseapy`
3. `python gsea_analysis.py`

## Lisans

Bu kod, TÜBİTAK 2204-A proje raporunun tekrarlanabilirliğini desteklemek amacıyla açık erişime sunulmuştur.
