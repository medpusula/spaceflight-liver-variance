# Uzay Ucusunun Fare Karacigerinde Gen Ifadesi Varyansi Uzerindeki Etkisi

TUBITAK 2204-A Lise Ogrencileri Arastirma Projeleri Yarismasi kapsaminda hazirlanmistir.

## Ne yapiyor bu kod?

NASA GeneLab'den alinan RNA-seq verisinde (OSD-379, RR-8/RRRM-1 gorevi, fare karacigeri),
uzay ucusunun (FLT) gen ifadesi varyansini yer kontrolune (GC) kiyasla artirip artirmadigini
test eder.

## Yontem

1. Dusuk ifadeli genler filtrelenir (ortalama < 10 okuma)
2. Veri log2(x+1) donusumune tabi tutulur
3. Her gen icin Levene testi (medyan-merkezli) ile varyans farki test edilir
4. Genom-capi toplu etki, 10.000 tekrarli bir permutasyon testiyle degerlendirilir
5. Benjamini-Hochberg (FDR) duzeltmesi uygulanir

## Nasil calistirilir?

1. GLDS-379_rna_seq_Normalized_Counts_rRNArm_GLbulkRNAseq.csv dosyasini
   NASA OSDR'den indirin: https://osdr.nasa.gov/bio/repo/data/studies/OSD-379
2. Dosyayi bu script ile ayni klasore koyun
3. `pip install pandas numpy scipy statsmodels`
4. `python variance_analysis.py`

## Veri kaynagi

NASA GeneLab Open Science Data Repository, OSD-379 (GLDS-379):
https://osdr.nasa.gov/bio/repo/data/studies/OSD-379

## Lisans

Bu kod acik erisimlidir, TUBITAK 2204-A proje raporunun tekrarlanabilirligini
desteklemek amaciyla paylasilmistir.
