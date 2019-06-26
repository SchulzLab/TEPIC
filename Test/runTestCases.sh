echo TestV1: Window 3kb - No annotation
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V1 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -w 3000
echo ""
echo TestV2: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V2 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000
echo ""
echo TestV3: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V3 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE
echo ""
echo TestV4: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V4 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u
echo ""
echo TestV5: Windows 3kb - Annotation - No Decay - Not Length Normalised - No Peak Features 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V5 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -l -u
echo ""
echo TestV6: Windows 3kb - Annotation - No Decay - Not Length Normalised - No Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V6 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -l -u -n 4 -x
echo ""
echo TestV7: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V7 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -n 4 -x
echo ""
echo TestV8: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V8 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4 -x
echo ""
echo TestV9: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V9 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -n 4 -x
echo ""
echo TestV10: Windows 3kb - Annotation - No Decay - Not Length Normalised - No Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V10 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -l -u -n 4
echo ""
echo TestV11: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V11 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4
echo ""
echo TestV12: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V12 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -n 4
echo ""
echo TestV13: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V13 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -n 4
echo ""
echo TestV14: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gzip
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V14 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -z
echo ""
echo TestV15: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - reduced peak set 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V15 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -f example_annotation.gtf
echo ""
echo TestV16: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V16 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y 
echo ""
echo TestV17: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body - reduced peak set 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V17 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y -f example_annotation.gtf
echo ""
echo TestV18:  Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body - reduced peak set - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V18 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y -f example_annotation.gtf -n 4
echo ""
echo TestV19:  Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body - reduced peak set - Scaling original 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V19 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y -f example_annotation.gtf -n 4 -x
echo ""
echo TestV20: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - Considering motif length
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V20 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -m ../PWMs/2.0/human_jaspar_hoc_kellis_Length.txt
echo ""
echo TestV21: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - Considering motif length - Chr prefix in the reference genome
bash ../Code/TEPIC.sh -g example_sequence_chr.fa -b example_regions.bed  -o Test_V21 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -m ../PWMs/2.0/human_jaspar_hoc_kellis_Length.txt -j
echo ""
echo TestV22: Windows 3kb - Annotation -Decay - Length Normalised - Peak Features - Compute discrete scoring using provided background regions
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V22 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -k background_regions.bed 
echo ""
echo TestV23: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features - Signal Feature - File with only peak features
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V23 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4 -q TRUE
echo ""
echo TestV24: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features - transcript based annotation
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V24 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4 -t
echo ""
echo TestV25: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features - transcript based annotation - reduced annotation
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V25 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM -f example_annotation.gtf  -a example_annotation.gtf  -w 3000 -e FALSE -n 4 -t
echo ""
echo TestV26: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features only - Chromatin conformation capture data 
bash ../Code/TEPIC.sh -g example_sequence.fa -b Toy_regions.bed  -o Test_V26 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM -a Toy_Annotation.gtf  -w 400 -h Toy_Loops.csv -s 500 -n 4
echo ""
