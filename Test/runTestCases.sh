echo TestV1: Window 3kb - No annotation
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V1 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -w 3000
echo ""
echo TestV2: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V2 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000
echo ""
echo TestV3: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V3 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE
echo ""
echo TestV4: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V4 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u
echo ""
echo TestV5: Windows 3kb - Annotation - No Decay - Not Length Normalised - No Peak Features 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V5 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -l -u
echo ""
echo TestV6: Windows 3kb - Annotation - No Decay - Not Length Normalised - No Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V6 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -l -u -n 4 -x
echo ""
echo TestV7: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V7 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -n 4 -x
echo ""
echo TestV8: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V8 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4 -x
echo ""
echo TestV9: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - Scaling original
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V9 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -n 4 -x
echo ""
echo TestV10: Windows 3kb - Annotation - No Decay - Not Length Normalised - No Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V10 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -l -u -n 4
echo ""
echo TestV11: Windows 3kb - Annotation - No Decay - Length Normalised - Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V11 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4
echo ""
echo TestV12: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V12 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -n 4
echo ""
echo TestV13: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V13 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -n 4
echo ""
echo TestV14: Windows 3kb - Annotation - Decay - Length Normalised - No Peak Features - Sparse Representation
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V14 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -s -u
echo ""
echo TestV15: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Sparse Representation
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V15 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -s -l -u
echo ""
echo TestV16: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Scaling original - Sparse Representation
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V16 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -n 4 -x -s
echo ""
echo TestV17: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gzip
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V17 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -z
echo ""
echo TestV18: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - reduced peak set 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V18 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -f example_annotation.gtf
echo ""
echo TestV19: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V19 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y 
echo ""
echo TestV20: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body - reduced peak set 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V20 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y -f example_annotation.gtf
echo ""
echo TestV21:  Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body - reduced peak set - Signal Feature
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V21 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y -f example_annotation.gtf -n 4
echo ""
echo TestV22:  Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - gene body - reduced peak set - Scaling original 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V22 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -y -f example_annotation.gtf -n 4 -x
echo ""
echo TestV23: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Scaling original - Sparse Representation - gzip 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V23 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -n 4 -x -s -z
echo ""
echo TestV24: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Signal Feature - Sparse Representation - gzip 
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V24 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -n 4 -s -z
echo ""
echo TestV25: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Signal Feature - Sparse Representation - bedGraph file - gzip
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V25 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -d example_regions.bg -s -z
echo ""
echo TestV26: Windows 3kb - Annotation - Decay - Not Length Normalised - No Peak Features - Original Scaling - Sparse Representation - bedGraph File - gzip
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V26 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -l -u -d example_regions.bg -s -z -x
echo ""
echo TestV27: Windows 3kb - Annotation - Decay - Length Normalised - Peak Features - Considering motif length
bash ../Code/TEPIC.sh -g example_sequence.fa -b example_regions.bed  -o Test_V27 -p ../PWMs/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -m ../PWMs/human_jaspar_hoc_kellis_Length.txt
echo ""