#!/bin/bash -e

help="Usage: ./TEPIC.sh [-g input fasta file] [-b bed file containing open chromatin regions] [-o prefix of outputfiles] [-p pwms] \n
Optional parameters:\n
[-c number of cores to use (default 1)]\n
[-d bedgraph file containing open chromatin signal, e.g. DNase1-seq]\n
[-a gene annotation file, required to generate the gene view]\n
[-n column in the -b file containg the average per base signal within a peak. If this option is used, the -d option must not be used.]\n
[-w size of the window to be considered to generate gene view (default 50000bp)]\n
[-e flag to be set if exponential decay should not be used]\n
[-l input Hi-C loopfile]\n
[-v size of the window to be considered to search for Hi-C loops around a TSS (default 25000bp)]\n
[-r defines the Hi-C resolution of the loops to use. Loops having the defined resolution should exist in the Hi-C file. r=All is also allowed, as it automatically searches for all available resolutions. (default 5000bp)]\n
[-j flag to be set if exponential decay should not be used in loops]\n
[-k flag to be set if scaling of open chromatin regions inside loop-sites using the loop-count value should be activated]\n
[-s sparse matrix representation output]"


#Initializing parameters
genome=""
regions=""
prefixP=""
cores=1
pwms=""
dnase=""
column=""
annotation=""
window=50000
decay="TRUE"
hicloops=""	# will be a normal fasta later on, TODO: integrate HiCCUPS pipeline
loopwindows=25000
resolution=5000
loopdecay="TRUE"
loopcountscaling="FALSE"
sparsity=0


#Parsing command line
while getopts "g:b:o:p:c:d:a:n:w:e:l:v:r:j:k:sh" o;
do
                    case $o in
                    g)                  genome=$OPTARG;;
                    b)                  regions=$OPTARG;;
                    o)                  prefixP=$OPTARG;;
                    p)                  pwms=$OPTARG;;
                    c)                  cores=$OPTARG;;
                    d)                  dnase=$OPTARG;;
                    a)                  annotation=$OPTARG;;
                    n)                  column=$OPTARG;;
                    w)                  window=$OPTARG;;
                    e)                  decay="FALSE";;
                    l)                  hicloops=$OPTARG;;
                    v)					loopwindows=$OPTARG;;
                    r)					resolution=$OPTARG;;
                    j)					loopdecay="FALSE";;
                    k)					loopcountscaling="TRUE";;
                    s)					sparsity=$OPTARG;;
					h)					echo -e $help
										exit 1;;
                    [?])                echo -e $help
                                        exit 1;;
                    esac
done

if [ $OPTIND -eq 1 ] ;
then
    echo -e $help
    exit 1;
fi

if [ -z "$genome" ] ;
then
	echo Reference genome must be specified using the -g parameter
	exit 1;
fi

if [ -z "$regions" ] ;
then
	echo Open chromatin regions must be specified using the -b parameter
	exit 1;
fi

if [ -z "$prefixP" ] ;
then
	echo Prefix of output files must be specified using the -o parameter
	exit 1;
fi

if [ -z "$pwms" ] ;
then
	echo PWMs must be specified using the -p parameter
	exit 1;
fi

if [ -n "$dnase" ] && [ -n "$column" ]
then
	echo The options -d and -n are mutually exclusive
	exit 1;
fi

re='^[0-9]+$'
if ! [[ $cores =~ $re ]] ;
then
   echo "Error: Number of cores is not a number" >&2;
   exit 1;
fi

if ! [[ $window =~ $re ]] ;
then
   echo "Error: Window size is not a number" >&2;
   exit 1;
fi

if [ $window -gt 500000 ] ;
then
   echo "WARNING: We recommend to use a smaller value for the TSS-window. Proceeding anyway...";
fi

if ! [[ $loopwindows =~ $re ]] ;
then
   echo "Error: Loopwindow size is not a number" >&2;
   exit 1;
fi

if [ $loopwindows -gt 500000 ] ;
then
   echo "WARNING: We recommend to use a smaller value for Loop-window. Proceeding anyway...";
fi

if ! [[ $resolution =~ $re ]] ;
then
	if [ "$x" == "valid" ];
   echo "Error: Hi-C resolution is not a number" >&2;
   exit 1;
fi

if ! [[ $sparsity =~ $re ]] ;
then
   echo "Error: Sparsity is not a number" >&2;
   exit 1;
fi


#Building prefix
d=$(date +%D)
d=`echo $d | sed 's/\//\_/g'`
t=$(date +%T | sed 's/:/_/g')

prefix=$prefixP"_TEPIC_"${d}"_"${t}
filteredRegions=`echo $regions | awk -F ".bed" '{print $1}'`


#Generating name of the fasta file containing the overlapping regions
openRegionSequences=${prefix}.OpenChromatin.fasta
metadatafile=${prefix}.amd.tsv


#Creating metadata file
touch $metadatafile
echo "[Description]" >> $metadatafile
echo "process	TEPICv1" >> $metadatafile
echo "run_by_user	"$USER >> $metadatafile
echo "date	"$d >> $metadatafile
echo "time	"$t >> $metadatafile
echo "analysis_id	"$prefix >> $metadatafile
echo "" >> $metadatafile
echo "[Inputs]" >> $metadatafile
echo "region_file	"$regions >> $metadatafile
echo "" >> $metadatafile
echo "[References]" >> $metadatafile
echo "genome_reference	"$genome >> $metadatafile
echo "pwms	"$pwms >> $metadatafile
if [ -n "$genome_annotation" ];
then
echo "genome_annotation	"$annotation>> $metadatafile
fi
if [ -n "$dnase" ];
then 
	echo "signale_file	"$dnase >> $metadatafile
fi
echo "" >> $metadatafile
echo "[Outputs]" >> $metadatafile
echo "regions_filtered	"$filteredRegions >> $metadatafile
echo "affinity_file_peak_view	"$prefix"_Affinity.txt" >> $metadatafile
echo "affinity_file_gene_view_filtered	"$prefix"_Affinity_Gene_View_Filtered.txt" >> $metadatafile
if [ -n "$dnase" ];
then 
	echo "signal_scaling_factors	"$prefix"_Peak_Coverage.txt" >> $metadatafile
	echo "scaled affinity_peak_view	"$prefix"_Scaled_Affinity.txt" >> $metadatafile
	echo "scaled_affinity_gene_view_filtered	"$prefix"_Scaled_Affinity_Gene_View_Filtered.txt" >> $metadatafile
fi
if [ -n "$column" ];
then 
	echo "signal_scaling_factors	"$prefix"_NOME_average.txt" >> $metadatafile
	echo "scaled affinity_peak_view	"$prefix"_Scaled_Affinity.txt" >> $metadatafile
	echo "scaled_affinity_gene_view_filtered	"$prefix"_Scaled_Affinity_Gene_View_Filtered.txt" >> $metadatafile
fi

echo "" >> $metadatafile
echo "[Parameters]" >> $metadatafile
echo "SampleID:	"$prefixP >> $metadatafile
echo "cores	"$cores >> $metadatafile
if [ -n "$column" ];
then
	echo "column	"$column >> $metadatafile
fi
if [ -n "$annotation" ];
then
echo "window	"$window >> $metadatafile
echo "decay	"$decay >> $metadatafile
echo "Hi-C file "$hicloops >> $metadatafile
echo "sparsity	"$sparsity >> $metadatafile
fi
echo "" >> $metadatafile
echo "[Metrics]" >> $metadatafile
numReg=`wc -l $regions | cut -f 1 -d " "`
echo "Number of analysed regions	"$numReg >> $metadatafile
numMat=`grep ">" $pwms | wc -l`
echo "Number of considered pwms	"$numMat >> $metadatafile 

#Preprocessing
echo "Preprocessing region file"
python removeInvalidGenomicPositions.py $regions
sort -s -V -k1,1 -k2,2 -k3,3 ${filteredRegions}_Filtered_Regions.bed > ${filteredRegions}_sorted.bed
rm ${filteredRegions}_Filtered_Regions.bed
echo "Runnig bedtools"
#Run bedtools to get a fasta file containing the sequence data for predicted open chromatin regions contained in the bedfile
bedtools getfasta -fi $genome -bed ${filteredRegions}_sorted.bed -fo $openRegionSequences


echo "Converting invalid characters"
#Remove R and Y from the sequence
python convertInvalidCharacterstoN.py $openRegionSequences $prefixP-FilteredSequences.fa


#Use TRAP to compute transcription factor affinities to the above extracted sequences
affinity=${prefix}_Affinity.txt
echo "Starting TRAP"
R3script TRAP.R3script $prefixP-FilteredSequences.fa ${affinity}_temp $cores $pwms


#Computing DNase Coverage in Peak regions
if [ -n "$dnase" ];
then 
	sort -s -V -k1,1 -k2,2 -k3,3 $dnase > ${dnase}_sorted
	python computeDNaseCoverage.py ${dnase}_sorted ${filteredRegions}_sorted.bed > ${prefix}_Peak_Coverage.txt
	rm ${dnase}_sorted
	python scaleAffinity.py ${prefix}_Peak_Coverage.txt ${affinity}_temp > ${prefix}_Scaled_Affinity_temp.txt
fi

if [ -n "${column}" ] ;
then
	cut -f ${column} ${filteredRegions}_sorted.bed > ${prefix}_NOMe_average.txt
	python scaleAffinity.py ${prefix}_NOMe_average.txt ${affinity}_temp > ${prefix}_Scaled_Affinity_temp.txt
fi	

echo "Filter regions that could not be annotated."
python filterInvalidRegions.py ${affinity}_temp $affinity
if [ -n "$dnase" ] ||  [ -n "$column" ];
then
	python filterInvalidRegions.py ${prefix}_Scaled_Affinity_temp.txt ${prefix}_Scaled_Affinity.txt
	rm ${prefix}_Scaled_Affinity_temp.txt
fi


#If an annotation file is provied, the gene view is generated
if [ -n "$annotation" ]; 
then
	echo "Generating gene scores"
	if [ -n "$dnase" ] ||  [ -n "$column" ];
	then
		if [ -n "$hicloops"];
		then
			python annotateTSS.py ${annotation} ${affinity}  "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--signalScale" ${prefix}_Scaled_Affinity.txt "--loopfile" $hicloops "--loopwindows" $loopwindows "--resolution" $resolution "--loopdecay" $loopdecay "--loopcountscaling" $loopcountscaling "--sparseRep" $sparsity
		else
			echo "Configuration: Scaled affinites, no Hi-C data."
			python annotateTSS.py ${annotation} ${affinity}  "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--signalScale" ${prefix}_Scaled_Affinity.txt "--sparseRep" $sparsity
		fi
	else
		if [ -n "$hicloops"];
		then
			python annotateTSS.py ${annotation} ${affinity}  "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--loopfile" $hicloops "--loopwindows" $loopwindows "--resolution" $resolution "--loopdecay" $loopdecay "--loopcountscaling" $loopcountscaling "--sparseRep" $sparsity
		else
			python annotateTSS.py ${annotation} ${affinity}  "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--sparseRep" $sparsity
		fi
	fi
	
	
	#Creating files containing only genes for which TF predictions are available
	echo "Filter genes for which no information is available."
	if [ "$decay" == "TRUE" ];
	then
		if [ -n "$dnase" ] || [ -n "$column" ];
		then
			python filterGeneView.py ${prefix}_Decay_Scaled_Affinity_Gene_View.txt
			rm ${prefix}_Decay_Scaled_Affinity_Gene_View.txt
		fi
			python filterGeneView.py ${prefix}_Decay_Affinity_Gene_View.txt
			rm ${prefix}_Decay_Affinity_Gene_View.txt
		
	else
		if [ -n "$dnase" ] ||  [ -n "$column" ];
		then
			python filterGeneView.py ${prefix}_Scaled_Affinity_Gene_View.txt
			rm ${prefix}_Scaled_Affinity_Gene_View.txt
		fi
		python filterGeneView.py ${prefix}_Affinity_Gene_View.txt
		rm ${prefix}_Affinity_Gene_View.txt
	fi
fi


#Delete temporary files generated for TRAP
rm $openRegionSequences
rm ${prefixP}-FilteredSequences.fa
rm ${affinity}_temp
