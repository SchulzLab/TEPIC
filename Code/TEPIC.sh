#!/bin/bash -e

help="Usage: ./TEPIC.sh [-g input fasta file] [-b bed file containing open chromatin regions] [-o prefix of outputfiles] [-p pwms] \n
Optional parameters:\n
[-c number of cores to use (default 1)]\n
[-d bedgraph file containing open chromatin signal, e.g. DNase1-seq]\n
[-a gene annotation file, required to generate the gene view]\n
[-n column in the -b file containg the average per base signal within a peak. If this option is used, the -d option must not be used.]\n
[-w size of the window to be considered to generate gene view (default 50000bp)]\n
[-e flag to be set if exponential decay should not be used]\n
[-s sparse matrix representation]"

#Initialising parameters
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
sparsity=0

#Parsing command line
while getopts "g:b:o:c:p:d:n:a:w:s:eh" o;
do
                    case $o in
                    g)                  genome=$OPTARG;;
                    b)                  regions=$OPTARG;;
                    o)                  prefixP=$OPTARG;;
                    c)                  cores=$OPTARG;;
                    p)                  pwms=$OPTARG;;
                    d)                  dnase=$OPTARG;;
                    n)                  column=$OPTARG;;
                    a)                  annotation=$OPTARG;;
                    w)                  window=$OPTARG;;
                    e)                  decay="FALSE";;
	s)	sparsity=$OPTARG;;
	h)	echo -e $help
		exit1;;
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

d=$(date +%D)
d=`echo $d | sed 's/\//\_/g'`
t=$(date +%T | sed 's/:/_/g')

prefix=$prefixP"_TEPIC_"${d}"_"${t}
filteredRegions=`echo $regions | awk -F ".bed" '{print $1}'`
#Generating name of the fasta file containing the overlapping regions
openRegionSequences=${prefix}.OpenChromatin.fasta
metadatafile=${prefix}.amd.tsv
#Create metadata file
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
echo "sparsity	"$sparsity >> $metadatafile
fi
echo "" >> $metadatafile
echo "[Metrics]" >> $metadatafile
numReg=`wc -l $regions | cut -f 1 -d " "`
echo "Number of analysed regions	"$numReg >> $metadatafile
numMat=`grep ">" $pwms | wc -l`
echo "Number of considered pwms	"$numMat >> $metadatafile 


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
		python annotateTSS.py ${annotation} ${affinity}  "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--signalScale" ${prefix}_Scaled_Affinity.txt "--sparseRep" $sparsity
	else
		python annotateTSS.py ${annotation} ${affinity}  "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--sparseRep" $sparsity
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
