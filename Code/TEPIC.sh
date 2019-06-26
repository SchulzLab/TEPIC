#!/bin/bash -e

help="TEPIC version 2.1\n
Usage: ./TEPIC.sh [-g input fasta file (in RefSeq format, without \"chr\" prefix] [-b bed file containing open chromatin regions] [-o prefix of outputfiles] [-p psems] \n
Optional parameters:\n
[-c number of cores to use (default 1)]\n
[-d bedgraph file containing open chromatin signal, e.g. DNase1-seq, or Histone-Mark signal]\n
[-n column in the -b file containg the average per base signal within a peak. If this option is used, the -d option must not be used]\n
[-a gene annotation file, required to generate the gene view]\n
[-w size of the window to be considered to generate gene view (default 50kb)]\n
[-f annotate only DNase peaks that are within a window specified by the -w option around all genes contained in the gene annotation file specified by this option]\n
[-e indicating whether exponential decay should be used (default TRUE)]\n
[-l flag to be set if affinities should not be normalised by peak length]\n
[-u flag to be set if peak features for peak length and peak counts should not be generated]\n 
[-x if -d or -n is used together with this flag, the original (Decay-)Scaling formulation of TEPIC is used to compute gene-TF scores]\n
[-m path to a tab delimited file containing the length of the used PSEMs]\n
[-y flag to be set if the entire gene body should be screened for TF binding. The search window is extended by a region half of the size that is specified by the -w option upstream of the genes 5' TSS]\n
[-z flag indicating that the output of TEPIC should be zipped]\n
[-k path to a set of background sequences that should be used to compute to generate a binary score for TF binding. Mutually exclusive to the -r option]\n
[-r path to a 2bit representation of the reference genome, required to generate a binary score for TF binding. The binary score is generated in addition to the standard affinity values. Mutually exclusive to the -k option]\n
[-v p-value cut off used to determine a cut off to derive a binary score for TF binding (default 0.05)]\n
[-i minutes that should be spend at most per chromosome to find matching random regions (default 3)]\n
[-j flag indicating that the reference genome contains a chr prefix]\n
[-t flag indicating that the annotation should be transcript and not gene based]\n
[-h a loop list file containing chromatin contacts]\n
[-s size of the loop window used around a genes promoter to link chromatin loops to genes (default 5000)]\n
[-q parameter to be set if only peak features should be computed (default FALSE)]\n"

#Initialising parameters
genome=""
randomGenome=""
regions=""
prefixP=""
cores=1
pwms=""
dnase=""
column=""
annotation=""
window=50000
decay="TRUE"
sparsity="FALSE"
geneBody="FALSE"
lengthNorm="TRUE"
peakFeatures="TRUE"
motifLength=""
working_dir=$(cd $(dirname "$0") && pwd -P)
originalScaling="FALSE"
zip="FALSE"
pvalue="0.05"
minutes=3
chrPrefix="FALSE"
backgroundRegions=""
onlyPeakFeatures="FALSE"
transcripts="FALSE"
loopWindow=5000
loopList=""
#Parsing command line
while getopts "g:b:o:c:p:d:n:a:w:f:m:e:r:v:k:i:q:h:s:yluhxzjt" o;
do
case $o in
	g)	genome=$OPTARG;;
	b)	regions=$OPTARG;;
	o)	prefixP=$OPTARG;;
	c)	cores=$OPTARG;;
	p)	pwms=$OPTARG;;
	d)	dnase=$OPTARG;;
	n)	column=$OPTARG;;
	a)	annotation=$OPTARG;;
	w)	window=$OPTARG;;
	f)	filter=$OPTARG;;
	m)	motifLength=$OPTARG;;
	y)	geneBody="TRUE";;
	e)	decay=$OPTARG;;
	l)	lengthNorm="FALSE";;
	u)	peakFeatures="FALSE";;
	x)	originalScaling="TRUE";;
	r)	randomGenome=$OPTARG;;
	z)	zip="TRUE";;
	v)	pvalue=$OPTARG;;
	i)	minutes=$OPTARG;;
	j)	chrPrefix="TRUE";;
	t)	transcripts="TRUE";;
	k)	backgroundRegions=$OPTARG;;
	q)	onlyPeakFeatures=$OPTARG;;
	h)	loopList=$OPTARG;;
	s)	loopWindow=$OPTARG;;
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

if [ -n "$dnase" ] && [ -n "$column" ] ;
then
	echo The options -d and -n are mutually exclusive
	exit 1;
fi

if [ -n "$backgroundRegions" ] && [ -n "$randomGenome" ] ;
then
	echo The options -k and -r are mutually exclusive
	exit 1;
fi

if [ -n "$randomGenome" ] ;
then
	sparsity="TRUE"
fi

if [ -n "$backgroundRegions" ] ;
then
	sparsity="TRUE"
fi

d=$(date +%D)
d=`echo $d | sed 's/\//\_/g'`
t=$(date +%H:%M:%S:%N | sed 's/:/_/g')

prefix=$prefixP"_TEPIC_"${d}"_"${t}
filteredRegions=$prefix"_canidate_binding_regions"
#Generating name of the fasta file containing the overlapping regions
openRegionSequences=${prefix}.OpenChromatin.fasta
metadatafile=${prefix}.amd.tsv
#Create metadata file
touch $metadatafile
echo "[Description]" >> $metadatafile
echo "process	TEPICv2.1" >> $metadatafile
echo "run_by_user	"$USER >> $metadatafile
echo "date	"$d >> $metadatafile
echo "time	"$t >> $metadatafile
echo "analysis_id	"$prefix >> $metadatafile
echo "" >> $metadatafile
echo "[Inputs]" >> $metadatafile
echo "region_file	"$regions >> $metadatafile
if [ -n "$filter" ];
then
	echo "gene_filter_file "$filter >> $metadatafile
fi
if [ -n "$backgroundRegions" ];
then
	echo "background_regions "$backgroundRegions >> $metadatafile
fi
echo "" >> $metadatafile
echo "[References]" >> $metadatafile
echo "genome_reference	"$genome >> $metadatafile
if [ -n "$randomGenome" ];
then
echo "2bit_reference_genome	"$randomGenome >> $metadatafile
fi
echo "pwms	"$pwms >> $metadatafile
if [ -n "$annotation" ];
then
echo "genome_annotation	"$annotation>> $metadatafile
fi
if [ -n "$dnase" ];
then 
	echo "signale_file	"$dnase >> $metadatafile
fi
if [ -n "$column" ];
then 
	echo "signale_column	"$column >> $metadatafile
fi
echo "" >> $metadatafile
echo "[Outputs]" >> $metadatafile
echo "affinity_file_peak_view	"$prefix"_Affinity.txt" >> $metadatafile

decayout=""
if [ "$decay" == "TRUE" ];
then
	decayout="_Decay"
fi

if [ "$peakFeatures" == "TRUE" ];
then
	echo "affinity_gene_view_peak_features "${prefix}${decayout}"_Peak_Based_Features_Affinity_Gene_View_Filtered.txt" >> $metadatafile
else
	echo "affinity_file_gene_view_filtered	"${prefix}${decayout}"_Affinity_Gene_View_Filtered.txt" >> $metadatafile
fi

if [ -n "$dnase" ];
then 
	echo "signal_scaling_factors	"$prefix"_Peak_Coverage.txt" >> $metadatafile
	if [ "$originalScaling" == "FALSE" ];
	then
		echo "affinity_gene_view_peak_features_signal "${prefix}${decayout}"_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt" >> $metadatafile
	fi
	if [ "$originalScaling" == "TRUE" ];
	then
		echo "scaled affinity_peak_view	"$prefix"_Scaled_Affinity.txt" >> $metadatafile
		echo "scaled_affinity_gene_view_filtered	"${prefix}${decayout}"_Scaled_Affinity_Gene_View_Filtered.txt" >> $metadatafile
	fi
fi
if [ -n "$column" ];
then
	echo "signal_scaling_factors	"$prefix"_Peak_Coverage.txt" >> $metadatafile
	if [ "$originalScaling" == "FALSE" ];
	then
		echo "affinity_gene_view_peak_features_signal "${prefix}${decayout}"_Three_Peak_Based_Features_Affinity_Gene_View_Filtered.txt" >> $metadatafile
	fi
	if [ "$originalScaling" == "TRUE" ];
	then
		echo "scaled affinity_peak_view	"$prefix"_Scaled_Affinity.txt" >> $metadatafile
		echo "scaled_affinity_gene_view_filtered	"${prefix}${decayout}"_Scaled_Affinity_Gene_View_Filtered.txt" >> $metadatafile
	fi
fi
echo "" >> $metadatafile
echo "[Parameters]" >> $metadatafile
echo "original_scaling	"$originalScaling >> $metadatafile
echo "SampleID	"$prefixP >> $metadatafile
echo "cores	"$cores >> $metadatafile
echo "Chr prefix	"$chrPrefix >> $metadatafile
echo "gzip	"$zip >> $metadatafile
if [ -n "$annotation" ];
then
	echo "window	"$window >> $metadatafile
	echo "decay	"$decay >> $metadatafile
	echo "sparsity	"$sparsity >> $metadatafile
	echo "genebody	"$geneBody >> $metadatafile
	echo "length_normalisation	"$lengthNorm >> $metadatafile
	echo "peak_features	"$peakFeatures >> $metadatafile
	echo "transcript based annotation	"$transcripts >> $metadatafile
	if [ -n "$motifLength" ];
	then
		echo "motif_length	"$motifLength >> $metadatafile
	fi
	if [ -n "$loopList" ];
	then
		echo "loop_list	"$loopList >> $metadatafile
		echo "loop_window	"$loopWindow >> $metadatafile
	fi
fi
if [ -n "$randomGenome" ];
then
  echo "minutes	"$minutes >> $metadatafile
	echo "p-value	"$pvalue >> $metadatafile
	echo "sparse_affinity_gene_view	"${prefix}${decayout}"_Sparse_Affinity_Gene_View_Filtered.txt" >> $metadatafile
fi
echo "" >> $metadatafile
echo "[Metrics]" >> $metadatafile
numReg=`wc -l $regions | cut -f 1 -d " "`
echo "Number of provided regions	"$numReg >> $metadatafile
numMat=`grep ">" $pwms | wc -l`
echo "Number of considered pwms	"$numMat >> $metadatafile 


echo "Preprocessing region file: Removing chr prefix, sorting regions and removing duplicats"
sed 's/chr//g' $regions >  ${filteredRegions}_Filtered_Regions.bed
sort -s -k1,1 -k2,2 -k3,3 ${filteredRegions}_Filtered_Regions.bed | uniq > ${filteredRegions}_sorted.bed
rm ${filteredRegions}_Filtered_Regions.bed

if [ -n "$backgroundRegions" ];
then
echo "Preprocessing background file"
sed 's/chr//g' $backgroundRegions >  ${filteredRegions}_Filtered_Regions_background.bed
sort -s -k1,1 -k2,2 -k3,3 ${filteredRegions}_Filtered_Regions_background.bed | uniq > ${prefix}_Random_Regions.bed
rm ${filteredRegions}_Filtered_Regions_background.bed
fi 

if [ -n "$filter" ];
then
echo "Filter total peak set"
python ${working_dir}/generateIntersectionWindows.py ${filter} ${window} ${geneBody} ${transcripts} > ${prefix}_gene_windows.temp.bed 
bedtools intersect -b ${prefix}_gene_windows.temp.bed -a ${filteredRegions}_sorted.bed -u > ${filteredRegions}_temp.bed
mv ${filteredRegions}_temp.bed ${filteredRegions}_sorted.bed
rm ${prefix}_gene_windows.temp.bed 
fi

#Generating random genomic regions and computing TF affinities
if [ -n "$randomGenome" ];
then
echo "Generating random genomic regions"
python ${working_dir}/findBackground.py -i ${filteredRegions}_sorted.bed -g ${randomGenome} -o ${prefix}_Random_Regions.bed -w ${cores} --time-out ${minutes}
fi

if [ ${chrPrefix} == "TRUE" ];
then
	echo "Adapting chr prefix in bed files for intersection with the reference genome"
	awk '{print "chr"$1"\t"$2"\t"$3}' ${filteredRegions}_sorted.bed > ${prefix}_tempFasta.bed
	getFastaRegion=${prefix}_tempFasta.bed
	if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ]; 
		then
		awk '{print "chr"$1"\t"$2"\t"$3}' ${prefix}_Random_Regions.bed > ${prefix}_tempFastaRandom.bed
		getFastaRegionRandom=${prefix}_tempFastaRandom.bed
		fi
else
	getFastaRegion=${filteredRegions}_sorted.bed
	getFastaRegionRandom=${prefix}_Random_Regions.bed
fi

echo "Runnig bedtools"
#Run bedtools to get a fasta file containing the sequence data for predicted open chromatin regions contained in the bedfile
bedtools getfasta -fi $genome -bed ${getFastaRegion} -fo $openRegionSequences
if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
then
bedtools getfasta -fi $genome -bed ${getFastaRegionRandom} -fo ${openRegionSequences}_Random
rm ${prefix}_Random_Regions.bed
fi
if [ ${chrPrefix} == "TRUE" ];
then
	rm ${prefix}_tempFasta.bed
	if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
	then
		rm ${prefix}_tempFastaRandom.bed	
	fi
fi

echo "Converting invalid characters"
#Remove R and Y from the sequence
python ${working_dir}/convertInvalidCharacterstoN.py $openRegionSequences $prefix-FilteredSequences.fa
rm $openRegionSequences

if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
then
python ${working_dir}/convertInvalidCharacterstoN.py ${openRegionSequences}_Random $prefix-Random-FilteredSequences.fa
rm ${openRegionSequences}_Random
fi

#Use TRAP to compute transcription factor affinities to the above extracted sequences
affinity=${prefix}_Affinity.txt

echo "Starting TRAP"
${working_dir}/TRAPmulti $pwms ${prefix}-FilteredSequences.fa $cores > ${affinity}_temp 
rm ${prefix}-FilteredSequences.fa
if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
then
${working_dir}/TRAPmulti $pwms ${prefix}-Random-FilteredSequences.fa $cores > ${affinity}_Random 
rm ${prefix}-Random-FilteredSequences.fa
fi

if [ ${chrPrefix} == "TRUE" ];
then
	sed -i 's/chr//g' ${affinity}_temp
	if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
	then
		sed -i 's/chr//g' ${affinity}_Random
	fi
fi



if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
then
echo "Discretising TF affinities"
Rscript ${working_dir}/IdentifyCut-Offs.R ${affinity}_Random ${affinity}_temp ${prefix}_Filtered_Affinities_temp.txt ${pvalue}
#rm ${affinity}_Random
fi

#Computing DNase Coverage in Peak regions
if [ -n "$dnase" ] ;
then 
	sort -s -k1,1 -k2,2 -k3,3 $dnase > ${dnase}_sorted
	python ${working_dir}/computeDNaseCoverage.py ${dnase}_sorted ${filteredRegions}_sorted.bed > ${prefix}_Peak_Coverage.txt
	rm ${dnase}_sorted
	if [ "$originalScaling" == "TRUE" ];
	then
		python ${working_dir}/scaleAffinity.py --is-sorted -s ${prefix}_Peak_Coverage.txt -a ${affinity}_temp > ${prefix}_Scaled_Affinity_temp.txt
		if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
		then
			python ${working_dir}/scaleAffinity.py --is-sorted -s ${prefix}_Peak_Coverage.txt -a ${prefix}_Filtered_Affinities_temp.txt > ${prefix}_Filtered_Scaled_Affinity_temp.txt
		fi
	fi
fi

if [ -n "${column}" ] ;
then
	if [ "$originalScaling" == "TRUE" ] ;
	then
		python ${working_dir}/scaleAffinity.py --is-sorted --scale-col ${column} -s ${filteredRegions}_sorted.bed -a ${affinity}_temp > ${prefix}_Scaled_Affinity_temp.txt
		if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
		then
			python ${working_dir}/scaleAffinity.py --is-sorted --scale-col ${column} -s ${filteredRegions}_sorted.bed -a ${prefix}_Filtered_Affinities_temp.txt > ${prefix}_Filtered_Scaled_Affinity_temp.txt
		fi
	else
		cut -f 1,2,3,${column} ${filteredRegions}_sorted.bed > ${prefix}_Peak_Coverage.txt
	fi
fi	
rm ${filteredRegions}_sorted.bed

#Removing regions that could not be annotated
echo "Filter regions that could not be annotated"
python ${working_dir}/filterInvalidRegions.py ${affinity}_temp $affinity
if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ] ;
then
	python ${working_dir}/filterInvalidRegions.py ${prefix}_Filtered_Affinities_temp.txt ${prefix}_Thresholded_Affinities.txt
	rm ${prefix}_Filtered_Affinities_temp.txt
fi
rm ${affinity}_temp

if [ -n "$dnase" ] ||  [ -n "$column" ] ;
then
	if [ "$originalScaling" == "TRUE" ] ;
	then
		python ${working_dir}/filterInvalidRegions.py ${prefix}_Scaled_Affinity_temp.txt ${prefix}_Scaled_Affinity.txt
		rm ${prefix}_Scaled_Affinity_temp.txt
		if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
		then
			python ${working_dir}/filterInvalidRegions.py ${prefix}_Filtered_Scaled_Affinity_temp.txt ${prefix}_Thresholded_Scaled_Affinity.txt
			rm ${prefix}_Filtered_Scaled_Affinity_temp.txt
		fi
	fi
fi


#If an annotation file is provied, the gene view is generated
if [ -n "$annotation" ]; 
then
	echo "Generating gene scores"
	if [ -n $loopList ]
	then
		if [ -n "$dnase" ] ||  [ -n "$column" ] ;
		then
			if [ "$originalScaling" == "FALSE" ] ;
			then
				python ${working_dir}/annotateTSS.py ${annotation} ${affinity} "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--peakCoverage" ${prefix}_Peak_Coverage.txt "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--additionalPeakFeatures" ${peakFeatures}  "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts} "--lwindows" ${loopWindow} "--conformationData" ${loopList}
			else
				python ${working_dir}/annotateTSS.py ${annotation} ${affinity} "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay  "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--signalScale" ${prefix}_Scaled_Affinity.txt "--additionalPeakFeatures" ${peakFeatures}  "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts} "--lwindows" ${loopWindow} "--conformationData" ${loopList}
			fi
		else
			python ${working_dir}/annotateTSS.py ${annotation} ${affinity} "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--geneBody" $geneBody "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--additionalPeakFeatures" ${peakFeatures}  "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts} "--lwindows" ${loopWindow} "--conformationData" ${loopList}
		fi 
	else
		if [ -n "$dnase" ] ||  [ -n "$column" ] ;
		then
			if [ "$originalScaling" == "FALSE" ] ;
			then
				python ${working_dir}/annotateTSS.py ${annotation} ${affinity} "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--peakCoverage" ${prefix}_Peak_Coverage.txt "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--additionalPeakFeatures" ${peakFeatures}  "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts}
			else
				python ${working_dir}/annotateTSS.py ${annotation} ${affinity} "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay  "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--signalScale" ${prefix}_Scaled_Affinity.txt "--additionalPeakFeatures" ${peakFeatures}  "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts}
			fi
		else
			python ${working_dir}/annotateTSS.py ${annotation} ${affinity} "--geneViewAffinity" ${prefix}_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--geneBody" $geneBody "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--additionalPeakFeatures" ${peakFeatures}  "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts}
		fi
	fi 

	if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
	then
		if [ -n "$dnase" ] ||  [ -n "$column" ] ;
			then
			if [ "$originalScaling" == "FALSE" ] ;
				then
				python ${working_dir}/annotateTSS.py ${annotation} ${prefix}_Thresholded_Affinities.txt "--geneViewAffinity" ${prefix}_Thresholded_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--peakCoverage" ${prefix}_Peak_Coverage.txt "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--additionalPeakFeatures" ${peakFeatures} "--sparseRep" $sparsity  "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts}
				else
				python ${working_dir}/annotateTSS.py ${annotation} ${prefix}_Thresholded_Affinities.txt "--geneViewAffinity" ${prefix}_Thresholded_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--sparseRep" $sparsity "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--signalScale" ${prefix}_Thresholded_Scaled_Affinity.txt "--additionalPeakFeatures" ${peakFeatures} "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts}
			fi
		else
		python ${working_dir}/annotateTSS.py ${annotation} ${prefix}_Thresholded_Affinities.txt "--geneViewAffinity" ${prefix}_Thresholded_Affinity_Gene_View.txt "--windows" $window "--decay" $decay "--geneBody" $geneBody "--geneBody" ${geneBody} "--normaliseLength" ${lengthNorm} "--motifLength" ${motifLength} "--additionalPeakFeatures" ${peakFeatures} "--sparseRep" $sparsity "--onlyPeakFeatures" ${onlyPeakFeatures} "--transcript" ${transcripts}
		fi
	fi
	#Creating files containing only genes for which TF predictions are available
	echo "Filter genes that could not be annotated"
	if [ "$decay" == "TRUE" ];
	then
		if [ -n "$dnase" ] || [ -n "$column" ];
		then
			if [ "$peakFeatures"  == "TRUE" ];
			then
				if [ "$originalScaling" == "FALSE" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Decay_Three_Peak_Based_Features_Affinity_Gene_View.txt
					rm ${prefix}_Decay_Three_Peak_Based_Features_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Decay_Three_Peak_Based_Features_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Decay_Three_Peak_Based_Features_Affinity_Gene_View.txt
					fi
				else
					python ${working_dir}/filterGeneView.py ${prefix}_Decay_Scaled_Peak_Features_Affinity_Gene_View.txt
					rm ${prefix}_Decay_Scaled_Peak_Features_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Decay_Scaled_Peak_Features_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Decay_Scaled_Peak_Features_Affinity_Gene_View.txt
					fi
				fi
			else
				if [ "$originalScaling" == "TRUE" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Decay_Scaled_Affinity_Gene_View.txt
					rm ${prefix}_Decay_Scaled_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Decay_Scaled_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Decay_Scaled_Affinity_Gene_View.txt
					fi
				else
					python ${working_dir}/filterGeneView.py ${prefix}_Decay_Signal_Feature_Affinity_Gene_View.txt
					rm ${prefix}_Decay_Signal_Feature_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Decay_Signal_Feature_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Decay_Signal_Feature_Affinity_Gene_View.txt
					fi
				fi
			fi
		fi
			if [ "$peakFeatures" == "TRUE" ];
			then
				python ${working_dir}/filterGeneView.py ${prefix}_Decay_Peak_Features_Affinity_Gene_View.txt
				rm ${prefix}_Decay_Peak_Features_Affinity_Gene_View.txt
				if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Decay_Peak_Features_Affinity_Gene_View.txt
					rm ${prefix}_Thresholded_Decay_Peak_Features_Affinity_Gene_View.txt
				fi
			else
				python ${working_dir}/filterGeneView.py ${prefix}_Decay_Affinity_Gene_View.txt
				rm ${prefix}_Decay_Affinity_Gene_View.txt
				if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Decay_Affinity_Gene_View.txt
					rm ${prefix}_Thresholded_Decay_Affinity_Gene_View.txt
				fi
			fi		
	else
		if [ -n "$dnase" ] || [ -n "$column" ];
		then
			if [ "$peakFeatures" == "TRUE" ];
			then
				if [ "$originalScaling" == "FALSE" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Three_Peak_Based_Features_Affinity_Gene_View.txt
					rm ${prefix}_Three_Peak_Based_Features_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Three_Peak_Based_Features_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Three_Peak_Based_Features_Affinity_Gene_View.txt
					fi
				else
					python ${working_dir}/filterGeneView.py ${prefix}_Scaled_Peak_Features_Affinity_Gene_View.txt
					rm ${prefix}_Scaled_Peak_Features_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Scaled_Peak_Features_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Scaled_Peak_Features_Affinity_Gene_View.txt
					fi
				fi
			else
				if [ "$originalScaling" == "TRUE" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Scaled_Affinity_Gene_View.txt
					rm ${prefix}_Scaled_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Scaled_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Scaled_Affinity_Gene_View.txt
					fi
				else
					python ${working_dir}/filterGeneView.py ${prefix}_Signal_Feature_Affinity_Gene_View.txt
					rm ${prefix}_Signal_Feature_Affinity_Gene_View.txt
					if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
					then
						python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Signal_Feature_Affinity_Gene_View.txt
						rm ${prefix}_Thresholded_Signal_Feature_Affinity_Gene_View.txt
					fi
				fi
			fi
		fi
			if [ "$peakFeatures" == "TRUE" ]
			then
				python ${working_dir}/filterGeneView.py ${prefix}_Peak_Features_Affinity_Gene_View.txt
				rm ${prefix}_Peak_Features_Affinity_Gene_View.txt
				if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Peak_Features_Affinity_Gene_View.txt
					rm ${prefix}_Thresholded_Peak_Features_Affinity_Gene_View.txt
				fi
			else
				python ${working_dir}/filterGeneView.py ${prefix}_Affinity_Gene_View.txt
				rm ${prefix}_Affinity_Gene_View.txt
				if [ -n "$randomGenome" ] || [ -n "$backgroundRegions" ];
				then
					python ${working_dir}/filterGeneView.py ${prefix}_Thresholded_Affinity_Gene_View.txt
					rm ${prefix}_Thresholded_Affinity_Gene_View.txt
				fi
			fi
	fi
	if [ "$zip" == "TRUE" ];
	then
		gzip ${prefix}*Affinity_Gene_View_Filtered.txt
		if [ "$sparsity" == "TRUE" ];
		then
			gzip ${prefix}*Affinity_Gene_View.txt
			gzip ${prefix}_Thresholded_Affinities.txt
		fi
		gzip ${affinity}
		if [ -n "$dnase" ] || [ -n "$column" ];
		then
			if [ "$originalScaling" == "TRUE" ];
			then
				gzip ${prefix}_Scaled_Affinity.txt
				if [ -z "$column" ];
				then
					gzip ${prefix}_Peak_Coverage.txt
				fi
			else
				gzip ${prefix}_Peak_Coverage.txt
			fi
		fi
	fi
fi

