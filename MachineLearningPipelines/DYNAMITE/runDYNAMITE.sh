#! /bin/bash
#Loading arguments from configFile differentialAnalysis.cfg
echo Loading arguments from configFile differentialAnalysis.cfg
source DYNAMITE.cfg

#This part runs the necessary scripts to perform a differential analysis using TEPIC affinity ratios
#DO NOT PERFORM CHANGES HERE!
if [ -n "$coverageFile" ] && [ -n "$coverageColumn" ] ;
then
     echo The parameters coverageFile and coverageColumn are mutualy exclusive
     exit 1;
fi

mkdir $outputDirectory
mkdir $outputDirectory"/Affinities"

if [ -n "$existing_TEPIC_Results_Group1" ];
then
	echo "Using already computed TF affinities"
else
mkdir $outputDirectory"/Affinities/group1"
mkdir $outputDirectory"/Affinities/group2"

echo "Running TEPIC for all files of group 1"
counter=0
for file in $open_regions_Group1
do
	prefix=$(basename $file)
	if [ "$peakFeatures" == "TRUE" ];
	then
		if [ -z "$coverage_Files_Group1" ] && [ -z "$coverage_Column_Group1" ] ;
		then
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay
		fi
		if [ -n "$coverage_Files_Group1" ] ;
		then
			cfile=${coverage_Files_Group1[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cfile  -f $geneAnnotation -e $decay
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cfile  -f $geneAnnotation -e $decay
			((counter++))
		fi
		if [ -n "$coverage_Column_Group1" ] ;
		then
			cCol=${coverage_Column_Group1[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay
			((counter++))
		fi
	else
		if [ -z "$coverage_Files_Group1" ] && [ -z "$coverage_Column_Group1" ] ;
		then
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay -u
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay -u
		fi
		if [ -n "$coverage_Files_Group1" ] ;
		then
			cfile=${coverage_Files_Group1[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cfile -f $geneAnnotation -e $decay -u
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cfile  -f $geneAnnotation -e $decay -u
			((counter++))
		fi
		if [ -n "$coverage_Column_Group1" ] ;
		then
			cCol=${coverage_Column_Group1[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay -u
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group1/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay -u
			((counter++))
		fi
	fi
done

echo "Running TEPIC for all files of group 2"
counter=0
for file in $open_regions_Group2
do
	prefix=$(basename $file)
	if [ "$peakFeatures" == "TRUE" ];
	then
		if [ -z "$coverage_Files_Group2" ] && [ -z "$coverage_Column_Group2" ] ;
		then
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay
		fi
		if [ -n "$coverage_Files_Group2" ] ;
		then
			cFile=${coverage_Files_Group2[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cFile -f $geneAnnotation -e $decay
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cFile -f $geneAnnotation -e $decay
			((counter++))
		fi
		if [ -n "$coverage_Column_Group2" ] ;
		then
			cCol=${coverage_Column_Group2[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay
			((counter++))
		fi
	else
		if [ -z "$coverage_Files_Group2" ] && [ -z "$coverage_Column_Group2" ] ;
		then
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay -u
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -f $geneAnnotation -e $decay -u
		fi
		if [ -n "$coverage_Files_Group2" ] ;
		then
			cFile=${coverage_Files_Group2[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cFile -f $geneAnnotation -e $decay -u
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -d $cFile -f $geneAnnotation -e $decay -u
			((counter++))
		fi
		if [ -n "$coverage_Column_Group2" ] ;
		then
			cCol=${coverage_Column_Group2[${counter}]}
			echo bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay -u
			bash ${path}/Code/TEPIC.sh -g $referenceGenome -b $file -o $outputDirectory/Affinities/group2/$prefix -p $pwm -c $cores_TEPIC -a $geneAnnotation -w $window -n $cCol -f $geneAnnotation -e $decay -u
			((counter++))
		fi
	fi
done
fi
echo "Postprocessing TF affinities"
mkdir $outputDirectory"/Affinities/Mean/"
mkdir $outputDirectory"/Affinities/Ratio/"
if [ -z $existing_TEPIC_Results_Group1 ];
then
	if [ -z "$coverage_Files_Group1" ] && [ -z "$coverage_Column_Group1" ] && [ -z "$coverage_Files_Group2" ] && [ -z "$coverage_Column_Group2" ] ; 
	then
		echo python ${scriptPath}/computeMeanRatioTFAffinities.py $outputDirectory/Affinities/group1/ $outputDirectory/Affinities/group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt"  "False" $peakFeatures
		python ${scriptPath}/computeMeanRatioTFAffinities.py $outputDirectory/Affinities/group1/ $outputDirectory/Affinities/group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" "False" $peakFeatures
	else
		echo python ${scriptPath}/computeMeanRatioTFAffinities.py $outputDirectory/Affinities/group1/ $outputDirectory/Affinities/group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" "True" $peakFeatures
		python ${scriptPath}/computeMeanRatioTFAffinities.py $outputDirectory/Affinities/group1/ $outputDirectory/Affinities/group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" "True" $peakFeatures
	fi
	else
	if [ -z "$coverage_Files_Group1" ] && [ -z "$coverage_Column_Group1" ] && [ -z "$coverage_Files_Group2" ] && [ -z "$coverage_Column_Group2" ] ;
	then
		echo python ${scriptPath}/computeMeanRatioTFAffinities.py $existing_TEPIC_Results_Group1/ $existing_TEPIC_Results_Group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" "False" $peakFeatures
		python ${scriptPath}/computeMeanRatioTFAffinities.py $existing_TEPIC_Results_Group1/ $existing_TEPIC_Results_Group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" "False" $peakFeatures
	else
		echo python ${scriptPath}/computeMeanRatioTFAffinities.py $existing_TEPIC_Results_Group1/ $existing_TEPIC_Results_Group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" "True" $peakFeatures
		python ${scriptPath}/computeMeanRatioTFAffinities.py $existing_TEPIC_Results_Group1/ $existing_TEPIC_Results_Group2/ $outputDirectory"/Affinities/Mean/Mean_Affinities_group1.txt" $outputDirectory"/Affinities/Mean/Mean_Affinities_group2.txt" $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" "True" $peakFeatures
	fi
fi

echo "Combining TF affinities with differential gene expression data"
mkdir $outputDirectory"/IntegratedData"
mkdir $outputDirectory"/IntegratedData/Log2"
echo python ${scriptPath}/integrateData.py $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" $differential_Gene_Expression_Data $outputDirectory"/IntegratedData/Log2/Integrated_Data_Log2_Quotient.txt"
python ${scriptPath}/integrateData.py $outputDirectory"/Affinities/Ratio/Ratio_Affinities_group1_vs_group2.txt" $differential_Gene_Expression_Data $outputDirectory"/IntegratedData/Log2/Integrated_Data_Log2_Quotient.txt"

mkdir $outputDirectory"/IntegratedData/Binary"
echo Rscript ${scriptPath}/prepareForClassificiation.R  $outputDirectory"/IntegratedData/Log2/Integrated_Data_Log2_Quotient.txt" $outputDirectory"/IntegratedData/Binary/Integrated_Data_For_Classification.txt"
Rscript ${scriptPath}/prepareForClassificiation.R  $outputDirectory"/IntegratedData/Log2/Integrated_Data_Log2_Quotient.txt" $outputDirectory"/IntegratedData/Binary/Integrated_Data_For_Classification.txt"

echo "Starting differential learning. Please wait..."
echo Rscript ${scriptPath}/DYNAMITE.R --dataDir=$outputDirectory"/IntegratedData/Binary/" --outDir=$outputDirectory"/Learning_Results/" --out_var="Expression" --Ofolds=$outerCV --Ifolds=$innerCV --performance=${performance} --alpha=$alpha_Step_Size --cores=$coresR
Rscript ${scriptPath}/DYNAMITE.R --dataDir=$outputDirectory"/IntegratedData/Binary/" --outDir=$outputDirectory"/Learning_Results/" --out_var="Expression" --Ofolds=${outerCV} --Ifolds=${innerCV} --performance=${performance} --alpha=${alpha_Step_Size} --cores=${coresR}

if [ "$randomise" == "TRUE" ];
then
mkdir $outputDirectory"/Learning_Results_Random/"
echo Rscript ${scriptPath}/DYNAMITE.R --dataDir=$outputDirectory"/IntegratedData/Binary/" --outDir=$outputDirectory"/Learning_Results_Random/" --out_var="Expression" --Ofolds=$outerCV --Ifolds=$innerCV --performance=${performance} --alpha=$alpha_Step_Size --cores=$coresR --randomise=$randomise
Rscript ${scriptPath}/DYNAMITE.R --dataDir=$outputDirectory"/IntegratedData/Binary/" --outDir=$outputDirectory"/Learning_Results_Random/" --out_var="Expression" --Ofolds=${outerCV} --Ifolds=${innerCV} --performance=${performance} --alpha=${alpha_Step_Size} --cores=${coresR} --randomise=${randomise}
fi
