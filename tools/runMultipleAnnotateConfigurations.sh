#!/bin/bash
# Runs the annotateTSS script with all combinations of defined parameters spawning at most processed as configured. Waits for all processed to be finished.

############################ CONFIG PART - Set your parameters here
maxprocessesinparallel=5

annotationfile="../Example/example_annotation.gtf"
affinityfile="../Example/example_affinities.txt"
loopfile="../Example/example_loopfile.txt"

geneViewAffinity=''

decay=("True")
hicresolutions=(5000 10000 25000) 	# optional: add "all"
tsswindows=(1000 2000 3000)
loopwindows=(5000 10000 25000 50000)
loopdecay=("True" "False")
									#TODO: also consider usemiddle and signaleScale options

############################ END CONFIG PART

runningprocesses=()
currentlyrunning=0
id=0

echo ${#decay[@]}
echo ${#hicresolutions[@]}
echo ${#tsswindows[@]}
echo ${#loopwindows[@]}
echo ${#loopdecay[@]}
return 0

for dec in "${decay[@]}":
	do
		for res in "${hicresolutions[@]}"
			do
				for twindow in "${tsswindows[@]}"
					do
						for lwindow in "${loopwindows[@]}"
							do
								for ldec in "${loopdecay[@]}"
									do	
										# build new combination and run it
										echo 'Spawned new annotateTSS instance'
										python ../Code/annotateTSS.py ${annotationfile} ${affinityfile}  "--geneViewAffinity" ${geneViewAffinity}_${id}_Affinity_Gene_View.txt "--windows" $twindow "--decay" $dec "--loopfile" $loopfile "--loopwindows" $lwindow "--resolution" $res "--loopdecay" $ldec &
										runningprocesses[currentlyRunning]=$!
										echo $!
										((id++))
										((currentlyrunning++))
										if [ "$currentlyrunning" -ge "$maxprocessesinparallel" ]
											then
												# wait for all processes until they terminated
												echo 'Waiting for spawned instances to finish...'
												wait ${runningprocesses[*]}
												currentlyrunning=0
												runningprocesses=()
												# sleep for some seconds to ensure a clear output buffer
												sleep 2
										fi
								done
						done
				done
		done
done

if [ ${#runningprocesses[@]} -ne 0 ]
	then
		# wait for remaining processes until they terminated
		wait ${runningprocesses[*]}
		currentlyrunning=0
		runningprocesses=()
fi

echo 'Finished all runs!'
