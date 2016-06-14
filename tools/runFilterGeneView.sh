#!/bin/bash

prefix="/Path/To/Your/Folder/"

for i in $(seq 0 71)
do
	python ../Code/filterGeneView.py $prefix"_"$i"_Scaled_Decay_Affinity_Gene_View.txt"
done
