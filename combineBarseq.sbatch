#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=12G
#SBATCH --err=combineBarseq.sbatch.err
#SBATCH --out=combineBarseq.sbatch.out
#SBATCH --job-name=combineBarseq.sbatch
#SBATCH --partition=caslake
mkdir combine
list=$(
for file in $(find /project/blekhman/jjcolgan/barseq/bin/g/lane1/MultiCounts -name '*.codes'); do
	echo $file
done
)
perl combineBarSeq.pl /project/blekhman/jjcolgan/barseq/bin/g/combine/ pool ${list}