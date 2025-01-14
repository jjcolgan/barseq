#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=16G
#SBATCH --err=designSortedLibrary.err
#SBATCH --out=designSortedLibrary.out
#SBATCH --job-name=designSortedLibrary
#SBATCH --partition=blekhman

python3 $RBP/DesignSortedLibrary.py3 --input /project/blekhman/jjcolgan/barseq/bin/g/lane1/mappedInserts/tncounts \
--genes /project/blekhman/jjcolgan/barseq/bin/g/bBreve/genes.tab \
--Nmin 10 --output /project/blekhman/jjcolgan/barseq/bin/g/bBreve/lookuptable