#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=16G
#SBATCH --err=designRandomPool.err
#SBATCH --out=designRandomPool.out
#SBATCH --job-name=designSortedLibrary
#SBATCH --partition=blekhman

perl DesignRandomPool.pl -pool pool -genes genes.tab \
'F1_tnpsites.txt' \
'F2_tnpsites.txt' \
'F3_tnpsites.txt' \
'F4_tnpsites.txt' \
'F5_tnpsites.txt'