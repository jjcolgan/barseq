#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=16G
#SBATCH --err=barcodeinclusion.err
#SBATCH --out=barcodeinclusion.out
#SBATCH --job-name=barcodeinclusion
#SBATCH --partition=blekhman

python3 /home/jjcolgan/bethaBarseqRepo/barcodeinclusion.py -k /project/blekhman/jjcolgan/barseq/bin/key.txt \
--sequencemodel 5:20 CAGCGTACG AGAGACCTCG \
-o /project/blekhman/jjcolgan/barseq/bin/g/lane1/MultiCounts/inclusionTable
