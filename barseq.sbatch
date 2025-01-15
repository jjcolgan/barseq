#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=16G
#SBATCH --err=barseq.err
#SBATCH --out=basrseq.out
#SBATCH --job-name=barseq
#SBATCH --partition=caslake

perl BarseqR.perl -org bBreve