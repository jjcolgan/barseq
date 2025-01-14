#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=16G
#SBATCH --err=resolveBarcodePosition.err
#SBATCH --out=resolveBarcodePosition.out
#SBATCH --job-name=barcodeinclusion
#SBATCH --partition=blekhman

matlab -nodesktop -nosplash -r \
"addpath('/home/jjcolgan/resolveBarcodePostion/scripts'); \
resolve_barcode_position('/project/blekhman/jjcolgan/barseq/bin/g/lane1/MultiCounts/inclusionTable', \
'/project/blekhman/jjcolgan/barseq/bin/g/bBreve/lookuptable', \
'/project/blekhman/jjcolgan/barseq/bin/g/bBreve/genes.tab', 'positiontable'); \
exit;"
