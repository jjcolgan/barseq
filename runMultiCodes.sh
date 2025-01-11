#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=12G
#SBATCH --err=multiCodes.err
#SBATCH --out=multiCodes.err
#SBATCH --job-name=multicodes
#SBATCH --partition=caslake
#Change lane to which ever run is being analyzed
mkdir g/lane1/MultiCounts/
mkdir g/lane1/MultiCodeOut
for sample in `cat samples.txt`
do
  pigz -d g/lane1/${sample}_L001_R1_001.fastq.gz
  perl MultiCodes.pl -out g/lane1/MultiCounts/${sample} \
  -index ${sample} \
  < bin/g/lane1/${sample}_L001_R1_001.fastq
  pigz g/lane1/${sample}_L001_R1_001.fastq
done
