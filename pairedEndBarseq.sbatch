#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=32
#SBATCH --mem=48G
#SBATCH --err=pairedbarseq.err
#SBATCH --out=pairedbarseq.out
#SBATCH --job-name=fullbarseq
#SBATCH --partition=blekhman

module load R
module load parallel
mkdir g/paired
cat samples.txt | parallel -j 4 '
  pigz -d g/full/{}_R1_001.fastq.gz;
  pigz -d g/full/{}_R1_002.fastq.gz;
  cat g/paired/{}*gz > g/paired/{}BothReads.fastq
  perl MultiCodes.pl -out g/paired/{}BothReads \
  -index {} \
  < g/full/{}_R1_001.fastq \
  > g/full/{}.out 2> g/full/{}.err;
  pigz g/full/{}_R1_001.fastq
  pigz -d g/full/{}_R1_002.fastq
  g/paired/{}BothReads.fastq
'

list=$(
for file in $(find g/paired/ -name '*.codes'); do
	echo $file
done
)
perl combineBarSeq.pl g/paired/combine pool ${list}

perl BarSeqR.pl -org paired