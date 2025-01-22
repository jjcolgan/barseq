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
cat samples.txt | parallel -j 4 '
  pigz -d /project/blekhman/jjcolgan/barseq/00_RAW/20250107_LH00315_0348_A225WL3LT1-EC-CG-2s-pls12/FastX/{}_R1_001.fastq.gz;
  pigz -d /project/blekhman/jjcolgan/barseq/00_RAW/20250107_LH00315_0348_A225WL3LT1-EC-CG-2s-pls12/FastX/{}_R1_002.fastq.gz;
  cat /project/blekhman/jjcolgan/barseq/00_RAW/20250107_LH00315_0348_A225WL3LT1-EC-CG-2s-pls12/FastX/{}*gz > g/paired/{}BothReads.fastq
  perl MultiCodes.pl -out g/paired/{}BothReads \
  -index {} \
  < g/full/{}_R1_001.fastq \
  > g/full/{}.out 2> g/full/{}.err;
  pigz g/full/{}_R1_001.fastq
  pigz g/full/{}_R1_002.fastq
  pigz g/paired/{}BothReads.fastq
'

list=$(
for file in $(find g/paired/ -name '*.codes'); do
	echo $file
done
)

perl combineBarSeq.pl g/paired/combine pool ${list}

perl BarSeqR.pl -org paired