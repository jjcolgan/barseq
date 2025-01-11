#adapted from Megan's barseq_analysis_share.sh
mkdir /project/blekhman/jjcolgan/barseq/bin/g/lane1/MultiCounts/
perl MultiCodes.pl -out /project/blekhman/jjcolgan/barseq/bin/g/lane1/MultiCounts/SAMPLE -index INDEX < \
/project/blekhman/jjcolgan/barseq/bin/g/lane1/EC-CG-2s-pl1-2481-Co_S18_L001_R1_001.fastq.gz
