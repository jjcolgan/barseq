'Wrapper script for writing bash code to pool barseq days'

library(tidyverse)

metadata = read_tsv('fullbarseqMeta.txt')

tissueDay = metadata%>%
  select(tissueDay)%>%
  distinct()%>%
  .$tissueDay

for (t in tissueDay){
  samples=metadata %>%
    filter(tissueDay == t)%>%
    .$sample
  for (s in samples){
    print(paste0('cat ',s, '*.fastq.gz >> ', t, '.fastq.gz'))
  }
}
