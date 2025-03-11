library(tidyverse)
library(Maaslin2)
fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$day = factor(metadata$day, levels = c('t0', 'day1', 'day3', 'day7', 'day14'))

colnames(fitnessScores) <- sub("setA", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("_.*", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("CO$", "Co", colnames(fitnessScores))
colnames(fitnessScores) <- sub("DJ$", "Dj", colnames(fitnessScores))

colonMeta = metadata%>%
  filter(tissue == 'colon')%>%
  column_to_rownames('sample')

djMeta = metadata %>%
  filter(tissue == 'dj')%>%
  column_to_rownames('sample')

maaslinIn=fitnessScores %>%
  select(-c(sysName, desc))%>%
  column_to_rownames('locusId')
# variance filter was primarily added to only focus on samples where there
# actually a dramatic change in abundance and reduce the number of tests
# performed
Maaslin2(input = maaslinIn,
         input_metadata = colonMeta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = 1,
         fixed_effects = c('day','pfClusters', 'millionBases'),
         random_effects = ,
         reference = 'day,day1;day3;day7;day14',
         output = 'colonLogRatiosMaaslin2'



)
sigResMaaslin=read_tsv('colonLogRatiosMaaslin2/significant_results.tsv')

sigResMaaslin %>%
  filter(metadata == 'day')%>%
  filter(qval < .05)%>%
  select(feature)%>%
  distinct()%>%
  nrow()

Maaslin2(input = maaslinIn,
         input_metadata = djMeta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = .5,
         fixed_effects = c('day','pfClusters', 'millionBases'),
         random_effects = ,
         output = 'djLogRatiosMaaslin2'



)