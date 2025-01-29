library(tidyverse)
library(ComplexHeatmap)

fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
metadata = read_tsv('fullbarseqMeta.txt')

colnames(fitnessScores) <- sub("setA", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("_.*", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("CO$", "Co", colnames(fitnessScores))
colnames(fitnessScores) <- sub("DJ$", "Dj", colnames(fitnessScores))


fitnessScoresLong=fitnessScores %>%
  select(-c(desc,
            sysName))%>%
  pivot_longer(cols = c(2:43), names_to = 'sample', values_to = 'logRatios')%>%
  left_join(metadata, by = 'sample')

fitnessScoresLong %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios')%>%
  column_to_rownames('locusId')%>%
  Heatmap(show_column_names = T, show_row_names = F)

fitnessScoresLong %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios')%>%
  column_to_rownames('locusId')%>%
  scale()%>%
  Heatmap(show_column_names = T, show_row_names = F)