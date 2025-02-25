library(tidyverse)
library(ComplexHeatmap)

tScores = read_tsv('barseqAdjustedParams/fit_t.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$day = factor(metadata$day, levels = c('t0', 'day1', 'day3', 'day7', 'day14'))

colnames(tScores) <- sub("setA", "", colnames(tScores))
colnames(tScores) <- sub("_.*", "", colnames(tScores))
colnames(tScores) <- sub("CO$", "Co", colnames(tScores))
colnames(tScores) <- sub("DJ$", "Dj", colnames(tScores))


tScoresLong=tScores %>%
  select(-c(desc,
            sysName))%>%
  pivot_longer(cols = c(2:43), names_to = 'sample', values_to = 'logRatios')%>%
  left_join(metadata, by = 'sample')

heatMapSamples = tScoresLong %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  select(sample)%>%
  distinct()%>%
  .$sample

heatmapMeta = metadata %>%
  filter(sample %in% heatMapSamples) %>%
  mutate(mouseTissueDay = paste(mouse, tissue, day, sep = '-')) %>%
  select(mouseTissueDay, tissue, day) %>%
  column_to_rownames('mouseTissueDay')

heatmapMatrix = tScoresLong %>%
  mutate(mouseTissueDay = paste(mouse, tissue, day, sep = '-')) %>%
  filter(tissue != 'T0') %>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios') %>%
  column_to_rownames('locusId')

heatmapMeta = heatmapMeta[colnames(heatmapMatrix), , drop = FALSE]



ha=HeatmapAnnotation(df = heatmapMeta, col = list(tissue = c('dj' = 'khaki1', 'colon' = 'khaki4'), day = c('day1' = 'olivedrab1',
                                                                                                           'day3'='olivedrab2',
                                                                                                           'day7' = 'olivedrab3',
                                                                                                           'day14' = 'olivedrab4')))

tScoresLong %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios')%>%
  column_to_rownames('locusId')%>%
  Heatmap(show_column_names = F, show_row_names = F, bottom_annotation = ha, name ='t-score')

tScoresLong %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios')%>%
  column_to_rownames('locusId')%>%
  scale()%>%
  Heatmap(show_column_names = T, show_row_names = F)

# just going to do this for genes with strong effects at some point
strongFitnessEffects=read_tsv('barseqAdjustedParams/strong.tab')

tScoresLong %>%
  filter(locusId %in% strongFitnessEffects$locusId) %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios')%>%
  column_to_rownames('locusId')%>%
  Heatmap(show_column_names = T, show_row_names = F)

heatMapSamples = tScoresLong %>%
  filter(locusId %in% strongFitnessEffects$locusId) %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  select(sample)%>%
  distinct()%>%
  .$sample

heatmapMeta = metadata %>%
  filter(sample %in% heatMapSamples) %>%
  mutate(mouseTissueDay = paste(mouse, tissue, day, sep = '-')) %>%
  select(mouseTissueDay, tissue, day) %>%
  column_to_rownames('mouseTissueDay')

heatmapMatrix = tScoresLong %>%
  filter(locusId %in% strongFitnessEffects$locusId) %>%
  mutate(mouseTissueDay = paste(mouse, tissue, day, sep = '-')) %>%
  filter(tissue != 'T0') %>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios') %>%
  column_to_rownames('locusId')

heatmapMeta = heatmapMeta[colnames(heatmapMatrix), , drop = FALSE]



ha=HeatmapAnnotation(df = heatmapMeta, col = list(tissue = c('dj' = 'khaki1', 'colon' = 'khaki4'), day = c('day1' = 'olivedrab1',
                                                                                                           'day3'='olivedrab2',
                                                                                                           'day7' = 'olivedrab3',
                                                                                                           'day14' = 'olivedrab4')))

tScoresLong %>%
  filter(locusId %in% strongFitnessEffects$locusId) %>%
  mutate('mouseTissueDay' = paste(mouse,tissue,day, sep = '-'))%>%
  filter(tissue != 'T0')%>%
  pivot_wider(id_cols = locusId, names_from = mouseTissueDay, values_from = 'logRatios')%>%
  column_to_rownames('locusId')%>%
  Heatmap(show_column_names = F, show_row_names = F, bottom_annotation = ha, name ='Fitness\nscore')
