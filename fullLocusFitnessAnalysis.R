library(tidyverse)
library(ggpubr)
library(Maaslin2)

fullLocusFitness<-read_tsv('full/fit_t.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$numericDay = 0
metadata$numericDay[metadata$day== 'day1'] = 1
metadata$numericDay[metadata$day== 'day3'] = 3
metadata$numericDay[metadata$day== 'day 7'] = 7
metadata$numericDay[metadata$day== 'day 14'] = 14

metadata$numericDay = as.integer(metadata$numericDay)

fullLocusFitnessTransposed<-fullLocusFitness %>%
  pivot_longer(cols = 4:43, values_to = 'fitnessScore', names_to = 'sample') %>%
  left_join(metadata, by ='sample')

fullLocusFitnessScores<-fullLocusFitnessTransposed %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()
fullLocusFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = tissue,shape = day, label = sample))+
  geom_point()


fullLocusFitnessScoresColon<-fullLocusFitnessTransposed %>%
  filter(tissue == 'colon') %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(fullLocusFitnessScoresColon)
fullLocusFitnessScoresColon$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = day, label = sample))+
  geom_point()+
  stat_ellipse()

fullLocusFitnessScoresDj=fullLocusFitnessTransposed %>%
  filter(tissue == 'dj') %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(fullLocusFitnessScores)
fullLocusFitnessScoresDj$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = day, label = sample))+
  geom_point()


maaslinDataSI<-fullLocusFitnessTransposed %>%
  filter(tissue == 'dj') %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')

maaslinSIMeta=metadata %>%
  filter(tissue == 'dj') %>%
  column_to_rownames('sample')

fit_data = Maaslin2(
  input_data = maaslinDataSI,
  input_metadata = maaslinSIMeta,
  output = "siFitnessScores",
  normalization = 'none',
  transform = 'none',
  min_abundance = 0,
  min_prevalence = 0,
  fixed_effects = c("numericDay"))

maaslinDataCo<-fullLocusFitnessTransposed %>%
  filter(tissue == 'colon') %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')

maaslinCoMeta=metadata %>%
  filter(tissue == 'colon') %>%
  column_to_rownames('sample')

fit_data = Maaslin2(
  input_data = maaslinDataCo,
  input_metadata = maaslinCoMeta,
  output = "CoFitnessScores",
  normalization = 'none',
  transform = 'none',
  min_abundance = 0,
  min_prevalence = 0,
  fixed_effects = c("numericDay"))

fullLocusFitness %>%
  select(locusId, desc)
coMaaslinRes = read_tsv('CoFitnessScores/all_results.tsv')
coMaaslinResWithFunctions<-coMaaslinRes %>%
  rename('locusId'=feature)%>%
  left_join(select(fullLocusFitness,locusId, desc), by = 'locusId')%>%
  as.data.frame()

coMaaslinResWithFunctions %>%
  ggplot(aes(x = coef,
             y = -log10(qval)))+
  geom_point()+
  geom_hline(yintercept = -log10(.05))

siMaaslinRes = read_tsv('siFitnessScores/all_results.tsv')
siMaaslinResWithFunctions<-siMaaslinRes %>%
  rename('locusId'=feature)%>%
  left_join(select(fullLocusFitness,locusId, desc), by = 'locusId')%>%
  as.data.frame()

siMaaslinResWithFunctions %>%
  ggplot(aes(x = coef,
             y = -log10(qval)))+
  geom_point()+
  geom_hline(yintercept = -log10(.05))

siSigNegatveFitness = siMaaslinResWithFunctions %>%
  filter(coef < 0)

coSigNegatveFitness = coMaaslinResWithFunctions %>%
  filter(coef < 0)


intersect(sigNegatveFitness$locusId,
          coSigNegatveFitness$locusId)%>%
  length()