library(tidyverse)
library(ggpubr)
library(Maaslin2)
library(ComplexHeatmap)
library(vegan)

fullLocusFitness<-read_tsv('full/fit_t.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$numericDay = 0
metadata$numericDay[metadata$day== 'day1'] = 1
metadata$numericDay[metadata$day== 'day3'] = 3
metadata$numericDay[metadata$day== 'day 7'] = 7
metadata$numericDay[metadata$day== 'day 14'] = 14

genes<-read_tsv('full/genes')
genes$begin = as.integer(genes$begin)
functions <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/bBreveFunctions.tsv')

kofams=functions %>%
  filter(source == 'KOfam') %>%
  dplyr::select(gene_callers_id,
                `function`,
                accession) %>%
  rename('kofamFunction' = `function`,
         'kofamAccession' = accession)

pfams=functions %>%
  filter(source == 'Pfam') %>%
  dplyr::select(gene_callers_id,
                `function`,
                accession) %>%
  rename('pfamFunction' = `function`,
         'pfamAcession' = accession)

CAZyme=functions %>%
  filter(source == 'CAZyme') %>%
  dplyr::select(gene_callers_id,
                `function`,
                accession) %>%
  rename('CAZymeFunction' = `function`,
         'cazymeAcession' =  accession)
geneCallers <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/bBreveGeneCalls.tsv')
geneCallers<-geneCallers %>%
  rename('begin' = start,
         'end' = stop)
geneCallers$begin = geneCallers$begin+1

geneCallersMerge=geneCallers %>%
  dplyr::select(begin, gene_callers_id, aa_sequence)

geneCallersMerge= geneCallersMerge %>%
  left_join(kofams, by = 'gene_callers_id') %>%
  left_join(pfams,  by = 'gene_callers_id')


genes = genes %>%
  left_join(geneCallersMerge, by = 'begin')

write_tsv(genes, 'genesWithAnvioAnnotations.tsv')

metadata$numericDay = as.integer(metadata$numericDay)

metadata=metadata%>%
  mutate(mouseDay = paste(mouse, day, sep = '-'))

metadata=metadata%>%
  mutate(mouseDayTissue = paste(mouse, day, tissue,sep = '-'))


fullLocusFitnessTransposed<-fullLocusFitness %>%
  pivot_longer(cols = 4:45, values_to = 'fitnessScore', names_to = 'sample') %>%
  left_join(metadata, by ='sample')

fullLocusFitnessScores<-fullLocusFitnessTransposed %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(fullLocusFitnessScores)
adonisData<-fullLocusFitnessTransposed %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore)

adonisMetadata<-metadata %>%
  filter(sample %in%adonisData$sample)
adonisData = adonisData %>%
  column_to_rownames('sample')
adonis2(adonisData~day+tissue, data= adonisMetadata, method = 'euc')
fullLocusFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = tissue,shape = day, label = mouseDayTissue))+
  geom_point()+
  labs(title = 'Gene fitness PCA SI vs LI',
       x = 'PC1 46.71% var',
       caption = 'p < .001, PERMANOVA\nR^2 = 0.38887',
       y = 'PC2 20.64% var')

fullLocusFitnessScoresNoT0<-fullLocusFitnessTransposed %>%
  filter(tissue == 'dj' | tissue == 'colon')%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()

fullLocusFitnessScoresNoT0$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  mutate(tissueDay = paste(tissue, day, sep ='-'))%>%
  ggplot(aes(x = PC1, y = PC2,
             col = tissueDay,
             label = sample))+
  geom_point()+
  stat_ellipse()

fullLocusFitnessTransposed %>%
  pivot_wider(names_from = mouseDayTissue,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  as.matrix()%>%
  Heatmap(show_row_names = F)

fullLocusFitnessTransposed %>%
  pivot_wider(names_from = mouseDayTissue,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  scale()%>%
  as.matrix()%>%
  Heatmap(show_row_names = F)

adonisDataColon = fullLocusFitnessTransposed %>%
  filter(tissue == 'colon' | tissue == 'T0') %>%
  pivot_wider(names_from = locusId,id_cols = sample,
              values_from =  fitnessScore)

adonisMetadataColon = metadata %>%
  filter(sample %in%adonisDataColon$sample)

adonisDataColon = adonisDataColon %>%
  column_to_rownames('sample')
adonis2(adonisDataColon ~ day, data = adonisMetadataColon, method = 'euc')
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
  labs(title = 'Gene fitness PCA Colon',
       x = 'PC1 52.72% var',
       caption = 'p < .001, PERMANOVA\nR^2 = 0.47774',
       y = 'PC2 14.36% var')

fullLocusFitnessTransposed %>%
  pivot_wider(names_from = mouseDayTissue,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  as.matrix()
  Heatmap(show_row_names = F)


fullLocusFitnessScoresDj=fullLocusFitnessTransposed %>%
  filter(tissue == 'dj') %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(fullLocusFitnessScoresDj)

adonisDataSi = fullLocusFitnessTransposed %>%
  filter(tissue == 'dj') %>%
  pivot_wider(names_from = locusId,id_cols = sample,
              values_from =  fitnessScore)

adonisMetadataSi = metadata %>%
  filter(sample %in%adonisDataSi$sample)

adonisDataColon = adonisDataColon %>%
  column_to_rownames('sample')
adonis2(adonisDataColon ~ day, data = adonisMetadataColon, method = 'euc')

fullLocusFitnessScoresDj$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = day, label = sample))+
  geom_point()+
  labs(title = 'Gene fitness PCA SI',
       x = 'PC1 54.62% var',
       caption = 'p = 0.004, PERMANOVA\nR^2 =0.47774',
       y = 'PC2 16.02% var')

fullLocusFitnessTransposed %>%
  filter(tissue == 'dj' | tissue == 'T0') %>%
  pivot_wider(names_from = mouseDay,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  as.matrix()%>%
  Heatmap(show_row_names = F)


fullLocusFitnessTransposed %>%
  filter(tissue == 'dj') %>%
  pivot_wider(names_from = mouseDay,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  as.matrix()%>%
  Heatmap(show_row_names = F)

maaslinDataSI<-fullLocusFitnessTransposed %>%
  filter(tissue == 'dj') %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')

maaslinSIMeta=metadata %>%
  filter(tissue == 'dj') %>%
  column_to_rownames('sample')

# fit_data = Maaslin2(
#   input_data = maaslinDataSI,
#   input_metadata = maaslinSIMeta,
#   output = "siFitnessScores",
#   normalization = 'none',
#   transform = 'none',
#   min_abundance = 0,
#   min_prevalence = 0,
#   fixed_effects = c("numericDay"),
#   reference = c('day,day1'))

# fit_data = Maaslin2(
#   input_data = maaslinDataSI,
#   input_metadata = maaslinSIMeta,
#   output = "siFitnessScoresCategoricalData",
#   normalization = 'none',
#   transform = 'none',
#   min_abundance = 0,
#   min_prevalence = 0,
#   fixed_effects = c("day"),
#   reference = c('day,day1'))

maaslinDataCo<-fullLocusFitnessTransposed %>%
  filter(tissue == 'colon') %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')

maaslinCoMeta=metadata %>%
  filter(tissue == 'colon') %>%
  column_to_rownames('sample')

# fit_data = Maaslin2(
#   input_data = maaslinDataCo,
#   input_metadata = maaslinCoMeta,
#   output = "CoFitnessScores",
#   normalization = 'none',
#   transform = 'none',
#   min_abundance = 0,
#   min_prevalence = 0,
#   fixed_effects = c("numericDay"))

fullLocusFitness %>%
  select(locusId, desc)
coMaaslinRes = read_tsv('CoFitnessScores/all_results.tsv')
coMaaslinResWithFunctions<-coMaaslinRes %>%
  rename('locusId'=feature)%>%
  left_join(select(fullLocusFitness,locusId, desc), by = 'locusId')%>%
  as.data.frame()

coMaaslinResWithFunctions=coMaaslinResWithFunctions %>%
  merge(genes, by = 'locusId')

coMaaslinResWithFunctions %>%
  filter(qval < .05) %>%
  arrange(qval)%>%
  view()

coMaaslinResWithFunctions$col <- NA
coMaaslinResWithFunctions$col[coMaaslinResWithFunctions$coef > 0 & coMaaslinResWithFunctions$qval < .05] <- 'Positive correlation'
coMaaslinResWithFunctions$col[coMaaslinResWithFunctions$coef < 0 & coMaaslinResWithFunctions$qval < .05] <- 'Negative correlation'
coMaaslinResWithFunctions %>%
  ggplot(aes(x = coef,
             col = col,
             y = -log10(qval)))+
  geom_hline(yintercept = -log10(.05))+
  geom_point(alpha = .25,)+
  labs(col = 'Impact of mutation',
       alpha = '',
       title = 'Mutations are assoicated with different fitness\
       outcomes temporally in the colon')


siMaaslinRes = read_tsv('siFitnessScores/all_results.tsv')
siMaaslinResWithFunctions<-siMaaslinRes %>%
  rename('locusId'=feature)%>%
  left_join(select(fullLocusFitness,locusId, desc), by = 'locusId')%>%
  as.data.frame()
siMaaslinResWithFunctions$col <- NA
siMaaslinResWithFunctions$col[siMaaslinResWithFunctions$coef > 0 & siMaaslinResWithFunctions$qval < .05] <- 'Positive correlation'
siMaaslinResWithFunctions$col[siMaaslinResWithFunctions$coef < 0 & siMaaslinResWithFunctions$qval < .05] <- 'Negative correlation'
siMaaslinResWithFunctions %>%
  ggplot(aes(x = coef,
             col = col,
             y = -log10(qval)))+
  geom_point(alpha = .25,)+
  geom_hline(yintercept = -log10(.05))+
  labs(col = 'Impact of mutation',
       alpha = '',
       title = 'Mutations are assoicated with different fitness\
       outcomes temporally in the SI')

siSigNegatveFitness = siMaaslinResWithFunctions %>%
  filter(coef < 0,
         qval < .05)

coSigNegatveFitness = coMaaslinResWithFunctions %>%
  filter(coef < 0,
         qval <.1)

#319 low fitness mutants shared between conditions over time
intersect(siSigNegatveFitness$locusId,
          coSigNegatveFitness$locusId)%>%
  length()

sharedLowFitness<-intersect(siSigNegatveFitness$locusId,
          coSigNegatveFitness$locusId)

siSigNegatveFitness %>%
  filter(locusId %in% sharedLowFitness)

siSigNegatveFitness %>%
  filter(!locusId %in% sharedLowFitness)

coSigNegatveFitness %>%
  filter(!locusId %in% sharedLowFitness)

day7Scores<-fullLocusFitnessTransposed %>%
  filter(day =='day 7')%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(day7Scores)
day7Scores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC2, y = PC3, col = tissue, shape = day,label = sample))+
  geom_point()