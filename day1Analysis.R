library(tidyverse)
library(ggpubr)
library(Maaslin2)
library(ComplexHeatmap)
library(vegan)

'This function generates a contingency table with columns "positive" and "negative".
The function should be used when testing to see if there is a difference in the fitness
effects regardless of their strength. Thus weak fitness scores |fitnessScore| < 2.5,
will be tested.'
genContingecyTableNegativeVersusPositiveSelection = function (geneMat){
  #input should be the full geneMat with metadata included.
  tissues=geneMat %>%
    dplyr::select(tissue) %>%
    distinct() %>%
    .$tissue
  contab = data.frame('positive' = integer(),
                      'neutral' = integer())
  for( i in 1:length(tissues)){
    #find samples with negative effect
    nPositive = geneMat %>%
      filter(tissue == tissues[i]) %>%
      filter(fitnessScore > 0) %>%
      nrow()
    #find samples with positive effect
    nNegative = geneMat %>%
      filter(tissue == tissues[i]) %>%
      filter(fitnessScore < 0) %>%
      nrow()
    temp = data.frame('positive' = nPositive,
                      'negative' = nNegative)
    row.names(temp) = tissues[i]
    contab = rbind(contab,
                   temp)
  }
  return(contab)
}

genContingecyTableNeutralVersusNegative = function (geneMat){
  #input should be the full geneMat with metadata included.
  tissues=geneMat %>%
    dplyr::select(tissue) %>%
    distinct() %>%
    .$tissue
  contab = data.frame('negative' = integer(),
                      'neutral' = integer())
  for( i in 1:length(tissues)){
    #find samples with negative effect
    nZero = geneMat %>%
      filter(tissue == tissues[i]) %>%
      filter(fitnessScore == 0) %>%
      nrow()
    #find samples with positive effect
    nNegative = geneMat %>%
      filter(tissue == tissues[i]) %>%
      filter(fitnessScore < 0) %>%
      nrow()
    temp = data.frame('negative' = nNegative,
                      'neutral' = nZero)
    row.names(temp) = tissues[i]
    contab = rbind(contab,
                   temp)
  }
  return(contab)
}

genContingecyTableNeutralVersusPositive = function (geneMat){
  #input should be the full geneMat with metadata included.
  tissues=geneMat %>%
    dplyr::select(tissue) %>%
    distinct() %>%
    .$tissue
  contab = data.frame('positive' = integer(),
                      'neutral' = integer())
  for( i in 1:length(tissues)){
    #find samples with negative effect
    nZero = geneMat %>%
      filter(tissue == tissues[i]) %>%
      filter(fitnessScore == 0) %>%
      nrow()
    #find samples with positive effect
    nPositive = geneMat %>%
      filter(tissue == tissues[i]) %>%
      filter(fitnessScore > 0) %>%
      nrow()
    temp = data.frame('positive' = nPositive,
                      'neutral' = nZero)
    row.names(temp) = tissues[i]
    contab = rbind(contab,
                   temp)
  }
  return(contab)
}

fullLocusFitness<-read_tsv('full/fit_t.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$numericDay = 0
metadata$numericDay[metadata$day== 'day1'] = 1
metadata$numericDay = as.integer(metadata$numericDay)
metadata=metadata%>%
  mutate(mouseDay = paste(mouse, day, sep = '-'))
metadata=metadata%>%
  mutate(mouseDayTissue = paste(mouse, day, tissue,sep = '-'))
day1Meta=metadata %>%
  filter(day == 'day1' | day == 't0')

genes=read_tsv('genesWithAnvioAnnotations.tsv')

fullLocusFitness<-read_tsv('full/fit_t.tab')

fullLocusFitnessTransposed<-fullLocusFitness %>%
  pivot_longer(cols = 4:45, values_to = 'fitnessScore', names_to = 'sample') %>%
  left_join(day1Meta, by ='sample')%>%
  filter(sample %in% day1Meta$sample)

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
adonis2(adonisData~tissue, data= adonisMetadata, method = 'euc')

fullLocusFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = tissue,shape = day, label = sample))+
  geom_point()+
  labs(title = 'Gene fitness PCA day 1 with T0',
       x = 'PC1 51.1% var',
       caption = 'p = .58, PERMANOVA\nR^2 = 0.12842',
       y = 'PC2 19.96% var')

fullLocusFitnessTransposedNoT0<-fullLocusFitnessTransposed %>%
  filter(day == 'day1')

fullLocusFitnessScoresNoT0<-fullLocusFitnessTransposedNoT0 %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()

summary(fullLocusFitnessScoresNoT0)
adonisData<-fullLocusFitnessTransposedNoT0 %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore)
adonisMetadata<-metadata %>%
  filter(sample %in%adonisData$sample)
adonisData = adonisData %>%
  column_to_rownames('sample')
adonis2(adonisData~tissue, data= adonisMetadata, method = 'euc')

fullLocusFitnessScoresNoT0$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = tissue,shape = day, label = sample))+
  geom_point()+
  labs(title = 'Gene fitness PCA day 1 no t0',
       x = 'PC1 53.97% var',
       caption = 'p < .01, PERMANOVA\nR^2 = 0.52705',
       y = 'PC2 12.14% var')


fullLocusFitnessTransposedNoT0%>%
  pivot_wider(names_from = mouseDayTissue,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  scale()%>%
  as.matrix()%>%
  Heatmap(show_row_names = F)

fullLocusFitnessTransposedNoT0%>%
  pivot_wider(names_from = mouseDayTissue,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  as.matrix()%>%
  Heatmap(show_row_names = F)


meanFitnessScores=fullLocusFitnessTransposedNoT0 %>%
  group_by(tissue, locusId)%>%
  summarise('meanFitnessScore' = mean(fitnessScore))

meanFitnessScores %>%
  pivot_wider(id_cols = 'locusId',
              names_from = 'tissue',
              values_from = 'meanFitnessScore') %>%
  ggplot(aes(x = dj,
             y = colon))+
  geom_point()+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  geom_vline(xintercept = 1.5)+
  geom_vline(xintercept = -1.5)


maaslinData = fullLocusFitnessTransposedNoT0%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')
maaslinMeta  = adonisMetadata %>%
  column_to_rownames('sample')

fit_data = Maaslin2(
  input_data = maaslinData,
  input_metadata = maaslinMeta,
  output = "day1DjVcolonGeneLocus",
  normalization = 'none',
  transform = 'none',
  min_abundance = 0,
  min_prevalence = 0,
  fixed_effects = c("tissue"))

maaslin2Res = read_tsv(file = 'day1DjVcolonGeneLocus/all_results.tsv')

maaslin2Res=maaslin2Res %>%
  rename(locusId = 'feature') %>%
  left_join(genes, by = 'locusId')

enrichedFitnessColonSig<-maaslin2Res %>%
  filter(qval < .05,
         coef < 0)

enrichedFitnessDjSig<-maaslin2Res %>%
  filter(qval < .05,
         coef > 0)

#fishers exact test for negative versus positive fitness effects
'This is probably going to require some more intelligent filtering. I also think in general,
the effect sizes are quite small, so will be hard to test.'
locusIds=fullLocusFitnessTransposedNoT0 %>%
  dplyr::select(locusId) %>%
  distinct() %>%
  .$locusId
output = data.frame('gene' = character(),
                    'pvalue' = double())
for (i in 1:length(locusIds)){
  genesOfInterest=fullLocusFitnessTransposedNoT0 %>%
    filter(locusId == locusIds[i])
  fisherResults=genesOfInterest %>%
    genContingecyTableNegativeVersusPositiveSelection()%>%
    fisher.test()

  temp = data.frame('gene' = locusIds[i],
                    'pvalue' = fisherResults$p.value)
  output = rbind(output,
                  temp)

}
output %>%
  filter(pvalue < .01)
output$padjust = p.adjust(output$pvalue, method = 'BH')
summary(output$padjust)

'Fishers test neutral versus positive selection.
Does not look like there is anything that is neutral in one population and beneifical in another.'
output = data.frame('gene' = character(),
                    'pvalue' = double())
for (i in 1:length(locusIds)){
  genesOfInterest=fullLocusFitnessTransposedNoT0 %>%
    filter(locusId == locusIds[i])
  fisherResults=genesOfInterest %>%
    genContingecyTableNeutralVersusPositive()%>%
    fisher.test()

  temp = data.frame('gene' = locusIds[i],
                    'pvalue' = fisherResults$p.value)
  output = rbind(output,
                 temp)

}


'Fishers negative versus neutral-
same as above.'
output = data.frame('gene' = character(),
                    'pvalue' = double())
for (i in 1:length(locusIds)){
  genesOfInterest=fullLocusFitnessTransposedNoT0 %>%
    filter(locusId == locusIds[i])
  fisherResults=genesOfInterest %>%
    genContingecyTableNeutralVersusNegative()%>%
    fisher.test()

  temp = data.frame('gene' = locusIds[i],
                    'pvalue' = fisherResults$p.value)
  output = rbind(output,
                 temp)

}