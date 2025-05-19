library(tidyverse)
library(lme4)
library(lmerTest)

fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$day = factor(metadata$day, levels = c('t0', 'day1', 'day3', 'day7', 'day14'))
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

metadata$phase = NA
metadata$phase[metadata$dayNumeric < 7] ='early'
metadata$phase[metadata$dayNumeric >= 7] ='late'

metadata$cage = as.factor(metadata$cage)

keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction)%>%
  distinct()
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

lmIn=fitnessScores %>%
  select(-c(sysName, desc))%>%
  column_to_rownames('locusId')

loci = rownames(lmIn)
'lm being used to compare is fitness ~ tissue+lane+millionBases'
'tissue+lane+millionBases+(tissue|day/cage) results in 13% of features having a decreased aic,
 significantly better performance for 24% of features, and a median F score of 1.26'

'fitness ~ tissue+lane+millionBases+(1|day/cage/mouse) results in aic decrease for 9.8% of features,
better anova for 24%, and median f stat of 1.26'
modelComparisons = data.frame('gene' = character(),
                              'deltaAic'= double(),
                              'pvalue' = double(),
                              'F' = double())

for (l in 1:length(loci)) {
  # Extract the input data for the specific gene
  input = lmIn[rownames(lmIn) == loci[l], ]
  input = input %>%
    pivot_longer(cols = c(1:ncol(.)),
                 names_to = 'sample',
                 values_to = 'fitness') %>%
    merge(metadata, by = 'sample') %>%
    filter(tissue != 'T0')

  # Convert tissue to a factor and set up tissueNumeric
  input$tissue = as.factor(input$tissue)
  input$millionBases = scale(input$millionBases)

  glm_best = glm(data = input,
                 formula = fitness ~ tissue+lane+millionBases)

  lmer_fit <- lmer(
    fitness ~ tissue + lane + millionBases + day + (1|day:cage/mouse),
    data = input
  )


  aicRes=AIC(glm_best, lmer_fit)
  anovaRes=anova(glm_best, lmer_fit)
  previousModAic=aicRes[1,2]
  newModAic=aicRes[2,2]
  deltaAic = newModAic - previousModAic
  p = anovaRes[2,6]
  f = anovaRes[2,5]
  newComparison = data.frame('gene' = loci[l],
                             'deltaAic'= deltaAic,
                             'pvalue' = p,
                             'f' = f)
  modelComparisons = rbind(modelComparisons,
                           newComparison)
}

aicIncreasesAbove2=modelComparisons %>%
  filter(deltaAic <= 2 )%>%
  nrow()

print(paste0('Percent aic decreases over 2: ',aicIncreasesAbove2/nrow(modelComparisons)))

betterPerformanceAnova=modelComparisons %>%
  filter(pvalue < .05)%>%
  nrow()


print(paste0('Percent significant better performace anova: ',betterPerformanceAnova/nrow(modelComparisons)))

print(paste0('Median F stat ', median(modelComparisons$f, na.rm = T)))

linearModelRes = read_tsv('linear models/lmTissueComparisionsAllTimePoints/linearModelRes/linearModelRes.tsv')

mixedEffectsModelRes = read_tsv('linear models/lmTissueComparisionsAllTimePoints/mixedEffectsModelRes/mixedEffectsModelRes.tsv')

sigMixedEffectsModel = mixedEffectsModelRes%>%
  filter(padj < .1)

sigLinearModel = linearModelRes %>%
  filter(padj < .1)

sigResOnlyInMixedEffects=sigMixedEffectsModel%>%
  filter(! gene %in% sigLinearModel$gene)

'Looks good so far, lets recheck the functions of these genes. Revise the temporal model as well to incorperate these
effects'
for (s in sigResOnlyInMixedEffects$gene){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(metadata, by = 'sample')%>%
    filter(tissue != 'T0',
           locusId == s)%>%
    ggplot(aes(x = tissue, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigResOnlyInMixedEffects$coef[sigResOnlyInMixedEffects$gene==s],
                          '\np-adj :', sigResOnlyInMixedEffects$padj[sigResOnlyInMixedEffects$gene==s]))
  plot(p)
}