library(lme4)
library(lmerTest)
library(tidyverse)
library(car)

'There is a better way to choose mutants which are interesting than trying to pick a variance filter that seems
right. Need to add a function calculate the fold change between groups and choose the genes in a similar manner to the
RNA seq. '

run_lm = function(inputMatrix){
  loci = rownames(inputMatrix)
  modelStatsCategoricalDay <- data.frame(gene = character(0),
                                         pvalue = numeric(0),
                                         coef = numeric(length = 0),
                                         absoluteDelta = numeric(length = 0)
  )
  for (l in 1:length(loci)) {
    input <- inputMatrix[rownames(inputMatrix) == loci[l], ] %>%
      pivot_longer(cols = c(1:ncol(.)),
                   names_to = 'sample',
                   values_to = 'fitness') %>%
      merge(colonMeta, by = 'sample') %>%
      filter(tissue != 'T0')

    input$millionBases <- scale(input$millionBases)

    lmer_fit <- lmer(
      fitness ~ day + percentPerfectBarcode + millionBases + (1|day:cage),
      data = input
    )

    coefs <- summary(lmer_fit)$coefficients
    pvalue <- coefs[2, "Pr(>|t|)"]
    coef = coefs[2,1]
    absoluteDelta = calculateMeanDeltaFitness(input)

    modelStatsCategoricalDay <- rbind(modelStatsCategoricalDay, data.frame(gene = loci[l],
                                                                           pvalue = pvalue,
                                                                           coef = coef,
                                                                           absoluteDelta = absoluteDelta
    ))

  }
  return(modelStatsCategoricalDay)
}

calculateMeanDeltaFitness = function(inputDf){
  groupMeans = inputDf %>%
    group_by(day) %>%
    summarise(groupMean = mean(fitness))
  absoluteDelta = abs(groupMeans$groupMean[1]-groupMeans$groupMean[2])
  return(absoluteDelta)
}


get_module <- function(accession) {
  keggout <- keggGet(accession)
  if (!is.null(keggout[[1]]$BRITE)) {
    brite <- keggout[[1]]$BRITE
    level1 <- if (length(brite) >= 2) as.character(brite[2]) else NA_character_
    level2 <- if (length(brite) >= 3) as.character(brite[3]) else NA_character_
    return(c(level1, level2))
  } else {
    return(c(NA_character_, NA_character_))
  }
  return
}

annotate_modules <- function(df) {
  df_no_na <- df %>% filter(!is.na(kofamAccession)) %>%
    mutate(moduleLevel1 = NA_character_, moduleLevel2 = NA_character_)

  with_progress({
    p <- progressor(steps = nrow(df_no_na))
    for (i in seq_len(nrow(df_no_na))) {
      p()
      keggOut <- get_module(df_no_na$kofamAccession[i])
      df_no_na$moduleLevel1[i] <- keggOut[1]
      df_no_na$moduleLevel2[i] <- keggOut[2]
    }
  })

  df_annotated <- df %>%
    left_join(df_no_na %>% select(kofamAccession, moduleLevel1, moduleLevel2),
              by = "kofamAccession")

  return(df_annotated)
}

'Need to test and make sure it is actually better than the basic LM'

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
  filter(tissue == 'colon')

passVarainceFilter=fitnessScores%>%
  select(-c(sysName,
            desc))%>%
  pivot_longer(cols = 2:ncol(.),
               values_to = 'fitnessScores',
               names_to = 'sample')%>%
  filter(sample %in% colonMeta$sample)%>%
  group_by(locusId)%>%
  summarise(locusVariance = var(fitnessScores))%>%
  filter(locusVariance > 0.5)


loci = passVarainceFilter$locusId

lmIn = fitnessScores%>%
  filter(locusId %in% loci)%>%
  select(-c(sysName,
            desc))%>%
  column_to_rownames('locusId')

modelStatsNumericDay <- data.frame(gene = character(0), pvalue = numeric(0), coef = numeric(length = 0))
for (l in 1:length(loci)) {
  # Extract the input data for the specific gene
  input = lmIn[rownames(lmIn) == loci[l], ]

  # Reshape the data and merge with metadata
  input = input %>%
    pivot_longer(cols = c(1:ncol(.)),
                 names_to = 'sample',
                 values_to = 'fitness') %>%
    merge(colonMeta, by = 'sample') %>%
    filter(tissue != 'T0')
  input$tissue = as.factor(input$tissue)
  input$millionBases = scale(input$millionBases)
  # Fit the model and capture the summary
  lmer_fit <- lmer(
    fitness ~ dayNumeric + percentPerfectBarcode + millionBases + (1|cage),
    data = input
  )

  lmer_summary <- summary(lmer_fit)

  coefs <- summary(lmer_fit)$coefficients

  pvalue <- coefs[2, "Pr(>|t|)"]

  coef =coefs[2,1]
  modelStatsNumericDay <- rbind(modelStatsNumericDay, data.frame(gene = loci[l], pvalue = pvalue, coef = coef))
}
#write_tsv(modelStats, 'linear models/lmTissueComparisionsAllTimePoints/mixedEffectsModelRes/fitnessCageNumericDayAndPercentPerfectBarcodeAsPredictorMouseAsRandomModelStats.tsv')

modelStatsNumericDay$padj = p.adjust(modelStatsNumericDay$pvalue, method = 'fdr')

sigNumericDay=modelStatsNumericDay%>%
  filter(padj < .1)

'Model is working well'
modelStatsNumericDay%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Dist of raw p-values numeric day colon')

sigNumericDay = sigNumericDay %>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in unique(sigNumericDay$locusId[1:100])){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(colonMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigNumericDay$coef[sigNumericDay$locusId==s],
                          '\np-adj :', sigNumericDay$padj[sigNumericDay$locusId==s]))
  plot(p)
}


day1Day3Samples=colonMeta%>%
  filter(day %in% c('day1', 'day3'))%>%
  .$sample

day1Day3Fitness=lmIn[colnames(lmIn)%in% day1Day3Samples]

day1Day3LmRes=run_lm(day1Day3Fitness)

day1Day3LmRes$padj = p.adjust(day1Day3LmRes$pvalue, method = 'fdr')

day1Day3LmRes%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 50)

day1Day3LmRes%>%
  ggplot(aes(x = absoluteDelta,
             y = -log10(padj)))+
  geom_point()

day1Day3LmRes%>%
  filter(padj < .1)

sigDay1Day3= day1Day3LmRes %>%
  filter(padj < .1,
         absoluteDelta> .5)%>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in sigDay1Day3$locusId){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(colonMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           day %in% c('day1', 'day3'),
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigDay1Day3$coef[sigDay1Day3$locusId==s],
                          '\np-adj :', sigDay1Day3$padj[sigDay1Day3$locusId==s]))
  plot(p)
}

day3Day7Samples=colonMeta%>%
  filter(day %in% c('day3', 'day7'))%>%
  .$sample

day3Day7Fitness=lmIn[colnames(lmIn)%in% day3Day7Samples]

day3Day7LmRes=run_lm(day3Day7Fitness)

day3Day7LmRes$padj = p.adjust(day3Day7LmRes$pvalue, method = 'fdr')

day3Day7LmRes%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 50)

day3Day7LmRes%>%
  ggplot(aes(x = absoluteDelta,
             y = -log10(padj)))+
  geom_point()

day3Day7LmRes%>%
  filter(padj < .1)

sigDay3Day7= day3Day7LmRes %>%
  filter(padj < .1,
         absoluteDelta > .75)%>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in sigDay3Day7$locusId[1:100]){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(colonMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           day %in% c('day3', 'day7'),
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigDay3Day7$coef[sigDay3Day7$locusId==s],
                          '\np-adj :', sigDay3Day7$padj[sigDay3Day7$locusId==s]))
  plot(p)
}

day7Day14Samples=colonMeta%>%
  filter(day %in% c('day7', 'day14'))%>%
  .$sample

day7Day14Fitness=lmIn[colnames(lmIn)%in% day7Day14Samples]

day7Day14LmRes=run_lm(day7Day14Fitness)

day7Day14LmRes$padj = p.adjust(day7Day14LmRes$pvalue, method = 'fdr')

day7Day14LmRes%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 50)

day7Day14LmRes%>%
  filter(padj < .1)

sigDay7Day14= day7Day14LmRes %>%
  filter(padj < .1,
         absoluteDelta > .75)%>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in sigDay7Day14$locusId[1:100]){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(colonMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           day %in% c('day7', 'day14'),
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigDay7Day14$coef[sigDay7Day14$locusId==s],
                          '\np-adj :', sigDay7Day14$padj[sigDay7Day14$locusId==s]))
  plot(p)
}

sigCategoricalDay=sigDay1Day3%>%
  rbind(sigDay3Day7)%>%
  rbind(sigDay7Day14)

sigCategoricalDay %>%
  select(-absoluteDelta)%>%
  rbind(sigNumericDay)%>%
  write_tsv('linear models/mixedEffectsColonOnly/sigRes.tsv')