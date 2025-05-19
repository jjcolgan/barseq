library(lme4)
library(lmerTest)
library(tidyverse)
library(car)

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

djMeta = metadata%>%
  filter(tissue == 'dj')


fitnessScores%>%
  select(-c(sysName,
            desc))%>%
  pivot_longer(cols = 2:ncol(.),
               values_to = 'fitnessScores',
               names_to = 'sample')%>%
  filter(sample %in% djMeta$sample)%>%
  group_by(locusId)%>%
  summarise(locusVariance = var(fitnessScores))%>%
  ggplot(aes(x =locusVariance))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept= .3)

passVarainceFilter=fitnessScores%>%
  select(-c(sysName,
            desc))%>%
  pivot_longer(cols = 2:ncol(.),
               values_to = 'fitnessScores',
               names_to = 'sample')%>%
  filter(sample %in% djMeta$sample)%>%
  group_by(locusId)%>%
  summarise(locusVariance = var(fitnessScores))%>%
  filter(locusVariance > 0)


loci = passVarainceFilter$locusId

lmIn = fitnessScores%>%
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
    merge(djMeta, by = 'sample') %>%
    filter(tissue != 'T0')
  input$tissue = as.factor(input$tissue)
  input$millionBases = scale(input$millionBases)
  # Fit the model and capture the summary
  lmer_fit <- lmer(
    fitness ~ dayNumeric + percentPerfectBarcode + millionBases +  (1|day:cage),
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
  labs(title = 'Dist of raw p-values numeric day DJ')

sigNumericDay = sigNumericDay %>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in unique(sigNumericDay$locusId[1:100])){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(djMeta, by = 'sample')%>%
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


modelStatsCategoricalDay <- data.frame(gene = character(0), pvalue = numeric(0), coef = numeric(length = 0))
for (l in 1:length(loci)) {
  input <- lmIn[rownames(lmIn) == loci[l], ] %>%
    pivot_longer(cols = c(1:ncol(.)),
                 names_to = 'sample',
                 values_to = 'fitness') %>%
    merge(djMeta, by = 'sample') %>%
    filter(tissue != 'T0')

  input$tissue <- as.factor(input$tissue)
  input$millionBases <- scale(input$millionBases)

  lmer_fit <- lmer(
    fitness ~ day + percentPerfectBarcode + millionBases + (1|day:cage),
    data = input
  )

  coefs <- summary(lmer_fit)$coefficients

  # Find rows for all 'day' levels
  day_rows <- grep("^day", rownames(coefs))

  if (length(day_rows) > 0) {
    for (d in day_rows) {
      modelStatsCategoricalDay <- rbind(
        modelStatsCategoricalDay,
        data.frame(
          gene   = loci[l],
          term   = rownames(coefs)[d],
          coef   = coefs[d, "Estimate"],
          pvalue = coefs[d, "Pr(>|t|)"]
        )
      )
    }
  }
}
#write_tsv(modelStats, 'linear models/lmTissueComparisionsAllTimePoints/mixedEffectsModelRes/fitnessCageCategoricalDayAndPercentPerfectBarcodeAsPredictorMouseAsRandomModelStats.tsv')

modelStatsCategoricalDay$padj = p.adjust(modelStatsCategoricalDay$pvalue, method = 'fdr')

sigCategoricalDay=modelStatsCategoricalDay%>%
  filter(padj < .1)

'Model is working well'
modelStatsCategoricalDay%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Dist of raw p-values categorical day DJ')

sigCategoricalDay = sigCategoricalDay %>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in unique(sigCategoricalDay$locusId[1:100])){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(djMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigCategoricalDay$coef[sigCategoricalDay$locusId==s],
                          '\np-adj :', sigCategoricalDay$padj[sigCategoricalDay$locusId==s]))
  plot(p)
}

djFitnessScores=fitnessScores%>%
  select(-c(sysName, desc))%>%
  filter(locusId %in% sigCategoricalDay$locusId | locusId %in% sigNumericDay$locusId)%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% djMeta$sample)

write_tsv(djFitnessScores, file = '/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/linear models/mixedEffectsJejunumOnly/mixedEffectsRes.tsv')