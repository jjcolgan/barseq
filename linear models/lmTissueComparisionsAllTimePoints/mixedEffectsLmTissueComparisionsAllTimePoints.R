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

'starting with the best linear model co-variates, which are just numeric day and percent perfect barcode.'
'It makes sense to have mouse as a nested effect of cage, since only certain mice appear in certain cages'
modelStats <- data.frame(gene = character(0), pvalue = numeric(0), coef = numeric(length = 0))
for (l in 1:length(loci)) {
  # Extract the input data for the specific gene
  input = lmIn[rownames(lmIn) == loci[l], ]

  # Reshape the data and merge with metadata
  input = input %>%
    pivot_longer(cols = c(1:ncol(.)),
                 names_to = 'sample',
                 values_to = 'fitness') %>%
    merge(metadata, by = 'sample') %>%
    filter(tissue != 'T0')
  input$tissue = as.factor(input$tissue)
  input$millionBases = scale(input$millionBases)
  # Fit the model and capture the summary
  lmer_fit <- lmer(
    fitness ~ tissue + lane + millionBases + day + (1|day:cage/mouse),
    data = input
  )

  lmer_summary <- summary(lmer_fit)

  coefs <- summary(lmer_fit)$coefficients

  pvalue <- coefs[2, "Pr(>|t|)"]

  coef =coefs[2,1]
  modelStats <- rbind(modelStats, data.frame(gene = loci[l], pvalue = pvalue, coef = coef))
}
#write_tsv(modelStats, 'linear models/lmTissueComparisionsAllTimePoints/mixedEffectsModelRes/fitnessCageNumericDayAndPercentPerfectBarcodeAsPredictorMouseAsRandomModelStats.tsv')

modelStats$padj = p.adjust(modelStats$pvalue, method = 'fdr')


 modelStats%>%
   ggplot(aes(x = pvalue))+
   geom_histogram(bins = 100)+
   labs(title = 'Histogram of unadjusted pvalues')
sig = modelStats %>%
  filter(padj < .1)

tissueScatterIn=fitnessScores%>%
  pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')%>%
  filter(tissue != 'T0')%>%
  group_by(locusId, tissue)%>%
  summarise(meanTissueFitness = mean(fitnessScore))%>%
  pivot_wider(id_cols = 'locusId',
              names_from = 'tissue',
              values_from = meanTissueFitness)
tissueScatterIn$col = 'P-adjust > .1'
tissueScatterIn$col[tissueScatterIn$locusId %in% sig$gene]= 'P-adjust < 0.1'

tissueScatterIn%>%
  ggplot(aes(x = dj,
             col = col,
             y = colon))+
  geom_point(alpha = .15)+
  geom_abline()+
  labs(col = 'Significance')


plot(lmer_fit)

qqnorm(residuals(lmer_fit))

qqline(residuals(lmer_fit))

sig = sig %>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in unique(sig$locusId)){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(metadata, by = 'sample')%>%
    filter(tissue != 'T0',
           locusId == s)%>%
    ggplot(aes(x = tissue, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sig$coef[sig$locusId==s],
                          '\np-adj :', sig$padj[sig$locusId==s]))
  plot(p)
}


modelStats %>%
  ggplot(aes(y = -log10(padj),
             x = coef))+
  geom_point()

modelStats%>%
  write_tsv('linear models/lmTissueComparisionsAllTimePoints/mixedEffectsModelRes/mixedEffectsModelRes.tsv')

higherMeanDj=sig %>%
  filter(padj < .1,
         coef > 0)

higherMeanColon=sig %>%
  filter(padj < .1,
         coef < 0)

higherMeanDj=annotate_modules(higherMeanDj)

higherMeanDj %>%
  mutate(moduleLevel1 = str_remove(moduleLevel1, "^\\s*[[:digit:]]+\\s+"), # Remove leading spaces & number
         moduleLevel2 = str_remove(moduleLevel2, "^\\s*[[:digit:]]+\\s+"),
         moduleLevel1 = wrap_format(20)(moduleLevel1),  # Wrap text for readability
         moduleLevel1 = fct_reorder(moduleLevel1, .x = moduleLevel1, .fun = function(x) length(x), .desc = TRUE)) %>%  # Order by count
  ggplot(aes(x = moduleLevel1, fill = moduleLevel2)) +
  geom_bar() +
  coord_flip() +  # Flip bars for better readability
  theme_minimal()+
  labs(x = 'Count',
       y = 'Brite level 1',
       fill = 'Brite level 2')

higherMeanColon=annotate_modules(higherMeanColon)

higherMeanColon %>%
  mutate(moduleLevel1 = str_remove(moduleLevel1, "^\\s*[[:digit:]]+\\s+"), # Remove leading spaces & number
         moduleLevel2 = str_remove(moduleLevel2, "^\\s*[[:digit:]]+\\s+"),
         moduleLevel1 = wrap_format(20)(moduleLevel1),  # Wrap text for readability
         moduleLevel1 = fct_reorder(moduleLevel1, .x = moduleLevel1, .fun = function(x) length(x), .desc = TRUE)) %>%  # Order by count
  ggplot(aes(x = moduleLevel1, fill = moduleLevel2)) +
  geom_bar() +
  coord_flip() +  # Flip bars for better readability
  theme_minimal()+
  labs(x = 'Count',
       y = 'Brite level 1',
       fill = 'Brite level 2')