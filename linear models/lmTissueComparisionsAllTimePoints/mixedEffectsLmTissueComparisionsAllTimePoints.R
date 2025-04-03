library(lme4)
library(lmerTest)
library(tidyverse)
library(car)

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
modelStats <- data.frame(gene = character(0), pvalue = numeric(0), aic = numeric(length = 0))
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


  # Fit the model and capture the summary
  lmer_fit <- lmer(
    fitness ~ tissue+lane+millionBases+(1|day/cage/mouse),
    data = input
  )

  lmer_summary <- summary(lmer_fit)

  coefs <- summary(lmer_fit)$coefficients
  pvalue <- coefs[2, "Pr(>|t|)"]


  aic <- AIC(lmer_fit)

  modelStats <- rbind(modelStats, data.frame(gene = loci[l], pvalue = pvalue, aic = aic))
}
#write_tsv(modelStats, 'linear models/lmTissueComparisionsAllTimePoints/mixedEffectsModelRes/fitnessCageNumericDayAndPercentPerfectBarcodeAsPredictorMouseAsRandomModelStats.tsv')

modelStats$padj = p.adjust(modelStats$pvalue, method = 'fdr')


 modelStats%>%
   ggplot(aes(x = pvalue))+
   geom_histogram(bins = 150)

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

tissueScatterIn$col = NA
tissueScatterIn$col[tissueScatterIn$locusId %in% sig$gene]= 'Significant'

tissueScatterIn%>%
  ggplot(aes(x = dj,
             col = col,
             y = colon))+
  geom_point(alpha = .15)+
  geom_abline()

plot(lmer_fit)

qqnorm(residuals(lmer_fit))

qqline(residuals(lmer_fit))

sig = sig %>%
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


