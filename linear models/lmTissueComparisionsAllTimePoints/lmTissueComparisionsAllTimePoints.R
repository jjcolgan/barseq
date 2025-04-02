library(tidyverse)
library(lme4)
library(car)
library(KEGGREST)
library(progressr)
library(stringr)
library(scales)

simpleQQPlot = function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)),
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
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

'Adding percentPerfectBarcode decreases all genes by 2%, percent lower significant of 41% and median F of
2.70'

'adding categorical day results in 80% decrease in aic, 58% lower significance by anova and median f of 3.85'

'adding numeric day results in 100% decrease over 2 in aic, 55% lower signficance by anova and median f of 5.48'

'million bases does not really increase anything, does result in 100% decrease in aic over 2'

'adding cage results in aic deacreases over 2 for 65% of genes, significantly decreased anova for 47% and median F
of 2.14'

'adding mouse results in aic decrease over 2 for 55% of genes, signficantly decreased anove for 47% and median f stat
of 1.49'

'adding lane decreases aic over 2 for all features, results in 59% of features showing better signifiance, and has a median f of 6.69'

'fitness~tissue + lane is the best model'

'Just being careful, not evidence that fitness~tissue + percentPerfectBarcode is any better than fitness~tissue + dayNumeric'

'adding percentPerfectBarcode decreases aic by 2 or more for all features, is signifcantly better for ~5% of features
and results in a median F stat of .54'

'i think tissue, lane, millionBases is the best model'

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

  glm_previousBest = glm(data = input, formula = fitness ~ tissue+lane+millionBases)
  glm_new = glm(data = input, formula = fitness ~ tissue+lane+pfClusters)

  vif(glm_previousBest)
  vif(glm_new)

  aicRes=AIC(glm_previousBest, glm_new)
  anovaRes=anova(glm_previousBest, glm_new)
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

'lets run it and look at the results '
output = data.frame('gene' = character(),
                    'pvalue'=double(),
                    'coef' = double())
for (l in 1:length(loci)) {
  # Extract the input data for the specific gene
  input = lmIn[rownames(lmIn) == loci[l], ]
  input = input %>%
    pivot_longer(cols = c(1:ncol(.)),
                 names_to = 'sample',
                 values_to = 'fitness') %>%
    merge(metadata, by = 'sample') %>%
    filter(tissue != 'T0')

  # Convert tissue to a factor
  input$tissue = as.factor(input$tissue)

  glm_best = glm(data = input, formula = fitness ~ tissue+lane+millionBases)
  glmRes=summary(glm_best)
  coef =glmRes$coefficients[2,1]
  p = glmRes$coefficients[2,4]

  output = rbind(output,
                 data.frame('gene' = loci[l],
                            'pvalue'=p,
                            'coef' = coef))
}

output$padj = p.adjust(output$pvalue, method = 'fdr')
summary(output$padj)

sig=output%>%
  filter(padj < .1)%>%
  arrange(padj)



output %>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Histogram of unadjusted pvalues')+
  theme_bw()

output %>%
  ggplot(aes(x = pvalue))+
  geom_density()

simpleQQPlot(output$pvalue)

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
tissueScatterIn$col[tissueScatterIn$locusId %in% sig$gene]= 'P-adjust < .1'

tissueScatterIn%>%
  ggplot(aes(x = dj,
             col = col,
             y = colon))+
  geom_point(alpha = .15)+
  geom_abline()+
  theme_bw()+
  labs(col = 'Significance')


sig = sig %>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')

view(sig)

for (s in sig$locusId){
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

hist(output$pvalue)

higherMedianFitnessColon = sig%>%
  filter(coef < 0)%>%
  arrange(desc(-log10(padj)),
          desc(abs(coef)))

higherMedianFitnessDj = sig%>%
  filter(coef>0)%>%
  arrange(desc(-log10(padj)),
          desc(coef))

view(higherMedianFitnessColon)
view(higherMedianFitnessDj)

higherMedianFitnessDj=annotate_modules(higherMedianFitnessDj)

higherMedianFitnessDj %>%
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




higherMedianFitnessColon=annotate_modules(higherMedianFitnessColon)
higherMedianFitnessColon %>%
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

title=higherMedianFitnessDj%>%
  filter(locusId == 'BBR_RS13450')%>%
  .$kofamFunction

fitnessScores%>%
  pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')%>%
  filter(tissue != 'T0',
         locusId == 'BBR_RS13450')%>%
  ggplot(aes(x = tissue, y = fitnessScore))+
  geom_boxplot(outliers = F)+
  geom_point(aes(col = day))+
  labs(title = title,
       caption = paste0('Coef: ', sig$coef[sig$locusId=='BBR_RS13450'],
                        '\np-adj :', sig$padj[sig$locusId=='BBR_RS13450']))+
  theme_bw()

title=higherMedianFitnessDj%>%
  filter(locusId == 'BBR_RS15455')%>%
  .$kofamFunction

fitnessScores%>%
  pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')%>%
  filter(tissue != 'T0',
         locusId == 'BBR_RS15455')%>%
  ggplot(aes(x = tissue, y = fitnessScore))+
  geom_boxplot(outliers = F)+
  geom_point(aes(col = day))+
  labs(title = title,
       caption = paste0('Coef: ', sig$coef[sig$locusId=='BBR_RS15455'],
                        '\np-adj :', sig$padj[sig$locusId=='BBR_RS15455']))+
  theme_bw()

title=higherMedianFitnessDj%>%
  filter(locusId == 'BBR_RS10755')%>%
  .$kofamFunction

fitnessScores%>%
  pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')%>%
  filter(tissue != 'T0',
         locusId == 'BBR_RS10755')%>%
  ggplot(aes(x = tissue, y = fitnessScore))+
  geom_boxplot(outliers = F)+
  geom_point(aes(col = day))+
  labs(title = title,
       caption = paste0('Coef: ', sig$coef[sig$locusId=='BBR_RS10755'],
                        '\np-adj :', sig$padj[sig$locusId=='BBR_RS10755']))+
  theme_bw()

title=higherMedianFitnessColon%>%
  filter(locusId == 'BBR_RS16470')%>%
  .$kofamFunction

fitnessScores%>%
  pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')%>%
  filter(tissue != 'T0',
         locusId == 'BBR_RS16470')%>%
  ggplot(aes(x = tissue, y = fitnessScore))+
  geom_boxplot(outliers = F)+
  geom_point(aes(col = day))+
  labs(title = title,
       caption = paste0('Coef: ', sig$coef[sig$locusId=='BBR_RS16470'],
                        '\np-adj :', sig$padj[sig$locusId=='BBR_RS16470']))+
  theme_bw()

title='PTS sugar transporter subunit IIA'

fitnessScores%>%
  pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')%>%
  filter(tissue != 'T0',
         locusId == 'BBR_RS17590')%>%
  ggplot(aes(x = tissue, y = fitnessScore))+
  geom_boxplot(outliers = F)+
  geom_point(aes(col = day))+
  labs(title = title,
       caption = paste0('Coef: ', sig$coef[sig$locusId=='BBR_RS17590'],
                        '\np-adj :', sig$padj[sig$locusId=='BBR_RS17590']))+
  theme_bw()
write_tsv(output, file = 'linear models/lmTissueComparisionsAllTimePoints/linearModelRes/sigResFirstModel.tsv')