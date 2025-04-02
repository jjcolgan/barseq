library(tidyverse)
'Trying to implement a variance filter for the model to increase the significant results. '

modelRes = read_tsv('linear models/lmTissueComparisionsAllTimePoints/linearModelRes/sigResFirstModel.tsv')

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

fullData=fitnessScores%>%
  pivot_longer(cols = c(4: ncol(fitnessScores)),
               names_to = 'sample',
               values_to = 'fitnessScore')%>%
  left_join(metadata,
            by = 'sample')

geneVariance=fullData%>%
  filter(tissue != 'T0')%>%
  group_by(locusId)%>%
  summarise(geneVariance = var(fitnessScore))

'Might be able to chop out genes with a variance less than .25'
geneVariance%>%
  rename(gene = locusId)%>%
  left_join(modelRes,
            by = 'gene')%>%
  ggplot(aes(x = geneVariance,
             fill = padj < .1))+
  geom_histogram(binwidth = .125)

geneMeanFitness=fullData%>%
  filter(tissue != 'T0')%>%
  group_by(locusId)%>%
  summarise(meanFitnessScore = mean(fitnessScore))

'Can probably chop out genes with an average fitness score less than .125'
geneMeanFitness%>%
  rename(gene = locusId)%>%
  left_join(modelRes,
            by = 'gene')%>%
  ggplot(aes(x = abs(meanFitnessScore),
             fill = padj < .1))+
  geom_histogram(binwidth = .125)

'Lets test it out and see how things change'
'Variance .125 is not good, only raises the number of significant genes
to 89'

'var .25 and .3 only yield 98 sig genes'
"var . 5 yields 93 genes"
genesToTest = geneVariance %>%
  filter(geneVariance >= .125)

'This filters ~ 400 genes'
print(nrow(genesToTest))

lmIn=fitnessScores %>%
  select(-c(sysName, desc))%>%
  column_to_rownames('locusId')

loci = genesToTest%>%
  .$locusId

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

'Filter increases the number of significant genes from 84 to 98'

sig=output%>%
  filter(padj < .1)%>%
  arrange(padj)

output %>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Histogram of unadjusted pvalues')+
  theme_bw()

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
  labs(col = 'Significance',
       title = 'Variance filtered')

'mean greater than or equal 2 .25 yields 129 significant results'
'In my mind retaining genes with a mean fitness score greater than .5 makes sense as these are going to be what
has been dubbed a "significant" phenotype previously. This needs to be done in a group specific manner however.
I think I am not testing genes that should be examined. I want to combine this with a variance filter.'
genesToTest = geneMeanFitness %>%
  filter(abs(meanFitnessScore) > .5)

'This filters ~ 700 genes'
print(nrow(genesToTest))

lmIn=fitnessScores %>%
  select(-c(sysName, desc))%>%
  column_to_rownames('locusId')

loci = genesToTest%>%
   .$locusId

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

'Filter increases the number of significant genes from 84 to 102'
sig=output%>%
  filter(padj < .1)%>%
  arrange(padj)

output %>%
  ggplot(aes(x = pvalue))+
  geom_histogram(binwidth = .01)+
  labs(title = 'Histogram of unadjusted pvalues')+
  theme_bw()

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
  filter(locusId %in% genesToTest$locusId)%>%
  ggplot(aes(x = dj,
             col = col,
             y = colon))+
  geom_point(alpha = .15)+
  geom_abline()+
  theme_bw()+
  labs(col = 'Significance',
       title = 'Mean filtered results')

geneMeanFitnessByTissuePass=fullData%>%
  filter(tissue != 'T0')%>%
  group_by(locusId, tissue)%>%
  summarise(meanFitnessScore = mean(abs(fitnessScore)))%>%
  filter(meanFitnessScore > .5)%>%
  select(locusId)%>%
  distinct()

variancePass=geneVariance %>%
  filter(geneVariance >= .2)


loci = intersect(variancePass$locusId, geneMeanFitnessByTissuePass$locusId)

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
  filter(padj <= .1)%>%
  arrange(padj)

output %>%
  ggplot(aes(x = pvalue))+
  geom_histogram(binwidth = .01)+
  labs(title = 'Histogram of unadjusted pvalues')+
  theme_bw()

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
  labs(col = 'Significance',
       title = 'Mean filtered results and variance filtered results')

for (s in sig$gene){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(metadata, by = 'sample')%>%
    filter(tissue != 'T0',
           locusId == s)%>%
    ggplot(aes(x = tissue, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sig$coef[sig$gene==s],
                          '\np-adj :', sig$padj[sig$gene==s]))
  plot(p)
}




