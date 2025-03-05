library(ComplexHeatmap)
library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyverse)
library(gaston)
#first try
'Initally ran on the output of DeSeq2 clusetering containing all samples (conventionalized and GF)
and all strong phenotypes identified by FEBA. This resulted in well over a million pairwise comparisons,
with only 5 passing FDR correction using estimated p-values. Exact pvalues from cor.test, has more but many of
the p-values are 0, which likely represents a rounding error. Going to retry, using sig degs from just barseq comparisons.

Using just the genes from the barseq comparisons reduces the number of significant genes from ~3,000 to 355. DEGs need to have L2FC greater than 1
in this case. Reducing that threshold to .5 brings the number to about 1000.  32 genes significant below FDR .1. Going to add a filter for mutants
must be present in at least 2 samples in a group to be considered. This leaves about 16 mutants with only 16 significant differences.

Going back to the input set of 355 genes. Even with the reduction in the number of mutants there are 256,000 tests this yields 91,235 tests and 9 significant correlations. Straight up
does not work. '

simpleQQPlot = function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)),
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
}
spearmanRes=read_tsv('sigColonicBarseqGenesStrongMutantSpearmanCorrelationEstimatedPvalues.tsv')
fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')

fitnessScores %>%
  colnames()

metadataRNA = read.csv('metadataRnaSeq.csv')%>%
  as.data.frame()

sampleLibraryRna = metadataRNA %>%
  filter(Tissue == 'Colon'& Treatment == 'Bar-seq')%>%
  filter(sample != '2493Co')%>%
  select(sample, library)

expression = read_tsv('sigGenesColonicBarseq.tsv')
expression=expression %>%
  rename_with(~ sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", .))

expression=expression %>%
  column_to_rownames('gene')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('library')%>%
  merge(sampleLibraryRna, by = 'library') %>%
  select(-library)%>%
  column_to_rownames('sample')
annotatedMutants = read_tsv('genesWithAnvioAnnotations.tsv')

fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
colnames(fitnessScores) <- sub("setA", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("_.*", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("CO$", "Co", colnames(fitnessScores))
colnames(fitnessScores) <- sub("DJ$", "Dj", colnames(fitnessScores))
fitnessScores=fitnessScores %>%
  select(-c(sysName, desc))%>%
  as.data.frame()%>%
  column_to_rownames('locusId')%>%
  t()

fitnessMeta = read_tsv('fullbarseqMeta.txt')
fitnessMeta=fitnessMeta%>%
  filter(tissue == 'colon')
fitnessMeta=fitnessMeta %>%
  mutate(mergeCol = paste0(mouse, 'Co'))

fitnessMeta = fitnessMeta %>%
  distinct(mergeCol, .keep_all = TRUE)


expression= expression[rownames(expression) %in% fitnessMeta$mergeCol, , drop = FALSE]

fitnessMeta=fitnessMeta %>%
  filter(mergeCol %in% row.names(expression))
fitnessScores=fitnessScores %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(fitnessMeta[c('sample', 'mergeCol')], by = 'sample')

expression= expression[rownames(expression) %in% fitnessScores$mergeCol, , drop = FALSE]
fitnessScores = fitnessScores %>%
  select(-sample)%>%
  column_to_rownames('mergeCol')


fitnessScores <- fitnessScores[order(rownames(fitnessScores)), , drop = FALSE]
expression <- expression[order(rownames(expression)), , drop = FALSE]



expressionLong=expression %>%
  rownames_to_column('Sample')%>%
  pivot_longer(cols = c(2:356), names_to = 'gene', values_to = 'rnaCounts')

fitnessLong=fitnessScores %>%
  rownames_to_column('Sample')%>%
  pivot_longer(cols = c(2:1676), names_to = 'mutant', values_to = 'fitness')

fulldata=expressionLong %>%
  rename(mergeCol = Sample)%>%
  left_join(fitnessMeta, by ="mergeCol")%>%
  rename(Sample = mergeCol)%>%
  left_join(fitnessLong, by = 'Sample')


spearmanRes$bh = p.adjust(spearmanRes$P_Value, method = 'BH')
sigRes=spearmanRes%>%
  filter(bh < .1)

spearmanRes$P_Value %>%
   hist()

for (i in 1:nrow(sigRes)){
  input = sigRes[i,]
  geneOfInterest = input$Gene[1]
  mutantOfInterest = input$Mutant[1]
  p=fulldata %>%
    filter(mutant == mutantOfInterest,
           gene == geneOfInterest)%>%
    ggplot(aes(x = rnaCounts, group = gene,col = day, y = fitness))+
    geom_point()+
    geom_smooth(method = 'lm')+
    labs( x = geneOfInterest,
          y =mutantOfInterest)
  plot(p)



}

geneSymbols=bitr(spearmanRes$Gene, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = org.Mm.eg.db)

spearmanRes=geneSymbols %>%
rename(Gene = ENSEMBL)%>%
  left_join(spearmanRes, by = 'Gene')
spearmanRes%>%
filter(Gene %in% sigRes$Gene,
       Mutant %in% sigRes$Mutant)%>%
  pivot_wider(names_from = 'Mutant', id_cols = 'Gene',
              values_from = 'Spearman_Correlation', )%>%
  column_to_rownames("Gene")%>%
  Heatmap(show_column_names = F, show_row_names = F)


spearmanRes %>%
  rename(locusId = Mutant)%>%
  left_join(annotatedMutants, by = 'locusId')%>%
  filter(padjust < .1)%>%
  select(c(locusId,Gene,SYMBOL, desc, kofamAccession,kofamFunction,padjust, Spearman_Correlation))%>%
  distinct()%>%
  view()
qqplot.pvalues(p = spearmanRes$P_Value)


