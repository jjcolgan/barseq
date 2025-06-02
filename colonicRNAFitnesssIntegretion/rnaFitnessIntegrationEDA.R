library(clusterProfiler)
library(KEGGREST)
library(ggplot2)
library(purrr)
library(progressr)
library(ComplexUpset)
library(org.Mm.eg.db)
library(tidyverse)

'Something weird is going on here, changing back to maaslin since the LMER model is basically saying everything
is significant'

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

spearmanRes = read_tsv('colonicRNAFitnesssIntegretion/spearmanRes/spearmanRes.tsv')
expression = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/sigGeneTab.tsv')
fitnessMeta = read_tsv('fullbarseqMeta.txt')
annotations = read_tsv('genesWithAnvioAnnotations.tsv')
allGenes = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')

keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction)%>%
  distinct()
metadataRNA = read.csv('metadataRnaSeq.csv')%>%
  as.data.frame()

sampleLibraryRna = metadataRNA %>%
  filter(Tissue == 'Colon'& Treatment == 'Bar-seq')%>%
  select(sample, library)
'Clean up column names'
expression=expression %>%
  rename_with(~ sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", .))

expression= expression %>%
  rename('2490Co'='2490Co-d')%>%
  column_to_rownames('SYMBOL')%>%
  t()%>%
  as.data.frame()

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

signficantSpearmanRes=spearmanRes%>%
  filter(padjust < .1)

'Generally more positive relationships than negative'
spearmanRes%>%
  ggplot(aes(x = Spearman_Correlation))+
  geom_histogram(bins=100)+
  labs(title = 'Histogram of all spearman values')

spearmanRes%>%
  filter(padjust < .1)%>%
  ggplot(aes(x = Spearman_Correlation))+
  geom_histogram(bins = 100)+
  labs(title = 'Histogram of signficant spearman values')

top20Genes=signficantSpearmanRes%>%
  group_by(Gene)%>%
  summarise('nSignficantMutants' = n())%>%
  arrange(desc(nSignficantMutants))%>%
  head(10)%>%
  .$Gene

top20GenesMutantBriteCategories=signficantSpearmanRes%>%
  filter(Gene %in% top20Genes)%>%
  select(Mutant)%>%
  distinct()%>%
  rename(locusId = Mutant)%>%
  left_join(annotations)%>%
  select(locusId, kofamAccession)%>%
  distinct()%>%
  annotate_modules()

numberOfMutantAssoicationsTop20genes=signficantSpearmanRes%>%
  group_by(Gene)%>%
  summarise('nSignficantMutants' = n())%>%
  ungroup()

p=signficantSpearmanRes%>%
  filter(#Gene %in% top20Genes,
    padjust < .1)%>%
  rename(locusId = Mutant)%>%
  left_join(top20GenesMutantBriteCategories, by = 'locusId')%>%
  left_join(numberOfMutantAssoicationsTop20genes,
            by = "Gene")%>%
  ggplot(aes(x = fct_reorder(Gene, -nSignficantMutants),
             fill = moduleLevel1))+
  geom_bar()

ggsave(filename = '/Users/johnjamescolgan/Desktop/testPlot.pdf', p)

p=signficantSpearmanRes%>%
  filter(#Gene %in% top20Genes,
    padjust < .1)%>%
  rename(locusId = Mutant)%>%
  left_join(top20GenesMutantBriteCategories, by = 'locusId')%>%
  filter(moduleLevel1 == ' 09100 Metabolism')%>%
  left_join(numberOfMutantAssoicationsTop20genes,
            by = "Gene")%>%
  ggplot(aes(x = fct_reorder(Gene, -nSignficantMutants),
             fill = moduleLevel2))+
  geom_bar()


ggsave(filename = '/Users/johnjamescolgan/Desktop/testPlotMetabolismOnly.pdf', p)

signficantSpearmanRes%>%
  group_by(Gene)%>%
  summarise('nSignficantMutants' = n())%>%
  ggplot(aes(x = nSignficantMutants,
             fill = Gene %in% top20Genes))+
  geom_histogram()

signficantSpearmanRes%>%
  group_by(Gene)%>%
  summarise('nSignficantMutants' = n())%>%
  arrange(desc(nSignficantMutants))%>%
  head(20)%>%
  write_tsv('colonicRNAFitnesssIntegretion/top20GenesWithMostAssoications.tsv')


signficantSpearmanRes%>%
  group_by(Mutant)%>%
  summarise('nSignficantGenes' = n())%>%
  ggplot(aes(x = nSignficantGenes))+
  geom_histogram()

signficantSpearmanRes%>%
  group_by(Mutant)%>%
  summarise('nSignficantGenes' = n())%>%
  arrange(desc(nSignficantGenes))%>%
  head(10)%>%
  rename(locusId = Mutant)%>%
  left_join(annotations,
            by = 'locusId')%>%
  select(locusId, nSignficantGenes, desc, kofamFunction)%>%
  distinct()%>%
  write_tsv('colonicRNAFitnesssIntegretion/top20MutantsWithMostAssoications.tsv')

signficantSpearmanRes%>%
  group_by(Mutant)%>%
  summarise('nSignficantGenes' = n())%>%
  arrange(nSignficantGenes)

mutantsTotest = signficantSpearmanRes%>%
  group_by(Mutant)%>%
  summarise('nSignficantGenes' = n())%>%
  arrange(desc(nSignficantGenes))%>%
  head(10)%>%
  .$Mutant

output = data.frame('ID' = character(),
                    'Description' = character(),
                    'GeneRatio' = numeric(),
                    'BgRatio' = numeric(),
                    'RichFactor' = numeric(),
                    'FoldEnrichment' = numeric(),
                    'zScore' = numeric(),
                    'pvalue' = numeric(),
                    'p.adjust' = numeric(),
                    'qvalue' = numeric(),
                    'geneID' = character(),
                    'Count' = numeric(),
                    'mutant' = character())

for (m in mutantsTotest){
  genes=signficantSpearmanRes%>%
    filter(Mutant == m)%>%
    .$Gene
  goRes=enrichGO(gene = genes,
                 ont = 'BP',
                 keyType = 'SYMBOL',
                 universe = allGenes$SYMBOL,
                 pvalueCutoff = .05,
                 org.Mm.eg.db,
                 minGSSize = 25,
                 maxGSSize = 300)%>%
    clusterProfiler::simplify()%>%
    as.data.frame()
  if(nrow(goRes)< 1){
    goRes = data.frame('ID' = 'none',
                       'Description' = 'none',
                       'GeneRatio' = 0,
                       'BgRatio' = 0,
                       'RichFactor' = 0,
                       'FoldEnrichment' = 0,
                       'zScore' = 0,
                       'pvalue' = 0,
                       'p.adjust' = 0,
                       'qvalue' = 0,
                       'geneID' = '',
                       'Count' = 0)
  }
  goRes$mutant = m
  output = rbind(goRes, output)
}

output %>%
  filter(Description != 'none',
         Count > 15) %>%
  ggplot(aes(y = fct_infreq(Description))) +
  geom_bar() +
  labs(y = "GO term", x = "Number of mutant assoications") +
  theme_minimal()+
  labs(title = 'GO term assoicated with the significant mutants in the colon background all tested genes')

output = data.frame('ID' = character(),
                    'Description' = character(),
                    'GeneRatio' = numeric(),
                    'BgRatio' = numeric(),
                    'RichFactor' = numeric(),
                    'FoldEnrichment' = numeric(),
                    'zScore' = numeric(),
                    'pvalue' = numeric(),
                    'p.adjust' = numeric(),
                    'qvalue' = numeric(),
                    'geneID' = character(),
                    'Count' = numeric(),
                    'mutant' = character())

for (m in mutantsTotest){
  genes=signficantSpearmanRes%>%
    filter(Mutant == m)%>%
    .$Gene
  goRes=enrichGO(gene = genes,
                 ont = 'BP',
                 keyType = 'SYMBOL',
                 universe = unique(spearmanRes$Gene),
                 pvalueCutoff = .05,
                 org.Mm.eg.db,
                 minGSSize = 25,
                 maxGSSize = 300)%>%
    clusterProfiler::simplify()%>%
    as.data.frame()
  if(nrow(goRes)< 1){
    goRes = data.frame('ID' = 'none',
                       'Description' = 'none',
                       'GeneRatio' = 0,
                       'BgRatio' = 0,
                       'RichFactor' = 0,
                       'FoldEnrichment' = 0,
                       'zScore' = 0,
                       'pvalue' = 0,
                       'p.adjust' = 0,
                       'qvalue' = 0,
                       'geneID' = '',
                       'Count' = 0)
  }
  goRes$mutant = m
  output = rbind(goRes, output)
}

output %>%
  filter(Description != 'none',
         Count > 15) %>%
  ggplot(aes(y = fct_infreq(Description))) +
  geom_bar() +
  labs(y = "GO term", x = "Number of mutant assoications") +
  theme_minimal()+
  labs(title = 'GO term assoicated with the top 10 mutants by number of gene associations background degs')

plotIn = spearmanRes%>%
  arrange(padjust)%>%
  head(20)

plotIn=expression%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('Gene')%>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = 'sample',
               values_to = 'vsdCounts'
  )%>%
  merge(plotIn, by = 'Gene')

plotIn=fitnessScores%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('Mutant')%>%
  pivot_longer(cols = c(2:ncol(.)),
               names_to = 'sample',
               values_to = 'FitnessScores'
  )%>%
  merge(plotIn,
        by = c('Mutant',
               'sample'))

plotIn  %>%
  rename(mergeCol = sample)%>%
  left_join(fitnessMeta, by = "mergeCol")%>%
  ggplot(aes(x = vsdCounts,
             col = day,
             group = tissue,
             y = FitnessScores))+
  geom_smooth(method = 'lm')+
  geom_point(alpha = .5)+
  facet_wrap(~Mutant+Gene,
             scales = 'free_y',
             nrow = 4,)

plotIn  %>%
  rename(mergeCol = sample)%>%
  left_join(fitnessMeta, by = "mergeCol")%>%
  ggplot(aes(y = vsdCounts,
             col = day,
             group = tissue,
             x = FitnessScores))+
  geom_smooth(method = 'lm')+
  geom_point(alpha = .5)+
  facet_wrap(~Mutant+Gene,
             scales = 'free_y',
             nrow = 4,)
# p=cbind(expressionTemp, fitnesstemp) %>%
#   rownames_to_column('mergeCol')%>%
#   left_join(fitnessMeta, by = 'mergeCol')%>%
#   ggplot(aes(x = (.data[[gene]]),
#              y = .data[[mutant]],
#              group = tissue,
#              col = day))+
#   geom_point() +
#   geom_smooth(method = 'lm')+
#   labs(x = paste0('log2', gene), y = mutant)
# plot(p)


