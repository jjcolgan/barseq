library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(KEGGREST)
library(viridis)
library(hrbrthemes)
setwd('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq')
strongFitnessEffects=read_tsv('barseqAdjustedParams/strong.tab')
strongFitnessEffects %>%
  select(locusId)%>%
  distinct()%>%
  nrow()
metadata = read_tsv('fullbarseqMeta.txt')
metadata = metadata %>%
  mutate('tissueDayMouse'=paste(tissueDay, mouse, sep = '-'))
quality=read_tsv('barseqAdjustedParams/fit_quality.tab')
quality$name <- sub("setA", "", quality$name)
quality$name <- sub("_.*", "", quality$name)
quality$name <- sub("CO$", "Co", quality$name)
quality$name <- sub("DJ$", "Dj", quality$name)
metadata = quality %>%
  rename('sample'= name)%>%
  left_join(metadata)

#metadata = metadata %>%
# filter(cor12 >= .2)

# clean up sample names
strongFitnessEffects$name <- sub("setA", "", strongFitnessEffects$name)
strongFitnessEffects$name <- sub("_.*", "", strongFitnessEffects$name)
strongFitnessEffects$name <- sub("CO$", "Co", strongFitnessEffects$name)
strongFitnessEffects$name <- sub("DJ$", "Dj", strongFitnessEffects$name)
strongFitnessEffects%>%
  select(name)%>%
  distinct()%>%
  nrow()

strongFitnessEffects=strongFitnessEffects %>%
  rename('sample'= name)

strongFitnessEffects=strongFitnessEffects%>%
  left_join(metadata, by = 'sample')

strongFitnessEffects%>%
  select(sample)%>%
  distinct()%>%
  nrow()

col_fun = colorRamp2(c(-10, 0, 5), c("blue", "white", "red"))
heatmapIn=strongFitnessEffects%>%
  filter(sample %in% metadata$sample)%>%
  pivot_wider(names_from = tissueDayMouse,id_cols = locusId, values_from = lrn)%>%
  column_to_rownames('locusId')
heatmapIn[is.na(heatmapIn)]=0
heatmapIn <- heatmapIn[rowSums(heatmapIn != 0) >= 3, ]
Heatmap(heatmapIn, show_row_names = F, col = col_fun, name = 'fitness scores')

functions = read_tsv('genesWithAnvioAnnotations.tsv')

sigGenes =strongFitnessEffects %>%
  left_join(functions, by = 'locusId')%>%
  filter(is.na(kofamAccession)!= T)%>%
  select(kofamAccession)%>%
  distinct()


sigGenes
sigGenes$brite <-NA
for (g in 1:nrow(sigGenes)){
  temp <- sigGenes %>%
    filter(kofamAccession == sigGenes$kofamAccession[g])
  keggout<- keggGet(temp$kofamAccession)
  b<-keggout[[1]]$BRITE[3]
  sigGenes[g, "brite"] <- as.character(b)
}

strongFitnessEffectsWithFunctions=functions %>%
  left_join(sigGenes, by = 'kofamAccession')%>%
  merge(strongFitnessEffects, by = 'locusId')

strongFitnessEffectsWithFunctions$selection = 'Negative'
strongFitnessEffectsWithFunctions[strongFitnessEffectsWithFunctions$lrn >0,]$selection = 'Positive'
strongFitnessEffectsWithFunctions[is.na(strongFitnessEffectsWithFunctions$kofamAccession)== T, ]$kofamAccession = 'Unannotated'
strongFitnessEffectsWithFunctions[is.na(strongFitnessEffectsWithFunctions$brite)== T, ]$brite = 'Unannotated'
strongFitnessEffectsWithFunctions[strongFitnessEffectsWithFunctions$brite== 'Unannotated', ]$brite = NA
strongFitnessEffectsWithFunctions$brite <- gsub("^\\s*\\d+\\s*", "", strongFitnessEffectsWithFunctions$brite)
strongFitnessEffectsWithFunctions$brite <- gsub("Protein families:", "Unclassified:", strongFitnessEffectsWithFunctions$brite)


strongFitnessEffectsWithFunctions$day = factor(strongFitnessEffectsWithFunctions$day,
                                               levels = c('day1', 'day3','day7','day14'))

# this is hard to read and IDK how i feel about it
strongFitnessEffectsWithFunctions %>%
  filter(selection == 'Negative')%>%
  ggplot(aes(x = tissue, fill = brite)) +
  geom_bar(position = "fill") +  # Normalize bar heights to 1
  facet_wrap(~day, nrow = 2) +
  ylab("Proportion")+theme_classic()