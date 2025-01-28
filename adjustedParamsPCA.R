library(tidyverse)
library(vegan)
library(ape)
library(ComplexHeatmap)
setwd('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq')
strongFitnessEffects=read_tsv('barseqAdjustedParams/strong.tab')
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
  left_join(metadata) %>%
  filter(cor12 >= .25)

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

#make PCoA based on jaccard distance, by looking at genes under positive selection first
strongFitnessEffects$jaccardPositiveIn = 0
strongFitnessEffects$jaccardPositiveIn[strongFitnessEffects$lrn > 2] = 1

positiveJaccardIn=strongFitnessEffects%>%
  filter(sample %in% metadata$sample)%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from = jaccardPositiveIn)%>%
  column_to_rownames('sample')

#filter invariant genes and samples
positiveJaccardIn[is.na(positiveJaccardIn)]= 0
column_variances <- apply(positiveJaccardIn, 2, var, na.rm = TRUE)
positiveJaccardIn <- positiveJaccardIn[, column_variances > 0]
row_variances <- apply(positiveJaccardIn, 1, var, na.rm = TRUE)
positiveJaccardIn <- positiveJaccardIn[row_variances > 0, ]

jaccardPositive=vegdist(positiveJaccardIn, method = 'jaccard', binary = T)
jaccardPositiveOut=cmdscale(jaccardPositive, eig = T)
jaccardPositiveOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = tissue,
             y = V2))+

  geom_point()
jaccardPositiveOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = day,
             y = V2))+
  geom_point()

jaccardPositiveOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = day,
             shape = tissue,
             y = V2))+
  geom_point()
# literally removes all of the day ones, trying again with negative selection

#make PCoA based on jaccard distance, by looking at genes under negative selection first
strongFitnessEffects$jaccardNegative= 0
strongFitnessEffects$jaccardNegative[strongFitnessEffects$lrn <= -2] = 1

negativeJaccardIn=strongFitnessEffects%>%
  filter(sample %in% metadata$sample)%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from = jaccardNegative)%>%
  column_to_rownames('sample')

#filter invariant genes and samples
negativeJaccardIn[is.na(negativeJaccardIn)]= 0
column_variances <- apply(negativeJaccardIn, 2, var, na.rm = TRUE)
negativeJaccardIn <- negativeJaccardIn[, column_variances > 0]
row_variances <- apply(negativeJaccardIn, 1, var, na.rm = TRUE)
negativeJaccardIn <- negativeJaccardIn[row_variances > 0, ]
#filter low preveleance
colsums <- apply(negativeJaccardIn, 2, sum, na.rm = TRUE)
negativeJaccardIn <- negativeJaccardIn[, colsums >= 3]

jaccardNegative=vegdist(negativeJaccardIn, method = 'jaccard', binary = T)
jaccardNegativeOut=cmdscale(jaccardNegative, eig = T, k =4)
jaccardNegativeOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = tissue,
             y = V2))+
  
  geom_point()+
  stat_ellipse()
jaccardNegativeOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = day,
             y = V2))+
  geom_point()+
  stat_ellipse()

jaccardNegativeOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = day,
             shape = tissue,
             y = V2))+
  geom_point()

jaccardNegativeOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = tissueDay,
             y = V2))+
  geom_point()+
  stat_ellipse()

jaccardNegativeOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V2,
             col = day,
             shape = tissue,
             y = V3))+
  geom_point()

jaccardNegativeOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V3,
             col = day,
             shape = tissue,
             y = V4))+
  geom_point()

strongFitnessEffects$jaccardAbsolute= 0
strongFitnessEffects$jaccardAbsolute[abs(strongFitnessEffects$lrn)>= 2] = 1

absoluteJaccardIn=strongFitnessEffects%>%
  filter(sample %in% metadata$sample)%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from = jaccardAbsolute)%>%
  column_to_rownames('sample')

#filter invariant genes and samples
absoluteJaccardIn[is.na(absoluteJaccardIn)]= 0
column_variances <- apply(absoluteJaccardIn, 2, var, na.rm = TRUE)
absoluteJaccardIn <- absoluteJaccardIn[, column_variances > 0]
row_variances <- apply(absoluteJaccardIn, 1, var, na.rm = TRUE)
absoluteJaccardIn <- absoluteJaccardIn[row_variances > 0, ]

absoluteJaccardIn=vegdist(absoluteJaccardIn, method = 'jaccard', binary = T)
jaccardAbsoluteOut=cmdscale(absoluteJaccardIn, eig = T, k =4)
jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = tissue,
             y = V2))+
  
  geom_point()
jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = day,
             y = V2))+
  geom_point()

jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = day,
             shape = tissue,
             y = V2))+
  geom_point()

jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V2,
             col = day,
             shape = tissue,
             y = V3))+
  geom_point()


strongFitnessEffects$jaccardAbsolute[abs(strongFitnessEffects$lrn)>= 5] = 1

absoluteJaccardIn=strongFitnessEffects%>%
  filter(sample %in% metadata$sample)%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from = jaccardAbsolute)%>%
  column_to_rownames('sample')

#filter invariant genes and samples
absoluteJaccardIn[is.na(absoluteJaccardIn)]= 0
column_variances <- apply(absoluteJaccardIn, 2, var, na.rm = TRUE)
absoluteJaccardIn <- absoluteJaccardIn[, column_variances > 0]
row_variances <- apply(absoluteJaccardIn, 1, var, na.rm = TRUE)
absoluteJaccardIn <- absoluteJaccardIn[row_variances > 0, ]

#filter low detection genes, lets do sum greater than or equal to 3
colsums = apply(absoluteJaccardIn, 2, sum, na.rm = TRUE)
absoluteJaccardIn <- absoluteJaccardIn[, colsums >=3]

absoluteJaccardIn=vegdist(absoluteJaccardIn, method = 'jaccard', binary = T)
jaccardAbsoluteOut=cmdscale(absoluteJaccardIn, eig = T, k =4)
jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = tissue,
             y = V2))+
  
  geom_point()

jaccardAbsoluteOut%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = tissue,
             y = V2))+

  geom_point()
jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = day,
             y = V2))+
  geom_point()

jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = u,
             shape = tissue,
             y = V2))+
  geom_point()

jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = mad12,
             y = V2))+
  geom_point()

jaccardAbsoluteOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = cor12,
             y = V2))+
  geom_point()


heatmapIn=strongFitnessEffects%>%
  filter(sample %in% metadata$sample)%>%
  pivot_wider(names_from = tissueDayMouse,id_cols = locusId, values_from = lrn)%>%
  column_to_rownames('locusId')
heatmapIn[is.na(heatmapIn)]=0
Heatmap(heatmapIn, show_row_names = F)
