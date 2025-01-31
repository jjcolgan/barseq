library(vegan)
library(tidyverse)

fitness = read_tsv('barseqAdjustedParams/fit_logratios.tab')
colnames(fitness) <- sub("setA", "", colnames(fitness))
colnames(fitness) <- sub("_.*", "", colnames(fitness))
colnames(fitness) <- sub("CO$", "Co", colnames(fitness))
colnames(fitness) <- sub("DJ$", "Dj", colnames(fitness))

metadata = read_tsv('fullbarseqMeta.txt')

noTo=metadata %>%
  filter(tissue != 'T0')%>%
  .$sample

pcoaIn=fitness %>%
  select(-c(desc,
            sysName))%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  vegdist(method = 'euc')%>%
  cmdscale(k =4, eig = T)

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V1,
             y = V2))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V1,
             col = tissue,
             y = V2))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V1,
             col = day,
             y = V2))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V2,
             y = V3))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V2,
             col = tissue,
             y = V3))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V2,
             col = day,
             y = V3))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V3,
             y = V4))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V3,
             col = tissue,
             y = V4))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V3,
             col = day,
             y = V4))+
  geom_point()

pcoaIn=fitness %>%
  select(-c(desc,
            sysName))%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% noTo)%>%
  column_to_rownames('sample')%>%
  vegdist(method = 'euc')%>%
  cmdscale(k =4, eig = T)

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V1,
             y = V2))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V1,
             col = tissue,
             y = V2))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V1,
             col = day,
             y = V2))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V2,
             y = V3))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V2,
             col = tissue,
             y = V3))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V2,
             col = day,
             y = V3))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V3,
             y = V4))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V3,
             col = tissue,
             y = V4))+
  geom_point()

pcoaIn$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by ='sample')%>%
  ggplot(aes(x = V3,
             col = day,
             y = V4))+
  geom_point()

eucDistMat=fitness %>%
  select(-c(desc,
            sysName))%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  vegdist(method = 'euc')%>%
  as.matrix()%>%
  as.data.frame()

eucDistLong=eucDistMat%>%
  rownames_to_column('sample1')%>%
  pivot_longer(cols = 2:43, names_to = 'sample2', values_to = 'eucDist')

eucDistLong=eucDistLong%>%
  rename('sample'= sample1)%>%
  merge(select(metadata, c('sample', 'tissueDay','mouse', 'tissue')), by = 'sample')%>%
  rename('sample1'=sample,
         'tissueDaySample1' = tissueDay,
         'tissue1' = tissue,
         'sample1Mouse' = mouse)%>%
  rename('sample'= sample2)%>%
  merge(select(metadata, c('sample', 'tissueDay','mouse','tissue')), by = 'sample')%>%
  rename('sample2'=sample,
         'tissueDaySample2' = tissueDay,
         'tissue2' = tissue,
         'sample2Mouse' = mouse)

eucDistLong=eucDistLong %>%
  mutate(tissueComparison = paste(tissue1, tissue2, sep = '-'))

eucDistLongNoDups <- eucDistLong %>%
  # Create a new column with the sorted pair of sample1 and sample2
  mutate(sorted_pair = paste(pmin(sample1, sample2), pmax(sample1, sample2), sep = "_")) %>%
  # Remove duplicate rows based on the sorted_pair column
  distinct(sorted_pair, .keep_all = TRUE) %>%
  # Remove the temporary sorted_pair column
  select(-sorted_pair)
eucDistLongNoDups[eucDistLongNoDups$tissueComparison == 'dj-colon', ]$tissueComparison = 'colon-dj'
eucDistLongNoDups%>%
  filter(tissue1 != 'T0',
         tissue2 != 'T0',
         sample1 != sample2)%>%
  ggplot(aes(x = tissueComparison, y = eucDist))+
  geom_violin()+
  ggpubr::stat_compare_means(comparisons = list(c('colon-colon','colon-dj'),
                                                c('colon-colon','dj-dj'),
                                                c('colon-dj','dj-dj')))

eucDistLongNoDups=eucDistLongNoDups %>%
  mutate(tissueDayComparison = paste(tissueDaySample1, tissueDaySample2, sep = '-'))
eucDistLongNoDups$tissueDayComparison <- sapply(strsplit(eucDistLongNoDups$tissueDayComparison, "-"), function(x) paste(sort(x), collapse = "-"))
eucDistLongNoDups$tissueDayComparison%>%
  unique()

eucDistLongNoDups %>%
  filter(tissue1 == 'dj',
         tissue2 == 'dj',
         sample1 != sample2) %>%
  ggplot(aes(x = reorder(tissueDayComparison, eucDist, FUN = median),
             y = eucDist)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0, 100)

eucDistLongNoDups %>%
  filter(tissue1 == 'colon',
         tissue2 == 'colon',
         sample1 != sample2) %>%
  ggplot(aes(x = reorder(tissueDayComparison, eucDist, FUN = median),
             y = eucDist)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0, 100)
