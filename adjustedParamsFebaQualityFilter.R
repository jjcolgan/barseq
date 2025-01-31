library(tidyverse)

fitness = read_tsv('barseqAdjustedParams/fit_logratios.tab')
colnames(fitness) <- sub("setA", "", colnames(fitness))
colnames(fitness) <- sub("_.*", "", colnames(fitness))
colnames(fitness) <- sub("CO$", "Co", colnames(fitness))
colnames(fitness) <- sub("DJ$", "Dj", colnames(fitness))

metadata = read_tsv('fullbarseqMeta.txt')

quality=read_tsv('barseqAdjustedParams/fit_quality.tab')
quality$name <- sub("setA", "", quality$name)
quality$name <- sub("_.*", "", quality$name)
quality$name <- sub("CO$", "Co", quality$name)
quality$name <- sub("DJ$", "Dj", quality$name)
quality = quality %>%
  rename('sample'= name)


metadata = metadata %>%
  left_join(quality, by = 'sample')%>%
  filter(u == T)

pcaOut=fitness %>%
  select(-c(desc,
            sysName))%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% metadata$sample)%>%
  column_to_rownames('sample')%>%
  prcomp(center = TRUE)

summary(pcaOut)

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissue))+
  geom_point()+
  labs(title = 'Default filter',
       x = 'PC1 50.01%',
       y = 'PC2 17.37%')

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = day))+
  geom_point()+
  labs(title = 'Default filter',
       x = 'PC1 50.01%',
       y = 'PC2 17.37%')

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = tissue,
             col = day))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = tissue,
             group = mouse,
             col = mouse))+
  geom_point()+
  geom_line()+
  theme(legend.position = 'none')

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = tissue,
             group = mouse,
             col = day))+
  geom_point()+
  geom_line()+
  theme(legend.position = 'none')

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = log10(millionBases)))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = meanQualityScore))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = gccor))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = adjcor))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = mad12))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = cor12))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = tissue))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = day))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             shape = tissue,
             col = day))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             shape = tissue,
             group = mouse,
             col = mouse))+
  geom_point()+
  geom_line()+
  theme(legend.position = 'none')

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             shape = tissue,
             group = mouse,
             col = day))+
  geom_point()+
  geom_line()+
  theme(legend.position = 'none')

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = log10(millionBases)))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = meanQualityScore))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = gccor))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = adjcor))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = mad12))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = cor12))+
  geom_point()