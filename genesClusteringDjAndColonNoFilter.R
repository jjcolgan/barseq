library(tidyverse)
library(ComplexHeatmap)

fullgenesFit = read_tsv('fullGenesTab.tsv')

fitnessScoresMat=fullgenesFit %>%
  select(sample, locusId, fitnessScore) %>%
  distinct()%>%
  pivot_wider(names_from = 'mouseDayTissue',id_cols = 'locusId', values_from = 'fitnessScore')

functions =fullgenesFit %>%
  select(locusId, name, desc,kofamAccession,kofamAccession, pfamAcession, pfamFunction) %>%
  distinct()

fitnessScoresMat%>%
  column_to_rownames('locusId')%>%
  Heatmap(show_row_names = F)

wss <- (nrow(fitnessScoresMat) - 1) * sum(apply(column_to_rownames(fitnessScoresMat,'locusId'), 2, var))
# Loop through cluster numbers from 2 to 20 to calculate WSS for each k-means solution
for (i in 2:20) {
  # Run k-means clustering with 'i' clusters and store the WSS for the solution
  wss[i] <- sum(kmeans(column_to_rownames(fitnessScoresMat,'locusId'),
                       centers = i)$withinss)
}
plot(1:20, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")

clusters=fitnessScoresMat%>%
  column_to_rownames('locusId')%>%
  kmeans(centers = 5)

clusterMembership=clusters$cluster%>%
  as.data.frame()

colnames(clusterMembership) = c('cluster')

clusterMembership=clusterMembership %>%
  rownames_to_column('locusId')

clusterMembership %>%
  group_by(cluster)%>%
  summarise('clusterSize'= n())

cluster1=clusterMembership %>%
  filter(cluster == 1) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster1$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  labs(title = 'Cluster 1 - Strong positive impact')

cluster2=clusterMembership %>%
  filter(cluster == 2) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster2$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  theme(legend.position = 'none')+
  labs(title = 'cluster 2 - weak fitness impacts')

cluster3=clusterMembership %>%
  filter(cluster == 3) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster3$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  theme(legend.position = 'none')+
  labs(title = 'Cluster 3 - Strong positive early, weakly negative late')

cluster4=clusterMembership %>%
  filter(cluster == 4) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster4$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  theme(legend.position = 'none')+
  labs('Cluster 4 - strong negative colonic')

cluster5=clusterMembership %>%
  filter(cluster == 5) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster5$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  theme(legend.position = 'none')+
  labs(title = 'cluster 5 very strong negative colonic')