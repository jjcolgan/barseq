fullLocusFitness<-read_tsv('full/fit_t.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$numericDay = 0
metadata$numericDay[metadata$day== 'day 14'] = 14
metadata$numericDay = as.integer(metadata$numericDay)
metadata=metadata%>%
  mutate(mouseDay = paste(mouse, day, sep = '-'))
metadata=metadata%>%
  mutate(mouseDayTissue = paste(mouse, day, tissue,sep = '-'))
day1Meta=metadata %>%
  filter(day == 'day 14' | day == 't0')

genes=read_tsv('genesWithAnvioAnnotations.tsv')

fullLocusFitness<-read_tsv('full/fit_t.tab')

fullLocusFitnessTransposed<-fullLocusFitness %>%
  pivot_longer(cols = 4:45, values_to = 'fitnessScore', names_to = 'sample') %>%
  left_join(day1Meta, by ='sample')%>%
  filter(sample %in% day1Meta$sample)

fullLocusFitnessScores<-fullLocusFitnessTransposed %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()

summary(fullLocusFitnessScores)
adonisData<-fullLocusFitnessTransposed %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore)
adonisMetadata<-metadata %>%
  filter(sample %in%adonisData$sample)
adonisData = adonisData %>%
  column_to_rownames('sample')
adonis2(adonisData~tissue, data= adonisMetadata, method = 'euc')

fullLocusFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = tissue,shape = day, label = sample))+
  geom_point()+
  labs(title = 'Gene fitness PCA day 14 with T0',
       x = 'PC1 68.65% var',
       caption = 'p = .927, PERMANOVA\nR^2 = .005873',
       y = 'PC2 10.52% var')

fullLocusFitnessTransposedNoT0<-fullLocusFitnessTransposed %>%
  filter(day == 'day 14')

fullLocusFitnessScoresNoT0<-fullLocusFitnessTransposedNoT0 %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')%>%
  prcomp()

summary(fullLocusFitnessScoresNoT0)
adonisData<-fullLocusFitnessTransposedNoT0 %>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore)
adonisMetadata<-metadata %>%
  filter(sample %in%adonisData$sample)
adonisData = adonisData %>%
  column_to_rownames('sample')
adonis2(adonisData~tissue, data= adonisMetadata, method = 'euc')

fullLocusFitnessScoresNoT0$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1, y = PC2, col = tissue,shape = day, label = sample))+
  geom_point()+
  labs(title = 'Gene fitness PCA day 14 no t0',
       x = 'PC1 57.36% var',
       caption = 'p < .01, PERMANOVA\nR^2 = 0.4651',
       y = 'PC2 14.93% var')


fullLocusFitnessTransposedNoT0%>%
  pivot_wider(names_from = mouseDayTissue,id_cols = locusId, values_from =  fitnessScore) %>%
  column_to_rownames('locusId')%>%
  scale()%>%
  as.matrix()%>%
  Heatmap(show_row_names = F)


meanFitnessScores=fullLocusFitnessTransposedNoT0 %>%
  group_by(tissue, locusId)%>%
  summarise('meanFitnessScore' = mean(fitnessScore))

meanFitnessScores %>%
  pivot_wider(id_cols = 'locusId',
              names_from = 'tissue',
              values_from = 'meanFitnessScore') %>%
  ggplot(aes(x = dj,
             y = colon))+
  geom_point()+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  geom_vline(xintercept = 1.5)+
  geom_vline(xintercept = -1.5)


maaslinData = fullLocusFitnessTransposedNoT0%>%
  pivot_wider(names_from = locusId,id_cols = sample, values_from =  fitnessScore) %>%
  column_to_rownames('sample')
maaslinMeta  = adonisMetadata %>%
  column_to_rownames('sample')

fit_data = Maaslin2(
  input_data = maaslinData,
  input_metadata = maaslinMeta,
  output = "day14DjVcolonGeneLocus",
  normalization = 'none',
  transform = 'none',
  min_abundance = 0,
  min_prevalence = 0,
  fixed_effects = c("tissue"))

maaslin2Res = read_tsv(file = 'day14DjVcolonGeneLocus/all_results.tsv')

maaslin2Res=maaslin2Res %>%
  rename(locusId = 'feature') %>%
  left_join(genes, by = 'locusId')

enrichedFitnessColonSig<-maaslin2Res %>%
  filter(qval < .05,
         coef < 0)

enrichedFitnessDjSig<-maaslin2Res %>%
  filter(qval < .05,
         coef > 0)