library(tidyverse)
library(Maaslin2)
library(ggrepel)
fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$day = factor(metadata$day, levels = c('t0', 'day1', 'day3', 'day7', 'day14'))
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

metadata$phase = NA
metadata$phase[metadata$dayNumeric < 7] ='early'
metadata$phase[metadata$dayNumeric >= 7] ='late'

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

maaslinIn=fitnessScores %>%
  select(-c(sysName, desc))%>%
  column_to_rownames('locusId')
# variance filter was primarily added to only focus on samples where there
# actually a dramatic change in abundance and reduce the number of tests
# performed
# Maaslin2(input = maaslinIn,
#          input_metadata = colonMeta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = .5,
#          fixed_effects = c('day','pfClusters', 'millionBases'),
#          random_effects = ,
#          reference = 'day,day1;day3;day7;day14',
#          output = 'colonLogRatiosMaaslin2'
#
#
#
# )
sigResMaaslin=read_tsv('colonLogRatiosMaaslin2/significant_results.tsv')

sigResMaaslin %>%
  filter(metadata == 'day')%>%
  filter(qval < .05)%>%
  select(feature)%>%
  distinct()%>%
  nrow()

# Maaslin2(input = maaslinIn,
#          input_metadata = colonMeta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = .5,
#          fixed_effects = c('dayNumeric','pfClusters', 'millionBases'),
#          random_effects = ,
#          reference = '',
#          output = 'colonLogRatiosMaaslin2NumericDay'
# )

sigResMaaslinNumericDayColon=read_tsv('colonLogRatiosMaaslin2NumericDay/significant_results.tsv')
resMaaslinNumericDayColon=read_tsv('colonLogRatiosMaaslin2NumericDay/all_results.tsv')
sigResMaaslinNumericDayColon %>%
  filter(metadata == 'dayNumeric')%>%
  filter(qval < .05)%>%
  select(feature)%>%
  distinct()%>%
  nrow()

resMaaslinNumericDayColon %>%
  filter(metadata == 'dayNumeric')%>%
  ggplot(aes(x = coef,
             y= -log10(qval)))+
  geom_point()

highNegativeCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantPositiveCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

label_df <- bind_rows(
  highSignificantNegativeCoefMutants,
  highSignificantPositiveCoefMutants,
  highNegativeCoefMutants,
  highPositiveCoefMutants
) %>%
  select(locusId, kofamFunction) %>%
  distinct() %>%
  mutate(kofamFunction = replace_na(kofamFunction, "Unannotated"))

# Initialize label column to NA
resMaaslinNumericDayColon$label <- NA

# Join annotations into resMaaslinNumericDayColon
resMaaslinNumericDayColon <- resMaaslinNumericDayColon %>%
  left_join(label_df, by = c("feature" = "locusId")) %>%
  mutate(label = if_else(!is.na(kofamFunction), kofamFunction, label)) %>%
  select(-kofamFunction)

resMaaslinNumericDayColon %>%
  filter(metadata == 'dayNumeric')%>%
  ggplot(aes(x = coef,
             label = label,
             y= -log10(qval)))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 2.)




# Maaslin2(input = maaslinIn,
#          input_metadata = djMeta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = .5,
#          fixed_effects = c('day','pfClusters', 'millionBases'),
#          random_effects = ,
#          output = 'djLogRatiosMaaslin2'
# )
# Maaslin2(input = maaslinIn,
#          input_metadata = djMeta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = .5,
#          fixed_effects = c('dayNumeric','pfClusters', 'millionBases'),
#          random_effects = ,
#          reference = '',
#          output = 'djLogRatiosMaaslin2NumericDay'
# )

sigResMaaslinNumericDayDj=read_tsv('djLogRatiosMaaslin2NumericDay/significant_results.tsv')
resMaaslinNumericDayDj=read_tsv('djLogRatiosMaaslin2NumericDay/all_results.tsv')
highNegativeCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantPositiveCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

label_df <- bind_rows(
  highSignificantNegativeCoefMutants,
  highSignificantPositiveCoefMutants,
  highNegativeCoefMutants,
  highPositiveCoefMutants
) %>%
  select(locusId, kofamFunction) %>%
  distinct() %>%
  mutate(kofamFunction = replace_na(kofamFunction, "Unannotated"))

# Initialize label column to NA
resMaaslinNumericDayDj$label <- NA
resMaaslinNumericDayDj <- resMaaslinNumericDayDj %>%
  left_join(label_df, by = c("feature" = "locusId")) %>%
  mutate(label = if_else(!is.na(kofamFunction), kofamFunction, label)) %>%
  select(-kofamFunction)

resMaaslinNumericDayDj %>%
  filter(metadata == 'dayNumeric')%>%
  ggplot(aes(x = coef,
             label = label,
             y= -log10(qval)))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 2.)+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness')

meta= metadata%>%
  filter(tissue != 'T0')%>%
  column_to_rownames('sample')
'This did not really work. '
 # Maaslin2(input = maaslinIn,
 #          input_metadata = meta,
 #          transform = 'none',
 #          normalization = 'none',
 #          min_prevalence = 0,
 #          min_abundance = 0,
 #          min_variance = .5,
 #          fixed_effects = c('tissue','pfClusters', 'millionBases'),
 #          random_effects = ,
 #          output = 'tissueLogRatiosMaaslin2'
 # )

resMaaslinTissue=read_tsv('tissueLogRatiosMaaslin2/all_results.tsv')
highNegativeCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantPositiveCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

label_df <- bind_rows(
  highSignificantNegativeCoefMutants,
  highSignificantPositiveCoefMutants,
  highNegativeCoefMutants,
  highPositiveCoefMutants
) %>%
  select(locusId, kofamFunction) %>%
  distinct() %>%
  mutate(kofamFunction = replace_na(kofamFunction, "Unannotated"))

# Initialize label column to NA
resMaaslinTissue$label <- NA
resMaaslinTissue <- resMaaslinTissue %>%
  left_join(label_df, by = c("feature" = "locusId")) %>%
  mutate(label = if_else(!is.na(kofamFunction), kofamFunction, label)) %>%
  select(-kofamFunction)

resMaaslinTissue %>%
  filter(metadata == 'tissue')%>%
  ggplot(aes(x = coef,
             label = label,
             y= -log10(qval)))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 2.)+
  geom_hline(yintercept = -log10(.1))+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness')

meta= metadata%>%
  filter(tissue != 'T0',
         day == 'day1')%>%
  column_to_rownames('sample')
# Maaslin2(input = maaslinIn,
#          input_metadata = meta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = .5,
#          fixed_effects = c('tissue','pfClusters', 'millionBases'),
#          random_effects = ,
#          output = 'tissueDay1LogRatiosMaaslin2'
# )

resMaaslinTissueDay1=read_tsv('tissueDay1LogRatiosMaaslin2/all_results.tsv')
highNegativeCoefMutants=resMaaslinTissueDay1%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=resMaaslinTissueDay1%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=resMaaslinTissueDay1%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantPositiveCoefMutants=resMaaslinTissueDay1%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

label_df <- bind_rows(
  highSignificantNegativeCoefMutants,
  highSignificantPositiveCoefMutants,
  highNegativeCoefMutants,
  highPositiveCoefMutants
) %>%
  select(locusId, kofamFunction) %>%
  distinct() %>%
  mutate(kofamFunction = replace_na(kofamFunction, "Unannotated"))

# Initialize label column to NA
resMaaslinTissueDay1$label <- NA
resMaaslinTissueDay1 <- resMaaslinTissueDay1 %>%
  left_join(label_df, by = c("feature" = "locusId")) %>%
  mutate(label = if_else(!is.na(kofamFunction), kofamFunction, label)) %>%
  select(-kofamFunction)

resMaaslinTissueDay1 %>%
  filter(metadata == 'tissue')%>%
  ggplot(aes(x = coef,
             label = label,
             y= -log10(qval)))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 2.)+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness')

meta= metadata%>%
  filter(tissue != 'T0',
         day == 'day3')%>%
  column_to_rownames('sample')
 Maaslin2(input = maaslinIn,
          input_metadata = meta,
          transform = 'none',
          normalization = 'none',
          min_prevalence = 0,
          min_abundance = 0,
          min_variance = 0,
          fixed_effects = c('tissue', 'millionBases'),
          random_effects = ,
          output = 'tissueDay3LogRatiosMaaslin2'
 )

meta= metadata%>%
  filter(tissue != 'T0',
         day == 'day7')%>%
  column_to_rownames('sample')
Maaslin2(input = maaslinIn,
         input_metadata = meta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = .0,
         fixed_effects = c('tissue', 'millionBases'),
         random_effects = ,
         output = 'tissueDay7LogRatiosMaaslin2'
)
meta= metadata%>%
  filter(tissue != 'T0',
         day == 'day14')%>%
  column_to_rownames('sample')
Maaslin2(input = maaslinIn,
         input_metadata = meta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = 0,
         fixed_effects = c('tissue','millionBases'),
         random_effects = ,
         output = 'tissueDay14LogRatiosMaaslin2'
)
resMaaslinTissueDay14=read_tsv('tissueDay14LogRatiosMaaslin2/all_results.tsv')

highNegativeCoefMutants=resMaaslinTissueDay14%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=resMaaslinTissueDay14%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=resMaaslinTissueDay14%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantPositiveCoefMutants=resMaaslinTissueDay14%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .05)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

label_df <- bind_rows(
  highSignificantNegativeCoefMutants,
  highSignificantPositiveCoefMutants,
  highNegativeCoefMutants,
  highPositiveCoefMutants
) %>%
  select(locusId, kofamFunction) %>%
  distinct() %>%
  mutate(kofamFunction = replace_na(kofamFunction, "Unannotated"))

# Initialize label column to NA
resMaaslinTissueDay14$label <- NA
resMaaslinTissueDay14 <- resMaaslinTissueDay14 %>%
  left_join(label_df, by = c("feature" = "locusId")) %>%
  mutate(label = if_else(!is.na(kofamFunction), kofamFunction, label)) %>%
  select(-kofamFunction)

resMaaslinTissueDay14 %>%
  filter(metadata == 'tissue')%>%
  ggplot(aes(x = coef,
             label = label,
             y= -log10(qval)))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 2.)+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness')

meta= metadata%>%
  filter(tissue != 'T0',
         phase == 'early')%>%
  column_to_rownames('sample')
Maaslin2(input = maaslinIn,
         input_metadata = meta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = 0,
         fixed_effects = c('tissue','millionBases'),
         random_effects = ,
         output = 'tissueEarlyPhaseLogRatiosMaaslin2'
)

meta= metadata%>%
  filter(tissue != 'T0',
         phase == 'late')%>%
  column_to_rownames('sample')
Maaslin2(input = maaslinIn,
         input_metadata = meta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = 0,
         fixed_effects = c('tissue','millionBases'),
         random_effects = ,
         output = 'tissueLatePhaseLogRatiosMaaslin2'
)