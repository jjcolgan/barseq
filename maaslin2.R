library(tidyverse)
library(Maaslin2)
library(ggrepel)
library(pwr)
library(KEGGREST)
library(progressr)
library(stringr)

get_module <- function(accession) {
  keggout <- keggGet(accession)
  if (!is.null(keggout[[1]]$BRITE)) {
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
Maaslin2(input = maaslinIn,
         input_metadata = colonMeta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = .0,
         fixed_effects = c('day','percentPerfectBarcode', 'millionBases'),
         random_effects =c('cage') ,
         reference = 'day,day1;day3;day7;day14',
         output = 'colonLogRatiosMaaslin2'
)
maaslinRes = read_tsv('colonLogRatiosMaaslin2/all_results.tsv')

maaslinRes %>%
  ggplot(aes(x = pval))+
  geom_histogram(bins = 100)

sigResMaaslin=read_tsv('colonLogRatiosMaaslin2/significant_results.tsv')

sigResMaaslin %>%
  filter(metadata == 'day')%>%
  filter(qval < .1)%>%
  select(feature)%>%
  distinct()%>%
  nrow()

  Maaslin2(input = maaslinIn,
          input_metadata = colonMeta,
           transform = 'none',
           normalization = 'none',
           min_prevalence = 0,
           min_abundance = 0,
           min_variance = .0,
           fixed_effects = c('dayNumeric','percentPerfectBarcode', 'millionBases'),
           random_effects =c('cage') ,
           reference = '',
           output = 'colonLogRatiosMaaslin2NumericDay'
 )

sigResMaaslinNumericDayColon=read_tsv('colonLogRatiosMaaslin2NumericDay/significant_results.tsv')
resMaaslinNumericDayColon=read_tsv('colonLogRatiosMaaslin2NumericDay/all_results.tsv')

resMaaslinNumericDayColon%>%
  ggplot(aes(x = pval))+
  geom_histogram(bins = 100)

sigResMaaslinNumericDayColon %>%
  filter(metadata == 'dayNumeric')%>%
  filter(qval < .1)%>%
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
  filter(qval < .1)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantPositiveCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
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
             col = qval < .1,
             y= -log10(qval)))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 1.5)+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness Co',
       col = 'P-adjust < .1')+
  theme_bw()


briteAnnotatedTemporalSignificantColon=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  annotate_modules()

briteAnnotatedTemporalSignificantColon$selection = 'Negative'
briteAnnotatedTemporalSignificantColon$selection[briteAnnotatedTemporalSignificantColon$coef < 0] = 'Positive'

briteAnnotatedTemporalSignificantColon%>%
  mutate(moduleLevel1 = str_remove(moduleLevel1, "^\\s*[[:digit:]]+\\s+"), # Remove leading spaces & number
         moduleLevel2 = str_remove(moduleLevel2, "^\\s*[[:digit:]]+\\s+"),
         moduleLevel1 = wrap_format(20)(moduleLevel1),  # Wrap text for readability
         moduleLevel1 = fct_reorder(moduleLevel1, .x = moduleLevel1, .fun = function(x) length(x), .desc = TRUE)) %>%  # Order by count
  ggplot(aes(x = moduleLevel1, fill = moduleLevel2)) +
  geom_bar() +
  theme_minimal()+
  coord_flip()+
  facet_wrap(~selection)+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness Co',
       x = 'Count',
       y = 'Brite level 1',
       fill = 'Brite level 2')

 # Maaslin2(input = maaslinIn,
 #         input_metadata = djMeta,
 #         transform = 'none',
 #          normalization = 'none',
 #          min_prevalence = 0,
 #         min_abundance = 0,
 #          min_variance = .5,
 #          fixed_effects = c('day','percentPerfectBarcode', 'millionBases'),
 #          random_effects = ,
 #          output = 'djLogRatiosMaaslin2'
 # )
Maaslin2(input = maaslinIn,
         input_metadata = djMeta,
         transform = 'none',
         normalization = 'none',
         min_prevalence = 0,
         min_abundance = 0,
         min_variance = .0,
         fixed_effects = c('dayNumeric','percentPerfectBarcode', 'millionBases'),
         random_effects = ,
         reference = '',
         output = 'djLogRatiosMaaslin2NumericDay'
)

sigResMaaslinNumericDayDj=read_tsv('djLogRatiosMaaslin2NumericDay/significant_results.tsv')

sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  ggplot(aes(x=pval))+
  geom_histogram(bins = 50)

resMaaslinNumericDayDj=read_tsv('djLogRatiosMaaslin2NumericDay/all_results.tsv')
highNegativeCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantPositiveCoefMutants=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
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

briteAnnotatedTemporalSignificantDj=sigResMaaslinNumericDayDj%>%
  filter(metadata == 'dayNumeric')%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  annotate_modules()

briteAnnotatedTemporalSignificantDj$selection = 'Positive'
briteAnnotatedTemporalSignificantDj$selection[briteAnnotatedTemporalSignificantDj$coef < 0] = 'Negative'

briteAnnotatedTemporalSignificantDj%>%
  mutate(moduleLevel1 = str_remove(moduleLevel1, "^\\s*[[:digit:]]+\\s+"), # Remove leading spaces & number
         moduleLevel2 = str_remove(moduleLevel2, "^\\s*[[:digit:]]+\\s+"),
         moduleLevel1 = wrap_format(20)(moduleLevel1),  # Wrap text for readability
         moduleLevel1 = fct_reorder(moduleLevel1, .x = moduleLevel1, .fun = function(x) length(x), .desc = TRUE)) %>%  # Order by count
  ggplot(aes(x = moduleLevel1, fill = moduleLevel2)) +
  geom_bar() +
  coord_flip() +  # Flip bars for better readability
  theme_minimal()+
  facet_wrap(~selection)+
  labs(x = 'Count',
       y = 'Brite level 1',
       fill = 'Brite level 2')+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness DJ')




resMaaslinNumericDayDj %>%
  filter(metadata == 'dayNumeric')%>%
  ggplot(aes(x = coef,
             label = label,
             y= -log10(qval),
             col = qval < .1))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 1.5)+
  labs(title = 'Genes assoicated with temporal increases or decreases in fitness DJ')+
  theme_bw()+
  labs(col = 'P-adjust < .1')

meta= metadata%>%
  filter(tissue != 'T0')%>%
  column_to_rownames('sample')
'adding day makes this look better'
 # Maaslin2(input = maaslinIn,
 #           input_metadata = meta,
 #           transform = 'none',
 #           normalization = 'none',
 #           min_prevalence = 0,
 #           min_abundance = 0,
 #           min_variance =.5,
 #           fixed_effects = c('tissue','numericDay','pfClusters', 'millionBases', 'cage'),
 #           reference = '',
 #           random_effects = ,
 #           output = 'tissueLogRatiosMaaslin2'
 #  )

resMaaslinTissue=read_tsv('tissueLogRatiosMaaslin2/all_results.tsv')
highNegativeCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highPositiveCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)



highSignificantPositiveCoefMutants=resMaaslinTissue%>%
  filter(metadata == 'tissue')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
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
             col = qval < .1,
             y= -log10(qval)))+
  geom_point(alpha = .25)+
  geom_text_repel(size = 2.)+
  geom_hline(yintercept = -log10(.1))+
  labs(title = 'Genes assoicated with tissue')+
  labs(col = 'P-adjust < .1')

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
 #          min_variance = 0.,
 #          fixed_effects = c('tissue','numericDay', 'millionBases', 'cage'),
 #          random_effects = ,
 #          output = 'tissueDay1LogRatiosMaaslin2'
 # )

resMaaslinTissueDay1=read_tsv('tissueDay1LogRatiosMaaslin2/all_results.tsv')
highNegativeCoefMutants=resMaaslinTissueDay1%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=resMaaslinTissueDay1%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
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
 # Maaslin2(input = maaslinIn,
 #          input_metadata = meta,
 #          transform = 'none',
 #          normalization = 'none',
 #          min_prevalence = 0,
 #          min_abundance = 0,
 #          min_variance = ,
 #          fixed_effects = c('tissue','numericDay', 'millionBases', 'cage'),
 #          random_effects = ,
 #          output = 'tissueDay3LogRatiosMaaslin2'
 # )

meta= metadata%>%
  filter(tissue != 'T0',
         day == 'day7')%>%
  column_to_rownames('sample')
# Maaslin2(input = maaslinIn,
#          input_metadata = meta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = 0,
#          fixed_effects = c('tissue','numericDay','millionBases', 'cage'),
#          random_effects = '',
#          output = 'tissueDay7LogRatiosMaaslin2'
# )
meta= metadata%>%
  filter(tissue != 'T0',
         day == 'day14')%>%
  column_to_rownames('sample')

# Maaslin2(input = maaslinIn,
#          input_metadata = meta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = 1,
#          fixed_effects = c('tissue','numericDay', 'millionBases', 'cage'),
#          random_effects = ,
#          output = 'tissueDay14LogRatiosMaaslin2',
#          max_significance = .1
# )
resMaaslinTissueDay14=read_tsv('tissueDay14LogRatiosMaaslin2/all_results.tsv')

highNegativeCoefMutants=resMaaslinTissueDay14%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(coef)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')%>%
  head(5)

highSignificantNegativeCoefMutants=resMaaslinTissueDay14%>%
  filter(metadata == 'tissue')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
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
  filter(qval < .1)%>%
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
  #geom_text_repel(size = 2.)+
  geom_hline(yintercept = -log10(.1))+
  labs(title = 'Genes assoicated SI or LI day 14')


'power calculation for f31'

powah=maaslinIn %>%
  rownames_to_column('mutant')%>%
  filter(mutant == 'BBR_RS14230' )%>%
  pivot_longer(cols = (2:ncol(.)), names_to = 'sample', values_to = 'fitness')%>%
  merge(rownames_to_column(meta, 'sample'), by = 'sample')

meanDj = powah %>%
  filter(tissue == 'dj')%>%
  .$fitness%>%
  mean()

meanColon = powah %>%
  filter(tissue == 'colon')%>%
  .$fitness%>%
  mean()

meanDif = meanDj - meanColon

sdDj = powah %>%
  filter(tissue == 'dj')%>%
  .$fitness%>%
  sd()

sdColon = powah %>%
  filter(tissue == 'colon')%>%
  .$fitness%>%
  sd()

nDj = powah %>%
  filter(tissue == 'dj')%>%
  nrow()

nColon = powah %>%
  filter(tissue == 'colon')%>%
  nrow()

pooled_sd <- sqrt(((nColon - 1)*sdColon^2 + (nDj - 1)*sdDj^2) / (nColon + nDj - 2))
cohens_d <- meanDif / pooled_sd
pwr.t.test(
  d = cohens_d,
  power = 0.8,         # desired power
  sig.level = 0.05,    # significance threshold
  type = "two.sample"
)

meta= metadata%>%
  filter(tissue != 'T0',
         phase == 'early')%>%
  column_to_rownames('sample')
# Maaslin2(input = maaslinIn,
#          input_metadata = meta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = 0,
#          fixed_effects = c('tissue','millionBases'),
#          random_effects = ,
#          output = 'tissueEarlyPhaseLogRatiosMaaslin2'
# )

meta= metadata%>%
  filter(tissue != 'T0',
         phase == 'late')%>%
  column_to_rownames('sample')
# Maaslin2(input = maaslinIn,
#          input_metadata = meta,
#          transform = 'none',
#          normalization = 'none',
#          min_prevalence = 0,
#          min_abundance = 0,
#          min_variance = 1,
#          fixed_effects = c('tissue','numericDay','millionBases','pfClusters','cage'),
#          random_effects = ,
#          output = 'tissueLatePhaseLogRatiosMaaslin2'
# )