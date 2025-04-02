library(tidyverse)
library(KEGGREST)
library(progressr)

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
  filter(tissue == 'colon')

sigResMaaslinNumericDayColon=read_tsv('colonLogRatiosMaaslin2NumericDay/significant_results.tsv')
resMaaslinNumericDayColon=read_tsv('colonLogRatiosMaaslin2NumericDay/all_results.tsv')

sigResMaaslinNumericDayColon%>%
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
  left_join(keggs, by = 'locusId')

highSignificantNegativeCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef < 0)%>%
  filter(qval < .1)%>%
  arrange(qval)%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')
highPositiveCoefMutants=sigResMaaslinNumericDayColon%>%
  filter(metadata == 'dayNumeric')%>%
  filter(coef > 0)%>%
  filter(qval < .1)%>%
  arrange(desc(coef))%>%
  rename(locusId = feature)%>%
  left_join(keggs, by = 'locusId')

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


briteAnnotatedTemporalSignificantColon %>%
  filter(locusId %in% highSignificantNegativeCoefMutants$locusId)%>%
  arrange(desc(-log(qval)))%>%
  view()


plotIn=fitnessScores%>%
  pivot_longer(cols = 4:ncol(fitnessScores), names_to = 'sample', values_to = 'fitness')%>%
  merge(colonMeta, by = 'sample')

plotIn %>%
  filter(locusId == 'BBR_RS15210')%>%
  ggplot(aes(x = dayNumeric,
             y = fitness))+
  geom_point()+
  geom_smooth(method = lm)

title=briteAnnotatedTemporalSignificantColon%>%
  filter(locusId == 'BBR_RS15210')%>%
  .$kofamFunction

plotIn %>%
  filter(locusId == 'BBR_RS15210')%>%
  left_join(briteAnnotatedTemporalSignificantColon, by = "locusId")%>%
  ggplot(aes(x = dayNumeric,
             col = day,
             y = fitness))+
  geom_jitter()+
  labs(x = "day",
       title = title,
       col = 'Day',
       x ='Day',
       y = 'Fitness Score')+
  theme_bw()

title=briteAnnotatedTemporalSignificantColon%>%
  filter(locusId == 'BBR_RS10300')%>%
  .$kofamFunction
plotIn %>%
  filter(locusId == 'BBR_RS10300')%>%
  left_join(briteAnnotatedTemporalSignificantColon, by = "locusId")%>%
  ggplot(aes(x = dayNumeric,
             col = day,
             y = fitness))+
  geom_jitter()+
  labs(x = "day",
       title = title,
       col = 'Day',
       x ='Day',
       y = 'Fitness Score')+
  theme_bw()



title=briteAnnotatedTemporalSignificantColon%>%
  filter(locusId == 'BBR_RS16470')%>%
  .$kofamFunction
plotIn %>%
  filter(locusId == 'BBR_RS16470')%>%
  left_join(briteAnnotatedTemporalSignificantColon, by = "locusId")%>%
  ggplot(aes(x = dayNumeric,
             col = day,
             y = fitness))+
  geom_jitter()+
  labs(x = "day",
       title = title,
       col = 'Day',
       x ='Day',
       y = 'Fitness Score')+
  theme_bw()
