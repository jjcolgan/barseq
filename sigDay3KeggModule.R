library(tidyverse)
library(KEGGREST)
library(progressr)

# Set up the progress handler
progressr::handlers("txtprogressbar")

# KEGG module fetcher
get_module <- function(accession) {
  keggout <- keggGet(accession)
  if (!is.null(keggout[[1]]$BRITE)) {
    return(as.character(keggout[[1]]$BRITE[1]))
  } else {
    return(NA_character_)
  }
}

# Helper function to annotate one table
annotate_modules <- function(df) {
  df$module <- NA
  df_no_na <- df %>% filter(!is.na(kofamAccession))

  with_progress({
    p <- progressor(steps = nrow(df_no_na))
    for (i in seq_len(nrow(df_no_na))) {
      p()
      df_no_na$module[i] <- get_module(df_no_na$kofamAccession[i])
    }
  })

  df[df$kofamAccession %in% df_no_na$kofamAccession, "module"] <- df_no_na$module
  return(df)
}

# Load annotations and datasets
keggAnno <- read_tsv('genesWithAnvioAnnotations.tsv') %>%
  select(locusId, kofamAccession) %>%
  distinct()

day3SigNegative <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day3SigNegative.tsv')

negativeDay3Ranked <- day3SigNegative %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

negativeDay3Ranked <- annotate_modules(negativeDay3Ranked)

day3SigPositive <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day3SigPositive.tsv')

positiveDay3Ranked <- day3SigPositive %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

positiveDay3Ranked <- annotate_modules(positiveDay3Ranked)


# Repeat for each dataset
day3SigPositiveDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjPositiveDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigPositiveDj <- annotate_modules(day3SigPositiveDj)


day3SigNegativeDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjNegativeDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigNegativeDj <- annotate_modules(day3SigNegativeDj)


day3SigPositiveColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonPositiveDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigPositiveColon <- annotate_modules(day3SigPositiveColon)


day3SigNegativeColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonNegativeDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigNegativeColon <- annotate_modules(day3SigNegativeColon)


# Example summary
day3SigNegativeColon %>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day3ColonNegativeModulesRanked.tsv')

day3SigPositiveColon %>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day3ColonPositiveModulesRanked.tsv')

day3SigNegativeDj%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day3DjNegativeModulesRanked.tsv')

day3SigPositiveDj%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day3DjPositiveModulesRanked.tsv')

negativeDay3Ranked%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('negativeDay3ModulesRanked.tsv')

positiveDay3Ranked%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('positiveDay3ModulesRanked.tsv')