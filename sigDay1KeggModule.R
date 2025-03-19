
library(tidyverse)
library(KEGGREST)
library(progressr)

# Set up the progress handler
progressr::handlers("txtprogressbar")

# KEGG module fetcher
get_module <- function(accession) {
  keggout <- keggGet(accession)
  if (!is.null(keggout[[1]]$BRITE)) {
    return(as.character(keggout[[1]]$BRITE[2]))
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

day1SigNegative <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day1SigNegative.tsv')

negativeDay1Ranked <- day1SigNegative %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

negativeDay1Ranked <- annotate_modules(negativeDay1Ranked)

day1SigPositive <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day1SigPositive.tsv')

positiveDay1Ranked <- day1SigPositive %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

positiveDay1Ranked <- annotate_modules(positiveDay1Ranked)


# Repeat for each dataset
day1SigPositiveDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjPositiveDay1.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day1SigPositiveDj <- annotate_modules(day1SigPositiveDj)


day1SigNegativeDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjNegativeDay1.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day1SigNegativeDj <- annotate_modules(day1SigNegativeDj)


day1SigPositiveColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonPositiveDay1.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day1SigPositiveColon <- annotate_modules(day1SigPositiveColon)


day1SigNegativeColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonNegativeDay1.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day1SigNegativeColon <- annotate_modules(day1SigNegativeColon)


# Example summary
day1SigNegativeColon %>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day1ColonNegativeModulesRanked.tsv')

day1SigPositiveColon %>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day1ColonPositiveModulesRanked.tsv')

day1SigNegativeDj%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day1DjNegativeModulesRanked.tsv')

day1SigPositiveDj%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day1DjPositiveModulesRanked.tsv')

negativeDay1Ranked%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('negativeDay1ModulesRanked.tsv')

positiveDay1Ranked%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('positiveDay1ModulesRanked.tsv')