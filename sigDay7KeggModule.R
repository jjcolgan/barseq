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

day7SigNegative <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day7SigNegative.tsv')

negativeDay7Ranked <- day7SigNegative %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

negativeDay7Ranked <- annotate_modules(negativeDay7Ranked)

day7SigPositive <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day7SigPositive.tsv')

positiveDay7Ranked <- day7SigPositive %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

positiveDay7Ranked <- annotate_modules(positiveDay7Ranked)


# Repeat for each dataset
day7SigPositiveDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjPositiveDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigPositiveDj <- annotate_modules(day7SigPositiveDj)


day7SigNegativeDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjNegativeDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigNegativeDj <- annotate_modules(day7SigNegativeDj)


day7SigPositiveColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonPositiveDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigPositiveColon <- annotate_modules(day7SigPositiveColon)


day7SigNegativeColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonNegativeDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigNegativeColon <- annotate_modules(day7SigNegativeColon)


# Example summary
day7SigNegativeColon %>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day7ColonNegativeModulesRanked.tsv')

day7SigPositiveColon %>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day7ColonPositiveModulesRanked.tsv')

day7SigNegativeDj%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day7DjNegativeModulesRanked.tsv')

day7SigPositiveDj%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('day7DjPositiveModulesRanked.tsv')

negativeDay7Ranked%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('negativeDay7ModulesRanked.tsv')

positiveDay7Ranked%>%
  group_by(module) %>%
  summarise(n = n()) %>%
  arrange(desc(n))%>%
  write_tsv('positiveDay7ModulesRanked.tsv')