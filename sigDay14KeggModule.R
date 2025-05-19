library(tidyverse)
library(KEGGREST)
library(progressr)

# Set up the progress handler
progressr::handlers("txtprogressbar")

# KEGG module fetcher
# KEGG module fetcher
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

# Load annotations and datasets
keggAnno <- read_tsv('genesWithAnvioAnnotations.tsv') %>%
  select(locusId, kofamAccession) %>%
  distinct()

day14SigNegative <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day14SigNegative.tsv')

negativeDay14Ranked <- day14SigNegative %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

negativeDay14Ranked <- annotate_modules(negativeDay14Ranked)

day14SigPositive <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day14SigPositive.tsv')

positiveDay14Ranked <- day14SigPositive %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

positiveDay14Ranked <- annotate_modules(positiveDay14Ranked)


# Repeat for each dataset
day14SigPositiveDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjPositiveDay14.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day14SigPositiveDj <- annotate_modules(day14SigPositiveDj)


day14SigNegativeDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjNegativeDay14.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day14SigNegativeDj <- annotate_modules(day14SigNegativeDj)


day14SigPositiveColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonPositiveDay14.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day14SigPositiveColon <- annotate_modules(day14SigPositiveColon)


day14SigNegativeColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonNegativeDay14.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day14SigNegativeColon <- annotate_modules(day14SigNegativeColon)


# Example summary
day14SigNegativeColon %>%
  write_tsv('day14ColonNegativeModulesRanked.tsv')

day14SigPositiveColon %>%
  write_tsv('day14ColonPositiveModulesRanked.tsv')

day14SigNegativeDj%>%
  write_tsv('day14DjNegativeModulesRanked.tsv')

day14SigPositiveDj%>%
  write_tsv('day14DjPositiveModulesRanked.tsv')

negativeDay14Ranked%>%
  write_tsv('negativeDay14ModulesRanked.tsv')

positiveDay14Ranked%>%
  write_tsv('positiveDay14ModulesRanked.tsv')