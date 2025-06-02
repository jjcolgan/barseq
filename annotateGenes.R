library(KEGGREST)
library(tidyverse)

# function to get indent level based on spaces at start of string
get_level <- function(x) {
  n_spaces <- str_length(str_extract(x, "^\\s*"))
  n_spaces / 1  # 1 space = 1 level, or adjust if your indent is 2/3/4 spaces
}

# function to parse KEGG BRITE output into a tidy dataframe
parseBrite <- function(keggOutput) {
  brite <- keggOutput[[1]]$BRITE

  brite_df <- tibble(raw = brite) %>%
    mutate(
      briteLevel = get_level(raw),
      code  = str_extract(raw, "^\\s*(\\d{5}|K\\d{5})"),
      KO    = str_detect(raw, "^\\s*K\\d{5}"),
      briteDescription = str_trim(str_remove(raw, "^\\s*(\\d{5}|K\\d{5})\\s*"))
    )%>%
    .[1:5,]

  max_level <- max(brite_df$briteLevel)
  level_cols <- paste0("brite_", 0:max_level)
  current_levels <- rep(NA_character_, max_level + 1)

  # add empty columns first
  brite_df <- brite_df %>%
    mutate(!!!setNames(rep(list(NA_character_), length(level_cols)), level_cols))

  result <- brite_df %>%
    rowwise() %>%
    mutate(
    {
      current_levels[briteLevel + 1] <<- briteDescription
      if (briteLevel < max_level) {
        current_levels[(briteLevel + 2):(max_level + 1)] <<- NA_character_
      }
      across(all_of(level_cols), ~current_levels[as.integer(str_extract(cur_column(), "\\d+")) + 1])
    }
    ) %>%
    ungroup()

  result %>%
    filter(KO) %>%
    select(code, briteDescription, all_of(level_cols))
}

# example KEGG annotation lookup on a gene
genes <- read_tsv(file = 'genesWithAnvioAnnotations.tsv')

library(KEGGREST)  # For fetching KEGG data
library(tidyverse) # For data manipulation
library(progress)  # For displaying progress bar in loops

# Function: get_level
# Purpose: Calculate indentation level of a BRITE line based on leading spaces
# Input: a string (one line of BRITE output)
# Output: integer indicating level (number of leading spaces)
get_level <- function(x) {
  n_spaces <- str_length(str_extract(x, "^\\s*")) # Count leading spaces
  n_spaces / 1  # Each space counts as one level
}

# Function: parseBrite
# Purpose: Parse BRITE hierarchy output from KEGG and create a tidy dataframe with hierarchy levels as columns
# Input: keggOutput - a list object returned by keggGet() for one KO accession
# Output: tibble with KO codes, description, and hierarchy columns (brite_0, brite_1, ...)
parseBrite <- function(keggOutput) {
  brite <- keggOutput[[1]]$BRITE  # Extract BRITE section

  # If BRITE section is missing, return NULL
  if (is.null(brite)) return(NULL)

  # Create initial tibble with raw BRITE lines and extract useful parts
  brite_df <- tibble(raw = brite) %>%
    mutate(
      briteLevel = get_level(raw), # Calculate indentation level
      code  = str_extract(raw, "^\\s*(\\d{5}|K\\d{5})"), # Extract numeric or KO code
      KO    = str_detect(raw, "^\\s*K\\d{5}"),          # Identify if line is a KO entry
      briteDescription = str_trim(str_remove(raw, "^\\s*(\\d{5}|K\\d{5})\\s*")) # Clean description text
    )

  max_level <- max(brite_df$briteLevel)                 # Maximum indentation level found
  level_cols <- paste0("brite_", 0:max_level)            # Create column names for each level
  current_levels <- rep(NA_character_, max_level + 1)    # Initialize vector to track current path in hierarchy

  # Initialize the hierarchy columns as NA
  brite_df <- brite_df %>%
    mutate(!!!setNames(rep(list(NA_character_), length(level_cols)), level_cols))

  # Fill hierarchy columns by tracking the current category at each level
  result <- brite_df %>%
    rowwise() %>%
    mutate(
    {
      # Update current level category
      current_levels[briteLevel + 1] <<- briteDescription

      # Reset deeper levels when going up the hierarchy
      if (briteLevel < max_level) {
        current_levels[(briteLevel + 2):(max_level + 1)] <<- NA_character_
      }

      # Assign current hierarchy path to respective columns
      across(all_of(level_cols), ~current_levels[as.integer(str_extract(cur_column(), "\\d+")) + 1])
    }
    ) %>%
    ungroup()

  # Return only rows corresponding to KO entries, selecting relevant columns
  result %>%
    filter(KO) %>%
    select(code, briteDescription, all_of(level_cols))
}

# Load your genes data with KO accessions and locus IDs
genes <- read_tsv(file = 'genesWithAnvioAnnotations.tsv')

# Initialize list to store parsed BRITE data for each gene
briteResults <- list()

# Initialize a progress bar to track processing progress
pb <- progress_bar$new(
  total = nrow(genes),                # Total steps = number of genes
  format = "[:bar] :current/:total (:percent) | eta: :eta", # Progress bar format
  clear = FALSE,
  width = 60
)

# Loop over each gene to fetch KEGG info and parse BRITE hierarchies
for (i in seq_len(nrow(genes))) {
  pb$tick()   # Update progress bar on each iteration

  temp <- genes[i,]  # Current gene row

  # Skip if KO accession is missing or empty
  if (is.na(temp$kofamAccession) | temp$kofamAccession == "") next

  # Safely attempt to fetch KEGG data; skip on error
  testOutput <- tryCatch(keggGet(temp$kofamAccession), error = function(e) NULL)

  # Skip if no valid KEGG output
  if (is.null(testOutput[[1]])) next

  # Extract primary and secondary pathways if available
  pathways <- testOutput[[1]]$PATHWAY
  primaryPathway <- pathways[1]
  secondaryPathway <- pathways[2]

  # Parse BRITE hierarchy for the current KO
  b <- parseBrite(testOutput)

  # If parsing was successful, add pathway and locus ID info and store
  if (!is.null(b)) {
    b <- b %>%
      mutate(
        primaryPathway = primaryPathway,
        secondaryPathway = secondaryPathway,
        locusId = temp$locusId
      )
    briteResults[[length(briteResults) + 1]] <- b
  }
}

# Combine all parsed BRITE results into one dataframe
briteOutDF <- bind_rows(briteResults)

# Join the BRITE data back to the original genes dataframe by locusId
final <- genes %>%
  left_join(briteOutDF, by = 'locusId')

# View the final merged dataframe with BRITE hierarchy and pathway info
view(final)

