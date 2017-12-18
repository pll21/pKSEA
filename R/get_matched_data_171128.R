
#' Filtering data to matched predictions
#'
#'
#' This function reformats summary statistic phosphoproteomicdata to single observations for each phosphorylation site,
#' duplicating other fields for multiple sites on the same peptide. Next, it attempts
#' to find predictions for each phosphorylation site in the provided database. It returns
#' observations (phosphorylation sites) for which a prediction is detected in the database,
#' matching based on HUGO gene name and phosphorylated residue.
#'
#' @param datafull Statistical summary data with an entry for each phosphopeptide. Required columns:
#' GN = gene name identifier that will be matched with prediction database, Peptide = unique peptide identifier
#' (for example, sequence with modifications), Phosphosites = comma-separated phosphorylation sites (eg. "T102,S105"),
#' pval= pairwise test p-value, fc= mean fold change, t= pairwise test t-statistic. pval and fc are used for results
#' reporting only, all others are important for database searching, calculation, and permutation testing.
#'
#' @param predictionDB Input database whose prediction scores will be used for calculations. Required columns:
#' substrate_name= name of substrate corresponding to GN in datafull, kinase_id = identifiers for kinase predictors,
#' position= phosphorylated residue number, score = numeric score for strength of prediction.
#'
#' @export
#'

##Input file annotated with Phosphosites, pval, fc, t (t-score),
## GN (HUGO Gene Name), and Peptide (sequence)

get_matched_data <- function(datafull, predictionDB){

  #Create separate entry for each phosphosite
  multiplePhos <- datafull[grep(",",datafull$Phosphosites),]
  residues <- as.character(multiplePhos$Phosphosites)
  residues <- as.matrix(residues, ncol = 1)
  split <- strsplit(residues, split = ",")
  numsites <- sapply(split, length)

  singlePhos = data.frame(GN = rep(multiplePhos$GN, numsites))
  singlePhos$GN <- rep(multiplePhos$GN, numsites)
  singlePhos$Peptide <- rep(multiplePhos$Peptide, numsites)
  singlePhos$Phosphosites <- unlist(split)
  singlePhos$pval <- rep(multiplePhos$pval, numsites)
  singlePhos$fc <- rep(multiplePhos$fc, numsites)
  # singlePhos$zscore <- rep(multiplePhos$zscore, numsites) #if applicable
  singlePhos$t <- rep(multiplePhos$t, numsites) #if applicable

  # New dataframe with an observation for each phosphosite on each peptide
  # Remove all entries with multiple phosphosites
  newfull <- datafull[!grepl(",", datafull$Phosphosites),]
  newfull <- newfull[,c("Peptide", "pval", "GN", "fc",
                        "t",
                        # "zscore",
                        "Phosphosites")]

  # Rebind all native single phosphosite entries with converted multiple phosphosite entries
  newfull <- rbind(newfull, singlePhos)

  #split number, residue, to facilitate matching
  newfull$pos <- as.numeric(gsub("[A-Z]", "", newfull$Phosphosites))
  newfull$res <- gsub("[0-9]", "", newfull$Phosphosites)

  #merge get pertinent rows
  newfull_filt <- merge(newfull, predictionDB, by.x = c("GN", "pos"), by.y= c("substrate_name", "position"))

  newfull_filt$siteid <- paste(newfull_filt$GN, newfull_filt$Phosphosites, sep = "_")
  # newfull_filt$siteid <- factor(newfull_filt$siteid)
  # newfull_filt$Phosphosites <- NULL
  return(newfull_filt)
}
