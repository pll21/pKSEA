
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
#' @examples
#' #Read in example summary statistics dataset from csv
#' summarydata_ex <- read.csv(system.file("extdata", "example_data1.csv", package="pKSEA"))
#'
#' #Get matched data using predictions from NetworKIN
#' matched_data_ex <- get_matched_data(summarydata_ex, NetworKINPred_db)

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


#' Runs pKSEA analysis on a dataset result from get_matched_data.
#'
#' Calculates score contributions from summary statistics (tscore) and prediction scores, and sums contribution scores
#' by kinase to calculate raw kinase activity change scores (KAC scores). Performs permutation test on summary statistic
#' data to assess significance of kinase activity change scores, and reports significance as a percentile score
#' (pKSEA significance score).
#'
#'
#' @param matched_data data after filtering against predictions (results from get_matched_data())
#' @param n_permutations number of mutations to perform (default 1000)
#' @param seed seed used for permutation testing
#' @param kin_ens_table optional table for inclusion of matched ensembl ids for kinases, with columns: ens = ensembl id,
#' kinases = kinase_id as otherwise used
#'
#' @export
#' @examples
#' #Read in example summary statistics dataset from csv
#' summarydata_ex <- read.csv(system.file("extdata", "example_data1.csv", package="pKSEA"))
#'
#' #Get matched data using predictions from NetworKIN
#' matched_data_ex <- get_matched_data(summarydata_ex, NetworKINPred_db)
#'
#' #Perform single run of pKSEA analysis
#' single_run_results_ex <- run_on_matched(matched_data_ex, n_permutations = 10)
#'


run_on_matched <- function (matched_data,
                            n_permutations= 1000,
                            seed = 123,
                            kin_ens_table = NULL){

  matched_calcs <- calc_contribution(matched_data)
  scores <- getscores(matched_calcs)
  substrate_list <- getsubs(matched_calcs)
  perms <- permtest(matched_calcs, n_permutations, seed)
  perm_results <- perc.permutation(scores, perms)
  perm_results <- perm_results[order(perm_results$permutationScore),]

  if(!is.null(kin_ens_table)){
  perm_results$ens <- kin_ens_table$ens[match(row.names(perm_results), kin_ens_table$kinases)]
  }

  return(perm_results)
}

#' Running analysis runs on known substrates, predicted substrates, and both.
#'
#' Performs up to three run_on_matched() runs on summary-prediction matcheddata from \code{get_matched_data()},
#' returning permutation significance score results.
#' If a KSEA database is provided for filtering and comparison, one full analysis will be performed on all
#' phosphosites, one on data with all known kinase substrates removed according to the provided KSEA database,
#' and one on known kinase substrates only.
#'
#' @usage compare(matched_data, predictionDB, kseadb, ...)
#' @param matched_data File path to summary statistic phosphoproteomics CSV data file
#' with an entry for each phosphopeptide. Required data file columns:
#' GN = gene name identifier that will be matched with prediction database, Peptide = unique peptide identifier
#' (for example, sequence with modifications), Phosphosites = comma-separated phosphorylation sites (eg. "T102,S105"),
#' pval= pairwise test p-value, fc= mean fold change, t= pairwise test t-statistic. pval and fc are used for results
#' reporting only, all others are important for database searching, calculation, and permutation testing.
#'
#' @param predictionDB Input database whose prediction scores will be used for calculations. Required columns:
#' substrate_name= name of substrate corresponding to GN in summary_data, kinase_id = identifiers for kinase predictors,
#' position= phosphorylated residue number, score = numeric score for strength of prediction.
#' @param kseadb Optional KSEA database for filtering purposes. Containing substrate gene name "SUB_GENE"
#' and phosphorylated residue "SUB_MOD_RSD" in standard form (ie. T302).
#' @param ... optional parameters to be passed on to downstream functions, including (default):
#' n_permutations (1000), seed (123), kin_ens_table (NULL). See \code{\link{run_on_matched}} for details.
#'
#' @export
#'
#' @examples
#' #Read in example summary statistics dataset from csv
#' summarydata_ex <- read.csv(system.file("extdata", "example_data1.csv", package="pKSEA"))
#'
#' #Get matched data using predictions from NetworKIN
#' matched_data_ex <- get_matched_data(summarydata_ex, NetworKINPred_db)
#'
#' #Perform comparative analysis using provided KSEAdb as filter
#' \dontrun{
#' compare_results_ex <- compare(matched_data_ex, kseadb = KSEAdb, n_permutations = 10)
#' }

compare <- function(matched_data,
                    predictionDB,
                    kseadb= NULL,
                    ...){
  results <- list()

  results$full <- run_on_matched(matched_data, ...)

  if(!is.null(kseadb)){
    matched_KSEAfiltered <- KSEAfilter(matched_data = matched_data, kseadb = kseadb)
    results$noksea <- run_on_matched(matched_KSEAfiltered, ...)

    matched_KSEAonly <- KSEAfilter(matched_data = matched_data, reverse = T , kseadb = kseadb)
    results$kseaonly <- run_on_matched(matched_KSEAonly, ...)
  }


  return(results)
}


#' Running pKSEA::compare() on multiple files
#'
#' For running compare() on multiple CSV data files in the same directory and for writing results to a folder in the
#' designated data directory. Can receive various arguments to be passed on to downstream functions. Writes to tempdir()
#' unless \code{outputpath} variable is specified by user (argument passed on to \code{\link{results_write}}).
#'
#' @importFrom utils read.csv
#' @usage batchrun(summaryfiledir, commonfilestring = ".csv",
#' predictionDB, results_folder = NULL, ...)
#' @param summaryfiledir Directory containing summary statistic CSV files. Required data file columns:
#' GN = gene name identifier that will be matched with prediction database, Peptide = unique peptide identifier
#' (for example, sequence with modifications), Phosphosites = comma-separated phosphorylation sites (eg. "T102,S105"),
#' pval= pairwise test p-value, fc= mean fold change, t= pairwise test t-statistic. pval and fc are used for results
#' reporting only, all others are important for database searching, calculation, and permutation testing.
#' @param commonfilestring Common string identifying all files to be included in analysis
#' @param predictionDB Input database whose prediction scores will be used for calculations. Required columns:
#' substrate_name= name of substrate corresponding to GN in summary_data, kinase_id = identifiers for kinase predictors,
#' position= phosphorylated residue number, score = numeric score for strength of prediction.
#' @param results_folder if desired, a single output folder. Else each run performed on each file
#' will have a separate output folder identified by run initiation time.
#' @param ... parameters to be passed on to downstream functions, including(default): outputpath (tempdir())
#' n_permutations (1000), seed (123), kseadb (NULL), kin_ens_table (NULL).
#' See \code{\link{run_on_matched}}, \code{\link{compare}} for details.
#' @inheritParams run_on_matched
#'
#' @export
#' @examples
#' #point to data directory that contains summary .csv files
#' datapath <- system.file("extdata", package = "pKSEA")
#'
#' #run batchrun function to analyze all files in that folder, with options
#' batchrun(datapath, predictionDB=NetworKINPred_db, kseadb = KSEAdb, n_permutations = 5)





#batch run on all files in summaryfiledir with a common file string label

batchrun <- function(summaryfiledir, commonfilestring = ".csv", predictionDB, results_folder = NULL, ...){

  sum_manifest <- list.files(summaryfiledir, pattern = commonfilestring)

  for (i in 1:length(sum_manifest)){

    summary_data <- read.csv(file.path(summaryfiledir, sum_manifest[i]))

    matched_data <- do.call(get_matched_data, c(list(summary_data, predictionDB)))


    run_results <- compare(matched_data, predictionDB, ...)


    results_write(run_results, outputname = gsub(pattern = ".csv", replacement = "", sum_manifest[i]), singlefolder = results_folder)

  }
}
