
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

#' Output writing of pKSEA compare() results
#'
#' Output only:
#' uses results from compare(), outputs up to three files labeled full.csv and no_ksea.csv and ksea_only.csv
#' appended to an output name (KSEA-filtered results only if KSEA database was provided to compare()).
#' @importFrom utils write.csv
#' @param full_ksea.results results from compare() including full and optional KSEA excluded and exclusive results
#' @param outputpath parent directory for output
#' @param outputname file name of output
#' @param singlefolder if desired, name of output folder within parent directory. Default is separate folders
#' for each compare() run
#'
#' @keywords internal
#' @export
#'
#'

results_write <- function(full_ksea.results, outputpath, outputname, singlefolder = NULL){
  if(!is.null(singlefolder)){
    outfolder <- file.path(outputpath, singlefolder)
    ifelse(!dir.exists(outfolder), dir.create(outfolder), F)
    write.csv(full_ksea.results$full, file = file.path(outfolder, paste(outputname, "full.csv")),
              quote = F, row.names = T)
    if(any(names(full_ksea.results) == "noksea")){
      write.csv(full_ksea.results$noksea, file = file.path(outfolder, paste(outputname, "no_ksea.csv")),
                quote = F, row.names = T)
      write.csv(full_ksea.results$kseaonly, file = file.path(outfolder, paste(outputname, "ksea_only.csv")),
                quote = F, row.names = T)
    }
  } else {
    runlabel <- mk_runlabel(parentdir = outputpath, customsuffix = outputname)
    write.csv(full_ksea.results$full, file = file.path(outputpath, runlabel, paste(outputname, "full.csv")),
              quote = F, row.names = T)
    if(any(names(full_ksea.results) == "noksea")){
      write.csv(full_ksea.results$noksea, file = file.path(outputpath, runlabel, paste(outputname, "noksea.csv")),
                quote = F, row.names = T)
      write.csv(full_ksea.results$kseaonly, file = file.path(outfolder, paste(outputname, "ksea only.csv")),
                quote = F, row.names = T)
    }
  }
}

#' mk_runlabel()
#'
#' Utility function for generating new identifiers for each run, labeled by time run was initiated and
#' custom suffix
#'
#' @param parentdir parent directory
#' @param customsuffix additional suffix to run identifier
#' @keywords internal
#' @export
#'

mk_runlabel <- function(parentdir= getwd(), customsuffix){
  runlabel <- paste(format(Sys.time(), "%F %H-%M"), customsuffix)
  ifelse(!dir.exists(file.path(parentdir, runlabel)), dir.create(file.path(parentdir, runlabel)), F)
  return(runlabel)
}

#' Running pKSEA::compare() on multiple files
#'
#' For running compare() on multiple CSV data files in the same directory and for writing results to a folder in the
#' designated data directory. Can receive various arguments to be passed on to downstream functions.
#'
#' @importFrom utils read.csv
#' @importFrom utils write.csv
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
#' @param ... parameters to be passed on to downstream functions, including(default):
#' n_permutations (1000), seed (123), kseadb (NULL), kin_ens_table (NULL).
#' See \code{\link{run_on_matched}}, \code{\link{compare}} for details.
#' @inheritParams run_on_matched
#'
#' @export


#batch run on all files in summaryfiledir with a common file string label

batchrun <- function(summaryfiledir, commonfilestring = ".csv", predictionDB, results_folder = NULL, ...){

  sum_manifest <- list.files(summaryfiledir, pattern = commonfilestring)

  for (i in 1:length(sum_manifest)){

    summary_data <- read.csv(file.path(summaryfiledir, sum_manifest[i]))
    matched_data <- get_matched_data(summary_data, predictionDB, ...)


    run_results <- compare(matched_data, predictionDB, ...)

    if(!is.null(results_folder)){
      results_write(run_results, summaryfiledir, gsub(pattern = ".csv", replacement = "", sum_manifest[i]), singlefolder = results_folder)
    } else{
      results_write(run_results, summaryfiledir, gsub(pattern = ".csv", replacement = "", sum_manifest[i]))
    }
  }
}
