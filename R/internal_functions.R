
#' Calculate score contributions by phosphorylation site
#'
#'
#' @param matched_data Input
#' @export
#' @keywords internal
#'

calc_contribution <- function(matched_data){
  # matched_data$contribution <- matched_data$zscore * log(matched_data$score)
  # matched_data$contribution <- matched_data$t *  log(matched_data$score)
  # matched_data$contribution <- -log(matched_data$pval) * log(matched_data$score, base = loglevel)
  matched_data$contribution <- matched_data$t * matched_data$score
  return(matched_data)
}

#' Sum score contributions for each kinase across all phosphopeptides
#'
#' @param matched_data Input with calculated contributions
#' @export
#' @keywords internal
#'
getscores <- function(matched_data){
  #calculate contributions
  matched_data <- calc_contribution(matched_data)

  #initialize data table from kinases represented in prediction dataframe
  kin_act_ch <- data.frame(row.names = levels(factor(matched_data$kinase_id)))

  #calculate totals
  for(i in 1:length(levels(factor(matched_data$kinase_id)))){
    kin_act_ch[i,1] <- sum(matched_data$contribution[matched_data$kinase_id == rownames(kin_act_ch)[i]])
  }

  colnames(kin_act_ch)<- "rawKAC"
  return(kin_act_ch)
}

#' Extract summary table with pertinent columns related to included substrates
#'
#' @param matched_data Input with calculated contributions
#' @export
#' @keywords internal
#'
#'
getsubs <- function(matched_data){
  #contributions
  matched_data <- calc_contribution(matched_data)

  substrates_out <- matched_data[, c("GN",
                                     "Phosphosites",
                                     "kinase_id",
                                     "score",
                                      "t",
                                      "pval", "fc", "contribution")]
  substrates_out <- substrates_out[order(substrates_out$kinase_id, substrates_out$contribution),]
  rownames(substrates_out) <- NULL
  return(substrates_out)
}


#' Perform permutation test
#'
#' Returns a table that has permuted the relationship between phosphopeptides and summary statistics
#' (ie. fold change, t-score)
#'
#' @param matched_data Input with calculated contributions
#' @param perms Number of permutations to run, default = 1000
#' @export
#' @keywords internal
#'
#'

permtest <- function(matched_data, perms = 1000, seed = 1){
  permtable <- data.frame(row.names = levels(factor(matched_data$kinase_id)))

  peptidelist <- matched_data[!duplicated(matched_data$Peptide), c("Peptide", "GN", "kinase_id", "siteid", "fc",
                                                                  # "zscore",
                                                                    "t",
                                                                    "pval")]

  set.seed(seed)
  for (i in 1:perms){
    #randomly reassign siteid calls
    randompep <- sample(length(levels(factor(matched_data$Peptide))), replace = F)
    randomframe <- matched_data
    randomframe$Peptide <- factor(randomframe$Peptide)
    randomframe[,c("t","pval")] <- randomframe[match(as.numeric(randomframe$Peptide), randompep),c("t", "pval")]

    currentperm <- getscores(randomframe)

    permtable <- cbind(permtable, currentperm)
  }

  return(permtable)
}

#' Obtain percentile rank comparing a single value to set
#'
#'
#' @param set Set of values to which given value will be compared
#' @param value Value for which percentile score will be calculated
#' @export
#' @keywords internal
#'
perc.rank <- function (set, value){
  length(set[set <=value])/length(set)*100
}

#' Get percentile ranks across permutations
#'
#' @param results Results of kinase scoring
#' @param permutations Results of permutations
#' @export
#' @keywords internal

perc.permutation <- function (results, permutations){
  scored <- as.data.frame(results)
  for (i in 1:nrow(results)){
    scored$permutationScore[i] <- perc.rank(permutations[i,], results[i,])
  }
  colnames(scored)[1] <- "rawKAC"
  return(scored)
}

#' Filter matched data to remove positive IDs from KSEA
#'
#' @param matched_data Results of get_matched_data function
#' @param kseadb KSEA database containing substrate gene name "SUB_GENE" and phosphorylated residue
#' "SUB_MOD_RSD" in standard form (ie. T302).
#' @export
#' @keywords internal


KSEAfilter <- function(matched_data, kseadb, reverse= F){

  #Generate siteid
  kseadb$siteid <- paste(kseadb$SUB_GENE, kseadb$SUB_MOD_RSD, sep = "_")


  #extract elements of netinput_filt that do not have corresponding KSEA info
  if(reverse== F) {
    matched_filtered <- subset(matched_data, !(matched_data$siteid %in% kseadb$siteid))
  } else {
    matched_filtered <- subset(matched_data, (matched_data$siteid %in% kseadb$siteid))
  }
  return(matched_filtered)
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
