# installation scripts:
#
# install splatter from source (e.g.):
# splatter_path <- "Your_path/splatter"
# install.packages(splatter_path, repos=NULL, type="source")
#
# install splatter from github repo (e.g.):
# require(devtools)
# devtools::install_github("YinuoJin/splatter")

suppressMessages(library(splatter))
suppressMessages(library(scater))
suppressMessages(library(proxyC))
suppressMessages(library(reshape2))
suppressMessages(library(msgr))


###############################################
#     Helper functions
###############################################

#' read real data with single cell-type (read dataframe without row & col names)
#' 
#' @param path input data path
#' @param fname input file name
#' 
#' @return matrix of real data counts
readSingleRealExpr <- function(path, fname) {
  df <- t(read.csv(paste0(path, fname), header=FALSE)) 
  counts <- as.matrix(df)
  
  return(counts)
} 

#' read real data with multiple cell-type (read dataframe with row names indicating cell type)
#' 
#' @param path input data path
#' @param fname input file name
#' 
#' @return matrix of real data counts
readGroupRealExpr <- function(path, fname) {
  df <- t(read.csv(paste0(path, fname), row.names=NULL, header=TRUE))
  counts <- df[2:nrow(df),]
  colnames(counts) <- df[1,]
  
  # convert real data observation matrix to numeric type
  mode(counts) = "numeric"
  
  return(counts)
}


#' print statistics of simulated observations
#' 
#' @param counts observation matrix 
#' @param simType type of simulation ("diploid", "CNV")
printSimStats <- function(counts, simType) {
  print(paste0("Stats of simulations type: ", simType))
  
  # sparsity level
  print(paste0("Sparsity level:", sum(counts == 0) / length(counts)))
  
  # mean & variance
  print(paste0("Expression mean:", mean(counts)))
  print(paste0("Expression variance:", var(counts[counts >= 0])))
}


#' simulate diploid observatoins for single cell type
#' 
#' @param counts matrix of real data counts
#' @param n_genes number of genes to simulate
#' @param n_cells number of cells to simulate
#' 
#' @return list of simulated expressions: $sims: simulated SingleCellExperiment object
#'                                        $observations: simulated diploid observations
simSingleDiploidExpr <- function(counts, n_genes, n_cells) {
  
  params <- splatEstimate(SingleCellExperiment(counts))  # initialize SplatParams object
  lib_loc <- getParams(params, "lib.loc")$lib.loc
  
  # update params for diploid simulations
  params <- setParams(params, update = list(
    nGenes = n_genes,
    batchCells = n_cells,
    out.prob = 0, # tmp, force outlier probability to be zero
    seed = sample(1:1000000, 1),
    lib.loc = lib_lic + log(n_genes / nrow(counts))
  ))
  
  # simulation
  sim <- splatSimultate(params = params)
  observations <- counts(diploid_sim)
  
  # return simulation object (SplatSimulate) and observations
  expr_list <- list("sims" = sim, "observations" = observations)
  
  return(expr_list)
}


#' simulate diploid observations for multiple cell types
#' 
#' @param counts matrix of real data counts
#' @param n_genes number of genes for simulation
#' @param n_cells_per_group number of cells in each group for simulation
#' @param n_clusters number of groups (cell types) for simulation
#' 
#' @return list of simulated expressions: $sims: list of simulated SingleCellExperiment object, one for each group 
#'                                        $observations: matrix of simulated diploid observations
simGroupDiploidExpr <- function(counts, n_genes, n_cells_per_group, n_clusters) {
  
  group_names <- unique(colnames(counts))[1:n_clusters] # name of each group indicating each specific cell type
  counts_per_group <- lapply(group_names, function(x) counts[, colnames(counts) == x])
  params_list <- lapply(counts_per_group, function(x) splatEstimate(SingleCellExperiment(x)))  # list of SplatParams object, one for each group
  
  sims <- lapply(params_list, function(params) {
    lib_loc <- getParams(params, "lib.loc")$lib.loc
    
    params <- setParams(params, update = list(
      nGenes = n_genes,
      batchCells = n_cells_per_group, 
      out.prob = 0,  # tmp: force outlier prob = 0
      seed = sample(1:1000000, 1),
      lib.loc = lib_loc + log(n_genes / params@nGenes)
    ))
    
    sim <- splatSimulate(params = params, verbose = FALSE)
    
    return(sim)
  })
  
  # vertically stack diploid observations for each group
  observations <- do.call(cbind, lapply(sims, function(sim) counts(sim)))

  # return list of simulation objects (SplatSimulate) and observations
  expr_list <- list("sims" = sims, "observations" = observations)
  
  return(expr_list)
} 

#' simulate CNV observations for single cell type
#' 
#' @param diploid_sim SingleCellExperiment object for corresponding diploid observations
#' @param n_genes number of genes for simulation
#' @param n_cells_per_subclone number of cells in each subclone for simulation
#' @param subclone_cn_states matrix of unique copy number states for each subclone
#' 
#' @return matrix of simulated CNV observations
simSingleCnvExpr <- function(diploid_sim, n_genes, n_cells_per_subclone, subclone_cn_states) {
  # base gene means  & params for corresponding diploid observations 
  base_gene_means <- rowData(diploid_sim)$BaseGeneMean
  params <- diploid_sim@metadata$Params
  
  # simulate  CNV observations iteratively for each subclone
  expressions_counts <- list()
  for (i in seq_along(nCellsPerSubclone)) {
    params <- setParams(params, update = list(seed = sample(1:1000000, 1)))
    sim <- splatSimulate(params = params,
                         baseGeneMeans = base_gene_means,
                         batchCells = n_cells_per_subclone[i],
                         copyNumStates = subclone_cn_states[i,],
                         alpha = 1,
                         verbose=FALSE)
    
    expressions_counts[[i]] <- counts(sim)
  }
  
  # merge CNV observations matrix by vertically stack observations from each subclone
  observations <- do.call(cbind, expression_counts)
  
  return(observations)
}

#' simulate CNV observations for multiple cell types
#' 
#' @param diploid_sims list of SingleCellExperiment object for corresponding diploid observations, one for each group
#' @param n_genes number of genes for simulation
#' @param n_cells_per_subclone number of cells in each subclone for simulation
#' @param subclone_cn_states matrix of unique copy number states for each subclone
#' @param n_clusters number of groups
#' 
#' return matrix of simulated CNV observations
simGroupCnvExpr <- function(diploid_sims, n_genes, n_cells_per_subclone, subclone_cn_states, n_clusters) {
  
  n_subclones <- nrow(subclone_cn_states)
  
  params_list <- lapply(diploid_sims, function(sim) sim@metadata$Params)
  base_gene_means_list <- do.call(rbind, lapply(diploid_sims, function(sim) rowData(sim)$BaseGeneMean))
  
  index=1
  expression_counts <- list()
  for (i in c(1:n_clusters)) {
    for (j in c(1:n_subclones)) {
      params <- setParams(params_list[[i]], update = list(seed = sample(1:1000000, 1)))
      sim <- splatSimulate(params = params,
                           baseGeneMeans = base_gene_means_list[i,],
                           batchCells = n_cells_per_subclone[j],
                           copyNumStates = true_cn_states[j,],
                           alpha=1,
                           verbose=FALSE)
      expression_counts[[index]] <- counts(sim)
      index = index + 1
    }
  }
  
  # merge CNV observations matrix by vertically stack observations from each subclone & each group
  observations <- do.call(cbind, expression_counts)
  
  return(observations)
}

###############################################
#     Main script
###############################################

# parse arguments
argReader = commandArgs(trailingOnly=TRUE)
real_data_path <- argReader[1]  # path to input real data
cn_path <- argReader[2] # path to simulation data & picasso output 
sim_groups <- as.integer(argReader[3])  # whether to simulate groups (multiple cell types)
n_clusters <- as.integer(argReader[4])  # number of clusters / groups / cell types to simulate

# import ground truth copy number states from PICASSO
true_cn_states <- as.matrix(read.csv(paste0(cn_path, "true_states.txt"), header=FALSE, sep=" "))
true_cn_states <- unique(true_cn_states, MARGIN=1)  # keep one unique path per subclone 
n_genes <- dim(true_cn_states)[2]

# import number of cells/subclone from PICASSO
n_cells_per_subclone <- t(as.matrix(read.csv(paste0(cn_path, "num_cells.txt"), header=FALSE, sep=" ")))  # vector of num cells per in each subclone
n_cells_per_group <- sum(n_cells_per_subclone)

# --- Simulate diploid observations --- #
if (!sim_groups) {
  message("Simulating diploid observations for single cell type...")
  
  real_data_counts <- readSingleRealExpr(real_data_path, "MET31_diploid_observations.csv")  # counts from real data
  diploid_sim_results <- simSingleDiploidExpr(counts = real_data_counts, 
                                              n_genes = n_genes,
                                              n_cells = n_cells_per_group)
} else {
  message("Simulating diploid observations for multiple cell types... ")
  
  real_data_counts <- readGroupRealExpr(real_data_path, "normal_tissue_cleaned_data.csv")  # counts from real data
  diploid_sim_results <- simGroupDiploidExpr(counts = real_data_counts,
                                             n_genes = n_genes,
                                             n_cells_per_group = n_cells_per_group,
                                             n_clusters = n_clusters)
}

diploid_observations = diploid_sim_results$observations
diploid_sims <- diploid_sim_results$sims  # single SplatSimulate object if simulate single group, list of SplatSimulate object if simualte multiple groups

# --- Simulate CNV observations --- #
if (!sim_groups) {
  message("Simulating CNV observations for single cell type...")
  
  # simSingleCnvExpr <- function(diploid_sim, n_genes, n_cells_per_subclone, subclone_cn_states)
  observations <- simSingleCnvExpr(diploid_sim = diploid_sims,
                                   n_genes = n_genes,
                                   n_cells_per_subclone = n_cells_per_subclone,
                                   subclone_cn_states = true_cn_states)
  
} else {
  message("Simulating CNV observations for multiple cell types...")
  
  observations <- simGroupCnvExpr(diploid_sims = diploid_sims,
                                  n_genes = n_genes,
                                  n_cells_per_subclone = n_cells_per_subclone,
                                  subclone_cn_states = true_cn_states,
                                  n_clusters = n_clusters)
}


# --- write to output --- #
write.table(as.data.frame(diploid_observations), file = paste0(cn_path, "splat_diploid_observations.csv"),
            quote = F,
            eol = "\n",
            col.names = F,
            row.names = F,
            append = FALSE)

write.table(as.data.frame(observations), file = paste0(cn_path, "splat_observations.csv"), 
            quote = F, 
            eol = "\n", 
            col.names = F,
            row.names = F,
            append = FALSE)

# --- print simulation stats --- #
cat("----------------------------\n")
print(paste0('Dimensions of diploid data: ', dim(diploid_observations)))
printSimStats(diploid_observations, "Diploid reference observations")   # diploid observation
cat("\n\n")
print(paste0('Dimensions of altered data: ', dim(observations)))
printSimStats(observations, "Observations")  # regular observations
cat("----------------------------\n")

message("Splatter simulations finished...")