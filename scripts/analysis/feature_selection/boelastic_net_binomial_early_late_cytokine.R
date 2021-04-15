#!/usr/bin/env Rscript

# run bootstrapped elastic net/lasso on normalized expression counts

library("optparse")
library("glmnet")

options(stringsAsFactors = FALSE)

args_list <- list(
  make_option(c("-i", "--input_exp"), type = "character", help = "RData file containing expression matrix of normalized unscaled counts"),
  make_option(c("-l", "--lambda_seq"), type = "character", help = "RData file containing numeric vector of lambda sequence"),
  make_option(c("-a", "--alpha"), type = "numeric", help = "lasso: alpha = 1, ridge: alpha = 0, elastic net: alpha = 0.5"),
  make_option(c("-m", "--metadata"), type = "character", help = "csv file w/ cell meta data"),
  make_option(c("-b", "--num_bootstrap"), type = "integer", help = "number of bootstraps to perform"),
  make_option(c("-o", "--out_dir"), type = "character", help = "path to where output files should be written"),
  make_option(c("-n", "--analysis_name"), type = "character", help = "name of analysis used to name output files"),
  make_option(c("-r", "--replace"), type = "logical", default = TRUE, help = "TRUE: Sample w/ replacement; FALSE: Sample w/out replacement"),
  make_option(c("-s", "--subsample_fraction"), type = "numeric", default = 1, help = "subsample fraction; set to 1 when sampling w/ replacement")
)

args_parser <- OptionParser(option_list = args_list)
args <- parse_args(args_parser)

load(args$input_exp)  # loads cardio -> sparse matrix of normalized counts
load(args$lambda_seq) # loads lambda_seq -> numeric vector containing sequence of lambdas

cell_metadata <- read.csv(args$metadata)
cardio <- t(cardio)

cytokine_cells <- cell_metadata$protocol_final == "C"
cardio <- cardio[cytokine_cells,]

diff_day <- cell_metadata$day[cytokine_cells]

diff_day <- gsub('D12|D14', "early", diff_day)
diff_day <- gsub('D24|D26', "late", diff_day)
names(diff_day) <- cell_metadata$X[cytokine_cells]

diff_day <- as.factor(diff_day)

all_cells <- rownames(cardio)

bootstraps <- list()

get_cell_time <- function(cell_id, time_factor) {
  return(time_factor[cell_id])
}

for (i in 1:args$num_bootstrap) {
  print(i)
  cells <- sample(all_cells, size = floor(length(all_cells) * args$subsample_fraction), replace = args$replace)
  current_mat <- cardio[cells,]
  current_day <- sapply(cells, get_cell_time, time_factor = diff_day, simplify = TRUE)
  
  boots_fit <- glmnet(current_mat, current_day, family = "binomial", alpha = args$alpha, lambda = lambda_seq)
  bootstraps[[i]] <- boots_fit
}

save(bootstraps, file = paste0(args$out_dir, args$analysis_name, ".RData"))

