#!/usr/bin/env Rscript

script.path <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}
script.dir = dirname(script.path())
args <- commandArgs(trailingOnly = T)
if (length(args) < 4) {
    cat("usage: archipelago-classify.R <NPCA> <NDA> <TARGET> <TRAINING> [TRAINING [TRAINING [...]]]\n")
    cat("\n")
    cat("where:\n")
    cat("\n")
    cat("<NPCA>          number of principle component axes to use\n")
    cat("                (specify 'NULL' for default of maximum).\n")
    cat("<NDA>           number of discriminant functions to use\n")
    cat("                (specify 'NULL' for default of 10).\n")
    cat("<TARGET>        is a CSV file containing the summary \n")
    cat("                statistics of the data to be classified.\n")
    cat("<TRAINING>      are CSV files containing the summary \n")
    cat("                statistics of the training data.\n")
    q()
}

n.pca = args[[1]]
if (n.pca == "NULL") {
    n.pca = NULL
} else {
    n.pca = as.numeric(n.pca)
}
n.da = args[[2]]
if (n.da == "NULL") {
    n.da = NULL
} else {
    n.da = as.numeric(n.da)
}
target.file = args[[3]]
training.files = args[c(-1,-2,-3)]

# cat(paste("n.pca = ", n.pca, "\n", sep=""))
# cat(paste("n.da = ", n.da, "\n", sep=""))
# cat(paste("target = ", target.file, "\n", sep=""))
# cat(paste("training = ", training.files, "\n", sep=""))

source(paste(script.dir, "..", "R", "analyze-dapc.R", sep="/"))
res = classify.data.from.files(
                               target.summary.stats.path=target.file,
                               training.summary.stats.paths=training.files,
                               n.pca=n.pca,
                               n.da=n.da)
write.csv(res, file="")

