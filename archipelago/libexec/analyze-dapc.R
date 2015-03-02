#!/usr/bin/env Rscript

# Imports {{{1
# ==============================================================================
suppressMessages(library(adegenet))
# }}}1

# Field Names Diagnostics: Discriminating between Model Name Fields, Parameter Fields and Summary Statistics Fields {{{1
# ==============================================================================

# These columns will be dropped from the training data set (they are
# typically parameters used to generate/simulate data).
NON.PREDICTOR.FIELD.NAMES <- c(
               "dispersal.model", # legacy
               "model.category",  # this is the key field
               "birth.rate",
               "death.rate",
               "extinction.rate",
               "dispersal.rate",
               "trait.transition.rate",
               "trait.evolution.rate",
               "num.focal.areas",
               "num.supplemental.areas",
               "ntax"
               )

# This column has the model label or category as the value.
CANDIDATE.GROUPING.FIELD.NAMES <- c( "model.category", "dispersal.model")
getGroupingFieldName <- function(summary.df) {
    fieldnames <- colnames(summary.df)
    for (field.name in CANDIDATE.GROUPING.FIELD.NAMES) {
        if (field.name %in% fieldnames) {
            return(field.name)
        }
    }
}
# }}}1

# Core Analytical Functions {{{1

# Given a data.frame, returns a list with two named elements:
#   `group`: data.frame
#       A data.frame consisting of a single-column, the grouping variable
#   `predictors`: data.frame
#       A data.frame consisting of (just) the predictor variables.
createGroupAndPredictors <- function(summary.df) {
    if (is.null(summary.df)) {
        return(NULL)
    }
    summary.df <- na.omit(summary.df)
    group <- summary.df[[getGroupingFieldName(summary.df)]]
    predictors <- summary.df[,!(names(summary.df) %in% NON.PREDICTOR.FIELD.NAMES)]
    rv <- list(
        group=group,
        predictors=predictors
        )
    rv
}

# Primary (back-end) workhorse function.
# Carries out the DAPC analysis, and packages the results.
calculateDAPC <- function(predictors, group, n.pca, n.da, verbose.on.insufficient.groups=F) {
    num.groups <- length(unique(group))
    if (num.groups < 2) {
        if (verbose.on.insufficient.groups) {
            warning(paste("Aborting: Only", num.groups, "groups:", group, "\n"))
        }
        return(NULL)
    }
    dapc.result <- dapc(predictors, group, n.pca=n.pca, n.da=n.da)
    var.contr <- data.frame(var=rownames(dapc.result$var.contr),
                        LD1=as.vector(dapc.result$var.contr)
                        )
    var.contr <- var.contr[order(-var.contr$LD1),]
    model.prefs <- data.frame(
                             group=group,                        # correct model
                             assign=dapc.result$assign,          # assigned model
                             # pp.model=dapc.result$posterior    # expands to `pp.model.{group}` for each value in `group`
                             data.frame(dapc.result$posterior)   # expands to one column for each model, with posterior probability of assignment to that model as value
                             )

    misassigns <- model.prefs[model.prefs$group != model.prefs$assign,1]
    correct.assigns <- model.prefs[model.prefs$group == model.prefs$assign,1]
    stopifnot(length(correct.assigns) + length(misassigns) == nrow(model.prefs))

    misassigns.prop <- prop.table(table(misassigns))
    mean.prop.wrong.model.preferred <- length(misassigns) / nrow(model.prefs)
    correct.assigns.prop <- prop.table(table(correct.assigns))
    true.model.proportion.correctly.assigned <- length(correct.assigns) / nrow(model.prefs)

    pps.of.correct.model <- c()
    for (model.name in unique(model.prefs$group)) {
        pps.of.correct.model <- c(pps.of.correct.model, model.prefs[model.prefs$group == model.name, model.name])
    }

    true.model.posterior.mean <- sum(pps.of.correct.model) / length(pps.of.correct.model)

    rv <- list(
              dapc.result=dapc.result,
              var.contr=var.contr,
              model.prefs=model.prefs,
              # true.model.posterior.mean=true.model.posterior.mean,
              # mean.count.correct.model.preferred=mean.count.correct.model.preferred,
              # true.model.proportion.correctly.assigned=true.model.proportion.correctly.assigned,
              correct.assigns=correct.assigns,
              correct.assigns.prop=correct.assigns.prop,
              misassigns=misassigns,
              misassigns.prop=misassigns.prop,
              mean.prop.wrong.model.preferred=mean.prop.wrong.model.preferred,
              true.model.proportion.correctly.assigned=true.model.proportion.correctly.assigned,
              true.model.posterior.mean=true.model.posterior.mean
              )
    rv
}

# Front end for analysis: (optionally) filters data, constructs groups and
# predictors, carries out DAPC analysis, and returns results.
analyzeDAPC <- function(summary.df, n.pca, n.da, verbose.on.insufficient.groups=NULL) {
    x <- createGroupAndPredictors(summary.df=summary.df)
    if (is.null(x)) {
        return(NULL)
    }
    rv <- calculateDAPC(
            predictors=x$predictors,
            group=x$group,
            n.pca=n.pca,
            n.da=n.da,
            verbose.on.insufficient.groups=verbose.on.insufficient.groups)
}

# Carries out multiple DAPC analysis with different numbers of PC axes
# retained, and selects a number of principal component axes to retain based on
# maximizing the proportion of correct classifications when the resulting DAPC
# function is reapplied to the training data. If `penalization.weight` > 0,
# then a penalty factor will be applied for each addition PC axis retained.
#
# - predictors           : data.frame of predictors
# - group                : data.frame of groups (model name)
# - penalization.weight  : penalization.weight  (0 = no penalty)
# - n.da                 : number of discriminant analysis axes to
#                          retain (NULL = one less than
#                          the number of groups)
# - n.pca.values         : vector of number of axes to try (if not given,
#                          1:MAX
# - verbose              : whether or not to dump progress
optimizeNumPCAxes <- function(
        predictors,
        group,
        penalization.weight=1.0,
        n.da=NULL,
        n.pca.values=NULL,
        verbose=F
        ) {
    if (is.null(n.da)) {
        n.da <- length(levels(group)) - 1
    }
    if (is.null(n.pca.values)) {
        n.pca.values <- 1:ncol(predictors)
    }
    max.n.pca.values <- max(n.pca.values)
    best.score <- NULL
    optima <- list(n.pca=0, n.da=0, true.model.proportion.correctly.assigned=0, best.score=0)
    raw.scores <- c()
    scores <- c()
    saved.n.pca.values <- c()
    saved.n.da.values <- c()
    for (n.pca in n.pca.values) {
        dapc.results <- calculateDAPC(
                           predictors=predictors,
                           group=group,
                           n.pca=n.pca,
                           n.da=n.da)
        raw.score <- dapc.results$true.model.proportion.correctly.assigned
        if (!is.null(raw.score)) {
            saved.n.pca.values <- c(saved.n.pca.values, n.pca)
            saved.n.da.values <- c(saved.n.da.values, n.da)
            raw.scores <- c(raw.scores, raw.score)
            score <-  (penalization.weight * (n.pca/max.n.pca.values)) - (raw.score)
            scores <- c(scores, score)
            if (verbose) {
                # cat(paste("[current: ", optima$n.pca, " (score: ", best.score, ")] ", n.pca, ": raw score <- ", raw.score, ", penalized score <- ", score, " (", 2*n.pca, "-", log(2 * raw.score), ")\n", sep=""))
                cat(paste("[current: ", optima$n.pca, " (raw: ", optima$true.model.proportion.correctly.assigned, ", penalized: ", best.score, ")] ", n.pca, ": raw score: ", raw.score, ", penalized score: ", score, "\n", sep=""))
            }
            if (is.null(best.score) || score < best.score) {
                best.score <- score
                # optima <- list(n.pca=n.pca, n.da=n.da)
                optima <- list(n.pca=n.pca, n.da=n.da, true.model.proportion.correctly.assigned=raw.score, best.score=best.score)
            }
        }
    }
    details <- list(n.pca=saved.n.pca.values, n.da=saved.n.da.values, true.model.proportion.correctly.assigned=raw.scores, score=scores)
    optima$details <- details
    return(optima)
}

# Carries out multiple DAPC analysis with different numbers of PC axes
# retained, and selects number of principal component axes to retain based on
# criteria.
#
# - predictors           : data.frame of predictors
# - group                : data.frame of groups (model name)
# - penalization.weight  : penalization.weight  (0 = no penalty)
# - n.da                 : number of discriminant analysis axes to
#                          retain (NULL = one less than
#                          the number of groups)
# - n.pca.values         : vector of number of axes to try (if not given,
#                          1:MAX
# - verbose              : whether or not to dump progress
optimizeNumPCAxesForDataFrame <- function(
        summary.df,
        penalization.weight=1.0,
        n.da=NULL,
        n.pca.values=NULL,
        verbose=F
        ) {
    groups.and.predictors <- createGroupAndPredictors(summary.df=summary.df)
    if (is.null(groups.and.predictors)) {
        return(NULL)
    }
    group <- groups.and.predictors$group
    predictors <- groups.and.predictors$predictors
    return(optimizeNumPCAxes(
                         predictors=predictors,
                         group=group,
                         penalization.weight=penalization.weight,
                         n.da=n.da,
                         n.pca.values=n.pca.values,
                         verbose=verbose,
                         ))
}

# Classifies target data.
#
# Constructs a DAPC function based on training data, and applies it to the
# target data to classify the generating model.
#
# - target.summary.stats
#       data.frame of summary statistics calculated on empirical (or other)
#       data to be classified
# - training.summary.stats
#       data.frame of training data to be used to construct DAPC classification
#       function
# - n.pca
#       Set of number of principal component axes to be retained for the
#       analysis.
#
#       This can be a numeric value specifying this directly or
#       a string:
#
#           - 'all'        : use maximum number of PC's available
#           - 'optimize'   : Optimize: try out different numbers of PC's and
#                            pick the one that yields the highest proportion
#                            of correct classifications when the resulting
#                            DAPC function is reapplied to the training data,
#                            with an a penalty factor (set by
#                            `penalization.weight`) to penalize over-fitting.
#
# - n.da
#       Number of discriminant analysis axes to retain (NULL = one less than
#       the number of groups)
#
# - n.pca.optimization.penalty.weight     : weight f0 = no penalty)
classifyData <- function(target.summary.stats,
                         training.summary.stats,
                         n.pca,
                         n.da=NULL,
                         n.pca.optimization.penalty.weight=1.0
                         ) {
    training.data <- createGroupAndPredictors(training.summary.stats)
    predictors <- training.data$predictors
    group <- training.data$group

    if (is.null(n.da)) {
        n.da <- length(levels(group)) - 1
    }
    if (!is.numeric(n.da)) {
        stop(paste("ERROR: Number of discriminant axes is not a number: '", n.da, "'", sep=""))
    } else if (n.da < 1) {
        stop(paste("ERROR: Number of discriminant axes retained is < 1: ", n.da, sep=""))
    }
    if (n.pca == "all") {
        n.pca <- ncol(training.data$predictors)
    } else if (n.pca == "optimize") {
        optima <- optimizeNumPCAxes(
                              predictors=predictors,
                              group=group,
                              penalization.weight=n.pca.optimization.penalty.weight,
                              n.da=n.da,
                              verbose=F,
                              )
        n.pca <- optima$n.pca
    } else if (!is.numeric(n.pca)) {
        stop(paste("ERROR: Number of principal component axes retained is not a number: '", n.pca, "'", sep=""))
    } else if (n.pca < 1) {
        stop(paste("ERROR: Number of principal component axes retained is < 1: '", n.pca, "'", sep=""))
    } else {
        n.pca <- n.pca
    }
    trained.model <- calculateDAPC(
            predictors,
            group,
            n.pca=n.pca,
            n.da=n.da)
    target.data <- createGroupAndPredictors(target.summary.stats)
    pred.sup <- predict.dapc(trained.model$dapc.result, newdata=target.data$predictors)
    results <- data.frame(pred.sup)
    results$n.pca <- n.pca
    results$n.da <- n.da
    results
}

# Classifies target data.
#
# Constructs a DAPC function based on training data, and applies it to the
# target data to classify the generating model.
#
# - target.summary.stats.path
#       Path to summary statistics calculated on empirical (or other) data to be classified
# - training.summary.stats.paths
#       `list` of one or more paths to training data to be used to construct DAPC classification function
# - n.pca
#       Set of number of principal component axes to be retained for the
#       analysis.
#
#       This can be a numeric value specifying this directly or
#       a string:
#
#           - 'all'        : use maximum number of PC's available
#           - 'optimize'   : Optimize: try out different numbers of PC's and
#                            pick the one that yields the highest proportion
#                            of correct classifications when the resulting
#                            DAPC function is reapplied to the training data,
#                            with an a penalty factor (set by
#                            `penalization.weight`) to penalize over-fitting.
#
# - n.da
#       Number of discriminant analysis axes to retain (NULL = one less than
#       the number of groups)
#
# - n.pca.optimization.penalty.weight     : weight f0 = no penalty)
classifyDataFromFiles <- function(
        target.summary.stats.path,
        training.summary.stats.paths,
        n.pca,
        n.da=NULL,
        n.pca.optimization.penalty.weight=1.0,
        output.path=NULL) {
    target.summary.stats <- read.csv(target.summary.stats.path, header<-T)
    training.summary.stats.sets <- lapply(training.summary.stats.paths, read.csv, header=T)
    training.summary.stats.merged <- Reduce(function(x,y){rbind(x,y)}, training.summary.stats.sets)
    results <- classifyData(
                            target.summary.stats=target.summary.stats,
                            training.summary.stats=training.summary.stats.merged,
                            n.pca=n.pca,
                            n.da=n.da,
                            n.pca.optimization.penalty.weight=n.pca.optimization.penalty.weight
                            )
    if (!is.null(output.path)) {
        write.csv(results, output.path, row.names=FALSE)
    }
    results
}

# }}}1

# Pre-Analysis Data Preparation and Diagnostics {{{1

# Reports the levels/values in each of non-predictor fields.
reportDataRegimes <- function(summary.df) {
    for (field.name in NON.PREDICTOR.FIELD.NAMES) {
        if (field.name %in% colnames(summary.df)) {
            cat(field.name, ": ", sort(unique(summary.df[,field.name])), "\n")
        }
    }
}

# Reports number of predictors in a data.frame
numPredictors <-function(summary.df) {
    x <- createGroupAndPredictors(summary.df=summary.df)
    if (is.null(x)) {
        return(NULL)
    } else {
        return(ncol(x$predictors))
    }
}

# Filters out data
filterData <- function(summary.df,
                        filter.for.birth.rate=NULL,
                        filter.for.death.rate=NULL,
                        filter.for.dispersal.rate=NULL,
                        filter.for.trait.transition.rate=NULL) {
    summary.df.copy <- summary.df
    if (!is.null(filter.for.birth.rate)) {
        summary.df.copy <- subset(summary.df.copy, birth.rate==filter.for.birth.rate)
    }
    if (!is.null(filter.for.death.rate)) {
        summary.df.copy <- subset(summary.df.copy, death.rate==filter.for.death.rate)
    }
    if (!is.null(filter.for.trait.transition.rate)) {
        summary.df.copy <- subset(summary.df.copy, trait.transition.rate==filter.for.trait.transition.rate)
    }
    if (!is.null(filter.for.dispersal.rate)) {
        summary.df.copy <- subset(summary.df.copy, dispersal.rate==filter.for.dispersal.rate)
    }
    # stopifnot(nrow(summary.df.copy) > 0)
    if (nrow(summary.df.copy) > 0) {
        return(summary.df.copy)
    } else {
        return(NULL)
    }
}

# }}}1

