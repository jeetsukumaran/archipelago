#!/usr/bin/env Rscript

suppressMessages(library(adegenet))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

# These columns will be dropped from the training data set (they are
# typically parameters used to generate/simulate data).
NON.PREDICTOR.FIELD.NAMES <- c(
               "dispersal.model",
               "model.category",
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

# Reports the levels/values in each of non-predictor fields.
reportDataRegimes <- function(summary.df) {
    for (field.name in NON.PREDICTOR.FIELD.NAMES) {
        if (field.name %in% colnames(summary.df)) {
            cat(field.name, ": ", sort(unique(summary.df[,field.name])), "\n")
        }
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

# Given a data.frame, returns a list with two named elements:
# `group`: data.frame
#   A data.frame consisting of a single-column, the grouping variable
# `predictors`: data.frame
#   A data.frame consisting of (just) the predictor variables.
createGroupAndPredictors <- function(summary.df,
                                       filter.for.birth.rate=NULL,
                                       filter.for.death.rate=NULL,
                                       filter.for.dispersal.rate=NULL,
                                       filter.for.trait.transition.rate=NULL) {
    source.df <- filterData(summary.df=summary.df,
                            filter.for.birth.rate=filter.for.birth.rate,
                            filter.for.death.rate=filter.for.death.rate,
                            filter.for.dispersal.rate=filter.for.dispersal.rate,
                            filter.for.trait.transition.rate=filter.for.trait.transition.rate)
    if (is.null(source.df)) {
        return(NULL)
    }
    source.df <- na.omit(source.df)
    group <- source.df[[getGroupingFieldName(summary.df)]]
    predictors <- source.df[,!(names(source.df) %in% NON.PREDICTOR.FIELD.NAMES)]
    rv <- list(
        group=group,
        predictors=predictors
        )
    rv
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

# Carries out DAPC analysis, and assesses
calculateTrueModelProportionCorrectlyAssigned <- function(predictors, group=group, n.pca, par, n.da) {
    x <- calculateDAPC(
                       predictors=predictors,
                       group=group,
                       n.pca=n.pca,
                       n.da=n.da)
    return(x$true.model.proportion.correctly.assigned)
}

# Front end for analysis: (optionally) filters data, constructs groups and
# predictors, carries out DAPC analysis, and returns results.
analyzeDAPC <- function(
                        summary.df,
                        n.pca,
                        n.da,
                        filter.for.birth.rate=NULL,
                        filter.for.death.rate=NULL,
                        filter.for.dispersal.rate=NULL,
                        filter.for.trait.transition.rate=NULL,
                        verbose.on.insufficient.groups=NULL) {
    x <- createGroupAndPredictors(summary.df=summary.df,
                                     filter.for.birth.rate=filter.for.birth.rate,
                                     filter.for.death.rate=filter.for.death.rate,
                                     filter.for.dispersal.rate=filter.for.dispersal.rate,
                                     filter.for.trait.transition.rate=filter.for.trait.transition.rate
                                    )
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

# Predictor assessment.
# Returns efficacy of summary stastitics
assessPredictorPerformance <- function(summary.df,
                                        filter.for.birth.rate=NULL,
                                        filter.for.death.rate=NULL,
                                        filter.for.dispersal.rate=NULL,
                                        filter.for.trait.transition.rate=NULL) {
    groups.and.predictors <- createGroupAndPredictors(summary.df=summary.df,
                                    filter.for.birth.rate=filter.for.birth.rate,
                                    filter.for.death.rate=filter.for.death.rate,
                                    filter.for.dispersal.rate=filter.for.dispersal.rate,
                                    filter.for.trait.transition.rate=filter.for.trait.transition.rate)
    if (is.null(groups.and.predictors)) {
        return(NULL)
    }

    group <- groups.and.predictors$group
    predictors <- groups.and.predictors$predictors
    result <- data.frame()
    for (n.pca in 2:ncol(predictors)) {
        n.da <- length(levels(group)) - 1
        # for (n.da in 1:ncol(predictors)) {
            x <- calculateDAPC(predictors, group, n.pca, n.da)
            cat(paste(n.pca,
                        n.da,
                        x$true.model.posterior.mean,
                        x$true.model.proportion.correctly.assigned,
                        "\n",
                        sep="\t\t"
                        ))
            subresult <- data.frame(
                        n.pca=n.pca,
                        n.da=n.da,
                        true.model.posterior.mean=x$true.model.posterior.mean,
                        true.model.proportion.correctly.assigned=x$true.model.proportion.correctly.assigned,
                        correct.assigns=as.list(x$correct.assigns.prop),
                        misassigns=as.list(x$misassigns.prop)
                        )
            result <- rbind(result, subresult)
        # }
    }
    result
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
        raw.score <- calculateTrueModelProportionCorrectlyAssigned(
                predictors=predictors,
                group=group,
                n.pca=n.pca,
                n.da=n.da)
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


# performance.df - data.frame with the following columns:
#
#   'true.model.proportion.correctly.assigned'
#   'true.model.posterior.mean'
#   'birth.rate'
#   'death.rate'
#   'dispersal.rate'
#   'trait.transition.rate'
#
plotPerformanceOverParameterSpace <- function(performance.df) {
    breaks = c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    # breaks = c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)
    # breaks = c(0.0, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.86, 0.9,  0.95, 0.99, 1.0)
    f1 <- cut(
            performance.df[["true.model.proportion.correctly.assigned"]],
            breaks=breaks,
            right=F,
            )
    performance.df$true.model.proportion.correctly.assigned.factor <- factor(f1, levels=levels(f1))
    f2 <- cut(
            performance.df[["true.model.posterior.mean"]],
            breaks=breaks,
            right=F,
            )
    performance.df$true.model.posterior.mean.factor <- factor(f2, levels=levels(f2))
    performance.df$birth.rate.factor <- factor(paste("b = ", performance.df$birth.rate, sep=""))
    performance.df$death.rate.factor <- factor(paste("e = ",
                                                        ifelse(performance.df$death.rate==0, "0", paste(format(round(performance.df$death.rate/performance.df$birth.rate, 2), nsmall=2), " x b", sep="")),
                                                        sep=""
                                                        ))
    p <- ggplot(performance.df, aes(trait.transition.rate, dispersal.rate))
    p <- p + scale_x_log10() + scale_y_log10()

    p <- p + geom_tile(aes(fill=true.model.proportion.correctly.assigned.factor))

    # palette = "Greys"
    # palette = "YlGn"
    palette = "RdYlBu"
    if (length(breaks) > 9) {
        color.fn <- colorRampPalette(brewer.pal(9, palette))
        p <- p + scale_fill_manual(values=color.fn(length(breaks)), name="")
    } else {
        # p <- p + scale_fill_brewer(type="seq", palette=palette, name="")
        p <- p + scale_fill_brewer(palette=palette, name="")
    }

    if (length(levels(performance.df$death.rate.factor)) > 1 && length(levels(performance.df$death.rate.factor)) > 1) {
        # p <- p + facet_grid(birth.rate.factor ~ death.rate.factor, labeller= label_parsed)
        p <- p + facet_grid(birth.rate.factor ~ death.rate.factor)
    } else if (length(levels(performance.df$death.rate.factor)) > 1) {
        p <- p + facet_wrap(~death.rate.factor)
    } else if (length(levels(performance.df$birth.rate.factor)) > 1) {
        p <- p + facet_wrap(~birth.rate.factor)
    }

    p <- p + theme(
                  panel.border=element_rect(size=0.5, fill=NA, color="black", linetype='dotted'),
                  legend.position="bottom"
                  # strip.background=element_rect(colour='black',  size=1)
                  )

    p
}
# performance.df - data.frame with the following columns:
#
#   'true.model.proportion.correctly.assigned'
#   'true.model.posterior.mean'
#   'birth.rate'
#   'death.rate'
#   'dispersal.rate'
#   'trait.transition.rate'
#
plotPerformanceOverParameterSpaceScaledtoDiversificationRate <- function(performance.df) {

    breaks = c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    f1 <- cut(performance.df[["true.model.proportion.correctly.assigned"]], breaks=breaks, right=F)
    performance.df$true.model.proportion.correctly.assigned.factor <- factor(f1, levels=levels(f1))
    f2 <- cut(performance.df[["true.model.posterior.mean"]], breaks=breaks, right=F)
    performance.df$true.model.posterior.mean.factor <- factor(f2, levels=levels(f2))

    performance.df$diversification.rate = performance.df$birth.rate - performance.df$death.rate
    performance.df$scaled.trait.transition.rate = performance.df$trait.transition.rate / performance.df$diversification.rate
    performance.df$scaled.dispersal.rate = performance.df$dispersal.rate / performance.df$diversification.rate

    p <- ggplot(performance.df, aes(scaled.trait.transition.rate, scaled.dispersal.rate))
    p <- p + scale_x_log10(breaks=c(0.01, 0.1, 1,10,100)) + scale_y_log10(breaks=c(0.01, 0.1, 1,10,100))
    p <- p + geom_point(aes(fill=true.model.proportion.correctly.assigned.factor), pch=21, size=3)

    # palette = "Greys"
    # palette = "YlGn"
    palette = "RdYlBu"
    if (length(breaks) > 9) {
        color.fn <- colorRampPalette(brewer.pal(9, palette))
        p <- p + scale_fill_manual(values=color.fn(length(breaks)), name="")
    } else {
        # p <- p + scale_fill_brewer(type="seq", palette=palette, name="")
        p <- p + scale_fill_brewer(palette=palette, name="")
    }

    p <- p + theme(
                  legend.position="bottom"
                  # strip.background=element_rect(colour='black',  size=1)
                  )

    p
}

# plot space, `performance.df` is a data.frame returned by
# `analyzePerformanceOverDiscretizedParameterSpace`, either directly or as loaded from a file.
plotPerformanceOverParameterSpace.old <- function(
                                         performance.df,
                                         plot.type="tile",
                                         characterization.schema="color-by-proportion-preferred",
                                         signficance.threshold=0.95
                                         ) {

    f1 <- cut(
            performance.df[["true.model.proportion.correctly.assigned"]],
            breaks=c(0.0, 0.5, signficance.threshold, 1.0),
            right=F,
            )
    performance.df$true.model.proportion.correctly.assigned.factor <- factor(f1, levels=rev(levels(f1)))
    f2 <- cut(
            performance.df[["true.model.posterior.mean"]],
            breaks=c(0.0, 0.5, signficance.threshold, 1.0),
            right=F,
            )
    performance.df$true.model.posterior.mean.factor <- factor(f2, levels=rev(levels(f2)))

    sweet.spot.data <- factor(ifelse(performance.df$true.model.proportion.correctly.assigned >= signficance.threshold & performance.df$true.model.posterior.mean >= signficance.threshold,
                             "yes",
                             ifelse(performance.df$true.model.proportion.correctly.assigned < signficance.threshold & performance.df$true.model.posterior.mean < signficance.threshold,
                                    "no", "partial")
                             ))
    performance.df$sweet.spot <- factor(sweet.spot.data, levels=rev(levels(sweet.spot.data)))

    performance.df$birth.rate.factor <- factor(paste("b <- ", performance.df$birth.rate, sep=""))

    # performance.df$birth.rate.factor <- factor(paste("b=", performance.df$birth.rate, sep=""))
    # performance.df$death.rate.factor <- factor(performance.df$death.rate/performance.df$birth.rate)
    ### requires labeler="parsed"
    # performance.df$death.rate.factor <- factor(factor(ifelse(performance.df$death.rate==0, "0", paste("frac(b,", performance.df$birth.rate/performance.df$death.rate,")", sep="")))
    performance.df$death.rate.factor <- factor(paste("e <- ",
                                                        ifelse(performance.df$death.rate==0, "0", paste(format(round(performance.df$death.rate/performance.df$birth.rate, 2), nsmall=2), " x b", sep="")),
                                                        sep=""
                                                        ))

    p <- ggplot(performance.df, aes(trait.transition.rate, dispersal.rate))
    p <- p + scale_x_log10() + scale_y_log10()
    posterior_legend_title <- "Mean Posterior of True Model"
    prop_legend_title <- "Mean Proportion True Model Preferred"
    sweet_spot_legend_title <- ">0.95 Success"

    if (characterization.schema == "color-by-success") {
        color.settings <- c("blue","dodgerblue", "orange") # third color for intermediate
    } else {
        color.settings <- c("blue", "dodgerblue", "orange")
    }

    if (plot.type == "scatter") {
        if (characterization.schema == "color-by-success") {
            p <- p + geom_point(aes(color=sweet.spot))
            p <- p + scale_color_manual(values=color.settings, name=sweet_spot_legend_title)
            # p <- p + scale_color_manual(values=c("dodgerblue", "orange", "red"), breaks=c("yes","partial","no"))
        } else  {
            if (characterization.schema == "color-by-posterior") {
                p <- p + geom_point(aes(
                                    fill=true.model.posterior.mean.factor,
                                    shape=true.model.proportion.correctly.assigned.factor
                                    ))
                fill_legend <- posterior_legend_title
                shape_legend <- prop_legend_title
            } else if (characterization.schema == "color-by-proportion-preferred") {
                p <- p + geom_point(aes(
                                    fill=true.model.proportion.correctly.assigned.factor,
                                    shape=true.model.posterior.mean.factor,
                                    ),
                                # size=2.5 # here, instead of in aes() because it is not a mapping
                                )
                fill_legend <- prop_legend_title
                shape_legend <- posterior_legend_title
            } else {
                stop(paste("Unrecognized characterization schema:'", characterization.schema, "'"))
            }
            p <- p + scale_shape_manual(values=c(24, 25), name=shape_legend)
            p <- p + scale_fill_manual(values=color.settings, name=fill_legend)
            # override to allow for fill-support shape
            # see: http://stackoverflow.com/questions/12488905/why-wont-the-ggplot2-legend-combine-manual-fill-and-scale-values
            # see: https://cloud.github.com/downloads/hadley/ggplot2/guide-col.pdf
            p <- p + guides(fill=guide_legend(override.aes <- list(shape <- 21)))
            # p <- p + scale_size(guide="none") # only needed if size is a mapping
        }
    } else if (plot.type == "tile") {

        # ggplot(df,aes(x=Var1,y=Var2,fill=factor(grp),alpha= z)) +
        #     geom_tile() +
        #     scale_fill_manual(values= c('red','blue'))
        if (characterization.schema == "color-by-success") {
            p <- p + geom_tile(aes(fill=sweet.spot))
            fill_legend <- sweet_spot_legend_title
        } else if (characterization.schema == "color-by-posterior") {
            p <- p + geom_tile(aes(fill=true.model.posterior.mean.factor))
            fill_legend <- posterior_legend_title
        } else if (characterization.schema == "color-by-proportion-preferred") {
            p <- p + geom_tile(aes(fill=true.model.proportion.correctly.assigned.factor))
            fill_legend <- prop_legend_title
        } else {
            stop(paste("Unrecognized characterization schema:'", characterization.schema, "'"))
        }
        p <- p + scale_fill_manual(values=color.settings, name=fill_legend)
    } else {
        stop(paste("Unrecognized plot type:'", plot.type, "'"))
    }
    p <- p + theme(
                  panel.border=element_rect(size=0.5, fill=NA, color="black", linetype='dotted'),
                  legend.position="bottom"
                  # strip.background=element_rect(colour='black',  size=1)
                  )
    p <- p + labs(x="Trait Transition Rate", y="Dispersal Rate")
    if (length(levels(performance.df$death.rate.factor)) > 1 && length(levels(performance.df$death.rate.factor)) > 1) {
        # p <- p + facet_grid(birth.rate.factor ~ death.rate.factor, labeller= label_parsed)
        p <- p + facet_grid(birth.rate.factor ~ death.rate.factor)
    } else if (length(levels(performance.df$death.rate.factor)) > 1) {
        p <- p + facet_wrap(~death.rate.factor)
    } else if (length(levels(performance.df$birth.rate.factor)) > 1) {
        p <- p + facet_wrap(~birth.rate.factor)
    }
    p
}

analyzePerformanceOverDiscretizedParameterSpace <- function(summary.df, n.pca, n.da, verbose=NULL) {
    birth.rates <- sort(unique(summary.df[,"birth.rate"]))
    if ("death.rate" %in% colnames(summary.df)) {
        death.rates <- sort(unique(summary.df[,"death.rate"]))
    } else {
        death.rates <- sort(unique(summary.df[,"extinction.rate"]))
    }
    dispersal.rates <- sort(unique(summary.df[,"dispersal.rate"]))
    trait.transition.rates <- sort(unique(summary.df[,"trait.transition.rate"]))
    result <- data.frame()
    for (birth.rate in birth.rates) {
        for (death.rate in death.rates) {
            for (dispersal.rate in dispersal.rates) {
                for (trait.transition.rate in trait.transition.rates) {
                    x <- analyzeDAPC(summary.df=summary.df,
                                     n.pca=n.pca,
                                     n.da=n.da,
                                     filter.for.birth.rate=birth.rate,
                                     filter.for.death.rate=death.rate,
                                     filter.for.dispersal.rate=dispersal.rate,
                                     filter.for.trait.transition.rate=trait.transition.rate,
                                     verbose.on.insufficient.groups=F
                                     )
                    if (is.null(x)) {
                        warning(paste("NULL result:",
                                  "birth.rate=",
                                  birth.rate,
                                  "death.rate=",
                                  death.rate,
                                  "dispersal.rate=",
                                  dispersal.rate,
                                  "trait.transition.rate=",
                                  trait.transition.rate,
                                  "\n",
                                  sep=","
                                  ))
                    } else {
                        if (!is.null(verbose) && verbose) {
                            cat(paste(
                                        birth.rate,
                                        death.rate,
                                        dispersal.rate,
                                        trait.transition.rate,
                                        x$true.model.posterior.mean,
                                        x$true.model.proportion.correctly.assigned,
                                        "\n",
                                        sep="\t\t"
                                        ))
                        }
                        subresult <- data.frame(
                                            birth.rate=birth.rate,
                                            death.rate=death.rate,
                                            dispersal.rate=dispersal.rate,
                                            trait.transition.rate=trait.transition.rate,
                                            true.model.proportion.correctly.assigned=x$true.model.proportion.correctly.assigned,
                                            true.model.posterior.mean=x$true.model.posterior.mean,
                                            correct.assigns=as.list(x$correct.assigns.prop),
                                            misassigns=as.list(x$misassigns.prop)
                                            )
                        result <- rbind(result, subresult)
                    }
                }
            }
        }
    }
    result
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

