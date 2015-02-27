#!/usr/bin/env Rscript

suppressMessages(library(adegenet))
suppressMessages(library(ggplot2))

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
CANDIDATE.GROUPING.FIELD.NAMES = c( "model.category", "dispersal.model")
get.grouping.field.name = function(summary.df) {
    fieldnames = colnames(summary.df)
    for (field.name in CANDIDATE.GROUPING.FIELD.NAMES) {
        if (field.name %in% fieldnames) {
            return(field.name)
        }
    }
}

# Reports the levels/values in each of non-predictor fields.
data.regimes <- function(summary.df) {
    for (field.name in NON.PREDICTOR.FIELD.NAMES) {
        if (field.name %in% colnames(summary.df)) {
            cat(field.name, ": ", sort(unique(summary.df[,field.name])), "\n")
        }
    }
}

# Filters out data
filter.data <- function(summary.df,
                        filter.for.birth.rate=NULL,
                        filter.for.death.rate=NULL,
                        filter.for.dispersal.rate=NULL,
                        filter.for.trait.transition.rate=NULL) {
    summary.df.copy = summary.df
    if (!is.null(filter.for.birth.rate)) {
        summary.df.copy = subset(summary.df.copy, birth.rate==filter.for.birth.rate)
    }
    if (!is.null(filter.for.death.rate)) {
        summary.df.copy = subset(summary.df.copy, death.rate==filter.for.death.rate)
    }
    if (!is.null(filter.for.trait.transition.rate)) {
        summary.df.copy = subset(summary.df.copy, trait.transition.rate==filter.for.trait.transition.rate)
    }
    if (!is.null(filter.for.dispersal.rate)) {
        summary.df.copy = subset(summary.df.copy, dispersal.rate==filter.for.dispersal.rate)
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
create.group.and.predictors = function(summary.df,
                                       filter.for.birth.rate=NULL,
                                       filter.for.death.rate=NULL,
                                       filter.for.dispersal.rate=NULL,
                                       filter.for.trait.transition.rate=NULL) {
    source.df = filter.data(summary.df=summary.df,
                            filter.for.birth.rate=filter.for.birth.rate,
                            filter.for.death.rate=filter.for.death.rate,
                            filter.for.dispersal.rate=filter.for.dispersal.rate,
                            filter.for.trait.transition.rate=filter.for.trait.transition.rate)
    if (is.null(source.df)) {
        return(NULL)
    }
    source.df = na.omit(source.df)
    group = source.df[[get.grouping.field.name(summary.df)]]
    predictors = source.df[,!(names(source.df) %in% NON.PREDICTOR.FIELD.NAMES)]
    rv = list(
        group=group,
        predictors=predictors
        )
    rv
}

# Reports number of predictors in a data.frame
num.predictors <-function(summary.df) {
    x = create.group.and.predictors(summary.df=summary.df)
    if (is.null(x)) {
        return(NULL)
    } else {
        return(ncol(x$predictors))
    }
}

# Primary (back-end) workhorse function.
# Carries out the DAPC analysis, and packages the results.
calculate.dapc = function(predictors, group, n.pca, n.da, verbose.on.insufficient.groups=F) {
    num.groups = length(unique(group))
    if (num.groups < 2) {
        if (verbose.on.insufficient.groups) {
            warning(paste("Aborting: Only", num.groups, "groups:", group, "\n"))
        }
        return(NULL)
    }

    dapc.result = dapc(predictors, group, n.pca=n.pca, n.da=n.da)
    var.contr = data.frame(var=rownames(dapc.result$var.contr),
                        LD1=as.vector(dapc.result$var.contr)
                        )
    var.contr = var.contr[order(-var.contr$LD1),]
    model.prefs = data.frame(
                             group=group,                        # correct model
                             assign=dapc.result$assign,          # assigned model
                             # pp.model=dapc.result$posterior    # expands to `pp.model.{group}` for each value in `group`
                             data.frame(dapc.result$posterior)   # expands to one column for each model, with posterior probability of assignment to that model as value
                             )

    misassigns = model.prefs[model.prefs$group != model.prefs$assign,1]
    correct.assigns = model.prefs[model.prefs$group == model.prefs$assign,1]
    stopifnot(length(correct.assigns) + length(misassigns) == nrow(model.prefs))

    misassigns.prop = prop.table(table(misassigns))
    mean.prop.wrong.model.preferred = length(misassigns) / nrow(model.prefs)
    correct.assigns.prop = prop.table(table(correct.assigns))
    true.model.proportion.correctly.assigned = length(correct.assigns) / nrow(model.prefs)

    pps.of.correct.model = c()
    for (model.name in unique(model.prefs$group)) {
        pps.of.correct.model = c(pps.of.correct.model, model.prefs[model.prefs$group == model.name, model.name])
    }

    true.model.posterior.mean = sum(pps.of.correct.model) / length(pps.of.correct.model)

    rv = list(
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
calculate.true.model.proportion.correctly.assigned = function(predictors, group=group, n.pca, par, n.da) {
    x = calculate.dapc(
                       predictors=predictors,
                       group=group,
                       n.pca=n.pca,
                       n.da=n.da)
    return(x$true.model.proportion.correctly.assigned)
}

# Front end for analysis: (optionally) filters data, constructs groups and
# predictors, carries out DAPC analysis, and returns results.
analyze.dapc = function(
                        summary.df,
                        n.pca,
                        n.da,
                        filter.for.birth.rate=NULL,
                        filter.for.death.rate=NULL,
                        filter.for.dispersal.rate=NULL,
                        filter.for.trait.transition.rate=NULL,
                        verbose.on.insufficient.groups=NULL) {
    x = create.group.and.predictors(summary.df=summary.df,
                                     filter.for.birth.rate=filter.for.birth.rate,
                                     filter.for.death.rate=filter.for.death.rate,
                                     filter.for.dispersal.rate=filter.for.dispersal.rate,
                                     filter.for.trait.transition.rate=filter.for.trait.transition.rate
                                    )
    if (is.null(x)) {
        return(NULL)
    }

    rv = calculate.dapc(
            predictors=x$predictors,
            group=x$group,
            n.pca=n.pca,
            n.da=n.da,
            verbose.on.insufficient.groups=verbose.on.insufficient.groups)
}

# Predictor assessment.
# Returns efficacy of summary stastitics
assess.predictor.performance = function(summary.df,
                                        filter.for.birth.rate=NULL,
                                        filter.for.death.rate=NULL,
                                        filter.for.dispersal.rate=NULL,
                                        filter.for.trait.transition.rate=NULL) {
    groups.and.predictors = create.group.and.predictors(summary.df=summary.df,
                                    filter.for.birth.rate=filter.for.birth.rate,
                                    filter.for.death.rate=filter.for.death.rate,
                                    filter.for.dispersal.rate=filter.for.dispersal.rate,
                                    filter.for.trait.transition.rate=filter.for.trait.transition.rate)
    if (is.null(groups.and.predictors)) {
        return(NULL)
    }

    group = groups.and.predictors$group
    predictors = groups.and.predictors$predictors
    result = data.frame()
    for (n.pca in 2:ncol(predictors)) {
        n.da = 10
        # for (n.da in 1:ncol(predictors)) {
            x = calculate.dapc(predictors, group, n.pca, n.da)
            cat(paste(n.pca,
                        n.da,
                        x$true.model.posterior.mean,
                        x$true.model.proportion.correctly.assigned,
                        "\n",
                        sep="\t\t"
                        ))
            subresult = data.frame(
                        n.pca = n.pca,
                        n.da = n.da,
                        true.model.posterior.mean=x$true.model.posterior.mean,
                        true.model.proportion.correctly.assigned=x$true.model.proportion.correctly.assigned,
                        correct.assigns=as.list(x$correct.assigns.prop),
                        misassigns=as.list(x$misassigns.prop)
                        )
            result = rbind(result, subresult)
        # }
    }
    result
}

optimize.dapc.axes = function(predictors, group, penalized=T, verbose=F, n.pca.values=NULL, n.da.values=NULL) {
    if (is.null(n.pca.values)) {
        n.pca.values = 1:ncol(predictors)
    }
    if (is.null(n.da.values)) {
        n.da.values = 1:ncol(predictors)
    }
    best.score = NULL
    optima = list(n.pca=0, n.da=0)
    for (n.da in n.da.values) {
        for (n.pca in n.pca.values) {
            raw.score = calculate.true.model.proportion.correctly.assigned(
                   predictors=predictors,
                   group=group,
                   n.pca=n.pca,
                   n.da=n.da)
            if (!is.null(raw.score)) {
                if (penalized) {
                    # score =  (2 * n.pca) - log(2 * raw.score)
                    score =  (2 * n.pca) - log(2 * raw.score)
                } else {
                    score = - raw.score
                }
                if (verbose) {
                    # cat(paste("[current: ", optima$n.pca, " (score: ", best.score, ")] ", n.pca, ": raw score = ", raw.score, ", penalized score = ", score, " (", 2*n.pca, "-", log(2 * raw.score), ")\n", sep=""))
                    cat(paste("[current: ", optima$n.pca, " (score: ", best.score, ")] ", n.pca, ": raw score = ", raw.score, ", penalized score = ", score, "\n", sep=""))
                }
                if (is.null(best.score) || score < best.score) {
                    best.score = score
                    optima = list(n.pca=n.pca, n.da=n.da)
                }
            }
        }
    }
    return(optima)
}

optimize.dapc.axes.for.data.frame = function(summary.df, penalized=T, verbose=F, n.pca.values=NULL, n.da.values=NULL) {
    groups.and.predictors = create.group.and.predictors(summary.df=summary.df)
    if (is.null(groups.and.predictors)) {
        return(NULL)
    }
    group = groups.and.predictors$group
    predictors = groups.and.predictors$predictors
    return(optimize.dapc.axes(
                              predictors=predictors,
                              group=group,
                              penalized=penalized,
                              verbose=verbose,
                              n.pca.values=n.pca.values,
                              n.da.values=n.da.values))
}

# plot space, `parameter.space.df` is a data.frame returned by
# `analyze.parameter.space.discrete`, either directly or as loaded from a file.
plot.parameter.space.discrete = function(
                                         parameter.space.df,
                                         plot.type="scatter",
                                         characterization.schema="color-by-proportion-preferred",
                                         signficance.threshold=0.95
                                         ) {

    f1 = cut(
            parameter.space.df[["true.model.proportion.correctly.assigned"]],
            breaks=c(0.0, 0.5, signficance.threshold, 1.0),
            right=F,
            )
    parameter.space.df$true.model.proportion.correctly.assigned.factor = factor(f1, levels=rev(levels(f1)))
    f2 = cut(
            parameter.space.df[["true.model.posterior.mean"]],
            breaks=c(0.0, 0.5, signficance.threshold, 1.0),
            right=F,
            )
    parameter.space.df$true.model.posterior.mean.factor = factor(f2, levels=rev(levels(f2)))

    sweet.spot.data = factor(ifelse(parameter.space.df$true.model.proportion.correctly.assigned >= signficance.threshold & parameter.space.df$true.model.posterior.mean >= signficance.threshold,
                             "yes",
                             ifelse(parameter.space.df$true.model.proportion.correctly.assigned < signficance.threshold & parameter.space.df$true.model.posterior.mean < signficance.threshold,
                                    "no", "partial")
                             ))
    parameter.space.df$sweet.spot = factor(sweet.spot.data, levels=rev(levels(sweet.spot.data)))

    parameter.space.df$birth.rate.factor = factor(paste("b = ", parameter.space.df$birth.rate, sep=""))

    # parameter.space.df$birth.rate.factor = factor(paste("b=", parameter.space.df$birth.rate, sep=""))
    # parameter.space.df$death.rate.factor = factor(parameter.space.df$death.rate/parameter.space.df$birth.rate)
    ### requires labeler="parsed"
    # parameter.space.df$death.rate.factor = factor(factor(ifelse(parameter.space.df$death.rate==0, "0", paste("frac(b,", parameter.space.df$birth.rate/parameter.space.df$death.rate,")", sep="")))
    parameter.space.df$death.rate.factor = factor(paste("e = ",
                                                        ifelse(parameter.space.df$death.rate==0, "0", paste(format(round(parameter.space.df$death.rate/parameter.space.df$birth.rate, 2), nsmall=2), " x b", sep="")),
                                                        sep=""
                                                        ))

    p = ggplot(parameter.space.df, aes(trait.transition.rate, dispersal.rate))
    p = p + scale_x_log10() + scale_y_log10()
    posterior_legend_title = "Mean Posterior of True Model"
    prop_legend_title = "Mean Proportion True Model Preferred"
    sweet_spot_legend_title = ">0.95 Success"

    if (characterization.schema == "color-by-success") {
        color.settings = c("blue","dodgerblue", "orange") # third color for intermediate
    } else {
        color.settings = c("blue", "orange")
    }

    if (plot.type == "scatter") {
        if (characterization.schema == "color-by-success") {
            p = p + geom_point(aes(color=sweet.spot))
            p = p + scale_color_manual(values=color.settings, name=sweet_spot_legend_title)
            # p = p + scale_color_manual(values=c("dodgerblue", "orange", "red"), breaks=c("yes","partial","no"))
        } else  {
            if (characterization.schema == "color-by-posterior") {
                p = p + geom_point(aes(
                                    fill=true.model.posterior.mean.factor,
                                    shape=true.model.proportion.correctly.assigned.factor
                                    ))
                fill_legend = posterior_legend_title
                shape_legend = prop_legend_title
            } else if (characterization.schema == "color-by-proportion-preferred") {
                p = p + geom_point(aes(
                                    fill=true.model.proportion.correctly.assigned.factor,
                                    shape=true.model.posterior.mean.factor,
                                    ),
                                # size=2.5 # here, instead of in aes() because it is not a mapping
                                )
                fill_legend = prop_legend_title
                shape_legend = posterior_legend_title
            } else {
                stop(paste("Unrecognized characterization schema:'", characterization.schema, "'"))
            }
            p = p + scale_shape_manual(values=c(24, 25), name=shape_legend)
            p = p + scale_fill_manual(values=color.settings, name=fill_legend)
            # override to allow for fill-support shape
            # see: http://stackoverflow.com/questions/12488905/why-wont-the-ggplot2-legend-combine-manual-fill-and-scale-values
            # see: https://cloud.github.com/downloads/hadley/ggplot2/guide-col.pdf
            p = p + guides(fill=guide_legend(override.aes = list(shape = 21)))
            # p = p + scale_size(guide="none") # only needed if size is a mapping
        }
    } else if (plot.type == "tile") {

        # ggplot(df,aes(x = Var1,y = Var2,fill = factor(grp),alpha = z)) +
        #     geom_tile() +
        #     scale_fill_manual(values = c('red','blue'))
        if (characterization.schema == "color-by-success") {
            p = p + geom_tile(aes(fill=sweet.spot))
            fill_legend = sweet_spot_legend_title
        } else if (characterization.schema == "color-by-posterior") {
            p = p + geom_tile(aes(fill=true.model.posterior.mean.factor))
            fill_legend = posterior_legend_title
        } else if (characterization.schema == "color-by-proportion-preferred") {
            p = p + geom_tile(aes(fill=true.model.proportion.correctly.assigned.factor))
            fill_legend = prop_legend_title
        } else {
            stop(paste("Unrecognized characterization schema:'", characterization.schema, "'"))
        }
        p = p + scale_fill_manual(values=color.settings, name=fill_legend)
    } else {
        stop(paste("Unrecognized plot type:'", plot.type, "'"))
    }
    p = p + theme(
                  panel.border = element_rect(size=0.5, fill=NA, color="black", linetype='dotted'),
                  legend.position = "bottom"
                  # strip.background=element_rect(colour='black',  size=1)
                  )
    p = p + labs(x="Trait Transition Rate", y="Dispersal Rate")
    if (length(levels(parameter.space.df$death.rate.factor)) > 1) {
        # p = p + facet_grid(birth.rate.factor ~ death.rate.factor, labeller= label_parsed)
        p = p + facet_grid(birth.rate.factor ~ death.rate.factor)
    }
    p
}

analyze.parameter.space.discrete = function(summary.df, n.pca, n.da, verbose=NULL) {
    birth.rates = sort(unique(summary.df[,"birth.rate"]))
    if ("death.rate" %in% colnames(summary.df)) {
        death.rates = sort(unique(summary.df[,"death.rate"]))
    } else {
        death.rates = sort(unique(summary.df[,"extinction.rate"]))
    }
    dispersal.rates = sort(unique(summary.df[,"dispersal.rate"]))
    trait.transition.rates = sort(unique(summary.df[,"trait.transition.rate"]))
    result = data.frame()
    for (birth.rate in birth.rates) {
        for (death.rate in death.rates) {
            for (dispersal.rate in dispersal.rates) {
                for (trait.transition.rate in trait.transition.rates) {
                    x = analyze.dapc(summary.df=summary.df,
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
                        subresult = data.frame(
                                            birth.rate=birth.rate,
                                            death.rate=death.rate,
                                            dispersal.rate=dispersal.rate,
                                            trait.transition.rate=trait.transition.rate,
                                            true.model.proportion.correctly.assigned=x$true.model.proportion.correctly.assigned,
                                            true.model.posterior.mean=x$true.model.posterior.mean,
                                            correct.assigns=as.list(x$correct.assigns.prop),
                                            misassigns=as.list(x$misassigns.prop)
                                            )
                        result = rbind(result, subresult)
                    }
                }
            }
        }
    }
    result
}


# `target.summary.stats`
#   - data.frame of summary statistics calculated on empirical (or other) data
#     to be classified
# `training.summary.stats`
#   - data.frame of training data to be used to construct DAPC classification function
# `n.dapc.axes`
#   -   Set of number of axes retained in the principle component ('n.pca')
#       and discriminant analysis steps ('n.da').
#       This can be a string value:
#       - 'all'             : use maximum number of PC's available
#       - 'penalized-fit'   : use number of PC's that optimizes the proportion
#                               of correct classification balanced with number of
#                               PC's
#       - 'maximized-fit'   : use number of PC's that maximize the proportion
#                               of correct classification
#       Or a list with two named elements, 'n.pca' and 'n.da':
#       list(n.pca, n.da)   : use 'n.pca' for number of axes retained in
#                             principal component step and 'n.da' for number of
#                             axes retained in discriminant analysis step.
classify.data = function(target.summary.stats, training.summary.stats, n.dapc.axes) {
    training.data = create.group.and.predictors(training.summary.stats)
    predictors = training.data$predictors
    group = training.data$group
    if (n.dapc.axes == "ax"ll) {
    } else if (n.dapc.axes == "max") {
        n.pca = ncol(training.data$predictors)
        n.da = ncol(training.data$predictors)
    } else if ((n.dapc.axes == "penalized-fit") || (n.dapc.axes == "maximized-fit")) {
        if (n.dapc.axes == "penalized-fit") {
            penalized = T
        } else if (n.dapc.axes == "maximized-fit") {
            penalized = F
        }
        optima = optimize.dapc.axes(
                              predictors=predictors,
                              group=group,
                              penalized=penalized,
                              verbose=F)
        n.pca = optima$n.pca
        n.da = optima$n.da
    } else {
        n.pca = n.dapc.axes$n.pca
        n.da = n.dapc.axes$n.da
    }
    trained.model = calculate.dapc(
            predictors,
            group,
            n.pca=n.pca,
            n.da=n.da)
    target.data = create.group.and.predictors(target.summary.stats)
    pred.sup <- predict.dapc(trained.model$dapc.result, newdata=target.data$predictors)
    data.frame(pred.sup)
}

# `target.summary.stats.path`
#   -   Path to summary statistics calculated on empirical (or other) data to be classified
# `training.summary.stats.paths`
#   -   `list` of one or more paths to training data to be used to construct DAPC classification function
# `n.dapc.axes`
#   -   Set of number of axes retained in the principle component ('n.pca')
#       and discriminant analysis steps ('n.da').
#       This can be a string value:
#       - 'all'             : use maximum number of PC's available
#       - 'penalized-fit'   : use number of PC's that optimizes the proportion
#                               of correct classification balanced with number of
#                               PC's
#       - 'maximized-fit'   : use number of PC's that maximize the proportion
#                               of correct classification
#       Or a list with two named elements, 'n.pca' and 'n.da':
#       list(n.pca, n.da)   : use 'n.pca' for number of axes retained in
#                             principal component step and 'n.da' for number of
#                             axes retained in discriminant analysis step.
classify.data.from.files = function(
        target.summary.stats.path,
        training.summary.stats.paths,
        n.dapc.axes,
        output.path=NULL) {
    # training.summary.stats.paths <- list(...)
    # training.summary.stats <- list()
    # for (i in 1:length(training.summary.stats.paths)){
    #     training.summary.stats[[i]] <- read.csv(training.summary.stats.paths[[i]])
    # }
    target.summary.stats <- read.csv(target.summary.stats.path, header<-T)
    training.summary.stats.sets <- lapply(training.summary.stats.paths, read.csv, header=T)
    training.summary.stats.merged <- Reduce(function(x,y){rbind(x,y)}, training.summary.stats.sets)
    results = classify.data(
                            target.summary.stats=target.summary.stats,
                            training.summary.stats=training.summary.stats.merged,
                            n.dapc.axes=n.dapc.axes
                            )
    if (!is.null(output.path)) {
        write.csv(results, output.path, row.names=FALSE)
    }
    results
}

