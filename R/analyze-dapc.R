library(adegenet)
library(ggplot2)

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

num.predictors <-function(summary.df) {
    x = create.group.and.predictors(summary.df=summary.df)
    if (is.null(x)) {
        return(NULL)
    } else {
        return(ncol(x$predictors))
    }
}

# Primary workhorse function.
# Carries out the DAPC analysis, and packages the results.
# TODO: generalize from "constrained"/"unconstrained"
calculate.dapc = function(predictors, group, n.pca, n.da, verbose.on.insufficient.groups=NULL) {
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
    mean.prop.correct.model.preferred = length(correct.assigns) / nrow(model.prefs)

    pps.of.correct.model = c()
    for (model.name in unique(model.prefs$group)) {
        pps.of.correct.model = c(pps.of.correct.model, model.prefs[model.prefs$group == model.name, model.name])
    }

    mean.pp.of.correct.model = sum(pps.of.correct.model) / length(pps.of.correct.model)

    # mean.pp.of.correct.model = mean(c(model.prefs[model.prefs$group == "constrained",]$pp.model.constrained, model.prefs[model.prefs$group == "unconstrained",]$pp.model.unconstrained))
    # debug1 = nrow(model.prefs[model.prefs$group=="constrained" & model.prefs$pp.model.constrained>0.5,])
    # debug2 = nrow(model.prefs[model.prefs$group=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    # mean.count.correct.model.preferred = nrow(model.prefs[model.prefs$group=="constrained" & model.prefs$pp.model.constrained>0.5,]) + nrow(model.prefs[model.prefs$group=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    # mean.prop.correct.model.preferred = mean.count.correct.model.preferred / nrow(model.prefs)

    rv = list(
              dapc.result=dapc.result,
              var.contr=var.contr,
              model.prefs=model.prefs,
              # mean.pp.of.correct.model=mean.pp.of.correct.model,
              # mean.count.correct.model.preferred=mean.count.correct.model.preferred,
              # mean.prop.correct.model.preferred=mean.prop.correct.model.preferred,
              correct.assigns=correct.assigns,
              correct.assigns.prop=correct.assigns.prop,
              misassigns=misassigns,
              misassigns.prop=misassigns.prop,
              mean.prop.wrong.model.preferred=mean.prop.wrong.model.preferred,
              mean.prop.correct.model.preferred=mean.prop.correct.model.preferred,
              mean.pp.of.correct.model=mean.pp.of.correct.model
              )
    rv
}

# Predictor assessment.
# Returns efficacy of summary stastitics
assess.predictor.performance = function(summary.df,
                                        filter.for.birth.rate=NULL,
                                        filter.for.death.rate=NULL,
                                        filter.for.dispersal.rate=NULL,
                                        filter.for.trait.transition.rate=NULL) {
    x = create.group.and.predictors(summary.df=summary.df,
                                    filter.for.birth.rate=filter.for.birth.rate,
                                    filter.for.death.rate=filter.for.death.rate,
                                    filter.for.dispersal.rate=filter.for.dispersal.rate,
                                    filter.for.trait.transition.rate=filter.for.trait.transition.rate)
    if (is.null(x)) {
        return(NULL)
    }

    group = x$group
    predictors = x$predictors
    result = data.frame()
    for (n.pca in 2:ncol(predictors)) {
        n.da = 10
        # for (n.da in 1:ncol(predictors)) {
            x = calculate.dapc(predictors, group, n.pca, n.da)
            cat(paste(n.pca,
                        n.da,
                        x$mean.pp.of.correct.model,
                        x$mean.prop.correct.model.preferred,
                        "\n",
                        sep="\t\t"
                        ))
            subresult = data.frame(
                        n.pca = n.pca,
                        n.da = n.da,
                        mean.pp.of.correct.model=x$mean.pp.of.correct.model,
                        mean.prop.correct.model.preferred=x$mean.prop.correct.model.preferred,
                        correct.assigns=as.list(x$correct.assigns.prop),
                        misassigns=as.list(x$misassigns.prop)
                        )
            result = rbind(result, subresult)
        # }
    }
    result
}

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

# plot space, `parameter.space.df` is a data.frame returned by
# `analyze.parameter.space.discrete`, either directly or as loaded from a file.
plot.parameter.space.discrete = function(parameter.space.df, plot.type="scatter") {
    characterization_schema = "color-by-proportion-preferred"

    f1 = cut(
            parameter.space.df[["mean.prop.correct.model.preferred"]],
            breaks=c(0.0, 0.5, 0.9, 1.0),
            right=F,
            )
    parameter.space.df$mean.prop.correct.model.preferred.factor = factor(f1, levels=rev(levels(f1)))
    f2 = cut(
            parameter.space.df[["mean.pp.of.correct.model"]],
            breaks=c(0.0, 0.5, 0.9, 1.0),
            right=F,
            )
    parameter.space.df$mean.pp.of.correct.model.factor = factor(f2, levels=rev(levels(f2)))

    parameter.space.df$sweet.spot = factor(
                                           ifelse(parameter.space.df$mean.prop.correct.model.preferred >= 0.90 & parameter.space.df$mean.pp.of.correct.model >= 0.90,
                                                "yes",
                                                ifelse(parameter.space.df$mean.prop.correct.model.preferred < 0.90 & parameter.space.df$mean.pp.of.correct.model < 0.90,
                                                "no", "partial")
                                                )
                                           )
    parameter.space.df$birth.rate.factor = factor(parameter.space.df$birth.rate)
    # parameter.space.df$birth.rate.factor = factor(paste("b=", parameter.space.df$birth.rate, sep=""))

    # parameter.space.df$death.rate.factor = factor(parameter.space.df$death.rate/parameter.space.df$birth.rate)

    # parameter.space.df$death.rate.factor = factor(
    #         factor(ifelse(parameter.space.df$death.rate==0, "0", paste("frac(b,", parameter.space.df$birth.rate/parameter.space.df$death.rate,")", sep="")))

    parameter.space.df$death.rate.factor = factor(ifelse(parameter.space.df$death.rate==0, "0", paste(format(round(parameter.space.df$death.rate/parameter.space.df$birth.rate, 2), nsmall=2), "*", "b", sep="")))

    p = ggplot(parameter.space.df, aes(trait.transition.rate, dispersal.rate))
    p = p + scale_x_log10() + scale_y_log10()
    if (characterization_schema == "both") {
        p = p + geom_point(aes(color=sweet.spot))
        # p = p + scale_color_manual(values=c("dodgerblue", "orange", "red"), breaks=c("yes","partial","no"))
    } else  {
        posterior_legend_title = "Mean Posterior of True Model"
        prop_legend_title = "Mean Proportion True Model Preferred"
        if (characterization_schema == "color-by-posterior") {
            p = p + geom_point(aes(
                                fill=mean.pp.of.correct.model.factor,
                                shape=mean.prop.correct.model.preferred.factor
                                ))
            fill_legend = posterior_legend_title
            shape_legend = prop_legend_title
        } else if (characterization_schema == "color-by-proportion-preferred") {
            p = p + geom_point(aes(
                                fill=mean.prop.correct.model.preferred.factor,
                                shape=mean.pp.of.correct.model.factor,
                                ),
                               # size=2.5 # here, instead of in aes() because it is not a mapping
                               )
            fill_legend = prop_legend_title
            shape_legend = posterior_legend_title
        } else {
            stop(paste("Unrecognized characterization schema:'", characterization_schema, "'"))
        }
        p = p + scale_shape_manual(values=c(24, 25), name=shape_legend)
        p = p + scale_fill_manual(values=c("dodgerblue", "orange", "red"), name=fill_legend)
        # override to allow for fill-support shape
        # see: http://stackoverflow.com/questions/12488905/why-wont-the-ggplot2-legend-combine-manual-fill-and-scale-values
        # see: https://cloud.github.com/downloads/hadley/ggplot2/guide-col.pdf
        p = p + guides(fill=guide_legend(override.aes = list(shape = 21)))
        # p = p + scale_size(guide="none") # only needed if size is a mapping
    }
    p = p + theme(legend.position = "bottom")
    p = p + labs(x="Trait Transition Rate", y="Dispersal Rate")
    if (length(levels(parameter.space.df$death.rate.factor)) > 1) {
        p = p + facet_grid(birth.rate.factor ~ death.rate.factor, labeller= label_parsed)
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
                                        x$mean.pp.of.correct.model,
                                        x$mean.prop.correct.model.preferred,
                                        "\n",
                                        sep="\t\t"
                                        ))
                        }
                        subresult = data.frame(
                                            birth.rate=birth.rate,
                                            death.rate=death.rate,
                                            dispersal.rate=dispersal.rate,
                                            trait.transition.rate=trait.transition.rate,
                                            mean.prop.correct.model.preferred=x$mean.prop.correct.model.preferred,
                                            mean.pp.of.correct.model=x$mean.pp.of.correct.model,
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

classify.trees = function(path.to.trees.summary, path.to.simulation.summary) {
    summary.df <- read.csv(path.to.simulation.summary, header<-T)
    training.data = create.group.and.predictors(summary.df)
    trained.model = calculate.dapc(
            training.data$predictors,
            training.data$group,
            n.pca=ncol(training.data$predictors),
            n.da=10)
    trees.df <- read.csv(path.to.trees.summary, header<-T)
    to.classify.data = create.group.and.predictors(trees.df)
    pred.sup <- predict.dapc(trained.model$dapc.result, newdata=to.classify.data$predictors)
}

