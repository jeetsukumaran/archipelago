library(adegenet)

# These columns will be dropped from the training data set (they are
# typically parameters used to generate/simulate data).
NON.PREDICTOR.FIELD.NAMES <- c(
               "dispersal.model",
               "birth.rate",
               "death.rate",
               "extinction.rate",
               "dispersal.rate",
               "trait.transition.rate",
               "trait.evolution.rate",
               "ntax"
               )
# This column has the model label or category as the value.
GROUPING.FIELD.NAME = "dispersal.model"

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
                        filter.for.dispersal.rate=NULL,
                        filter.for.trait.transition.rate=NULL) {
    summary.df.copy = summary.df
    if (!is.null(filter.for.birth.rate)) {
        summary.df.copy = subset(summary.df.copy, birth.rate==filter.for.birth.rate)
    }
    if (!is.null(filter.for.trait.transition.rate)) {
        summary.df.copy = subset(summary.df.copy, trait.transition.rate==filter.for.trait.transition.rate)
    }
    if (!is.null(filter.for.dispersal.rate)) {
        summary.df.copy = subset(summary.df.copy, dispersal.rate==filter.for.dispersal.rate)
    }
    summary.df.copy
}

# Given a data.frame, returns a list with two named elements:
# `group`: data.frame
#   A data.frame consisting of a single-column, the grouping variable
# `predictors`: data.frame
#   A data.frame consisting of (just) the predictor variables.
create.group.and.predictors = function(summary.df,
                                       filter.for.birth.rate=NULL,
                                       filter.for.dispersal.rate=NULL,
                                       filter.for.trait.transition.rate=NULL) {
    source.df = filter.data(summary.df=summary.df,
                            filter.for.birth.rate=filter.for.birth.rate,
                            filter.for.dispersal.rate=filter.for.dispersal.rate,
                            filter.for.trait.transition.rate=filter.for.trait.transition.rate)
    source.df = na.omit(source.df)
    group = source.df[[GROUPING.FIELD.NAME]]
    predictors = source.df[,!(names(source.df) %in% NON.PREDICTOR.FIELD.NAMES)]
    rv = list(
        group=group,
        predictors=predictors
        )
    rv
}

# Primary workhorse function.
# Carries out the DAPC analysis, and packages the results.
# TODO: generalize from "constrained"/"unconstrained"
calculate.dapc = function(predictors, group, n.pca, n.da) {
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
                                        filter.for.dispersal.rate=NULL,
                                        filter.for.trait.transition.rate=NULL) {
    x = create.group.and.predictors.for.regime(summary.df=summary.df,
                                                filter.for.birth.rate=filter.for.birth.rate,
                                                filter.for.dispersal.rate=filter.for.dispersal.rate,
                                                filter.for.trait.transition.rate=filter.for.trait.transition.rate)
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
                        mean.prop.correct.model.preferred=x$mean.prop.correct.model.preferred
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
                        filter.for.dispersal.rate=NULL,
                        filter.for.trait.transition.rate=NULL) {
    x = create.group.and.predictors(summary.df=summary.df,
                                     filter.for.birth.rate=filter.for.birth.rate,
                                     filter.for.dispersal.rate=filter.for.dispersal.rate,
                                     filter.for.trait.transition.rate=filter.for.trait.transition.rate)
    group = x$group
    predictors = x$predictors
    rv = calculate.dapc(
            predictors=predictors,
            group=group,
            n.pca=n.pca,
            n.da=n.da)
    # rv
}

analyze.parameter.space.discrete = function(summary.df) {
    birth.rates = sort(unique(summary.df[,"birth.rate"]))
    dispersal.rates = sort(unique(summary.df[,"dispersal.rate"]))
    trait.transition.rates = sort(unique(summary.df[,"trait.transition.rate"]))
    result = data.frame()
    n.pca = 47
    n.da = 10
    for (birth.rate in birth.rates) {
        for (dispersal.rate in dispersal.rates) {
            for (trait.transition.rate in trait.transition.rates) {
                x = analyze.dapc(summary.df=summary.df,
                                n.pca=n.pca,
                                n.da=n.da,
                                filter.for.birth.rate=birth.rate,
                                filter.for.dispersal.rate=dispersal.rate,
                                filter.for.trait.transition.rate=trait.transition.rate
                                )
                i1 = dispersal.rate / trait.transition.rate
                cat(paste(
                            birth.rate=birth.rate,
                            dispersal.rate=dispersal.rate,
                            trait.transition.rate=trait.transition.rate,
                            i1,
                            x$mean.pp.of.correct.model,
                            x$mean.prop.correct.model.preferred,
                            "\n",
                            sep="\t\t"
                            ))
                subresult = data.frame(
                                       birth.rate=birth.rate,
                                       dispersal.rate=dispersal.rate,
                                       trait.transition.rate=trait.transition.rate,
                                       i1=i1,
                                       mean.pp.of.correct.model=x$mean.pp.of.correct.model,
                                       mean.prop.correct.model.preferred=x$mean.prop.correct.model.preferred
                                       )
                result = rbind(result, subresult)
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

