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
        cat(field.name, ": ", sort(unique(summary.df[,field.name])), "\n")
    }
}


filter.data <- function(summary.df,
                        filter.for.birth.rate=NULL,
                        filter.for.dispersal.rate=NULL,
                        filter.for.niche.evolution.prob=NULL) {
    summary.df.copy = summary.df
    if (!is.null(filter.for.birth.rate)) {
        summary.df.copy = subset(summary.df.copy, birth.rate==filter.for.birth.rate)
    }
    if (!is.null(filter.for.niche.evolution.prob)) {
        summary.df.copy = subset(summary.df.copy, niche.evolution.prob==filter.for.niche.evolution.prob)
    }
    if (!is.null(filter.for.dispersal.rate)) {
        summary.df.copy = subset(summary.df.copy, dispersal.rate==filter.for.dispersal.rate)
    }
    summary.df.copy
}

# Given a data.frame, returns a list with two named elements:
# `groups`: data.frame
#   A data.frame consisting of a single-column, the groupsing variable
# `predictors`: data.frame
#   A data.frame consisting of (just) the predictor variables.
create.groups.and.predictors = function(summary.df,
                                       filter.for.birth.rate=NULL,
                                       filter.for.dispersal.rate=NULL,
                                       filter.for.niche.evolution.prob=NULL) {
    source.df = filter.data(summary.df=summary.df,
                            filter.for.birth.rate=filter.for.birth.rate,
                            filter.for.dispersal.rate=filter.for.dispersal.rate,
                            filter.for.niche.evolution.prob=filter.for.niche.evolution.prob)
    source.df = na.omit(source.df)
    groups = source.df[[GROUPING.FIELD.NAME]]
    predictors = source.df[,!(names(source.df) %in% NON.PREDICTOR.FIELD.NAMES)]
    rv = list(
        groups=groups,
        predictors=predictors
        )
    rv
}

# Primary workhorse function.
# Carries out the DAPC analysis, and packages the results.
# TODO: generalize from "constrained"/"unconstrained"
calculate.dapc = function(predictors, groups, n.pca, n.da) {
    dapc.result = dapc(predictors, groups, n.pca=n.pca, n.da=n.da)
    var.contr = data.frame(var=rownames(dapc.result$var.contr),
                        LD1=as.vector(dapc.result$var.contr)
                        )
    var.contr = var.contr[order(-var.contr$LD1),]
    model.prefs = data.frame(groups=groups, pp.model=dapc.result$posterior)
    mean.pp.for.correct.model = mean(c(model.prefs[model.prefs$groups == "constrained",]$pp.model.constrained, model.prefs[model.prefs$groups == "unconstrained",]$pp.model.unconstrained))
    # debug1 = nrow(model.prefs[model.prefs$groups=="constrained" & model.prefs$pp.model.constrained>0.5,])
    # debug2 = nrow(model.prefs[model.prefs$groups=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    mean.count.correct.model.preferred = nrow(model.prefs[model.prefs$groups=="constrained" & model.prefs$pp.model.constrained>0.5,]) + nrow(model.prefs[model.prefs$groups=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    mean.prop.correct.model.preferred = mean.count.correct.model.preferred / nrow(model.prefs)
    rv = list(
              dapc.result=dapc.result,
              var.contr=var.contr,
              model.prefs=model.prefs,
              mean.pp.for.correct.model=mean.pp.for.correct.model,
              mean.count.correct.model.preferred=mean.count.correct.model.preferred,
              mean.prop.correct.model.preferred=mean.prop.correct.model.preferred
              )
    rv
}

# Predictor assessment.
# Returns efficacy of summary stastitics
assess.predictor.performance = function(summary.df,
                                        filter.for.birth.rate=NULL,
                                        filter.for.dispersal.rate=NULL,
                                        filter.for.niche.evolution.prob=NULL) {
    x = create.groups.and.predictors.for.regime(summary.df=summary.df,
                                                filter.for.birth.rate=filter.for.birth.rate,
                                                filter.for.dispersal.rate=filter.for.dispersal.rate,
                                                filter.for.niche.evolution.prob=filter.for.niche.evolution.prob)
    groups = x$groups
    predictors = x$predictors
    result = data.frame()
    for (n.pca in 2:ncol(predictors)) {
        n.da = 10
        # for (n.da in 1:ncol(predictors)) {
            x = calculate.dapc(predictors, groups, n.pca, n.da)
            cat(paste(n.pca,
                        n.da,
                        x$mean.pp.for.correct.model,
                        x$mean.prop.correct.model.preferred,
                        "\n",
                        sep="\t\t"
                        ))
            subresult = data.frame(
                        n.pca = n.pca,
                        n.da = n.da,
                        mean.pp.for.correct.model=x$mean.pp.for.correct.model,
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
                        filter.for.niche.evolution.prob=NULL) {
    x = create.groups.and.predictors(summary.df=summary.df,
                                     filter.for.birth.rate=filter.for.birth.rate,
                                     filter.for.dispersal.rate=filter.for.dispersal.rate,
                                     filter.for.niche.evolution.prob=filter.for.niche.evolution.prob)
    groups = x$groups
    predictors = x$predictors
    rv = calculate.dapc(
            predictors=predictors,
            groups=groups,
            n.pca=n.pca,
            n.da=n.da)
    # rv
}

classify.trees = function(path.to.trees.summary, path.to.simulation.summary) {
    summary.df <- read.csv(path.to.simulation.summary, header<-T)
    training.data = create.groups.and.predictors(summary.df)
    trained.model = calculate.dapc(
            training.data$predictors,
            training.data$groups,
            n.pca=ncol(training.data$predictors),
            n.da=10)
    trees.df <- read.csv(path.to.trees.summary, header<-T)
    to.classify.data = create.groups.and.predictors(trees.df)
    pred.sup <- predict.dapc(trained.model$dapc.result, newdata=to.classify.data$predictors)
}


run = function() {
    x1 = read.csv("data.csv", header=T)
    x2 = analyze.dapc(x1, 10, 10)
}
