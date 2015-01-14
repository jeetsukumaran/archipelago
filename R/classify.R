library(adegenet)
SIM.PARAMS <- c(
               "dispersal.model"
               )
EXPERIMENTAL.TO.DROP <- c(
        # "community.by.habitat.normalized.unweighted.mntd.obs.p.H0",
        # "community.by.habitat.normalized.unweighted.mntd.obs.p.H1",
        # "community.by.habitat.normalized.unweighted.mntd.obs.p.mean",
        # "community.by.habitat.normalized.unweighted.mntd.obs.p.var",
        # "community.by.habitat.normalized.unweighted.mpd.obs.p.H0",
        # "community.by.habitat.normalized.unweighted.mpd.obs.p.H1",
        # "community.by.habitat.normalized.unweighted.mpd.obs.p.mean",
        # "community.by.habitat.normalized.unweighted.mpd.obs.p.var",
        # "community.by.habitat.normalized.weighted.mntd.obs.p.H0",
        # "community.by.habitat.normalized.weighted.mntd.obs.p.H1",
        # "community.by.habitat.normalized.weighted.mntd.obs.p.mean",
        # "community.by.habitat.normalized.weighted.mntd.obs.p.var",
        # "community.by.habitat.normalized.weighted.mpd.obs.p.H0",
        # "community.by.habitat.normalized.weighted.mpd.obs.p.H1",
        # "community.by.habitat.normalized.weighted.mpd.obs.p.mean",
        # "community.by.habitat.normalized.weighted.mpd.obs.p.var",
        # "community.by.island.normalized.unweighted.mntd.obs.p.mean",
        # "community.by.island.normalized.unweighted.mntd.obs.p.var",
        # "community.by.island.normalized.unweighted.mpd.obs.p.mean",
        # "community.by.island.normalized.unweighted.mpd.obs.p.var",
        # "community.by.island.normalized.weighted.mntd.obs.p.mean",
        # "community.by.island.normalized.weighted.mntd.obs.p.var",
        # "community.by.island.normalized.weighted.mpd.obs.p.mean",
        # "community.by.island.normalized.weighted.mpd.obs.p.var"
)
RESULT.METADATA <- c(
                    )
COLS.TO.DROP <- c(SIM.PARAMS, RESULT.METADATA, EXPERIMENTAL.TO.DROP)

analyze.dapc = function(predictors, groups, n.pc, n.da) {
    dapc.result = dapc(predictors, groups, n.pc=n.pc, n.da=n.da)
    var.contr = data.frame(var=rownames(dapc.result$var.contr),
                        LD1=as.vector(dapc.result$var.contr)
                        )
    var.contr = var.contr[order(-var.contr$LD1),]
    model.prefs = data.frame(dispersal.model=groups,
                        pp.model=dapc.result$posterior)
    mean.pp.for.correct.model = mean(c(model.prefs[model.prefs$dispersal.model == "constrained",]$pp.model.constrained,
                                    model.prefs[model.prefs$dispersal.model == "unconstrained",]$pp.model.unconstrained))
    debug1 = nrow(model.prefs[model.prefs$dispersal.model=="constrained" & model.prefs$pp.model.constrained>0.5,])
    debug2 = nrow(model.prefs[model.prefs$dispersal.model=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    mean.count.correct.model.preferred = nrow(model.prefs[model.prefs$dispersal.model=="constrained" & model.prefs$pp.model.constrained>0.5,]) + nrow(model.prefs[model.prefs$dispersal.model=="unconstrained" & model.prefs$pp.model.unconstrained>0.5,])
    mean.prop.correct.model.preferred = mean.count.correct.model.preferred / nrow(model.prefs)
    rv = list(
              dapc.result=dapc.result,
              var.contr=var.contr,
              model.prefs=model.prefs,
              mean.pp.for.correct.model=mean.pp.for.correct.model,
              mean.count.correct.model.preferred=mean.count.correct.model.preferred,
              mean.prop.correct.model.preferred=mean.prop.correct.model.preferred,
              debug1=debug1,
              debug2=debug2
              )
    rv
}

create.groups.and.predictors = function(source.df) {
    source.df = na.omit(source.df)
    groups = source.df$dispersal.model
    predictors = source.df[,!(names(source.df) %in% COLS.TO.DROP)]
    rv = list(
        groups=groups,
        predictors=predictors
        )
    rv
}

classify.trees = function(path.to.trees.summary, path.to.simulation.summary) {
    summary.df <- read.csv(path.to.simulation.summary, header<-T)
    training.data = create.groups.and.predictors(summary.df)
    # print(nrow(training.data$predictors))
    # print(nrow(training.data$groups))
    # print(ncol(training.data$predictors))
    trained.model = analyze.dapc(
            training.data$predictors,
            training.data$groups,
            n.pc=ncol(training.data$predictors),
            n.da=10)

    trees.df <- read.csv(path.to.trees.summary, header<-T)
    to.classify.data = create.groups.and.predictors(trees.df)
    pred.sup <- predict.dapc(trained.model$dapc.result, newdata=to.classify.data$predictors)
}


run = function() {
    path.to.trees.summary = "/Users/jeet/Documents/Projects/Phyloinformatics/Archipelago/archipelago-studies/archipelago-study-west-indian-birds/reference/summaries/west_indies.habitat_simplified.summary.csv"
    path.to.simulation.summary = "/Users/jeet/Documents/Projects/Phyloinformatics/Archipelago/archipelago-studies/archipelago-study-west-indian-birds-analyses/x1/experiment/data/summary.csv"
    t = classify.trees(
            path.to.trees.summary=path.to.trees.summary,
            path.to.simulation.summary=path.to.simulation.summary)
    t
}
