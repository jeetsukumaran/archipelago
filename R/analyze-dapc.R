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

summary.df <- read.csv("data/summary.csv", header<-T)
# summary.b1 = summary.df[summary.df$birth.rate==0.03,]
# dispersal.model :  1 2
# birth.rate :  0.003 0.03
# death.rate :  0
# dispersal.rate :  0.03 0.15 0.3
# niche.evolution.prob :  0.03 0.15 0.5 1
# > qqq = summary.df[summary.df$birth.rate==0.03 & summary.df$dispersal.rate==0.3 & summary.df$niche.evolution.prob==1,]
# > x = create.data(qqq)
# > predictors = x$predictors
# > groups = x$groups

data.regimes <- function(summary.df) {
    for (sim.param in SIM.PARAMS) {
        cat(sim.param, ": ", sort(unique(summary.df[,sim.param])), "\n")
    }
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

assess.vars = function(predictors, groups) {
    result = data.frame()
    for (n.pc in 2:ncol(predictors)) {
        n.da = 10
        # for (n.da in 1:ncol(predictors)) {
            x = analyze.dapc(predictors, groups, n.pc, n.da)
            cat(paste(n.pc,
                        n.da,
                        x$mean.pp.for.correct.model,
                        x$mean.prop.correct.model.preferred,
                        "\n",
                        sep="\t\t"
                        ))
            subresult = data.frame(
                        n.pc = n.pc,
                        n.da = n.da,
                        mean.pp.for.correct.model=x$mean.pp.for.correct.model,
                        mean.prop.correct.model.preferred=x$mean.prop.correct.model.preferred
                        )
            result = rbind(result, subresult)
        # }
    }
    result
    #  z.pp = result[order(-result$mean.pp.for.correct.model,result$n.pc,result$n.da),]
    #  z.prop = result[order(-result$mean.prop.correct.model.preferred,result$n.pc,result$n.da),]
}

create.data.for.regime <- function(
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

create.groups.and.predictors.for.regime = function(
                                                   filter.for.birth.rate=NULL,
                                                   filter.for.dispersal.rate=NULL,
                                                   filter.for.niche.evolution.prob=NULL) {
    regime.summary = create.data.for.regime(
                                            filter.for.birth.rate=filter.for.birth.rate,
                                            filter.for.dispersal.rate=filter.for.dispersal.rate,
                                            filter.for.niche.evolution.prob=filter.for.niche.evolution.prob)
    x = create.groups.and.predictors(regime.summary)
    x
}

assess.vars.for.regime = function(
                   filter.for.birth.rate=NULL,
                   filter.for.dispersal.rate=NULL,
                   filter.for.niche.evolution.prob=NULL) {
    x = create.groups.and.predictors.for.regime(
                   filter.for.birth.rate=filter.for.birth.rate,
                   filter.for.dispersal.rate=filter.for.dispersal.rate,
                   filter.for.niche.evolution.prob=filter.for.niche.evolution.prob)
    groups = x$groups
    predictors = x$predictors
    rv = assess.vars(predictors, groups)
    rv
}

analyze.dapc.for.regime = function(
                   n.pc,
                   n.da,
                   filter.for.birth.rate=NULL,
                   filter.for.dispersal.rate=NULL,
                   filter.for.niche.evolution.prob=NULL) {
    x = create.groups.and.predictors.for.regime(
                   filter.for.birth.rate=filter.for.birth.rate,
                   filter.for.dispersal.rate=filter.for.dispersal.rate,
                   filter.for.niche.evolution.prob=filter.for.niche.evolution.prob)
    groups = x$groups
    predictors = x$predictors
    rv = analyze.dapc(predictors, groups, n.pc, n.da)
    rv
}
