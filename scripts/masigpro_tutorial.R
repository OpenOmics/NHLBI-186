
library('maSigPro')

data("data.abiotic")
data("edesign.abiotic")


design <- make.design.matrix(edesign.abiotic, degree = 2,
                             )
View(design$dis)
fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
fit$i
fit$FDR
fit$dat
tstep <- T.fit(fit, step.method = "two.ways.backward", alfa = 0.05)
?T.fit
sigs2 <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
names(sigs)
sigs$summary
sigs$sig.genes$Control$sig.profiles
see.genes(sigs2$sig.genes$ColdvsControl, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9)
see.genes(sigs$sig.genes$Control, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9)
STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, edesign = edesign.abiotic)
PlotGroups (STMDE66, edesign = edesign.abiotic, show.fit = T,dis = design$dis, groups.vector = design$groups.vector)
suma2Venn(sigs$summary[, c(2:4)])


#~~~~~~~
fit <- p.vector(exps.deseq2, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
fit$i
fit$FDR
fit$dat
tstep <- T.fit(fit, step.method = "two.ways.backward", alfa = 0.05)
?T.fit
sigs2 <- get.siggenes(tstep, rsq = 0.6, vars = "groups")




see.genes(sigs2$sig.genes$ACAT1vsControl, show.fit = T, dis =design$dis,groups.vector = design$groups.vector,
          cluster.method="kmeans" ,cluster.data = 1, k = 2)
sigs$sig.genes$ACAT1vsControl$
