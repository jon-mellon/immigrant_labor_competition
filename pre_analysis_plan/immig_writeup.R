## ----setup, cache=TRUE, echo=FALSE, message=FALSE, warning=FALSE---------
# setwd("C:/Users/Jon/Dropbox/immig_exp")
# rm(list = ls())
# setwd("C:/Users/Jon/Documents/immigrant_labor_competition/pre_analysis_plan")
Sys.setenv(TEXINPUTS=getwd(),
           BIBINPUTS=getwd(),
           BSTINPUTS=getwd())
library(MASS)
library(haven)
library(knitr)

# Function that fixes a common read problem with the Haven package
labelDataset <- function(data) {
  correctLabel <- function(x) {
    if(!is.null(attributes(x)$labels)) {
      class(attributes(x)$labels) <- typeof(x)
    }
    return(x)
  }
  for(i in colnames(data)) {
    data[, i] <- correctLabel(data[, i])
  }
  return(data)
}


### Function adjusted from ocME function in the erer package. The original function did not support interaction terms

olMFX <- function(w) {
  rev.dum <- T
  digits <- 3
  lev <- w$lev
  J <- length(lev)
  x.name <- attr(x = w$terms, which = "term.labels")
  # x2 <- w$model[, x.name]
  x2 <- w$model[, -1]
  ww <- paste("~ 1", paste("+", x.name, collapse = " "), collapse = " ")
  x <- model.matrix(as.formula(ww), data = x2)[, -1]
  x.bar <- as.matrix(colMeans(x))
  b.est <- as.matrix(coef(w))
  K <- nrow(b.est)
  xb <- t(x.bar) %*% b.est
  z <- c(-10^6, w$zeta, 10^6)
  pfun <- switch(w$method, probit = pnorm, logistic = plogis)
  dfun <- switch(w$method, probit = dnorm, logistic = dlogis)
  V2 <- vcov(w)
  V3 <- rbind(cbind(V2, 0, 0), 0, 0)
  ind <- c(1:K, nrow(V3) - 1, (K + 1):(K + J - 1), nrow(V3))
  V4 <- V3[ind, ]
  V5 <- V4[, ind]
  f.xb <- dfun(z[1:J] - xb) - dfun(z[2:(J + 1)] - xb)
  me <- b.est %*% matrix(data = f.xb, nrow = 1)
  colnames(me) <- paste("effect", lev, sep = ".")
  se <- matrix(0, nrow = K, ncol = J)
  for (j in 1:J) {
    u1 <- c(z[j] - xb)
    u2 <- c(z[j + 1] - xb)
    if (w$method == "probit") {
      s1 <- -u1
      s2 <- -u2
    }
    else {
      s1 <- 1 - 2 * pfun(u1)
      s2 <- 1 - 2 * pfun(u2)
    }
    d1 <- dfun(u1) * (diag(1, K, K) - s1 * (b.est %*% t(x.bar)))
    d2 <- -1 * dfun(u2) * (diag(1, K, K) - s2 * (b.est %*% 
                                                   t(x.bar)))
    q1 <- dfun(u1) * s1 * b.est
    q2 <- -1 * dfun(u2) * s2 * b.est
    dr <- cbind(d1 + d2, q1, q2)
    V <- V5[c(1:K, K + j, K + j + 1), c(1:K, K + j, K + j + 
                                          1)]
    cova <- dr %*% V %*% t(dr)
    se[, j] <- sqrt(diag(cova))
  }
  colnames(se) <- paste("SE", lev, sep = ".")
  rownames(se) <- colnames(x)
  if (rev.dum) {
    for (k in 1:K) {
      if (identical(sort(unique(x[, k])), c(0, 1))) {
        for (j in 1:J) {
          x.d1 <- x.bar
          x.d1[k, 1] <- 1
          x.d0 <- x.bar
          x.d0[k, 1] <- 0
          ua1 <- z[j] - t(x.d1) %*% b.est
          ub1 <- z[j + 1] - t(x.d1) %*% b.est
          ua0 <- z[j] - t(x.d0) %*% b.est
          ub0 <- z[j + 1] - t(x.d0) %*% b.est
          me[k, j] <- pfun(ub1) - pfun(ua1) - (pfun(ub0) - 
                                                 pfun(ua0))
          d1 <- (dfun(ua1) - dfun(ub1)) %*% t(x.d1) - 
            (dfun(ua0) - dfun(ub0)) %*% t(x.d0)
          q1 <- -dfun(ua1) + dfun(ua0)
          q2 <- dfun(ub1) - dfun(ub0)
          dr <- cbind(d1, q1, q2)
          V <- V5[c(1:K, K + j, K + j + 1), c(1:K, K + 
                                                j, K + j + 1)]
          se[k, j] <- sqrt(c(dr %*% V %*% t(dr)))
        }
      }
    }
  }
  t.value <- me/se
  p.value <- 2 * (1 - pt(abs(t.value), w$df.residual))
  out <- list()
  for (j in 1:J) {
    out[[j]] <- round(cbind(effect = me[, j], error = se[,j], t.value = t.value[, j], p.value = p.value[, j]), 
                      digits)
  }
  out[[J + 1]] <- round(me, digits)
  names(out) <- paste("ME", c(lev, "all"), sep = ".")
  result <- listn(w, out)
  class(result) <- "ocME"
  return(result)
}


immig <- read_spss("BES2016_wave7.sav")
bes.data <- read_dta("bes_immig.dta")
bes.data <- labelDataset(bes.data)
jobs <- read.csv("job_assignments_final.csv",
                 stringsAsFactors = FALSE)

# jobs[match(immig$id, jobs$id), c("job.1", "job.2")]==immig[, c("Job1", "Job2")]
orig.occs <- read.csv("occsmall2_coded.csv", stringsAsFactors = FALSE)
immig$own.job <- orig.occs[match(immig$id, orig.occs$id), "exptext"]

immig$dv <- as_factor(immig$immigExpDV)
# immig$dv.dk <- as_factor(immig$immigExpDV)
immig$dv[immig$dv=="Don't know"] <- NA
immig <- immig[!is.na(immig$own.job), ]
immig <- merge(immig, bes.data, by = "id", all.x = TRUE)

immig <- merge(immig, jobs, by= "id")
immig$treatment <- immig$own.job==immig$job.2
immig$immigManipCheck[immig$immigManipCheck==997] <- NA
immig$immigManipCheck2[immig$immigManipCheck2==997] <- NA
rm(bes.data)

manip.check.table <- table(immig$immigManipCheck, immig$treatment)
colnames(manip.check.table) <- c("Control groups", "Own-job Treatment")
manip.check.table <- round(prop.table(manip.check.table, 2) * 100, 1)

library(Hmisc, quietly = TRUE)

# chisq.test(y = immig$immigManipCheck, x = immig$treatment)


low.skill.jobs <- c("waiters", "drivers", "receptionists",
                    "shop assistants", "carers", "cleaners")
high.skill.jobs <- c("senior managers", "doctors", 
                     "lawyers", "lecturers", "engineers",
                     "programmers")

immig$skill.1 <- NA
immig$skill.1[immig$job.1 %in% low.skill.jobs] <- "Low"
immig$skill.1[immig$job.1 %in% high.skill.jobs] <- "High"

immig$skill.2 <- NA
immig$skill.2[immig$job.2 %in% low.skill.jobs & !immig$treatment] <- "Low"
immig$skill.2[immig$job.2 %in% high.skill.jobs & !immig$treatment] <- "High"
immig$skill.2[immig$treatment] <- "treatment"

immig$group <- NA
immig$group[immig$skill.1=="Low" & immig$skill.2=="Low"] <- "Low"
immig$group[immig$skill.1=="High" & immig$skill.2=="High"] <- "High"
immig$group[(immig$skill.1=="High" & immig$skill.2=="Low") | 
              (immig$skill.1=="Low" & immig$skill.2=="High")] <- "Mixed"

immig$group[immig$skill.1=="Low" & immig$treatment] <- "low-treat"
immig$group[immig$skill.1=="High" & immig$treatment] <- "high-treat"

immig$edlevel1 <- NA
immig$edlevel1[immig$edlevel==0] <- 1
immig$edlevel1[immig$edlevel==1] <- 2
immig$edlevel1[immig$edlevel==2] <- 2
immig$edlevel1[immig$edlevel==3] <- 3
immig$edlevel1[immig$edlevel==4] <- 4
immig$edlevel1[immig$edlevel==5] <- 4
immig$edlevel1 <- factor(immig$edlevel1, labels = c("None", "GCSE", "A-Level", "Degree"))


immig$edlevel2 <- NA
immig$edlevel2[immig$edlevel1 %in% c("None", "GCSE")] <- "Low Ed"
immig$edlevel2[immig$edlevel1 %in% c("A-Level")] <- "Med Ed"
immig$edlevel2[immig$edlevel1 %in% c("Degree")] <- "High Ed"

library(ggplot2)
library(reshape)

immig$dv <- factor(immig$dv)
immig$edlevel2 <- as.factor(immig$edlevel2)

immig$dv <- factor(immig$dv)
immig$Treatment <- as.numeric(immig$treatment)

immig.replicate <- immig[immig$group %in% c("High", "Mixed", "Low"), ]
immig.replicate$dv.dk <- immig.replicate$dv



nice.names <- matrix(c("edlevel2Low Ed" , "Low education", 
                       "edlevel2Med Ed", "Medium education",
                       "skill.1Low", "Job 1: low skill", 
                       "skill.2Low", "Job 2: Low skill", 
                       "skill.2treatment", "Job 2: own job", 
                       "edlevel2Low Ed:skill.2Low", "Low ed. * low skill job", 
                       "edlevel2Med Ed:skill.2Low", "Medium ed. * low skill job", 
                       "edlevel2Low Ed:skill.2treatment", "Low ed. * own job", 
                       "edlevel2Med Ed:skill.2treatment", "Medium ed. * own job", 
                       "groupMixed", "Treatment: low/high", 
                       "groupHigh", "Treatment: both high", 
                       "edlevel2Low Ed:groupMixed", "Low ed. * low/high", 
                       "groupMixed:edlevel2Low Ed", "Low ed. * low/high", 
                       "edlevel2Med Ed:groupMixed", "Medium ed. * low/high",
                       "edlevel2Low Ed:groupHigh", "Low ed. * both high", 
                       "edlevel2Med Ed:groupHigh", "Medium ed. * both high",
                       "groupLow:edlevel2Low Ed" , "Low ed. * both low",
                       "groupLow:edlevel2Med Ed", "Medium ed. * both low",
                       "groupMixed:edlevel2Med Ed", "Medium ed. * low/high",
                       "edlevel2Low Ed:groupLow", "Low ed. * both low",
                       "edlevel2Med Ed:groupLow", "Medium ed. * both low",
                       "groupLow", "Treatment: both low", 
                       "grouphigh-treat", "Treatment: high/ own job",
                       "grouplow-treat", "Treatment: low/ own job") ,
                     byrow = TRUE, ncol = 2)
colnames(nice.names) <- c("model", "nice")
nice.names <- data.frame(nice.names, stringsAsFactors = FALSE)
nice.order <- c("Treatment: high/ own job",
                "Treatment: low/ own job",
                "Treatment: both high",
                "Treatment: low/high",
                "Treatment: both low",
                "Job 2: high skill",
                "Job 2: own job", 
                "Job 2: Low skill", 
                "Job 1: high skill",
                "Job 1: low skill", 
                "High education",
                "Medium education",
                "Low education",
                "Medium ed. * both high",
                "Low ed. * both high", 
                "Medium ed. * low/high",
                "Low ed. * low/high", 
                "Medium ed. * both low",
                "Low ed. * both low",
                "Medium ed. * own job",
                "Low ed. * own job", 
                "Medium ed. * low skill job",
                "Low ed. * low skill job")



makeNice <- function(x) {
  output <- nice.names$nice[match(x, nice.names$model)]
  return(output)
}
makeNiceFactor <- function(x) {
  output <- factor(x, levels = nice.order[length(nice.order):1])
  output <- factor(output)
  return(output)
}

immig.replicate$group <- factor(immig.replicate$group, 
                                levels = c("High", "Mixed", "Low"))


## ----assignments, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE----
all.assign <- table(immig$group)
assign.table <- data.frame(matrix(c("Low skilled/own occupation", "25", all.assign["low-treat"], 
                                    "High skilled/own occupation", "25", all.assign["high-treat"], 
                                    "Both high skilled", "12.5",  all.assign["High"], 
                                    "Both low skilled", "12.5", all.assign["Low"], 
                                    "Low skilled and high skilled", "25", all.assign["Mixed"]), byrow = TRUE, ncol = 3), stringsAsFactors = FALSE)

colnames(assign.table) <- c("Assignment", "% in limit", "Cases assigned")
assign.table$`% assigned` <- round(prop.table(as.numeric(assign.table$`Cases assigned`)) * 100, 1)
library(xtable)
xtable(assign.table, caption = "Experimental assignments", label = "table:assignment")

## ----mainLabMarketDisplay, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE----
library(pander)
library(texreg)
library(erer)
library(MASS)

immig.main.reg <- polr(data = immig, formula = dv~ edlevel2 + 
                         edlevel2:skill.2 + skill.1 + skill.2)
texreg(immig.main.reg, caption = "Ordered logistic regression predicting acceptance of immigrants. (Reference categories: high education, Job 1: high, skill, and Job 2: high skill)", 
       label = "table:mainLabMarket", 
       custom.coef.names = makeNice(names(immig.main.reg$coefficients)))

## ----mainLabMarketMFX, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE, fig.cap = "Effects of labor market competition treatment on probability of `strongly disagreeing' that more immigrants should be allowed to come to Britain"----
mfx.main <- olMFX(immig.main.reg)

mfx.df <- data.frame(mfx= mfx.main$out$`ME.Strongly disagree`[, "effect"], 
                     se = mfx.main$out$`ME.Strongly disagree`[, "error"])
mfx.df$u.ci <- mfx.df$mfx + (1.96 * mfx.df$se)
mfx.df$l.ci <- mfx.df$mfx - (1.96 * mfx.df$se)
mfx.df$var <- makeNice(rownames(mfx.df))
mfx.df <- rbind(mfx.df, NA, NA, NA)

mfx.df$var[is.na(mfx.df$var)] <- c("Job 1: high skill", "High education", 
                                   "Job 2: high skill")
mfx.df$var <- makeNiceFactor(mfx.df$var)
mfx.df$mfx[is.na(mfx.df$mfx)] <- 0

ggplot(data = mfx.df, aes(x = var, y = mfx)) + geom_point() + 
  geom_errorbar(aes(ymax = u.ci, ymin = l.ci), width = 0.25) +
  geom_hline(aes(yintercept=0)) +  theme_bw() + 
  xlab("") + ylab("Marginal effect of strongly disagreeing with more immigrants") + 
  coord_flip()

## ----manipCheckSimple, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE----

library(pander)
library(texreg)

immig.manip.compete <- polr(data = immig,
                            formula = factor(immigManipCheck)~ Treatment)
immig.manip.culture <- polr(data = immig, 
                            formula = factor(immigManipCheck2)~ Treatment)

texreg(list(immig.manip.compete, immig.manip.culture), 
       caption = "Manipulation checks. Ordered logistic regression predicting worry about job prospects.", 
       custom.model.names = c("Job prospects", "Cultural threat"), 
       label = "table:manipulation_checks")

## ----mainReplication, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE, fig.cap="Effects of immigrant skill levels on probability of `strongly disagreeing' that more immigrants should be allowed to come to Britain (excludes `own job' treatment group)"----
reg.replicate <- polr(data = immig.replicate, dv ~ edlevel2 + group + group:edlevel2)

replicate.to.plot <- olMFX(reg.replicate)

mfx.df.rep <- data.frame(mfx= replicate.to.plot$out$`ME.Strongly disagree`[, "effect"], 
                         se = replicate.to.plot$out$`ME.Strongly disagree`[, "error"])
mfx.df.rep$u.ci <- mfx.df.rep$mfx + (1.96 * mfx.df.rep$se)
mfx.df.rep$l.ci <- mfx.df.rep$mfx - (1.96 * mfx.df.rep$se)

mfx.df.rep$var <- makeNice(rownames(mfx.df.rep))
mfx.df.rep <- rbind(mfx.df.rep, NA, NA)

mfx.df.rep$var[is.na(mfx.df.rep$var)] <- c("High education", "Treatment: both high")
mfx.df.rep$mfx[is.na(mfx.df.rep$mfx)] <- 0
mfx.df.rep$var <- makeNiceFactor(mfx.df.rep$var)

ggplot(data = mfx.df.rep, aes(x = var, y = mfx)) + geom_point() + 
  geom_errorbar(aes(ymax = u.ci, ymin = l.ci), width = 0.25) +
  geom_hline(aes(yintercept=0)) +  theme_bw() + 
  xlab("") + 
  ylab("Marginal effect on 'strongly disagreeing' with more immigrants") +
  coord_flip()

## ----manipCheckCompete, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE, fig.cap="Effects of education and immigrant skill level on probability of respondent being `not at all worried' about their own job prospects (excludes `own job' treatment group)"----


# prop.table(table(immig.replicate$immigManipCheck, immig.replicate$edlevel2), 1)
# immig.replicate$group

large.labmarket.rep <- polr(data = immig.replicate, 
                            formula = factor(immigManipCheck) ~ group + edlevel2 + group:edlevel2)
manip.compet.mfx <- olMFX(large.labmarket.rep)
manip.compet.mfx <- data.frame(manip.compet.mfx$out$ME.1)

manip.compet.mfx$u.ci <- manip.compet.mfx[, "effect"] + (manip.compet.mfx[, "error"] * 1.96)
manip.compet.mfx$l.ci <- manip.compet.mfx[, "effect"] - (manip.compet.mfx[, "error"] * 1.96)
manip.compet.mfx$var <- makeNice(rownames(manip.compet.mfx))
manip.compet.mfx <- rbind(manip.compet.mfx, NA, NA)
manip.compet.mfx$var[is.na(manip.compet.mfx$var)] <- 
  c("High education", "Treatment: both high")
manip.compet.mfx$var <- makeNiceFactor(manip.compet.mfx$var)
manip.compet.mfx$effect[is.na(manip.compet.mfx$effect)] <- 0

ggplot(data = manip.compet.mfx, aes(x = var, y = effect)) + geom_point() + 
  geom_errorbar(aes(ymax = u.ci, ymin = l.ci), width = 0.25) +
  geom_hline(aes(yintercept=0)) +  theme_bw() + 
  xlab("") + ylab("Marginal effect on probability of respondent feeling 'not at all worried' about their job prospects") + coord_flip()

## ----replicateRegTable, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE----
texreg(reg.replicate, caption = "Ordered logistic regression predicting acceptance of immigrants in control groups", label = "table:mainreplication",
       custom.coef.names = makeNice(names(reg.replicate$coefficients)))


## ----manipCheckControlGroup, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE----

texreg(large.labmarket.rep, 
       label = "table:manipcheckreplicate",
       caption = "Ordered logistic regression predictors of economic threat",
       custom.coef.names = makeNice(names(large.labmarket.rep$coefficients)))

## ----balanceCheck, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE----
immig$immigSelf[immig$immigSelf==111] <- 0
immig$immigSelf[immig$immigSelf==99] <-NA

immig.att.balance <- polr(data = immig, factor(immigSelf)~group)

ed.balance <- chisq.test(x = immig$edlevel2, y = immig$group)


texreg(immig.att.balance, 
       label = "table:balanceCheckImmig",
       caption = "Ordered logistic regression predictors of immigration attitudes using random groups")

## ----hhFig1, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message = FALSE, fig.height=5, fig.cap="Bivariate relationship between education and anti-immigrant responses in the experiment (excludes own job experimental group)"----

immig.ed <- prop.table(table(immig.replicate$dv, immig.replicate$edlevel1), 2)
immig.ed<- melt(immig.ed)
colnames(immig.ed) <- c("More immigrants", "Education", "Fraction")

immig.ed$`More immigrants` <- factor(immig.ed$`More immigrants`, levels = levels(immig.replicate$dv))
immig.ed$Education <- factor(immig.ed$Education, levels = levels(immig.replicate$edlevel1))

ggplot(data = immig.ed, aes(x = Education, y = Fraction, group = `More immigrants`,
                            fill = `More immigrants`)) +  
  geom_bar(position = "dodge", colour = "black", 
           stat = "identity") + 
  theme_bw() + 
  scale_fill_grey(start = 0, end = 1, na.value = "blue") + 
  ylab("Fraction") + xlab("Education")


## ----hhFig2, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message = FALSE, fig.height=5, fig.cap = "Control group support for more immigrants by skill level of immigrants  (excludes own job experimental group)"----
### H&H figure 2 replication

counts.replicate <- table(immig.replicate$dv, immig.replicate$group)
counts.replicate <- prop.table(counts.replicate, 2)
counts.replicate <- melt(counts.replicate)

colnames(counts.replicate) <- c("More immigrants", "Treatment", "Count")
counts.replicate$Treatment <- factor(counts.replicate$Treatment, levels = c("Low", "Mixed", "High"))

counts.replicate$`More immigrants` <- factor(counts.replicate$`More immigrants`, levels = levels(immig.replicate$dv))

ggplot(counts.replicate, aes(x = Treatment, group = `More immigrants`, y = Count, fill = `More immigrants`)) +   geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  theme_bw() + 
  scale_fill_grey(start = 0, end = 1, na.value = "red")  + 
  ylab("Fraction") + xlab("Skill level (within control group)")

kable(manip.check.table, format = "latex", caption = "Level of concern about their own job for respondents in own job group and other groups")


## ----dk.check, cache=TRUE, eval=TRUE, dpi=100, echo=FALSE, results='asis', message=FALSE, warning=FALSE----

immig$dk.binary <- immig$immigExpDV==99
dk.pred <- glm(data = immig, dk.binary ~ group)
texreg(dk.pred, 
       label = "table:dk.pred",
       caption = "Logistic regression of reporting don't know to support for immigration across experimental groups", custom.coef.names = makeNice(names(dk.pred$coefficients)))

