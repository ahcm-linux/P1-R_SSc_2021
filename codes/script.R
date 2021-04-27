# PARAMETERS ----------------------------------------------------------------------------------------------------

max_perc_NA <- 0.05 # maximum percentage of missing values
max_perc_0 <- 0.8   # maximum percentage of zeros
min_perc_0 <- 0.2   # minimum percentage of zeros
AABS <- c("ENG", "VEGFB", "PDGFB", "MET", "F2RL1", "C5AR1", "CASR", "FLT1", "PIGF", "FGF1", "C3AR1", "F2R", 
          "CHRM4", "HGF",  "ADRB1", "CHRM1", "ADRB2", "CHRM2", "CHRM5", "EDNRB", "KDR", "VEGFA", "AGTR1", "EGFR",
          "EDNRA", "EGF", "CHRM3", "CXCR4", "CXCR3", "PDGFA", "YBX1") # selected set of aabs for analysis
col_sscf_1 <- "tomato2"     # color for represent sscf = 1 (diffuse SSc)
col_sscf_2 <- "dodgerblue2" # color for represent sscf = 2 (limited SSc)

# WARN: sex (1 = male, 2 = female); symptom (1 = present, 0 = absent); SSc (1 = diffuse, 2 = limited)


# EXTERNAL SOURCES ----------------------------------------------------------------------------------------------

source("codes/packages.R")
source("codes/functions.R")


# DATA ----------------------------------------------------------------------------------------------------------

# covariates
covariates_set1 <- read.csv("data/covariates_set1.csv", header = TRUE)
covariates_set2 <- read.csv("data/covariates_set2.csv", header = TRUE)

# symptoms as factors
symptoms <- read.csv("data/symptoms.csv", header = TRUE, colClasses = "factor")

# aabs
aabs_set1 <- read.csv("data/aabs_set1.csv", header = TRUE)
aabs_set2 <- read.csv("data/aabs_set2.csv", header = TRUE)
aabs_set3 <- read.csv("data/aabs_set3.csv", header = TRUE)
aabs_set4 <- read.csv("data/aabs_set4.csv", header = TRUE)


# DELETING DATA -------------------------------------------------------------------------------------------------

# identify patients with undefined SSCF
del <- which(covariates_set2$SSCF == 4)

# delete patients with undefined type of SSCF
symptoms <- symptoms[-del, ]
covariates_set1 <- covariates_set1[-del, ]
covariates_set2 <- covariates_set2[-del, ]

# delete patients with undefined type of SSCF and keep only aabs in AABS
aabs_set1 <- aabs_set1[-del, colnames(aabs_set1) %in% AABS]
aabs_set2 <- aabs_set2[-del, colnames(aabs_set2) %in% AABS]
aabs_set3 <- aabs_set3[-del, colnames(aabs_set3) %in% AABS]
aabs_set4 <- aabs_set4[-del, colnames(aabs_set4) %in% AABS]


# AABS DATA -----------------------------------------------------------------------------------------------------

# list of aabs data
aabs_list <- list(aabs_set1, aabs_set2, aabs_set3, aabs_set4)

# types of aabs
aabs_group <- unlist(lapply(1:length(aabs_list), function(.x) rep(.x, ncol(aabs_list[[.x]]))))

# aabs data
aabs <- do.call(cbind, aabs_list)


# SYMPTOMS SELECTION --------------------------------------------------------------------------------------------

# check variablity level in factor variables
prop_0 <- apply(symptoms, 2, function(.x) prop.table(table(.x))[1])
select_symptoms <- names(prop_0[prop_0 > min_perc_0 & prop_0 < max_perc_0])

# final set of symptoms
symptoms <- symptoms[, select_symptoms]


# COVARIATES SELECTION ------------------------------------------------------------------------------------------

# selected variables from covariates_set1
select_covar_set1 <- c("Age", "sex")

# selected variables from covariates_set2
n <- nrow(covariates_set2)
prop_na <- apply(covariates_set2, 2, function(.x) sum(is.na(.x)) / n)
select_covar_set2 <- colnames(covariates_set2)[prop_na < max_perc_NA]

# selected covariates after checking for level of missing data
select_covariates <- cbind(covariates_set1[, select_covar_set1], covariates_set2[, select_covar_set2])

# classify some of the covariates as factors
to_factor <- c(2, 5, 11:ncol(select_covariates))
for (.x in to_factor) select_covariates[, .x] <- factor(select_covariates[, .x])

# check variablity level in factor variables
prop_1 <- apply(select_covariates[, to_factor], 2, function(.x) prop.table(table(.x))[1])
select_vars_1 <- names(prop_1[prop_1 > min_perc_0 & prop_1 < max_perc_0])
if (!("sex" %in% select_vars_1)) select_vars_1 <- c(select_vars_1, "sex")

# final set of covariates
selec_vars_2 <- setdiff(colnames(select_covariates)[-to_factor], select_vars_1)
covariates <- select_covariates[, c(select_vars_1, selec_vars_2)]

# missing data imputation
covariates <- mice::complete(mice::mice(covariates, m = 1, method = "cart"))


# DATA TRANSFORMATION -------------------------------------------------------------------------------------------

# generalized log transformation
aabs_glog <- sapply(aabs, DescTools::LogSt)

# scaled covariates
covariates_scale <- covariates
covariates_scale[sapply(covariates, is.numeric)] <- scale(covariates[sapply(covariates, is.numeric)])


# DESCRIPTIVE STATISTICS ----------------------------------------------------------------------------------------

# number of variables and observations
n_symp <- ncol(symptoms)
n_aabs <- ncol(aabs)
n_vars <- n_symp + n_aabs
n_obs <- nrow(covariates)
n_obs_sscf_1 <- nrow(subset(covariates, SSCF == "1"))
n_obs_sscf_2 <- nrow(subset(covariates, SSCF == "2"))

# table
t_desc <- data.frame(n_symp = n_symp, n_aabs = n_aabs, n_obs_diffuseSSc = n_obs_sscf_1, n_obs_limitedSSc = n_obs_sscf_2)

## untransformed variables -------------------------------------------------------------------------------------

# symptoms: summary
summ_symp <- summary(symptoms)

# covariates: summary and normality checking for continuous variables
summ_covars <- summary(covariates)
mvn_covars <- mvn(covariates[sapply(covariates, is.numeric)])

# aabs: summary and normality checking
summ_aabs <- summary(aabs)
mvn_aabs <- mvn(aabs)

# symptoms as matrix
Msymptoms <- as.matrix(sapply(symptoms, function(.x) as.numeric(as.character(.x))))

# plot description of symptoms for each ssc form
df_prop_symp <- reshape2::melt(data.frame(Msymptoms, sscf = covariates$SSCF), id.vars = "sscf")
prop_symp <- plyr::ddply(df_prop_symp, ~ sscf + variable, dplyr::summarise, total_1 = sum(value))
prop_symp$perc <- ifelse(prop_symp$sscf == "1", prop_symp$total_1 / n_obs_sscf_1, prop_symp$total_1 / n_obs_sscf_2)

diff_prop_symp <- abs(plyr::arrange(subset(prop_symp, sscf == "1"), variable)$perc - plyr::arrange(subset(prop_symp, sscf == "2"), variable)$perc)
order_prop_symp <- arrange(subset(prop_symp, sscf == "1"), variable)$variable[rev(order(diff_prop_symp))]

prop_symp$variable <- factor(prop_symp$variable, levels = as.character(order_prop_symp), ordered = TRUE)

p_prop_symp <- ggplot2::ggplot(data = prop_symp, ggplot2::aes(x = variable, y = perc, group = sscf, fill = sscf)) +
  ggplot2::geom_bar(position = ggplot2::position_dodge(), stat = "identity", alpha = 0.6, colour = "black", width = 0.4) +
  ggplot2::scale_fill_manual("SSc: ", labels = c("1" = "Diffuse", "2" = "Limited"), values = c("1" = col_sscf_1, "2" = col_sscf_2)) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90)) +
  ggplot2::labs(x = "Symptoms", y = "Sample prevalence")

# plot description of aabs for each ssc form
df_aabs_sscf <- data.frame(aabs, sscf = covariates$SSCF)

median_aabs_sscf <- df_aabs_sscf %>% dplyr::group_by(sscf) %>% dplyr::summarise(dplyr::across(dplyr::everything(), median))
q75_aabs_sscf <- df_aabs_sscf %>% dplyr::group_by(sscf) %>% dplyr::summarise(dplyr::across(dplyr::everything(), ~quantile(.x, probs = 0.75)))
q25_aabs_sscf <- df_aabs_sscf %>% dplyr::group_by(sscf) %>% dplyr::summarise(dplyr::across(dplyr::everything(), ~quantile(.x, probs = 0.25)))

melt_median <- melt(median_aabs_sscf, id.vars = "sscf")
melt_q75 <- melt(q75_aabs_sscf, id.vars = "sscf")
melt_q25 <- melt(q25_aabs_sscf, id.vars = "sscf")

melt_aabs_sscf <- Reduce(function(...) merge(..., by = c("sscf", "variable")), list(melt_median, melt_q75, melt_q25))
colnames(melt_aabs_sscf)[3:5] <- c("q50", "q75", "q25")
melt_aabs_sscf <- plyr::arrange(melt_aabs_sscf, sscf, -q50, variable)
melt_aabs_sscf$variable <- factor(melt_aabs_sscf$variable, levels = subset(melt_aabs_sscf, sscf == "1")$variable, ordered = TRUE)

p_quantile_aabs_sscf <- ggplot2::ggplot(data = melt_aabs_sscf, ggplot2::aes(x = variable, y = q50, colour = sscf, group = sscf)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = q25, ymax = q75), alpha = 0.7, size = 0.8, width = 0.5, position = ggplot2::position_dodge(0.5)) +
  ggplot2::geom_point(alpha = 1, size = 1.3, position = ggplot2::position_dodge(0.5)) +
  ggplot2::scale_colour_manual("SSc: ", labels = c("1" = "Diffuse", "2" = "Limited"), values = c("1" = col_sscf_1, "2" = col_sscf_2)) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90)) +
  ggplot2::labs(x = "Autoantibodies", y = "Sample concentration")

## transformed variables ---------------------------------------------------------------------------------------

# scaled covariates: summary and normality checking for continuous variables
summ_covars_scale <- summary(covariates_scale)
mvn_covars_scale <- mvn(covariates_scale[sapply(covariates_scale, is.numeric)])

# aabs: summary and normality checking
summ_aabs_glog <- summary(aabs_glog)
mvn_aabs_glog <- mvn(aabs_glog)

# plot tansformed and untransformed aabs
df_aabs <- rbind(data.frame(aabs, trans = "0"), data.frame(aabs_glog, trans = "1"))
melt_df_aabs <- reshape2::melt(df_aabs, id.vars = "trans")

p_aabs_dist <- ggplot2::ggplot(data = melt_df_aabs, ggplot2::aes(x = variable, y = value)) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  ggplot2::facet_wrap(~ trans, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(c("0" = "Raw autoantibodies concentration", "1" = "Generalized logarithm of autoantibodies concentration"))) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "", y = "Value")


# NONPARAMETRIC TESTS -------------------------------------------------------------------------------------------

# tests for difference in proportions of symptoms between sscf
prop_test_symp <- do.call(rbind, lapply(colnames(symptoms), function(.x) {
  symp <- subset(prop_symp, variable == .x)
  total <- symp[, "total_1"]
  DescTools::BinomDiffCI(total[1], n_obs_sscf_1, total[2], n_obs_sscf_2, method = "ac", conf.level = .95)
}))
rownames(prop_test_symp) <- colnames(symptoms)

# plot difference of proportions tests
prop_test_symp <- data.frame(prop_test_symp, var = rownames(prop_test_symp))
df_prop_test_symp <- plyr::arrange(prop_test_symp, est)
df_prop_test_symp$symp <- factor(df_prop_test_symp$var, levels = df_prop_test_symp$var, ordered = TRUE)

p_prop_test_symp <- ggplot2::ggplot(data = df_prop_test_symp, ggplot2::aes(x = symp, ymin = lwr.ci, ymax = upr.ci)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7, size = 0.5) +
  ggplot2::geom_point(aes(y = est), shape = 20, size = 2) +
  ggplot2::geom_errorbar(alpha = 0.5, size = 0.5, width = 0.5) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Symptom", y = "CI for difference in prevalences between\nforms of SSc (diffuse SSc - limited SSc)")

# Mann-Whitney tests for differences in aabs between sscf
suppressWarnings({
  mann_whitney_test_aabs <- do.call(rbind, lapply(colnames(aabs), function(.x) {
    aab <- aabs[, .x]
    test <- stats::wilcox.test(aab ~ covariates$SSCF, correct = TRUE, paired = FALSE)
    data.frame(statistic = test$statistic, pvalue = test$p.value)
  }))
  rownames(mann_whitney_test_aabs) <- colnames(aabs)
})

# plot Mann-Whitney tests
mann_whitney_test_aabs$aabs <- factor(rownames(mann_whitney_test_aabs), levels = unique(melt_aabs_sscf$variable), ordered = TRUE)

p_MannWhitney_aabs <- ggplot2::ggplot(data = mann_whitney_test_aabs, ggplot2::aes(x = aabs, y = pvalue, colour = ifelse(pvalue < 0.1, "1", "0"))) +
  ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed", alpha = 0.7, size = 0.5) +
  ggplot2::geom_point(size = 2, shape = 20) +
  ggplot2::geom_segment(ggplot2::aes(xend = aabs, y = 0, yend = pvalue), size = 0.5) +
  ggplot2::scale_colour_manual("p-value: ", values = c("1" = "black", "0" = "grey"), labels = c("1" = "< 0.10", "0" = "\u2265 0.10")) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Autoantibody", y = "p-value of Mann-Whitney test for difference\nin aab concentration between forms of SSC")


# MIXED CANONICAL CORRELATION ANALYSIS --------------------------------------------------------------------------

# canonical correlation analysis
cca <- mixedCCA::mixedCCA(X1 = Msymptoms, X2 = aabs_glog, type1 = "binary", type2 = "continuous", BICtype = 1, trace = TRUE)

# canonical weights
w11 <- cca$w1
w12 <- cca$w2
cor_w11_w12 <- cca$cancor

# canonical scores
cca_scores <- data.frame(score_w11 = Msymptoms %*% w11, score_w12 = aabs_glog %*% w12)

# kendal correlation matrix
K <- cca$KendallR
K12 <- K[colnames(symptoms), colnames(aabs)]

# plot of the latent kendal correlation between symptoms and aabs
p_heatmap_kendal <- pheatmap::pheatmap(K12, clustering_method = "ward.D2", border_color = grey(0.4), fontsize = 7,
                                color = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"), silent = TRUE,
                                cutree_rows = 2, cutree_cols = 2, legend_breaks = seq(-1, 1, 0.1))

# plot canonical weights
df_cca <- data.frame(weight = rbind(w11, w12), cv = rep(c("0", "1"), c(nrow(w11), nrow(w12))),
                     variable = c(colnames(symptoms), colnames(aabs)))
df_cca <- plyr::arrange(df_cca, -abs(weight))
df_cca$variable <- factor(df_cca$variable, levels = df_cca$variable, ordered = TRUE)

p_canonical_weights <- ggplot2::ggplot(data = df_cca, ggplot2::aes(xend = variable, yend = weight)) +
  ggplot2::geom_segment(ggplot2::aes(x = variable, y = 0), alpha = 0.7, size = 0.5) +
  ggplot2::geom_point(ggplot2::aes(x = variable, y = weight), alpha = 0.8, size = 1, show.legend = FALSE) +
  ggplot2::facet_wrap(~ cv, scales = "free", ncol = 1, labeller = ggplot2::as_labeller(c("0" = "First Canonical Variate for Symptoms", "1" = "First Canonical Variate for Autoantibodies"))) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "", y = "Estimated canonical weight")

# table
t_cca <- data.frame(canonical_cor = round(cor_w11_w12, 3))


# HIGH DIMENTISIONAL LINEAR DISCRIMINANT ANALYSIS ---------------------------------------------------------------

# lda for high dimensional data
suppressWarnings(lda <- HiDimDA::Slda(data = aabs_glog, grouping = covariates$SSCF, VSelfunct = "none"))

# error of classification
pred_lda <- predict(lda, aabs_glog, grpcodes = levels(covariates$SSCF))
error_lda <- table(pred_lda$class, covariates$SSCF)

# plot lda scores
df_lda_scores <- data.frame(score = as.vector(pred_lda$Z), class = covariates$SSCF)

p_lda_scores <- ggplot2::ggplot(data = df_lda_scores, ggplot2::aes(x = score, colour = class, fill = class)) +
  ggplot2::geom_density(position = "identity", alpha = 0.1, size = 0.5) +
  ggplot2::scale_colour_manual("SSc: ", labels = c("1" = "Diffuse", "2" = "Limited"), values = c("1" = col_sscf_1, "2" = col_sscf_2)) +
  ggplot2::scale_fill_manual(values = c("1" = col_sscf_1, "2" = col_sscf_2), guide = FALSE) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::labs(x = "Patient score obtained through the linear discriminant function", y = "Density")

# plot lda weights
df_lda_weights <- data.frame(weight = as.vector(lda$scaling), variable = rownames(lda$scaling))
df_lda_weights <- arrange(df_lda_weights, -abs(weight))
df_lda_weights$variable <- factor(df_lda_weights$variable, levels = df_lda_weights$variable, ordered = TRUE)

ct_lda <- round(quantile(abs(as.vector(lda$scaling)), probs = 0.75), 1) # 75-th percentile of the absolute value of lda weights  

p_lda_weights <- ggplot2::ggplot(data = df_lda_weights, ggplot2::aes(xend = variable, yend = weight)) +
  ggplot2::geom_segment(ggplot2::aes(x = variable, y = 0), alpha = 0.7, size = 0.5) +
  ggplot2::geom_point(ggplot2::aes(x = variable, y = weight), alpha = 0.8, size = 1, show.legend = FALSE) +
  ggplot2::geom_hline(yintercept = c(-ct_lda, ct_lda), linetype = "dashed", size = 0.5, colour = "black", alpha = 0.6) +
  ggplot2::scale_y_continuous(breaks = seq(-10, 10, 1)) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.y = ggplot2::element_blank(),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Autoantibody", y = "Estimated linear discriminant weight")

# table
t_lda <- error_lda
colnames(t_lda) <- rownames(t_lda) <- c("Diffuse", "Limited")


# SYMPTOMS MODEL ------------------------------------------------------------------------------------------------

# select symptoms and aabs with nonzero canonical weights
symptoms_from_cca <- symptoms[, which(w11 != 0 | colnames(symptoms) %in% c("ILF"))]
aabs_from_cca <- aabs_glog[, which(w12 != 0)]

# data for regression
df_symp <- data.frame(ID = factor(1:n_obs), symptoms_from_cca)
df_covar <- data.frame(ID = df_symp$ID, aabs_from_cca, covariates[, c("Age", "sex", "SSCF")])
melt_symp <- reshape2::melt(df_symp, id.vars = "ID")
df_reg_symp <- merge(melt_symp, df_covar, by = "ID")
df_reg_symp$value <- factor(df_reg_symp$value)

# data for lasso regularization
vars_ri_symp <- colnames(df_reg_symp)[colnames(df_reg_symp) %in% colnames(aabs)]
df_ri_symp <- df_reg_symp[, vars_ri_symp]

# modeling probability of a symptom (using CXCR3 and SSCF as random coefficients)
set.seed(2021) # seed to replicate randomized residuals
fit_reg_symp <- gamlss::gamlss(value ~ log(Age) + sex + SSCF + CXCR3 + re(random =~ 1 + SSCF + CXCR3|variable) + re(random =~ 1|ID),
                               family = BI, data = df_reg_symp, trace = FALSE)

# model adequacy
box_test_symp <- Box.test(residuals(fit_reg_symp))
#plot(fit_reg_symp)
#rqres.plot(fit_reg_symp, plot = 'all', howmany = 100)
drop_test_symp <- drop1(fit_reg_symp)

# parameter estimates
invisible(utils::capture.output(summary_fit_reg_symp <- summary(fit_reg_symp)))

# smooth parameters and random effects
ranef_symp <- ranef(getSmo(fit_reg_symp, which = 1))
ranef_ID_symp <- ranef(getSmo(fit_reg_symp, which = 2))

(ranef_summary_symp <- summary(getSmo(fit_reg_symp, which = 1))) # print
(ranef_summary_ID_symp <- summary(getSmo(fit_reg_symp, which = 2))) # print

# model predictions
pred_data_symp_CXCR3 <- data.frame(
  Age = mean(df_reg_symp$Age),
  sex = "1",
  SSCF = rep(c("1", "2"), each = 100 * 2 * nrow(symptoms_from_cca)),
  CXCR3 = rep(seq(min(df_reg_symp$CXCR3), max(df_reg_symp$CXCR3), length = 100), 2 * nrow(symptoms_from_cca)),
  variable = rep(rep(colnames(symptoms_from_cca), each = 100), 2),
  ID = "1")
pred_symp_CXCR3 <- predict(fit_reg_symp, newdata = pred_data_symp_CXCR3, type = "response")

pred_data_symp_Age <- data.frame(
  Age = rep(seq(min(df_reg_symp$Age), max(df_reg_symp$Age), length = 100), 2 * nrow(symptoms_from_cca)),
  sex = "1",
  SSCF = rep(c("1", "2"), each = 100 * 2 * nrow(symptoms_from_cca)),
  CXCR3 = mean(df_reg_symp$CXCR3),
  variable = rep(rep(colnames(symptoms_from_cca), each = 100), 2),
  ID = "1")
pred_symp_Age <- predict(fit_reg_symp, newdata = pred_data_symp_Age, type = "response")

# plot predictions for CXCR3
df_pred_symp_CXCR3 <- data.frame(pred_data_symp_CXCR3, pred = pred_symp_CXCR3)
df_pred_symp_CXCR3$variable <- factor(df_pred_symp_CXCR3$variable, levels = sort(unique(df_pred_symp_CXCR3$variable)), ordered = TRUE)

p_pred_symp_CXCR3 <- ggplot2::ggplot(data= df_pred_symp_CXCR3, ggplot2::aes(x = CXCR3, y = pred, group = SSCF, colour = SSCF)) +
  ggplot2::geom_line(alpha = 0.7, size = 0.5) +
  ggplot2::geom_rug(sides = "tb", alpha = 0.7, size = 0.2, inherit.aes = FALSE, data = df_reg_symp, ggplot2::aes(x = CXCR3, y = ifelse(value == "0", min(pred_symp_CXCR3), max(pred_symp_CXCR3)), colour = SSCF)) +
  ggplot2::facet_wrap(~ variable) +
  ggplot2::scale_colour_manual("SSc: ", labels = c("1" = "Diffuse", "2" = "Limited"), values = c("1" = col_sscf_1, "2" = col_sscf_2)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Generalized logarithm of CXCR3 concentration", y = "Predicted probability of symptom occurrence")

# plot predictions for Age
df_pred_symp_Age <- data.frame(pred_data_symp_Age, pred = pred_symp_Age)
df_pred_symp_Age$variable <- factor(df_pred_symp_Age$variable, levels = sort(unique(df_pred_symp_Age$variable)), ordered = TRUE)

p_pred_symp_Age <- ggplot2::ggplot(data= df_pred_symp_Age, ggplot2::aes(x = Age, y = pred, group = SSCF, colour = SSCF)) +
  ggplot2::geom_line(alpha = 0.7, size = 0.5) +
  ggplot2::geom_rug(sides = "tb", alpha = 0.7, size = 0.2, inherit.aes = FALSE, data = df_reg_symp, ggplot2::aes(x = Age, y = ifelse(SSCF == "1", min(pred_symp_Age), max(pred_symp_Age)), colour = SSCF)) +
  ggplot2::facet_wrap(~ variable) +
  ggplot2::scale_colour_manual("SSC form: ", labels = c("1" = "Diffuse", "2" = "Limited"), values = c("1" = col_sscf_1, "2" = col_sscf_2)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Age", y = "Predicted probability of symptom occurrence")

# residuals plots
df_res_symp <- data.frame(res = residuals(fit_reg_symp), fit = fitted(fit_reg_symp), index = 1:nrow(df_reg_symp))
p_res_fit_symp <- ggplot2::ggplot(df_res_symp, ggplot2::aes(x = fit, y = res)) +
  ggplot2::geom_hline(yintercept = c(-2, 2), linetype = 'dashed', size = 0.5) +
  ggplot2::geom_point(shape = 20, size = 2, alpha = 0.5) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Quantile residual", x = "Fitted value")
p_res_index_symp <- ggplot2::ggplot(df_res_symp, ggplot2::aes(x = index, y = res)) +
  ggplot2::geom_hline(yintercept = c(-2, 2), linetype = 'dashed', size = 0.5) +
  ggplot2::geom_point(shape = 20, size = 2, alpha = 0.5) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Quantile residual", x = "Index of observation")
p_res_density_symp <- ggplot2::ggplot(df_res_symp, ggplot2::aes(x = res, fill = "1")) +
  ggplot2::geom_density(position = "identity", alpha = 0.1, size = 0.5, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = grey(0.4)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Density", x = "Quantile residual")
p_res_qq_symp <- ggplot2::ggplot(df_res_symp, ggplot2::aes(sample = res)) +
  ggplot2::stat_qq(size = 1, alpha = 0.5) + ggplot2::stat_qq_line() +
  ggplot2::scale_fill_manual(values = grey(0.4)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Sample quantile", x = "Theoretical quantile")
p_diagnostics_symp <- gridExtra::arrangeGrob(p_res_fit_symp, p_res_index_symp, p_res_density_symp, p_res_qq_symp, nrow = 2)

# tables
t_reg_symp <- as.data.frame(summary_fit_reg_symp)
t_reg_symp$`P(>LRT)` <- NA
last_rows_symp <- grep('^(re)', rownames(drop_test_symp))
t_last_rows_symp <- data.frame(NA, NA, NA, NA, drop_test_symp$`Pr(Chi)`[last_rows_symp])
colnames(t_last_rows_symp) <- colnames(t_reg_symp)
rownames(t_last_rows_symp) <- rownames(drop_test_symp)[last_rows_symp]
t_reg_symp <- round(rbind(t_reg_symp, t_last_rows_symp), 3)


# SSC FORM MODEL ------------------------------------------------------------------------------------------------

# select aabs with nonzero canonical weights
aabs_from_lda <- aabs_glog[, which(abs(lda$scaling) > ct_lda | colnames(aabs) %in% c("ADRB1"))]

# all aabs
abbs_cca_lda <- unique(c(colnames(aabs_from_cca), colnames(aabs_from_lda)))
aabs_all <- cbind(aabs_from_cca, aabs_from_lda)[, abbs_cca_lda]

# data for regression
df_reg_sscf <- data.frame(aabs_all, covariates[, c("Age", "sex", "SSCF")])

# data for lasso regularization
vars_ri_sscf <- colnames(df_reg_sscf)[colnames(df_reg_sscf) %in% colnames(aabs)]
df_ri_sscf <- df_reg_sscf[, vars_ri_sscf]

# reference level of SSCF
df_reg_sscf$SSCF <- stats::relevel(df_reg_sscf$SSCF, ref = "1")

# modeling probability of a ssc form
set.seed(2021) # seed to replicate randomized residuals
fit_reg_sscf <- gamlss::gamlss(SSCF ~ log(Age) + sex + ri(df_ri_sscf, x.vars = vars_ri_sscf, Lp = 1),
                               family = BI, data = df_reg_sscf, trace = FALSE)

# model adequacy
box_test_sscf <- Box.test(residuals(fit_reg_sscf))
#plot(fit_reg_sscf)
#rqres.plot(fit_reg_sscf, plot = 'all', howmany = 100)
drop_test_sscf <- drop1(fit_reg_sscf)

# parameter estimates
invisible(utils::capture.output(summary_fit_reg_sscf <- summary(fit_reg_sscf)))

# smooth parameters and random effects
coefs_aabs <- round(coef(getSmo(fit_reg_sscf, which = 1)), 3)
lambda_lasso <- getSmo(fit_reg_sscf, which = 1)$lambda

# model predictions
x_aabs <- do.call(cbind, lapply(df_reg_sscf[, c("ADRB2", "YBX1", "EDNRB", "FLT1")], function(.x) seq(min(.x), max(.x), length = 100)))
colnames(x_aabs) <- c("ADRB2", "YBX1", "EDNRB", "FLT1")

pred_data_sscf_aabs <- data.frame(
  intercept = 1,
  logAge = mean(log(df_reg_sscf$Age)),
  sex = 0,
  ADRB2 = c(rep(mean(df_reg_sscf$ADRB2), 0 * 100), x_aabs[, "ADRB2"], rep(mean(df_reg_sscf$ADRB2), 3 * 100)),
  AGTR1 = mean(df_reg_sscf$AGTR1),
  CHRM4 = mean(df_reg_sscf$CHRM4),
  CXCR3 = mean(df_reg_sscf$CXCR3),
  FGF1 = mean(df_reg_sscf$FGF1),
  YBX1 = c(rep(mean(df_reg_sscf$YBX1), 1 * 100), x_aabs[, "YBX1"], rep(mean(df_reg_sscf$YBX1), 2 * 100)),
  ADRB1 = mean(df_reg_sscf$ADRB1),
  CXCR4 = mean(df_reg_sscf$CXCR4),
  EDNRA = mean(df_reg_sscf$EDNRA),
  EDNRB = c(rep(mean(df_reg_sscf$EDNRB), 2 * 100), x_aabs[, "EDNRB"], rep(mean(df_reg_sscf$EDNRB), 1 * 100)),
  FLT1 =  c(rep(mean(df_reg_sscf$FLT1),  3 * 100), x_aabs[, "FLT1"],  rep(mean(df_reg_sscf$FLT1),  0 * 100)),
  PIGF = mean(df_reg_sscf$PIGF))
#pred_sscf_aabs <- DescTools::LogitInv(as.matrix(pred_data_sscf_aabs) %*% matrix(c(coef(fit_reg_sscf)[1:3], coefs_aabs), ncol = 1))
pred_sscf_aabs <- as.matrix(pred_data_sscf_aabs) %*% matrix(c(coef(fit_reg_sscf)[1:3], coefs_aabs), ncol = 1)

pred_data_sscf_Age <- data.frame(
  intercept = 1,
  logAge = seq(min(log(df_reg_sscf$Age)), max(log(df_reg_sscf$Age)), length = 100),
  sex = 0,
  ADRB2 = mean(df_reg_sscf$ADRB2),
  AGTR1 = mean(df_reg_sscf$AGTR1),
  CHRM4 = mean(df_reg_sscf$CHRM4),
  CXCR3 = mean(df_reg_sscf$CXCR3),
  FGF1 = mean(df_reg_sscf$FGF1),
  YBX1 = mean(df_reg_sscf$YBX1),
  ADRB1 = mean(df_reg_sscf$ADRB1),
  CXCR4 = mean(df_reg_sscf$CXCR4),
  EDNRA = mean(df_reg_sscf$EDNRA),
  EDNRB = mean(df_reg_sscf$EDNRB),
  FLT1 =  mean(df_reg_sscf$FLT1),
  PIGF = mean(df_reg_sscf$PIGF))
#pred_sscf_Age <- DescTools::LogitInv(as.matrix(pred_data_sscf_Age) %*% matrix(c(coef(fit_reg_sscf)[1:3], coefs_aabs), ncol = 1))
pred_sscf_Age <- as.matrix(pred_data_sscf_Age) %*% matrix(c(coef(fit_reg_sscf)[1:3], coefs_aabs), ncol = 1)

# plot prediction for aabs
melt_df_pred_sscf_aabs <- cbind(reshape2::melt(x_aabs), pred = pred_sscf_aabs)
melt_df_pred_sscf_aabs$Var2 <- factor(melt_df_pred_sscf_aabs$Var2, levels = c("EDNRB", "YBX1", "FLT1", "ADRB2"), ordered = TRUE)

df_obs_sscf_aabs <- df_reg_sscf[, c("ADRB2", "YBX1", "EDNRB", "FLT1", "SSCF")]
melt_df_obs_sscf_aabs <- reshape2::melt(df_obs_sscf_aabs, id.vars = "SSCF")
colnames(melt_df_obs_sscf_aabs)[2] <- "Var2"
melt_df_obs_sscf_aabs$Var2 <- factor(melt_df_obs_sscf_aabs$Var2, levels = c("EDNRB", "YBX1", "FLT1", "ADRB2"), ordered = TRUE)

p_pred_sscf_aabs <- ggplot2::ggplot(data= melt_df_pred_sscf_aabs, ggplot2::aes(x = value, y = pred, group = Var2)) +
  ggplot2::geom_line(alpha = 0.5, size = 0.5, linetype = "solid", colour = "black") +
  ggplot2::geom_rug(alpha = 0.4, sides = "tb", size = 0.2, inherit.aes = FALSE, data = melt_df_obs_sscf_aabs, ggplot2::aes(x = value, y = ifelse(SSCF == "1", min(pred_sscf_aabs), max(pred_sscf_aabs)), group = Var2)) +
  ggplot2::facet_wrap(~ Var2, ncol = 2, scales = "free_x") +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Generalized logarithm of autoantibody concentration", y = "Log odds of diffuse SSc against limited SSc")

# plot predictions for Age
df_pred_sscf_Age <- data.frame(logAge = pred_data_sscf_Age[, "logAge"], pred = pred_sscf_Age)

p_pred_sscf_Age <- ggplot2::ggplot(data= df_pred_sscf_Age, ggplot2::aes(x = logAge, y = pred)) +
  ggplot2::geom_line(alpha = 0.7, size = 0.5) +
  ggplot2::geom_rug(sides = "tb", alpha = 0.7, size = 0.2, inherit.aes = FALSE, data = df_reg_sscf, ggplot2::aes(x = log(Age), y = ifelse(SSCF == "1", min(pred_sscf_Age), max(pred_sscf_Age)))) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Age", y = "Log odds of diffuse SSc against limited SSc")

# residuals plots
df_res_sscf <- data.frame(res = residuals(fit_reg_sscf), fit = fitted(fit_reg_sscf), index = 1:nrow(df_reg_sscf))
p_res_fit_sscf <- ggplot2::ggplot(df_res_sscf, ggplot2::aes(x = fit, y = res)) +
  ggplot2::geom_hline(yintercept = c(-2, 2), linetype = 'dashed', size = 0.5) +
  ggplot2::geom_point(shape = 20, size = 2, alpha = 0.5) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Quantile residual", x = "Fitted value")
p_res_index_sscf <- ggplot2::ggplot(df_res_sscf, ggplot2::aes(x = index, y = res)) +
  ggplot2::geom_hline(yintercept = c(-2, 2), linetype = 'dashed', size = 0.5) +
  ggplot2::geom_point(shape = 20, size = 2, alpha = 0.5) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Quantile residual", x = "Index of observation")
p_res_density_sscf <- ggplot2::ggplot(df_res_sscf, ggplot2::aes(x = res, fill = "1")) +
  ggplot2::geom_density(position = "identity", alpha = 0.1, size = 0.5, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = grey(0.4)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Density", x = "Quantile residual")
p_res_qq_sscf <- ggplot2::ggplot(df_res_sscf, ggplot2::aes(sample = res)) +
  ggplot2::stat_qq(size = 1, alpha = 0.5) + ggplot2::stat_qq_line() +
  ggplot2::scale_fill_manual(values = grey(0.4)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::labs(y = "Sample quantile", x = "Theoretical quantile")
p_diagnostics_sscf <- gridExtra::arrangeGrob(p_res_fit_sscf, p_res_index_sscf, p_res_density_sscf, p_res_qq_sscf, nrow = 2)

# tables
t_reg_sscf <- as.data.frame(summary_fit_reg_sscf)
t_reg_sscf$`P(>LRT)` <- NA
last_rows_sscf <- data.frame(coefs_aabs, NA, NA, NA, drop_test_sscf$`Pr(Chi)`[grep('^(ri)', rownames(drop_test_sscf))])
colnames(last_rows_sscf) <- colnames(t_reg_sscf)
rownames(last_rows_sscf) <- names(coefs_aabs)
t_reg_sscf <- round(rbind(t_reg_sscf, last_rows_sscf), 3)


# FACTOR ANALYSIS -----------------------------------------------------------------------------------------------

# data for fa
aabs_fa <- aabs_all

# kaiser criterion for the number of latent factors
k_sscf_1 <- sum(eigen(cor(subset(aabs_fa, covariates$SSCF == "1")))$values > 1)
k_sscf_2 <- sum(eigen(cor(subset(aabs_fa, covariates$SSCF == "2")))$values > 1)
k <- max(k_sscf_1, k_sscf_2)

# exploratory fa for sscf = 1
fa_sscf_1 <- suppressWarnings(fa(subset(aabs_fa, covariates$SSCF == "1"), nfactors = k, rotate = "varimax"))

# exploratory fa for sscf = 2
fa_sscf_2 <- suppressWarnings(fa(subset(aabs_fa, covariates$SSCF == "2"), nfactors = k, rotate = "varimax"))

# models for confirmatory fa
cfa_model <- 'lf1 =~ 1*ADRB2+AGTR1+CHRM4+CXCR3+FGF1+YBX1+ADRB1+CXCR4+EDNRA+EDNRB+FLT1+PIGF
              lf2 =~ 0*ADRB2+1*AGTR1+CHRM4+CXCR3+FGF1+YBX1+ADRB1+CXCR4+EDNRA+EDNRB+FLT1+PIGF
              lf3 =~ 0*ADRB2+0*AGTR1+1*CHRM4+CXCR3+FGF1+YBX1+ADRB1+CXCR4+EDNRA+EDNRB+FLT1+PIGF
              lf1 ~~ 1*lf1
              lf2 ~~ 1*lf2
              lf3 ~~ 1*lf3
              lf1 ~~ 0*lf2
              lf1 ~~ 0*lf3
              lf2 ~~ 0*lf3'

# data for confirmatory fa
df_cfa <- data.frame(aabs_fa, SSCF = covariates$SSCF)

# confirmatory fa
cfa_model_0 <- lavaan::cfa(cfa_model, data = df_cfa, std.ov = TRUE, meanstructure = FALSE, group = "SSCF", estimator = "MLM")
cfa_model_l <- lavaan::cfa(cfa_model, data = df_cfa, std.ov = TRUE, meanstructure = FALSE, group = "SSCF", group.equal = "loadings", estimator = "MLM")
cfa_model_lr <- lavaan::cfa(cfa_model, data = df_cfa, std.ov = TRUE, meanstructure = FALSE, group = "SSCF", group.equal = c("loadings", "residuals"), estimator = "MLM")

# model comparison
fit_anova <- lavaan::lavTestLRT(cfa_model_0, cfa_model_l, cfa_model_lr)
fit_aic <- AIC(cfa_model_0, cfa_model_l, cfa_model_lr)
fit_baic <- BIC(cfa_model_0, cfa_model_l, cfa_model_lr)

# chosen cfa model
fit_cfa <- cfa_model_lr

summary(fit_cfa) # print

est <- lavaan::lavInspect(fit_cfa, what = "est")
est_std <- lavaan::lavInspect(fit_cfa, what = "std")

est_list <- lavaan::parameterEstimates(fit_cfa)
est_std_list <- lavaan::standardizedsolution(fit_cfa)

fit_measures <- lavaan::fitMeasures(fit_cfa)
df_fit_measures <- data.frame(value = round(fit_measures, 3))

rotated_loading_matrix <- fa.sort(loadings(GPArotation::Varimax(est_std[[1]]$lambda)))
commonality <- apply(rotated_loading_matrix, 1, function(.x) sum(.x^2))
uniqueness <- diag(est_std[[1]]$theta)[names(commonality)]

# plot cfa model
cfa_l <- rotated_loading_matrix
df_cfa_results <- data.frame(Loading = c(cfa_l), Var = rep(rownames(cfa_l), 3), Factor = rep(c("1", "2", "3"), each = nrow(cfa_l)))
df_cfa_results$Var <- factor(df_cfa_results$Var, levels = rownames(cfa_l), ordered = TRUE)

p_cfa <- ggplot2::ggplot(df_cfa_results, ggplot2::aes(x = Var, y = abs(Loading), fill = Loading)) + 
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_hline(yintercept = 0.75, linetype = "dashed", colour = "black", size = 0.5) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ Factor, nrow = 1, labeller = ggplot2::as_labeller(c("1" = "Factor 1", "2" = "Factor 2", "3" = "Factor 3"))) +
  ggplot2::scale_fill_gradient2(name = "Loading: ", high = "blue", mid = "white", low = "red", midpoint = 0, guide = FALSE) +
  ggplot2::ylab("Loading Strength") +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Autoantibody", y = "Factor loading estimate")


# Mahalanobis distance
ind_sscf_1 <- as.matrix(subset(df_cfa, SSCF == "1")[, -which(colnames(df_cfa) == "SSCF")])
ind_sscf_2 <- as.matrix(subset(df_cfa, SSCF == "2")[, -which(colnames(df_cfa) == "SSCF")])

MD_sscf_1 <- sapply(1:n_obs_sscf_1, function(.x) {
  Sigma <- solve((est[[1]]$lambda %*% est[[1]]$psi %*% t(est[[1]]$lambda)) + est[[1]]$theta)
  y <- matrix(ind_sscf_1[.x, ], ncol = 1)
  MD <- t(y) %*% Sigma %*% y
}) 

MD_sscf_2 <- sapply(1:n_obs_sscf_2, function(.x) {
  Sigma <- solve((est[[1]]$lambda %*% est[[1]]$psi %*% t(est[[1]]$lambda)) + est[[1]]$theta)
  y <- matrix(ind_sscf_2[.x, ], ncol = 1)
  MD <- t(y) %*% Sigma %*% y
})

MD <- data.frame(md = c(MD_sscf_1, MD_sscf_2), sscf = rep(c("1", "2"), c(n_obs_sscf_1, n_obs_sscf_2)), ID = c(1:n_obs_sscf_1, 1:n_obs_sscf_2))

p_MD <- ggplot2::ggplot(data = MD, ggplot2::aes(x = ID, xend = ID, y = 0, yend = md, group = sscf)) +
  ggplot2::geom_hline(yintercept = qchisq(0.95, df = ncol(aabs_fa)), linetype = "dashed", colour = "black", alpha = 0.7, size = 0.5) +
  ggplot2::geom_segment(alpha = 0.5, size = 0.5) +
  ggplot2::facet_grid(~ sscf, scale = "free_x", labeller = ggplot2::as_labeller(c("1" = "Diffuse SSC", "2" = "Limited SSC"))) +
  ggplot2::theme_light(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 legend.box.background = ggplot2::element_rect(colour = "black", fill = NA),
                 panel.grid.minor.x = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold", colour = "black"),
                 strip.background = ggplot2::element_rect(fill = grey(0.9), colour = "black")) +
  ggplot2::labs(x = "Index of patient", y = "Mahalanobis distance")

# bootstrap analysis for model adequacy
set.seed(2021) # seed for the bootstrap simulation

fit_cfa_chisq <- lavaan::fitMeasures(fit_cfa, "chisq")
T_boot_chisq <- lavaan::bootstrapLavaan(fit_cfa, R = 1000, type = "bollen.stine", FUN = fitMeasures, fit.measures = "chisq")
pvalue_boot_chisq <- length(which(T_boot_chisq > fit_cfa_chisq)) / length(T_boot_chisq)

# tables
t_cfa <- est_std_list
t_cfa_lrt_pvalue <- pvalue_boot_chisq


# SAVE PLOTS -----------------------------------------------------------------------------------------------------

# pdf
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA1_Sample_Prevalence-Symptoms", ".pdf"), plot = p_prop_symp, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA3_Sample_aabs-distribution", ".pdf"), plot = p_quantile_aabs_sscf, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigS1_Sample_aabs-transformation", ".pdf"), plot = p_aabs_dist, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA2_Inference_Agresti-Coff-CIs", ".pdf"), plot = p_prop_test_symp, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA4_Inference_Mann-Whitney-test", ".pdf"), plot = p_MannWhitney_aabs, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA5_Sample_CCA-weights", ".pdf"), plot = p_canonical_weights, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA6_Sample_LDA-scores", ".pdf"), plot = p_lda_scores, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA7_Sample_LDA-weights", ".pdf"), plot = p_lda_weights, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA8_Inference_MGLM-CXCR3-SSc", ".pdf"), plot = p_pred_symp_CXCR3, width = 2 * 105, height = 1.7 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigS2_Inference_MGLM-Age-SSc", ".pdf"), plot = p_pred_symp_Age, width = 2 * 105, height = 1.7 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA9_Inference_PLR-aabs", ".pdf"), plot = p_pred_sscf_aabs, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigS3_Inference_PLR-Age", ".pdf"), plot = p_pred_sscf_Age, width = 1 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigA10_Inference_CFA-loadings", ".pdf"), plot = p_cfa, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigS6_Inference_CFA-diagnostics", ".pdf"), plot = p_MD, width = 1.5 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigS4_Inference_MGLM-diagnostics", ".pdf"), plot = p_diagnostics_symp, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/pdf/", "FigS5_Inference_PLR-diagnostics", ".pdf"), plot = p_diagnostics_sscf, width = 2 * 105, height = 1.5 * 74.25, units = "mm")

# png
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA1_Sample_Prevalence-Symptoms", ".png"), plot = p_prop_symp, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA3_Sample_aabs-distribution", ".png"), plot = p_quantile_aabs_sscf, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigS1_Sample_aabs-transformation", ".png"), plot = p_aabs_dist, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA2_Inference_Agresti-Coff-CIs", ".png"), plot = p_prop_test_symp, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA4_Inference_Mann-Whitney-test", ".png"), plot = p_MannWhitney_aabs, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA5_Sample_CCA-weights", ".png"), plot = p_canonical_weights, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA6_Sample_LDA-scores", ".png"), plot = p_lda_scores, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA7_Sample_LDA-weights", ".png"), plot = p_lda_weights, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA8_Inference_MGLM-CXCR3-SSc", ".png"), plot = p_pred_symp_CXCR3, width = 2 * 105, height = 1.7 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigS2_Inference_MGLM-Age-SSc", ".png"), plot = p_pred_symp_Age, width = 2 * 105, height = 1.7 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA9_Inference_PLR-aabs", ".png"), plot = p_pred_sscf_aabs, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigS3_Inference_PLR-Age", ".png"), plot = p_pred_sscf_Age, width = 1 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigA10_Inference_CFA-loadings", ".png"), plot = p_cfa, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigS6_Inference_CFA-diagnostics", ".png"), plot = p_MD, width = 1.5 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigS4_Inference_MGLM-diagnostics", ".png"), plot = p_diagnostics_symp, width = 2 * 105, height = 1.5 * 74.25, units = "mm")
ggplot2::ggsave(filename = paste0("./outputs/figures/png/", "FigS5_Inference_PLR-diagnostics", ".png"), plot = p_diagnostics_sscf, width = 2 * 105, height = 1.5 * 74.25, units = "mm")


# SAVE TABLES ----------------------------------------------------------------------------------------------------

write.table(t_desc, paste0("./outputs/tables/", "TabA1_Sample_description", ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(t_cca, paste0("./outputs/tables/", "TabA2_Sample_CCA-canonical-corr", ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(t_lda, paste0("./outputs/tables/", "TabA3_Sample_LDA-error-classification-table", ".txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.table(t_reg_symp, paste0("./outputs/tables/", "TabS1_Inference_MGLM-estimates", ".txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.table(t_reg_sscf, paste0("./outputs/tables/", "TabS2_Inference_PLR-estimates", ".txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.table(t_cfa, paste0("./outputs/tables/", "TabS3_Inference_CFA-estimates", ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(LRT_pvalue = t_cfa_lrt_pvalue), paste0("./outputs/tables/", "TabS4_Inference_CFA-LRT-pvalue", ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# SAVE R IMAGE (.RData) ------------------------------------------------------------------------------------------

save.image("script_objects_image.RData")