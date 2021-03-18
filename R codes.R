## Group 13

## Question 1
bmt.2 <- read.csv("~/Documents/MS-2nd-year/Winter/BIOST 537/Assignments/Group project/bmt-2.csv")
s.data<-with(data, Surv(tdfs, deltadfs))
surv<-survfit(s.data~1,data,conf.type = "log-log", conf.int=0.9)
surv
plot(surv, 
     xlab = "Time (in days)", ylab = "Survival probability", 
     main = "Figure 1. Kaplan-Meier estimates of disease-free survival")


## Question 2: Table 1

# Rename variables
data_names <- data %>% rename(`Patient Age` = age, `Donor Age` = donorage, Sex = male, `Donor Sex`= donormale, CMV = cmv, `Donor CMV` = donorcmv,
                              `Wait Time` = waittime, Methotrexate = mtx, FAB = fab,  Hospital = hospital, aGVHD = deltaa, Survival = deltadfs,
                              `Disease Group` = disgroup)

data_names$Sex <- ifelse(data_names$Sex == 1, "Male", "Female")
data_names$`Donor Sex` <- ifelse(data_names$`Donor Sex` == 1, "Male", "Female")
data_names$CMV <- ifelse(data_names$CMV == 1, "Positive", "Negative")
data_names$`Donor CMV` <- ifelse(data_names$`Donor CMV` == 1, "Positive", "Negative")
data_names$Methotrexate <- ifelse(data_names$Methotrexate == 1, "Yes", "No")
data_names$Survival <- ifelse(data_names$Survival == 1, "Dead or Relapsed", "Alive and Disease Free")
data_names$FAB <- ifelse(data_names$FAB == 1, "FAB Grade 4 or 5 and AML", "Otherwise")
data_names$Hospital[data_names$Hospital == 1] <- "Ohio State"
data_names$Hospital[data_names$Hospital == 2] <- "Alfred"
data_names$Hospital[data_names$Hospital == 3] <- "St. Vincent"
data_names$Hospital[data_names$Hospital == 4] <- "Hahnemann"
data_names$aGVHD <- ifelse(data_names$aGVHD == 1, "Yes", "No")
data_names$`Disease Group`[data_names$`Disease Group` == 1] <- "ALL"
data_names$`Disease Group`[data_names$`Disease Group` == 2] <- "AML Low Risk"
data_names$`Disease Group`[data_names$`Disease Group` == 3] <- "AML High Risk"

# Declare summary variables and desired summary
s_vars <-c("Patient Age" = "mean_se", "Sex" = "count_perc", "CMV" = "count_perc", "Wait Time" = "mean_se",
           "Methotrexate" = "count_perc", "Disease Group" = "count_perc", "Donor Age" = "mean_se",
           "Donor Sex" = "count_perc", "Donor CMV" = "count_perc")
# Make table
tbl_one <-make_table_one(data_names, summary_vars = s_vars,split_across = "FAB",num_rnd = 1)

knitr::kable(tbl_one, row.names = FALSE)
write.csv(tbl_one, "project_2_fab.csv")

## Make second table for disease groups
s_vars_dg <-c("Patient Age" = "mean_se", "Sex" = "count_perc","CMV" = "count_perc", "Wait Time" = "mean_se",
              "Methotrexate" = "count_perc", "FAB" = "count_perc", "Donor Age" = "mean_se",
              "Donor Sex" = "count_perc", "Donor CMV" = "count_perc")

tbl_one_dg <-make_table_one(data_names, summary_vars = s_vars_dg,split_across = "Disease Group",num_rnd = 1)
knitr::kable(tbl_one_dg, row.names = FALSE)
write.csv(tbl_one_dg, "project_2_dg.csv")

## Question 3
# Model for each variable
# Age
cox_age = coxph(s.data ~ age + cmv, data = data)
summary(cox_age)
# Sex
cox_sex = coxph(s.data ~ male + cmv, data = data)
summary(cox_sex)
# CMV
cox_cmv = coxph(s.data ~ cmv, data = data)
summary(cox_cmv)
# Donor Age
cox_donorage = coxph(s.data ~ donorage + cmv, data = data)
summary(cox_donorage)
# Donor Sex
cox_donorsex = coxph(s.data ~ donormale + cmv, data = data)
summary(cox_donorsex)
# Donor CMV
cox_donorcmv = coxph(s.data ~ donorcmv + cmv, data = data)
summary(cox_donorcmv)
# Wait time
cox_waittime = coxph(s.data ~ waittime + cmv, data = data)
summary(cox_waittime)
# Disease Group
cox_dg = coxph(s.data ~ as.factor(disgroup) + cmv, data = data)
summary(cox_dg)
# FAB
cox_fab = coxph(s.data ~ fab + cmv, data = data)
summary(cox_fab)
# MTX
cox_mtx = coxph(s.data ~ mtx + cmv, data = data)
summary(cox_mtx)
# Full model
cox = coxph(s.data ~ age + male + cmv + donorage + donormale + donorcmv + waittime + as.factor(disgroup) + fab +  mtx, data = data)
summary(cox)
# Combined significant factors
cox_fab_dg = coxph(s.data ~ as.factor(disgroup) + fab + cmv, data = data)
summary(cox_fab_dg)
# Test remaining factors
# Age
cox_fab_dg_age = coxph(s.data ~ as.factor(disgroup) + fab + age + cmv, data = data)
anova(cox_fab_dg_age, cox_fab_dg)
# Sex
cox_fab_dg_sex = coxph(s.data ~ as.factor(disgroup) + fab + male + cmv, data = data)
anova(cox_fab_dg_sex, cox_fab_dg)
# Donor Age
cox_fab_dg_dage = coxph(s.data ~ as.factor(disgroup) + fab + donorage + cmv, data = data)
anova(cox_fab_dg_dage, cox_fab_dg)
# Donor Sex
cox_fab_dg_dsex = coxph(s.data ~ as.factor(disgroup) + fab + donormale + cmv, data = data)
anova(cox_fab_dg_dsex, cox_fab_dg)
# Donor CMV
cox_fab_dg_dcmv = coxph(s.data ~ as.factor(disgroup) + fab + donorcmv + cmv, data = data)
anova(cox_fab_dg_dcmv, cox_fab_dg)
# Wait time
cox_fab_dg_wait = coxph(s.data ~ as.factor(disgroup) + fab + waittime + cmv, data = data)
anova(cox_fab_dg_wait, cox_fab_dg)
# Mtx
cox_fab_dg_mtx = coxph(s.data ~ as.factor(disgroup) + fab + mtx + cmv, data = data)
anova(cox_fab_dg_mtx, cox_fab_dg)

# 90% Confidence intervals
cox_coef <- summary(cox_fab_dg)$coef[,1]
cox_se <- summary(cox_fab_dg)$coef[,3]
exp(cox_coef - qnorm(.95)*cox_se)
exp(cox_coef + qnorm(.95)*cox_se)

## Table One Function
make_table_one <- function(data, summary_vars, 
                           split_across = NULL,
                           num_rnd = 1){
  require(magrittr)
  require(dplyr)
  options(dplyr.summarise.inform = FALSE)
  mean_se <- function(vec, sig_figs = 1){
    vals = paste0(formatC(mean(vec, na.rm = TRUE), 
                          digits = sig_figs, format = "f"), " (",
                  formatC(sd(vec, na.rm = TRUE),
                          digits = sig_figs, format = "f"), ")")
    return(data.frame("Description" = "Mean (SD)",
                      "num" = vals, stringsAsFactors = FALSE))
  }
  
  median_iqr <- function(vec, sig_figs = 1){
    vals = paste0(formatC(median(vec, na.rm = TRUE),
                          digits = sig_figs, format = "f"), " (",
                  paste0(formatC(
                    quantile(vec, probs = c(0.25, 0.75),
                             na.rm = TRUE),
                    digits = sig_figs, format = "f"),
                    collapse = "-"), ")")
    return(data.frame("Description" = "Median (IQR)",
                      "num" = vals, stringsAsFactors = FALSE))
  }
  
  mis_perec <- function(vec, sig_figs = 1){
    vals <- vec
    vals <- paste0(sum(is.na(vals)), " (",
                   formatC( 
                     100 * mean(is.na(vals)) ,
                     digits = sig_figs, format = "f"), "%)")
    return(data.frame("Description" = "Number Missing (%)",
                      "num" = vals, stringsAsFactors = FALSE))
  }
  
  count_perc <- function(vec, sig_figs = 1){
    df <- data.frame("value" = vec, stringsAsFactors = FALSE)
    if (class(df$value) == "character") {
      df$value[is.na(df$value)] <- "Missing"
    }
    df <- df %>% group_by(value, .drop = FALSE) %>%
      summarise(n = n()) %>% 
      mutate(n = n,
             perc = 100 * n/sum(n, na.rm = TRUE), 
             both = paste0(n, " (", 
                           formatC(perc, digits = sig_figs,
                                   format = "f"), "%)"))
    return(df %>% select("Description" = value, "num" = both))
  }
  
  guess_type <- function(data, vec){
    n_vec <- vec
    cnames <- names(vec)
    w_guess <- which(cnames == "")
    if (length(w_guess) < 1) {
      return(n_vec)
    }else{
      names(n_vec)[w_guess] <- n_vec[w_guess]
      for (n_idx in 1:length(w_guess)) {
        s_dat <- data[, vec[w_guess[n_idx]], drop = TRUE]
        if (length(unique(s_dat)) > 10 & 
            class(s_dat) %in% c("numeric", "integer")) {
          n_vec[w_guess[n_idx]] <- "cont"
        }else{
          n_vec[w_guess[n_idx]] <- "cat"
        }
      }
      return(n_vec)
      
    }
  }
  
  summary_funcs <- list(
    "cat" = count_perc,
    "cont" = mean_se,
    "count_perc" = count_perc,
    "mean_se" = mean_se,
    "median_iqr" = median_iqr,
    "missing_perc" = mis_perec
  )
  
  up_s_vars <- guess_type(data, summary_vars)
  sub_data <- data[, c(names(up_s_vars), split_across), drop = FALSE]
  if (!is.null(split_across)) {
    colnames(sub_data) <- c(colnames(sub_data)[-ncol(sub_data)],
                            "split_var")
  }
  which_cat <- grep("cat", up_s_vars)
  for (cat_idx in seq_along(which_cat)) {
    sub_data[, which_cat[cat_idx]] <- 
      as.factor(sub_data[, which_cat[cat_idx]])
  }
  if (!is.null(split_across)) {
    counts <- sub_data %>% group_by(split_var) %>% 
      summarise(n = n())
    counts$split_var <- as.character(counts$split_var)
    counts <- rbind(counts, 
                    cbind(split_var = "Overall", 
                          n = sum(counts$n))) %>%
      mutate(counts = paste0("(N=", n, ")"))
    counts <- as.data.frame(counts,
                            stringsAsFactors = FALSE )
  }else{
    counts <- nrow(sub_data)
    counts <- cbind(split_var = "Overall", 
                    n = counts) %>%
      mutate(counts = paste0("(N=", n, ")"))
  }
  
  new_counts <- data.frame(counts$split_var,
                           counts$counts, 
                           stringsAsFactors = FALSE) %>% t()
  colnames(new_counts) <- 1:ncol(new_counts) 
  fin_all <- cbind(
    "Variable" = " ", "Description" = " ",
    new_counts
  )
  
  colnames(fin_all)[-(1:2)] <- counts$split_var
  fin_all <- as.data.frame(fin_all, 
                           stringsAsFactors = FALSE)
  for (var_idx in 1:length(up_s_vars)) {
    sum_fn <- summary_funcs[[up_s_vars[var_idx]]]
    for (col_idx in 1:nrow(counts)) {
      if (col_idx == nrow(counts)) {
        smal_dat <- sub_data
      }else{
        smal_dat <- sub_data %>%
          filter(split_var == counts$split_var[col_idx])
      }
      out_unit <- smal_dat[, names(up_s_vars)[var_idx]] %>% 
        sum_fn(sig_figs = num_rnd)
      if (col_idx == 1) {
        fin_row <- out_unit
      }else{
        fin_row <- left_join(fin_row, out_unit, by = "Description")
      }
    }
    fin_row <- cbind(
      "Variable" = rep(names(up_s_vars)[var_idx], nrow(fin_row)),
      fin_row
    ) %>% as.data.frame(stringsAsFactors = FALSE)
    colnames(fin_row)[-(1:2)] <- counts$split_var
    fin_all <- rbind(fin_all, fin_row)
  }
  return(fin_all[-1, ])
}

## Question 4
## Effect of aGVHD on disease-free survival
bmt.tvc1<-tmerge(data1 = bmt.2,
                data2 = bmt.2,
                id=id,
                death=event(tdfs,deltadfs),
                postgvhd=tdc(ta))
bmt.tvc1[bmt.tvc1$id<5,]
surv.bmt.tvc1<-Surv(bmt.tvc1$tstart, bmt.tvc1$tstop, bmt.tvc1$deltadfs)
cox.bmt.tvc1<-coxph(surv.bmt.tvc1~postgvhd, data = bmt.tvc1)
cox.bmt.tvc1
exp(confint(cox.bmt.tvc1,level = 0.95))
cox.bmt.tvc1.adjust<-coxph(surv.bmt.tvc1~postgvhd+donormale+cmv+
                              as.factor(disgroup)+fab+mtx, data = bmt.tvc1)
cox.bmt.tvc1.adjust
exp(confint(cox.bmt.tvc1.adjust,level = 0.90))

## Effect of aGVHD on relapse
bmt.tvc2<-tmerge(data1 = bmt.2,
                 data2 = bmt.2,
                 id=id,
                 death=event(tdfs,deltar),
                 postgvhd=tdc(ta))
cox.bmt.tvc2<-coxph(surv.bmt.tvc2~postgvhd, data = bmt.tvc2)
cox.bmt.tvc2
exp(confint(cox.bmt.tvc2,level = 0.90))
cox.bmt.tvc2.adjust<-coxph(surv.bmt.tvc2~postgvhd+donormale+cmv+
                             as.factor(disgroup)+fab+mtx, data = bmt.tvc2)
cox.bmt.tvc2.adjust
exp(confint(cox.bmt.tvc2.adjust,level = 0.90))

## Question 5
bmt_gvhd <- subset(bmt,deltaa==1)
s.surv.gvhd <- with(bmt_gvhd,Surv(tdfs,deltadfs))

## full model
model5 <- coxph(s.surv.gvhd~age+donorage+male+donormale+cmv+donorcmv+waittime+
                  as.factor(disgroup)+fab+mtx,data=bmt_gvhd)
summary(model5)

## one variable in one model
coxph(s.surv.gvhd~age,data=bmt_gvhd)
coxph(s.surv.gvhd~donorage,data=bmt_gvhd)
coxph(s.surv.gvhd~male,data=bmt_gvhd)
coxph(s.surv.gvhd~donormale,data=bmt_gvhd)
coxph(s.surv.gvhd~cmv,data=bmt_gvhd)
coxph(s.surv.gvhd~donorcmv,data=bmt_gvhd)
coxph(s.surv.gvhd~waittime,data=bmt_gvhd)
coxph(s.surv.gvhd~as.factor(disgroup),data=bmt_gvhd)
coxph(s.surv.gvhd~fab,data=bmt_gvhd)
coxph(s.surv.gvhd~mtx,data=bmt_gvhd)

## combined model with variables p>0.1
m1<-coxph(s.surv.gvhd~donorage+as.factor(disgroup)+waittime+mtx,data=bmt_gvhd)
summary(m1)

## log-rank test for remaining variables
m2<-coxph(s.surv.gvhd~donorage+disgroup+waittime+mtx+age,data=bmt_gvhd)
anova(m1,m2)
m3<-coxph(s.surv.gvhd~donorage+disgroup+waittime+mtx+male,data=bmt_gvhd)
anova(m1,m3)
m4<-coxph(s.surv.gvhd~donorage+disgroup+waittime+mtx+donormale,data=bmt_gvhd)
anova(m1,m4)
m5<-coxph(s.surv.gvhd~donorage+disgroup+waittime+mtx+donorcmv,data=bmt_gvhd)
anova(m1,m5)
m6<-coxph(s.surv.gvhd~donorage+disgroup+waittime+mtx+donorcmv,data=bmt_gvhd)
anova(m1,m6)
m7<-coxph(s.surv.gvhd~donorage+disgroup+waittime+mtx+fab,data=bmt_gvhd)
anova(m1,m7)
m8<-coxph(s.surv.gvhd~donorage+disgroup+waittime+mtx+mtx,data=bmt_gvhd)
anova(m1,m8)

## Question 6
s.gvhd <- with(bmt,Surv(ta,deltaa))

## fit parametric models to check the distribution of survival time
fitweibull <- flexsurvreg(s.gvhd~1,data=bmt,dist = "weibull")
fitggamma <- flexsurvreg(s.gvhd~1,data=bmt,dist = "gengamma")
fitexp <- flexsurvreg(s.gvhd~1,data=bmt,dist = "exp")

## plot curves for non-paramatric and paramatric models
plot(survfit(s.gvhd~1,data=bmt),conf.int = FALSE,mark.time = FALSE,
     xlab = "time (in days)", ylab="survival probability",lwd=1.5)
lines(fitweibull,col=4,ci=FALSE,lwd=1.8,lty=3)
lines(fitggamma,col=2,ci=FALSE,lwd=1.8,lty=3)
lines(fitexp,col=3,ci=FALSE,lwd=1.8,lty=3)

legend("bottomright",legend = c("exponential distribution",
                                "Weibull distribution","generalized gamma distribution",
                                "nonparametric estimator"),fill = c(3,4,2,1))

# unadjusted model cox
model6.unadjust <- coxph(s.gvhd~mtx,data=bmt)
summary(model6.unadjust)
exp(confint.default(model6.unadjust,level = 0.90))

# adjusted model cox
model6.full <- coxph(s.gvhd~mtx+age+donorage+male+donormale+as.factor(disgroup),data=bmt)
model6.full

weibull.aft.adjust <- flexsurvreg(s.gvhd~mtx+age+donorage+male+donormale+as.factor(disgroup),data=bmt,dist = "weibull")
weibull.aft.adjust

## Question 7
## Effect of recovery of normal platelet levels on disease-free survival
bmt.tvc3<-tmerge(data1 = bmt.2,
                 data2 = bmt.2,
                 id=id,
                 death=event(tdfs,deltadfs),
                 postrnpl=tdc(tp))
surv.bmt.tvc3<-Surv(bmt.tvc3$tstart, bmt.tvc3$tstop, bmt.tvc3$deltadfs)
bmt.tvc3[bmt.tvc3$id<6,]
cox.bmt.tvc3<-coxph(surv.bmt.tvc3~postrnpl, data = bmt.tvc3)
cox.bmt.tvc3
exp(confint(cox.bmt.tvc3, level = 0.9))
cox.bmt.tvc3.adjust<-coxph(surv.bmt.tvc3~postrnpl+donormale+cmv+
                             as.factor(disgroup)+fab+mtx, data = bmt.tvc3)
cox.bmt.tvc3.adjust
exp(confint(cox.bmt.tvc3.adjust, level = 0.9))

cox.bmt.tvc3.curve<-survfit(surv.bmt.tvc3~postrnpl, data = bmt.tvc3)
## Effect of recovery of normal platelet levels on relapse
bmt.tvc4<-tmerge(data1 = bmt.2,
                 data2 = bmt.2,
                 id=id,
                 death=event(tdfs,deltar),
                 postrnpl=tdc(tp))
surv.bmt.tvc4<-Surv(bmt.tvc4$tstart, bmt.tvc4$tstop, bmt.tvc4$deltar)
cox.bmt.tvc4<-coxph(surv.bmt.tvc4~postrnpl, data = bmt.tvc4)
cox.bmt.tvc4
exp(confint(cox.bmt.tvc4,level = 0.9))
cox.bmt.tvc4.adjust<-coxph(surv.bmt.tvc4~postrnpl+donormale+cmv+
                             as.factor(disgroup)+fab+mtx, data = bmt.tvc4)
cox.bmt.tvc4.adjust
exp(confint(cox.bmt.tvc4.adjust,level = 0.9))
cox.bmt.tvc4.curve<-survfit(surv.bmt.tvc4~postrnpl, data = bmt.tvc4)
