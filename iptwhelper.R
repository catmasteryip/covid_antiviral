library(WeightIt)
library(dplyr)
library(tableone)
library(survival)
library(survey)
library(broom)

iptw_tableone = function(df, 
                         treatment, 
                         x_factors, 
                         path_of_unweighted_tableone,
                         path_of_weighted_tableone){
  #' IPTW and table one generation (unweighted and weighted), save to 
  #'
  #' @param df Dataframe. A dataframe that contains the treatment and x_factors
  #' @param treatment string. A string that specifies the treatment column
  #' @param x_factors c(string). A vector of strings that specfiy the covariate columns passed to IPTW
  #' @param path_of_unweighted_tableone string. 
  # A string that specifies the path to save unweighted table one in .xlsx
  #' @param path_of_weighted_tableone string. 
  # A string that specifies the path to save weighted table one in .xlsx
  #' @return Weightit object/output that contains the weights and data, and other parameters
  
  # weightit 
  function_call<-paste0("weight <- weightit(",treatment," ~ ",paste(x_factors, collapse = "+"), 
                        ", data = df, method = \"ps\", 
                estimand = \"ATT\")"
  )
  eval(parse(text = function_call))
  
  # unweighted table one
  tab1_unweighted = CreateTableOne(data = select(df, c(treatment, x_factors)), 
                                   strata = treatment, 
                                   test = F)
  
  tab1_unweighted_Mat <- print(tab1_unweighted, quote = FALSE, 
                               noSpaces = TRUE, printToggle = FALSE, smd = T)
  tab1_unweighted_Mat = as.data.frame(cbind(" " = rownames(tab1_unweighted_Mat),tab1_unweighted_Mat))
  ## Save to a xlsx
  writexl::write_xlsx(tab1_unweighted_Mat, path = path_of_unweighted_tableone)
  
  # IPTW-ed table one
  tab1_weighted = svyCreateTableOne(data = clus, strata = treatment, 
                                    test = F)
  tab1_weighted_Mat <- print(tab1_weighted, quote = FALSE, 
                             noSpaces = TRUE, printToggle = FALSE, smd = T)
  tab1_weighted_Mat = as.data.frame(cbind(" " = rownames(tab1_weighted_Mat),tab1_weighted_Mat))
  ## Save to a xlsx
  writexl::write_xlsx(tab1_weighted_Mat, path = path_of_weighted_tableone)
  
  return(weight)
}

iptw_regression_tidy = function(df, 
                                weight, 
                                outcome, 
                                treatment){
  #' IPTW-weighted regression to generate the exponentiated coeff/odds ratio given treatment to outcome
  #'
  #' @param df Dataframe. A dataframe that is used for iptw
  #' @param weight Weightit object/output. which contains the weights and data, and other parameters
  #' @param outcome String. A string that specifies the outcome column in df
  #' @param treatment String. A string that specifies the treatment column in df
  #' @return A DataFrame that contains the exponentiated coeff/odds ratio given treatment to outcome
  
  clus <- svydesign(id =~ 1, weights = weight$weights, data = df)
  res <- svyglm(paste0(outcome, " ~ ", treatment), design = clus,family = binomial)
  tidy_reg = tidy(res, exponentiate = T, conf.int = T, conf.level = .95)
  return(tidy_reg %>% filter(term != "(Intercept)"))
}

iptw_hr_one_col = function(df, 
                        x_factors,
                        event,
                        daystillevent,
                        treatment){
  #' IPTW-weighted univariate cox regression to generate the exponentiated hazard ratio given treatment to outcome
  #'
  #' @param df Dataframe. A dataframe that is for iptw and has iptw weight in column called "weight"
  #' @param x_factors c(String). A vector of strings for category segmentation
  #' @param event String. String specifying the event column
  #' @param daystillevent String. String specifying no. of days till event
  #' @param treatment String. A string that specifies the treatment column in df
  #' @return A DataFrame that contains the exponentiated hazard ratio given treatment to outcome, 
  #' segmented by each level of each factor in x_factors 
  
  iformula <- as.formula(sprintf("Surv(as.numeric(daystillevent), as.numeric(death28d)) ~ paxlovid"))
  hr.df = c()
  
  for(x in x_factors){
    
    hr.df.small = df %>%
      group_by(!!!syms(x)) %>%
      do(models = tidy(
        coxph(formula = iformula, 
              data = .,
              weights = weight), 
        exponentiate = T, conf.int = T, conf.level = .95))%>%
      unnest(models) %>%
      mutate(variable = paste0(x, !!!syms(x))) %>%
      select(-x) %>%
      relocate(variable, .before = term) %>%
      select(-term)
    hr.df = rbind(hr.df, hr.df.small)
  }
  
  hr.df
}
