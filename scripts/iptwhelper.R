library(WeightIt)
library(dplyr)
library(tableone)
library(survival)
library(survey)
library(broom)
library(boot)

iptw_tableone = function(df, 
                         treatment, 
                         x_factors, 
                         x_factors_tableone,
                         output = T,
                         path_of_unweighted_tableone = NULL,
                         path_of_weighted_tableone = NULL){
  #' IPTW and table one generation (unweighted and weighted), save to 
  #'
  #' @param df Dataframe. A dataframe that contains the treatment and x_factors
  #' @param treatment string. A string that specifies the treatment column
  #' @param x_factors c(string). A vector of strings that specfiy the covariate columns passed to IPTW
  #' @param x_factors_tableone c(string). A vector of strings that specify the variable columns shown on table one
  #' @param output (logical). A flag of whether needing tableone outputs
  #' @param path_of_unweighted_tableone string. 
  # A string that specifies the path to save unweighted table one in .xlsx
  #' @param path_of_weighted_tableone string. 
  # A string that specifies the path to save weighted table one in .xlsx
  #' @return Weightit object/output that contains the weights and data, and other parameters
  
  # weightit 
  if(length(unique(df %>% pull(c(treatment)))) > 2){
    function_call<-paste0("weight <- weightit(",treatment," ~ ",paste(x_factors, collapse = "+"), 
                          ", data = df, method = \"ps\", focal = \"none\",
                estimand = \"ATT\")"
    )
  }else{
    function_call<-paste0("weight <- weightit(",treatment," ~ ",paste(x_factors, collapse = "+"), 
                          ", data = df, method = \"cbps\", 
                estimand = \"ATT\")"
    )
  }
  
  eval(parse(text = function_call))
  # browser()
  if(output){
    # unweighted table one
    tab1_unweighted = CreateTableOne(vars = x_factors_tableone,
                                     data = df, 
                                     strata = treatment, 
                                     test = F)
    
    tab1_unweighted_Mat <- print(tab1_unweighted, quote = FALSE, 
                                 noSpaces = TRUE, printToggle = FALSE, smd = T)
    tab1_unweighted_Mat = as.data.frame(cbind(" " = rownames(tab1_unweighted_Mat),tab1_unweighted_Mat))
    ## Save to a xlsx
    writexl::write_xlsx(tab1_unweighted_Mat, path = path_of_unweighted_tableone)
    
    # IPTW-ed table one
    clus <- svydesign(ids = ~ 1, 
                      weights = ~ weight$weights,
                      data = select(df, c(treatment, x_factors_tableone)))
    tab1_weighted = svyCreateTableOne(vars = x_factors_tableone,
                                      data = clus, 
                                      strata = treatment, 
                                      test = F)
    tab1_weighted_Mat <- print(tab1_weighted, quote = FALSE, 
                               noSpaces = TRUE, printToggle = FALSE, smd = T)
    tab1_weighted_Mat = as.data.frame(cbind(" " = rownames(tab1_weighted_Mat),tab1_weighted_Mat))
    ## Save to a xlsx
    writexl::write_xlsx(tab1_weighted_Mat, path = path_of_weighted_tableone)
  }
  
  
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

iptw_hr_tidy = function(df, 
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
  
  iformula <- as.formula(paste("Surv(as.numeric(",daystillevent,"),as.numeric(", event,")) ~ ", treatment, collapse = ''))
  hr.df = c()
  
  for(x in x_factors){
    if(x == "overall"){
      # group by nothing 
      group_var = c()
    }else{
      # group by x, or segment
      group_var = c(x)
    }
    hr.df.small = df %>%
      group_by_at(group_var) %>%
      do(models = tidy(
        coxph(formula = iformula, 
              data = .,
              weights = weight), 
        exponentiate = T, conf.int = T, conf.level = .95))%>%
      unnest(models) 
    
    if(x == "overall"){
      hr.df.small = hr.df.small %>% 
        mutate(segment = "overall")
      
    }else{
      hr.df.small = hr.df.small %>% 
        mutate(segment = paste0(x, !!!syms(x))) %>%
        select(-x) 
    }
    # browser()
    hr.df.small = hr.df.small %>%
      relocate(segment, .before = term) %>%
      select(-term)
    hr.df = rbind(hr.df, hr.df.small)
  }
  
  hr.df
}


iptw_hr_boot = function(df, 
                        x_factors,
                        x_factors_weighting,
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
  
  iformula <- as.formula(paste("Surv(as.numeric(",daystillevent,"),as.numeric(", event,")) ~ ", treatment, collapse = ''))
  hr.df = c()
  
  for(x in x_factors){
    if(x == "overall"){
        group_var = c()
      }else{
        group_var = c(x)
      }
    hr.df.small = df %>%
      group_by_at(group_var) %>%
      nest() %>%
      mutate(boot.obj = map(.x = data,
                            .f = ~boot(data = .x,
                                       statistic = iptw.coxph.tidy,
                                       R = 200,
                                       stype = "i",
                                       parallel = "multicore",
                                       ncpus = 15,
                                       # coxph.tidy arguments
                                       formula = iformula,
                                       treatment = treatment,
                                       x_factors = x_factors_weighting
                                       ))
      ) %>%
      select(-data) %>%
      mutate(booted.ci = map(.x = boot.obj,
                             .f = ~booting.ci(.x))) %>%
      ungroup() %>%
      unnest_wider(booted.ci) %>%
      ungroup()
    if(x == "overall"){
      hr.df.small = hr.df.small %>% 
        mutate(segment = "overall")
      
    }else{
      hr.df.small = hr.df.small %>% 
        mutate(segment = paste0(x, !!!syms(x))) %>%
        select(-x) 
    }
    hr.df.small = hr.df.small %>%
      ungroup() %>%
      rename("estimate" = !!names(.[2]),
             "conf.low" = !!names(.[3]),
             "conf.high" = !!names(.[4])) %>%
      
      select(-boot.obj) 
    hr.df = rbind(hr.df, hr.df.small)
  }
  
  hr.df
}


# boostrap ci helpers

iptw.coxph.tidy = function(df, i, formula, x_factors, treatment){
  # browser()
  # re-iptw
  
  function_call<-paste0("psweight <- weightit(",treatment," ~ ",paste(x_factors, collapse = "+"), 
                        ", data = df[i,], method = \"ps\", 
                estimand = \"ATT\")"
  )
  eval(parse(text = function_call))
  # browser()
  df = df[i,] %>%
    mutate(weight = psweight$weights)
  
  coxph.obj = coxph(formula = formula, 
        data = df,
        weights = weight)
  
  results = tidy(coxph.obj, exponentiate = T, conf.int = T, conf.level = .95)
  results = results$estimate
  return(results)
}

booting.ci = function(boot.obj){
  booted.ci = boot.ci(boot.obj, conf = .95, type = c("perc"))
  estimate = booted.ci$t0
  ci = booted.ci$percent[4:5]
  # add col names to results for return
  return(c(estimate,ci))
}
par.boot.coxph = function(nested.df, iformula){
  return(nested.df %>%
    # partition(cluster) %>%
    mutate(boot.obj = map(.x = data,
                          .f = ~boot(data = .x,
                                     statistic = coxph.tidy,
                                     R = 200,
                                     stype = "i",
                                     # coxph.tidy arguments
                                     formula = iformula,
                                     parallel = "multicore",
                                     ncpus = 15))
    ))
}

