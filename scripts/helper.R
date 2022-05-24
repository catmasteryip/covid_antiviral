# it's a helper functions compilation
library(WeightIt)
library(dplyr)
library(tableone)

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
  tab1_unweighted = CreateTableOne(data = select(df, c(treatment, x_factors)), 
                                   strata = treatment, 
                                   test = F)
  # unweighted table one
  tab1_unweighted_Mat <- print(tab1_unweighted, quote = FALSE, 
                               noSpaces = TRUE, printToggle = FALSE, smd = T)
  tab1_unweighted_Mat = as.data.frame(cbind(" " = rownames(tab1_unweighted),tab1_unweighted))
  ## Save to a xlsx
  writexl::write_xlsx(tab1_unweighted_Mat, path = path_of_unweighted_tableone)
  
  # IPTW-ed table one
  tab1_weighted = svyCreateTableOne(data = clus, strata = treatment, 
                                    test = F)
  tab1_weighted_Mat <- print(tab1_weighted, quote = FALSE, 
                             noSpaces = TRUE, printToggle = FALSE, smd = T)
  tab1_weighted_Mat = as.data.frame(cbind(" " = rownames(tab1_weighted),tab1_weighted))
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
  
  clus <- svydesign(id =~ 1, weights = weight$weights, data = ip_propensity)
  res <- svyglm(paste0(outcome, " ~ ", treatment), design = clus,family = binomial)
  tidy_reg = tidy(res, exponentiate = T, conf.int = T, conf.level = .95)
  return(tidy_reg %>% filter(term != "(Intercept)"))
}


impute_meds = function(df){
  #' impute medication status in binary 
  #'
  #' @param df Dataframe. A dataframe that contains the Drug Name column from raw data
  #' @return Dataframe that mutated new cols of binary meds status
  df = df %>%
    # case-sensitive imputation of drug names
    mutate(
      # paxlovid
      paxlovid = ifelse(grepl("PAXLOVID", `Drug Name`), 1, 0),
      # aspirin
      aspirin = ifelse(grepl("ASPIRIN|ACETYLSALICYLIC ACID", `Drug Name`), 1, 0),
      # ACE inhibitors
      ace_inhibitors = ifelse(grepl("BENAZEPRIL|CAPTOPRIL|CLIZAPRIL|ENALAPRIL|ENALAPRILAT|FOSINOPRIL|LISINOPRIL|MOEXIPRIL|PERINDOPRIL|QUINAPRIL|RAMIPRIL|TRANDOLAPRIL", `Drug Name`), 1, 0),
      # beta blockers
      beta_blockers = ifelse(grepl("ACEBUTALOL|ATENOLOL|BETAXOLOL|BISOPROLOL|CARTEOLOL|CARVEDILOL|ESMOLOL|LABETALOL|LEVOBUNOLOL|METOPROLOL|NADOLOL|NEBIVOLOL|PINDOLOL|PROPANOLOL|SOTALOL|TIMOLOL", `Drug Name`), 1, 0),
      # calcium channel blockers
      calcium_channel_blockers = ifelse(grepl("AMLODIPINE|CLEVIDIPINE|FELODIPINE|FLUNARZINE|ISRADIPINE|LEVAMLODIPINE|NICARDIPINE|NIFEDIPINE", `Drug Name`), 1, 0),
      # statins
      statins = ifelse(grepl("ATORVASTATIN|SIMVASTATIN|FLUVASTATIN|LOVASTATIN|PITAVASTATIN|PRAVASTATIN|ROSUVASTATIN", `Drug Name`), 1, 0),
      # molnupiravir
      molnupiravir = ifelse(grepl("MOLNUPIRAVIR", `Drug Name`), 1, 0),
      
    ) %>% 
    mutate_at(c("paxlovid", "aspirin","ace_inhibitors","beta_blockers","calcium_channel_blockers", "statins","molnupiravir"), ~replace_na(., 0))
  return(df)
}
# severe covid related meds
impute_meds_severe_covid = function(df){
  #' impute medication status in binary 
  #'
  #' @param df Dataframe. A dataframe that contains the Drug Name column from raw data
  #' @return Dataframe that mutated new cols of binary meds status
  df = df %>%
    # case-sensitive imputation of drug names
    mutate(
      # baricitinib 
      baricitinib = ifelse(grepl("BARICITINIB", `Drug Name`), 1, 0)
    ) %>% 
    mutate_at(c("baricitinib"), ~replace_na(., 0))
  return(df)
}


# --------
impute_comorbs = function(df, diag_cols){
  #' impute comorbidities status in binary 
  #'
  #' @param df Dataframe. A dataframe that contains the icd9 code columns from raw data at diag_cols column indexes
  #' @param diag_cols. A vector of column indexes of icd9 code columns
  #' @return Dataframe that mutated new cols of binary comorbidity status with no icd9 code
  df = df %>%
  mutate(
  diabetes = ifelse(if_any(diag_cols, ~grepl("^250", .)), 1, 0),
  hypertension = ifelse(if_any(diag_cols,
                            ~between(., "401", "405.99")), 1, 0),
  stroke = ifelse(if_any(diag_cols,
                      ~(between(., "362.3", "362.49") |
                        between(., "430", "435.99"))), 1, 0),
  heart_failure = ifelse(if_any(diag_cols, ~grepl("^428", .)), 1, 0),
  # atrial fibrillation
  atrial_fib = ifelse(if_any(diag_cols, ~grepl("^427.31", .)), 1, 0),
  parkinsons = ifelse(if_any(diag_cols, ~grepl("^332.00", .)), 1, 0),
  # schizophrenia
  schizo = ifelse(if_any(diag_cols, ~grepl("^295", .)), 1, 0),
  # liver cirrhosis
  liver = ifelse(if_any(diag_cols,
                     ~(grepl("^571.5|^517.2|^546|^567.89|^789.59", .) |
                     between(., "456", "456.21") |
                     between(., "567", "567.21") |
                     between(., "572.2", "572.49"))), 1, 0),
  # depression
  depression = ifelse(if_any(diag_cols,
                    ~( grepl("^300.4|^269.5", .) |
                     between(., "269.2", "260.39") |
                     between(., "309", "311"))), 1, 0),
  # chronic kidney disease
  ckd = ifelse(if_any(diag_cols,
                     ~(grepl("^593.9|^592", .) |
                     between(., "583", "586.99"))), 1, 0),
  # rheumatoid arthritis
  arthritis = ifelse(if_any(diag_cols,
                         ~(grepl(., "^446.5|^725") |
                         between(., "714", "714.29") |
                         between(., "714.8", "718.89")))
                         , 1, 0),
  # obesity
  obesity =  ifelse(if_any(diag_cols,
                     ~between(., "277.70", "280.09")), 1, 0),
  # alcohol abuse
  alcohol_abuse = ifelse(if_any(diag_cols,
                             ~grepl("^303|^305", .)), 1, 0),
  
)
  # specify binary comorbity column index
  comorb_cols = 1:12 + max(diag_cols)
  # replace NA as 0
  df = df %>%
    mutate_at(comorb_cols, ~replace_na(., 0))
  # remove icd codes from df_diag
  df = df %>%
    select(-diag_cols)
  df
}
# -------------

# to be used to replace IP_data repeated calculation of HR
calculate_hr = function(df, 
                        covariates = c("HN Number",'age_group', 'Sex', 'obesity', 'diabetes', 'molnupiravir', 'paxlovid'),
                        eventdata,
                        weights = NULL){
  
  #' a helper function for univariate hazard ratio calculation and descriptive statistics table generation
  #'
  #' @param df Dataframe. A dataframe that contains the covariates
  #' @param covariates. [string]. vector of strings of col names that are covariates
  #' @param eventdata Dataframe. A dataframe that contains HN Number and time to event 
  #' with exactly mentioned column names
  #' @return (dataframe, dataframe). A list of hazard ratio results and table of descriptive stats
  #' 
  ip_icu = df %>%
    select(`HN Number`) %>%
    left_join(select(df, covariates)) %>%
    # left-join eventdata
    left_join(eventdata) %>%
    select(-`HN Number`)
  if("event" %in% colnames(eventdata)){
    # if there is event col
    ip_icu = ip_icu %>%
      # impute event, if NA, event is 0
      mutate(event = ifelse(is.na(event), 0, event))
  }else{
    # if there isnt event col
    ip_icu = ip_icu %>%
      # impute event, if timetoevent NA, event is 0
      mutate(event = ifelse(!is.na(timetoevent), 1, 0))
  }
  ip_icu = ip_icu %>%
    # replace timetoevent NA, factorise
    mutate(timetoevent = as.numeric(timetoevent),
           timetoevent = ifelse(is.na(timetoevent), 28, timetoevent),
           across(age_group:event, factor),
           age_group = ifelse(age_group != "under65", "1", "0"),
           Sex = ifelse(Sex == "M", "1", "0")) %>%
    rename(plus65 = age_group,
           SexM = Sex)
  
  # ip_icu %>% filter_all(any_vars(is.na(.)))
  
  # nrow(ip_icu)
  icu_hr = as.data.frame(c())
  # obesity is excluded due to non-convergence
  x_factors = colnames(ip_icu)[!colnames(ip_icu) %in% c('event', 'timetoevent', 'obesity', 'molnupiravir', 'paxlovid')]
  # iterate over cols
  for(x in x_factors){
    if(!is.null(weights)){
      # susbet weights
      df_hr = filter(ip_icu, eval(parse(text = paste0(x, "==1")))) %>% 
        mutate(id = 1:n()) %>% 
        select(-x) 
      coxfit_positive <- coxph(formula = Surv(as.numeric(timetoevent), as.numeric(event), type = "right") ~ paxlovid, 
                           data = df_hr,
                           weights = weights[df_hr$id],
                           model = T)
      coxfit_negative <- coxph(formula = Surv(as.numeric(timetoevent), as.numeric(event), type = "right") ~ paxlovid, 
                               data = df_hr,
                               weights = weights[df_hr$id],
                               model = T)
    }else{
      coxfit_positive <- coxph(formula = Surv(as.numeric(timetoevent), as.numeric(event), type = "right") ~ paxlovid, 
                           data = df_hr,
                           model = T)
      coxfit_negative <- coxph(formula = Surv(as.numeric(timetoevent), as.numeric(event), type = "right") ~ paxlovid, 
                               data = df_hr,
                               model = T)
    }
    # iformula <- as.formula(sprintf("Surv(as.numeric(timetoevent), as.numeric(event)) ~ %s", x))
    
    icu_hr = bind_rows(icu_hr, 
                       tidy(coxfit_positive, exponentiate = T, conf.int = T, conf.level = .95) %>%
                         mutate(xfactor = x,
                                cox_model = list(coxfit_dexa)),
                       tidy(coxfit_negative, exponentiate = T, conf.int = T, conf.level = .95) %>%
                         mutate(xfactor = x,
                                cox_model = list(coxfit_dexa)))
  }
  icu_hr 
  ip_tab_icu = CreateTableOne(data = select(ip_icu, -timetoevent), strata = c("event"), test = F)
  print(ip_tab_icu, smd = T)
  list(icu_hr, ip_tab_icu, ip_icu)
}

