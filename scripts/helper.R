# it's a helper functions compilation

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
                        eventdata){
  
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
    left_join(select(ip_tab1_df, covariates)) %>%
    # left-join eventdata
    left_join(eventdata) %>%
    select(-`HN Number`)
  ip_icu = ip_icu %>%
    # replace timetoevent NA, impute event, right-censor, factorise
    mutate(event = ifelse(is.na(timetoevent), 0, 1),
           timetoevent = as.numeric(timetoevent),
           timetoevent = ifelse(is.na(timetoevent), 28, timetoevent),
           across(age_group:event, factor),
           age_group = ifelse(age_group != "under65", "65plus", age_group))
  
  # ip_icu %>% filter_all(any_vars(is.na(.)))
  
  # nrow(ip_icu)
  icu_hr = as.data.frame(c())
  x_factors = colnames(ip_icu)[!colnames(ip_icu) %in% c('event', 'timetoevent')]
  # iterate over cols
  for(x in x_factors){
    
    iformula <- as.formula(sprintf("Surv(as.numeric(timetoevent), as.numeric(event)) ~ %s", x))
    coxfit_dexa <- coxph(formula = iformula, 
                         data = ip_icu)
    icu_hr = bind_rows(icu_hr, tidy(coxfit_dexa, exponentiate = T, conf.int = T, conf.level = .95))
  }
  icu_hr 
  ip_tab_icu = CreateTableOne(data = select(ip_icu, -timetoevent), strata = c("event"), test = F)
  print(ip_tab_icu, smd = T)
  list(icu_hr, ip_tab_icu)
}
