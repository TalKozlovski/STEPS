rm(list = ls())

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()


if(!require('data.table'))  install.packages('data.table')
if(!require('magrittr'))    install.packages('magrittr')
if(!require('dplyr'))       install.packages('dplyr')
if(!require('tidyr'))       install.packages('tidyr')
if(!require('MASS'))        install.packages('MASS')
if(!require('rpart'))       install.packages('rpart')
if(!require('foreach'))     install.packages('foreach')  ; library(foreach)
library(doRNG)
if(!require('doParallel'))  install.packages('doParallel')  ; library(doParallel)
if(!require('parallel'))    install.packages('parallel') ; library(parallel)
if(!require('purrr'))       install.packages('purrr') ; library(purrr)
if(!require('pROC'))        install.packages('pROC') ; library(pROC)
if(!require('plotROC'))     install.packages('plotROC') ; library(plotROC)



### Tree functions
if(!require('STEPS'))     devtools::install_github('TalKozlovski/STEPS/STEPS') ; library(STEPS)





#### Percentize a vector to valuess 0 to 1
percentize <- function(x) {
  min.x <- min(x)
  max.x <- max(x)
  out <- (x - min(x)) / (max.x -min.x)
  return(out)
}



#### Create X matrix from stress test
x_sim_fun <- function(ID_tab, levels.tab, pX) {
  tab <- levels.tab %>% 
    left_join(ID_tab)
  tab <- tab %>% 
    rowwise() %>%
    mutate("obs_severity" = case_when(sum(difficulty1, difficulty2) == 0 ~ logistic(severity01, k = 8, a = 1), 
                                      sum(difficulty1, difficulty2) == 1 ~ logistic(severity01, k = 8, a = 0.8), 
                                      sum(difficulty1, difficulty2) == 2 ~ logistic(severity01, k = 8, a = 0.6), 
                                      sum(difficulty1, difficulty2) == 3 ~ logistic(severity01, k = 8, a = 0.4), 
                                      sum(difficulty1, difficulty2) == 4 ~ logistic(severity01, k = 8, a = 0.2)))
  
  tab[, paste("x", c(1), sep = "")] <- tab$obs_severity 
  tab[, paste("x", c(2:pX), sep = "")] <- 0
  
  return(tab)
}



### Logistic function
logistic <- function(x, L = 1, k = 1, a = 0) {
  return(L / (1 + exp(-k * (x - a))))
}






SIM_ordering_setR_fun <- function(IDtimes.tab, 
                                  data) {


  # Create R groups names
  ordered_tab.R <- IDtimes.tab 
  # add index variable
  ordered_tab.R[, "index"] <- c(1:nrow(IDtimes.tab))
  ordered_tab.R[, "R.group"] <- as.character(apply(ordered_tab.R[, c("ID", "difficulty1", "difficulty2")],
                                                   1, paste, collapse = "_"))
  # aggregate R groups to summarise table
  ordered_tab.R %<>% 
    group_by(R.group) %>% 
    mutate("n.R" = n()) %>% 
    ungroup() 
  
  
  return(ordered_tab.R)
  
}



SIM_ordering_setW_fun <- function(IDtimes.tab,
                                  data, 
                                  setR) {
  tmp_data.W <- setR %>%
    dplyr::select(ID, difficulty1, difficulty2,
                  R.group, n.R)
  duped <- duplicated(tmp_data.W)
  tmp_data.W <- tmp_data.W[!duped, ]
  
  ### X is the better state and Y is the severe
  tmp_data.W1 <- tmp_data.W %>%
    full_join(tmp_data.W, by = c("ID", "difficulty1")) %>%
    filter((difficulty2.y > difficulty2.x)) %>% 
    mutate( difficulty1.y = difficulty1, 
            ID.y = ID)
  
  colnames(tmp_data.W1) <- c("ID.better", 
                             "difficulty1.better",
                             "difficulty2.better",
                             "R.group.better",
                             "n.R.better",
                             "difficulty2.severe", 
                             "R.group.severe",
                             "n.R.severe",
                             "difficulty1.severe",
                             "ID.severe")
  
  
  tmp_data.W2 <- tmp_data.W %>%
    full_join(tmp_data.W, by = c("ID", "difficulty2")) %>%
    filter((difficulty1.y > difficulty1.x)) %>%
    mutate(difficulty2.y = difficulty2,
           ID.y = ID)
  
  colnames(tmp_data.W2) <- c("ID.better",
                             "difficulty1.better",
                             "difficulty2.better",
                             "R.group.better",
                             "n.R.better",
                             "difficulty1.severe",
                             "R.group.severe",
                             "n.R.severe",
                             "difficulty2.severe",
                             "ID.severe")
  
  
  ordered_tab.W <- tmp_data.W1 %>% bind_rows(tmp_data.W2)

  ordered_tab.W <- ordered_tab.W %>%
    group_by(ID.better, ID.severe) %>%
    mutate("n.Wi" = n()) %>% 
    ungroup()
  
  return(ordered_tab.W)
}



ord_pairs_fun_oracle <- function(pair, sort.by1) {
  ## ord_pairs_fun sorts pairs of observations to the better (left) and severe (right) states 
  # Input: @ pair - a pair of two observations indexes
  #        @ sort.by  A vector with values to sort according to the given pair
  # Output: vector of length 2, sorted from better to severe disease states. 

  if(sort.by1[pair[1]] < sort.by1[pair[2]]) { 
    return(pair)
  }
  return(pair[c(2:1)])
  
}


ord_pairs_fun <- function(pair, sort.by1, sort.by2 = NULL, thresh = 1) {
  ## ord_pairs_fun sorts pairs of observations to the better (left) and severe (right) states 
  # Input: @ pair - a pair of two observations indexes
  #        @ sort.by  A vector with values to sort according to the given pair
  # Output: vector of length 2, sorted from better to severe disease states. 

  if((sort.by1[pair[1]] < (sort.by1[pair[2]] - thresh)) & (sort.by2[pair[1]] < (sort.by2[pair[2]] - thresh))) { 
    return(pair)
  }
  if(((sort.by1[pair[1]] - thresh)  > sort.by1[pair[2]]) & ((sort.by2[pair[1]] - thresh) > sort.by2[pair[2]])) { 
    return(pair[c(2:1)])
  } else {
    return(NULL) 
  }
    
}


orderIDfun <- function(data, oracle) {
  n.id <- unique(data$ID)
  pairs <- t(combn(x = c(1:length(n.id)), m = 2))
  if(oracle) {
    ord.pairs <- t(apply(pairs, 1, ord_pairs_fun_oracle, 
                         sort.by1 = data$severity01) )
  } else {
    ord.pairs <- apply(pairs, 1, ord_pairs_fun, 
                       sort.by1 = data$score1, 
                       sort.by2 = data$score2) 
    ord.pairs <- do.call('rbind', ord.pairs)
  }
  ord.pairs <- matrix(n.id[ord.pairs], ncol = 2)
  return(ord.pairs)
}

SIM_ordering_setB_fun <- function(obs_cov, 
                                  # IDtimes.tab, 
                                  # data, 
                                  setR.ord.pairs, 
                                  oracle) {
  
  
  # Ordering pairs according to observed severity, score 1 & score 2

  orderID <- orderIDfun(data = obs_cov, oracle = oracle)

  tmp.setR.ord.pairs <- setR.ord.pairs %>% 
    dplyr::select(ID, difficulty1, difficulty2, R.group, n.R) %>% 
    unique()
  
  orderID <- as.data.frame(orderID)
  colnames(orderID) <- c("ID", "ID.severe")
  
  full_tab1 <- orderID %>%
    # left_join(IDtimes.tab) %>% 
    left_join(tmp.setR.ord.pairs) %>%
    unique() %>% 
    rename(ID.better = ID, 
           ID = ID.severe,
           Group.better = R.group, 
           n.R.better = n.R
           )

  full_tab2 <- full_tab1 %>% 
    left_join(tmp.setR.ord.pairs) %>% 
    rename(ID.severe = ID, 
           Group.severe = R.group, 
           n.R.severe = n.R)

  full_tab3 <- full_tab2 %>%
    group_by(ID.severe, ID.better) %>%
    mutate("dij" = n())  %>% 
    ungroup()

  
  return(full_tab3)
}




SIM_sets_fun <- function(data, 
                         obs_cov,
                         # inds, 
                         oracle = F
) {
  
  IDtimes.tab <- data[, c("ID", "R", "difficulty1", "difficulty2")]
  # Set R: set of repetitions
  setR.ord.pairs  <- SIM_ordering_setR_fun(IDtimes.tab = IDtimes.tab,
                                           data = data)
  
  # Set W : within subject ordering
  setW.ord.pairs  <- SIM_ordering_setW_fun(IDtimes.tab = IDtimes.tab,
                                           data = data, 
                                           setR = setR.ord.pairs) 

  ## Set B: between subjects

  # Detailed Set B: between subjects
  setB.ord.pairs  <- SIM_ordering_setB_fun(obs_cov = obs_cov, 
                                           # IDtimes.tab = IDtimes.tab %>%
                                           #   dplyr::select(-R) %>% unique(),
                                           # data = data,
                                           setR.ord.pairs = setR.ord.pairs, 
                                           oracle = oracle)
  
  return(list("setB" = setB.ord.pairs, 
              "setR" = setR.ord.pairs, 
              "setW" = setW.ord.pairs))
  
}




lm_pred_score1_fun <- function(train_data, test_data) {
  mod <- lm(score1 ~. , train_data %>% dplyr::select(-ID, -severity, -obs_severity, -score2, -status, -name_var,
                                                     -difficulty1, -difficulty2))
  pred <- predict(mod, test_data)
  pred_tab <- data.frame(test_data[, c( "ID", "difficulty1", "difficulty2")], "lm_score1" = pred)
  return(pred_tab)
}

lm_pred_score2_fun <- function(train_data, test_data) {
  mod <- lm(score2 ~. , train_data %>% dplyr::select(-ID, -severity, -obs_severity, -score1, -status, -name_var,
                                                     -difficulty1, -difficulty2))
  pred <- predict(mod, test_data)
  pred_tab <- data.frame(test_data[, c( "ID", "difficulty1", "difficulty2")], "lm_score2" = pred)
  return(pred_tab)
}


addNoise <- function(data, p , noise) {
  noised_data <- data
  noised_data[, paste0("x", c(1:p))] <- noised_data[, grepl("x", colnames(noised_data))] + 
    rnorm(nrow(noised_data) * p, sd = sqrt(noise))
  return(noised_data)
}


divideTestTrain <- function(percTrain, data) {
  n <-  length(unique(data$ID))
  n.train <- floor(percTrain * n)
  n.test  <- n - n.train
  
  train.ind <- sort(sample(c(1:n), size = n.train, replace = F ))
  test.ind <- c(1:n)[!(c(1:n) %in% train.ind)]
  
  ## Divide the data to train and test
  dat_train <- data %>%
    filter(ID %in% train.ind)
  dat_test <- data %>%
    filter(ID %in% test.ind)
  
  return(list("Train" = dat_train, "Test" = dat_test, "train_index" = train.ind, "test_index" = test.ind))
}

lmDataPrep <- function(data) {
  lm.data <- data %>% 
    dplyr::select( -R) %>% 
    group_by(ID, severity, score1, score2, status, difficulty1, difficulty2 ) %>% 
    summarise_all(mean) %>% 
    ungroup()
  
  lm.data_list <- lm.data %>% 
    mutate("name_var" = paste0("M_", difficulty1, "_C_", difficulty2)) %>%
    group_split(difficulty1, difficulty2)
  
  
  lm.data_list %>%
    purrr::map(~pull(.,name_var)) %>% # Pull out Species variable
    purrr::map(~unique(.)) -> names(lm.data_list) 
  
  return(lm.data_list)
  
}




lm_pred_fun <- function(train_data, test_data, data) {
  
  lm_res_score1_list <- list(train_data = train_data,
                             test_data = test_data) %>% 
    purrr::pmap(lm_pred_score1_fun) 
  lm_res1 <- do.call('rbind', lm_res_score1_list)
  
  
  lm_res_score2_list <- list(train_data = train_data,
                             test_data = test_data) %>% 
    purrr::pmap(lm_pred_score2_fun) 
  lm_res2 <- do.call('rbind', lm_res_score2_list)
  
  lm_res <- data[, c("ID", "difficulty1", "difficulty2")] %>% 
    left_join(full_join(lm_res1, lm_res2))
  
  return(lm_res) 
}





Wpairs_conc <- function(setW, setR, Fx, by.vars, conc.name = "mean_concW") {
  Fx.mean <- tapply(X = Fx[setR$index],
                    INDEX = setR$R.group,
                    FUN = mean)
  tmp.W <- setW %>%
    left_join(data.frame("R.group.better" = names(Fx.mean),
                         "Fx.better.mean" = Fx.mean,
                         stringsAsFactors = F)) %>%
    left_join(data.frame("R.group.severe" = names(Fx.mean),
                         "Fx.severe.mean" = Fx.mean,
                         stringsAsFactors = F))
  Wconc <- mean(tmp.W$Fx.severe.mean > tmp.W$Fx.better.mean)
  return(Wconc)
  
}


Bpairs_conc <- function(setB, setR, Fx, by.vars, conc.name = "mean_concB") {
  Fx.mean <- tapply(X = Fx[setR$index],
                    INDEX = setR$R.group,
                    FUN = mean)
  tmp.B <- setB %>%
    left_join(data.frame("Group.better" = names(Fx.mean),
                         "Fx.better.mean" = Fx.mean,
                         stringsAsFactors = F)) %>%
    left_join(data.frame("Group.severe" = names(Fx.mean),
                         "Fx.severe.mean" = Fx.mean,
                         stringsAsFactors = F))
  
  Bconc <- tmp.B %>% 
    group_by_at(by.vars) %>% 
    summarise(!!conc.name :=  mean(Fx.severe.mean > Fx.better.mean)) %>% 
    ungroup()
  Bconc_mean <- mean(Bconc[, conc.name, drop = T], na.rm = T)
  return(Bconc_mean)
  
}


aucFun <- function(status, var) {
  AUC <- as.numeric(roc(status, var, direction = "<")$auc)
  return(AUC)
}



summariseIndex <- function(test_data, train_data) {
  
  
  ### Summarize STEPS predictions
  healthy.tab_STEPS_summarise <- train_data %>% 
    filter(severity01 < 0.8) %>% 
    group_by(Method, difficulty2, difficulty1) %>%
    summarise("Median_y_hat" = median(y_hat, na.rm = T), 
              "SD_scale_y_hat" = sd(y_hat, na.rm = T), 
              "n" = n())
  
  test_data_class <- test_data %>% 
    left_join(healthy.tab_STEPS_summarise) %>% 
    mutate("diff_from_median" = (y_hat - Median_y_hat)) %>% 
    mutate("positive_diff_from_median" = diff_from_median > 0,
           "greater_sd_median" = diff_from_median > (2*SD_scale_y_hat)) %>%
    group_by(Method, ID, 
             # status, 
             severity, severity01) %>% 
    summarise("mean_diff_from_median" = mean(diff_from_median), 
              "mean_trim0.05_from_median" = mean(diff_from_median, trim = 0.05), 
              "positive_mean_from_median" = sum(positive_diff_from_median * diff_from_median) / n(),
              "positive_SD_median" = sum(greater_sd_median * positive_diff_from_median * diff_from_median) / n(),
              "n" = n()) %>% 
    ungroup() %>% 
    pivot_longer(cols = c(mean_diff_from_median, mean_trim0.05_from_median, positive_mean_from_median, positive_SD_median), 
                 names_to = "index_type", values_to = "diff_values")
  return(test_data_class)
  
}




aggResults <- function(obs_cov, 
                       test_data, 
                       train_data,
                       ordered_sets, 
                       ordered_sets_oracle) {
  
  
  test_data_noR <- test_data %>% 
    group_by(ID, difficulty1, difficulty2, 
             # status, 
             severity, severity01 ) %>% 
    summarise("STEPS" = mean(STEPS), 
              "STEPS_oracle" = mean(STEPS_oracle)) %>% 
    ungroup()
  
  train_data_noR <- train_data %>% 
    group_by(ID, difficulty1, difficulty2,
             # status, 
             severity, severity01 ) %>% 
    summarise("STEPS" = mean(STEPS), 
              "STEPS_oracle" = mean(STEPS_oracle)) %>% 
    ungroup()
  
  
  test_data_noR <- test_data_noR %>% 
    pivot_longer(cols = c(STEPS, STEPS_oracle), 
                 names_to = "Method", 
                 values_to = "y_hat")
  train_data_noR <- train_data_noR %>% 
    pivot_longer(cols = c(STEPS, STEPS_oracle), 
                 names_to = "Method", 
                 values_to = "y_hat")
  
  
  test_data_summarized <- summariseIndex(test_data_noR, train_data_noR)
  
  ### Correlations withs everity
  correlations <- test_data_summarized %>%
    group_by(Method, index_type) %>% 
    summarise("cor_severity" = cor(diff_values, severity, method = "spearman")) %>% 
    ungroup() %>% 
    bind_rows(data.frame("Method" = c("score1", "score2"), 
                         "index_type" = c("", ""),
                         "cor_severity" = c(cor(obs_cov$score1, obs_cov$severity, method = "spearman"), 
                                            cor(obs_cov$score2, obs_cov$severity, method = "spearman"))))
  
  #### Classifications
  classifications <- test_data_summarized %>%
    mutate("status" = as.numeric(severity01 > 0.5)) %>%
    group_by(Method, index_type) %>% 
    summarise("AUC" = aucFun(status, diff_values)) %>% 
    ungroup() %>% 
    bind_rows(data.frame("Method" = c("score1", "score2"), 
                         "index_type" = c("", ""),
                         "AUC" = c(aucFun(as.numeric(obs_cov$severity01 > 0.5), obs_cov$score1), 
                                   aucFun(as.numeric(obs_cov$severity01 > 0.5), obs_cov$score2))))
  
  #### Concordance 
  concordances <- concResultsFun(ordered_sets = ordered_sets, ordered_sets_oracle = ordered_sets_oracle, STEPS_output = test_data)
  return(list("resCor" = correlations, "resClass" = classifications, "resConc" = concordances))
}


### Concordance
concResultsFun <- function(ordered_sets, ordered_sets_oracle, STEPS_output) {
  
  #### Conc between subjects
  ## STEPS 
  idComp.STEPS.observed <- Bpairs_conc(setB = ordered_sets[["setB"]], 
                                       setR = ordered_sets[["setR"]], 
                                       Fx = STEPS_output$STEPS, 
                                       by.vars = c("ID.better", "ID.severe"))
  idComp.STEPS.oracle <- Bpairs_conc(setB = ordered_sets_oracle[["setB"]], 
                                     setR = ordered_sets_oracle[["setR"]], 
                                     Fx = STEPS_output$STEPS_oracle, 
                                     by.vars = c("ID.better", "ID.severe"))
  
  STEPS.Rsd <- STEPS_output %>% 
    group_by(ID, difficulty1, difficulty2) %>% 
    summarise("sd_STEPS" = sd(STEPS), "sd_STEPS_oracle" = sd(STEPS_oracle)) %>% 
    ungroup() %>% 
    group_by(ID) %>% 
    summarise("averageSDperID_observed" = mean(sd_STEPS), 
              "averageSDperID_oracle" = mean(sd_STEPS_oracle)) %>% 
    ungroup()
  
  conc_tab <- data.frame(
    "concB.STEPS.observed" = idComp.STEPS.observed, 
    "concB.STEPS.oracle" = idComp.STEPS.oracle, 
    "concW.STEPS.observed" = Wpairs_conc(ordered_sets[["setW"]], 
                                         ordered_sets[["setR"]], 
                                         Fx = STEPS_output$STEPS), 
    "concW.STEPS.oracle" = Wpairs_conc(ordered_sets_oracle[["setW"]],
                                       ordered_sets_oracle[["setR"]],
                                       Fx = STEPS_output$STEPS_oracle), 
    "concR.STEPS.obserevd" = mean(STEPS.Rsd$averageSDperID_observed),
    "concR.STEPS.oracle"   = mean(STEPS.Rsd$averageSDperID_oracle)
    
    
  )
  
  return(conc_tab)
}




sim_fun <- function(sim_dat,
                    obs_cov,
                    n,
                    p,
                    Rho,
                    noise.level,
                    noise.ratio,
                    lambdaB,
                    lambdaW,
                    lambdaR, 
                    results.tab,
                    iter.num, 
                    ntrees
) {
  
  resultsMat.iter <- data.frame("Iteration" = sort(rep(c(1:iter.num), nrow(results.tab))), results.tab, row.names = NULL)
  output_list <- list(list())
  for(k in 1:nrow(resultsMat.iter)) {
    ### Create Px explanatory variables
    obs.cov <- obs_cov
    obs.cov$score1 <- obs.cov$score1 + rnorm(n, mean = 0, sd = noise.level)
    obs.cov$score2 <- obs.cov$score2 + rnorm(n, mean = 0, sd = noise.level)
    noised_sim_dat <- addNoise(data = sim_dat, p = p, noise = noise.level / noise.ratio)


    trainTest <- divideTestTrain(percTrain = 0.8, data = noised_sim_dat)
    sim_dat_train <- trainTest[["Train"]]
    sim_dat_test  <- trainTest[["Test" ]]

    obs_cov_train <- obs.cov %>% filter(ID %in% trainTest[["train_index"]])
    obs_cov_test  <- obs.cov %>% filter(ID %in% trainTest[["test_index"]])
    #### Replace to extract pairs that are relevant
    ordered_sets_train <- SIM_sets_fun(data = sim_dat_train, 
                                       obs_cov = obs_cov_train,
                                       # inds = trainTest[["train_index"]], 
                                       oracle = F)
    ordered_sets_train_oracle <- SIM_sets_fun(data = sim_dat_train, 
                                              obs_cov = obs_cov_train,
                                              # inds = trainTest[["train_index"]], 
                                              oracle = T)
    ordered_sets_test  <- SIM_sets_fun(data = sim_dat_test, 
                                       obs_cov = obs_cov_test,
                                       # inds = trainTest[["test_index"]], 
                                       oracle = F)
    ordered_sets_test_oracle <- SIM_sets_fun(data = sim_dat_test, 
                                             obs_cov = obs_cov_test,
                                             # inds = trainTest[["test_index"]], 
                                             oracle = T)

    # Define x variables to use in model
    x_vars <- c(paste0("x", c(1:p)))
    # Extract explanatory variables
    train.x <- sim_dat_train %>%
      dplyr::select(one_of(x_vars))
    test.x <- sim_dat_test %>%
      dplyr::select(one_of(x_vars)) 
    # Number of of pairs of participants compared in set O
    n.Bpairs.train <- nrow(ordered_sets_train[["setB"]] %>%
                             group_by(ID.better, ID.severe) %>%
                             summarise( n()))

    n.Bpairs.train_oracle <- nrow(ordered_sets_train_oracle[["setB"]] %>%
                             group_by(ID.better, ID.severe) %>%
                             summarise( n()))
   
    # test_data <- list("data_test" = test.x %>% dplyr::select(-score1, -score2))
    test_data <- list("data_test" = test.x)

    # STEPS oracle solution 
    STEPS.output       <- STEPS(data = train.x,
                                setB.ord.pairs = ordered_sets_train[["setB"]],
                                setW.ord.pairs = ordered_sets_train[["setW"]],
                                setR.ord.pairs = ordered_sets_train[["setR"]],
                                n.B = n.Bpairs.train,
                                lambdaB = lambdaB,
                                lambdaW = lambdaW,
                                lambdaR = lambdaR,
                                rpart_control = NULL,
                                rho = Rho,
                                shrinkage = 1,
                                stop_terms = list("iter_count" = ntrees,
                                                  "minimal_weight" = 10^(-4)),
                                test_data = test_data,
                                plot.weights = FALSE,
                                plot.iteration.trees = FALSE,
                                y_plot = NULL,
                                y_test_plot = NULL,
                                value_initialize = 0, 
                                h = 0.5
    )
    
    STEPS.y.hat       <- STEPS.output$test$test_predictions
    STEPS.y.hat_train <- STEPS.output$predictions
    
    
    #### STEPS with pairs constructed according to real severity
    STEPS.oracle.output       <- STEPS(data = train.x,
                                       # data = train.x %>% dplyr::select(-score1, -score2, -R, -difficulty1, -difficulty2),
                                       setB.ord.pairs = ordered_sets_train_oracle[["setB"]],
                                       setW.ord.pairs = ordered_sets_train_oracle[["setW"]],
                                       setR.ord.pairs = ordered_sets_train_oracle[["setR"]],
                                       n.B = n.Bpairs.train_oracle,
                                       lambdaB = lambdaB,
                                       lambdaW = lambdaW,
                                       lambdaR = lambdaR,
                                       rpart_control = NULL,
                                       rho = Rho,
                                       shrinkage = 1,
                                       stop_terms = list("iter_count" = ntrees,
                                                         "minimal_weight" = 10^(-4)),
                                       test_data = test_data,
                                       plot.weights = FALSE,
                                       plot.iteration.trees = FALSE,
                                       y_plot = NULL,
                                       y_test_plot = NULL,
                                       value_initialize = 0, 
                                       h = 0.5
    )
    STEPS.oracle.y.hat       <- STEPS.oracle.output$test$test_predictions
    STEPS.oracle.y.hat_train <- STEPS.oracle.output$predictions
    
    results_list_test <-   aggResults(obs_cov = obs_cov_test, 
                                    test_data = data.frame(sim_dat_test, 
                                                           "STEPS" = STEPS.y.hat, 
                                                           "STEPS_oracle" = STEPS.oracle.y.hat),
                                    train_data = data.frame(sim_dat_train, 
                                                            "STEPS" = STEPS.y.hat_train, 
                                                            "STEPS_oracle" = STEPS.oracle.y.hat_train), 
                                    ordered_sets = ordered_sets_test, 
                                    ordered_sets_oracle = ordered_sets_test_oracle)

    output_list[["cors_tab"]][[k]] <- data.frame(resultsMat.iter[k, ],results_list_test[["resCor"]])
    output_list[["concordance_tab"]][[k]] <- data.frame(resultsMat.iter[k, ],results_list_test[["resConc"]])
    output_list[["classification_tab"]][[k]] <- data.frame(resultsMat.iter[k, ],results_list_test[["resClass"]])
    
     
  }
  output_list_out <- list()
  output_list_out[["cors_tab"]]           <- do.call('rbind', output_list[["cors_tab"]])
  output_list_out[["concordance_tab"]]    <- do.call('rbind', output_list[["concordance_tab"]])
  output_list_out[["classification_tab"]] <- do.call('rbind', output_list[["classification_tab"]])
  return(output_list_out)
  
}





### Save RDS file 
save_RDS_file_wrap_fun <- function(arg, dir, file_name) {
  ## Input: 1. arg = R argument to save
  ##        2. dir = path to directory
  ##        3. file_name = name for the saved file without the file type extension
  full_file_name <-  paste(dir, file_name, ".rds", sep = "")
  while(!file.exists(full_file_name)){
    tryCatch(saveRDS(arg, full_file_name), 
             error = function(e) Sys.sleep(5)) 
  }
  return()
}






### Save CSV file 
write_CSV_file_wrap_fun <- function(arg, dir, file_name, row_names = F) {
  ## Input: 1. arg = R argument to save
  ##        2. dir = path to directory
  ##        3. file_name = name for the saved file without the file type extension
  full_file_name <-  paste(dir, file_name, 
                           # Sys.Date(),
                           ".csv", sep = "")
  while(!file.exists(full_file_name)){
    tryCatch(expr =  write.csv(arg, full_file_name, 
                               row.names = row_names), 
             error = function(e) Sys.sleep(5)) 
  }
  return()
}



comb <- function(x, ...) {
  
  cors_results <- lapply(x, function(l) l[["cors_tab"]])
  conc_results <- lapply(x, function(l) l[["concordance_tab"]])
  class_results <- lapply(x, function(l) l[["classification_tab"]])
  
  cors.results  <- do.call('rbind', cors_results) 
  conc.results  <- do.call('rbind', conc_results)
  class.results <- do.call('rbind', class_results)
  
  
  mapply(rbind, x, ..., SIMPLIFY = F)
}

cat(paste("#", capture.output(sessionInfo()), "\n", collapse =""))