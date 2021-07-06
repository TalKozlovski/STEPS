### Create basic dataset
## Number of subjects in each group

if(!require('magrittr'))      install.packages('magrittr')
if(!require('dplyr'))         install.packages('dplyr')
if(!require('heatmaply'))     install.packages('heatmaply')

generateIDsWithSeverity <- function(n.subjects) {
   ID_tab <- data.frame("ID" = c(1:n.subjects)  ## obs index
                       )  %>%
    mutate("severity" =  c(1:n.subjects)) %>%
    mutate("severity01" = percentize(severity))

  return(ID_tab)
}


generateDifficultyLevelPerID <- function(n.subjects) {
  difficulties_per_ID <- expand.grid("ID" = c(1:n.subjects),
                            "difficulty1" = c(0, 1, 2),
                            "R" = c(1, 2))
  return(difficulties_per_ID)
}


### Logistic function
logistic <- function(x, L = 1, k = 1, a = 0) {
  return(L / (1 + exp(-k * (x - a))))
}


#### Create X matrix from stress test
x_sim_fun <- function(ID_tab, levels.tab, pX) {
  tab <- levels.tab %>%
    left_join(ID_tab)
  # tab <- tab %>%
  #   rowwise() %>%
  #   mutate("obs_severity" = case_when(sum(difficulty1, difficulty2) == 0 ~ logistic(severity01, k = 8, a = 1),
  #                                     sum(difficulty1, difficulty2) == 1 ~ logistic(severity01, k = 8, a = 0.8),
  #                                     sum(difficulty1, difficulty2) == 2 ~ logistic(severity01, k = 8, a = 0.6),
  #                                     sum(difficulty1, difficulty2) == 3 ~ logistic(severity01, k = 8, a = 0.4),
  #                                     sum(difficulty1, difficulty2) == 4 ~ logistic(severity01, k = 8, a = 0.2)))
  tab <- tab %>%
    rowwise() %>%
    mutate("obs_severity" = case_when(sum(difficulty1) == 0 ~ logistic(severity01, k = 8, a = 1),
                                      sum(difficulty1) == 1 ~ logistic(severity01, k = 8, a = 0.6),
                                      sum(difficulty1) == 2 ~ logistic(severity01, k = 8, a = 0.2)))

  tab[, paste("x", c(1), sep = "")] <- tab$obs_severity
  tab[, paste("x", c(2:pX), sep = "")] <- 0

  return(tab)
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





SIM_sets_fun <- function(data,
                         obs_cov
) {

  IDtimes.tab <- data[, c("ID", "R", "difficulty1")]
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
                                           setR.ord.pairs = setR.ord.pairs)

  return(list("setB" = setB.ord.pairs,
              "setR" = setR.ord.pairs,
              "setW" = setW.ord.pairs))

}


ord_pairs_fun <- function(pair, sort.by1, thresh = 1) {
  ## ord_pairs_fun sorts pairs of observations to the better (left) and severe (right) states
  # Input: @ pair - a pair of two observations indexes
  #        @ sort.by  A vector with values to sort according to the given pair
  # Output: vector of length 2, sorted from better to severe disease states.

  if((sort.by1[pair[1]] < (sort.by1[pair[2]] - thresh)) ) {
    return(pair)
  }
  if(((sort.by1[pair[1]] - thresh)  > sort.by1[pair[2]]) ) {
    return(pair[c(2:1)])
  } else {
    return(NULL)
  }

}


orderIDfun <- function(data) {
  n.id <- unique(data$ID)
  pairs <- t(combn(x = c(1:length(n.id)), m = 2))
  if(nrow(pairs) == 1) {
    ord.pairs <- list(ord_pairs_fun(pairs,    sort.by1 = data$score1))
  } else {
    ord.pairs <- apply(pairs, 1, ord_pairs_fun,
                       sort.by1 = data$score1)
  }


  ord.pairs <- do.call('rbind', ord.pairs)

  ord.pairs <- matrix(n.id[ord.pairs], ncol = 2)
  return(ord.pairs)
}

SIM_ordering_setB_fun <- function(obs_cov,
                                  setR.ord.pairs) {


  # Ordering pairs according to observed severity, score 1 & score 2

  orderID <- orderIDfun(data = obs_cov)

  tmp.setR.ord.pairs <- setR.ord.pairs %>%
    dplyr::select(ID, difficulty1,  R.group, n.R) %>%
    unique()

  orderID <- as.data.frame(orderID)
  colnames(orderID) <- c("ID", "ID.severe")

  full_tab1 <- orderID %>%
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



SIM_ordering_setW_fun <- function(IDtimes.tab,
                                  data,
                                  setR) {
  tmp_data.W <- setR %>%
    dplyr::select(ID, difficulty1,
                  R.group, n.R)
  duped <- duplicated(tmp_data.W)
  tmp_data.W <- tmp_data.W[!duped, ]

  ### X is the better state and Y is the severe
  tmp_data.W1 <- tmp_data.W %>%
    full_join(tmp_data.W, by = c("ID")) %>%
    filter((difficulty1.y > difficulty1.x)) %>%
    mutate(ID.y = ID)

  colnames(tmp_data.W1) <- c("ID.better",
                             "difficulty1.better",
                             "R.group.better",
                             "n.R.better",
                             "difficulty1.severe",
                             "R.group.severe",
                             "n.R.severe",
                             "ID.severe")

  ordered_tab.W <- tmp_data.W1 %>%
    group_by(ID.better, ID.severe) %>%
    mutate("n.Wi" = n()) %>%
    ungroup()

  return(ordered_tab.W)
}





SIM_ordering_setR_fun <- function(IDtimes.tab,
                                  data) {


  # Create R groups names
  ordered_tab.R <- IDtimes.tab
  # add index variable
  ordered_tab.R[, "index"] <- c(1:nrow(IDtimes.tab))
  ordered_tab.R[, "R.group"] <- as.character(apply(ordered_tab.R[, c("ID", "difficulty1")],
                                                   1, paste, collapse = "_"))
  # aggregate R groups to summarise table
  ordered_tab.R <- ordered_tab.R %>%
    group_by(R.group) %>%
    mutate("n.R" = n()) %>%
    ungroup()


  return(ordered_tab.R)

}

STEPS_configuration <- function(n_subjects) {
  IDs_tab <- generateIDsWithSeverity(n.subjects = n_subjects)
  difficulties_per_ID <- generateDifficultyLevelPerID(n.subjects = n_subjects )
  obs_cov <-  data.frame(IDs_tab,
                         "score1" = IDs_tab$severity)
  return(list("IDs_tab" = IDs_tab, "difficulties_per_ID" = difficulties_per_ID, "obs_cov" = obs_cov))
}
