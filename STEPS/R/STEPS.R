#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom data.table ":=" "%like%" "%between%"
#' @export
magrittr::`%>%`


#' Sorts pairs of observation to better (left) and severe (right) states
#'
#' @param pair a vector of length 2 with two observations indices.
#' @param sort.by A numeric or integer vector at the same length of number of observations to use for ordering. Higher value indiated severe state.
#' @return The function returns the \code{pair} ordered from better (left) to severe (right) state according to \code{sort.by}
#' @details The function checks which of the obsevations given in \code{pair} have the higher
#' @examples
#' set.seed(263)
#` x <- runif(20, -2, 2)
#` y <- 2 + 0.5*x + x^2
#' pairs <- t(combn(x = c(1:20), m = 2))
#' head(pairs)
#' # STEPS::ord_pairs_fun(pairs[1, ], sort.by = y)
#'
.ord_pairs_fun <- function(pair, sort.by) {
  # Check requirement to the sort.by argument to be numeric or integer
  if(!(class(sort.by)  %in% c("numeric", "integer"))) {
    stop("sort.by parametr need to be numeric or integer")
  }
  # If the sort.by for index 1 is lower than that of index 2 the pair is returned
  if(sort.by[pair[1]] < sort.by[pair[2]]) {
    return(pair)
  }
  else { # The pair is returned in opposite order
    return(pair[c(2:1)])
  }
}

#' Huber loss
#'
#' @param a a number to calculate huber loss of
#' @param h a positive threshold for fine. Default is 0.5.
#' @return The function returns the huber loss value of given a/
#' @details This function is slightly different from the regular huber loss, as all values smaller than -h
#' @examples
#' x <- runif(1000, -5, 5)
#' (h_loss <- huber_loss(x, h = 0.5))
#' plot(h_loss ~ x)
#'@export
huber_loss <- function(a, h = 0.5) {
  out <- rep(0,  length(a))
  ind_2 <- which(abs(a) < h)
  ind_3 <- which(a >= h)
  out[ind_2] <-  (a[ind_2] + h)^2 / (4 * h)
  out[ind_3] <-  a[ind_3]
  return(out)

}

#' Initial ranking function
#'
#' @param x a data.frame or a matrix that holds one row per observation.
#'          each row represents an index of an observation to look for in the \code{pairs.tab}
#' @param ord.pairs a table with sorted pairs of observations.
#' @return The function counts for each observation how many times it was smaller than
#'          other observations and how many times it was greater.
#' @details As a good starting point for y values we will use the information we have on our sorted pairs.
#'          An observations will get a higher rank for starter if they appear a lot of times in a
#'         severe state comapring to others observations.
#' @examples
#' set.seed(6324)
#` \dontrun{x <- matrix(sort(runif(10, -2, 2)), ncol = 1)}
#` \dontrun{rank_init_fun(x, ord.pairs = matrix(c(1:20), ncol = 2))}
#' @export
rank_init_fun <- function(x, ord.pairs) {
  n <- nrow(x)
  scores <- as.vector(tapply(c(1:n), c(1:n), count_rank_fun, ord.pairs = ord.pairs))
  return(c(scores))
}

#' Initial ranking function
#'
#' @param ind observation index to look for in \code{ord.pairs}
#' @param ord.pairs a table with sorted pairs of observations.
#' @return The function counts for a specific observation how many times it was smaller than
#'          other observations and how many times it was greater.
#' @details As a good starting point for y values we will use the information we have on our sorted pairs.
#'          An observations will get a higher rank for starter if they appear a lot of times in a
#'         severe state comapring to others observations.
#' @examples
#' \dontrun{
#' set.seed(6324)
#` x <- matrix(sort(runif(10, -2, 2)), ncol = 1)
#` y <- 2*x + 0.5 + x^2
#' pairs <- t(combn(x = c(1:10), m = 2))
#' ord.pairs <- t(apply(pairs, 1, ord_pairs_fun, sort.by = y))
#' count_rank_fun(ind = 1, ord.pairs)
#' }
count_rank_fun <- function(ind, ord.pairs) {
  minus <- sum(ind == ord.pairs[, 1])
  plus <- sum(ind == ord.pairs[, 2])
  score <- plus - minus
  return( score)
}

# Linear modeling ---------------------------------------------------------



#' Linear Model Objective Function
#'
#' @param pars Linear model equation's parameters.
#' @param y.hat.fold a table with two columns, and a row per pair we know its ordering. Each row includes the predictions for each pair.
#' @param h sapce parameter
#' @return The function returns the value of the linear DSS model objective function.
#' @details The objective function is of the shape:
#' \deqn{\frac{1}{2}\norm{w} + \frac{1}{2}*\sum_{pairs\in B}L_{h}\left( 1 - \hat{g}(x_{i}) - \hat{g}(x_{j}) \right)}
#' @examples
#' @export
linear_obj_term <- function(pars, y.hat.fold, h = 0.5) {

  # out  <- 0.5 * (t(eq) %*% eq) + 0.5 * (sum(sapply(1 - (y.hat.fold[ , 2]  - y.hat.fold[ , 1]), huber_lossC, h = 0.5)))
  # out  <- 0.5 * (t(pars) %*% pars) +
  #   0.5 * (sum(sapply(1 - (y.hat.fold[ , 2]  - y.hat.fold[ , 1]), huber_loss, h = 0.5)))
  out  <- 0.5 * (t(pars) %*% pars) +
    0.5 * (sum(huber_loss(1 - (y.hat.fold[ , 2] - y.hat.fold[ , 1]), h = h)))
  return(drop(out))
}


#' Non-linear Model Objective Function
#' @param Fx y
#' @param setB set B
#' @param setW set W
#' @param setR set R
#' @param lambdaB lambdaB
#' @param lambdaW lambdaW
#' @param lambdaR lambdaR
#' @param h space parameter
#' @return The function returns the value of the linear DSS model objective function.
#' @details The objective function is of the shape:
#' \deqn{\frac{1}{2}\norm{w} + \frac{1}{2}*\sum_{pairs\in O}L_{h}\left( 1 - \hat{g}(x_{i}) - \hat{g}(x_{j}) \right)}
#' @examples print("write an example")
#' @importFrom stats optimise predict sd setNames var
#'
#' @export
non_linear_loss_fun <- function(Fx, setB, setW, setR,
                                lambdaB, lambdaW, lambdaR,
                                h) {
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

  b.term <- tmp.B %>%
    # group_by(ID.severe, ID.better, Group.better, Group.severe) %>%
    mutate("loss_ijd" = huber_loss(1 - (Fx.severe.mean - Fx.better.mean), h = h)) %>%
    group_by(ID.better, ID.severe) %>%
    summarise("mean_ij" = mean(loss_ijd))

  b.term <- lambdaB * mean(b.term$mean_ij)


  tmp.W <- setW %>%
    left_join(data.frame("R.group.better" = names(Fx.mean),
                         "Fx.better.mean" = Fx.mean,
                         stringsAsFactors = F)) %>%
    left_join(data.frame("R.group.severe" = names(Fx.mean),
                         "Fx.severe.mean" = Fx.mean,
                         stringsAsFactors = F))

  w.term <- tmp.W %>%
    mutate("loss_id1d2" = huber_loss(1 - (Fx.severe.mean - Fx.better.mean), h = h)) %>%
    group_by(ID.better, ID.severe) %>%
    summarise("mean_loss_i" = mean(loss_id1d2))
  w.term <- lambdaW * mean(w.term$mean_loss_i)

  Fx.var <- tapply(X = Fx[setR$index],
                   INDEX = setR$R.group,
                   FUN = var_error_term)
  tmp.R <-  setR %>%
    left_join(data.frame("R.group" = names(Fx.var),
                         "Fx.var.R.group" = Fx.var,
                         stringsAsFactors = F))

  mean.R <- mean(tmp.R$Fx.var.R.group, na.rm = T)
  mean.R[is.nan(mean.R)] <- 0
  tmp.R[is.na(tmp.R$Fx.var.R.group), "Fx.var.R.group"] <- mean.R

  r.term <- tmp.R %>%
    group_by(ID) %>%
    summarise("mean_loss_i" = mean(Fx.var.R.group))
  r.term <- lambdaR * mean(r.term$mean_loss_i)

  loss.out <- b.term + w.term + r.term
  return(loss.out)
}

#' Variance error term
#'
#' @param subject.Fx subject.Fx
#' @return The function returns variance error term
#' @details The
#' @examples var_error_term(c(1, 2, 3))
#' @export
var_error_term <- function(subject.Fx) {
  # subject.Fx <- Fx
  if(length(subject.Fx) < 2) {
    return(NA)
  } else {
    return(var(subject.Fx))
  }
}


#' Derivatibe of Loss Objective Function according to F(x)
#'
#' @param Fx Fx
#' @param n.B n.B
#' @param setB setB
#' @param setW setW
#' @param setR setR
#' @param lambdaB lambdaB
#' @param lambdaW lambdaW
#' @param lambdaR lambdaR
#' @param h h
#' @return The function returns the value of the linear DSS model objective function.
#' @details The objective function is of the shape:
#' @examples
#' @export
deriv_fun <- function(Fx,
                      n.B,
                      setB, setW, setR,
                      lambdaB, lambdaW, lambdaR,
                      h) {

  n.id <- length(unique(setR$ID))
  Fx.mean <- tapply(X = Fx[setR$index],
                    INDEX = setR$R.group,
                    FUN = mean)


  # tic("o1")
  tmp.B <- setB %>%
    left_join(data.frame("Group.better" = names(Fx.mean),
                         "Fx.better.mean" = Fx.mean,
                         stringsAsFactors = F)) %>%
    left_join(data.frame("Group.severe" = names(Fx.mean),
                         "Fx.severe.mean" = Fx.mean,
                         stringsAsFactors = F))

  tmp.B <- data.table(tmp.B)
  # setDT(tmp.B)
  tmp.B[, "grad_loss_ijd" := .((1 / dij) * grad_huber_loss(1 - (Fx.severe.mean - Fx.better.mean), h = h))]

  b.term <- setR
  b.term <- data.table(b.term)
  # setDT(b.term)
  means.better <- tmp.B[, sum(grad_loss_ijd * ( 1 / n.R.better)), by = Group.better][, setNames(V1, Group.better)]
  means.severe <- tmp.B[, sum(grad_loss_ijd * ( -1 / n.R.severe)), by = Group.severe][, setNames(V1, Group.severe)]
  b.term[, "b.term" := .( sum(0, (lambdaB / n.B) * ( means.severe[as.character(R.group)] +
                                                       means.better[as.character(R.group)]), na.rm = T)),
         by = index]
  tmp.W <- setDT(setW)

  tmp.W <- setW %>%
    left_join(data.frame("R.group.better" = names(Fx.mean),
                         "Fx.better.mean" = Fx.mean,
                         stringsAsFactors = F)) %>%
    left_join(data.frame("R.group.severe" = names(Fx.mean),
                         "Fx.severe.mean" = Fx.mean,
                         stringsAsFactors = F))

  setDT(tmp.W)

  tmp.W[, "grad_loss_id1d2" := .((1 / n.Wi) * grad_huber_loss(1 - (Fx.severe.mean - Fx.better.mean), h = h))]

  w.term <- setR
  setDT(w.term)
  Wmeans.better <- tmp.W[, sum(grad_loss_id1d2 * ( 1 / n.R.better)), by = R.group.better][, setNames(V1, R.group.better)]
  Wmeans.severe <- tmp.W[, sum(grad_loss_id1d2 * ( -1 / n.R.severe)), by = R.group.severe][, setNames(V1, R.group.severe)]
  w.term[, "w.term" := .( sum(0, (lambdaW / n.id) * ( Wmeans.severe[as.character(R.group)] +
                                                 Wmeans.better[as.character(R.group)]), na.rm = T)),
         by = index]
  # toc()

  r.term <- setR %>%
    left_join(data.frame("R.group" = names(Fx.mean),
                         "Fx.mean" = Fx.mean,
                         stringsAsFactors = F)) %>%
    left_join(setR %>% group_by(ID) %>% summarise("di" = length(unique(R.group)))) %>%
    mutate("Fx" = Fx[index]) %>%
    mutate("r.term" = (lambdaR / n.id) * (1 / di) * (1 / max(1, (n.R - 1 ))) * 2 * (Fx - Fx.mean))

  out.deriv <- r.term[, c("index", "r.term")] %>%
    left_join(b.term[, c("index", "b.term")]) %>%
    left_join( w.term[, c("index", "w.term")] ) %>%
    arrange(index) %>%
    mutate("deriv" = b.term + r.term + w.term )


  return(out.deriv$deriv)


}


#' Huber loss gradient
#'
#' @param a vectore of gradients
#' @param h space parameter
#' @return The huber loss gradient
#' @examples grad_huber_loss(c(-0.5, 0.2, 1, 2))
#' @export
grad_huber_loss <- function(a, h = 0.5) {
  out <- rep(0,  length(a))
  ind_2 <- which(abs(a) < h)
  ind_3 <- which(a > h)
  out[ind_2] <- (a[ind_2] + h)/ (2 * h)
  out[ind_3] <-  1
  return(out)
}


#' Stopping terms validation
#'
#' @param iter_count iter_count
#' @param stop_terms stop_terms
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' @export
stop_term_fun <- function(iter_count, stop_terms) {
  if(iter_count >= stop_terms[["iter_count"]]) {
    return(T)
  }
  return(F)

}

#' Objective function per tree iteration
#'
#' @param par par
#' @param Fx_iter Fx_iter
#' @param tree tree
#' @param tree_data tree_data
#' @param setB setB
#' @param setW setW
#' @param setR setR
#' @param lambdaB lambdaB
#' @param lambdaW lambdaW
#' @param lambdaR lambdaR
#' @param h h = 0.5
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' @export
obj_term_per_step <- function(par,
                              Fx_iter,
                              tree,
                              tree_data,
                              setB,
                              setW,
                              setR,
                              lambdaB, lambdaW, lambdaR,
                              h = 0.5) {
  tmp_pred <- Fx_iter + par * predict(tree, tree_data, type = "vector")
  val_error <- non_linear_loss_fun(Fx = tmp_pred,
                                   setB,
                                   setW,
                                   setR,
                                   lambdaB, lambdaW, lambdaR,
                                   h = h)
  return(val_error)
}

#' Non-linear DSS training function
#'
#' @param data data
#' @param setB.ord.pairs A table with 2 columns that holds the pairs ordering information.
#' On the first column are the observations of a better state, and on the second column are the worsen state.
#' The table has the observations indices in respect to the set O data.
#' @param setW.ord.pairs Set W
#' @param setR.ord.pairs Set R
#' @param n.B number
#' @param lambdaB lambs B
#' @param lambdaW lambda W
#' @param lambdaR lambda R
#' @param rpart_control rapart control
#' @param rho rho
#' @param shrinkage shrinkage
#' @param stop_terms stop terms
#' @param test_data test
#' @param plot.weights binary
#' @param plot.iteration.trees binary
#' @param y_plot vec
#' @param y_test_plot vec test
#' @param value_initialize default 0
#' @param h space parameter
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' @export
STEPS <- function(data,
                  setB.ord.pairs,
                  setW.ord.pairs = NULL,
                  setR.ord.pairs = NULL,
                  n.B,
                  lambdaB = 1,
                  lambdaW = 1,
                  lambdaR = 1,
                  rpart_control = NULL,
                  rho = 1,
                  shrinkage = 0.1,
                  stop_terms = list("iter_count" = 5,
                                    "minimal_weight" = 10^(-6)),
                  test_data = NULL,
                  plot.weights = FALSE,
                  plot.iteration.trees = FALSE,
                  y_plot,
                  y_test_plot,
                  value_initialize = 0,
                  h = 0.5
)
{


  if(!is.data.frame(data)) {
    # predict.rpart requires data.frame as input
    warning("data was coerced into data.frame")
    data <- as.data.frame(data, stringsAsFactors = F)
  }
  if(is.null(test_data)) {
    warning("Test data is the same as the training data")
    # If we so not have a test data, then the train dataset is used as test
    test_data[["data_test"]] <- data
  }

  if(!is.data.frame(test_data[["data_test"]])) {
    # predict.rpart requires data.frame as input
    warning("The test data set was coerced into data.frame")
    test_data[["data_test"]] <- as.data.frame(test_data[["data_test"]], stringsAsFactors = F)
  }

  ## Initialize arguments for NLDSS training
  # Initial y values - We start from the value 0
  y_init <- rep(value_initialize, nrow(data))
  # Initial y values are our first prediction
  Fx_iter <- y_init
  # Initial our predictions for the test data
  Fx_test_iter <- rep(value_initialize, nrow(test_data[["data_test"]]))
  # Calculate the error from first prediction to train and test
  error <- non_linear_loss_fun(Fx = Fx_iter,
                               setB = setB.ord.pairs,
                               setW = setW.ord.pairs,
                               setR = setR.ord.pairs,
                               lambdaB, lambdaW, lambdaR,
                               h = h)

  error_test <- 0

  if(is.null(rpart_control)) { # Define rpart control parameters if were not given by the suer
    rpart_control <- rpart.control(minsplit = 2,
                                   cp = 1e-8,
                                   maxcompete = 0,
                                   maxsurrogate = 0,
                                   usesurrogate = 0,
                                   xval = 0,
                                   surrogatestyle = 0,
                                   maxdepth = 5)
  }
  # Initial prediction table per iteration for the function to return
  Fx_out <- matrix(Fx_iter, ncol = 1)
  Fx_test_out <- matrix(Fx_test_iter, ncol = 1)
  # Iterations counting
  iter_count <- 0
  # Initial trees list and trees weights vector
  trees_list <- list(NA)
  trees_weights <- NULL
  stop_reason <- NA

  repeat{ # GBRT algorithm

    # Update iteration count
    iter_count <- iter_count + 1
    print(paste("Iteration number: ", iter_count))
    # Update target variable to be the negative gradient
    # tic("grad_fun")
    y_iter <- -1 * deriv_fun(Fx = Fx_iter,
                             n.B = n.B,
                             setB = setB.ord.pairs,
                             setW = setW.ord.pairs,
                             setR = setR.ord.pairs,
                             lambdaB, lambdaW, lambdaR,
                             h = h)

    # toc()
    # Merge target and explanatory variable to one data.table
    tree_data <- data.frame("y_iter" = c(y_iter), data, stringsAsFactors = F)
    # Build a shallow regression tree
    # tic("rpart tree")
    tree_iter <- rpart(y_iter ~. ,
                       data = tree_data,
                       method = "anova",
                       x = FALSE,
                       y = FALSE,
                       control = rpart_control)
    # toc()
    # Save the iteration's tree to output
    trees_list[[iter_count]] <- tree_iter

    ### find the step size should be taken (parallel to rho in GBRT algorithm)
    # tic("optimise")
    rho.opt <- optimise(f = obj_term_per_step,
                        interval = c(0, 1),
                        Fx_iter = Fx_iter,
                        tree = tree_iter,
                        tree_data = tree_data,
                        setB = setB.ord.pairs,
                        setW = setW.ord.pairs,
                        setR = setR.ord.pairs,
                        lambdaB,
                        lambdaW,
                        lambdaR,
                        h = h,
                        maximum = FALSE)
    rho <- rho.opt$minimum
    # toc()
    # Save step size for output
    trees_weights[iter_count + 1] <- rho
    # Update predictions
    Fx_iter <- Fx_iter +
      (shrinkage * rho * predict(tree_iter, data, type = "vector"))
    # Update test predictions
    # tic("predict")
    Fx_test_iter <- Fx_test_iter +
      (shrinkage * rho * predict(tree_iter, newdata = test_data[["data_test"]], type = "vector"))
    # toc()
    # if(plot.iteration.trees) {
    #   plot(Fx_iter ~ y_plot, main = paste("iter ", iter_count))
    #   # plot(Fx_test_iter ~ y_test_plot, main = paste("iter ", iter_count))
    # }

    # calculate the error of the updated pedictions
    # tic("error")
    error <- c(error, non_linear_loss_fun(Fx = Fx_iter,
                                          setB = setB.ord.pairs,
                                          setW = setW.ord.pairs,
                                          setR = setR.ord.pairs,
                                          lambdaB,
                                          lambdaW,
                                          lambdaR,
                                          h = h))
    # toc()
    # error_test <- c(error_test, non_linear_obj_term(Fx = Fx_test_iter,
    #                                                 subject.indices.tab = test_data[["subject.indices.tab_test"]],
    #                                                 setB.ord.pairs = test_data[["setB.ord.pairs_test"]],
    #                                                 lambdaR))
    # error_test <- 0
    # If the error value started to increase we break the loop
    if(abs(error[iter_count + 1] - error[iter_count]) < stop_terms[["minimal_weight"]]) {
      stop_reason <- "Index converged"
      break
    }
    # print(paste("Error value: ", round(error[iter_count + 1], 3)))
    # Save updated predictions for output
    Fx_out <- cbind(Fx_out, Fx_iter)
    Fx_test_out <- cbind(Fx_test_out, Fx_test_iter)
    if(iter_count >= stop_terms[["iter_count"]]) {
      stop_reason <- "Reached number of iterations"
      break
    }

  }
  # print(paste("stop reason: ", stop_reason))
  p <- ncol(Fx_out)
  colnames(Fx_out) <- paste("Fx", "Iteration", c(0:(p - 1)), sep = "_")
  colnames(Fx_test_out) <- paste("Fx", "test", "Iteration", c(0:(p - 1)), sep = "_")
  return(list("trees_info" = trees_list,
              "scale_parameters" = list("mean" = mean(Fx_out[, p]),
                                        "sd" = sd(Fx_out[, p])),
              "weights" = trees_weights,
              "prediction_error" = error,
              "predictions" = Fx_iter,
              "Fx_per_iteration" = Fx_out,
              "test" = list("test_predictions" = Fx_test_iter,
                            "Fx_test_per_iteration" = Fx_test_out,
                            "prediction_error_test" = error_test),
              "stop_reason" = stop_reason))


}



#' DSS prediction
#'
#' @param tree_model list of trees output. Must include arguments with the names "weights" an
#' "trees_info". even if it is just one model it should still be a list
#' @param newdata newdata set
#' @param shrinkage shrinakge rate, default is 0.1
#' @param rho if NULL than the weights per iteration are taken from tree_model$weights
#' @param summarise_fun if the models are from cross validation, or we have more than one,
#' then we can summarise the different predictions per observations by mean or emdian for example.
#' Default is set to mean.
#' @return A vector with predictions
#' @examples print("missing ecample")
#' @export
prediction.STEPS <- function(tree_model, newdata, shrinkage = 0.1, rho = NULL, summarise_fun = mean) {
  len_model <- length(tree_model)
  predictions <- matrix(0, ncol = len_model, nrow = nrow(newdata))
  scale_predictions <- matrix(0, ncol = len_model, nrow = nrow(newdata))
  if(is.null(rho)) {
    rhos <- lapply(tree_model, function(l) l$weights[-1])
  }
  trees_info <- lapply(tree_model, function(l) l$trees_info)
  scale_pars <- lapply(tree_model, function(l) l$scale_parameters)
  for(i in 1:len_model) {
    tmp.rhos <- rhos[[i]]
    tmp.trees <- trees_info[[i]]
    preds <- lapply(tmp.trees, predict, newdata = newdata)

    for(j in 1:length(preds)) {
      predictions[, i] <- predictions[ , i] + shrinkage * tmp.rhos[j] * preds[[j]]

    }
    scale_predictions[, i] <- (predictions[, i] - scale_pars[[i]][["mean"]]) / scale_pars[[i]][["sd"]]
  }
  final_pred <- apply(scale_predictions, 1, summarise_fun)
  return(final_pred)
}


