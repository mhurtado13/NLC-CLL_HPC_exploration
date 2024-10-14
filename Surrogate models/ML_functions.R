
##Functions from https://github.com/VeraPancaldiLab/CellTFusion

##Load libraries

libraries_set <- function(){
  suppressMessages(library("BiocManager"))
  suppressMessages(library("devtools"))
  suppressMessages(library("remotes"))
  suppressMessages(library("tidyr"))
  suppressMessages(library("dplyr"))
  suppressMessages(library("matrixStats"))
  suppressMessages(library("reshape2"))
  suppressMessages(library("purrr"))
  suppressMessages(library("tidygraph"))
  suppressMessages(library("stringr"))
  suppressMessages(library("tibble"))
  suppressMessages(library("gplots"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("AnnotationDbi"))
  suppressMessages(library("RColorBrewer"))
  suppressMessages(library("pheatmap"))
  suppressMessages(library("ggfortify"))
  suppressMessages(library("msigdbr"))
  suppressMessages(library("Hmisc"))
  suppressMessages(library("ggpubr"))
  suppressMessages(library("ggstatsplot"))
  suppressMessages(library("dendextend"))
  suppressMessages(library("stats"))
  suppressMessages(library("caret"))
  suppressMessages(library("pROC"))
  suppressMessages(library("MLeval"))
  suppressMessages(library("rms"))
  suppressMessages(library("igraph"))
  suppressMessages(library("uuid"))
  suppressMessages(library("parallel"))
  suppressMessages(library("factoextra"))
  suppressMessages(library("doParallel"))
  suppressMessages(library("foreach"))
  suppressMessages(library("kernlab")) ##For Gaussiann Process Regression
  suppressMessages(library("mgcv")) #For Radial Basis Function
  suppressMessages(library("MASS")) #For polynomial
  
}

libraries_set()

dir.create(file.path(getwd(), "Results"))

compute.boruta <- function(data, seed, fix = TRUE) {
  
  set.seed(seed)
  boruta_output <- Boruta(target ~ ., data = data, doTrace = 0)
    
  if (fix) {
    roughFixMod <- TentativeRoughFix(boruta_output)
    boruta_output <- roughFixMod
  }
    
  imps <- attStats(boruta_output)
  decision <- as.character(imps$decision)
  
  res <- imps %>%
    data.frame() %>%
    rownames_to_column("Variable") %>%
    dplyr::select(-decision)
  
  
  
  return(list(res, decision))
}

merge_boruta_results = function(importance_values, decisions, file_name, iterations, threshold, return = T){
  
  ### Construct matrix of importance
  combined_importance <- do.call(rbind, importance_values)
  combined_results_long <- combined_importance %>% #Matrix for plotting
    pivot_longer(cols = meanImp, names_to = "Measure", values_to = "Value")
  
  median_df <- combined_importance %>% #Calculate the median for each column, grouped by the variable name
    group_by(Variable) %>%
    dplyr::summarize(across(everything(), \(x) median(x, na.rm = TRUE)))
  
  ### Retrieve important and tentatives variables
  
  combined_results <- do.call(cbind, decisions)
  rownames(combined_results) = median_df$Variable
  decisions_summary <- apply(combined_results, 1, function(x) {
    table(factor(x, levels = c("Confirmed", "Tentative", "Rejected")))
  })
  confirmed_vars <- names(which(decisions_summary["Confirmed",] >= round(threshold*iterations))) 
  tentative_vars <- names(which(decisions_summary["Tentative",] >= round(threshold*iterations))) 
  
  # For plotting
  combined_results_long$Decision = "Rejected"
  combined_results_long$Decision[which(combined_results_long$Variable %in% confirmed_vars)] = "Confirmed"
  combined_results_long$Decision[which(combined_results_long$Variable %in% tentative_vars)] = "Tentative"
  
  mean_order <- median_df %>% #Extract the order of variables for plotting
    arrange(meanImp) %>%
    pull(Variable)
  
  # For result 
  median_df$Decision = "Rejected"
  median_df$Decision[which(median_df$Variable %in% confirmed_vars)] = "Confirmed"
  median_df$Decision[which(median_df$Variable %in% tentative_vars)] = "Tentative"
  
  # Plot variable importance boxplots
  if(return){
    pdf(paste0("Results/Boruta_variable_importance_", file_name, ".pdf"), width = 8, height = 12)
    print(ggplot(combined_results_long, aes(x = factor(Variable, levels = mean_order), y = Value, fill = Decision)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            coord_flip() +
            labs(x = "Features", y = "Importance", title = paste0("Variable Importance by Boruta after ", iterations, " bootstraps\n", file_name)) +
            scale_fill_manual(values = c("Confirmed" = "green", "Tentative" = "yellow", "Rejected" = "red")) +
            facet_wrap(~ Measure, scales = "free_y"))
    dev.off()
  }

  return(list(Confirmed = confirmed_vars, Tentative = tentative_vars, Matrix_Importance = median_df))
}

feature.selection.boruta <- function(data, iterations = NULL, fix, doParallel = F, workers=NULL, file_name = NULL, threshold = NULL, return) {
  if(doParallel){
    if(is.null(iterations) == T){
      stop("No iterations specified for running in parallel, please set a number. If you want to run feature selection once consider setting doParallel = F")
    }else{
      if(is.null(workers)==T){
        num_cores <- detectCores() - 1  
      }else{
        num_cores <- workers
      }
      
      cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
      doParallel::registerDoParallel(cl)
      
      message("Running ", iterations, " iterations of the Boruta algorithm using ", num_cores, " cores")
      # arg_list <- replicate(iterations, list(data, sample.int(100000, 1), fix), simplify = FALSE)
      # system.time({
      #   res <- mclapply(arg_list, function(x) {
      #     do.call(compute.boruta, x)
      #   }, mc.cores = num_cores)
      # }) 
      
      res <- foreach(seed = sample.int(100000, iterations)) %dopar% {
        
        source("src/environment_set.R") 
        
        tryCatch({
          # If successful, return the result and the seed
          list(result = compute.boruta(data, seed, fix), 
               error = NULL, 
               seed = seed)
        }, error = function(e) {
          # If an error occurs, return the error message and the seed for debugging
          list(result = NULL, error = e$message, seed = seed)
        })
      }
      
      parallel::stopCluster(cl)
      unregister_dopar() #Stop Dopar from running in the background
      
    }
    
    # Extract the first sublist of each element
    matrix_of_importance <- lapply(res, function(x) x[[1]])
    
    # Extract the second sublist of each element
    features_labels <- lapply(res, function(x) x[[2]])
    
    res = merge_boruta_results(matrix_of_importance, features_labels, file_name = file_name, iterations = iterations, threshold = threshold, return = return)
    
  }else{
    if(is.null(iterations) == T){
      stop("No iterations specified for running Boruta algorithm for feature selection")
    }else{
      message("Running ", iterations, " iterations of the Boruta algorithm")
      res = list()
      for (i in 1:iterations) {
        res[[i]] = compute.boruta(data, seed = sample.int(100000, 1), fix)
      }
      
      # Extract the first sublist of each element
      matrix_of_importance <- lapply(res, function(x) x[[1]])
      
      # Extract the second sublist of each element
      features_labels <- lapply(res, function(x) x[[2]])
      
      res = merge_boruta_results(matrix_of_importance, features_labels, file_name = file_name, iterations = iterations, threshold = threshold, return = return)
    }

  }
  
  return(res)
  
}

# Create folds for each repetition with different seeds
create_folds_for_repetitions <- function(data, k_folds, n_rep) {
  all_folds <- list()
  for (rep in 1:n_rep) {
    set.seed(sample.int(100000, 1)) # Change seed for each repetition
    folds <- createFolds(data$target, k = k_folds, returnTrain = TRUE, list = TRUE)
    all_folds[[rep]] <- folds
  }
  return(all_folds)
}

computed.boruta.kfolds = function(folds, data_model, boruta_iterations, fix_boruta, tentative, threshold, file_name){
  
  folds_threshold = 0.8*length(folds) 
  features_folds = list()
  
  for (i in 1:length(folds)) {
    message("Feature selection using Boruta...............................................................\n\n")
    training_set = data_model[folds[[i]],]
    res_boruta = feature.selection.boruta(training_set, iterations = boruta_iterations, fix = fix_boruta, thres = threshold, file_name = file_name, return = T)
    
    if(tentative == F){
      if(length(res_boruta$Confirmed) <= 1){
        features_folds[[i]] = list()
        message("\nNo features were confirmed in more than ", round(threshold*100) ,"% of the times for training in this specific fold.......................\n\n")
      }
      message("\nKeeping only features confirmed in more than", round(threshold*100) ,"% of the times for training in this specific fold......................\n\n")
      message("If you want to consider also tentative features, please specify tentative = T in the parameters.\n\n")
      features_folds[[i]] = res_boruta$Confirmed
    }else{
      sum_features = length(res_boruta$Confirmed) + length(res_boruta$Tentative)
      if(sum_features <= 1){
        features_folds[[i]] = list()
        message("\nNo features were confirmed in more than ", round(thresh*100) ,"% of the times for training in this specific fold.......................\n\n")
      }
      message("\nKeeping only features confirmed and tentative in more than", round(thresh*100) ,"% of the times for training in this specific fold............................\n\n")
      features_folds[[i]] = c(res_boruta$Confirmed, res_boruta$Tentative)
    }
  }
  
  all_features <- unlist(features_folds)
  feature_freq <- table(all_features)
  selected_features <- names(feature_freq[feature_freq >= folds_threshold])
  
  if(length(selected_features)<=1){
    message("No features selected meet the requirements. Try with different parameter values.")
  }else{
    return(selected_features)
  }
    
}

# polynomial_model = function(data, target){
#   model = lm(target ~ poly(., degree = 2), data = data)
#   return(model)
# }

polynomial_model <- list(
  type = "Regression",
  library = NULL,  # No need to load any external libraries
  loop = NULL,
  
  # Define the hyperparameter grid (polynomial degree)
  parameters = data.frame(
    parameter = 'degree', 
    class = 'numeric', 
    label = 'Polynomial Degree'
  ),
  
  # Grid function for tuning the polynomial degree
  grid = function(x, y, len = NULL, search = "grid") {
    data.frame(degree = seq(1, 5, by = 1))  # Tune degrees from 1 to 5
  },
  
  # Fitting function: applies polynomial regression
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    lm(y ~ poly(as.matrix(x), degree = param$degree))
  },
  
  # Predict function: predicts based on the model
  predict = function(modelFit, newdata, submodels = NULL) {
    predict(modelFit, as.matrix(newdata))
  },
  
  # We don't need probability predictions for regression
  prob = NULL
)

gaussian_process_model = function(data, target){
  model = gausspr(target ~ ., data = data)
  return(model)
}

rbf_model = function(data, target){
  model = gam(target ~ s(.), data = data)
  return(rbf_model)
}

compute_surrogate_models = function(train_data, k_folds, n_rep, file_name = NULL){
  cat("Training surrogate models...............................................................\n\n")
  
  ######### Machine Learning models
  metric <- "RMSE" #metric to use for selecting best methods (default: Accuracy -- for AUC see below and parameter must be equal to cv_metric = "AUC")
  
  ### Stratify K fold cross-validation 
  multifolds <- createMultiFolds(train_data[,'target'], k = k_folds, times = 1) #repeated folds
  trainControl <- trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=1, verboseIter = F, allowParallel = F, savePredictions=T)
  
  ################## Polynomial
  fit.polynomial <- train(target~., data = train_data, method = polynomial_model, metric = metric, trControl = trainControl) 
  predictions.polynomial <- data.frame(Polynomial = predict(polynomial_fit, newdata = train_data))
  polynomial_fit <- train(target ~ ., data = train_data, method = polynomial_model, metric = metric, trControl = trainControl, tuneLength = 2)
  ################## Gaussian
  fit.gaussian <- train(target~., data = train_data, method = gaussian_process_model, metric = metric,trControl = trainControl)
  predictions.gaussian = data.frame(Gaussian = predict(fit.gaussian, newdata = train_data))
  
  ################## Radial Basis Function RBF
  fit.rbf <- train(target~., data = train_data, method=rbf_model, metric=metric,trControl=trainControl)
  predictions.rbf = data.frame(RBF = predict(fit.rbf, newdata = train_data))
  
  ############################################################## Save models
  
  ensembleResults <- list(Polynomial = fit.polynomial,
                          Gaussian = fit.gaussian,
                          RBF = fit.rbf)
  
  
  model_predictions = list(Polynomial = predictions.polynomial,
                           Gaussian = predictions.gaussian,
                           RBF = predictions.rbf)
  
  
  rm(fit.polynomial, fit.gaussian, fit.rbf, multifolds)
  gc()
  
  model_predictions = do.call(cbind, model_predictions) #Join as data frame
  
  metrics = compute_cv_RMSE(ensembleResults, file_name)
  
  top_model = metrics[["Top_model"]]
  
  model = ensembleResults[[top_model]]
  
  cat("Best surrogate model found: ", top_model, "\n")
  
  cat("Returning model trained\n")
  
  output = list("Model" = model, "ML_Models" = ensembleResults)
  
  return(output)
  
}

compute.RMSE_k_fold_CV = function(train_data, k_folds, n_rep, stacking = F, file_name = NULL){
  
  cat("Training machine learning model...............................................................\n\n")
  
  ######### Machine Learning models
  metric <- "RMSE" #metric to use for selecting best methods (default: Accuracy -- for AUC see below and parameter must be equal to cv_metric = "AUC")
  
  ######### Stratify K fold cross-validation 
  #folds <- createFolds(train_data[,'target'], k = k_folds, returnTrain = T, list = T) #this for single folds
  multifolds <- createMultiFolds(train_data[,'target'], k = k_folds, times = n_rep) #repeated folds
  trainControl <- trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = F)
  
  ##################################################### ML models
  ################## Bagged CART
  fit.treebag <- train(target~., data = train_data, method = "treebag", metric = metric,trControl = trainControl) 
  predictions.bag <- data.frame(BAG = predict(fit.treebag, newdata = train_data))
  
  ################## RF
  fit.rf <- train(target~., data = train_data, method = "rf", metric = metric,trControl = trainControl)
  predictions.rf = data.frame(RF = predict(fit.rf, newdata = train_data))
  
  ################## LG - Logistic Regression
  fit.glm <- train(target~., data = train_data, method="glm", metric=metric,trControl=trainControl)
  predictions.glm = data.frame(GLM = predict(fit.glm, newdata = train_data))
  
  ################## GLMNET - Regularized Logistic Regression (Elastic net)
  fit.glmnet <- train(target~., data = train_data, method="glmnet", metric=metric,trControl=trainControl)
  predictions.glmnet = data.frame(GLMNET = predict(fit.glmnet, newdata = train_data))
  
  ################## KNN - k-Nearest Neighbors 
  fit.knn <- train(target~., data = train_data, method="knn", metric=metric,trControl=trainControl)
  predictions.knn = data.frame(KNN = predict(fit.knn, newdata = train_data))
  
  ################## CART - Classification and Regression Trees (CART), 
  fit.cart <- train(target~., data = train_data, method="rpart", metric=metric,trControl=trainControl)
  predictions.cart = data.frame(CART = predict(fit.cart, newdata = train_data))
  
  ################## Regularized Lasso
  fit.lasso <- train(target~., data = train_data, method="glmnet", metric=metric,trControl=trainControl, tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 20)))
  predictions.lasso = data.frame(LASSO = predict(fit.lasso, newdata = train_data))
  
  ################## Ridge regression
  fit.ridge <- train(target~., data = train_data, method="glmnet", metric=metric,trControl=trainControl, tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 1, length = 20)))
  predictions.ridge = data.frame(RIDGE = predict(fit.ridge, newdata = train_data))
  
  ################## Support Vector Machine with Radial Kernel
  fit.svm_radial <- train(target ~ ., data = train_data, method = "svmRadial", metric = metric, trControl = trainControl)
  predictions.svm_radial = data.frame(SVM_radial = predict(fit.svm_radial, newdata = train_data, type = "prob"))
  
  ################## Support Vector Machine with Linear Kernel
  fit.svm_linear <- train(target ~ ., data = train_data, method = "svmLinear", metric = metric, trControl = trainControl)
  predictions.svm_linear = data.frame(SVM_linear = predict(fit.svm_linear, newdata = train_data))
  
  ############################################################## Save models
  
  ensembleResults <- list(BAG = fit.treebag,
                          RF = fit.rf,
                          GLM = fit.glm,
                          KNN = fit.knn,
                          CART = fit.cart,
                          GLMNET = fit.glmnet,
                          LASSO = fit.lasso,
                          RIDGE = fit.ridge,
                          SVM_radial = fit.svm_radial,
                          SVM_linear = fit.svm_linear)
  
  
  model_predictions = list(BAG = predictions.bag,
                           RF = predictions.rf,
                           GLM = predictions.glm,
                           KNN = predictions.knn,
                           CART = predictions.cart,
                           GLMNET = predictions.glmnet,
                           LASSO = predictions.lasso,
                           RIDGE = predictions.ridge,
                           SVM_radial = predictions.svm_radial,
                           SVM_linear = predictions.svm_linear)
  

  rm(fit.treebag, fit.rf, fit.glm, fit.knn, fit.cart, fit.glmnet, fit.lasso, fit.ridge, fit.svm_radial, fit.svm_linear, multifolds)
  gc()
  
  model_predictions = do.call(cbind, model_predictions) #Join as data frame
  
  if(stacking){
    #Base models using ML models with best accuracy from each family
    base_models = compute_cv_RMSE(ensembleResults, base_models = T)
    
    features_predictions = model_predictions %>%
      t() %>%
      data.frame() %>%
      rownames_to_column("Models") %>%
      filter(grepl(paste0("\\b(", paste(base_models$Base_models, collapse = "|"), ")\\b"), Models)) %>%
      column_to_rownames("Models") %>%
      t() %>%
      data.frame()
    
    meta_features = cbind(features_predictions, "true_label" = train_data$target) 
    
    meta_learner <- train(true_label ~ ., data = meta_features, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression
    
    #Base models using ALL ML models 
    meta_features_all = cbind(model_predictions, "true_label" = train_data$target) 
    
    meta_learner_all <- train(true_label ~ ., data = meta_features_all, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression
    
    cat("Meta-learners ML model based on GLM\n")
    output = list("Meta_learners" = list("simple" = meta_learner, "all" = meta_learner_all), "Base_models" = base_models$Base_models, "ML_models" = ensembleResults)
    
  }else{
    metrics = compute_cv_RMSE(ensembleResults, file_name)
    
    top_model = metrics[["Top_model"]]
    
    model = ensembleResults[[top_model]]
    
    cat("Best ML model found: ", top_model, "\n")
    
    cat("Returning model trained\n")
    
    output = list("Model" = model, "ML_Models" = ensembleResults)
  }
  
  
  return(output)
  
}

compute.ML = function(samples, target_variable, partition, stack, seed, file_name = NULL){

  samples = samples %>%
    rename_with(~"target", all_of(target_variable)) #In case target variable is with a different name
  
  set.seed(seed)   
  
  # Do partition
  index = createDataPartition(samples[,"target"], times = 1, p = partition, list = FALSE) 
  
  # Train cohort
  training_set = samples[index,]
  
  # Test cohort
  testing_set = samples[-index,]
  
  ###############################################################################################################################################################################
  
  ####################### Training 5 kfolds and 100 repetitions 
  training = compute_surrogate_models(training_set, k_folds = 5, n_rep = 100, file_name = file_name)

  if(length(training)!=0){
    
    #Extract target variable
    target = testing_set %>%
      pull(target)
    #Remove target variable
    testing_set$target = NULL 
    
    ####################### Testing set
    if(stack){
      meta_learner = training[["Meta_learners"]]
      prediction = compute.prediction.stacked(meta_learner, testing_set, target, training[["ML_models"]], training[["Base_models"]])
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy from CV per partition
      prediction = compute.prediction(model, testing_set, target)
    }
    
    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]
    
    gc() #Clean garbage
    
    return(list(Predictions = predictions, Prediction_metrics = metrics))
  }else{  #No features are selected as predictive
    
    gc() #Clean garbage
    
    return(NULL)
  }
  
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


compute.bootstrap.ML = function(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition = 0.7, iterations, doFS = F, stack, workers = NULL, plots = F, file.name = NULL){
  
  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- detectCores() - 1
    }else{
      num_cores <- workers
    }
    
    #options(cluster.timeout = 300)
    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)
    
    #future::plan("multicore", workers = num_cores)
    #future::plan()
    
    message("Running ", iterations, " splits for training and test using ", num_cores, " cores")
    
    #List of arguments inputs
    #arg_list <- replicate(iterations, list(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition, sample.int(100000, 1)), simplify = FALSE)
    
    
    # Run in parallel with error handling
    # system.time({
    #   res <- mclapply(seq_along(arg_list), function(i) {
    #     args <- arg_list[[i]]
    #     result = tryCatch({
    #         list(result = do.call(compute.ML, args), error = NULL, data = args)
    #       }, error = function(e) {
    #         list(result = NULL, error = e$message, data = args)
    #       })
    #     rm(args) #Remove input arguments to free memory after each iteration
    #     invisible(gc()) #Do garbage collection
    #     return(result)
    #   }, mc.cores = num_cores)
    # })

    #options(future.globals.maxSize = 10 * 1024^3)  # Set limit to 10 GiB
    # res <- future.apply::future_lapply(sample.int(100000, iterations), function(seed) {
    #   compute.ML(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition, seed, universe)
    # }, future.seed = TRUE)
    
    # res <- future.apply::future_lapply(sample.int(100000, iterations), function(seed) {
    #   # Use tryCatch to handle any errors during the function call
    #   result <- tryCatch({
    #     # If successful, return the result and the seed
    #     list(result = compute.ML(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition, seed, universe), 
    #          error = NULL, 
    #          seed = seed)
    #   }, error = function(e) {
    #     # If an error occurs, return the error message and the seed for debugging
    #     list(result = NULL, error = e$message, seed = seed)
    #   })
    #   
    #   return(result)
    # }, future.seed = TRUE)
    
    res <- foreach(random.seed = sample.int(100000, iterations)) %dopar% {
      
      source("src/environment_set.R") 
      
      tryCatch({
        # If successful, return the result and the seed
        list(result = compute.ML(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition, stack = stack, feature.selection = doFS, seed = random.seed, return = plots, file_name = file.name), 
             error = NULL, 
             seed = random.seed)
      }, error = function(e) {
        # If an error occurs, return the error message and the seed for debugging
        list(result = NULL, error = e$message, seed = random.seed)
      })
    }
    
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background
    
  }
  
  # # Extract the first sublist of each element
  # matrix_of_importance <- lapply(res, function(x) x[[1]])
  # 
  # # Extract the second sublist of each element
  # features_labels <- lapply(res, function(x) x[[2]])
  # 
  # res = merge_boruta_results(matrix_of_importance, features_labels, file_name, iterations, threshold = 0.7)
  
  #future::plan(future::sequential)
  
  return(res)
}

compute_cv_RMSE = function(models, file_name = NULL, base_models = F, return = T){
  
  #Bind accuracy values from each model
  rmse = list()
  for (i in 1:length(models)){
    rmse[[i]] = models[[i]]$resample %>% 
      mutate(model = names(models)[i])
    names(rmse)[i] = names(models)[i]
  }
  rmse_data = do.call(rbind, rmse)
  
  #Retrieve top model based on accuracy
  res_rmse <- rmse_data %>%
    group_by(model) %>%
    summarise(RMSE = mean(RMSE)) %>%
    arrange(RMSE)  # Arrange in ascending order to get min RMSE
  
  top_model <- res_rmse %>%
    slice_head(n = 1) %>%
    pull(model)
  
  if(return){
    pdf(paste0("Results/RMSE_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot(rmse_data, aes(x = model, y = RMSE, fill = model)) +
           geom_boxplot() +
           labs(title = "Distribution of RMSE by Model",
                x = "Model",
                y = "RMSE") +
           theme_minimal() +
           theme(legend.position = "none"))
    dev.off()
  }
  
  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = "RMSE")
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("RMSE" = res_rmse, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("RMSE" = res_rmse, "Top_model" = top_model))
  }
  
}

compute_cv_AUC = function(models, file_name = NULL, base_models = F, return = T){
  
  #Bind AUC values from each model
  auc = list()
  for (i in 1:length(models)){
    auc[[i]] = models[[i]]$resample %>% #we use the resample matrix and not directly the results matrix as some have hyperparameters so we will need to define best on the tuned parameter (=more code) - resample matrix is made based on the best tuning
      mutate(model = names(models)[i])
    names(auc)[i] = names(models)[i]
  }
  auc_data = do.call(rbind, auc)

  #Retrieve top model based on accuracy
  res_auc <- auc_data %>%
    group_by(model) %>%
    summarise(AUC = mean(AUC))  %>%
    arrange(desc(AUC)) 
  
  top_model = res_auc %>%
    slice(1) %>%
    pull(model)
  
  if(return){
    pdf(paste0("Results/AUC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot(auc_data, aes(x = model, y = AUC, fill = model)) +
           geom_boxplot() +
           labs(title = "Distribution of AUC scores by Model",
                x = "Model",
                y = "AUC") +
           theme_minimal() +
           theme(legend.position = "none"))
    dev.off()
  }
  
  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = "AUC")
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("AUC" = res_auc, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("AUC" = res_auc, "Top_model" = top_model))
  }
  
}

choose_base_models = function(models, metric = "Accuracy"){
  
  #Bind metrics values from each model
  resample_df = list()
  for (i in 1:length(models)){
    resample_df[[i]] = models[[i]]$resample %>% 
      mutate(model = names(models)[i])
    names(resample_df)[i] = names(models)[i]
  }
  resample_df = do.call(rbind, resample_df)
  
  resample_df <- resample_df %>%
    mutate(Category = case_when(
      model %in% c("BAG", "C50", "CART", "RF") ~ "Tree-based Methods",
      model %in% c("GLM", "LDA", "GLMNET", "LASSO", "RIDGE") ~ "Linear Models",
      model %in% c("KNN", "SVM_linear", "SVM_radial") ~ "Instance-based Methods",
      TRUE ~ "Other"  # In case there are models not in the above lists
    ))
  
  if(metric == "Accuracy"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(Category = unique(Category)[1],    # Keeps the first occurrence of Category for each model
                Accuracy = mean(Accuracy, na.rm = TRUE), # Calculates the mean RMSE
                .groups = 'drop')  # Drops the grouping after summarization
    
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(Accuracy == max(Accuracy)) %>%
      ungroup()     
    
  }else if(metric == "AUC"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(Category = unique(Category)[1],    # Keeps the first occurrence of Category for each model
                AUC = mean(AUC, na.rm = TRUE), # Calculates the mean RMSE
                .groups = 'drop')  # Drops the grouping after summarization
    
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(AUC == max(AUC)) %>%
      ungroup()     
    
  }else if(metric == "RMSE"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(Category = unique(Category)[1],    # Keeps the first occurrence of Category for each model
                RMSE = mean(RMSE, na.rm = TRUE), # Calculates the mean RMSE
                .groups = 'drop')  # Drops the grouping after summarization
    
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(RMSE == min(RMSE)) %>% #Take the model with the min RMSE per category
      ungroup()     
  }
 
  #Retrieve top model based on accuracy/auc/rmse
  base_models <- groupped_df %>%
    pull(model)
  
  return(base_models)
}

calculate_auc_resample = function(obs, pred){

  prob_obs = data.frame("yes" = pred, "obs" = obs)
  
  prob_obs = prob_obs %>%
    arrange(desc(pred)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (obs == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold 
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           fpr = fp/sum(obs == 'no'),
           tpr = tp/sum(obs == 'yes'))

  auc_value = calculate_auc(prob_obs$fpr, prob_obs$tpr)
  
  return(auc_value)
}

get_sensitivity_specificity = function(predictions, observed, ml.model){
  prob_obs = bind_cols(predictions, observed = observed) 

  prob_obs = prob_obs %>%
    arrange(desc(yes)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (observed == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold 
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           sensitivity = tp/sum(observed == 'yes'),
           fpr = fp/sum(observed == 'no'),
           specificity = 1 - fpr) %>%
    select(sensitivity, specificity, fpr) %>%
    mutate(model = ml.model)
    
  starts_at_zero <- any(prob_obs$sensitivity == 0 & prob_obs$fpr == 0)
  
  ##Add dummy row if it doesnt start at 0
  if(!starts_at_zero){
    dummy_row <- data.frame(
      sensitivity = 0,
      specificity = 1,
      fpr = 0,
      model = ml.model
    )

    prob_obs = rbind(dummy_row, prob_obs)
  }
  
  prob_obs = prob_obs %>%
    mutate(Accuracy = calculate_accuracy(., observed),
           Precision = calculate_precision(., observed)) 
    
  
  return(prob_obs)
  
}

#Take sensitivities values based on values of specificities
get_sensitivity = function(x, data){
  data %>%
    filter(specificity - x >= 0)%>% #Take specificity values above threshold x
    top_n(sensitivity, n=1) %>% #Take highest sensitivity from that threshold
    mutate(specificity = x, fpr = 1-x) %>% #Define sensitivity based on the specified threshold
    distinct() #If multiple thresholds have same sensitivity values take only one
}

#Compute customized roc curve
get_roc_curve = function(pred.prob, observed, model, return_roc = T){
  
  sens_spec = get_sensitivity_specificity(pred.prob, observed, model) 
  auc = calculate_auc(sens_spec$fpr, sens_spec$sensitivity)
  
  #sens_spec = map_dfr(specificity, get_sensitivity, sens_spec)
  p = sens_spec %>% 
      group_by(model, specificity) %>%
      summarise(lquartile = quantile(sensitivity, prob =0.25),
                uquartile = quantile(sensitivity, prob = 0.75),
                sensitivity = median(sensitivity),
                .groups = "drop") %>%
      ggplot(aes(x=1-specificity, y=sensitivity, 
                 ymin=lquartile, ymax=uquartile)) +
      geom_ribbon(alpha=0.25, aes(fill=model)) +
      geom_step(aes(color=model)) +
      geom_abline(slope = 1, intercept = 0) +
      theme_classic() +
      theme(legend.position.inside = c(0.8, 0.2))+
      scale_color_discrete(name = "Model",
                           labels = paste(model, 
                                          "\nAUC = ", round(auc, 4))) +
      scale_fill_discrete(name = "Model",
                          labels = paste(model, 
                                         "\nAUC = ", round(auc, 4)))  +
      xlim(0, 1) +  # Ensure the x-axis is from 0 to 1
      ylim(0, 1)    # Ensure the y-axis is from 0 to 1
      
    
  if(return_roc){
    print(p)
  }

  return(list(sens_spec, auc))
}

calculate_auc <- function(fpr, tpr) {
  #tpr = sensitivity 
  
  # Sort by FPR to ensure trapezoidal rule is correctly applied
  # ordered <- order(fpr)
  # fpr <- fpr[ordered]
  # tpr <- tpr[ordered]
  # 
  
  auc <- 0
  for (i in 1:(length(fpr) - 1)) { #-1 to avoid NA cause last terms are TPR = 1 and FPR = 1
    # Trapezoidal rule: (TPR_i + TPR_{i+1}) / 2 * (FPR_{i+1} - FPR_i)
    auc <- auc + ((tpr[i+1] + tpr[i]) / 2) * (fpr[i+1] - fpr[i])
  }
  return(auc)
  
}

compute.prediction = function(model, test_data, target){
  
  cat("Predicting target variable using provided ML model")
  
  predict <- as.numeric(predict(model, test_data))

  prediction_res = evaluate_prediction(predict, target)

  return(list(Metrics = prediction_res, Predictions = prediction_res))
}

compute.prediction.stacked = function(super.learner, test_data, target, ml.models, base.models){
  
  #Learning from simple meta-learner
  base_predictions = list()
  for (i in 1:length(base.models)) {
    base_predictions[[i]] = predict(ml.models[[base.models[i]]], test_data)
    names(base_predictions)[i] = base.models[i]
  }
  
  base_predictions = do.call(cbind, base_predictions)
  
  prediction_simple = as.numeric(predict(super.learner[["simple"]], base_predictions)) 
  
  #Learning from simple meta-learner
  all_predictions = list()
  for (i in 1:length(ml.models)) {
    all_predictions[[i]] = predict(ml.models[[i]], test_data)
    names(all_predictions)[i] = names(ml.models)[i]
  }
  
  all_predictions = do.call(cbind, all_predictions)
  
  prediction_all = as.numeric(predict(super.learner[["all"]], all_predictions)) 
  
  #Metrics
  prediction_res_simple = evaluate_prediction(prediction_simple, target)
  prediction_res_all = evaluate_prediction(prediction_all, target)

  return(list(Metrics = list("Simple" = prediction_res_simple, "All" = prediction_res_all), Predictions = list("Simple" = prediction_simple, "All" = prediction_all)))    
  
}

evaluate_prediction = function(predicted, observed){

  ######Metrics
  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  
  # Calculate MAE
  mae <- mean(abs(observed - predicted))
  
  # Calculate R-squared
  ss_res <- sum((observed - predicted) ^ 2)
  ss_tot <- sum((observed - mean(observed)) ^ 2)
  r_squared <- 1 - (ss_res / ss_tot)
  
  metrics = data.frame(RMSE = rmse,
                       MAE = mae,
                       R_square = r_squared)
  
  data <- data.frame(
    Actual = observed,
    Predicted = predicted
  )
  
  # Fit a linear regression model
  model <- lm(Actual ~ Predicted, data = data)
  
  # Regression plot
  p1 <- ggplot(data, aes(x = Predicted, y = Actual)) +
    geom_point(color = 'blue', size = 2, alpha = 0.6) +  # Scatter plot
    geom_smooth(method = 'lm', color = 'red', se = FALSE, linetype = "dashed") +  # Regression line
    labs(title = 'Regression Plot', x = 'Predicted Values', y = 'Actual Values') +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    ) +
    annotate("text", x = max(data$Predicted) * 0.1, y = max(data$Actual) * 0.8,
             label = paste("RMSE:", round(rmse, 2)), color = "black", size = 4, hjust = 0) +
    annotate("text", x = max(data$Predicted) * 0.1, y = max(data$Actual) * 0.7,
             label = paste("R-squared:", round(r_squared, 3)), color = "black", size = 4, hjust = 0)
  
  print(p1)
  
  # Residuals analysis
  data$residuals <- residuals(model)
  
  p2 <- ggplot(data, aes(x = Predicted, y = residuals)) +
    geom_point(color = 'green', alpha = 0.6, size = 2) +  # Residuals scatter plot
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Horizontal line at y=0
    labs(title = 'Residuals Plot', 
         x = 'Predicted Values', 
         y = 'Residuals',
         subtitle = paste("Mean of Residuals =", round(mean(data$residuals), 2))) +  # Mean as subtitle
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  print(p2)
  
  return(metrics)
}


calculate_accuracy <- function(metrics, target) {
  sensitivity = metrics[,"sensitivity"]
  specificity = metrics[,"specificity"]
  total_positives = sum(target == "yes")
  total_negatives = sum(target == "no")
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  return(accuracy)
}

calculate_precision <- function(metrics, target) {
  sensitivity = metrics[,"sensitivity"]
  specificity = metrics[,"specificity"]
  total_positives = sum(target == "yes")
  total_negatives = sum(target == "no")
  
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN
  
  # Calculate Precision
  precision <- TP / (TP + FP)
  
  return(precision)
}

compute_AUC_hyperparameters_tuning = function(models, ml.predictions, training_df){
  
  for (i in 1:length(models)) {
    #Define model
    ml_model = models[[i]]
    
    if(names(models)[i] == "BAG"){ ################## Bagged CART
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% 
        group_by(Resample) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% #Calculate resamples AUC scores 
        ungroup() 
      
      ## Integrate AUCs into resamples matrix
      auc = c()
      for (i in 1:nrow(ml_model$resample)) {
        auc_val = ml_model$pred %>%
          filter(Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) %>%
        select(AUC, everything())
      
      ## Integrate average CV AUCs into results
      ml_model$results = ml_model$results %>%
        mutate(AUC = mean(ml_model$resample$AUC))
      
      ## Compute predictions based on tuned model
      predictions <- data.frame(predict(ml_model$finalModel, newdata = training_df, type = "prob")) %>% #Predictions using tuned model
        select(yes) %>%
        rename(BAG = yes) 
      
    }else if(names(models)[i] == "RF"){ ################## RF
  
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>%
        group_by(Resample, mtry) %>% #Parameters for tunning
        mutate(AUC = calculate_auc_resample(obs, yes)) %>%
        ungroup()
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(mtry) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = "mtry")
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(mtry = ml_model$results$mtry[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(mtry == as.numeric(ml_model$bestTune),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
      
      predictions = data.frame(predict(ml_model$finalModel, newdata = training_df, type = "prob"))[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(RF = yes) #Predictions of model (already ordered)
      
    }else if(names(models)[i] == "C50"){ ################## C5.0
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>%
        group_by(trials, model, winnow) %>% #Parameters for tunning
        mutate(AUC = calculate_auc_resample(obs, yes)) %>%
        ungroup()
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(trials, model, winnow) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = c("trials", "model", "winnow"))
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(trials = ml_model$results$trials[tune],
               model = ml_model$results$model[tune],
               winnow = ml_model$results$winnow[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(trials == as.numeric(ml_model$bestTune$trials),
                 model == as.character(ml_model$bestTune$model),
                 winnow == as.character(ml_model$bestTune$winnow),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
      
      predictions = data.frame(predict(ml_model$finalModel, newdata = training_df, type = "prob"))[,"yes", drop=F]  %>% 
        select(yes) %>%
        rename(C50 = yes)  #Predictions of model (already ordered)
      
    }else if(names(models)[i] == "GLM"){   ################## LG - Logistic Regression
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into resamples matrix
      auc = c()
      for (i in 1:nrow(ml_model$resample)) {
        auc_val = ml_model$pred %>%
          filter(Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        
        auc = c(auc, auc_val)
      }
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
      
      ## Integrate average CV AUCs into results
      ml_model$results = ml_model$results %>%
        mutate(AUC = mean(ml_model$resample$AUC))
      
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(GLM = yes)  #Predictions of model (already ordered)
  
    }else if(names(models)[i] == "LDA"){ ################## LDA - Linear Discriminate Analysis
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into resamples matrix
      auc = c()
      for (i in 1:nrow(ml_model$resample)) {
        auc_val = ml_model$pred %>%
          filter(Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        
        auc = c(auc, auc_val)
      }
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
      
      ## Integrate average CV AUCs into results
      ml_model$results = ml_model$results %>%
        mutate(AUC = mean(ml_model$resample$AUC))
      
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(LDA = yes)  #Predictions of model (already ordered)
      
    }else if(names(models)[i] == "GLMNET"){ ################## GLMNET - Regularized Logistic Regression (Elastic net)
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample, alpha, lambda) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(alpha, lambda) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = c("alpha", "lambda"))
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(alpha = ml_model$results$alpha[tune],
               lambda = ml_model$results$lambda[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(alpha == as.numeric(ml_model$bestTune$alpha),
                 lambda == as.numeric(ml_model$bestTune$lambda),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
      
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(GLMNET = yes)  #Predictions of model (already ordered)
      
  
    }else if(names(models)[i] == "KNN"){ ################## KNN - k-Nearest Neighbors
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample, k) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(k) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = "k")
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(k = ml_model$results$k[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(k == as.numeric(ml_model$bestTune$k),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc)
      
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(KNN = yes) #Predictions of model (already ordered)
      
    }else if(names(models)[i] == "CART"){ ################## CART - Classification and Regression Trees (CART)
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample, cp) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(cp) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = "cp")
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(cp = ml_model$results$cp[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(cp == as.numeric(ml_model$bestTune$cp),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
  
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(CART = yes)  #Predictions of model (already ordered)
      
    }else if(names(models)[i]=="LASSO"){ ################## Regularized Lasso
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample, alpha, lambda) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(alpha, lambda) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = c("alpha", "lambda"))
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(alpha = ml_model$results$alpha[tune],
               lambda = ml_model$results$lambda[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(alpha == as.numeric(ml_model$bestTune$alpha),
                 lambda == as.numeric(ml_model$bestTune$lambda),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
      
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(LASSO = yes)  #Predictions of model (already ordered)
      
    }else if(names(models)[i]=="RIDGE"){   ################## Ridge regression
  
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample, alpha, lambda) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(alpha, lambda) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = c("alpha", "lambda"))
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(alpha = ml_model$results$alpha[tune],
               lambda = ml_model$results$lambda[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(alpha == as.numeric(ml_model$bestTune$alpha),
                 lambda == as.numeric(ml_model$bestTune$lambda),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc) 
     
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(RIDGE = yes)  #Predictions of model (already ordered)
      
    }else if(names(models)[i] == "SVM_radial"){ ################## Support Vector Machine with Radial Kernel
      
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample, sigma, C) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(sigma, C) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = c("sigma", "C"))
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(sigma = ml_model$results$sigma[tune],
               C = ml_model$results$C[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(sigma == as.numeric(ml_model$bestTune$sigma),
                 C == as.numeric(ml_model$bestTune$C),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc)
      
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(SVM_radial = yes)  #Predictions of model (already ordered)
      
    }else if(names(models)[i] == "SVM_linear"){
      ## Integrate AUCs into prediction matrix
      ml_model$pred = ml_model$pred %>% #Calculate resamples AUC scores 
        group_by(Resample, C) %>%
        mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
        ungroup() 
      
      ## Integrate AUCs into results per parameter 
      auc_values = ml_model$pred %>%
        group_by(C) %>%
        summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
      
      ml_model[["results"]] <- ml_model[["results"]] %>%
        left_join(auc_values, by = "C")
      
      #Tuning parameter (select combination with top AUC)
      tune = which.max(ml_model$results$AUC)
      ml_model$bestTune = ml_model$bestTune %>%
        mutate(C = ml_model$results$C[tune])
      
      #Configure resamples to have the AUCs only using tuned parameter
      ml_model$resample = ml_model$resample[order(ml_model$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
      auc = c()
      for (i in 1:nrow(ml_model$resample)){
        auc_val = ml_model$pred %>%
          filter(C == as.numeric(ml_model$bestTune$C),
                 Resample == ml_model$resample$Resample[i]) %>%
          pull(AUC) %>% #AUC per resample is the same
          unique()
        auc = c(auc, auc_val)
      }
      
      ml_model$resample = ml_model$resample %>%
        mutate(AUC = auc)
      
      ################## Support Vector Machine with Linear Kernel
      predictions = predict(ml_model, newdata = training_df, type = "prob")[,"yes", drop=F]  %>%
        select(yes) %>%
        rename(SVM_linear = yes)  #Predictions of model (already ordered)
      
    }
    
    models[[i]] = ml_model
    ml.predictions[[i]] = predictions
  } 
  
  return(list(Models = models, ML_predictions = ml.predictions))
  
}
