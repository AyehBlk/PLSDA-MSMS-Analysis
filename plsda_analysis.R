# ============================================================================
# PLS-DA ANALYSIS FOR IDA MS/MS SPECTROMETRY DATA
# Partial Least Squares Discriminant Analysis
# ============================================================================

cat(paste(rep("=", 70), collapse=""), "\n")
cat("PLS-DA SUPERVISED CLUSTERING ANALYSIS\n")
cat("MS/MS IDA DATASET\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# ============================================================================
# SECTION 1: CREATE EXAMPLE MS/MS IDA DATASET
# ============================================================================

set.seed(123)

create_msms_data <- function(n_samples = 150, n_features = 50) {
  cat("Creating simulated MS/MS IDA dataset...\n")
  
  n_per_class <- n_samples / 3
  
  # Control group
  control <- matrix(rnorm(n_per_class * n_features, mean = 100, sd = 20),
                    nrow = n_per_class, ncol = n_features)
  
  # Treatment1 group (features 1-15 upregulated)
  treatment1 <- matrix(rnorm(n_per_class * n_features, mean = 100, sd = 20),
                       nrow = n_per_class, ncol = n_features)
  treatment1[, 1:15] <- treatment1[, 1:15] + 50
  
  # Treatment2 group (features 16-30 upregulated, 1-10 downregulated)
  treatment2 <- matrix(rnorm(n_per_class * n_features, mean = 100, sd = 20),
                       nrow = n_per_class, ncol = n_features)
  treatment2[, 16:30] <- treatment2[, 16:30] + 60
  treatment2[, 1:10] <- treatment2[, 1:10] - 30
  
  # Combine and ensure positive values
  data_matrix <- rbind(control, treatment1, treatment2)
  data_matrix <- abs(data_matrix + rnorm(length(data_matrix), 0, 5))
  
  # Create feature names
  mz_values <- seq(200, 1000, length.out = n_features)
  colnames(data_matrix) <- paste0("mz_", round(mz_values, 2))
  
  # Sample information
  sample_info <- data.frame(
    SampleID = paste0("Sample_", 1:n_samples),
    Class = factor(rep(c("Control", "Treatment1", "Treatment2"), 
                      each = n_per_class)),
    Batch = factor(rep(1:3, length.out = n_samples)),
    stringsAsFactors = FALSE
  )
  
  complete_data <- cbind(sample_info, data.frame(data_matrix))
  
  list(data = complete_data, 
       features = colnames(data_matrix),
       feature_matrix = data_matrix)
}

# Generate dataset
msms_data <- create_msms_data(n_samples = 150, n_features = 50)
data_full <- msms_data$data
features <- msms_data$feature_matrix
labels <- data_full$Class

cat("Dataset created:\n")
cat("- Samples:", nrow(features), "\n")
cat("- Features:", ncol(features), "\n")
cat("- Classes:", paste(levels(labels), collapse=", "), "\n")
print(table(labels))
cat("\n")


# ============================================================================
# SECTION 2: DATA PREPROCESSING
# ============================================================================

cat("Preprocessing data...\n")

# Log transformation (common for MS data)
features_log <- log2(features + 1)

# Mean centering and scaling
features_scaled <- scale(features_log, center = TRUE, scale = TRUE)

cat("- Log2 transformation applied\n")
cat("- Mean centering and unit variance scaling completed\n\n")


# ============================================================================
# SECTION 3: TRAIN-TEST SPLIT (STRATIFIED)
# ============================================================================

set.seed(456)
train_prop <- 0.7

train_idx <- c()
for(class_name in levels(labels)) {
  class_idx <- which(labels == class_name)
  n_train <- floor(length(class_idx) * train_prop)
  train_idx <- c(train_idx, sample(class_idx, n_train))
}

X_train <- features_scaled[train_idx, ]
y_train <- labels[train_idx]
X_test <- features_scaled[-train_idx, ]
y_test <- labels[-train_idx]

cat("Train-test split (70-30):\n")
cat("Training samples:", nrow(X_train), "\n")
cat("Testing samples:", nrow(X_test), "\n\n")


# ============================================================================
# SECTION 4: PLS-DA IMPLEMENTATION
# ============================================================================

cat("Building PLS-DA model...\n\n")

# Convert class labels to dummy matrix for PLS regression
create_dummy_matrix <- function(y) {
  classes <- levels(y)
  n <- length(y)
  Y <- matrix(0, nrow = n, ncol = length(classes))
  colnames(Y) <- classes
  
  for(i in 1:n) {
    Y[i, which(classes == y[i])] <- 1
  }
  
  return(Y)
}

# PLS algorithm (NIPALS algorithm)
pls_nipals <- function(X, Y, ncomp = 2) {
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # Storage for components
  T <- matrix(0, n, ncomp)  # X scores
  P <- matrix(0, p, ncomp)  # X loadings
  W <- matrix(0, p, ncomp)  # X weights
  U <- matrix(0, n, ncomp)  # Y scores
  Q <- matrix(0, q, ncomp)  # Y loadings
  C <- matrix(0, q, ncomp)  # Y weights
  
  # Regression coefficients
  B <- array(0, dim = c(p, q, ncomp))
  
  X_residual <- X
  Y_residual <- Y
  
  for(comp in 1:ncomp) {
    
    # Initialize u as first column of Y
    u <- Y_residual[, 1, drop = FALSE]
    
    # NIPALS iterations
    for(iter in 1:100) {
      # X weights
      w <- t(X_residual) %*% u
      w <- w / sqrt(sum(w^2))  # Normalize
      
      # X scores
      t_score <- X_residual %*% w
      
      # Y weights
      c_weight <- t(Y_residual) %*% t_score
      c_weight <- c_weight / sqrt(sum(c_weight^2))  # Normalize
      
      # Y scores
      u_new <- Y_residual %*% c_weight
      
      # Check convergence
      if(sum((u_new - u)^2) < 1e-10) break
      u <- u_new
    }
    
    # X loadings
    p_loading <- t(X_residual) %*% t_score / as.numeric(t(t_score) %*% t_score)
    
    # Y loadings
    q_loading <- t(Y_residual) %*% t_score / as.numeric(t(t_score) %*% t_score)
    
    # Store components
    T[, comp] <- t_score
    P[, comp] <- p_loading
    W[, comp] <- w
    U[, comp] <- u
    Q[, comp] <- q_loading
    C[, comp] <- c_weight
    
    # Regression coefficient for this component
    B[, , comp] <- w %*% t(q_loading)
    
    # Deflate X and Y
    X_residual <- X_residual - t_score %*% t(p_loading)
    Y_residual <- Y_residual - t_score %*% t(q_loading)
  }
  
  # Calculate explained variance
  TSS_X <- sum(X^2)
  RSS_X <- sum(X_residual^2)
  R2X <- (TSS_X - RSS_X) / TSS_X
  
  TSS_Y <- sum(Y^2)
  RSS_Y <- sum(Y_residual^2)
  R2Y <- (TSS_Y - RSS_Y) / TSS_Y
  
  list(
    T = T,           # X scores
    P = P,           # X loadings
    W = W,           # X weights
    U = U,           # Y scores
    Q = Q,           # Y loadings
    C = C,           # Y weights
    B = B,           # Regression coefficients
    R2X = R2X,       # Variance explained in X
    R2Y = R2Y        # Variance explained in Y
  )
}

# Train PLS-DA model
Y_train <- create_dummy_matrix(y_train)
Y_test <- create_dummy_matrix(y_test)

ncomp <- min(10, nrow(X_train) - 1, ncol(X_train))  # Number of components

plsda_model <- pls_nipals(X_train, Y_train, ncomp = ncomp)

cat("PLS-DA model built with", ncomp, "components\n")
cat("R2X (variance explained in X):", round(plsda_model$R2X, 4), "\n")
cat("R2Y (variance explained in Y):", round(plsda_model$R2Y, 4), "\n\n")


# ============================================================================
# SECTION 5: PREDICTIONS AND EVALUATION
# ============================================================================

cat("Making predictions...\n")

# Predict function for PLS-DA
predict_plsda <- function(model, X_new, ncomp_use = 2) {
  
  # Calculate scores for new data
  T_new <- X_new %*% model$W[, 1:ncomp_use, drop = FALSE]
  
  # Predict Y values
  Y_pred <- T_new %*% t(model$Q[, 1:ncomp_use, drop = FALSE])
  
  # Convert to class predictions (highest value)
  class_pred <- apply(Y_pred, 1, which.max)
  
  list(Y_pred = Y_pred, class_pred = class_pred, scores = T_new)
}

# Optimize number of components using cross-validation
cat("Performing cross-validation to select optimal components...\n")

cv_folds <- 5
cv_accuracy <- numeric(ncomp)

for(n in 1:ncomp) {
  fold_acc <- numeric(cv_folds)
  folds <- cut(seq(1, nrow(X_train)), breaks = cv_folds, labels = FALSE)
  
  for(f in 1:cv_folds) {
    test_idx <- which(folds == f)
    train_idx <- which(folds != f)
    
    # Train on CV fold
    cv_model <- pls_nipals(X_train[train_idx, ], 
                           Y_train[train_idx, ], 
                           ncomp = n)
    
    # Predict on validation fold
    cv_pred <- predict_plsda(cv_model, X_train[test_idx, ], ncomp_use = n)
    
    # Calculate accuracy
    pred_classes <- colnames(Y_train)[cv_pred$class_pred]
    actual_classes <- as.character(y_train[test_idx])
    fold_acc[f] <- mean(pred_classes == actual_classes)
  }
  
  cv_accuracy[n] <- mean(fold_acc)
}

optimal_ncomp <- which.max(cv_accuracy)
cat("Optimal number of components:", optimal_ncomp, "\n")
cat("Cross-validation accuracy:", round(cv_accuracy[optimal_ncomp], 4), "\n\n")

# Final predictions on test set
test_pred <- predict_plsda(plsda_model, X_test, ncomp_use = optimal_ncomp)
pred_classes <- colnames(Y_test)[test_pred$class_pred]
actual_classes <- as.character(y_test)

# Confusion matrix
conf_matrix <- table(Predicted = factor(pred_classes, levels = levels(y_test)), 
                     Actual = y_test)

accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)

cat(paste(rep("=", 70), collapse=""), "\n")
cat("TEST SET RESULTS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("Overall Accuracy:", round(accuracy, 4), "\n\n")

cat("Confusion Matrix:\n")
print(conf_matrix)
cat("\n")

# Per-class metrics
classes <- levels(y_test)
for(i in 1:length(classes)) {
  TP <- conf_matrix[i, i]
  FP <- sum(conf_matrix[i, ]) - TP
  FN <- sum(conf_matrix[, i]) - TP
  
  sensitivity <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  
  cat(sprintf("%s - Sensitivity: %.4f, Precision: %.4f\n", 
              classes[i], sensitivity, precision))
}
cat("\n")


# ============================================================================
# SECTION 6: VARIABLE IMPORTANCE (VIP SCORES)
# ============================================================================

cat("Calculating Variable Importance in Projection (VIP) scores...\n")

calculate_vip <- function(model, ncomp_use) {
  
  W <- model$W[, 1:ncomp_use, drop = FALSE]
  T <- model$T[, 1:ncomp_use, drop = FALSE]
  Q <- model$Q[, 1:ncomp_use, drop = FALSE]
  
  p <- nrow(W)
  
  # Sum of squares of Y explained by each component
  SSY <- colSums((T %*% t(Q))^2)
  
  # VIP calculation
  VIP <- sqrt(p * rowSums((W^2) %*% diag(SSY)) / sum(SSY))
  
  return(VIP)
}

vip_scores <- calculate_vip(plsda_model, optimal_ncomp)

vip_df <- data.frame(
  Feature = colnames(X_train),
  VIP = vip_scores
)
vip_df <- vip_df[order(-vip_df$VIP), ]

cat("\nTop 20 Most Important Features (VIP > 1 are important):\n")
print(head(vip_df, 20), row.names = FALSE)
cat("\n")

important_features <- vip_df$Feature[vip_df$VIP > 1]
cat("Number of important features (VIP > 1):", length(important_features), "\n\n")


# ============================================================================
# SECTION 7: VISUALIZATIONS
# ============================================================================

cat("Generating visualizations...\n")

# Plot 1: Scores Plot (Component 1 vs 2)
png("01_plsda_scores_plot.png", width = 1000, height = 800, res = 100)
par(mar = c(5, 5, 4, 8), xpd = TRUE)

colors <- c("Control" = "blue", "Treatment1" = "red", "Treatment2" = "green")
shapes <- c("Control" = 16, "Treatment1" = 17, "Treatment2" = 15)

# Calculate variance explained for each component
var_explained <- numeric(ncomp)
for(i in 1:ncomp) {
  var_explained[i] <- sum(plsda_model$T[, i]^2) / sum(X_train^2) * 100
}

plot(plsda_model$T[, 1], plsda_model$T[, 2],
     col = colors[as.character(y_train)],
     pch = shapes[as.character(y_train)],
     cex = 1.5,
     xlab = paste0("Component 1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("Component 2 (", round(var_explained[2], 1), "%)"),
     main = "PLS-DA Scores Plot - Training Set")

legend("topright", inset = c(-0.2, 0),
       legend = levels(y_train),
       col = colors[levels(y_train)],
       pch = shapes[levels(y_train)],
       cex = 1.2,
       title = "Class")

# Add 95% confidence ellipses
for(class in levels(y_train)) {
  idx <- which(y_train == class)
  if(length(idx) > 2) {
    scores_class <- plsda_model$T[idx, 1:2]
    
    # Calculate ellipse
    center <- colMeans(scores_class)
    cov_mat <- cov(scores_class)
    
    # Eigenvalues and eigenvectors
    eig <- eigen(cov_mat)
    
    # 95% confidence (chi-square with 2 df)
    radius <- sqrt(qchisq(0.95, df = 2))
    
    # Ellipse points
    theta <- seq(0, 2*pi, length.out = 100)
    ellipse <- matrix(0, 100, 2)
    for(i in 1:100) {
      point <- c(cos(theta[i]), sin(theta[i]))
      ellipse[i, ] <- center + radius * (eig$vectors %*% diag(sqrt(eig$values)) %*% point)
    }
    
    lines(ellipse[, 1], ellipse[, 2], 
          col = colors[class], lwd = 2, lty = 2)
  }
}

dev.off()


# Plot 2: Test Set Predictions
png("02_plsda_test_predictions.png", width = 1000, height = 800, res = 100)
par(mar = c(5, 5, 4, 8), xpd = TRUE)

plot(test_pred$scores[, 1], test_pred$scores[, 2],
     col = colors[actual_classes],
     pch = shapes[actual_classes],
     cex = 1.5,
     xlab = paste0("Component 1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("Component 2 (", round(var_explained[2], 1), "%)"),
     main = "PLS-DA Test Set Predictions")

# Mark misclassifications with circles
misclass_idx <- which(pred_classes != actual_classes)
if(length(misclass_idx) > 0) {
  points(test_pred$scores[misclass_idx, 1], 
         test_pred$scores[misclass_idx, 2],
         cex = 2.5, col = "black", pch = 1, lwd = 2)
}

legend("topright", inset = c(-0.2, 0),
       legend = c(levels(y_test), "Misclassified"),
       col = c(colors[levels(y_test)], "black"),
       pch = c(shapes[levels(y_test)], 1),
       cex = 1.2)

dev.off()


# Plot 3: VIP Scores
png("03_plsda_vip_scores.png", width = 1000, height = 800, res = 100)
par(mar = c(5, 12, 4, 2))

top_vip <- head(vip_df, 20)
barplot(top_vip$VIP[20:1], 
        names.arg = top_vip$Feature[20:1],
        horiz = TRUE,
        las = 1,
        col = ifelse(top_vip$VIP[20:1] > 1, "darkred", "steelblue"),
        xlab = "VIP Score",
        main = "Top 20 Features by Variable Importance")
abline(v = 1, lty = 2, lwd = 2, col = "red")
legend("bottomright", 
       legend = c("VIP > 1 (Important)", "VIP < 1"),
       fill = c("darkred", "steelblue"))
dev.off()


# Plot 4: Loadings Plot
png("04_plsda_loadings_plot.png", width = 1000, height = 800, res = 100)
par(mar = c(5, 5, 4, 2))

plot(plsda_model$P[, 1], plsda_model$P[, 2],
     pch = 16, col = "darkblue", cex = 1,
     xlab = "Component 1 Loadings",
     ylab = "Component 2 Loadings",
     main = "PLS-DA Loadings Plot")

# Highlight important features
important_idx <- which(vip_scores > 1)
points(plsda_model$P[important_idx, 1], 
       plsda_model$P[important_idx, 2],
       pch = 16, col = "red", cex = 1.5)

# Label top 10 features
top_10_idx <- order(-vip_scores)[1:10]
text(plsda_model$P[top_10_idx, 1], 
     plsda_model$P[top_10_idx, 2],
     labels = colnames(X_train)[top_10_idx],
     pos = 3, cex = 0.7, col = "darkred")

legend("topright",
       legend = c("All features", "VIP > 1", "Top 10"),
       col = c("darkblue", "red", "darkred"),
       pch = c(16, 16, NA),
       text.col = c("black", "black", "darkred"))

abline(h = 0, v = 0, lty = 2, col = "gray")
dev.off()


# Plot 5: Cross-Validation Results
png("05_plsda_cv_results.png", width = 800, height = 600, res = 100)
par(mar = c(5, 5, 4, 2))

plot(1:ncomp, cv_accuracy, type = "b",
     pch = 16, col = "darkblue", lwd = 2, cex = 1.5,
     xlab = "Number of Components",
     ylab = "Cross-Validation Accuracy",
     main = "PLS-DA: Model Selection via CV",
     ylim = c(0, 1))

points(optimal_ncomp, cv_accuracy[optimal_ncomp],
       pch = 16, col = "red", cex = 2)

text(optimal_ncomp, cv_accuracy[optimal_ncomp],
     labels = paste("Optimal:", optimal_ncomp),
     pos = 3, col = "red", font = 2)

grid()
dev.off()


# Plot 6: Confusion Matrix Heatmap
png("06_plsda_confusion_matrix.png", width = 800, height = 700, res = 100)
par(mar = c(5, 6, 4, 2))

image(1:ncol(conf_matrix), 1:nrow(conf_matrix), 
      t(conf_matrix[nrow(conf_matrix):1, ]),
      col = colorRampPalette(c("white", "lightblue", "darkblue"))(20),
      xlab = "Actual Class", ylab = "Predicted Class",
      main = "PLS-DA Confusion Matrix (Test Set)",
      axes = FALSE)

axis(1, at = 1:ncol(conf_matrix), labels = colnames(conf_matrix))
axis(2, at = 1:nrow(conf_matrix), labels = rev(rownames(conf_matrix)), las = 1)

# Add text labels
for(i in 1:nrow(conf_matrix)) {
  for(j in 1:ncol(conf_matrix)) {
    text(j, nrow(conf_matrix) - i + 1, conf_matrix[i, j], cex = 2)
  }
}

dev.off()


# Plot 7: Component Variance Explained
png("07_plsda_variance_explained.png", width = 800, height = 600, res = 100)
par(mar = c(5, 5, 4, 2))

barplot(var_explained, 
        names.arg = paste0("Comp", 1:ncomp),
        col = "steelblue",
        xlab = "Component",
        ylab = "Variance Explained (%)",
        main = "Variance Explained by Each Component")

# Add cumulative line
cumvar <- cumsum(var_explained)
lines(seq(0.7, by = 1.2, length.out = ncomp), 
      cumvar, col = "red", lwd = 2, type = "b", pch = 16)

legend("right",
       legend = c("Individual", "Cumulative"),
       fill = c("steelblue", NA),
       border = c("black", NA),
       lty = c(NA, 1),
       col = c(NA, "red"),
       lwd = c(NA, 2))

dev.off()


# ============================================================================
# SECTION 8: EXPORT RESULTS
# ============================================================================

cat("Exporting results...\n")

# Export predictions
results_df <- data.frame(
  SampleID = data_full$SampleID[-train_idx],
  Actual = actual_classes,
  Predicted = pred_classes,
  Correct = actual_classes == pred_classes,
  Component1_Score = test_pred$scores[, 1],
  Component2_Score = test_pred$scores[, 2]
)

write.csv(results_df, "plsda_predictions.csv", row.names = FALSE)

# Export VIP scores
write.csv(vip_df, "plsda_vip_scores.csv", row.names = FALSE)

# Export model summary
model_summary <- data.frame(
  Metric = c("Optimal_Components", "CV_Accuracy", "Test_Accuracy", 
             "R2X", "R2Y", "Important_Features"),
  Value = c(optimal_ncomp, 
            round(cv_accuracy[optimal_ncomp], 4),
            round(accuracy, 4),
            round(plsda_model$R2X, 4),
            round(plsda_model$R2Y, 4),
            length(important_features))
)

write.csv(model_summary, "plsda_model_summary.csv", row.names = FALSE)

# Export loadings
loadings_df <- data.frame(
  Feature = colnames(X_train),
  Component1_Loading = plsda_model$P[, 1],
  Component2_Loading = plsda_model$P[, 2],
  VIP = vip_scores
)

write.csv(loadings_df, "plsda_loadings.csv", row.names = FALSE)


# ============================================================================
# SECTION 9: FINAL SUMMARY
# ============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("PLS-DA ANALYSIS COMPLETED!\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("DATASET SUMMARY:\n")
cat("- Total samples:", nrow(features), "\n")
cat("- Features:", ncol(features), "\n")
cat("- Classes:", length(levels(labels)), "\n\n")

cat("MODEL PERFORMANCE:\n")
cat("- Optimal components:", optimal_ncomp, "\n")
cat("- Cross-validation accuracy:", round(cv_accuracy[optimal_ncomp], 4), "\n")
cat("- Test set accuracy:", round(accuracy, 4), "\n")
cat("- R2X (variance in X):", round(plsda_model$R2X, 4), "\n")
cat("- R2Y (variance in Y):", round(plsda_model$R2Y, 4), "\n")
cat("- Important features (VIP>1):", length(important_features), "\n\n")

cat("OUTPUT FILES GENERATED:\n")
cat("1. 01_plsda_scores_plot.png - Training set scores with ellipses\n")
cat("2. 02_plsda_test_predictions.png - Test set predictions\n")
cat("3. 03_plsda_vip_scores.png - Variable importance ranking\n")
cat("4. 04_plsda_loadings_plot.png - Feature loadings\n")
cat("5. 05_plsda_cv_results.png - Cross-validation results\n")
cat("6. 06_plsda_confusion_matrix.png - Prediction accuracy heatmap\n")
cat("7. 07_plsda_variance_explained.png - Component variance\n")
cat("8. plsda_predictions.csv - Test predictions with scores\n")
cat("9. plsda_vip_scores.csv - All VIP scores ranked\n")
cat("10. plsda_model_summary.csv - Model performance metrics\n")
cat("11. plsda_loadings.csv - Feature loadings and VIP\n")

cat("\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Analysis complete! All files saved to working directory.\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")