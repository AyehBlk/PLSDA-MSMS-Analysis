# PLS-DA Analysis for MS/MS IDA Data
## Complete Guide and Documentation

---

## Table of Contents

1. [Introduction](#introduction)
2. [What is PLS-DA?](#what-is-pls-da)
3. [Why Use PLS-DA for MS/MS Data?](#why-use-pls-da-for-msms-data)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Understanding the Workflow](#understanding-the-workflow)
7. [Output Files Explained](#output-files-explained)
8. [Interpreting Results](#interpreting-results)
9. [Using Your Own Data](#using-your-own-data)
10. [Parameter Tuning](#parameter-tuning)
11. [Best Practices](#best-practices)
12. [Troubleshooting](#troubleshooting)
13. [Statistical Background](#statistical-background)
14. [References](#references)

---

## Introduction

This package provides a complete implementation of **PLS-DA (Partial Least Squares Discriminant Analysis)** specifically designed for MS/MS (tandem mass spectrometry) IDA (Information Dependent Acquisition) data analysis.

### What You Get:
-  Complete PLS-DA implementation from scratch (no external packages needed)
-  Cross-validation for model optimization
-  Variable Importance in Projection (VIP) scores
-  7 publication-ready visualizations
-  4 detailed CSV output files
-  Example dataset for testing
-  Comprehensive documentation

### Key Features:
- **Dimensionality Reduction**: Reduces thousands of m/z features to interpretable components
- **Biomarker Discovery**: Identifies most important features via VIP scores
- **Classification**: Accurately classifies samples into treatment groups
- **Visualization**: Creates beautiful plots for publications
- **Interpretability**: Provides loadings and scores for biological insight

---

## What is PLS-DA?

### Simple Explanation:
PLS-DA is a supervised machine learning method that:
1. Finds patterns in your MS/MS data that **best separate** different groups
2. Reduces complex data to a few **meaningful components**
3. Identifies which m/z features are **most important** for classification
4. Creates **visual representations** of group separation

### Technical Definition:
PLS-DA combines:
- **PLS Regression** (Partial Least Squares): Finds latent variables that explain variance in both X (features) and Y (classes)
- **Discriminant Analysis**: Maximizes separation between predefined groups

### Think of it as:
Imagine you have 50 different measurements (m/z ratios) for each sample. PLS-DA finds the "super-measurements" (components) that combine these 50 measurements in the best way to tell your groups apart.

---

## Why Use PLS-DA for MS/MS Data?

### PLS-DA is Ideal for MS/MS Because:

#### 1. **Handles High Dimensionality**
- MS/MS data: 100s-1000s of m/z features
- PLS-DA: Works excellently with p >> n (more features than samples)

#### 2. **Deals with Multicollinearity**
- MS/MS features are often correlated (related metabolites)
- PLS-DA handles this naturally (unlike methods that assume independence)

#### 3. **Provides Interpretable Results**
- **Scores plots**: Visualize sample clustering
- **Loadings plots**: See which m/z ratios drive separation
- **VIP scores**: Rank biomarker candidates

#### 4. **Dimensionality Reduction**
- Reduces 1000s of features to 2-5 components
- Easy to visualize in 2D/3D
- Removes noise while retaining signal

#### 5. **Supervised Classification**
- Uses group information to maximize separation
- Better than unsupervised PCA for classification
- Provides class predictions with confidence

### Comparison with Other Methods:

| Method | Best For | Limitations |
|--------|----------|-------------|
| **PLS-DA** | High-dimensional, correlated data; biomarker discovery | Requires class labels; can overfit |
| **PCA** | Exploratory analysis, outlier detection | Unsupervised (ignores groups) |
| **LDA** | Low-dimensional data | Assumes normality; struggles with p > n |
| **Random Forest** | Non-linear relationships | Black box; less interpretable |
| **SVM** | Complex decision boundaries | No feature importance; hard to interpret |

---

## Installation

### Requirements:
- **R version**: 4.0.0 or higher
- **Packages**: Base R only (MASS package - comes with R)
- **RAM**: Minimum 4GB (8GB recommended)
- **Storage**: ~10MB for outputs

### Installation Steps:
```bash
# Step 1: Install R (if not already installed)
# Ubuntu/Debian:
sudo apt-get update
sudo apt-get install r-base

# macOS (with Homebrew):
brew install r

# Windows: Download from https://cran.r-project.org/

# Step 2: No additional packages needed!
# The script uses only base R functions
```

---

## Quick Start

### 1. Save the Script
Save the PLS-DA code as `plsda_analysis.R`

### 2. Run It
```r
# In R console:
source("plsda_analysis.R")

# Or from terminal:
Rscript plsda_analysis.R
```

### 3. Check Results
The script creates 11 files in your working directory:
- 7 PNG plots
- 4 CSV data files

### 4. View Main Results
```r
# View predictions
read.csv("plsda_predictions.csv")

# View important features
read.csv("plsda_vip_scores.csv")

# View model summary
read.csv("plsda_model_summary.csv")
```

That's it! The example runs in ~5-10 seconds.

---

## Understanding the Workflow

### Step-by-Step Process:
```
┌─────────────────────────────────────────────┐
│  1. DATA LOADING                            │
│  - Load MS/MS intensity data                │
│  - 150 samples × 50 m/z features            │
│  - 3 classes (Control, Treatment1, Treatment2) │
└────────────────┬────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────┐
│  2. PREPROCESSING                           │
│  - Log2 transformation (normalize intensities)│
│  - Mean centering (subtract mean)           │
│  - Unit variance scaling (divide by SD)     │
└────────────────┬────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────┐
│  3. TRAIN-TEST SPLIT                        │
│  - 70% training (105 samples)               │
│  - 30% testing (45 samples)                 │
│  - Stratified (maintains class proportions) │
└────────────────┬────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────┐
│  4. PLS-DA MODEL BUILDING                   │
│  - NIPALS algorithm (iterative)             │
│  - Extract latent components                │
│  - Calculate scores, loadings, weights      │
│  - Build regression model                   │
└────────────────┬────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────┐
│  5. CROSS-VALIDATION                        │
│  - 5-fold CV on training set                │
│  - Test 1-10 components                     │
│  - Select optimal number (highest accuracy) │
└────────────────┬────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────┐
│  6. TEST SET PREDICTION                     │
│  - Apply model to test samples              │
│  - Calculate scores for new samples         │
│  - Predict class memberships                │
│  - Generate confusion matrix                │
└────────────────┬────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────┐
│  7. VARIABLE IMPORTANCE                     │
│  - Calculate VIP scores                     │
│  - VIP > 1 = important features             │
│  - Rank all m/z features                    │
│  - Identify biomarker candidates            │
└────────────────┬────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────┐
│  8. VISUALIZATION & EXPORT                  │
│  - Generate 7 plots                         │
│  - Export 4 CSV files                       │
│  - Create summary report                    │
└─────────────────────────────────────────────┘
```

---

## Output Files Explained

### Visualization Files (PNG)

#### 1. **01_plsda_scores_plot.png**
**What it shows:**
- Each point = one sample from training set
- X-axis = Component 1 (most important)
- Y-axis = Component 2 (second most important)
- Colors = actual class labels
- Dashed ellipses = 95% confidence regions

**How to interpret:**
-  **Good**: Classes form separate clusters
-  **Excellent**: No overlap between ellipses
-  **Warning**: Overlapping points suggest similar samples
-  **Problem**: Complete overlap means classes are indistinguishable

**Example interpretation:**
```
If Control (blue) clusters in bottom-right,
Treatment1 (red) in top-right,
Treatment2 (green) in left,
→ Your treatments have distinct metabolic profiles!
```

#### 2. **02_plsda_test_predictions.png**
**What it shows:**
- Test set samples projected onto PLS-DA space
- Colors = actual classes
- Black circles = misclassified samples

**How to interpret:**
- Check if test samples fall within training ellipses
- Misclassified samples (circled) may be:
  - Outliers
  - Mislabeled
  - Representing biological variation

**Use this to:**
- Validate model performance
- Identify problematic samples
- Assess model generalization

#### 3. **03_plsda_vip_scores.png**
**What it shows:**
- Top 20 m/z features ranked by importance
- Red bars = VIP > 1 (important)
- Blue bars = VIP < 1 (less important)
- Vertical dashed line at VIP = 1

**How to interpret:**
- **VIP > 1**: Feature is important for classification
- **VIP > 1.5**: Very important, strong biomarker candidate
- **VIP < 1**: Not important, can potentially be removed

**Action items:**
```
VIP Score    Interpretation              Action
─────────────────────────────────────────────────
> 2.0        Extremely important         Priority biomarker
1.5-2.0      Very important              Validate experimentally
1.0-1.5      Important                   Include in model
0.8-1.0      Moderately important        Consider for follow-up
< 0.8        Not important               Can exclude
```

#### 4. **04_plsda_loadings_plot.png**
**What it shows:**
- Each point = one m/z feature
- X-axis = Component 1 loading
- Y-axis = Component 2 loading
- Red points = VIP > 1 (important features)
- Labels = Top 10 most important features

**How to interpret:**
- **Distance from origin**: Importance (farther = more important)
- **Direction**: Which component the feature contributes to
- **Clustering**: Groups of related features
- **Quadrants**: Features that vary together

**Biological insight:**
```
Top-right quadrant: Features high in Treatment1
Bottom-left quadrant: Features high in Control
Top-left quadrant: Features high in Treatment2 (if 3 classes)
```

#### 5. **05_plsda_cv_results.png**
**What it shows:**
- X-axis = Number of components used
- Y-axis = Cross-validation accuracy
- Blue line = CV accuracy trend
- Red point = Optimal number of components

**How to interpret:**
- **Rising then plateau**: Good (model captures signal, then stabilizes)
- **Rising then dropping**: Overfitting (too many components)
- **Flat line**: No useful signal (or preprocessing issue)

**Optimal components:**
- Usually 2-5 components
- More is not always better
- Balance accuracy vs. interpretability

#### 6. **06_plsda_confusion_matrix.png**
**What it shows:**
- Heatmap of predictions vs. actual classes
- Rows = Predicted classes
- Columns = Actual classes
- Numbers = Count of samples
- Darker blue = More samples

**How to interpret:**
```
Perfect model:     Problem model:
  C  T1 T2           C  T1 T2
C 15  0  0        C  10  3  2
T1 0 15  0        T1  2 10  3
T2 0  0 15        T2  3  2 10

Diagonal = correct     Off-diagonal = errors
```

**Metrics from confusion matrix:**
- **Accuracy** = Sum of diagonal / Total
- **Sensitivity** = Diagonal / Column sum (per class)
- **Precision** = Diagonal / Row sum (per class)

#### 7. **07_plsda_variance_explained.png**
**What it shows:**
- Blue bars = Variance explained by each component
- Red line = Cumulative variance explained

**How to interpret:**
- **Component 1**: Usually explains most (30-50%)
- **Components 2-3**: Add meaningful information (10-20% each)
- **Later components**: Often just noise (<5% each)

**Target:**
- First 2-3 components should explain >60% of variance
- If <40%, may need better preprocessing or more samples

---

### Data Files (CSV)

#### 1. **plsda_predictions.csv**
**Columns:**
- `SampleID`: Sample identifier
- `Actual`: True class label
- `Predicted`: Model's prediction
- `Correct`: TRUE/FALSE (was prediction correct?)
- `Component1_Score`: Score on component 1
- `Component2_Score`: Score on component 2

**Use for:**
- Identifying misclassified samples
- Quality control
- Manual review of predictions
- Integration with sample metadata

**Example:**
```csv
SampleID,Actual,Predicted,Correct,Component1_Score,Component2_Score
Sample_3,Control,Control,TRUE,-2.34,1.56
Sample_56,Treatment1,Treatment1,TRUE,3.21,-0.89
Sample_107,Treatment2,Treatment1,FALSE,2.45,-1.12
```

#### 2. **plsda_vip_scores.csv**
**Columns:**
- `Feature`: m/z ratio name
- `VIP`: Variable Importance in Projection score

**Use for:**
- Biomarker discovery
- Feature selection for targeted analysis
- Pathway analysis input
- Follow-up experimental validation

**How to use:**
```r
# Load VIP scores
vip <- read.csv("plsda_vip_scores.csv")

# Get important features
important <- vip[vip$VIP > 1, ]

# Get top 10
top10 <- head(vip, 10)

# Export for pathway analysis
write.csv(important, "biomarkers_for_pathway_analysis.csv")
```

#### 3. **plsda_model_summary.csv**
**Columns:**
- `Metric`: Performance metric name
- `Value`: Metric value

**Contains:**
- `Optimal_Components`: Best number of components from CV
- `CV_Accuracy`: Cross-validation accuracy (training set)
- `Test_Accuracy`: Held-out test set accuracy
- `R2X`: Variance explained in X (features)
- `R2Y`: Variance explained in Y (classes)
- `Important_Features`: Number of features with VIP > 1

**Interpret values:**
```
Metric              Good         Excellent      Warning
──────────────────────────────────────────────────────
Test_Accuracy       >0.80        >0.90         <0.70
R2X                 >0.50        >0.70         <0.40
R2Y                 >0.80        >0.95         <0.60
Important_Features  5-20         10-30         >50 or <3
```

#### 4. **plsda_loadings.csv**
**Columns:**
- `Feature`: m/z ratio name
- `Component1_Loading`: Loading on component 1
- `Component2_Loading`: Loading on component 2
- `VIP`: Variable importance score

**Use for:**
- Understanding which features drive each component
- Identifying feature patterns
- Biological interpretation
- Creating custom plots

**Interpretation:**
```r
# Features with high positive loading on Comp1
high_comp1 <- loadings[loadings$Component1_Loading > 0.1, ]

# Features with high negative loading (opposite direction)
low_comp1 <- loadings[loadings$Component1_Loading < -0.1, ]

# These features separate groups along component 1
```

---

## Interpreting Results

### Step 1: Check Model Performance
```r
# Read model summary
summary <- read.csv("plsda_model_summary.csv")
print(summary)
```

**What to look for:**
-  Test accuracy > 0.80 (good)
-  Test accuracy > 0.90 (excellent)
-  Test accuracy < 0.70 (needs improvement)
-  Large gap between CV and test accuracy (overfitting)

**If accuracy is low:**
1. Collect more samples (minimum 20 per class)
2. Improve preprocessing (remove outliers)
3. Check if classes are truly separable
4. Consider if labels are correct

### Step 2: Identify Biomarkers
```r
# Read VIP scores
vip <- read.csv("plsda_vip_scores.csv")

# Get important features
biomarkers <- vip[vip$VIP > 1, ]
print(head(biomarkers, 20))
```

**Prioritization:**
```
Priority    VIP Range    Action
────────────────────────────────────────
High        > 2.0        Validate immediately
Medium      1.5-2.0      Include in panel
Low         1.0-1.5      Monitor in future
Exclude     < 1.0        Not discriminative
```

**Next steps:**
1. Literature search for known biomarkers
2. Pathway analysis (KEGG, MetaboAnalyst)
3. Experimental validation (targeted MS/MS)
4. Clinical validation (independent cohort)

### Step 3: Visualize Separation

**Look at scores plot (01_plsda_scores_plot.png):**
```
Good Separation:              Poor Separation:

  T1 ●●●                       ●○●○
      ●●                       ○●○●
         C ○○○                ●○●○●
           ○○                 ○●●○
  T2 ■■■                      ○●■○
     ■■■                      

Classes distinct              Classes overlap
```

**Biological interpretation:**
- Separated clusters → Distinct metabolic/proteomic profiles
- Overlapping clusters → Similar profiles (biological similarity)
- Outliers → Check for:
  - Sample preparation issues
  - Biological outliers (disease subtypes)
  - Labeling errors

### Step 4: Understand Feature Contributions

**Look at loadings plot (04_plsda_loadings_plot.png):**
```
Component 2 ↑
            │
    Q2      │      Q1
  ●    ●   │   ●     ●
     ●     │       ●
──────────┼──────────→ Component 1
     ●     │    ●
  ●    ●   │   ●    ●
    Q3      │      Q4
            │
```

**Quadrant interpretation:**
- **Q1** (top-right): High in one treatment
- **Q2** (top-left): High in another treatment
- **Q3** (bottom-left): Opposite to Q1
- **Q4** (bottom-right): Opposite to Q2

### Step 5: Validate Predictions
```r
# Check predictions
pred <- read.csv("plsda_predictions.csv")

# Find misclassifications
errors <- pred[pred$Correct == FALSE, ]
print(errors)

# Calculate error rate
error_rate <- sum(!pred$Correct) / nrow(pred)
print(paste("Error rate:", round(error_rate * 100, 2), "%"))
```

**Investigate misclassifications:**
- Are they systematic (same classes confused)?
- Are they biological (intermediate phenotypes)?
- Are they technical (batch effects)?

---

## Using Your Own Data

### Data Format Requirements

Your CSV file should look like this:
```csv
SampleID,Class,Batch,mz_200.00,mz_250.00,mz_300.00,...
S001,Control,1,1234.5,2345.6,3456.7,...
S002,Control,1,1345.6,2456.7,3567.8,...
S003,Treatment,2,5678.9,6789.0,7890.1,...
...
```

**Required columns:**
- **Class**: Group labels (e.g., Control, Treatment1, Treatment2)
  - Must be text or factor
  - Minimum 2 classes, maximum recommended: 5-6
  - Should have meaningful names

- **Feature columns**: m/z intensity values
  - Must be numeric
  - Should be named consistently (e.g., mz_XXX.XX)
  - Intensities should be positive

**Optional columns:**
- SampleID: Unique identifier for each sample
- Batch: Batch number (for quality control)
- Any other metadata

### Step-by-Step Integration

#### 1. Load Your Data
```r
# Load your data
my_data <- read.csv("path/to/your/msms_data.csv")

# Check structure
str(my_data)
head(my_data)

# Check dimensions
cat("Samples:", nrow(my_data), "\n")
cat("Columns:", ncol(my_data), "\n")
```

#### 2. Extract Features and Labels
```r
# Option A: If features start with "mz_"
feature_cols <- grep("^mz_", colnames(my_data), value = TRUE)
features <- as.matrix(my_data[, feature_cols])

# Option B: If features are numeric columns after metadata
metadata_cols <- c("SampleID", "Class", "Batch")  # adjust as needed
feature_cols <- setdiff(colnames(my_data), metadata_cols)
features <- as.matrix(my_data[, feature_cols])

# Extract labels
labels <- factor(my_data$Class)

# Verify
cat("Feature matrix:", dim(features), "\n")
cat("Classes:", levels(labels), "\n")
print(table(labels))
```

#### 3. Quality Check
```r
# Check for missing values
missing <- sum(is.na(features))
cat("Missing values:", missing, "\n")

# If missing values exist:
if(missing > 0) {
  # Option 1: Remove features with >20% missing
  missing_prop <- colSums(is.na(features)) / nrow(features)
  features <- features[, missing_prop < 0.2]
  
  # Option 2: Remove samples with >20% missing
  missing_prop <- rowSums(is.na(features)) / ncol(features)
  keep_samples <- missing_prop < 0.2
  features <- features[keep_samples, ]
  labels <- labels[keep_samples]
  
  # Option 3: Impute with median
  for(i in 1:ncol(features)) {
    features[is.na(features[, i]), i] <- median(features[, i], na.rm = TRUE)
  }
}

# Check for zero or negative values
negative <- sum(features <= 0, na.rm = TRUE)
cat("Zero/negative values:", negative, "\n")

# Replace zeros with small value for log transformation
features[features <= 0] <- min(features[features > 0], na.rm = TRUE) / 2

# Check for outliers (optional)
outliers <- apply(features, 1, function(x) {
  sum(abs(scale(x)) > 3, na.rm = TRUE)
})
cat("Samples with potential outliers:", sum(outliers > 5), "\n")
```

#### 4. Preprocess
```r
# Log transformation
features_log <- log2(features + 1)

# Check if transformation worked
par(mfrow = c(1, 2))
hist(features[, 1], main = "Before Log", xlab = "Intensity")
hist(features_log[, 1], main = "After Log", xlab = "Log2 Intensity")

# Mean centering and scaling
features_scaled <- scale(features_log, center = TRUE, scale = TRUE)

# Verify scaling
cat("Mean of first feature:", round(mean(features_scaled[, 1]), 6), "\n")  # Should be ~0
cat("SD of first feature:", round(sd(features_scaled[, 1]), 6), "\n")      # Should be ~1
```

#### 5. Run Analysis

Now use `features_scaled` and `labels` in the PLS-DA script starting from **Section 3: Train-Test Split**.
```r
# Continue with the script from Section 3
# The rest of the code remains the same!
```

### Example: Complete Integration
```r
# ========== YOUR DATA LOADING ==========
my_data <- read.csv("my_msms_data.csv")

# Extract
feature_cols <- grep("^mz_", colnames(my_data), value = TRUE)
features <- as.matrix(my_data[, feature_cols])
labels <- factor(my_data$Class)

# Quality check
features[features <= 0] <- min(features[features > 0]) / 2
features <- features[complete.cases(features), ]
labels <- labels[complete.cases(features)]

# Preprocess
features_log <- log2(features + 1)
features_scaled <- scale(features_log, center = TRUE, scale = TRUE)

# ========== NOW USE THE SCRIPT ==========
# Copy sections 3-9 from the PLS-DA script
# They will work with your data!
```

---

## Parameter Tuning

### 1. Number of Components

**Default**: Automatic selection via CV

**Manual selection:**
```r
# If you want to force a specific number
optimal_ncomp <- 3  # Change this

# Re-run prediction with chosen number
test_pred <- predict_plsda(plsda_model, X_test, ncomp_use = optimal_ncomp)
```

**Guidelines:**
- **2 components**: Easiest to visualize, usually sufficient
- **3-5 components**: Better accuracy, more complex
- **>5 components**: Risk of overfitting

**How to choose:**
1. Look at CV plot (05_plsda_cv_results.png)
2. Choose where accuracy plateaus
3. Prefer fewer components for interpretability

### 2. Train-Test Split Ratio

**Default**: 70-30

**Change it:**
```r
train_prop <- 0.8  # Use 80% for training, 20% for testing
```

**Guidelines:**
- **Small dataset (<100)**: Use 80-20 or cross-validation only
- **Medium dataset (100-500)**: Use 70-30
- **Large dataset (>500)**: Use 60-40 or 50-50

### 3. Preprocessing Options

#### A. Different Transformations
```r
# Option 1: Log10 instead of Log2
features_log <- log10(features + 1)

# Option 2: Square root (for count data)
features_sqrt <- sqrt(features)

# Option 3: No transformation (if already normalized)
features_log <- features
```

#### B. Different Scaling Methods
```r
# Option 1: Pareto scaling (metabolomics standard)
features_scaled <- scale(features_log) / sqrt(apply(features_log, 2, sd))

# Option 2: Range scaling (0-1)
features_scaled <- apply(features_log, 2, function(x) (x - min(x)) / (max(x) - min(x)))

# Option 3: No scaling (already normalized)
features_scaled <- scale(features_log, center = TRUE, scale = FALSE)
```

#### C. Feature Filtering
```r
# Remove low-variance features
feature_var <- apply(features_scaled, 2, var)
features_scaled <- features_scaled[, feature_var > 0.1]

# Remove highly correlated features (>0.95)
cor_matrix <- cor(features_scaled)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.95)  # requires caret package
features_scaled <- features_scaled[, -high_cor]
```

### 4. Cross-Validation Strategy

**Default**: 5-fold CV

**Change it:**
```r
# More folds (more conservative)
cv_folds <- 10

# Leave-one-out CV (very small datasets)
cv_folds <- nrow(X_train)  # LOOCV
```

### 5. VIP Threshold

**Default**: VIP > 1

**Adjust threshold:**
```r
# More stringent (fewer, but stronger biomarkers)
important_features <- vip_df$Feature[vip_df$VIP > 1.5]

# More lenient (more potential biomarkers)
important_features <- vip_df$Feature[vip_df$VIP > 0.8]
```

---

## Best Practices

### 1. Sample Size Guidelines

**Minimum Requirements:**
- **Absolute minimum**: 15 samples per class
- **Recommended**: 30-50 samples per class
- **Ideal**: 100+ samples per class

**For different numbers of classes:**
```
2 classes:  30 samples (15 per class) minimum
3 classes:  60 samples (20 per class) minimum
4 classes:  80 samples (20 per class) minimum
5+ classes: 100+ samples (20+ per class) minimum
```

**With small samples:**
- Use leave-one-out cross-validation
- Report confidence intervals
- Validate on independent cohort
- Use fewer components (2-3 max)

### 2. Feature Selection

**When to filter features:**
-  >500 features with <100 samples
-  Many low-variance features
-  Highly correlated features (>0.95)
-  Features with >20% missing values

**How to filter:**
```r
# Step 1: Remove features with low variance
feature_var <- apply(features_scaled, 2, var)
keep_features <- feature_var > quantile(feature_var, 0.1)  # Keep top 90%

# Step 2: Remove highly correlated
cor_matrix <- cor(features_scaled[, keep_features])
# (use correlation threshold or hierarchical clustering)

# Step 3: Use in model
features_filtered <- features_scaled[, keep_features]
```

** Warning**: Don't use test set information for feature selection!

### 3. Handling Batch Effects

**Check for batch effects:**
```r
# Visual inspection
pca <- prcomp(features_scaled)
plot(pca$x[, 1:2], col = as.numeric(my_data$Batch), pch = 16)
```

**If batch effects exist:**
```r
# Option 1: Include batch as blocking factor
# (requires advanced implementation)

# Option 2: Use ComBat correction
# library(sva)  # Bioconductor package
# features_corrected <- ComBat(dat = t(features_scaled), batch = my_data$Batch)

# Option 3: Model batch as covariate
# (requires modification of PLS-DA algorithm)
```

### 4. Validation Strategy

**Three levels of validation:**

**Level 1: Internal Validation (Minimum)**
```r
# Cross-validation on training set
# Already included in the script!
```

**Level 2: Test Set Validation (Good)**
```r
# Held-out test set
# Already included in the script!
```

**Level 3: External Validation (Best)**
```r
# Independent cohort
# Run your trained model on completely new data
# This is the gold standard for publication
```

**Recommended approach:**
1. Split your data: 70% train, 30% test
2. Optimize model on training set only
3. Evaluate once on test set
4. If possible, validate on independent cohort

### 5. Reproducibility Checklist

 **Set random seed**
```r
set.seed(123)  # Use same seed for reproducible results
```

 **Document preprocessing steps**
```r
# Keep notes on:
# - Transformation used (log2, log10, sqrt)
# - Scaling method (unit variance, pareto, range)
# - Missing value handling
# - Outlier removal criteria
```

 **Save preprocessing parameters**
```r
# Save for applying to new data
preprocessing_params <- list(
  transformation = "log2",
  scaling = "unit_variance",
  center_values = colMeans(features_log),
  scale_values = apply(features_log, 2, sd)
)
saveRDS(preprocessing_params, "preprocessing_params.rds")
```

 **Record model parameters**
```r
# Already saved in plsda_model_summary.csv!
```

### 6. Publication-Ready Results

**What to report:**
1. **Sample size**: Training and test set sizes
2. **Preprocessing**: Exact steps taken
3. **CV strategy**: Number of folds, repetitions
4. **Optimal components**: How selected
5. **Performance metrics**: Accuracy, sensitivity, specificity per class
6. **Important features**: Top 10-20 with VIP scores
7. **Validation**: Independent cohort results if available

**Figures to include:**
1. Scores plot (with 95% ellipses)
2. VIP scores (top 20)
3. Confusion matrix
4. Cross-validation results
5. (Optional) Loadings plot

**Table to include:**
```
Model Performance Summary

Metric                      Value
────────────────────────────────
Samples (train/test)        105/45
Features                    50
Optimal components          3
CV accuracy                 0.96
Test accuracy               0.93
R²X                         0.72
R²Y                         0.94
Important features (VIP>1)  15

Per-class Performance (Test Set)
Class        Sensitivity  Specificity  Precision
─────────────────────────────────────────────────
Control      0.93         0.97         0.93
Treatment1   0.93         0.97         0.93
Treatment2   0.93         0.97         0.93
```

---

## Troubleshooting

### Problem 1: Low Accuracy (<70%)

**Possible causes:**
1. Classes are not separable (biological similarity)
2. Insufficient samples
3. Too much noise in data
4. Incorrect preprocessing

**Solutions:**
```r
# Check 1: Are classes similar?
# Look at scores plot - if completely overlapping, classes may be too similar

# Check 2: Sample size
table(labels)  # Should have at least 15-20 per class

# Check 3: Check for outliers
pca <- prcomp(features_scaled)
plot(pca$x[, 1:2])  # Look for extreme outliers

# Check 4: Try different preprocessing
# See "Parameter Tuning" section above
```

### Problem 2: Overfitting (CV accuracy >> Test accuracy)

**Example**: CV = 0.95, Test = 0.65

**Causes:**
- Too many components
- Too few samples
- Data leakage

**Solutions:**
```r
# Solution 1: Use fewer components
optimal_ncomp <- 2  # Force 2 components

# Solution 2: Collect more data
# No code solution - need more samples!

# Solution 3: Check for data leakage
# Make sure test samples not used in any way during training
```

### Problem 3: Perfect Accuracy (100%)

**This is suspicious!**

**Check for:**
```r
# 1. Data leakage (test samples in training)
intersect(rownames(X_train), rownames(X_test))  # Should be empty

# 2. Information leakage (e.g., SampleID encoded as feature)
# Remove non-feature columns before analysis

# 3. Overfitting (usually shows on test set)
# If test = 100%, you may have data leakage

# 4. Too easy problem (perfect separation)
# Look at scores plot - widely separated?
```

### Problem 4: Error Messages

**Error: "Singular matrix"**
```r
# Cause: Features are perfectly correlated
# Solution: Remove highly correlated features

cor_matrix <- cor(features_scaled)
high_cor_pairs <- which(abs(cor_matrix) > 0.99 & cor_matrix != 1, arr.ind = TRUE)
# Remove one from each pair
```

**Error: "Non-numeric data"**
```r
# Cause: Class labels in feature matrix
# Solution: Check that features are numeric

str(features)  # All should say "num"
features <- as.matrix(features)  # Force to numeric matrix
```

**Error: "Missing values not allowed"**
```r
# Check for NAs
sum(is.na(features))

# Remove or impute
features <- na.omit(features)  # Remove rows with NA

# OR impute with median
for(i in 1:ncol(features)) {
  features[is.na(features[, i]), i] <- median(features[, i], na.rm = TRUE)
}
```

### Problem 5: Unexpected VIP Scores

**All VIP < 1:**
```r
# Possible causes:
# 1. Too many components (diluted importance)
# 2. No discriminative features
# 3. Scaling issue

# Try:
# - Use fewer components
# - Check preprocessing
# - Verify classes are different
```

**One feature VIP >> others:**
```r
# This feature dominates classification
# Check if it's:
# 1. A true strong biomarker (good!)
# 2. An artifact or technical error
# 3. Sample ID leaked into features

# Verify by checking loadings plot
```

### Problem 6: Can't Visualize Results

**Plots not opening:**
```r
# Make sure you're in correct directory
getwd()
setwd("path/to/output/folder")

# Check if files were created
list.files(pattern = "*.png")

# Try opening manually
png_files <- list.files(pattern = "*.png")
# Open in image viewer
```

**Poor quality plots:**
```r
# Increase resolution
png("plot.png", width = 1200, height = 1000, res = 150)  # Higher res
# ... plotting code ...
dev.off()
```

---

## Statistical Background

### PLS-DA Theory

#### What PLS-DA Does:

1. **Finds latent variables** that:
   - Explain variance in X (features)
   - Predict Y (class membership)
   - Maximize covariance between X and Y

2. **Creates components** that are:
   - Linear combinations of original features
   - Orthogonal to each other (uncorrelated)
   - Ordered by importance (Component 1 > Component 2 > ...)

3. **Uses components** for:
   - Dimensionality reduction (1000 features → 2-5 components)
   - Classification (predict class from components)
   - Interpretation (understand which features matter)

#### Mathematical Formulation:

**Goal**: Find weight vector **w** such that:
```
t = Xw  (scores)
maximize: cov(t, y)
subject to: ||w|| = 1
```

Where:
- X = feature matrix (n × p)
- y = class labels (n × 1)
- t = scores (n × 1)
- w = weights (p × 1)

**NIPALS Algorithm** (used in this script):
```
1. Initialize u = y
2. Repeat until convergence:
   a. w = X'u / ||X'u||        (X weights)
   b. t = Xw                   (X scores)
   c. c = Y't / ||Y't||        (Y weights)
   d. u = Yc                   (Y scores)
3. Deflate: X = X - tp', Y = Y - tq'
4. Repeat for next component
```

#### Variable Importance in Projection (VIP):

VIP measures contribution of each feature to the model:
```
VIP_j = sqrt(p × Σ(w²_jk × SSY_k) / Σ(SSY_k))
```

Where:
- p = number of features
- w_jk = weight of feature j in component k
- SSY_k = sum of squares of Y explained by component k

**Interpretation:**
- VIP > 1: Feature more important than average
- VIP ≈ 1: Average importance
- VIP < 1: Less important than average

### Cross-Validation Theory

**Why cross-validate?**
- Prevents overfitting
- Estimates generalization performance
- Selects optimal hyperparameters

**K-fold CV process:**
```
1. Split data into K folds
2. For fold i = 1 to K:
   a. Train on K-1 folds
   b. Test on fold i
   c. Record accuracy
3. Average K accuracies
```

**Leave-One-Out CV (LOOCV):**
- K = n (number of samples)
- Most unbiased
- Computationally expensive
- High variance (one bad sample affects result)

**Stratified CV:**
- Maintains class proportions in each fold
- Recommended for imbalanced data
- Used in this script

### Model Evaluation Metrics

**Confusion Matrix:**
```
                Actual
              C    T1   T2
Predicted C  [TP] [FP] [FP]
         T1  [FN] [TP] [FP]
         T2  [FN] [FN] [TP]
```

**Metrics:**

1. **Accuracy** = (TP + TN) / Total
   - Overall correctness
   - Can be misleading with imbalanced classes

2. **Sensitivity (Recall)** = TP / (TP + FN)
   - How many actual positives were caught
   - Important when missing positives is costly

3. **Specificity** = TN / (TN + FP)
   - How many actual negatives were correctly identified
   - Important when false alarms are costly

4. **Precision** = TP / (TP + FP)
   - How many predicted positives are actually positive
   - Important when acting on predictions is costly

5. **F1-Score** = 2 × (Precision × Recall) / (Precision + Recall)
   - Harmonic mean of precision and recall
   - Balances both metrics

**Multi-class metrics:**
- **Macro-average**: Average across classes (treats all classes equally)
- **Micro-average**: Aggregate TP/FP/FN then calculate (favors large classes)
- **Weighted average**: Weight by class size

### Variance Explained

**R²X** (Variance explained in X):
```
R²X = 1 - (RSS_X / TSS_X)
```
- RSS_X = Residual sum of squares in X after model
- TSS_X = Total sum of squares in X
- Measures how well model captures feature variance

**R²Y** (Variance explained in Y):
```
R²Y = 1 - (RSS_Y / TSS_Y)
```
- Measures how well model predicts class membership
- Should be high (>0.8) for good classification

**Q²** (Predictive ability - from CV):
```
Q² = 1 - (PRESS / TSS)
```
- PRESS = Prediction error sum of squares
- Estimates model's predictive power on new data
- Q² < 0: Model no better than random
- Q² > 0.5: Good predictive ability

---

## References

### Key Papers on PLS-DA

1. **Original PLS method:**
   - Wold, S., Sjöström, M., & Eriksson, L. (2001). "PLS-regression: a basic tool of chemometrics." *Chemometrics and Intelligent Laboratory Systems*, 58(2), 109-130.

2. **PLS-DA for metabolomics:**
   - Barker, M., & Rayens, W. (2003). "Partial least squares for discrimination." *Journal of Chemometrics*, 17(3), 166-173.

3. **Variable importance:**
   - Wold, S., Johansson, E., & Cocchi, M. (1993). "PLS - partial least squares projections to latent structures." *3D QSAR in Drug Design*, 1, 523-550.

4. **Validation strategies:**
   - Westerhuis, J. A., et al. (2008). "Assessment of PLSDA cross validation." *Metabolomics*, 4(1), 81-89.

### MS/MS Data Analysis

5. **Preprocessing:**
   - Broadhurst, D., et al. (2018). "Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays." *Metabolomics*, 14(6), 72.

6. **Biomarker discovery:**
   - Smilde, A. K., et al. (2005). "ANOVA-simultaneous component analysis (ASCA): a new tool for analyzing designed metabolomics data." *Bioinformatics*, 21(13), 3043-3048.

### Online Resources

7. **MetaboAnalyst** (web-based metabolomics analysis):
   - https://www.metaboanalyst.ca/

8. **KEGG Pathway Database**:
   - https://www.genome.jp/kegg/pathway.html

9. **PLS-DA Tutorial**:
   - Eriksson, L., et al. (2013). "Multi- and Megavariate Data Analysis." Umetrics Academy.

### R Packages (if you want to extend)

10. **mixOmics** (PLS-DA with additional features):
```r
    install.packages("BiocManager")
    BiocManager::install("mixOmics")
```

11. **ropls** (PLS-DA optimized for metabolomics):
```r
    BiocManager::install("ropls")
```

12. **caret** (machine learning framework):
```r
    install.packages("caret")
```

---

## Frequently Asked Questions (FAQ)

### Q1: PLS-DA vs PCA - When to use which?

**Use PCA when:**
- Exploratory analysis (don't know groups yet)
- Outlier detection
- Data quality check
- Unsupervised clustering

**Use PLS-DA when:**
- Known groups (classes)
- Want to classify new samples
- Want biomarker discovery
- Want maximum class separation

**Key difference:** PCA ignores class labels; PLS-DA uses them.

### Q2: How many components should I use?

**Guidelines:**
- Start with 2 (easiest to visualize)
- Use CV to optimize (script does this automatically)
- Rarely need more than 5
- More components ≠ better (risk overfitting)

**Signs you have too many:**
- CV accuracy decreases
- Test accuracy much lower than training
- Variance explained per component <5%

### Q3: What if my classes are imbalanced?

**Example:** 100 controls, 20 patients

**Solutions:**
```r
# Option 1: Stratified sampling (already in script)
# Ensures train/test split maintains proportions

# Option 2: Undersample majority class
minority_n <- min(table(labels))
balanced_idx <- c()
for(class in levels(labels)) {
  class_idx <- which(labels == class)
  balanced_idx <- c(balanced_idx, sample(class_idx, minority_n))
}

# Option 3: Oversample minority class (with caution!)
# Not recommended without SMOTE or similar methods

# Option 4: Use weighted cross-validation
# (requires modification of CV function)
```

### Q4: Can PLS-DA handle >3 classes?

**Yes!** PLS-DA works with any number of classes.

**Practical limits:**
- 2-5 classes: Excellent
- 6-10 classes: Good (need more samples)
- >10 classes: Difficult (need many components)

**With many classes:**
- Need more samples per class
- Visualization becomes harder (use first 2-3 components)
- Consider hierarchical approach (group similar classes)

### Q5: What VIP threshold should I use?

**Standard threshold: VIP > 1**

But you can adjust based on goals:
```
Threshold    Use Case
────────────────────────────────────────
VIP > 2.0    Very conservative, top biomarkers only
VIP > 1.5    Standard conservative
VIP > 1.0    Standard (recommended)
VIP > 0.8    Liberal, exploratory analysis
VIP > 0.5    Very liberal, hypothesis generation
```

**Consider:**
- Sample size (smaller n → higher threshold)
- Cost of validation (expensive → higher threshold)
- Prior knowledge (known biomarkers → may have lower VIP)

### Q6: How do I apply my model to completely new data?
```r
# Step 1: Load your saved model components
# (You'll need to save these after training)

# Step 2: Load new data
new_data <- read.csv("new_samples.csv")
new_features <- as.matrix(new_data[, feature_cols])

# Step 3: Apply SAME preprocessing
new_log <- log2(new_features + 1)
new_scaled <- scale(new_log)  # Uses same parameters!

# Step 4: Predict
new_pred <- predict_plsda(plsda_model, new_scaled, ncomp_use = optimal_ncomp)

# Step 5: Extract predictions
new_classes <- colnames(Y_train)[new_pred$class_pred]
```

**Critical:** Use EXACT same preprocessing as training data!

### Q7: My VIP scores are all < 1. Is this bad?

**Not necessarily!**

VIP scores are relative. If all are <1:
- Could mean features contribute equally
- Could mean need more components
- Could mean weak signal overall

**Check:**
1. Is test accuracy still good? (If yes, you're fine!)
2. Are some VIP close to 1? (Those are still useful)
3. Does model explain variance well? (R²X, R²Y)

**Action:**
- Report top 20 features by VIP anyway
- Look at loadings plot for interpretation
- Consider if preprocessing is appropriate

### Q8: Can I use PLS-DA for regression (continuous outcome)?

**Short answer:** Use PLS regression instead.

**Long answer:**
- PLS-DA is specifically for classification (discrete classes)
- For continuous Y (e.g., drug concentration, survival time):
  - Use standard PLS regression
  - Remove dummy matrix creation
  - Predict continuous values directly

**Modification needed in script:**
```r
# Instead of creating dummy matrix:
Y_train <- as.matrix(continuous_outcome_train)

# Prediction returns continuous value:
Y_pred <- plsda_model$T %*% t(plsda_model$Q)
```

### Q9: How do I cite this analysis?

**For the method:**
```
Partial Least Squares Discriminant Analysis (PLS-DA) was performed
following the methodology of Barker & Rayens (2003). Cross-validation
was used to determine the optimal number of components. Variable 
importance in projection (VIP) scores were calculated to identify
discriminative features (VIP > 1).
```

**For the software:**
```
Analysis was performed using R (version X.X.X) with custom
implementation of the NIPALS algorithm for PLS-DA.
```

### Q10: What's the difference between PLS-DA and OPLS-DA?

**PLS-DA:**
- All components contribute to discrimination
- Orthogonal components
- Simpler interpretation

**OPLS-DA (Orthogonal PLS-DA):**
- Separates predictive and orthogonal variation
- Component 1 = discrimination only
- Other components = within-class variation
- Slightly better interpretation
- Same predictive performance

**When to use OPLS-DA:**
- Want cleaner scores plots
- Focus on single discriminant component
- Have strong within-class variation

**Note:** This script implements PLS-DA. For OPLS-DA, use the `ropls` R package.

---

## Glossary

**Biomarker**: A measurable indicator of biological state or condition. In MS/MS, typically an m/z feature that differs between groups.

**Component (Latent Variable)**: Linear combination of original features that captures variance in data. Created by PLS algorithm.

**Cross-Validation (CV)**: Technique to assess model performance by splitting data into training and validation sets multiple times.

**Discriminant Analysis**: Statistical method to classify samples into predefined groups.

**Dummy Matrix**: Binary matrix encoding class membership (1 = belongs to class, 0 = doesn't).

**IDA (Information Dependent Acquisition)**: MS/MS data acquisition mode that selects precursor ions based on initial scan.

**Loadings**: Coefficients describing contribution of each original feature to a component.

**m/z (mass-to-charge ratio)**: Fundamental measurement in mass spectrometry.

**NIPALS (Nonlinear Iterative Partial Least Squares)**: Algorithm for calculating PLS components.

**Overfitting**: Model performs well on training data but poorly on new data due to learning noise.

**PLS (Partial Least Squares)**: Regression method that finds latent variables maximizing covariance between X and Y.

**R²X**: Proportion of variance in feature matrix explained by model.

**R²Y**: Proportion of variance in class labels explained by model.

**Scores**: Coordinates of samples in component space.

**Stratified Sampling**: Sampling that maintains class proportions in train and test sets.

**VIP (Variable Importance in Projection)**: Metric quantifying contribution of each feature to PLS model.

**Weights**: Coefficients used to calculate scores from original features.

---

## Appendix: Advanced Topics

### A. Permutation Testing

**Purpose**: Assess statistical significance of model
```r
# Permutation test
n_perm <- 1000
perm_accuracy <- numeric(n_perm)

for(i in 1:n_perm) {
  # Permute labels
  perm_labels <- sample(y_train)
  
  # Train model
  Y_perm <- create_dummy_matrix(perm_labels)
  perm_model <- pls_nipals(X_train, Y_perm, ncomp = optimal_ncomp)
  
  # Test
  perm_pred <- predict_plsda(perm_model, X_test, ncomp_use = optimal_ncomp)
  perm_classes <- colnames(Y_test)[perm_pred$class_pred]
  perm_accuracy[i] <- mean(perm_classes == y_test)
}

# P-value
p_value <- sum(perm_accuracy >= accuracy) / n_perm
cat("Permutation p-value:", p_value, "\n")

# Plot
hist(perm_accuracy, main = "Permutation Test",
     xlab = "Permuted Accuracy")
abline(v = accuracy, col = "red", lwd = 2)
```

### B. Confidence Intervals via Bootstrap
```r
# Bootstrap confidence intervals for accuracy
n_boot <- 1000
boot_accuracy <- numeric(n_boot)

for(i in 1:n_boot) {
  # Sample with replacement
  boot_idx <- sample(1:nrow(X_test), replace = TRUE)
  
  # Calculate accuracy on bootstrap sample
  boot_pred <- test_pred$class_pred[boot_idx]
  boot_actual <- y_test[boot_idx]
  boot_classes <- colnames(Y_test)[boot_pred]
  boot_accuracy[i] <- mean(boot_classes == as.character(boot_actual))
}

# 95% CI
ci_lower <- quantile(boot_accuracy, 0.025)
ci_upper <- quantile(boot_accuracy, 0.975)

cat("Accuracy: ", round(accuracy, 3), "\n")
cat("95% CI: [", round(ci_lower, 3), ", ", round(ci_upper, 3), "]\n")
```

### C. Feature Selection via VIP
```r
# Iterative feature selection
vip_threshold <- 1.0
features_selected <- features_scaled

repeat {
  # Train model
  model_temp <- pls_nipals(features_selected[train_idx, ], 
                           Y_train, 
                           ncomp = optimal_ncomp)
  
  # Calculate VIP
  vip_temp <- calculate_vip(model_temp, optimal_ncomp)
  
  # Check if any VIP < threshold
  if(all(vip_temp >= vip_threshold)) break
  
  # Remove features with VIP < threshold
  keep_features <- vip_temp >= vip_threshold
  features_selected <- features_selected[, keep_features]
  
  cat("Features remaining:", ncol(features_selected), "\n")
}

# Final model with selected features
final_model <- pls_nipals(features_selected[train_idx, ], 
                          Y_train, 
                          ncomp = optimal_ncomp)
```

### D. Hierarchical PLS-DA

**For complex multi-class problems:**
```r
# Step 1: Group 1 vs Groups 2+3
model_level1 <- # ... binary PLS-DA ...

# Step 2: Group 2 vs Group 3 (within Groups 2+3)
model_level2 <- # ... binary PLS-DA ...

# Hierarchical prediction
# First classify at level 1, then at level 2
```

### E. Integration with Pathway Analysis
```r
# Extract important features
important_mz <- vip_df$Feature[vip_df$VIP > 1]

# Map to metabolite IDs (requires database)
# Then perform pathway enrichment

# Example output format for MetaboAnalyst
output <- data.frame(
  mz = important_mz,
  VIP = vip_df$VIP[vip_df$VIP > 1]
)
write.csv(output, "for_pathway_analysis.csv", row.names = FALSE)
```

---

## Version History

**Version 1.0** (2025-10-30)
- Initial release
- Complete PLS-DA implementation
- Cross-validation
- VIP scores
- 7 visualizations
- 4 data exports
- Comprehensive documentation

---

## License

This code is provided free for research and educational purposes.

**You are free to:**
- Use for academic research
- Modify for your needs
- Share with colleagues

**Please:**
- Cite appropriately in publications
- Share improvements with community
- Report bugs or issues

---

## Contact & Support

**For questions about:**
- **The method**: See references section
- **The code**: Review comments in script
- **Your specific data**: See "Using Your Own Data" section
- **Troubleshooting**: See "Troubleshooting" section

**Recommended resources:**
- R-help mailing list
- Stack Overflow (tag: [r] [pls])
- Metabolomics forums
- Your local bioinformatics core

---

## Acknowledgments

This implementation is based on:
- Original PLS algorithm by Wold et al.
- NIPALS algorithm for efficient computation
- Best practices from metabolomics community
- Feedback from MS/MS data analysts

Special thanks to the R community for creating an excellent platform for statistical analysis.

---

**End of README**

---

*Last updated: 2025-10-30*
*Version: 1.0*
*Compatible with: R 4.0.0+*
