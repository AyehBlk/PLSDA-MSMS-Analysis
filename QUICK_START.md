# Quick Start Guide

## Getting Started with This Repository

### Prerequisites
- R version 4.0.0 or higher
- No additional packages required (uses base R only)

### Running the Analysis

1. **Clone or download this repository**
```bash
git clone https://github.com/[your-username]/PLSDA-MSMS-Analysis.git
cd PLSDA-MSMS-Analysis
```

2. **Run the script**
```r
# Option 1: In R console
source("plsda_analysis.R")

# Option 2: From command line
Rscript plsda_analysis.R
```

3. **Check your outputs**
The script will generate 11 files in your current directory:
- 7 visualization files (PNG format)
- 4 data files (CSV format)

### Expected Runtime
- **Example data**: ~5-10 seconds
- **Your own data**: Depends on sample size and features

### Using Your Own Data

See the comprehensive [README.md](README.md) for detailed instructions on:
- Data format requirements
- Preprocessing options
- Parameter tuning
- Troubleshooting

### Output Files

**Visualizations:**
1. `01_plsda_scores_plot.png` - Sample clustering
2. `02_plsda_test_predictions.png` - Test set validation
3. `03_plsda_vip_scores.png` - Feature importance
4. `04_plsda_loadings_plot.png` - Feature contributions
5. `05_plsda_cv_results.png` - Model optimization
6. `06_plsda_confusion_matrix.png` - Classification accuracy
7. `07_plsda_variance_explained.png` - Component variance

**Data Files:**
1. `plsda_predictions.csv` - Sample predictions and scores
2. `plsda_vip_scores.csv` - All features ranked by importance
3. `plsda_model_summary.csv` - Model performance metrics
4. `plsda_loadings.csv` - Feature loadings and VIP scores

### Next Steps

1. Review the comprehensive [README.md](README.md) for in-depth documentation
2. Modify the script for your own MS/MS data
3. Explore parameter tuning options
4. See troubleshooting section if you encounter issues

### Getting Help

- Check the [Troubleshooting section](README.md#troubleshooting) in the README
- Review the [FAQ](README.md#frequently-asked-questions-faq)
- Open an issue on GitHub for bugs or feature requests

### Citation

If you use this code in your research, please cite:
```
[Your Name]. (2025). PLS-DA Analysis for MS/MS IDA Data. 
GitHub repository: https://github.com/[your-username]/PLSDA-MSMS-Analysis
```

And cite the original PLS-DA methodology:
```
Barker, M., & Rayens, W. (2003). Partial least squares for discrimination. 
Journal of Chemometrics, 17(3), 166-173.
```
