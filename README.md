# PLS-DA Analysis for MS/MS IDA Data

![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-active-success)
![Contributions](https://img.shields.io/badge/contributions-welcome-brightgreen)

**Complete implementation of Partial Least Squares Discriminant Analysis (PLS-DA) for tandem mass spectrometry data analysis.**

![PLS-DA Overview](https://img.shields.io/badge/Method-PLS--DA-orange) ![Application](https://img.shields.io/badge/Application-Metabolomics-purple) ![Language](https://img.shields.io/badge/Language-R-276DC3)

---

## 🎯 Overview

This repository provides a comprehensive, production-ready implementation of PLS-DA specifically designed for MS/MS (Mass Spectrometry) IDA (Information Dependent Acquisition) data analysis. Perfect for metabolomics, proteomics, and biomarker discovery workflows.

### ✨ Key Features

- ✅ **Complete PLS-DA implementation** using NIPALS algorithm
- ✅ **Zero dependencies** - uses only base R
- ✅ **Cross-validation** for optimal model selection
- ✅ **Variable Importance (VIP)** scores for biomarker discovery
- ✅ **7 publication-ready visualizations**
- ✅ **4 detailed CSV exports** for further analysis
- ✅ **Comprehensive documentation** with examples
- ✅ **Example dataset included** for immediate testing

---

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/[your-username]/PLSDA-MSMS-Analysis.git
cd PLSDA-MSMS-Analysis
```

### Run the Example

```r
# In R console
source("plsda_analysis.R")

# Or from command line
Rscript plsda_analysis.R
```

**That's it!** The script will:
1. Generate example MS/MS data (150 samples, 50 features, 3 classes)
2. Perform PLS-DA analysis with cross-validation
3. Create 7 plots and 4 data files
4. Display comprehensive results summary

⏱️ **Runtime:** ~5-10 seconds for example data

---

## 📊 Example Output

The analysis automatically generates:

### Visualizations (PNG)
1. **Scores Plot** - Sample clustering and group separation
2. **Test Predictions** - Model validation on held-out data
3. **VIP Scores** - Top 20 most important features
4. **Loadings Plot** - Feature contribution to separation
5. **CV Results** - Model optimization curve
6. **Confusion Matrix** - Classification performance heatmap
7. **Variance Explained** - Component importance

### Data Files (CSV)
1. **Predictions** - Sample classifications with scores
2. **VIP Scores** - All features ranked by importance
3. **Model Summary** - Performance metrics (accuracy, R²X, R²Y)
4. **Loadings** - Feature contributions to each component

---

## 📖 Documentation

### For Beginners
- **[QUICK_START.md](QUICK_START.md)** - Get running in 5 minutes
- **Sections 1-5 of README_FULL.md** - Basic concepts and usage

### For Advanced Users
- **[README_FULL.md](README_FULL.md)** - Complete 2000+ line documentation including:
  - Statistical theory and mathematical background
  - Using your own data (step-by-step guide)
  - Parameter tuning and optimization
  - Troubleshooting common issues
  - Best practices for publication
  - Advanced topics (permutation tests, bootstrap, etc.)

---

## 🔬 Use Cases

Perfect for:
- **Metabolomics** - Identify metabolite biomarkers between groups
- **Proteomics** - Discover differentially expressed proteins
- **Lipidomics** - Classify samples based on lipid profiles
- **Biomarker Discovery** - Rank features by discriminative power
- **Quality Control** - Validate analytical methods
- **Method Development** - Optimize sample preparation protocols

---

## 💡 Why This Implementation?

### Advantages of PLS-DA for MS/MS Data

| Feature | Benefit |
|---------|---------|
| **Handles high dimensionality** | Works with 1000s of m/z features |
| **Manages multicollinearity** | Correlated features (common in MS data) |
| **Supervised classification** | Uses group labels for maximum separation |
| **Interpretable results** | VIP scores, loadings, and scores plots |
| **Robust to noise** | Dimensionality reduction filters noise |
| **No external dependencies** | Pure R implementation |

### Why Not Just Use PCA?

- **PCA** is unsupervised (ignores group labels)
- **PLS-DA** maximizes separation between known groups
- **PLS-DA** provides biomarker rankings (VIP scores)
- **PLS-DA** is designed for classification tasks

---

## 📋 Requirements

- **R** version ≥ 4.0.0
- **No additional packages required!**
- **Memory:** Minimum 4GB RAM (8GB recommended)
- **Storage:** ~10MB for outputs

---

## 🎓 Methodology

This implementation uses:
- **NIPALS algorithm** for PLS component extraction
- **Stratified train-test split** (70-30 default)
- **K-fold cross-validation** for component optimization
- **VIP scores** for feature importance ranking
- **Dummy matrix encoding** for multi-class problems

See [README_FULL.md](README_FULL.md) Section 13 for complete statistical background.

---

## 📝 Using Your Own Data

Your CSV should have this structure:

```csv
SampleID,Class,Batch,mz_200.00,mz_250.00,mz_300.00,...
S001,Control,1,1234.5,2345.6,3456.7,...
S002,Treatment,1,5678.9,6789.0,7890.1,...
```

See [README_FULL.md](README_FULL.md) Section 9 for detailed integration guide.

**Quick integration:**
```r
# Load your data
my_data <- read.csv("your_msms_data.csv")

# Extract features and labels
feature_cols <- grep("^mz_", colnames(my_data), value = TRUE)
features <- as.matrix(my_data[, feature_cols])
labels <- factor(my_data$Class)

# Continue with script from Section 2 (Preprocessing)
```

---

## 🤝 Contributing

Contributions are welcome! Areas for improvement:
- Additional preprocessing methods
- Support for more data formats
- Integration with pathway analysis tools
- Additional visualization options
- Performance optimizations

Please open an issue or submit a pull request.

---

## 📚 Citation

If you use this code in your research, please cite:

**This repository:**
```
AyebBlk. (2025). PLS-DA Analysis for MS/MS IDA Data. 
GitHub: https://github.com/[your-username]/PLSDA-MSMS-Analysis
```

**Original PLS-DA method:**
```
Barker, M., & Rayens, W. (2003). Partial least squares for discrimination. 
Journal of Chemometrics, 17(3), 166-173.
```

**NIPALS algorithm:**
```
Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: 
a basic tool of chemometrics. Chemometrics and Intelligent Laboratory 
Systems, 58(2), 109-130.
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**You are free to:**
- ✅ Use for academic research
- ✅ Use for commercial projects
- ✅ Modify and distribute
- ✅ Include in your own projects

**Just include the license and attribution!**

---

## 🐛 Troubleshooting

### Common Issues

**Problem:** Low accuracy (<70%)
- Check if classes are actually separable (biological question)
- Ensure sufficient samples (minimum 20 per class recommended)
- Try different preprocessing methods

**Problem:** Error "Singular matrix"
- Features are too highly correlated
- Remove features with correlation >0.95

**Problem:** All VIP scores <1
- May need fewer components
- Check if preprocessing is appropriate
- Verify classes are actually different

See [README_FULL.md](README_FULL.md) Section 12 for complete troubleshooting guide.

---

## 📞 Support

- 📖 Check the [comprehensive documentation](README_FULL.md)
- ❓ Read the [FAQ](README_FULL.md#frequently-asked-questions-faq)
- 🐛 Open an [issue](https://github.com/[your-username]/PLSDA-MSMS-Analysis/issues) for bugs
- 💬 Start a [discussion](https://github.com/[your-username]/PLSDA-MSMS-Analysis/discussions) for questions

---

## 🌟 Star History

If you find this useful, please consider giving it a star! ⭐

---

## 🔗 Related Resources

- [MetaboAnalyst](https://www.metaboanalyst.ca/) - Web-based metabolomics analysis
- [KEGG](https://www.genome.jp/kegg/) - Pathway database
- [mixOmics R package](http://mixomics.org/) - Extended multivariate methods
- [ropls Bioconductor](https://bioconductor.org/packages/ropls/) - Alternative PLS-DA implementation

---

## 👤 Author

**Ayeh Bolouki**
- GitHub: [@AyehBlk](https://github.com/your-username)
- Role: Computational Biologist / Bioinformatician

---

## 🎯 Project Status

✅ **Active Development** - Maintained and open to contributions

**Current Version:** 1.0
**Last Updated:** October 2025

---

<p align="center">
 Made with ❤️ - Let's make free science for everybody around the world.
</p>

<p align="center">
  <sub>If this helped your research, consider citing it in your publications!</sub>
</p>
