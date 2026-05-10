# Hindgut Microbiome and Growth Traits in FUNAAB-Alpha Chickens

This project investigates the relationship between hindgut microbial composition and key growth performance traits in **FUNAAB-Alpha chickens**.

## Project Overview
The analysis explores how various bacterial genera influence physical growth metrics, using linear regression models and multivariate analysis (PCA).

## Growth Traits Analyzed
- **BW**: Body Weight
- **WL**: Wing Length
- **SL**: Shank Length
- **BL**: Body Length
- **TL**: Thigh Length
- **KL**: Keel Length

## Analysis Pipeline
1. **Data Preprocessing**: Subsetting microbial abundance data and growth metrics.
2. **Regression Modeling**: Identifying significant microbial predictors for each trait.
3. **PCA**: Visualizing microbiome composition differences across Breeds and Sex.
4. **Relative Abundance**: Normalizing data for comparative analysis.

## Requirements
- R (version 4.0+)
- Packages: `dplyr`, `ggplot2`, `ggfortify`
