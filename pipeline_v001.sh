#!/bin/bash
# 1. Data preprocessing: In this step, we will prepare the input data by filtering, cleaning, and formatting the data.
# Data preprocessing
# input: genotypes, phenotype, and covariates
# output: cleaned, filtered, and formatted data

# Filter genotypes
plink --file input_genotypes --maf 0.05 --geno 0.05 --hwe 1e-6 --make-bed --out filtered_genotypes

# Filter phenotype and covariates
awk '{print $1,$2,$3,$4,$5,$6}' input_phenotype > filtered_phenotype
awk '{print $1,$2,$3,$4,$5,$6}' input_covariates > filtered_covariates

# Merge genotypes, phenotype, and covariates
plink --bfile filtered_genotypes --pheno filtered_phenotype --covar filtered_covariates --make-bed --out merged_data

# 2. GWAS analysis: In this step, we will perform genome-wide association analysis to identify the genetic variants associated with the phenotype of interest.
# GWAS analysis
# input: merged data
# output: GWAS results

# Run GWAS
plink --bfile merged_data --linear --pheno filtered_phenotype --covar filtered_covariates --out GWAS_results

# 3. Polygenic risk score calculation: In this step, we will calculate the polygenic risk score based on the GWAS results and the weights for each SNP.
# Polygenic risk score calculation
# input: GWAS results, weights
# output: polygenic risk score

# Extract SNPs and weights
awk '{print $2,$5}' GWAS_results.assoc.linear > snps_weights.txt

# Calculate polygenic risk score
python calculate_PRS.py --snps_weights snps_weights.txt --genotypes filtered_genotypes --out polygenic_risk_score.txt
