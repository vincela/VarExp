# Loading data
data("GWAS")
data("COHORT")

# Computation of the genotype correlation matrix (may take some time)
S = getGenoCorMatrix(GWAS$CHR, GWAS$POS, GWAS$A0, "EUR")

# Parameters of the outcome and the exposure in the pooled sample
parsY = calculateContParams(COHORT$pheno1_N, COHORT$pheno1_Mean, COHORT$pheno1_SD)
parsE = calculateExpoParams(COHORT, "pheno1", "expo1")

# Standardizing betas
std_betaG = standardizeBeta(GWAS$MAIN_EFFECT, GWAS$INT_EFFECT, GWAS$FREQ_A0, parsE[1], parsE[2], "G")
std_betaI = standardizeBeta(GWAS$MAIN_EFFECT, GWAS$INT_EFFECT, GWAS$FREQ_A0, parsE[1], parsE[2], "I")

# Estimation of the fraction of variance explained
fracG = calculateVarFrac_v2(std_betaG, std_betaI, S, parsY[2], sum(COHORT$pheno1_N), "G")
fracI = calculateVarFrac_v2(std_betaG, std_betaI, S, parsY[2], sum(COHORT$pheno1_N), "I")
fracJ = calculateVarFrac_v2(std_betaG, std_betaI, S, parsY[2], sum(COHORT$pheno1_N), "J")

print(c(fracG,fracI,fracJ))
