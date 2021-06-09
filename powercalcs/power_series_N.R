# R functions for power calculations of GWAS studies, under varying MAF and beta parameters. 
# Suitable for classical (i.e. single-SNP single-trait) GWAS studies using linear regression models, i.e for quantitative traits.

# Modification of https://github.com/kaustubhad/gwas-power with the following addition: 
# Can also parameterize the Linkage disequillibrium (r2) between the causal SNP and the genotyped SNP

# Uses formulae for power calculation presented in Appendix A of
# Visscher PM, Wray NR, Zhang Q, et al. 10 Years of GWAS Discovery: Biology, Function, and Translation. Am J Hum Genet 2017;101(1):5-22. doi: 10.1016/j.ajhg.2017.06.005.
# Assumes that additional covariates, if any, are uncorrelated to the SNP.
# in such cases, if the SNP is severely stratified, these formulae will not be applicable, and a simulation-based approach might be preferrable.
# Calculates power by calculating the non-centrality parameter (NCP) of the chi-square statistic under the alternative hypothesis.
# NCP is the degree to which the null distribution is wrong. 

# Parameters:
# n		=    Sample size (has to be >= 0)
# qsq	=    Fraction of trait variance explained by the SNP. Must be in range [0,1]. 
# beta	= Effect size of the SNP on the trait, in SD units. Only square of beta is used, so result is symmetric in +ve and -ve beta values. For simplicity, preferable to have all beta to be +ve.
# maf	= Minor allele frequency of the SNP (between 0 and 0.5).
# het	= Heterozygous genotype frequency. Equal to 2*maf*(1-maf) under Hardy-Weinberg Equilibrium. Must be in range [0,1].  
# pval	= P-value threshold for significance. Common values are 5E-8 (for genome-wide significance) or 1E-5 (for suggestive significance). 

# Use as: source("power_calcs.R")

# Tested under R version 3.4.1.

library(reshape2)
library(ggplot2)
is.scalar <- function(v) { is.numeric(v) && length(v)==1; }
# Accessory function: Check whether the input is numeric scalar.
# Used by other functions for power calculation in checking input consistency.
# The user doesn't need to use this function directly.


power_n_hsq <- function(n, qsq, pval = 5E-8) {
  # Calculates power based on a range of values of n and qsq. pval has to be fixed.
  # Parameters as defined above.
  # n and hsq can be vectors.
  # pval a single numeric value. Default is the common genome-wide significance threshold 5E-8.
  # Returns pow, a matrix containing power values, with n along rows and qsq along columns.
  # Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if qsq=1.
  # Example call: pow = power_n_hsq(n = (1:5)*1000, qsq = (1:10)/100, pval=5E-8)
  
  if (missing(n)) { stop("Parameter n not found."); }
  if (missing(qsq)) { stop("Parameter qsq not found."); }
  if ( !( is.vector(n) && is.numeric(n) ) ) { stop("Parameter n not a numeric vector."); }
  if ( !( is.vector(qsq) && is.numeric(qsq) ) ) { stop("Parameter qsq not a numeric vector."); }
  if ( !( is.scalar(pval) ) ) { stop("Parameter pval not a numeric scalar."); }
  if ( any( n<0 | !is.finite(n) ) ) { stop("Parameter n has unacceptable values."); }
  if ( any( qsq<0 | qsq>1 | !is.finite(qsq) ) ) { stop("Parameter qsq has unacceptable values."); }
  if ( pval<=0 | pval>=1 | !is.finite(pval) ) { stop("Parameter pval has unacceptable value."); }
  
  th=qchisq(pval,df=1,lower.tail=F); # Significance threshold for chi-square, corresponding to P-value threshold
  
  pow=matrix(NA,nrow=length(n),ncol=length(qsq)); # Pre-allocalte matrix to store power values
  rownames(pow)=n; colnames(pow)=qsq; # Assign row and column names
  
  for (i in 1:length(n)) {
    for (j in 1:length(qsq)) {
      ncp=n[i]*qsq[j]/(1-qsq[j]); # Calculate NCP parameter
      if (is.finite(ncp) && ncp>=0) { pow[i,j]=pchisq(th,df=1,lower.tail=F,ncp=ncp); } # Calculate power when NCP is finite and >=0
    }}
  
  return(pow); # Return calculated power matrix
}

## 1. Genome-wide significance (alpha = 5e-8) and UK Biobank-size cohort
pow <- power_n_hsq(n = (1:5)*1000, qsq = (1:10)/100, pval=5E-8)
m <- melt(pow)
names(m) <- c("n", "qsq", "power")
ggplot(m, aes(x = n, y = qsq, fill = power, label = round(power, digits = 2))) + geom_tile(aes(fill = power)) + 
  scale_fill_distiller(palette = "YlGnBu") + geom_text(aes(label = round (power, digits = 2))) + 
  labs (x = "sample size", y = "Frac variance explained", title =  bquote(paste(alpha,"=5e-8,",'',sep='')))  + 
  scale_x_continuous (breaks = unique(m$n)) + scale_y_continuous (breaks = unique(m$qsq))
ggsave(plot = last_plot(), filename = '/home/sorokin/projects/power/power_calcs_ppp_2.png',dpi = 330, units = "in", height = 5, width = 10)
                                     