library(ukbb)
library(tidyverse)
library(data.table)
library(ggplot2)

rm(list=ls())

##Setup the UK Biobank version.
setup_ukbb("ukbb_v4")

ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
       return(y)
}

args = commandArgs(trailingOnly = TRUE)

if (length(args > 0)) {
  phenotype = args[1]
  fam_file = args[2]
  outdir = args[3]
} else {
  phenotype = "aac_ml_mean_20200720"
  fam_file = "/data/ukbb/uk-biobank/aging_18448/genotypes/raw/v2/exomes/ukb23155_c1_b0_v1_s200632.fam"
  outdir = "~/ukbb_project/aortic_calcif/model_data/aging_v5/pheno"
}

fam = fread(fam_file, header = FALSE)
setkey(fam, 'V2')
message("Size of exome cohort: ", nrow(fam))

##Get European ancestry indivs n=409624
eur = get_field(22006) %>%
  select(individual_id) 
eur = data.table(eur)
setkey(eur, "individual_id")
fam = fam[fam$V2 %in% eur$individual_id,]
message("Size of European exome cohort: ", nrow(fam))

genetic_sex = get_named_field(22001, "genetic_sex") %>%
  mutate(genetic_sex = 2 - as.numeric(genetic_sex)) %>%
  filter(individual_id %in% fam$V2)
message("Number of nonmissing genetic sex values:", nrow(genetic_sex))

##Age at the relevant assessment
#if (phenotype_class %in% c("imaging", "brain")) {
#  age = get_phenotype("imaging_age", "covariates") %>%
#    rename(age=imaging_age)
#  center = get_phenotype("imaging_center", "covariates") %>%
#    rename(center=imaging_center)
#} else if (phenotype_class == "baseline") {
#  age = get_phenotype("age", "covariates")
#  center = get_phenotype("baseline_center", "covariates") %>%
#    rename(center=baseline_center)
#} 

## Baseline age 
age = get_phenotype("age", subdir="covariates") %>%
  filter(individual_id %in% fam$V2)
message("Number of nonmissing baseline age values:", nrow(age))

## Age at the relevant assessment ##TODO: Not available for DEXA scan participants
age = get_phenotype("imaging_age", subdir="covariates") %>%
  rename(age=imaging_age) %>%
  filter(individual_id %in% fam$V2)
message("Number of nonmissing imaging age values:", nrow(age))


##PCs needs its own processing step because its stored in melted table
pcs_raw  = get_field(22009)
pcs = pcs_raw[c('individual_id', 'array_id', 'value')] %>% 
  pivot_wider(names_from='array_id', values_from='value')
names(pcs) = c('individual_id',paste("PC", seq(1:40), sep = "_"))
pcs = pcs %>%
  filter(individual_id %in% fam$V2)
message("Number of nonmissing PCs:", nrow(pcs))

covar = age %>%
  #inner_join(center, by="individual_id") %>% #Only for a brain imaging subset.. 
  inner_join(genetic_sex, by = "individual_id") %>%
  inner_join(pcs, by="individual_id") 
message("Size of covariate dataframe: ", nrow(covar))

# ##Quantitative phenotypes
# aac  = get_phenotype("aac_ml_mean_20200722") 

# covar = covar %>% 
#   inner_join(aac,by='individual_id') 
# dim(covar)

dat = tibble(data.frame(individual_id = covar$individual_id))
p = "aac_ml_mean_20200722"
temp = get_phenotype(p)
standardized_name = sprintf("%s.std", p)
rank_norm_name = sprintf("%s.rank", p)
log_name = sprintf("%s.log", p)
ihs_name = sprintf("%s.ihs", p)
resid_name = sprintf("%s.rnres", p )

temp = temp %>%
filter (individual_id %in% covar$individual_id)
dat = temp %>%
  mutate (standardized = scale(temp[[p]])) %>%
  mutate (ranknormed = ukbb::rank_norm(temp[[p]])) %>%
  mutate (log_transfo = log(temp[[p]])) %>%
  mutate (ihs_transfo = ihs(temp[[p]])) %>%
  rename(!!standardized_name := standardized) %>%
  rename(!!rank_norm_name := ranknormed) %>%
  rename(!!log_name := log_transfo) %>%
  rename(!!ihs_name := ihs_transfo) %>%
  inner_join(covar, by = 'individual_id')
# Plotting code
message("Phenotype: ", p)
message("Size of phenotypes dataframe: ", nrow(unique(dat)))
ggplot(dat, aes_string(x = p)) + geom_histogram() + theme_classic()
ggsave(plot=last_plot(),
  filename = sprintf("%s/%s.png", outdir, p), width = 5, height = 5, dpi = 330)
ggplot(dat, aes_string(x = standardized_name)) + geom_histogram() + theme_classic()
ggsave(plot=last_plot(),
  filename = sprintf("%s/%s.std.png", outdir, p), width = 5, height = 5, dpi = 330)
ggplot(dat, aes_string(x = rank_norm_name)) + geom_histogram() + theme_classic()
ggsave(plot=last_plot(),
  filename = sprintf("%s/%s.rank.png", outdir, p), width = 5, height = 5, dpi = 330)
ggplot(dat, aes_string(x = log_name)) + geom_histogram() + theme_classic()
ggsave(plot=last_plot(),
  filename = sprintf("%s/%s.log.png", outdir, p), width = 5, height = 5, dpi = 330)
ggplot(dat, aes_string(x = ihs_name)) + geom_histogram() + theme_classic()
ggsave(plot=last_plot(),
  filename = sprintf("%s/%s.ihs.png", outdir, p), width = 5, height = 5, dpi = 330)
dat = dat %>%
mutate (resid_transfo = ukbb::rank_norm(resid(lm( data = dat, formula = sprintf ("%s ~ age + genetic_sex + PC_1 + PC_2 + PC_3 + PC_4
  + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10", p), family = "linear")))) %>%
rename(!!resid_name:= resid_transfo)
ggplot(dat, aes_string(x = resid_name)) + geom_histogram() + theme_classic()
ggsave(plot=last_plot(), filename = sprintf("%s/%s.rnres.png", outdir, p), width = 5, height = 5, dpi = 330)

write.table(dat, sprintf("%s/%s.tr2.txt", outdir, p), sep = '\t', quote = F, row.names = F)
write.table(dat[,c(1,1)], sprintf("%s/%s.tr2.idvs", outdir, p), sep = '\t', quote = F, row.names = F)
