# Quality control for SNPs and sample 
# Written by Lucia Trastulla -- email: lucia_trastulla@psych.mpg.de

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(doParallel))

parser <- ArgumentParser(description = "Quality control of genomestudio analysis for samples and SNPs" )

parser$add_argument("--SNP_table_file", type = "character", help = ".txt file, SNP table  from genomestudio")
parser$add_argument("--Full_table_file", type = "character", help = ".txt file, Full data table from genomestudio")
parser$add_argument("--Sample_table_file", type = "character", help = ".txt file, sample table from genomestudio")
parser$add_argument("--PAR_file", type = "character", help = ".txt file, PAR region for haploid chr depends on the chr37 reference")
parser$add_argument("--fold_fun", type = "character", help = "folder with the functions to be loaded")
parser$add_argument("--clust_sep_as", type = "double", default = 0.3, help = "cluster separation parameter for autosomal snps")
# parser$add_argument("--call_freq_as", type = "double", default = 0.95, help = "call frequency parameter for autosomal snps")
parser$add_argument("--ABR_mean_as", type = "double", default = 0.2, help = "AB R mean paramter for autosomal SNPs")
parser$add_argument("--AAR_mean_as", type = "double", default = 0.2, help = "AA R mean paramter for autosomal SNPs")
parser$add_argument("--BBR_mean_as", type = "double", default = 0.2, help = "BB R mean paramter for autosomal SNPs")
parser$add_argument("--ABT_mean_low_as", type = "double", default = 0.1, help = "AB T mean paramter for autosomal SNPs (lower bound)")
parser$add_argument("--ABT_mean_up_as", type = "double", default = 0.9, help = "AB T mean paramter for autosomal SNPs (upper bound)")
parser$add_argument("--het_low_as", type = "double", default = -0.9, help = "Het excess paramter for autosomal SNPs (lower bound)")
parser$add_argument("--het_up_as", type = "double", default = 0.9, help = "Het excess paramter for autosomal SNPs (upper bound)")
parser$add_argument("--MAF_as", type = "double", default = 0, help = "MAF parameter for autosomal SNPs")
parser$add_argument("--AAT_mean_as", type = "double", default = 0.3, help = "AA T mean paramter for autosomal SNPs")
parser$add_argument("--AAT_dev_as", type = "double", default = 0.06, help = "AA T dev paramter for autosomal SNPs")
parser$add_argument("--BBT_mean_as", type = "double", default = 0.7, help = "BB T mean paramter for autosomal SNPs")
parser$add_argument("--BBT_dev_as", type = "double", default = 0.06, help = "BB T dev paramter for autosomal SNPs")
#parser$add_argument("--ABfreq_hp", type = "double", default = 0.2, help = "AB freq paramter for X (males) SNPs")
parser$add_argument("--male_frac", type = "double", default = 0.5, help = "percentage of males samples with GT AB and ABfreq>ABfreq_hp")
parser$add_argument("--R_hpY", type = "double", default = 0.2, help = "R parameter for Y (females) SNPs")
parser$add_argument("--female_frac", type = "double", default = 0.5, help = "percentage of females samples with R>R_hpY and call freq > call_freq_hpY")
# parser$add_argument("--call_freq_hpY_male", type = "double", default = 0.95, help = "call frequency parameter for Y snps in males")
parser$add_argument("--manifest_name", type="character", help = "name of the manifest file used") # GSA.24v1.0_A6
parser$add_argument("--ncores", type="integer", help = "n cores for the parallelization")
parser$add_argument("--outf", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
#############################


# tables
fold_fun <- args$fold_fun
Full_table_file <- args$Full_table_file
Sample_table_file <- args$Sample_table_file
SNP_table_file <- args$SNP_table_file
PAR_file <- args$PAR_file
outFile <- args$outf

# parameters
# autosomal chr
clust_sep_as <- args$clust_sep_as
# call_freq_as <- args$call_freq_as
ABR_mean_as <- args$ABR_mean_as
AAR_mean_as <- args$AAR_mean_as
BBR_mean_as <- args$BBR_mean_as
ABT_mean_low_as <- args$ABT_mean_low_as
ABT_mean_up_as <- args$ABT_mean_up_as
het_low_as <- args$het_low_as
het_up_as <- args$het_up_as
MAF_as <- args$MAF_as
AAT_mean_as <- args$AAT_mean_as
AAT_dev_as <- args$AAT_dev_as
BBT_mean_as <- args$BBT_mean_as
BBT_dev_as <- args$BBT_dev_as

# haploid chr
#ABfreq_hp <- args$ABfreq_hp
male_frac <- args$male_frac
R_hpY <- args$R_hpY
#call_freq_hpY <- args$call_freq_hpY
female_frac <- args$female_frac
# call_freq_hpY_male <- args$call_freq_hpY_male
manifest_name <- args$manifest_name
ncores <- args$ncores

#################
# load the function
source(paste(fold_fun, 'gt_functions.R', sep = ""))
#################

# read the tables
Sample_table <- read.table(Sample_table_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE, quote = "", dec=",")
SNP_table <- read.csv(SNP_table_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE, quote = "", dec=",")
Full_table <- read.table(Full_table_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE, quote = "", dec=",")

# used for haploid SNPs
PAR_table <- read.table(PAR_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE, dec=".", skip = 1)
nsamples <- nrow(Sample_table)

#### produce a txt file which include the cll rate for each sample before and after the quality control
df_sample_callrate <- data.frame(sample_ID = Sample_table$Sample.ID)
df_sample_callrate$call_rate <- Sample_table$Call.Rate
# number od samples with a call rate > 0.98
nsamples_QC <- length(which(Sample_table$Call.Rate > 0.98))

####### removal of snps for autosomal chr ######
# exclude X, Y (not pseudo autosomal regions) MT and 0 chr
id_chrX <- which(SNP_table$Chr == 'X')
id_chrY <- which(SNP_table$Chr == 'Y')
id_chrMT <- which(SNP_table$Chr == 'MT')
id_chr0 <- which(SNP_table$Chr == 0)
id_chrXY <- which(SNP_table$Chr == 'XY')

# exclude only not pseudo autosomal region
SNP_table_X <- SNP_table[id_chrX,]
chrX_SNPname <- SNP_table_X$Name

SNP_table_Y <- SNP_table[id_chrY,]
chrY_SNPname <- SNP_table_Y$Name

chrX_notPAR_SNPname <- PAR_table$V1[which(PAR_table$V1 %in% chrX_SNPname & PAR_table$V2 == 0)]
chrY_notPAR_SNPname <- PAR_table$V1[which(PAR_table$V1 %in% chrY_SNPname & PAR_table$V2 == 0)]

id_chrX_notPAR <- which(SNP_table$Name %in% chrX_notPAR_SNPname)
id_chrY_notPAR <- which(SNP_table$Name %in% chrY_notPAR_SNPname)

# exclude chr0 and MT and XY from the further analysis
# chrXY pseudo-autosomal regions
SNP_table_as <- SNP_table[-c(id_chrX_notPAR, id_chrY_notPAR, id_chrMT, id_chr0),]
# SNP_table_as <- SNP_table[-c(id_chrX_notPAR, id_chrY_notPAR, id_chrMT, id_chr0, id_chrXY),]
# SNP_table_as <- SNP_table[-c(id_chrX, id_chrY, id_chrMT, id_chr0),]

# create a df explaining how many snps are excluded for each step
step_name <- c('initial_nSNPs', 'clust_sep', 'call_freq', 'AB_R_mean', 'AA_R_mean', 'BB_R_mean', 'AB_T_mean', 'het', 'MAF', 'AA_T_mean', 
'AA_T_dev', 'BB_T_mean', 'BB_T_dev', 'tot_unique_as', 'ABfreqX_male', 'callfreqY_female', 'callfreqY_male', 'tot_unique_hp', 'chr_0_MT',
'tot', 'perc')
#step_name <- c('initial_nSNPs', 'clust_sep', 'call_freq', 'AB_R_mean', 'AA_R_mean', 'BB_R_mean', 'AB_T_mean', 'het', 'MAF', 'AA_T_mean',  
#'AA_T_dev', 'BB_T_mean', 'BB_T_dev', 'tot_unique_as', 'ABfreqX_male', 'callfreqY_female', 'callfreqY_male', 'tot_unique_hp', 'chr_0_MT_XY', 
#'tot', 'perc')

df_excl_SNP <- data.frame(step = step_name)
vect_lensnps <- nrow(SNP_table)

# cluster sep
id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.Cluster.Sep'))
id_clustsep <- which(SNP_table_as[,id_name] <= clust_sep_as)
vect_lensnps <- c(vect_lensnps, length(id_clustsep))

# call freq
# to compute the call freq use only the number of samples that pass the quality control
call_freq_as <- thr_SNPscallrate_fun(nsamples_QC)
id_callfreq <- which(SNP_table_as$Call.Freq <= call_freq_as)
vect_lensnps <- c(vect_lensnps, length(id_callfreq))
# AB R mean
id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.AB.R.Mean'))
id_ABRmean <- which(SNP_table_as[,id_name] <= ABR_mean_as)
vect_lensnps <- c(vect_lensnps, length(id_ABRmean))
# AA R mean
id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.AA.R.Mean'))
id_AARmean <- which(SNP_table_as[,id_name] <= AAR_mean_as)
vect_lensnps <- c(vect_lensnps, length(id_AARmean))
# BB R mean
id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.BB.R.Mean'))
id_BBRmean <- which(SNP_table_as[,id_name] <= BBR_mean_as)
vect_lensnps <- c(vect_lensnps, length(id_BBRmean))
# AB T mean
id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.AB.T.Mean'))
id_ABTmean <- which(SNP_table_as[,id_name] < ABT_mean_low_as | SNP_table_as[,id_name] > ABT_mean_up_as)
vect_lensnps <- c(vect_lensnps, length(id_ABTmean))
# het excess
id_het <- which(SNP_table_as$Het.Excess < het_low_as | SNP_table_as$Het.Excess > het_up_as)
vect_lensnps <- c(vect_lensnps, length(id_het))
# MAF
id_MAF <- which(SNP_table_as$AB.Freq == 0 & SNP_table_as$Minor.Freq > MAF_as)
vect_lensnps <- c(vect_lensnps, length(id_MAF))
# AA and BB freq
id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.AA.T.Mean'))
id_AATmean <- which(SNP_table_as$AA.Freq == 1 & SNP_table_as[, id_name] > AAT_mean_as)
vect_lensnps <- c(vect_lensnps, length(id_AATmean))

id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.AA.T.Dev'))
id_AATdev <- which(SNP_table_as$AA.Freq == 1 & SNP_table_as[,id_name] > AAT_dev_as)
vect_lensnps <- c(vect_lensnps, length(id_AATdev))

id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.BB.T.Mean'))
id_BBTmean <- which(SNP_table_as$BB.Freq == 1 & SNP_table_as[,id_name] < BBT_mean_as)
vect_lensnps <- c(vect_lensnps, length(id_BBTmean))

id_name <-  which(colnames(SNP_table_as) == paste0(manifest_name, '.bpm.BB.T.Dev'))
id_BBTdev <- which(SNP_table_as$BB.Freq == 1 & SNP_table_as[,id_name] > BBT_dev_as)
vect_lensnps <- c(vect_lensnps, length(id_BBTdev))

# delete all the SNPs selected
id_SNPs_as <- unique(c(id_clustsep,id_callfreq, id_ABRmean, id_AARmean, id_BBRmean, id_ABTmean, id_het, id_MAF, id_AATmean, id_AATdev, id_BBTmean, id_BBTdev))
vect_lensnps <- c(vect_lensnps, length(id_SNPs_as))
# SNP_table_as_filt <- SNP_table_as[-id_SNPs_as,]



#### infer sex if not known
# theta, R, LogR, BAF columns are in the wrong format (charcater, transform "," with ".")
# Full_table contains Gtype, score, Theta, R, BAF, LRR
# Full_table_old <- Full_table

id_column <- as.vector(sapply(0:(nsamples-1), function(x) (13:16)+(6*x)))

for(i in 1:length(id_column)){
  
  Full_table[,id_column[i]] <- as.numeric(gsub(",", ".", gsub("\\.", "", Full_table[,id_column[i]])))
  
}
# NOTE: na produced where n.def is present in the original Full Data Table

# check the chr X, if the majority of Gtype if A or B then is a male
Full_table_X <- Full_table[which(Full_table$Chr == 'X'),]
id_column_GT <- as.vector(sapply(0:(nsamples-1), function(x) 11 + (6*x)))

# summarize
GT_chrX <- sapply(id_column_GT, function(x) table(Full_table_X[,x]))

# if the structure is saved in list form, change the format in a matrix
if(is.list(GT_chrX)){
  
  temp <- matrix(nrow = 4, ncol = nsamples, 0)
  GT_type <- c('AA', 'AB', 'BB', 'NC')
  rownames(temp) <- GT_type
  
  for(i in 1:nsamples){
  
        id <- which(GT_type %in% names(GT_chrX[[i]]) )
        temp[id,i] <- GT_chrX[[i]]
      
      }
 
  GT_chrX <- temp   
}
  

# compute percentage of AB
colnames(GT_chrX) <- colnames(Full_table_X[,id_column_GT])
AB_chrX_perc <- apply(GT_chrX,2, function(x) x[2]/nrow(Full_table_X))

# male and female detected from X
id_male_X <- which(AB_chrX_perc <= 0.01)
id_female_X <- which(AB_chrX_perc > 0.01)


# check also the Y chromosome
Full_table_Y <- Full_table[which(Full_table$Chr == 'Y'),]
# summarize
GT_chrY <- sapply(id_column_GT, function(x) table(Full_table_Y[,x]))

# if the structure is saved in list form, change the format in a matrix
if(is.list(GT_chrY)){
  
  temp <- matrix(nrow = 4, ncol = nsamples, 0)
  GT_type <- c('AA', 'AB', 'BB', 'NC')
  rownames(temp) <- GT_type
  
  for(i in 1:nsamples){
    
    id <- which(GT_type %in% names(GT_chrY[[i]]) )
    temp[id,i] <- GT_chrY[[i]]
    
  }
  
  GT_chrY <- temp   
}

colnames(GT_chrY) <- colnames(Full_table_Y[,id_column_GT])
NC_chrY_perc <- apply(GT_chrY,2, function(x) x[4]/nrow(Full_table_Y))

# male and female detected from Y
id_male_Y <- which(NC_chrY_perc < 0.7)
id_female_Y <- which(NC_chrY_perc >= 0.7)

id_male <- list(id_male_X, id_male_Y)
id_female <- list(id_female_X, id_female_Y)

# three 'type' of gender: male, female and undefined (if the info from X and Y don't match)
# id_big <- which.max(sapply(id_male, length))
# id_small <- setdiff(1:2,id_big)

id_male <- intersect(id_male[[1]], id_male[[2]])
id_female <- intersect(id_female[[1]], id_female[[2]])
id_undef <- setdiff(1:nsamples, c(id_male, id_female))

# sample name from the full data table 
sample_name_FT <- lapply(colnames(GT_chrX), function(x) strsplit(x, "")[[1]])
sample_name_FT <- lapply(sample_name_FT, function(x) x[-c(1,(length(x)-6+1):length(x))])
sample_name_FT <- sapply(sample_name_FT, function(x) paste(x, sep="", collapse=""))

# add the gender in the loaded Sample_table
if(length(id_undef)>0){
  id_undef_ST <- sapply(sample_name_FT[id_undef], function(x) which(x == Sample_table$Sample.ID))
  Sample_table$Gender[id_undef_ST] <- "Undef"
}
if(length(id_male)>0){
  id_male_ST <- sapply(sample_name_FT[id_male], function(x) which(x == Sample_table$Sample.ID))
  Sample_table$Gender[id_male_ST] <- "Male"
}
if(length(id_female)>0){
  id_female_ST <- sapply(sample_name_FT[id_female], function(x) which(x == Sample_table$Sample.ID))
  Sample_table$Gender[id_female_ST] <- "Female"
}

# update the new table
df_sample_callrate$gender <- Sample_table$Gender


#### removal of SNPs for the Haploid chr
# use the PAR file
#### chr X
nmale <- length(id_male)
nfemale <- length(id_female)
nundef <- length(id_undef)


registerDoParallel(cores = ncores)

bool_SNP_to_exclude <- foreach(i=1:length(chrX_notPAR_SNPname))%dopar%{

  print(i)

  id <- which(Full_table$Name == chrX_notPAR_SNPname[i])
  GT_id <- Full_table[id,id_column_GT]
  GT_male <- GT_id[id_male]

  (length(which(GT_male == 'AB'))/nmale) > male_frac

}

id_SNP_to_exclude <- which(unlist(bool_SNP_to_exclude))
name_SNP_toexclude_X <- chrX_notPAR_SNPname[id_SNP_to_exclude]

vect_lensnps <- c(vect_lensnps, length(name_SNP_toexclude_X))


#### chr Y

id_column_R <- as.vector(sapply(0:(nsamples-1), function(x) 14 + (6*x)))

if(length(id_female)>0){
  registerDoParallel(cores = ncores)
  bool_SNP_to_exclude <- foreach(i=1:length(chrY_notPAR_SNPname))%dopar%{
  
    print(i)
  
    id <- which(Full_table$Name == chrY_notPAR_SNPname[i])
    GT_id <- Full_table[id,id_column_GT]
    GT_female <- GT_id[id_female]
    R_id <- Full_table[id,id_column_R]
    R_female <- R_id[id_female]
  
  
  
    (length(which(GT_female != 'NC' | R_female > R_hpY))/nfemale) > female_frac
  
  }
  id_SNP_to_exclude <- which(unlist(bool_SNP_to_exclude))
  name_SNP_toexclude_Y_female <- chrY_notPAR_SNPname[id_SNP_to_exclude]
  
}else{
  name_SNP_toexclude_Y_female <- c()
}


vect_lensnps <- c(vect_lensnps, length(name_SNP_toexclude_Y_female))


SNP_table_Y_notPAR <- SNP_table[which(SNP_table$Name %in% chrY_notPAR_SNPname),]
call_freq_hpY_male <- thr_SNPscallrate_fun(nmale)
name_SNP_toexclude_Y_male <- SNP_table_Y_notPAR$Name[which(SNP_table_Y_notPAR$X..Calls/nmale < call_freq_hpY_male)]
vect_lensnps <- c(vect_lensnps, length(name_SNP_toexclude_Y_male))

#### save the final name for the SNP to be excluded
name_SNP_toexclude_hp <- unique(c(name_SNP_toexclude_X, name_SNP_toexclude_Y_male, name_SNP_toexclude_Y_female))
vect_lensnps <- c(vect_lensnps, length(name_SNP_toexclude_hp))

# exclude also 0 and MT chromosomes
# name_SNP_toexclude_chr <- c(SNP_table$Name[id_chr0], SNP_table$Name[id_chrMT], SNP_table$Name[id_chrXY])
name_SNP_toexclude_chr <- c(SNP_table$Name[id_chr0], SNP_table$Name[id_chrMT])
vect_lensnps <- c(vect_lensnps, length(name_SNP_toexclude_chr))

name_SNP_toexclude <- c(SNP_table_as$Name[id_SNPs_as], name_SNP_toexclude_hp, name_SNP_toexclude_chr)
vect_lensnps <- c(vect_lensnps, length(name_SNP_toexclude))

print(paste('n. excluded SNPs:', length(name_SNP_toexclude), '/', nrow(Full_table)))

vect_lensnps <- c(vect_lensnps,length(name_SNP_toexclude)/nrow(Full_table))
df_excl_SNP$n_SNPs <- vect_lensnps

##### write a table with the info for each removal of SNPs #####
write.table(x = df_excl_SNP, file = paste(outFile, 'info_QC.txt', sep = ''), dec = '.', sep = '\t', quote = FALSE, row.names = FALSE)


##### write filtered tables #####
# save just the info needed for the Full_table: BAF, LRR, GT
id_ft <- which(Full_table$Name %in% name_SNP_toexclude)

# save the correct number of columns for the Full table
id_column_BAF <- as.vector(sapply(0:(nsamples-1), function(x) 15 + (6*x)))
id_column_LRR <- as.vector(sapply(0:(nsamples-1), function(x) 16 + (6*x)))
id_initial_info <- 1:10

id_columns <- sort(c(id_initial_info, id_column_GT, id_column_BAF, id_column_LRR))
Full_table_filt <- Full_table[-id_ft,id_columns]

id_st <- which(SNP_table$Name %in% name_SNP_toexclude)
SNP_table_filt <- SNP_table[-id_st,]

# recompute the call freq for each sample:
# NOTE: for female exclude the chr Y PAR == 0 
# for the undefined, consider them as male
n_SNPs <- nrow(Full_table_filt)

id_column_GT_new <-  as.vector(sapply(0:(nsamples-1), function(x) 11 + (3*x)))

if(length(id_female)>0){
  GT_male_undef <- sapply(id_column_GT_new[-id_female], function(x) table(Full_table_filt[,x]))
  colnames(GT_male_undef) <- colnames(Full_table_filt[,id_column_GT_new[-id_female]])
  call_freq_filt_male_undef <- apply(GT_male_undef, 2, function(x) sum(x[1:3])/n_SNPs)
}else{
  GT_male_undef <- sapply(id_column_GT_new, function(x) table(Full_table_filt[,x]))
  colnames(GT_male_undef) <- colnames(Full_table_filt[,id_column_GT_new])
  call_freq_filt_male_undef <- apply(GT_male_undef, 2, function(x) sum(x[1:3])/n_SNPs)
}

id_chrY <- which(Full_table_filt$Chr == 'Y')
id_notPAR <- which(PAR_table$V2 == 0)
chrY_notPAR_filt <- Full_table_filt$Name[id_chrY][which(Full_table_filt$Name[id_chrY] %in% PAR_table$V1[id_notPAR])]

id <- which(Full_table_filt$Name %in% chrY_notPAR_filt)
Full_table_filt_female <- Full_table_filt[-id,] 

n_SNPs_female <- nrow(Full_table_filt_female)

call_freq_filt <- rep(0, nsamples)
if(length(id_female)>0){
  GT_female <- sapply(id_column_GT_new[id_female], function(x) table(Full_table_filt_female[,x]))
  colnames(GT_female) <- colnames(Full_table_filt_female[,id_column_GT_new[id_female]])
  call_freq_filt_female <- apply(GT_female, 2, function(x) sum(x[1:3])/n_SNPs_female)
  call_freq_filt[id_female_ST] <- call_freq_filt_female
  call_freq_filt[-id_female_ST] <- call_freq_filt_male_undef  
}else{
  call_freq_filt <- call_freq_filt_male_undef  
}

df_sample_callrate$call_rate_filt <- call_freq_filt

write.table(x = df_sample_callrate, file = paste(outFile, 'Samples_Table_filt.txt', sep = ''), dec = '.', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(x = SNP_table_filt, file = paste(outFile, 'SNP_Table_filt.txt', sep = ''), dec = '.', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(x = Full_table_filt, file = paste(outFile, 'Full_Data_Table_filt.txt', sep = ''), dec = '.', sep = '\t', quote = FALSE, row.names = FALSE)


################################## check the genotype (used then for the match) #################################
# id_column_GT <- as.vector(sapply(0:(nsamples-1), function(x) 11 + (6*x)))
GT_samples <- Full_table_filt[, id_column_GT_new]

# create a matrix for all the possible combination
registerDoParallel(cores = ncores)

GT_comp <- foreach(i = 1:(nsamples-1))%dopar%{
  # for(i in 1:N){
  print(i)
  id_in <- rep(i, nsamples-i)
  id_fin <- (i+1):nsamples
  
  # s_m_vect <- c(s_m_vect, mapply(function(x,y) frobenious_dotprod(X[[x]],X[[y]],K), x = id_in, y = id_fin))
  mapply(function(x,y) compareGT(GT_s1 = GT_samples[,x], GT_s2 = GT_samples[,y]), x = id_in, y = id_fin)
}

GT_comp_vect <- unlist(GT_comp)

GT_comp <- matrix(0, nrow=nsamples, ncol=nsamples) 
GT_comp[lower.tri(GT_comp, diag = FALSE)] <- GT_comp_vect 
GT_comp <- GT_comp + t(GT_comp) - diag(diag(GT_comp)) 
diag(GT_comp) <- 1

# add samples name
colnames(GT_comp) <- sample_name_FT
rownames(GT_comp) <- sample_name_FT

write.table(x = GT_comp, file = paste(outFile, 'GT_comp.txt', sep = ''), dec = '.', sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)


################################## save LRR and BAF and GT from the Full_table_filt #################################
# extract LRR and BAF values for each sample, reorder position
id_column_BAF_new <- 12 + 3 * (0:(nsamples-1))
id_column_LRR_new <- 13 + 3 * (0:(nsamples-1))

### LRR
LRR_table <- data.frame(CHROM = Full_table_filt$Chr, POS = Full_table_filt$Position, ID = Full_table_filt$Name, stringsAsFactors = FALSE)
# LRR_table$CHROM <-  Full_table_filt$Chr
# LRR_table$POS <- Full_table_filt$Position

# rename X and Y with 23 and 24
LRR_table$CHROM[which(LRR_table$CHROM == "X")] <- '23'
LRR_table$CHROM[which(LRR_table$CHROM == "Y")] <- '24'

LRR_table_sample <- Full_table_filt[, id_column_LRR_new]
colnames(LRR_table_sample) <- sample_name_FT
LRR_table <- cbind(LRR_table, LRR_table_sample)
LRR_table_temp <- LRR_table

### BAF
BAF_table <- LRR_table_temp[,1:3]
BAF_table <- cbind(BAF_table, Full_table_filt[, id_column_BAF_new])
colnames(BAF_table) <- colnames(LRR_table_temp)
BAF_table_temp <- BAF_table

### GT
GT_table <- LRR_table_temp[,1:3]
GT_table <- cbind(GT_table, Full_table_filt[, id_column_GT_new])
colnames(GT_table) <- colnames(LRR_table_temp)
GT_table_temp <- GT_table

LRR_table <- vector(mode = 'list', length = length(unique(LRR_table_temp$CHROM)))
BAF_table <- vector(mode = 'list', length = length(unique(LRR_table_temp$CHROM)))
GT_table <- vector(mode = 'list', length = length(unique(LRR_table_temp$CHROM)))

for(id  in 1:length(unique(LRR_table_temp$CHROM))){
  
  print(id)
  chr <- unique(LRR_table_temp$CHROM)[id]
  
  LRR_chr <- LRR_table_temp[LRR_table_temp$CHROM == chr,]
  BAF_chr <- BAF_table_temp[BAF_table_temp$CHROM == chr,]
  GT_chr <- GT_table_temp[GT_table_temp$CHROM == chr,]
  
  ord_CHR <- order(LRR_chr$POS)
  
  LRR_table[[id]] <- LRR_chr[ord_CHR,]
  BAF_table[[id]] <- BAF_chr[ord_CHR,]
  GT_table[[id]] <-  GT_chr[ord_CHR,]
  
}

LRR_table <- do.call(rbind, LRR_table)
BAF_table <- do.call(rbind, BAF_table)
GT_table <- do.call(rbind, GT_table)

# replace NA with dot


### save
write.table(x = LRR_table, file = paste(outFile, 'LRR_table.txt', sep = ""), sep = '\t', quote = FALSE, dec = '.', row.names = FALSE, col.names = TRUE, na = '.')
write.table(x = BAF_table, file = paste(outFile,'BAF_table.txt', sep = ""), sep = '\t', quote = FALSE, dec = '.', row.names = FALSE, col.names = TRUE, na = '.')
write.table(x = GT_table, file = paste(outFile,'GT_table.txt', sep = ""), sep = '\t', quote = FALSE, dec = '.', row.names = FALSE, col.names = TRUE, na = '.')






