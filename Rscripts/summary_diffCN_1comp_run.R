# Compute the summary report for different CN detected in each comparison
# Use the summary.tab file produced from each comparison
# Find the centromere position, only one dat.tab file is needed
# filter out CN of poor quality
# Written by Lucia Trastulla -- email: lucia_trastulla@psych.mpg.de


suppressPackageStartupMessages(library(argparse))
##############################
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('gridExtra'))
suppressPackageStartupMessages(library('grid'))
##############################

parser_plot <- ArgumentParser(description="cnv summary (differnce in predicted CN) for all the comparisons in a sample")

parser_plot$add_argument("--inputdir", type = "character", help="input folder")
parser_plot$add_argument("--fold_fun", type = "character", help="folder for the input function")
parser_plot$add_argument("--post", type = "character", nargs = '*', help=" vector of genome studio name identifyng the post repreogramming lines")
parser_plot$add_argument("--pre", type = "character",  help="genome studio name for the pre reprogramming lines")
parser_plot$add_argument("--ann_file", type = "character",  help="annotation for the samples")
parser_plot$add_argument("--sample_name", type = "character", help="vector of name identifying the sample (patient)")
parser_plot$add_argument("--qs_thr", type = "double", default = 2, help="quality score threshold")
parser_plot$add_argument("--nSites_thr", type = "integer", default = 10,  help="number of sites thersold for CN1")
parser_plot$add_argument("--nHet_thr", type = "integer", default = 10, help="number of sites thersold for CN3")
parser_plot$add_argument("--CN_len_thr", type = "integer", default = 200, help="CN length thr in kb")
parser_plot$add_argument("--sex", type="character", default = "M", help = "sex of the sample (patient)")
parser_plot$add_argument("--width_plot", type="double", default = 10)
parser_plot$add_argument("--outf", type="character", help = "Output file ")

args <- parser_plot$parse_args()

inputdir <- args$inputdir
fold_fun <- args$fold_fun
post <- args$post 
pre <- args$pre
ann_file <- args$ann_file
qs_thr <- args$qs_thr
nSites_thr <- args$nSites_thr
nHet_thr <- args$nHet_thr
CN_len_thr <- args$CN_len_thr
sample_name <- args$sample_name
sex_sample <- args$sex
width_plot <- args$width_plot
outFile <- args$outf

# set heigth_plot based on the number of post
#heigth_plot <- 9+(length(post)-1)
# only 1 sample 
heigth_plot <- 4

#### load function
source(paste(fold_fun, 'plot_summary_cnv_function.R', sep = ""))
####

setwd(sprintf('%s/outdir_%s_%s/',inputdir, sample_name, pre)) 

comp_dir <- sapply(post, function(x) sprintf('output_%s/', x))
cnv_file <- sapply(comp_dir, function(x) paste(x, 'summary.tab', sep = ""))

M <- 10^6


cnv_table_list <- vector(mode = 'list', length = length(post))
cnv_filt <- vector(mode = 'list', length = length(post))

for(i in 1:length(post)){
  
  cnv_table_list[[i]] <- read.csv(cnv_file[[i]], header = TRUE, sep = '\t', skip=3, stringsAsFactors = FALSE)
  # exclude the first column (contain only RG)
  cnv_table_list[[i]] <- cnv_table_list[[i]][,-1]
 
  # insert colnames
  colnames(cnv_table_list[[i]]) <- c('Chr', 'Start', 'End', sprintf('CN_post_%s', post[i]),  sprintf('CN_pre_%s', pre), 'QS', 'nSites_post', 'nHets_post', 'nSites_pre', 'nHets_pre')
  # add length column
  cnv_table_list[[i]]$Len <- cnv_table_list[[i]][,3] - cnv_table_list[[i]][,2]
  
  # report the id of the rows to delete:
  # quality score
  id_qs <- which(cnv_table_list[[i]]$QS < qs_thr)
  # length 
  id_len <- which(cnv_table_list[[i]]$Len <= CN_len_thr*1000)
  
  #### Haploipd chr
  if(sex_sample == 'F'){
    
    id_chrY <- which(cnv_table_list[[i]]$Chr == 24)
    chr_lim <- 24
    id_Hap <- id_chrY
    
  }else{
    
    chr_lim <- 23
    # number of sites for delition in pre and post
    id_nH_pre_XY <- which(cnv_table_list[[i]]$nHets_pre < nHet_thr & cnv_table_list[[i]][,4] > 1 & cnv_table_list[[i]]$Chr >= chr_lim)
    id_nH_post_XY <- which(cnv_table_list[[i]]$nHets_post < nHet_thr & cnv_table_list[[i]][,5] > 1 & cnv_table_list[[i]]$Chr >= chr_lim)
    
    id_Hap <- unique(c(id_nH_pre_XY, id_nH_post_XY))
  }
  
  #### Autosomal chr 
  # number of sites for delition in pre and post
  id_nS_pre <- which(cnv_table_list[[i]]$nSites_pre < nSites_thr & cnv_table_list[[i]][,4] == 1 & cnv_table_list[[i]]$Chr < chr_lim)
  id_nH_pre <- which(cnv_table_list[[i]]$nHets_pre < nHet_thr & cnv_table_list[[i]][,4] == 3 & cnv_table_list[[i]]$Chr < chr_lim)
  
  id_nS_post <- which(cnv_table_list[[i]]$nSites_post < nSites_thr & cnv_table_list[[i]][,5] == 1 & cnv_table_list[[i]]$Chr < chr_lim)
  id_nH_post <- which(cnv_table_list[[i]]$nHets_post < nHet_thr & cnv_table_list[[i]][,5] == 3 & cnv_table_list[[i]]$Chr < chr_lim)
  
  ### total
  id_tot <- unique(c(id_qs, id_len,id_Hap, id_nS_pre,id_nH_pre, id_nS_post,id_nH_post))
  cnv_filt[[i]] <- cnv_table_list[[i]][-id_tot,]
  
  # save the original table but add the column Filt
  cnv_table_list[[i]]$Filt <- FALSE
  cnv_table_list[[i]]$Filt[id_tot] <- TRUE
  
  # add the column with the different CN (TRUE if the CN is the same, FALSE otherwise)
  cnv_filt[[i]]$diff_CN <- cnv_filt[[i]][,4] == cnv_filt[[i]][,5]
  
  write.table(x = cnv_table_list[[i]], file = paste(outFile, sprintf('summary_pair_%s_%s.tab', sample_name, post[i]), sep = ""), col.names = TRUE, quote = FALSE, row.names = FALSE, sep ='\t')
  
}
  
cnv_diff <- lapply(cnv_filt, function(x) x[!x$diff_CN,])

# create the file (df)
# read the annotation file
ann_table <- read.csv(ann_file, stringsAsFactors = FALSE, header = TRUE)

# sample_name_ann <- strsplit(sample_name, '_')[[1]][1]
complete_id <- sapply(1:nrow(ann_table), function(x) paste(ann_table$SentrixBarcode_A[x],  ann_table$SentrixPosition_A[x], sep = "_"))
pre_name <- ann_table$ID[which(complete_id ==  pre)]
sample_name_ann <- ann_table$Sample_name_new[which(complete_id ==  pre)]

post_name <- sapply(post, function(x) ann_table$ID[which(complete_id ==  x)])

df <- data.frame(Sample_name = rep(sample_name_ann, 25*length(post)), pre_ID = rep(pre_name, 25*length(post)), pre_Chip_ID = rep(pre, 25*length(post)))
df$post_ID <- as.vector(sapply(post_name, function(x) rep(x, 25)))
df$post_Chip_ID <- as.vector(sapply(post, function(x) rep(x, 25)))
df$Chr <- rep(c(1:22, 'X', 'Y', 'tot'),length(post))
df$Ndiff_200kb <- rep(0, 25*length(post))
df$Ndiff_1Mb <- rep(0, 25*length(post))


for(i in 1:length(post)){
  
  id <- which(!cnv_filt[[i]]$diff_CN)
  
  if(length(id)>0){
   
    chr_diff <- sort(unique(cnv_filt[[i]]$Chr[id]))
    n_diff200 <- sapply(chr_diff, function(x) length(which(cnv_filt[[i]]$Chr[id] == x))) 
    n_diff1 <- sapply(chr_diff, function(x) length(which(cnv_filt[[i]]$Chr[id] == x & cnv_filt[[i]]$Len[id] > 10^6)))
    
    df$Ndiff_200kb[25*(i-1) + chr_diff] <- n_diff200
    df$Ndiff_200kb[25*(i-1) + 25] <- sum(n_diff200)
    
    df$Ndiff_1Mb[25*(i-1) + chr_diff] <- n_diff1
    df$Ndiff_1Mb[25*(i-1) + 25] <- sum(n_diff1)
    
    
      
  }
}


# write table
write.table(x = df, file = paste(outFile, sprintf('%s_CNV_diff.txt', sample_name), sep = ''), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)


# make the plot for each post
df_plot <- rbind(df[,1:6], df[,1:6])
df_plot$Ndiff = c(df$Ndiff_200kb, df$Ndiff_1Mb)
df_plot$thr = c(rep('0.2 Mb', nrow(df)), rep('1 Mb', nrow(df)))
df_plot$Chr <- factor(df_plot$Chr, levels = c(1:22, 'X', 'Y', 'tot'))

plot_CN <- ggplot(data = df_plot, aes(x = Chr, y = Ndiff))+
        geom_bar(stat = "identity", aes(fill = thr), colour = "black",  position=position_dodge())+
  theme_bw()+
  ggtitle(sprintf('Different CN wrt %s', pre_name))+
  ylab("N. of diff CN")+xlab("Chr")+
  facet_wrap(~post_ID, ncol = 1)+
  scale_fill_discrete(name = "length >")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=14), plot.title = element_text(hjust = 0, size = 17))

ggsave(filename = paste(outFile, sprintf('%s_CNV_diff.pdf', sample_name), sep = ''), plot = plot_CN, width = width_plot, height = heigth_plot)
ggsave(filename = paste(outFile, sprintf('%s_CNV_diff.png', sample_name), sep = ''), plot = plot_CN, width = width_plot, height = heigth_plot)



# modify the single summary files
post_file <- mapply(function(x,y) paste(x, 'summary.', y, '.tab', sep = ""), x = comp_dir, y = post)
pre_file <- sapply(comp_dir, function(x) paste(x, 'summary.', pre, '.tab', sep = ""))


# create a function that filter the tab files
filter_cnv <- function(file, sex_sample, qs_thr, nSites_thr, nHet_thr, CN_len_thr){
  
  cnv_table <- read.csv(file, header = FALSE, sep = '\t', skip=0, stringsAsFactors = FALSE)
  cnv_table <- cnv_table[-1,-1]
  colnames(cnv_table) <-  c('Chr', 'Start', 'End', 'CN',  'QS', 'nSites', 'nHets')
  cnv_table <- data.frame(apply(cnv_table, 2, as.numeric))
  
  cnv_table$Len <- cnv_table[,3] - cnv_table[,2]
  
  # report the id of the rows to delete:
  # quality score
  id_qs <- which(cnv_table$QS < qs_thr)
  # length 
  id_len <- which(cnv_table$Len <= CN_len_thr*1000)
  
  #### Haploipd chr
  if(sex_sample == 'F'){
    
    id_chrY <- which(cnv_table$Chr == 24)
    chr_lim <- 24
    id_Hap <- id_chrY
    
  }else{
    
    chr_lim <- 23
    # number of sites for delition in pre and post
    id_nH_XY <- which(cnv_table$nHets < nHet_thr & cnv_table[,4] > 1 & cnv_table$Chr >= chr_lim)
    
    id_Hap <- id_nH_XY
  }
  
  #### Autosomal chr 
  # number of sites for delition in pre and post
  id_nS <- which(cnv_table$nSites < nSites_thr & cnv_table[,4] == 1 & cnv_table$Chr < chr_lim)
  id_nH <- which(cnv_table$nHets < nHet_thr & cnv_table[,4] == 3 & cnv_table$Chr < chr_lim)
  
  ### total
  id_tot <- unique(c(id_qs, id_len,id_Hap, id_nS, id_nH))
  
  cnv_table$Filt <- FALSE
  cnv_table$Filt[sort(id_tot)] <- TRUE
  
  return(cnv_table)
  
}


post_filt <- lapply(post_file, function(x) filter_cnv(file = x , sex_sample = sex_sample, qs_thr = qs_thr, nSites_thr = nSites_thr, nHet_thr = nHet_thr, CN_len_thr = CN_len_thr))
pre_filt <- lapply(pre_file, function(x) filter_cnv(file = x , sex_sample = sex_sample, qs_thr = qs_thr, nSites_thr = nSites_thr, nHet_thr = nHet_thr, CN_len_thr = CN_len_thr))

# save the tables
for(i in 1:length(post)){
  
  write.table(x = post_filt[[i]], file = paste(comp_dir[[i]], sprintf('summary_QC.%s.tab', post[i]), sep = ""), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
  
  write.table(x = pre_filt[[i]], file = paste(comp_dir[[i]], sprintf('summary_QC.%s.tab', pre), sep = ""), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
  
}


