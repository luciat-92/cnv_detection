# Compute the heatmap for the GT amtch and update the annotation file with the GT match
# Written by Lucia Trastulla -- email: lucia_trastulla@psych.mpg.de

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library('ComplexHeatmap'))
suppressPackageStartupMessages(library('RColorBrewer'))
suppressPackageStartupMessages(library('circlize'))

parser <- ArgumentParser(description = "Compute GT match between samples and update annotation file" )

parser$add_argument("--ann_file", type = "character", help = ".csv file, annotation file manually edited")
parser$add_argument("--comp_GT_file", type = "character", help = ".txt file, compatibility matrix already computed")
parser$add_argument("--Samples_table_file", type = "character", help = ".txt file, sample table after quality control")
parser$add_argument("--outf", type="character", help = "Output file [basename only]")

#######
args <- parser$parse_args()
ann_file <- args$ann_file
Samples_table_file <- args$Samples_table_file
comp_GT_file <- args$comp_GT_file
outFile <- args$outf
#######

comp_GT <- read.table(comp_GT_file, header = TRUE, stringsAsFactors = FALSE, sep = '\t', dec = '.')
Sample_table_filt <- read.table(Samples_table_file, header = TRUE, stringsAsFactors = FALSE, sep = '\t', dec = '.')
ann_samples <- read.csv(ann_file, header = TRUE, stringsAsFactors = FALSE)

# save the different samples considered
samples_name <- unique(ann_samples$Sample_name)
P <- length(samples_name)
colnames(comp_GT) <- rownames(comp_GT)

# filter the GT considering only the samples in the annotation file:
sample_name_ann <- sapply(1:nrow(ann_samples), function(x) paste(ann_samples$SentrixBarcode_A[x], ann_samples$SentrixPosition_A[x], sep = "_"))
id <- which(rownames(comp_GT) %in% sample_name_ann)
comp_GT <- comp_GT[id,id]

diss <- as.dist(1-comp_GT)
hclust_sample <- hclust(d = diss, method = 'ward.D')



##### set colours #####
coul <- brewer.pal(9, "BuPu") 
coul <- colorRampPalette(coul)(10)
coul[1] <- 'white'
nbins <- 10
######################

# create the colour palette
col_name <- brewer.pal(n = 8, name = 'Paired')
col_name <- colorRampPalette(col_name)(P)
if(P <= length(col_name)){

  col_name <- col_name[1:P]

}else{

  n_rep <- ceiling(P/8)
  col_name <- rep(col_name, n_rep)
  col_name <- col_name[1:P]

}

col_match <- col_name
names(col_match) <- samples_name


ord_samples <- sapply(names(comp_GT) ,function(i) which(i  == sample_name_ann))
df_samples <- data.frame(sample_id = names(comp_GT), Sample_name = ann_samples$Sample_name[ord_samples])
df_samples <- data.frame(Sample_name = ann_samples$Sample_name[ord_samples])


ha = HeatmapAnnotation(df = df_samples, col = list(Sample_name = col_match),  annotation_legend_param = list(Sample_name = list(ncol = ceiling(P/45), title = "Samples", title_gp = gpar(fontsize = 15, fontface = "bold"), title_position = "topcenter",                                                       
                                                                                                                            grid_height = unit(8, "mm"), gap = unit(4, "mm"), labels_gp = gpar(fontsize = 11))))


hm <- Heatmap(comp_GT, name = "GT match",  cluster_rows = hclust_sample, cluster_columns = hclust_sample, row_dend_reorder = FALSE, column_dend_reorder = FALSE, top_annotation = ha,
              show_row_names = FALSE, show_column_names = TRUE, show_row_dend = FALSE, show_column_dend = TRUE,
                 column_title_gp = gpar(fontsize = 25, fontface = "bold"), col = colorRamp2(seq(min(comp_GT), 1, length = nbins), coul), heatmap_legend_param = list(color_bar = "continuous", legend_height = unit(10, "cm"), labels_gp = gpar(fontsize = 17), title_gp = gpar(fontsize = 15, fontface = "bold")))


padding = unit.c(unit(9, "mm"), unit(0.5, "cm"), unit(c(3, 3), "mm"))

png(paste(outFile, 'hm_GT.png', sep = ""), width = 4000+(ceiling(P/45)*200), height = 4000, res = 300, pointsize = 12)
draw(hm, padding = padding)
dev.off()

pdf(paste(outFile, 'hm_GT.pdf', sep = ""), width = 15, height = 15)
draw(hm, padding = padding)
dev.off()

#### check for mismatch
# same element should have a match higher than 96%, if a sample is of poor quality labelled in a different way
#############
list.flatten = function (x, use.names = TRUE, classes = "ANY") 
{
  len <- sum(rapply(x, function(x) 1L, classes = classes))
  y <- vector("list", len)
  i <- 0L
  items <- rapply(x, function(x) {
    i <<- i + 1L
    y[[i]] <<- x
    TRUE
  }, classes = classes)
  if (use.names && !is.null(nm <- names(items))) 
    names(y) <- nm
  y
}

#############

gr_names <- NULL
logic_vect <- rep(TRUE, nrow(comp_GT))
names_samples <- rownames(comp_GT)

for(i in 1:nrow(comp_GT)){
  
  if(logic_vect[i]){
    
    id <- which(comp_GT[i,] > 0.96)
    id <- id[which(logic_vect[id])]
    # print(names_samples[id])
    gr_names <- list.flatten(list(gr_names, names_samples[id]))
    logic_vect[id] <- F 
    
  }
  
}

# create a file that match the orignal sample annotations and the one obtained from GT match
# if they are concordant use the original name, otherwise reassign
ann_samples$sample_complete_id <- sample_name_ann 
ann_samples$Sample_name_new <- rep(0, nrow(ann_samples))

count <- 0

for(i in 1:length(gr_names)){
  
  id <- which(sample_name_ann %in% gr_names[[i]])
  cond1 <- length(unique(ann_samples$Sample_name[id]))==1 # check if the matching one have the same name
  cond2 <- length(which(ann_samples$Sample_name %in% ann_samples$Sample_name[id])) == length(id) # check if the sample considered are labelled with other ones
  cond <- cond1 & cond2 # both must be TRUE
  
  if(cond){ann_samples$Sample_name_new[id] <- ann_samples$Sample_name[id]}
  
  else{
    count <- count+1
    ann_samples$Sample_name_new[id] <- sprintf('MISMATCHED_SAMPLE_%i', count)}

}



##############
# add gender fro Sample table
id <- sapply(ann_samples$sample_complete_id, function(x) which(Sample_table_filt$sample_ID == x))
gender <- Sample_table_filt$gender[id]
# NOTE: keep undefined info, chane manually if necessary
gender[which(gender == 'Male')] <- 'M'
gender[which(gender == 'Female')] <- 'F'
gender[which(gender == 'Undef')] <- 'U'
ann_samples$Gender <- gender
##############


############################
samples_name_new <- unique(ann_samples$Sample_name_new)
P <- length(unique(ann_samples$Sample_name_new))
# make a new plot with the correct annotation
# create the colour palette
col_name <- brewer.pal(n = 8, name = 'Paired')
col_name <- colorRampPalette(col_name)(P)
if(P <= length(col_name)){
  
  col_name <- col_name[1:P]
  
}else{
  
  n_rep <- ceiling(P/8)
  col_name <- rep(col_name, n_rep)
  col_name <- col_name[1:P]
  
}

col_match <- col_name
names(col_match) <- samples_name_new


ord_samples <- sapply(names(comp_GT) ,function(i) which(i  == sample_name_ann))
df_samples <- data.frame(sample_id = names(comp_GT), Sample_name = ann_samples$Sample_name_new[ord_samples])
df_samples <- data.frame(Sample_name = ann_samples$Sample_name_new[ord_samples])


ha = HeatmapAnnotation(df = df_samples, col = list(Sample_name = col_match),  annotation_legend_param = list(Sample_name = list(ncol = ceiling(P/45), title = "Samples", title_gp = gpar(fontsize = 15, fontface = "bold"), title_position = "topcenter",                                                       
                                                                                                                                grid_height = unit(8, "mm"), gap = unit(4, "mm"), labels_gp = gpar(fontsize = 11))))


hm <- Heatmap(comp_GT, name = "GT match",  cluster_rows = hclust_sample, cluster_columns = hclust_sample, row_dend_reorder = FALSE, column_dend_reorder = FALSE, top_annotation = ha,
              show_row_names = FALSE, show_column_names = TRUE, show_row_dend = FALSE, show_column_dend = TRUE,
              column_title_gp = gpar(fontsize = 25, fontface = "bold"), col = colorRamp2(seq(min(comp_GT), 1, length = nbins), coul), heatmap_legend_param = list(color_bar = "continuous", legend_height = unit(10, "cm"), labels_gp = gpar(fontsize = 17), title_gp = gpar(fontsize = 15, fontface = "bold")))


padding = unit.c(unit(9, "mm"), unit(0.5, "cm"), unit(c(3, 3), "mm"))

png(paste(outFile, 'hm_GT_correct.png', sep = ""), width = 4000+(ceiling(P/45)*200), height = 4000, res = 300, pointsize = 12)
draw(hm, padding = padding)
dev.off()

pdf(paste(outFile, 'hm_GT_correct.pdf', sep = ""), width = 15, height = 15)
draw(hm, padding = padding)
dev.off()


# save the table
write.table(x = ann_samples, file = paste(outFile, 'annotation_file_GTmatch.csv', sep = ''), dec = '.', sep = ',', quote = FALSE, row.names = FALSE)








