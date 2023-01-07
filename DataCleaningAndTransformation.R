# Applying WGCNA to Human Alopecia Areata Skin Biopsy Samples
# Public GEO data: GSE68801 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68801)
# WGCNA: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559

# AA Sample Summary:
# AAP (patchy type disease, n=22 patients; lesional samples n=22, nonlesional samples n=20)
# AAP.T (transient patchy type disease, disease duration less than 1 year, n=6; lesional samples n=6, nonlesional samples n=6)
# AU (alopecia universalis, n=23 patients)
# AT (alopecia totalis, n=9 patients)
# Normal (healthy controls, n=36 patients)

geodir = "/opt/projects/../AA_transcriptomics/GSE68801/expression_data/"
datadir = "/opt/projects/../AA_transcriptomics/GSE68801/WGCNA/Data_WGCNA/"

#=====================================================================================
# Data Cleaning
#=====================================================================================
setwd(geodir)

norm_data = read.table("GSE68801.normalised_data.tsv",sep="\t",header=T)
exp = norm_data %>% relocate(ensembl_gene_id , .before = ProbeName) # move the gene ID col forward 
norm_data_w_pheno_data_samples_as_rows = read.table("GSE68801.norm_data_w_pheno_data_samples_as_rows.tsv",sep="\t",header=T)
clinics = norm_data_w_pheno_data_samples_as_rows[,c(1:11)] # only the first 11 cols are useful info here
clinics$name = paste(paste("s",clinics$sample,sep = ""), clinics$disease_status,sep = "_") # rename sample names, adding on disease state info
clinics$name = paste(clinics$name, clinics$skintype,sep = "_") # refine sample names, adding on skin type info (L, lesional vs. NL, nonlesional)

# combine clinical info for samples of interest
clinicsA = clinics[!grepl('AAP.T', clinics$disease_status),]
clinicsA = clinicsA[!grepl('N', clinicsA$disease_status),]
dim(clinicsA) # 63 samples left after removing AAP.T and N samples
clinicsA = clinicsA %>% relocate(name , .before = sample)
clinicsAAP.L = clinicsA[grepl('AAP_L', clinicsA$name),]
clinicsAAP.NL = clinicsA[grepl('AAP_NL', clinicsA$name),]

# create clinical info table for Normal samples
clinicsN = clinics[grepl('N', clinics$disease_status),]
clinicsN = clinicsN %>% relocate(name , .before = sample)

colnames(exp) = c("ensembl_gene_id","ProbeName",clinics$name) # rename col names of exp table to match clinical info table
geneanno = exp[,c(1,2)] # extract gene annotation info for future use

#=====================================================================================
# Select Data to Feed into WGCNA Algorithm
#=====================================================================================
# expression info table for Normal samples
expN = cbind( exp[,c(1,2)], exp[grepl('Normal', colnames(exp))]) #33
# expression info table for all AA samples
expA = cbind( exp[,c(1,2)], 
              exp[grepl('AAP_L', colnames(exp))],
              exp[grepl('AAP_NL', colnames(exp))],
              exp[grepl('AU', colnames(exp))],
              exp[grepl('AT', colnames(exp))]) #63
# expression info table for all of the samples of interest (all AA + N)
expA.N = cbind( exp[,c(1,2)],
                exp[grepl('AAP_L', colnames(exp))],
                exp[grepl('AAP_NL', colnames(exp))],
                exp[grepl('AU', colnames(exp))],
                exp[grepl('AT', colnames(exp))],
                exp[grepl('Normal', colnames(exp))]) # 96
expAAP.L = cbind( exp[,c(1,2)], exp[grepl('AAP_L', colnames(exp))]) #18
expAAP.NL = cbind( exp[,c(1,2)], exp[grepl('AAP_NL', colnames(exp))]) #17

#=====================================================================================
# Save tables for WGCNA Algorithm
#=====================================================================================
setwd(datadir)

# write exp tables
# write.table(expA, file='GSE68801_WGCNA_expA.tsv', quote=FALSE, col.names = TRUE, row.names = TRUE, sep = "\t",na = "")
aa_list = list(expA, expAAP.L, expAAP.NL, expA.N , expN) # create a list of 5 exp tables to save
names(aa_list) = c('expA', 'expAAP.L', 'expAAP.NL', 'expA.N', 'expN') # name the data frames
# save each new data frame as an individual .tsv file based on its name
lapply(1:length(my_list), function(i) write.csv(aa_list[[i]], 
                                                file = paste0("GSE68801_WGCNA_", names(my_list[i]), ".tsv"),
                                                quote=FALSE, col.names = TRUE, row.names = TRUE, sep = "\t",na = ""))

# write clinical info tables
write.table(clinicsA, file='GSE68801_WGCNA_clinicsA.tsv', quote=FALSE, col.names = TRUE, row.names = TRUE, sep = "\t",na = "")
write.table(clinicsN, file='GSE68801_WGCNA_clinicsN.tsv', quote=FALSE, col.names = TRUE, row.names = TRUE, sep = "\t",na = "")
# write gene annotation table
write.table(geneanno, file='GSE68801_WGCNA_geneanno.tsv', quote=FALSE, col.names = TRUE, row.names = TRUE, sep = "\t",na = "")

