#### Script to filter data for pipeline report ####
### Proteins to keep

keep <- read.table("data/selected_proteins_final.txt")[,1]
keep_peptides <- read.table("data/selected_peptides_final.txt")[,1]

# sa1 <- readRDS("~/polybox/SemPro_Picotti/euler_output/significant_all_ATP_k15.rds")
# peps <- sa1$scores %>% filter(Protein %in% keep) %>% pull(Peptide) %>% unique()
# write.table(peps,"data/selected_peptides_final.txt",quote = FALSE,row.names = FALSE,col.names=FALSE)

### Spectronaut output
Spectronaut <- read.csv("~/polybox/SemPro_Picotti/001_ATP/018/ATP/20211116_085130_018_atp_Report.tsv",sep = "\t")

f_Spectronaut <- Spectronaut %>% filter(PG.ProteinAccessions %in% keep) 

### Python output

# FT
MLL_ATP_FT <- read.csv("~/polybox/SemPro_Picotti/001_ATP/018/ATP/python_output/MLL_ATP_FT.csv")
solution_ATP_FT <- read.csv("~/polybox/SemPro_Picotti/001_ATP/018/ATP/python_output/solution_ATP_FT.csv")

f_solution_ATP_FT <- solution_ATP_FT %>% filter(peptide %in% keep_peptides) 
f_MLL_ATP_FT <- MLL_ATP_FT %>% filter(peptide %in% keep_peptides)

# HT
MLL_ATP_HT <- read.csv("~/polybox/SemPro_Picotti/001_ATP/018/ATP/python_output/MLL_ATP_HT.csv")
solution_ATP_HT <- read.csv("~/polybox/SemPro_Picotti/001_ATP/018/ATP/python_output/solution_ATP_HT.csv")

f_solution_ATP_HT <- solution_ATP_HT %>% filter(peptide %in% keep_peptides)
f_MLL_ATP_HT <- MLL_ATP_HT %>% filter(peptide %in% keep_peptides)

write.csv(f_solution_ATP_FT,"~/polybox/SemPro_Picotti/006_PipelineReport/Rmarkdown/data/ATP/python_output/solution_ATP_FT.csv",row.names = FALSE)

### Clustering

clusters <- readRDS("~/polybox/SemPro_Picotti/euler_output/clusters_ATP_k15.rds")

f_clusters <- lapply(clusters[2:3],function(x){x %>% filter(Peptide %in% keep_peptides)})
f_clusters$plot <- clusters$plot
f_clusters$manual_order <- clusters$manual_order

### all_results

all_results <- readRDS("~/polybox/SemPro_Picotti/euler_output/combined_results_ATP.rds")
f_all_results <- lapply(all_results$ATP,function(x){x %>% filter(peptide %in% keep_peptides)})

### significant_all

sa <- readRDS("~/polybox/SemPro_Picotti/euler_output/significant_all_ATP_k15.rds")
f_sa <- list()
f_sa$scores <- sa$scores %>% filter(Peptide %in% keep_peptides)
f_sa$long <- sa$long %>% filter(peptide %in% keep_peptides)
f_sa$AA_level <- sa$AA_level %>% filter(Peptide %in% keep_peptides)
f_sa$Protein_level <- sa$Protein_level %>% filter(Protein %in% keep)
#### save output

save(f_clusters,f_MLL_ATP_FT,f_MLL_ATP_HT,f_solution_ATP_FT,f_solution_ATP_HT,f_Spectronaut,f_sa,f_all_results,
                             file = "~/polybox/SemPro_Picotti/006_PipelineReport/Rmarkdown/data/pipeline_report_data.RData")


# t <- readRDS("~/polybox/SemPro_Picotti/euler_output/significant_all_ATP_k15.rds")

