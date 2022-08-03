# Python import function
# Author: Monika Pepelnjak
# Description: The function takes the import from Spectronout and transforms the data for python fitting approach
# The function also saves the files in the described path

#Python analysis of thermal stabilization
# The script firstly filters the data - for now only complete cases, proteotypic peptides, removes contaminants and iRT peptides
# In the next step peptides are being scaled. 


Spectronaut_to_python <- function(sm, data,annotation_file, common_path , include_replicates = c(1,2,3)) {
  Annotation <- read.csv(annotation_file)  %>%
    mutate(run = paste(.$condition, .$temperature, .$replicate, sep="_")) %>%
    filter(condition %in% c("Control",sm))
  
  colnames(data)[grepl("StrippedSequence", colnames(data))] <- "EG.StrippedSequence"
  colnames(data)[grepl("IsProteotypic", colnames(data))] <- "EG.IsProteotypic"
  
  FT_peptides <- data$EG.StrippedSequence[data$PEP.DigestType....Trypsin.P. == "Specific"] %>% unique()
  HT_peptides <- data$EG.StrippedSequence[data$PEP.DigestType....Trypsin.P. %in% c( "Specific-N", "Specific-C")] %>% unique()
  
  Annotated <- data %>% 
    #join with the annotation file
    plyr::join(., Annotation, by="R.FileName") %>%
    filter(PEP.Quantity > 100) %>%
    filter(EG.IsProteotypic == "True",
           #remove the iRT peptides
           PG.FastaFiles == "200331_ecoli", 
           replicate %in% include_replicates)
  
  #only keep the peptides with full profiles
  Annotated_full <- Annotated %>%
    reshape2::dcast(EG.StrippedSequence ~ run, value.var = "PEP.Quantity", fun.aggregate = function(x) mean(x, na.rm=T)) %>% 
    reshape2::melt() %>% 
    `colnames<-`(c("EG.StrippedSequence", "run", "PEP.Quantity")) %>%
    plyr::join(., Annotation, by="run")
  
  Annotated_full %<>%
    group_by(EG.StrippedSequence) %>%
    mutate(number_missing = sum(is.na(PEP.Quantity)))
  
  Annotated_full %<>%
    mutate(which_missing = ifelse(is.na(PEP.Quantity), R.FileName, NA))
  
  #python export file
  hundred_peptides <-Annotated_full$EG.StrippedSequence[Annotated_full$number_missing < 5 ] %>% 
    as.character() %>% 
    unique() 
  
  python_export <- Annotated_full %>% 
    ungroup() %>%
    filter(.$EG.StrippedSequence %in% hundred_peptides) #%>% 
  
  python_export <- python_export %>% 
    ungroup() %>% 
    dplyr::select(-which_missing) %>%
    group_by(EG.StrippedSequence, condition, replicate) %>%
    na.omit() %>%
    #scaling the data pcr:
    mutate( pcr_minmax = (PEP.Quantity - min(PEP.Quantity))/(max(PEP.Quantity)-min(PEP.Quantity))) %>%
    #  mutate( pc_minmax = (2^PEP.Quantity - min(2^PEP.Quantity))/(max(2^PEP.Quantity)-min(2^PEP.Quantity))) %>%
    ungroup() 
  
  export_table <- python_export %>%
    dplyr::select(EG.StrippedSequence, condition, temperature, pcr_minmax) %>% 
    `colnames<-`(c("uniqueID", "condition", "x", "y"))
  

  # create output directory python_import
  path_output <- paste(common_path,sm,"python_import",sep="/")

  if (!dir.exists(path_output)) {
  	dir.create(path_output,showWarnings = FALSE,recursive = TRUE)
  }

  # save tables 
  write.csv(export_table[export_table$uniqueID %in% FT_peptides, ], paste(path_output,"/",sm,"_FT.csv", sep=""))
  write.csv(export_table[export_table$uniqueID %in% HT_peptides, ], paste(path_output,"/",sm,"_HT.csv", sep=""))

  print(paste("Succesfully exported", sm, "files. The files can be found in", path_output))
  
}






