### Function to load the python data and organise it with a list
### Peptide matching functions (by position)

python_to_list <- function(sm,data,out_dir) {
  
  # generate filepaths
  mll_FT <- paste0(out_dir,"/",sm,"/python_output/MLL_",sm,"_FT.csv",sep = "")
  sol_FT <- paste0(out_dir,"/",sm,"/python_output/solution_",sm,"_FT.csv",sep = "")
  mll_HT <- paste0(out_dir,"/",sm,"/python_output/MLL_",sm,"_HT.csv",sep = "")
  sol_HT <- paste0(out_dir,"/",sm,"/python_output/solution_",sm,"_HT.csv",sep = "")


  #defined_list <- list()
  sm_list <- list()
  Import <- data %>% 
    dplyr::select(PG.ProteinAccessions, PEP.StrippedSequence) %>%
    `colnames<-`(c("protein", "peptide")) %>% distinct()
  
  #path <- paste(common_path, "0" , unique(path_file$Experiment[path_file$sm == sm]), "/", sm, "/", "python_output/", sep="" )
  HT_import_fit <- read.csv(sol_HT)
  HT_import_fit %<>% plyr::join(., Import, by="peptide") %>%
    mutate(trypticity = "HT")
  
  FT_import_fit <- read.csv(sol_FT)
  FT_import_fit %<>% plyr::join(., Import, by="peptide") %>%
    mutate(trypticity = "FT")
  
  sm_list[[sm]][["python_fit"]] <- rbind(HT_import_fit, FT_import_fit)
  
  HT_import_score <- read.csv(mll_HT)
  HT_import_score %<>% plyr::join(., Import, by="peptide") %>%
    mutate(trypticity = "HT")
  
  FT_import_score <- read.csv(mll_FT)
  FT_import_score %<>% plyr::join(., Import, by="peptide",match = "first") %>%
    mutate(trypticity = "FT")
  
  sm_list[[sm]][["python_score"]] <- rbind(HT_import_score, FT_import_score)
  
  return(sm_list)
}

