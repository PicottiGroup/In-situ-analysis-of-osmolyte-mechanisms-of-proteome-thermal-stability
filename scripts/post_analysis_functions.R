## Author: Monika Pepelnjak
## Date: 14.08.2020
## Post-python analysis
# python post analysis functions

library(seqinr)
library(magrittr)
library(dplyr)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(ggpubr)

plot_individual_peptides <- function(all_results, small_molecule, control = "Control", peptide_set) {
  
  peptide_set <- all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$peptide %in% peptide_set] %>% unique()
  
  split_peptides <- ceiling(length(peptide_set)/12)

    for(i in 1:split_peptides) {
    test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model")
    
    test$condition <- factor(test$condition, levels = c(small_molecule, control))
    levels(test$peptide) <- peptide_set
    
    sm_color <- c(Control = "#8D8D92", 
                  Proline = "#3C4285", 
                  Betaine = "#2C938B",
                  Trehalose = "#A84343",
                  Glucose = "#6B1818",
                  TMAO_2 = "#6B1818",
                  PK_0.01 = "#A84343",
                  TMAO = "#2C8AA4", 
                  TMAO1 = "#2C8AA4", 
                  ATP = "#3C4285", 
                  Low_PK = "#4AB698", 
                  Glycerol_2 = "#87667F")

    plot3 <- ggplot(test) +  
    geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
    scale_color_manual("Condition",values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F"),label = c("ATP","Control")) +
    scale_fill_manual("Condition",values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F"),label = c("ATP","Control"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide) +
      theme_minimal() +
      theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(), legend.position = "none") +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
      xlab("Temperature") +
      ylab("Scaled intensity")
  return(plot3)}

}

plot_single_peptide <- function(all_results, small_molecule, control = "Control", pep) {
  
  peptide_set <- all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$peptide %in% pep] %>% unique()
  
  test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide == pep, ] 
    
  test$condition <- factor(test$condition, levels = c(small_molecule, control))
  levels(test$peptide) <- pep
    
  sm_color <- c(Control = "#8D8D92", 
                Proline = "#3C4285", 
                Betaine = "#2C938B",
                Trehalose = "#A84343",
                Glucose = "#6B1818",
                TMAO_2 = "#6B1818",
                PK_0.01 = "#A84343",
                TMAO = "#2C8AA4", 
                TMAO1 = "#2C8AA4", 
                ATP = "#3C4285", 
                Low_PK = "#4AB698", 
                Glycerol_2 = "#87667F")
    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
      scale_color_manual("Condition",values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F","black"),label = c("ATP","Control","full model")) +
      scale_fill_manual("Condition",values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F","black"),label = c("ATP","Control","full model"))+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      theme_minimal() +
      theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(), legend.position = "none",
            axis.text = element_text(size = 18),
            text = element_text(size = 20)) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
      xlab("Temperature") +
      ylab("Scaled intensity")
    return(plot3)
}
  

plot_peptides_from_protein <- function(all_results, small_molecule, uniprot_code, n_split = 12) {
#  uniprot_code <- "P0A6Y8"
#  small_molecule <- "ATP"
  peptide_set <- significant_all_aascore[[small_molecule]][["scores"]][significant_all_aascore[[small_molecule]][["scores"]]$Protein == uniprot_code,] %>%
    .[order(.$start),] %>%
    mutate(peptide = factor(Peptide, levels=unique(.$Peptide)))
  
  split_peptides <- ceiling(length(peptide_set$peptide)/n_split)

    for(i in 1:split_peptides) {
    test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide_set$peptide[((n_split*(i-1))+1) : ((n_split*(i-1))+n_split)], ]
    
    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition)) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition)) +
      scale_color_manual(values = c("#00059C", "#8F8F8F", "#767791")) +
      scale_fill_manual(values = c("#00059C", "#8F8F8F", "#767791"))+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      facet_wrap(~peptide)
    print(plot3)}
  
}

plot_peptide_color_area <- function(all_results, small_molecule,
                                    control = "Control",peptide,area = "binding",
                                    col_area = "red",legend_pos = "none") {
  
  peptide_data <- all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide,] %>%
    unique() %>% 
    filter(condition != "full model")
  
  peptide_data_fitted.control <- peptide_data %>% filter(condition == "Control" & type == "fitted")
  peptide_data_fitted.sm <- peptide_data %>% filter(condition == "ATP"& type == "fitted")
  
  # binding
  if(area == "binding" | area == "Binding"){
    stop <- 1:which((peptide_data_fitted.sm$conflik_lower < peptide_data_fitted.control$conflik_upper) == TRUE)[1]  
  }
  # stabilisation
  if((area == "Stabilisation" | area == "stabilisation") & peptide_data$trypticity[1] == "FT"){
    stop <- which((peptide_data_fitted.sm$conflik_lower > peptide_data_fitted.control$conflik_upper) == TRUE)
  }
  if((area == "Stabilisation" | area == "stabilisation") & peptide_data$trypticity[1] == "HT"){
    stop <- which((peptide_data_fitted.control$conflik_lower > peptide_data_fitted.sm$conflik_upper) == TRUE)
  }
  # destabilisation
  if((area == "Destabilisation" | area == "destabilisation")){
    stop <- which((peptide_data_fitted.control$conflik_lower > peptide_data_fitted.sm$conflik_upper) == TRUE)
  }
  # aggregation
  if(area == "Aggregation" | area == "aggregation"){
    stop <- which((peptide_data_fitted.sm$conflik_lower > peptide_data_fitted.control$conflik_upper) == TRUE)
  }
  
  lower <- peptide_data_fitted.control$y[stop[1]:stop[length(stop)]]
  upper <- peptide_data_fitted.sm$y[stop[1]:stop[length(stop)]]
  
  data_ribbon <- data.frame(t = peptide_data_fitted.sm$t[stop[1]:stop[length(stop)]],
                            lower = lower,
                            upper = upper)
  
  peptide_data$condition <- factor(peptide_data$condition, levels = c(small_molecule, control))
  
  sm_color <- c(Control = "#8D8D92",
                Proline = "#3C4285", 
                Betaine = "#2C938B",
                Trehalose = "#A84343",
                Glucose = "#6B1818",
                TMAO_2 = "#6B1818",
                PK_0.01 = "#A84343",
                TMAO = "#2C8AA4", 
                TMAO1 = "#2C8AA4", 
                ATP = "#3C4285", 
                Low_PK = "#4AB698", 
                Glycerol_2 = "#87667F")
  
  plot <- ggplot(peptide_data) +  
    # color effect area
    geom_ribbon(data = data_ribbon,
                aes(x = t,ymin=lower,ymax = upper),alpha = 0.2,fill = col_area) +
    # CI
    geom_ribbon(data = peptide_data[peptide_data$type == "fitted",],
                aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.3) + 
    # data points
    geom_point(data = peptide_data[peptide_data$type == "measured",],
               aes(y=y, x=t, col=condition), size=2) +
    # fitted curves
    geom_line(data = peptide_data[peptide_data$type == "fitted",],
              aes(y=y, x=t, col=condition), lwd=1.2) +
    scale_color_manual("Condition",
                       values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F"),
                       label = c(small_molecule,control)) +
    scale_fill_manual("Condition",
                      values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F"),
                      label = c(small_molecule,control))+
    theme_minimal() +
    theme(axis.line = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), 
          legend.position = legend_pos) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
    xlab("Temperature") +
    ylab("Scaled intensity")
  return(plot)
}

plot_CI_difference <- function(pep){
  data <- significant_all$ATP$long %>% filter(peptide == pep) 
  plot1 <- data%>% 
    ggplot(.,aes(x = t,y = differences,col = goal)) + 
    geom_point(size = 2) + 
    geom_point(inherit.aes=FALSE,data = subset(data,apex == TRUE),
               aes(x = t,y = differences,fill = "red",size = 4),shape = 4,show.legend = FALSE) +
    theme_classic() +  
    ylab(expression(paste("Difference of \nnon-overlapping CI"))) + 
    xlab("Temperature") + 
    theme(text = element_text(size = 10),legend.position = "none")
  return(plot1)
}

plot_idv_peptide <- function(all_results, small_molecule, pep) {
  
  peptide_set <- all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$peptide %in% pep] %>% unique()
  
  test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide == pep, ] %>%
    filter(condition == small_molecule)
  
  test$condition <- factor(test$condition, levels = c(small_molecule))
  levels(test$peptide) <- pep
  
  sm_color <- c(Control = "#8D8D92", 
                Proline = "#3C4285", 
                Betaine = "#2C938B",
                Trehalose = "#A84343",
                Glucose = "#6B1818",
                TMAO_2 = "#6B1818",
                PK_0.01 = "#A84343",
                TMAO = "#2C8AA4", 
                TMAO1 = "#2C8AA4", 
                ATP = "#3C4285", 
                Low_PK = "#4AB698", 
                Glycerol_2 = "#87667F")
  
  test_fitted <- test %>% filter(type == "fitted")
  data_arrow <- test_fitted[test_fitted$y == min(test_fitted$y) | test_fitted$y == max(test_fitted$y),]
  plot3 <- ggplot(test) +  
    geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
    scale_color_manual("Condition",values = c(sm_color[names(sm_color) == small_molecule]),label = c(small_molecule)) +
    scale_fill_manual("Condition",values = c(sm_color[names(sm_color) == small_molecule]),label = c(small_molecule))+
    theme_minimal() +
    theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.position = "none") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
    geom_segment(x = 46.75,y=1.014879,xend = 46.75,yend=-0.02676764,arrow = arrow(length= unit(0.03,"npc"),ends = "both"),color = "green") +
    geom_segment(x = 45.75,xend = 47.75,y = 1.014879,yend = 1.014879,linetype = "dashed",color = "grey") +
    geom_segment(x = 45.75,xend = 47.75,y = -0.02676764,yend = -0.02676764,linetype = "dashed",color = "grey") +
    geom_segment(x = 65.64,xend = 66.64,y = -0.02676764,yend = -0.02676764,linetype = "dashed",color = "grey") +
    geom_segment(x = 37,xend = 37,y = 0.942107685,yend = 0.005768004,arrow = arrow(length= unit(0.03,"npc"),ends = "both"),col = "red") +
    geom_segment(x = 36,xend = 38,y = 0.005768004,yend = 0.005768004,linetype = "dashed",color = "red") +
    geom_segment(x = 75,xend = 77,y = 0.005768004,yend = 0.005768004,linetype = "dashed",color = "red") +
    xlab("Temperature") +
    ylab("Scaled intensity")
  return(plot3)
}

plot_individual_peptides_full_model <- function(all_results, small_molecule, control = "Control", peptide_set) {
  
  peptide_set <- all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$peptide %in% peptide_set] %>% unique()
  
  split_peptides <- ceiling(length(peptide_set)/12)
  
  for(i in 1:split_peptides) {
    test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ]
      #filter(condition == "full model")
    
    test$condition <- factor(test$condition, levels =c(small_molecule, control,"full model"))
    levels(test$peptide) <- peptide_set
    
    sm_color <- c(Control = "#8D8D92", 
                  Proline = "#3C4285", 
                  Betaine = "#2C938B",
                  Trehalose = "#A84343",
                  Glucose = "#6B1818",
                  TMAO_2 = "#6B1818",
                  PK_0.01 = "#A84343",
                  TMAO = "#2C8AA4", 
                  TMAO1 = "#2C8AA4", 
                  ATP = "#3C4285", 
                  Low_PK = "#4AB698", 
                  Glycerol_2 = "#87667F")
    
fig <-     ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
      scale_color_manual("Condition",values = c("#3C4285", "black","#8F8F8F"),label = c("ATP","Control","null model")) +
      scale_fill_manual("Condition",values = c("#3C4285", "black","#8F8F8F"),label = c("ATP","Control","null model"))+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      
      theme_minimal() +
      theme(axis.line = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 18),
            panel.grid.major = element_blank()) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
      xlab("Temperature") +
      ylab("Scaled intensity")
    return(plot3)}
  
}

