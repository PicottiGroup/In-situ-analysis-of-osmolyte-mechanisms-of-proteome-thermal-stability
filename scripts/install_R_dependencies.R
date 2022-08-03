# necessary packages
packages <- c("plotly","ggplot2",
	"magrittr","dplyr",
	"data.table","ggpubr",
	"spatstat","seqinr",
	"fclust") 

# check for packages, which are not yet installed
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]


# install packages
if(length(new.packages)!= 0){
	print("Installing following packages:")
	new.packages
	install.packages(new.packages)
} else{
	print("Packages are already installed")
}

