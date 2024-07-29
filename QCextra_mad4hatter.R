#############################
#### Mad4hatter QC extra ####
#############################

library(ggplot2)
library(tidyverse)
library(data.table)
library(forcats)
library(ggrepel)
library(grid)
library(writexl)
library(dplyr)
library(readxl)
library(ggpubr) # for regression lines
library(RColorBrewer)
library(reshape)
library(gridExtra)

# Set working directory
setwd("~/Escritorio/00-GenMoz/03-Data_analysis/QC_EXTRA/MDACD_NX02/")

# Read input files, clean and merge by SampleID
# Input 1: Sample Number of Loci with 100 Reads or More per Pool
ampliok <- read.table("MDACD_NX02_Amplicons.txt", header = TRUE)
ampliok <- ampliok %>%
  mutate(key = substr(SampleID, 2, 15),
         Pool_all = rowSums(across(c(2, 4, 6)))) %>%
  select(SampleID, 2, 4, 6, 7, 8)
# Input 2:SEQ excel file, it needs library information for every sample
sheet <- read_excel("MDACD_NextSEQ02_100724.xlsx", sheet="Balancing") %>%
select(1, 3, 4, 8, 10) %>%
  mutate(`Sample ID` = gsub("\\.", "_", `Sample ID`),
         key = `Sample ID`)
# Merge data
df1 <- merge(sheet, ampliok, by = "key") %>%
  mutate(Library = as.factor(Library))

#################################################################################################
# Parasitemia and libraries performance
#################################################################################################
# Density plots
# Parasitemia
d1 <- ggplot(data = df1, aes(x = Parasitemia, color = Library)) +
  geom_density() +
  labs(color = "Library plate") +# Change legend title
  theme_light() +
  ggtitle("MDACD_NX02") +
  theme(
    legend.title = element_text(size = 8),       
    legend.text = element_text(size = 6),        
    axis.title = element_text(size = 10),        
    axis.text = element_text(size = 8),          
  )
# Performance
d2 <- ggplot(data = df1, aes(x = Pool_all, color = Library)) +
  geom_density() +
  labs(color = "Library plate") +# Change legend title
  theme_light() +
  xlab("Loci > 100 reads") +
  theme(
    legend.title = element_text(size = 8),       
    legend.text = element_text(size = 6),        
    axis.title = element_text(size = 10),        
    axis.text = element_text(size = 8),          
  )
# Geom point
d3 <- ggplot(df1, aes(x = Parasitemia, y = Pool_all, color = Library)) +
  geom_point() +
  labs(color = "Library plate") +
  theme_light() +
  ylab("Loci >100 reads") +
  theme(
    legend.title = element_text(size = 8),       
    legend.text = element_text(size = 6),        
    axis.title = element_text(size = 10),        
    axis.text = element_text(size = 8),          
    plot.title = element_text(size = 10)         # Plot title size
  )
# Boxplot
d4 <- ggplot(df1, aes(x = Library, y = Pool_all)) +
  geom_boxplot() +
  geom_point(shape=21, position=position_jitter(0.2), aes(fill = Parasitemia), size =1) +
  scale_fill_gradientn(
    colours = colorRampPalette(c('#d64138', '#FBF719', '#5fc435'))(100)
  ) +
  xlab("Library plate") +
  ylab("Loci >100 reads") +
  theme_light() +
  theme(
    legend.title = element_text(size = 8),       
    legend.text = element_text(size = 6),        
    axis.title = element_text(size = 10),        
    axis.text = element_text(size = 8),          
  )

# Export and combine plots in a PDF
d5 <- grid.arrange(d1, d2, d3, d4, nrow = 2, ncol =2)
ggsave(filename = "MDACD_NX02_parasitemia_performance.pdf", plot = d5, width = 10, height = 10)

#################################################################################################
# Heatmaps in 96 well plates
#################################################################################################

# Prepare input
df2 <- df1 %>% pivot_longer(cols = starts_with("Pool"), names_to = 'Primer_pools', values_to = 'POkAmp')
# Function to extract letters
extract_letters <- function(x) {
  str_extract(x, "[A-Za-z]+")
}
# Function to extract numbers
extract_numbers <- function(x) {
  str_extract(x, "\\d+")
}
# Extract letters and numbers from Well column
df3 <- df2 %>%
  mutate(
    Row = as.factor(extract_letters(Well)),
    Column = as.factor(extract_numbers(Well))
  )
# Function to extract subsets of library experiments from the sequencing run
crear_subsets <- function(data, columna, sufijo) {
  valores_unicos <- unique(data[[columna]])
  for (valor in valores_unicos) {
    subset_name <- paste0(valor, "_", sufijo)
    subset_data <- subset(data, data[[columna]] == valor)
    assign(subset_name, subset_data, envir = .GlobalEnv)
  }
}
# Function to plot heatmap of performance in the plate
heatmap_data <- function(data, title, low_color, high_color) {
  ggplot(data, aes(x = Column, y = Row, fill = POkAmp)) +
    geom_tile(color = "black", lwd = 0.2, linetype = 1) +
    scale_y_discrete(limits = rev) +
    scale_x_discrete(limits = 1:12) +
    coord_fixed() +
    scale_fill_gradient(low = low_color, high = high_color) +
    ggtitle(title) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 10)
    )
}
# Run heatmap_data function for every library plate/experiment
lib_plate <- "4"
Lib <- df3[df3$Library %in% lib_plate,]
crear_subsets(Lib, "Primer_pools", lib_plate)
p1 <- heatmap_data(Pool_all_4, "MDACD_NX02_LIB04", "#070738", "steelblue2")
p2 <- heatmap_data(Pool_1A_4, "Pool1A", "goldenrod4", "goldenrod1")
p3 <- heatmap_data(Pool_1B_4, "Pool5", "darkgreen", "darkolivegreen2")
p4 <- heatmap_data(Pool_2_4, "Pool2", "darkorange4", "darkorange")
# Save every plot in a different object (p5, p6, p7, p8...)
p8 <- grid.arrange(p1, p2, p3, p4, nrow = 2, ncol =2)
# Export as a PDF file
ggsave(filename = "MDACD_NX02_Heatmap_LIB04.pdf", plot = p8, width = 10, height = 10)

#########################################################################################
# COMBINE AND EXPORT ALL PLOTS #
#########################################################################################
pdf("MDACD_NX02_QC.pdf")

# List of plot objects to arrange
plots <- list(d5, p5, p6, p7, p8)
# Loop through the list of plots in batches and arrange them

for (plot in plots) {
  grid.arrange(plot)
}

# Close the PDF device
dev.off()
