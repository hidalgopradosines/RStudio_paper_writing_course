# RStudio paper writing course - excercises
# Gaspar Jekely, 2024

# sourcing, installing and loading packages -------------------------------

# installing packages
install.packages("tidyverse")

# source multiple packages from one file
source("analysis/scripts/packages_and_functions.R") #we are sourcing/loading all the packages we need during the analysis

#check your working dir (should be the .rproject dir)
getwd()

# in case you used setwd() in a script (you should not) you can find your Rproject directory with here()
here::here()

#test, but you should not do this
setwd("~")
getwd()

#get back to Rproj dir
setwd(here::here())
getwd()
# all good again

list.files()
list.files("analysis/data")
#loading packages
library(tidyverse)
library(png)


# share session info ------------------------------------------------------

#save session info and Rstudio version info for reproducibility
sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt") #captures the output in a file

# load data ---------------------------------------------------------------

data_Jose <- readxl::read_excel("analysis/data/data - José - March 2024.xlsx")

# alternatively, use the rio package with import() to automatically recognise format and import
data_Jose <- rio::import("analysis/data/data - José - March 2024.xlsx")

head(data_Jose)
glimpse(data_Jose)
str(data_Jose)
summary(data_Jose)
data_Jose


# rio: A Swiss-Army Knife for Data I/O  -----------------------------------

#export as csv
rio::export(data_Jose, "analysis/data/data_Jose_March2024.csv")
#export as compressed csv
export(data_Jose, "analysis/data/data_Jose_March2024.csv.zip")

#convert file formats
rio::convert("analysis/data/data - José - March 2024.xlsx",
"analysis/data/data_Jose_March2024.csv")


# another data file -------------------------------------------------------

data_Ashwini <- readxl::read_excel("analysis/data/24_01_22qpcr_1 - ash pal.xls")
#or with rio
data_Ashwini <- rio::import("analysis/data/24_01_22qpcr_1 - ash pal.xls")

head(data_Ashwini)
glimpse(data_Ashwini)
str(data_Ashwini)
summary(data_Ashwini)


# a tidy dataset ----------------------------------------------------------

head(iris)
vignette("tibble")

# overwrite gene names ----------------------------------------------------

Gene <- c(rep("Nanog", 15), rep("oct4", 15), rep("sox2", 15),
  rep("Nestin", 15), rep("pax6", 15), rep("Foxg1", 15),
  rep("GAPDH", 15))
Gene

data_Ashwini$Gene <- Gene
head(data_Ashwini)

# Select only relevant columns and clean up names -------------

data_Ashwini_sel <- data_Ashwini %>%
  select(1:6) %>%
  janitor::clean_names()
data_Ashwini_sel

# Add mean and SD columns with group_by() and mutate() --------

data_Ashwini_sel_M_SD <- data_Ashwini_sel %>%
  group_by(gene, days) %>%
  mutate(mean2dct = mean(x2_dct)) %>%
  mutate(sd2dct = sd(x2_dct))
data_Ashwini_sel_M_SD

# Change data type -----------

data_Ashwini_sel_M_SD <- data_Ashwini_sel_M_SD %>%
  mutate(ct_value = as.double(ct_value))
data_Ashwini_sel_M_SD


# tidying data ------------------------------------------------------------

data_Syn <- read_csv("analysis/data/a-Syn-Data.csv")
data_Syn

# rename
data_Syn_clean <- data_Syn  %>%
  rename_with(~ gsub("_", "-", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("...", "_", .x, fixed = TRUE))

tb_syn <- data_Syn_clean |>
  pivot_longer(matches("aSyn"), 
               names_to = c("chemistry", "sample"), 
               names_sep = "_",
               values_to = "fluorescence")

plot_syn <- tb_syn %>%
  ggplot(aes(x = Time, y = fluorescence, color = chemistry)) +
  geom_smooth(
    method = 'loess', span = 0.1
    ) +
  theme_minimal()
plot_syn

plot_syn +
  annotate("segment", x = 20, xend = 50, y = 1, yend = 1, linewidth = 1)+
  annotate("text", x = 34, y = 300, label = "30 sec", size = 3)

# Aesthetics, plot types and themes ------------
head(iris)
iris %>%  
  ggplot(aes(x = Petal.Length, y = Sepal.Width, color = Species)) +
  geom_boxplot(notch = TRUE) + #with the notch we already do a statistical test. If the notches do not overlap, they are statistically different
  theme_minimal()
#also possible to write it like this:
#ggplot(iris, aes(x = Petal.Length, y = Sepal.Width, color = Species)) +
#  geom_boxplot(notch = TRUE) +
#  theme_minimal()

iris %>%
  ggplot(aes(x = Sepal.Width)) +
  geom_histogram()

iris %>%  
  ggplot(aes(
    x = Petal.Length, y = Sepal.Width, 
    color = Species, size = Sepal.Length)
    ) +
  geom_point() +
  facet_wrap(~Species) #separate it into different plots for each species. If we leave this out, we plot all species in a single plot


# Plot data - Jose ---------------
data_Jose
plot_Jose1 <- data_Jose %>%
  ggplot(aes(x = genotype, y = length, fill = factor(Treatment, level=c('Control', 'ABA', 'Sulfate')), na.rm = TRUE)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#cccccc")) +
  guides(fill = guide_legend(title = "Treatment")) +
  coord_flip() + #flips the plot from vertical to horizontal
  scale_y_log10() #to have log10 scale in the y axis

plot_Jose1  

plot_Jose2 <- data_Jose %>%
  ggplot(aes(x = genotype, y = length, fill = factor(Treatment, level=c('Control', 'ABA', 'Sulfate')), na.rm = TRUE)) +
  geom_violin() +
  geom_point( position=position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.5, size = 0.4) +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#aaaaaa", "#dddddd")) +
  guides(fill = guide_legend(title = "Treatment")) +
  coord_flip() +
  scale_y_log10() 
plot_Jose2


# Plot data - Ashwini ----------------------

data_Ashwini_sel_M_SD %>%
  group_by(gene) %>%
  ggplot(aes(x = days, y = dt_ct, fill = gene )) +
  geom_boxplot()

data_Ashwini_sel_M_SD %>%
  ggplot(aes(x = ct_value)) +
  geom_histogram()

plot_Ashwini_ct <- data_Ashwini_sel_M_SD %>%
  group_by(gene) %>%
  ggplot(aes(x = gene, y = ct_value, fill = gene )) +
  geom_boxplot(na.rm = TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

plot_Ashwini_ct

# Plot the Synuclein data ----------
tb_syn
plot_syn <- tb_syn %>%
  ggplot(aes(x = Time, y = fluorescence, color = fluorescence)) +
  geom_smooth(method = 'loess') +
  theme_minimal()
plot_syn

# Read and preview data 3 -------------- #This data is super messy

data_Anchel <- rio::import("analysis/data/240323 CIN Exp278 reporter assay - Anchel.xlsx")
head(data_Anchel)
glimpse(data_Anchel)
str(data_Anchel)
summary(data_Anchel)

# Save tidy data as source data for the plot/figure/paper -----------

write_csv2(data_Ashwini_sel_M_SD, "manuscript/source_data/FigureX_Ashwini_source_data1.csv")
#now if we go inside the project to manuscrip --> source data, you can find the file. You should do this for every plot. Good journals ask for the source data associated to a plot

# check
read_csv2("manuscript/source_data/FigureX_Ashwini_source_data1.csv")

# Format plots with predefined complete ggplot2 themes ------------ #Good to use the same theme for all plots in your publication 
plot_Jose1
plot_Jose1 +
  theme_dark()
plot_Jose1 +
  theme_bw()
plot_Jose1 +
  theme_linedraw()

plot_Jose2 +
  theme_classic()
plot_Jose2 +
  theme_minimal()
plot_Jose2 +
  theme_light()


# Format plots with a common custom theme() -------------

args(theme) #list of all things that you can modify in your plot (e.g. font, spacing, size of letters...); everything about the graphical/aesthetic aspect

theme_plots <- theme_minimal() + #we can modify the theme_mininmal() with + theme(different arguments of things we want to modify). We can save all the modifications/settings we like under the variable theme_plots. Then, we can add + theme_plots to all our plots so that they are all the same (below there are some examples)
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.size = unit(7, "mm"),
    legend.title.position = "top",
    legend.background = element_rect(color = "grey"),
    plot.title.position = "panel"
  )

plot_Ashwini_ct <- plot_Ashwini_ct +
  theme_plots
plot_Ashwini_ct

plot_Jose1 <- plot_Jose1 +
  theme_plots
plot_Jose1

plot_Jose1 <- plot_Jose1 +
  theme_plots +
  axis.text.x(element.text(angle=90))
plot_Jose1

plot_Jose2 <- plot_Jose2 +
  theme_plots
plot_Jose2

plot_syn <- plot_syn +
  theme_plots
plot_syn

# Optional - save plots as png--------------- #png have a fixed size. Better to use ggplot directly to assemble the different plots into a figure

ggsave( "analysis/pictures/plot_Jose1a.png",
  limitsize = FALSE,
  units = c("px"), plot_Jose1,
  width = 3400, height = 1600 #background will be transparent
  )

# save in a different size
ggsave( "analysis/pictures/plot_Jose1b.png",
  limitsize = FALSE,
  units = c("px"), plot_Jose2,
  width = 2400, height = 2000, bg = "white" # to have the background white
  )

ggsave(
  "analysis/pictures/synuclein_plot.png", plot_syn, 
  bg = "white"
  )

# Assemble figure with cowplot and patchwork --------------

#read images

img1 <- readPNG("analysis/pictures/plot_Jose1a.png")
img2 <- readPNG("analysis/pictures/plot_Jose1b.png")

#convert to panels
panel_JoseA <- ggdraw() + draw_image(img1) #converts image into a panel
panel_JoseB <- ggdraw() + draw_image(img2)

#define layout with textual representation #define the layout for the figure
layout <- "
AB
CD"
#other option would be
layout <- "
ABCD"
#other options are
layout <- "
ABCDEF
GHHHHI
JHHHHK
LHHHHM
NOPQRS"
layout <- "
ABCDEF
AHHHHI
JHHHHK
LHHHHM
NOPQRS"
#you can define empty space with #
#if you run out of alphabetic letters, you can use capital and small letters
#check cheatsheet/tutorial of patchwork

#assemble multipanel figure based on layout
Figure_Jose <- panel_JoseA + panel_JoseB + plot_Jose1 + plot_Jose2 +
  plot_layout(design = layout, heights = c(1, 1)) + #heights of the pannels
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain')) #tags of the panels
#you can also add insets to your figure 

layout <- "
AB
##
CD"
Figure_Jose <- panel_JoseA + panel_JoseB + plot_Jose1 + plot_Jose2 +
  plot_layout(design = layout, heights = c(1, 0.05, 1)) + #heights of the panels
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))


#save figure as png and pdf
ggsave(
  "manuscript/figures/Figure_Jose.png", limitsize = FALSE, 
  units = c("px"), Figure_Jose, width = 4000, height = 1600,
  bg = "white"
  ) #you cantry to save the image with different dimensions. This is quite critical for your final figure

ggsave(
  "manuscript/figures/Figure_Jose.pdf", limitsize = FALSE, 
  units = c("px"), Figure_Jose, width = 4000, height = 1600
  )

image_read("manuscript/figures/Figure_Jose.png")

#it is also possible to instead of first saving the plots as png and then read them, use the variable you saved the plot to to make the panels directly
#because although we used the same theme for all the plots, when we saved them as png we saved them with different dimensions, and then it is very difficult to adjust the exact same font size for all plots 
#we can also save all plots with the same theme and same dimensions?
Figure_Jose <- plot_Jose1 + plot_Jose2 +
  plot_layout(design = layout, heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))



# Annotating a ggplot object ----------------

plot_syn_ann <- plot_syn +
  annotate("segment", x = 20, xend = 50, y = 1, yend = 1, linewidth = 1)+
  annotate("text", x = 34, y = 300, label = "30 sec", size = 3)
plot_syn_ann

plot_Jose1 +
  annotate("segment", x = 1, xend = 2, y = 1, linewidth = 1)

# Annotating an image ------------

#read images
#first you do some image processing in imageJ and the final image we import here
img_INNOS <- magick::image_read("analysis/pictures/INNOS_synapses.png")

#define arrow endpoints 
arrow <- data.frame(x1 = 0.95, x2 = 0.95, y1 = 0.8, y2 = 0.9)

#add text labels
panel_INNOS <- ggdraw() + 
  draw_image(img_INNOS) +
  draw_label("INNOS", x = 0.5, y = 0.99, size = 10) +
  draw_label("NS plexus", x = 0.485, y = 0.59, size = 8) +
  draw_label("outgoing", x = 0.9, y = 0.45, size = 10, color='#E29F00') +
  draw_label("incoming", x = 0.89, y = 0.5, size = 10, color='#0072B2') +
  draw_label("D", x = 0.95, y = 0.93, size = 6) +
  draw_label("V", x = 0.95, y = 0.77, size = 6) +
  draw_label("*", x = 0.5, y = 0.29, color='black',size = 18,fontface='plain') +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = arrow, 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.2)

#define layout
layout <- "AB"

#assemble multipanel figure based on layout
Figure_INNOS <- plot_syn_ann + panel_INNOS +
  plot_layout(design = layout, widths = c(2, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

#save figure as png
ggsave(
  "manuscript/figures/Figure_INNOS.png", limitsize = FALSE,
  units = c("px"), Figure_INNOS, 
  width = 3000, height = 1000,
  bg = "white"
  )

#save figure as pdf
ggsave(
  "manuscript/figures/Figure_IHC.pdf", limitsize = FALSE, 
  units = c("px"), Figure_INNOS, width = 3000, height = 1000
  )

image_read("manuscript/figures/Figure_IHC.png")


# Adding consistent scale bars -----------------

#read images and make annotated panel
panel_NOS2d_HCR <- ggdraw() + draw_image(readPNG("analysis/pictures/HCR-IHC_51_AP_NOS_actub_56um.png")) +
  draw_label("in situ HCR", x = 0.3, y = 0.99, size = 10) + #title of the panel
  draw_label("NOS", x = 0.12, y = 0.9, color="magenta", size = 11, fontface="italic") + #this adds a label
  draw_label("acTub", x = 0.36, y = 0.9, color="green", size = 11, fontface="plain") + #this adds another label
  draw_line(x = c(0.1, 0.46), y = c(0.08, 0.08), color = "white", size = 0.5) + #this makes the scalebar
  draw_label(expression(paste("20 ", mu, "m")), x = 0.28, y = 0.11, color = "white", size = 8) #label of the scalebar
  #expression and mu is to get the micro letter (greek) 
panel_NIT_HCR <- ggdraw() + draw_image(readPNG("analysis/pictures/HCR_72_AP_NIT_94um.png")) +
  draw_label("transgene + IHC", x = 0.5, y = 0.99, size = 10) +
  draw_label("NOSp::palmi-3xHA", x = 0.34, y = 0.9, color="magenta", size = 10, fontface="plain") +
  draw_label("acTub", x = 0.8, y = 0.9, color="green", size = 10, fontface="plain") +
  draw_line(x = c(0.1, 0.31), y = c(0.08, 0.08), color = "white", size = 0.5) 
panel_NIT_HCR

#magick package allows to work with TIFs instead of pngs
?magick

# introduce gaps in layout --------------

layout <- "A#B"

#assemble multipanel figure based on layout
Figure_scalebars <- panel_NOS2d_HCR + panel_NIT_HCR +
  plot_layout(design = layout, widths = c(1, 0.01, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

#save figure as png
ggsave(
  "manuscript/figures/Figure_scalebars.png",
  units = c("px"), Figure_scalebars, 
  width = 1700, height = 940, bg = "white" #the size of the figure is the width of the two confocal images (800 pixel each) plus a little gap and the height of the two confocal images (800 pixel each) plus the text
  )

ggsave(
  "manuscript/figures/Figure_scalebars.pdf",
  units = c("px"), Figure_scalebars, 
  width = 1700, height = 940
  )
image_read("manuscript/figures/Figure_scalebars.png")
image_read("manuscript/figures/Figure_scalebars.pdf")

# Fine-tuning figure size and gaps ----------------

# Excercises
# - save the figure in different sizes 
ggsave(
  "manuscript/figures/Figure_scalebars.png",
  units = c("px"), Figure_scalebars, 
  width = 2000, height = 2000
)
# - introduce gap with # into layout, also need to define width of gap as say 0.05
layout <- "A#B"
Figure_scalebars <- panel_NOS2d_HCR + panel_NIT_HCR +
  plot_layout(design = layout, widths = c(1, 0.05, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave(
  "manuscript/figures/Figure_scalebars.png",
  units = c("px"), Figure_scalebars, 
  width = 1700, height = 940
)
# - change position of scalebar and scalebar legend
panel_NOS2d_HCR <- ggdraw() + draw_image(readPNG("analysis/pictures/HCR-IHC_51_AP_NOS_actub_56um.png")) +
  draw_label("in situ HCR", x = 0.3, y = 0.99, size = 10) + 
  draw_label("NOS", x = 0.12, y = 0.9, color="magenta", size = 11, fontface="italic") + #this adds a label
  draw_label("acTub", x = 0.36, y = 0.9, color="green", size = 11, fontface="plain") + #this adds another label
  draw_line(x = c(0.1, 0.46), y = c(1, 0.08), color = "white", size = 0.5) + #this makes the scalebar
  draw_label(expression(paste("20 ", mu, "m")), x = 0.1, y = 0.5, color = "white", size = 8) #label of the scalebar

panel_NIT_HCR <- ggdraw() + draw_image(readPNG("analysis/pictures/HCR_72_AP_NIT_94um.png")) +
  draw_label("transgene + IHC", x = 0.5, y = 0.99, size = 10) +
  draw_label("NOSp::palmi-3xHA", x = 0.34, y = 0.9, color="magenta", size = 10, fontface="plain") +
  draw_label("acTub", x = 0.8, y = 0.9, color="green", size = 10, fontface="plain") +
  draw_line(x = c(0.1, 0.31), y = c(0.5, 0.5), color = "white", size = 0.5) 

Figure_scalebars <- panel_NOS2d_HCR + panel_NIT_HCR +
  plot_layout(design = layout, widths = c(1, 0.05, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))
ggsave(
  "manuscript/figures/Figure_scalebars.png",
  units = c("px"), Figure_scalebars, 
  width = 1700, height = 940
)
#no gap in layout
layout1 <- "AB"

#assemble multipanel figure based on layout
Figure_scalebars <- panel_NOS2d_HCR + panel_NIT_HCR +
  plot_layout(design = layout1, widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

#save figure as png
ggsave(
  "manuscript/figures/Figure_scalebars_no_gap.png", 
  limitsize = FALSE,
  units = c("px"), Figure_scalebars, 
  width = 1700, height = 940,
  bg = "white"
  )

image_read("manuscript/figures/Figure_scalebars_no_gap.png")

#introduce gap in layout
layout2 <- "A#B"

#assemble multipanel figure based on layout
Figure_scalebars <- panel_NOS2d_HCR + panel_NIT_HCR +
  plot_layout(design = layout2, widths = c(1, 0.03, 1)) +
  plot_annotation(tag_levels = list(c('A','','B',''))) & 
  theme(plot.tag = element_text(size = 12, face='plain'))

#save figure as png
ggsave(
  "manuscript/figures/Figure_scalebars_gap.png", 
  limitsize = FALSE,
  units = c("px"), Figure_scalebars, 
  width = 1700, height = 940,
  bg = "white"
  )

image_read("manuscript/figures/Figure_scalebars_gap.png")

# More complex figure layouts --------------

#read images and make annotated panel
panel_Platy <- ggdraw() + draw_image(readPNG("analysis/pictures/Platynereis_SEM_inverted_nolabel.png"))
panel_NOS <- ggdraw() + draw_image(readPNG("analysis/pictures/HCR-IHC_51_AP_NOS_actub_56um.png"))
panel_FVRI <- ggdraw() + draw_image(readPNG("analysis/pictures/FVRIa_rhoPhall_31h_200um.png"))
panel_Jose <- ggdraw() + draw_image(readPNG("analysis/pictures/plot_Jose1b.png"))
panel_INNOS <- ggdraw() + draw_image(readPNG("analysis/pictures/INNOS_synapses.png"))
panel_NIT <- ggdraw() + draw_image(readPNG("analysis/pictures/IHC_55_AP_NITGC2_actub_61um.png"))
panel_DAF <- ggdraw() + draw_image(readPNG("analysis/pictures/DAFFM.png"))
panel_model <- ggdraw() + draw_image(readPNG("analysis/pictures/Magnitude_model_cPRC.png"))

#introduce gap in layout
layout <- "
AAAABBBBCCCC
AAAABBBBDDDD
############
EEEFFFGGGHHH
"

#assemble multipanel figure based on layout
Figure_complex <- panel_Platy + panel_FVRI +  panel_NOS + 
  panel_NIT +
  panel_INNOS + panel_Jose + panel_DAF +
  panel_model +
  plot_layout(design = layout, heights = c(1, 1, 0.05, 2)) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

#save figure as png
ggsave(
  "manuscript/figures/Figure_complex.png",
  units = c("px"), Figure_complex, 
  width = 2600, height = 1700, bg = "white"
  )

image_read("manuscript/figures/Figure_complex.png")


# Image saved with defined resolution (dpi) -------------------------------


#save figure as png at 300 dpi
ggsave(
  "manuscript/figures/Figure_complex_300dpi.png",
  units = c("cm"), Figure_complex, 
  width = 22, height = 13, dpi = 300, bg = "white"
  )
image_read("manuscript/figures/Figure_complex_300dpi.png")


# Read tif files -------------
install.packages("magick") #also allows to read svg files
?magick
?image_read()

# Statistical comparisons -------------


data_Jose %>%
  ggplot(aes(x = genotype, y = length, fill = factor(Treatment, level=c('Control', 'ABA', 'Sulfate')), na.rm = TRUE)) +
  geom_violin() +
  geom_point( position=position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.5, size = 0.4) +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#aaaaaa", "#dddddd")) +
  guides(fill = guide_legend(title = "Treatment")) +
  coord_flip() +
  scale_y_log10()


