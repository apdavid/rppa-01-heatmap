library(limma)
library(lumi)
library(gplots)
library(readxl)
library(tidyverse)
library(pvclust)

crossref <- read.csv("paths.csv")

# read in RPPA data saved as csv from excel
# CHANGE THIS TO RUN LOCALLY
# allData <- read.csv("~/BTSync/MASS_SPEC/CrazyRPPA/RPPA_EGFR_IGF1R.csv")

allData <- read.csv("RPPA_EGFR_IGF1R.csv")

# grab only "targets" data. Limma uses the language of "targets" and "expression matrix". "targets"
# refers to annotation data that describes the experimental design. The "expression matrix" refers
# to the expression data for which the linear model will be fit. 

targets <- allData[-c(1:7),c(1:8)] %>%
  setNames(t(allData[7,c(1:8)]))

targets <- targets[,c(3:8)]

# limma relies heavily on factors for experimental design, so levels must be explictly set

# grabbing expression matrix data here
rppa.raw <- as.matrix(allData[-c(1:7),c(12:169)])
# limma expects an expression matrix to have genes as rows and samples as columns, so we take the transverse
colnames(rppa.raw) <- t(allData[7,c(12:169)])
rppa.trnsv <- t(rppa.raw)
class(rppa.trnsv) <- "numeric"

# this will induce NAs ("ND" is not able to be cast as a numeric)

# "ND"  means unsure of value 
# "0" means the analyte is truly not there0

# limma expects log transformed data
rppa.BC <- log2(rppa.trnsv)

# turn matrix into object that limma needs
row.names(rppa.BC)[150] <- paste(row.names(rppa.BC)[150], "ALT", " ")

# Rather than trying to run all analysis at once, we are running two independent paired analyses
# One to extract the consistent response to stimulation
# One to extract the consistent response to treatment with either OSI or BMS

# Make a targets data frame with data for the treatment effects
treat_targets <- targets %>%
  # filter() assumes "AND" behavior when you use ","
  # Any data we grab must satisfy ALL filters specified in the command
  filter(`Cell Line` != "SCC25GR1", # remove unwanted cell line
         Stimulus %in% c("IGF"), # grab stimulated data only 
         `Inhibitor - EGFR` == "0", # grab data that wasn't treated with EGFR inhibitor
         `Inhibitor - IGF1R` %in% c("0", "BMS", "OSI")) %>% # grab IGF1R inhibitor data
  # As mentioned earlier, limma needs factors to behave properly, doubly so for paired analysis
  mutate(Treatment = factor(`Inhibitor - IGF1R`, levels = c("0", "BMS", "OSI")),
         `Cell Line` = factor(`Cell Line`, levels = unique(`Cell Line`))
  )

# Make a targets data frame with data for the stimulation effect, much as the same as above
stim_targets <- targets %>%
  filter(`Cell Line` != "SCC25GR1",
         Stimulus %in% c("0", "IGF"),
         `Inhibitor - EGFR` == "0",
         `Inhibitor - IGF1R` == "0") %>%
  mutate(Stimulus = factor(Stimulus, levels = c("0", "IGF")),
         `Cell Line` = factor(`Cell Line`, levels = unique(`Cell Line`))
  )

# subset appropriate experimental data for each analysis, which returns the index that satisfy the constraints
treat_rppa <- rppa.BC[,which((targets$`Cell Line` != "SCC25GR1" &
                               targets$Stimulus %in% c("IGF") &
                               targets$`Inhibitor - EGFR` == "0" &
                               targets$`Inhibitor - IGF1R` %in% c("0", "BMS", "OSI")) == T)]

stim_rppa <- rppa.BC[,which((targets$`Cell Line` != "SCC25GR1" &
                               targets$Stimulus %in% c("0", "IGF") &
                               targets$`Inhibitor - EGFR` == "0" &
                               targets$`Inhibitor - IGF1R` == "0") == T)]

# Paired analysis for limma (see limma user manual section 9.4 for details)
stim_design <- model.matrix(~stim_targets$`Cell Line` + stim_targets$Stimulus)
stim_fit <- lmFit(stim_rppa, stim_design)
stim_fit2 <- eBayes(stim_fit)
stim_stats <- topTable(stim_fit2, coef = "stim_targets$StimulusIGF", number = Inf)

treat_design <- model.matrix(~treat_targets$`Cell Line` + treat_targets$Treatment)
treat_fit <- lmFit(treat_rppa, treat_design)
treat_fit2 <- eBayes(treat_fit)

# grab BMS results, then OSI stats. Usually we do this with a loop for contrasts, but that
# is overkill for just two contrasts
BMS_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentBMS", number = Inf)
OSI_stats <- topTable(treat_fit2, coef = "treat_targets$TreatmentOSI", number = Inf)

# combine the _stats results together for plotting
allStatsData <- bind_rows(stim_stats %>% rownames_to_column("Protein"),
                          BMS_stats %>% rownames_to_column("Protein"),
                          OSI_stats %>% rownames_to_column("Protein"),
                          .id = "Treatment") %>%
  left_join(data.frame("Treatment" = c("1", "2", "3"),
                       "Condition" = c("StimVsUnstim",
                                       "StimBMS_Vs_Stim",
                                       "StimOSI_Vs_Stim")))

allStatsData <- bind_rows(stim_stats %>% rownames_to_column("Protein"),
                          BMS_stats %>% rownames_to_column("Protein"),
                          OSI_stats %>% rownames_to_column("Protein"),
                          .id = "Treatment") %>%
  left_join(data.frame("Treatment" = c("1", "2", "3"),
                       "Condition" = c("Uninhibited",
                                       "BMS",
                                       "OSI")))
allStatsData$Condition <- relevel(allStatsData$Condition, "Uninhibited")
# filter to significant genes, then pass into ggplot for heatmap
allStatsData %>% filter(adj.P.Val <=0.05) %>% write_csv("FDR 5 Percent.csv")
# allStatsData <- allStatsData %>% filter(adj.P.Val <=0.05)
allStatsData <- merge(allStatsData, crossref, by = "Protein", all.x = TRUE)

allStatsData %>%
  filter(adj.P.Val <=0.05) %>% # Change this cutoff if you want a different FDR (adj.P.Val) 
  # now passed into ggplot, notice how the ggplot "+" operator is analogous to "%>%" operator for tidy
  ggplot(aes(Condition, fct_rev(Protein))) +
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(name = expression("log"[2] * "(Fold Change)"), low = "royalblue4", high = "red", mid = "white") + 
  theme_bw() + 
  scale_y_discrete(name = "") +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", legend.justification = c("center", "bottom")) + 
  geom_text(aes(label = "*"), data = allStatsData %>% filter(adj.P.Val<=0.01), vjust = 0.75, hjust = 0.5) + facet_grid(Pathway ~ ., scales = "free_y", space = "free")
# ggsave saves the last plot made, but you can pass it saved plots (i.e. p <- ggplot(...))
# ggsave(p, fileName, etc...)

# ggsave("2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Heat-Map APD 5-PCT New.pdf",
#        width = 8.5,
#        height = 11)
ggsave("Version 3 bottom.png",
       width = 5.5,
       height = 11, dpi = 300, device = "png")




### 10

# ggsave("~/BTSync/MASS_SPEC/CrazyRPPA/2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Heat-Map APD.pdf",
#        width = 8.5,
#        height = 11)

# save stats to csv file, openable by excel and R
# allStatsData %>% write_csv("~/BTSync/MASS_SPEC/CrazyRPPA/2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Results-Table.csv")

allStatsData %>% write_csv("2017-Dec-01-OSI-Or-BMS-Vs-IGF-Stim-Paired-By-Cell-Line-Results-Table.csv")

allStatsData %>% filter(adj.P.Val <=0.05) %>% write_csv("FDR 5 Percent.csv")

source("FunctionsFromBrian.R")
myClustResult <- pvclust(allStatsData, method.hclust="average", method.dist="correlation", parallel=TRUE)

plot(myClustResult) # dendogram with p values
#add rectangles around groups highly supported by the data
pvrect(myClustResult, alpha=.1)