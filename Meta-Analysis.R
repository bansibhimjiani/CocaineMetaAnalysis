#Pipeline Aligning Results 
#Bansi Bhimjiani
#July 28 2026

############

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

############
library(dplyr)
library(plyr)

if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}

0

############

#Download the necessary functions from our Github repository into your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_AligningDEResults.R

############

#Reading in the functions:

source("Function_AligningDEResults (1).R")

###########

#Aligning the mouse datasets with each other:

#Example Usage:

ListOfMouseDEResults <-list (DEResults_GSE43410_529724, DEResults_GSE18751_516219, DEResults_GSE63748_536281, DEResults_GSE9134_560143, DEResults_GSE76868_550679, DEResults_GSE276453_590389, DEResults_GSE141520_642255, DEResults_GSE110344_492905, DEResults_GSE110344_492906)

str(ListOfMouseDEResults[[3]])

AligningMouseDatasets(ListOfMouseDEResults)

str(Mouse_MetaAnalysis_FoldChanges)

ListOfRatDEResults <-list (DEResults_GSE169665_579445, DEResults_GSE66349_524663, DEResults_GSE169666_578302, DEResults_GSE88736_512558)

AligningRatDatasets(ListOfRatDEResults)

str(ListOfRatDEResults[[3]])

str(Rat_MetaAnalysis_FoldChanges)


#################

#Code for aligning the rat and mice results:

#This code isn't nicely functionalized yet
#It also assumes that there are mouse datasets
#It will break if there are only rat datasets - this needs to be fixed.

################

#First: What are gene orthologs?

# Homology refers to biological features including genes and their products that are descended from a feature present in a common ancestor.

# Homologous genes become separated in evolution in two different ways: separation of two populations with the ancestral gene into two species or gene duplication of the ancestral gene within a lineage:

### Genes separated by speciation are called orthologs.
### Genes separated by gene duplication events are called paralogs.

#This definition came from NCBI (https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/ortho_para.html)

#We have the ortholog database that we downloaded from Jackson Labs on April 25, 2024
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20240425.R"

MouseVsRat_NCBI_Entrez<-read.csv("HOM_MouseVsRat_EntrezEnsemblAgree_NoMultimapped_20260511.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

str(MouseVsRat_NCBI_Entrez)

#We want to join this ortholog database to our mouse results (Log2FC and SV):

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	28920 obs. of  30 variables:

Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	28920 obs. of  30 variables:


#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	36792 obs. of  32 variables:

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)
#'data.frame':	36792 obs. of  32 variables:



###############################

#Not all of these have annotation...
#Those would be the genes that have Entrez IDs that aren't 1:1 with Ensembl
#This matters more this year because we will be joining across datasets with different annotation

sum(is.na(MetaAnalysis_FoldChanges$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_FoldChanges$Rat_Ensembl))
#[1] 18577

MetaAnalysis_FoldChanges<-MetaAnalysis_FoldChanges[(is.na(MetaAnalysis_FoldChanges$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_FoldChanges$Rat_Ensembl))==FALSE,]

dim(MetaAnalysis_FoldChanges)
#[1] 22695    48

MetaAnalysis_SV<-MetaAnalysis_SV[(is.na(MetaAnalysis_SV$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_SV$Rat_Ensembl))==FALSE,]

dim(MetaAnalysis_SV)
#[1] 22695    48

#For simplicity's sake for labeling later charts, let's make a combo annotation with both mouse and rat gene symbol:
MetaAnalysis_FoldChanges$MouseRat_GeneSymbol<-paste(MetaAnalysis_FoldChanges$Mouse_Symbol, MetaAnalysis_FoldChanges$Rat_Symbol, sep="_")

MetaAnalysis_SV$MouseVsRat_GeneSymbol<-paste(MetaAnalysis_SV$Mouse_Symbol, MetaAnalysis_SV$Rat_Symbol, sep="_")


#######################

#We should probably spend some time renaming the comparisons included in our MetaAnalysis_FoldChanges and MetaAnalysis_SV objects now...

colnames(MetaAnalysis_FoldChanges)
#  [1] "DB.Class.Key"                                                                                                                                                                                                                            
#  [2] "Mouse_Common.Organism.Name"                                                                                                                                                                                                              
#  [3] "Mouse_NCBI.Taxon.ID"                                                                                                                                                                                                                     
#  [4] "Mouse_Symbol"                                                                                                                                                                                                                            
#  [5] "Mouse_EntrezGene.ID"                                                                                                                                                                                                                     
#  [6] "Mouse_Mouse.MGI.ID"                                                                                                                                                                                                                      
#  [7] "Mouse_HGNC.ID"                                                                                                                                                                                                                           
#  [8] "Mouse_OMIM.Gene.ID"                                                                                                                                                                                                                      
#  [9] "Mouse_Genetic.Location"                                                                                                                                                                                                                  
# [10] "Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38."                                                                                                                                                                                  
# [11] "Mouse_Name"                                                                                                                                                                                                                              
# [12] "Mouse_Synonyms"                                                                                                                                                                                                                          
# [13] "Mouse_ENSEMBLGene.ID"                                                                                                                                                                                                                    
# [14] "Ensembl_Entrez"                                                                                                                                                                                                                          
# [15] "Rat_Common.Organism.Name"                                                                                                                                                                                                                
# [16] "Rat_NCBI.Taxon.ID"                                                                                                                                                                                                                       
# [17] "Rat_Symbol"                                                                                                                                                                                                                              
# [18] "Rat_EntrezGene.ID"                                                                                                                                                                                                                       
# [19] "Rat_Mouse.MGI.ID"                                                                                                                                                                                                                        
# [20] "Rat_HGNC.ID"                                                                                                                                                                                                                             
# [21] "Rat_OMIM.Gene.ID"                                                                                                                                                                                                                        
# [22] "Rat_Genetic.Location"                                                                                                                                                                                                                    
# [23] "Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38."                                                                                                                                                                                    
# [24] "Rat_Name"                                                                                                                                                                                                                                
# [25] "Rat_Synonyms"                                                                                                                                                                                                                            
# [26] "Rat_Ensembl"                                                                                                                                                                                                                             
# [27] "Rat_DBSymbol"                                                                                                                                                                                                                            
# [28] "Rat_DBName"                                                                                                                                                                                                                              
# [29] "GSE43410_cocaine"                                                                                                                                                                                                                        
# [30] "GSE18751_cocaine.delivered.at.dose.1.x.20.mg.kg"                                                                                                                                                                                         
# [31] "GSE18751_cocaine.delivered.at.dose.8.x.20.mg.kg"                                                                                                                                                                                         
# [32] "GSE18751_cocaine.delivered.at.dose.7.x.20.mg.kg...1.x.20.mg.kg.challenge.after.1.week"                                                                                                                                                   
# [33] "GSE63748_cocaine"                                                                                                                                                                                                                        
# [34] "GSE9134_cocaine"                                                                                                                                                                                                                         
# [35] "GSE76868_cocaine"                                                                                                                                                                                                                        
# [36] "GSE276453_cocaine"                                                                                                                                                                                                                       
# [37] "GSE141520_cocaine.delivered.at.dose.10.mg.kg"                                                                                                                                                                                            
# [38] "GSE110344_cocaine.dependence"                                                                                                                                                                                                            
# [39] "GSE110344_cocaine"                                                                                                                                                                                                                       
# [40] "GSE169665_impaired.conditioned.place.preference.behavior"                                                                                                                                                                                
# [41] "GSE169665_enhanced.conditioned.place.preference.behavior"                                                                                                                                                                                
# [42] "GSE66349_cocaine.with.withdrawal.and.delivered.for.duration.1.d..cocaine.with.withdrawal.and.delivered.for.duration.1.d"                                                                                                                 
# [43] "GSE66349_cocaine.with.withdrawal.and.delivered.for.duration.30.d..cocaine.with.withdrawal.and.delivered.for.duration.30.d"                                                                                                               
# [44] "GSE66349_cocaine.with.withdrawal..and.delivered.for.duration.1.d..extinction.test..cocaine.with.withdrawal..and.delivered.for.duration.1.d..extinction.test..cocaine.with.withdrawal..and.delivered.for.duration.1.d..extinction.test"   
# [45] "GSE66349_cocaine.with.withdrawal..and.delivered.for.duration.30.d..extinction.test..cocaine.with.withdrawal..and.delivered.for.duration.30.d..extinction.test..cocaine.with.withdrawal..and.delivered.for.duration.30.d..extinction.test"
# [46] "GSE169666_cocaine.dependence"                                                                                                                                                                                                            
# [47] "GSE169666_resistant.to.towards.cocaine.dependence"                                                                                                                                                                                       
# [48] "GSE88736_cocaine"                                                                                                                                                                                                                        
# [49] "MouseRat_GeneSymbol"  


write.csv(colnames(MetaAnalysis_FoldChanges), "Colnames_MetaAnalysisFoldChanges.csv")




###################

#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:

#Example scatter plot comparing two datasets:


#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(
  as.matrix(MetaAnalysis_FoldChanges[,29:48]),
  use="pairwise.complete.obs",
  method="spearman"
)
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(
  cor(
    as.matrix(MetaAnalysis_FoldChanges[,29:48]),
    use="pairwise.complete.obs",
    method="spearman"
  )
)


str(MetaAnalysis_FoldChanges[,29:48])












################
#Example code showing the pipeline for running a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include heterogeneity statistics, publication bias statistics, and robustness statistics
#Updated version: March 10, 2026
#Updated again: July 28, 2026 (to make the code more generalizable for the full 2026 cohort)

######################

#Grabbing input data and setting the working directory:

######################

#Installing and loading relevant code packages:

if (!require("metafor", quietly = TRUE)){
  install.packages("metafor")
}

library(metafor)

if (!require("plyr", quietly = TRUE)){
  install.packages("plyr")
}

library(plyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")                                 

library(multtest)

######################

#Download the two necessary functions from our Github repository:

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_RunBasicMetaAnalysis.R

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_FalseDiscoveryCorrection.R

######################


#Read in the functions:

source("Function_RunBasicMetaAnalysis.R")

source("Function_FalseDiscoveryCorrection.R")

source("Function_MakeForestPlots.R")

source("Function_VolcanoPlot.R")

######################

#Example Usage:

#Figure out which columns contain differential expression output:
colnames(MetaAnalysis_FoldChanges)

#Then make a variable that is a vector of with the column numbers containing the DE output:
#Example code for making a variable representing the columns 
Columns_DE<-c(29:48)

#Then make a variable that is the shorthand annotation that we will use to label the genes in our output, e.g. MouseRat_GeneSymbol:
Column_GeneName<-49

NumberOfComparisons=19
CutOffForNAs=9
#I have 6 comparisons
#2 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV, Columns_DE, Column_GeneName)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

#Take a peek at the output:

str(metaOutput)
head(metaOutput)
tail(metaOutput)

write.csv(metaOutput, "metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")

colnames(metaOutput)

write.csv(influence_dfbs, "influence_dfbs.csv")
write.csv(influence_cookd, "influence_cookd.csv" )
write.csv(influence_TF, "influence_TF.csv")

###############

#FDR correction:

FalseDiscoveryCorrection(metaOutput, MetaAnalysis_Annotation)
library(multtest)

###############

#Writing out the meta-analysis input for record-keeping:

write.csv(MetaAnalysis_FoldChanges_ForMeta, "MetaAnalysis_FoldChanges_ForMeta.csv")
write.csv(MetaAnalysis_SV_ForMeta, "MetaAnalysis_SV_ForMeta.csv")

###########
# Set axis limits and plot height
Xaxis_LowerAndUpperBound <- c(-7, 7)
HeightInInches <- 5

# Top 15 genes by p-value
MakeForestPlots(metaOutputFDR_annotated, "Gm5617_C8h11orf71", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Ripply2_Ripply2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Fam246a_LOC120095536", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Gimap9_Gimap9", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Kirrel2_Kirrel2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Tfcp2l1_Tfcp2l1", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Oxld1_Oxld1", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Arpp19_Arpp19", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Inafm1_NA", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Fkbp11_NA", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Bola1_Bola1", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Dnajb3_Dnajb3", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Fosb_Fosb", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Znrd2_Znrd2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Pacrg_Pacrg", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)

# Additional recommended genes
MakeForestPlots(metaOutputFDR_annotated, "Rorb_Rorb", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)
MakeForestPlots(metaOutputFDR_annotated, "Plekhf2_Plekhf2", Columns_DE, Xaxis_LowerAndUpperBound, HeightInInches)


###############
# Make a volcano plot
###############

DE_Results <- metaOutputFDR_OrderbyPval

VariableOfInterest <- "Cocaine"

CoefficientCol <- "Log2FC_estimate"

PvalueCol <- "pval"

FDRCol <- "FDR"

VolcanoPlot(
  DE_Results,
  VariableOfInterest,
  CoefficientCol,
  PvalueCol,
  FDRCol
)

###############
# Heatmap of correlations across studies
###############

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

# Calculate the Spearman correlation matrix across study/contrast Log2FC columns
Log2FC_CorMatrix <- cor(
  as.matrix(MetaAnalysis_FoldChanges[, Columns_DE, drop = FALSE]),
  use = "pairwise.complete.obs",
  method = "spearman"
)

# Save the heatmap
pheatmap(
  Log2FC_CorMatrix,
  filename = "Heatmap_CorMatrix_AllStudies.pdf",
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  breaks = seq(-1, 1, length.out = 101),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  width = 8.5,
  height = 11
)
