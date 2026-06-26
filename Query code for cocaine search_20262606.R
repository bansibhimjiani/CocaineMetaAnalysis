#This is R code for a practice search in GEO
#2026-06-22
#Bansi Bhimjiani 

######################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

#######################



searchGEO

MyQueryTerms <- '(cocain*[All Fields] OR psychostimulant*[All Fields] OR stimulant*[All Fields]) AND ("accumbens"[All Fields] OR hippocamp*[All Fields] OR "dentate gyrus"[All Fields] OR CA1[All Fields] OR CA2[All Fields] OR CA3[All Fields] OR CA4[All Fields] OR "CA field"[All Fields] OR subiculum[All Fields] OR fimbria[All Fields] OR "cornu ammonis"[All Fields]) AND ("Mus musculus"[ORGN] OR "Rattus norvegicus"[ORGN]) AND ("Expression profiling by high throughput sequencing"[DataSet Type] OR "Expression profiling by array"[DataSet Type]) AND "gse"[Filter]'


QueryResults <- searchGEO(MyQueryTerms)

# This will show you an overview of the identified GEO Records:

str(QueryResults)

#Adding columns to hold the additional metadata to our Query Results object:

QueryResults$Citation<-character(length=nrow(QueryResults))
QueryResults$PMID<-character(length=nrow(QueryResults))
QueryResults$Contributor<-character(length=nrow(QueryResults))
QueryResults$Date<-character(length=nrow(QueryResults))
QueryResults$Abstract<-character(length=nrow(QueryResults))

#Looping over each of the identified GEO records and extracting the desired metadata:


for(i in c(1:nrow(QueryResults))){
  
  gse_raw <- getGEO(QueryResults$`Series Accession`[i], GSEMatrix=FALSE)
  
  QueryResults$Citation[i] <- paste(Meta(gse_raw)$citation, collapse=" ")
  
  QueryResults$PMID[i] <- paste(Meta(gse_raw)$pubmed_id, collapse=" ")
  
  QueryResults$Contributor[i] <- paste(Meta(gse_raw)$contributor, collapse = " ")
  
  QueryResults$Date[i] <- paste(Meta(gse_raw)$submission_date, collapse= " ")
  
  QueryResults$Abstract[i] <- Meta(gse_raw)$summary
  
  rm(gse_raw)
}

#Getting an overview of our Query Result object with its new additions:

str(QueryResults)

#Adding empty columns to hold additional information that we will find while reviewing the dataset records:

QueryResults$Tissue<-character(length=nrow(QueryResults))
QueryResults$DevelopmentalStage<-character(length=nrow(QueryResults))
QueryResults$ManipulatedVariables<-character(length=nrow(QueryResults))
QueryResults$Notes<-character(length=nrow(QueryResults))

#Phrases that indicate the variable manipulated:
# Subjects were treated with ___ and rna-sequencing performed after testing behaviorally for...
# Subjects were divided into groups and one group experienced...
# Subjects received one of two interventions...

#Lets add some empty columns for taking inclusion/exclusion notes:

QueryResults$ManipulationUnrelatedToTopic<-character(length=nrow(QueryResults))
QueryResults$WrongTissue<-character(length=nrow(QueryResults))
QueryResults$NotBulkDissection_ParticularCellTypeOrSubRegion<-character(length=nrow(QueryResults))
QueryResults$IncorrectDevelopmentalStage<-character(length=nrow(QueryResults))
QueryResults$NotFullTranscriptome<-character(length=nrow(QueryResults))
QueryResults$MetadataIssues_MissingInfo_Retracted_Duplicated<-character(length=nrow(QueryResults))

QueryResults$Excluded<-character(length=nrow(QueryResults))
QueryResults$WhyExcluded<-character(length=nrow(QueryResults))

#Output the query results as a comma-separated variable file:

write.csv(QueryResults, "QueryResults.csv")
getwd()

