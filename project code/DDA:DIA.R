#########################################
#Installation and initialization

#Should install Xcode(if Mac) or Rtools(if Windows) to help things run smoother

#if asked to update packages or to install packages from source, reply yes to all

#HELPFUL: If this is the first time running, run this section by itself, then proceed to the rest of the workflow, this will make things run smoother

Dependencies=dget("functions/Dependencies.R")
Dependencies()

#Spectroanut data for Msstats
MSSU=read_csv("data/peptide_level_Ecoli83972wt-NHS-CC_labeling_07122021_Report.csv")
comparison1=matrix(c(1,-1,0,0,0,0,0,0,0,0,0,0),nrow=1)
comparison2=matrix(c(1,0,-1,0,0,0,0,0,0,0,0,0),nrow=1)
comparison3=matrix(c(1,0,0,-1,0,0,0,0,0,0,0,0),nrow=1)
comparison4=matrix(c(0,0,0,0,1,-1,0,0,0,0,0,0),nrow=1)
comparison5=matrix(c(0,0,0,0,1,0,-1,0,0,0,0,0),nrow=1)
comparison6=matrix(c(0,0,0,0,1,0,0,-1,0,0,0,0),nrow=1)
comparison7=matrix(c(0,0,0,0,0,0,0,0,1,-1,0,0),nrow=1)
comparison8=matrix(c(0,0,0,0,0,0,0,0,1,0,-1,0),nrow=1)
comparison9=matrix(c(0,0,0,0,0,0,0,0,1,0,0,-1),nrow=1)
comparison=rbind(comparison1,
                 comparison2,
                 comparison3,
                 comparison4,
                 comparison5,
                 comparison6,
                 comparison7,
                 comparison8,
                 comparison9)
colnames(comparison)=c(unique(MSSU$R.Condition)[5],
                       unique(MSSU$R.Condition)[11],
                       unique(MSSU$R.Condition)[6],
                       unique(MSSU$R.Condition)[12],
                       unique(MSSU$R.Condition)[3],
                       unique(MSSU$R.Condition)[9],
                       unique(MSSU$R.Condition)[4],
                       unique(MSSU$R.Condition)[10],
                       unique(MSSU$R.Condition)[1],
                       unique(MSSU$R.Condition)[7],
                       unique(MSSU$R.Condition)[2],
                       unique(MSSU$R.Condition)[8])
rownames(comparison)=c(unique(MSSU$R.Condition)[11],
                       unique(MSSU$R.Condition)[6],
                       unique(MSSU$R.Condition)[12],
                       unique(MSSU$R.Condition)[9],
                       unique(MSSU$R.Condition)[4],
                       unique(MSSU$R.Condition)[10],
                       unique(MSSU$R.Condition)[7],
                       unique(MSSU$R.Condition)[12],
                       unique(MSSU$R.Condition)[8])

Spectronaut_MSstats=dget("functions/Spectronaut_MSstats.R")
MSstat=Spectronaut_MSstats(MSSU,comparison)

#Limma Data

#Read-in
peptide=read_csv("data/peptidePivotTable_Ecoli83972wt-NHS-CC_labeling_07122021_Report.csv")
SampleInfo=read_csv("data/RunSummary_Ecoli83972wt-NHS-CC_labeling_07122021_Report.csv")

#Formatting
#Formatting is user defined and will differ between input file formats
Names=paste(peptide$PG.ProteinAccessions,
            peptide$PG.Genes,
            peptide$PEP.StrippedSequence,
            peptide$EG.ModifiedSequence,
            sep=":")
PEP.Quantity=peptide[,grep("PEP.Quantity",names(peptide))]
names(PEP.Quantity)=SampleInfo$R.Condition
PEP_Quan=cbind(Names,PEP.Quantity)

PEP_Quan$Names[duplicated(PEP_Quan$Names)]=paste(PEP_Quan$Names[duplicated(PEP_Quan$Names)],"A",sep=":")
PEP_Quan$Names[duplicated(PEP_Quan$Names)]=paste(PEP_Quan$Names[duplicated(PEP_Quan$Names)],"B",sep=":")

#Expreimetnal design matrix
Time=SampleInfo$Time
Treat=as.factor(SampleInfo$Treat)
Enrichment=as.factor(SampleInfo$Enrichment)
design=model.matrix(~Time+Treat+Enrichment)

#Instantiate Limma function
Limma=dget("functions/Limma.R")

#Limma(data, SampleInfo, index, IDs, design)
#data: formatted data set
#SampleInfo: small table describing which conditions are present in which samples
#index: indeces of columns in data that contains the value being quantified
#IDs: column name from data which contains identifiers for the rows of data
#design: the experimental design matrix for the comparison/fit being tested, 
#defualts to using every column of SampleInfo individually
#p: cutoff for adjusted p-values to be considered significant, defualts to p<=0.05
LM=Limma(data=PEP_Quan,
         SampleInfo=SampleInfo[,c(1,3,4,5,6)],
         index=c(2:13),
         IDs="Names",
         design=design,
         p=0.05)

#topGo Data

#Formatting
#this comes from significant results found by Limma
Tab=read_csv("Out/Fit.csv")
Tab$X1=gsub("\\:.*","",Tab$X1)
Tab$X1=gsub("\\;.*","",Tab$X1)


#Gene ontologies Ids

#replace taxId with the relevant species Id,
#use availableUniprotSpecies(pattern,n=Inf) to search Ids, with pattern allowing for filtering of results, and n being the maximum number of resutls to return
#defualt Id is 9606 which is for human species
#will take a while
index=UniProt.ws(taxId=9606)

#Instantiate topGO function
topGO=dget("functions/topGO.R")

#topGO(data, index)
#data: data set
#index: the index from UniProt that contains the Id relations
#column: column of the data set withnthe Id
topGO=topGO(Tab,index)

#GeneWalk
annot=select(index,Tab$X1,columns="ENSEMBL",keytype="UNIPROTKB")
ID_MAP=list()
for(i in 1:dim(annot)[1]){
  a=str_split(annot[i,2],"; ")
  names(a)=annot$UNIPROTKB[i]
  ID_MAP=append(ID_MAP,a)
}

#Exports for GeneWalk
write.table(as.factor(names(inverseList(ID_MAP))),
            "Out/ENS.csv",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
write.csv(data.frame(project="test",
                     genes="Out/ENS.csv",
                     id_type="ensembl_id",
                     nproc=4),"Out/ARGS.csv")

#Run GeneWalk
py_run_file("functions/Gwalk.py",convert=FALSE)

#Removal of excess functions and objects for cleaner environment
rm(Dependencies,a,index,annot,ID_MAP,i,Limma,topGO)

#Garbage collection for better memory usage
gc()

#To see results for any function, call the object name in the console, 
#e.g. to look at the limma results, call LM