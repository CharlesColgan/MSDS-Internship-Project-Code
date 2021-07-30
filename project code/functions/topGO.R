function(Table, Index){
  #Mark significants
  Table$sig=0
  for(i in 1:dim(Table)[1]){
    if(Table$adj.P.Val[i]<0.05){
      Table$sig[i]=1
    }
  }
  
  #Form gene list
  scores=Table$sig
  names(scores)=Table$X1
  scores=as.factor(scores)
  
  #Create mapping
  annot=select(Index,Table$X1,columns="GO-ID",keytypes="UNIPROTKB")
  gene2GO=list()
  for(i in 1:dim(annot)[1]){
    a=str_split(annot[i,2],"; ")
    names(a)=annot$UNIPROTKB[i]
    gene2GO=append(gene2GO,a)
  }
  
  #initialise topGO object
  GOdat=new("topGOdata",
            description="Simple session",
            ontology="BP",
            allGenes=scores,
            nodeSize=10,
            annot=annFUN.gene2GO,
            gene2GO=gene2GO)
  
  #Perform tests
  resultFisher=runTest(GOdat,
                       algorithm="classic",
                       statistic="fisher")
  resultFisher.elim=runTest(GOdat,
                            algorithm="elim",
                            statistic="fisher")
  resultKS=runTest(GOdat,
                   algorithm="classic",
                   statistic="ks")
  resultKS.elim=runTest(GOdat,
                        algorithm="elim",
                        statistic="ks")

  #Tabulate and compare tests
  allRes=GenTable(GOdat,
                  classicFisher=resultFisher,
                  elimFisher=resultFisher.elim,
                  classicKS=resultKS,
                  elimKS=resultKS.elim,
                  orderBy="elimFisher",
                  ranksOf="classicFisher",
                  topNodes=100)
  
  #P-values for KS on classic and elim
  pValue.classic.KS=score(resultKS)
  pValue.elim.KS=score(resultKS.elim)[names(pValue.classic.KS)]
  
  #Plot elements
  gstat=termStat(GOdat, names(pValue.classic.KS))
  gSize=gstat$Annotated/max(gstat$Annotated)*4
  colMap=function(x){
    .col=rep(rev(heat.colors(length(unique(x)))),time=table(x))
    return(.col[match(1:length(x), order(x))])
  }
  gCol=colMap(gstat$Significant)
  
  #Scatter Plot
  SCP=plot(pValue.classic.KS,
       pValue.elim.KS,
       xlab="p-value classic",
       ylab="p-value elim",
       pch=19,
       cex=gSize,
       col=gCol)
  
  #GO subgraph
  #based on p-values and the n most significant terms
  SSN=showSigOfNodes(GOdat,
                     termsP.value=score(resultKS.elim),
                     firstSigNodes=5,
                     useInfo='all')
  
  listGO=list("topGOobject"=GOdat,
            "SCP"=allRes,
            "Subgraph"=SSN,
            "Genes"=names(gene2GO),
            "GO Terms"=names(inverseList(gene2GO)))
  
  write.csv(allRes,"Out/topGo_results.csv")
  
  return(listGO)
}
