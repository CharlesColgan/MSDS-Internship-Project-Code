function(data, SampleInfo, index, IDs, design=NULL, p=0.05){
  ##Processing to ExpressionSet
  
  ##Filter out NAs by removal
  
  #Create new instance of base data
  proteins=data
  
  #remove proteins with at least 1 missing value
  proteins=proteins[complete.cases(proteins[,index]),]
  
  #Processing to ExpressionSet
  exprs=as.matrix(proteins[,index])
  
  exprs=normalizeMedianValues(exprs)
  
  rownames(exprs)=as.vector(as.matrix(proteins[,grep(IDs,colnames(proteins))]))
  
  feats=as.data.frame(proteins[,-index])
  
  rownames(feats)=as.vector(as.matrix(feats[,grep(IDs,colnames(proteins))]))
  
  feats[,grep(IDs,colnames(proteins))]=NULL
  
  feats=as(feats,"AnnotatedDataFrame")
  
  Proteins=ExpressionSet(exprs)
  
  ##Design for Comparison between Experiments
  
  if(is.null(design)==TRUE){
    SampleInfo=as.data.frame(unclass(SampleInfo))
    
    for(i in 1:dim(SampleInfo)[2]){
      SampleInfo[,i]=as.factor(SampleInfo[,i])
    }
    
    #Model Matrix
    design=model.matrix(~.,data=SampleInfo)
  }
  
  ##Fit and Bayes for filtered data
  fit=lmFit(Proteins,design,method="robust")
  
  Bayes_fit=eBayes(fit,
                   proportion=0.01,
                   stdev.coef.lim=c(0.1,4),
                   trend=TRUE,
                   robust=TRUE)
  
  resultsB=decideTests(Bayes_fit)
  
  filt=summary(resultsB)
  
  Tab1=topTable(Bayes_fit,number=dim(exprs(Proteins))[1])
  
  write.csv(Tab1,"Out/Fit.csv")
  
  #Summaries of the regulation of the proteins for each factor, given in percent
  filt_sum=round(filt/dim(exprs)[1],3)*100
  
  ##Filtering data for only proteins with adjusted Pvalue<0.05 using Benjamin-Hochberg method
  Sig=Tab1[Tab1$adj.P.Val<=p,]
  
  #Q-Q plot with fisher combined pvalue for data set with removed NAs
  #return(plotp(Tab1$adj.P.Val))
  #return(title(paste("Combined P=",signif(sumlog(Tab1$adj.P.Val)$p,4))))
  
  #List for return
  list=list("Summary"=filt,
            "regulation"=filt_sum,
            "Significants"=Sig)
  
  return(list)
}
