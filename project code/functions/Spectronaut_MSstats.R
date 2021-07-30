function(data, comparison){
  #Cleaning and formatiing from Spectronaut format to Msstats format
  MSSU_Clean=SpectronauttoMSstatsFormat(data,log_file_path="MSstat_Logs/MSstats_format_log.log")
  
  #=====================
  #Function: dataProcess
  #pre-processing data quality control of MS format to be used later
  #Defaults to log transform of base 2 and equalized median normalization,
  #
  quantData=dataProcess(MSSU_Clean,log_file_path="MSStat_Logs/MSstats_dataProcess_log.log")
  
  # finaly use the groupcomparson function
  resultMultiComparisons<-groupComparison(contrast.matrix=comparison,
                                          data=quantData,
                                          log_file_path="MSStat_Logs/MSstats_groupComparison_log.log")
  Comp_Sult=resultMultiComparisons$ComparisonResult
  
  ## turn groupComparision into a data frame
  GCF=data.frame(Comp_Sult)
  #Filter out NA values for log2FC
  GCF=GCF[complete.cases(GCF$log2FC),]
  
  ###############################################
  #Explore possible cutoffs values for significant proteins
  proteins=1:nrow(GCF)
  log2FC=sort(GCF$log2FC)
  #Qvalues=sort(GCF$adj.pvalue)
  #expl=data.frame(proteins,log2FC,Qvalues)
  #ggplot(expl,aes(x=proteins,y=log2FC))+geom_point()+geom_hline(yintercept=1,color="red")+geom_hline(yintercept=-1,color="red")
  #ggplot(expl,aes(x=proteins,y=Qvalues))+geom_point()+geom_hline(yintercept=0.05,color="red")
  
  #Filtering based on chosen cutoffs
  GCF=GCF[abs(GCF$log2FC)>1&GCF$adj.pvalue<0.05,]
  rownames(GCF)=NULL
  write.csv(GCF,"Out/Spectronaut_Significant.csv")

  #####################################
  
  #=====================
  # Function: dataProcessPlots
  # visualization 
  # takes about 15mins to load for each 
  
  #Profile Plot
  dataProcessPlots(quantData,
                   type="ProfilePlot",
                   which.Protein=as.vector(GCF$Protein),
                   ylimUp=35,
                   address="Out/Spec_")
  
  #QC Plot
  dataProcessPlots(quantData,
                   type="QCPlot",
                   which.Protein=as.vector(GCF$Protein),
                   ylimUp=35,
                   address="Out/Spec_")  
  
  #Condition Plot
  dataProcessPlots(quantData,
                   type="ConditionPlot",
                   which.Protein=as.vector(GCF$Protein),
                   address="Out/Spec_")
  
  # extra functions/ other functions the package provides 
  # Function: groupComparisonPlots
  # visualization for testing results
  
  # Visualization 1: Volcano plot
  # Both FDR cutoff = 0.05 and fold change cutoff = 1.5
  groupComparisonPlots(data=GCF,
                       type="VolcanoPlot",
                       ylimUp=70,
                       FCcutoff=1.5,
                       ProteinName=FALSE,
                       which.Protein=as.vector(GCF$Protein),
                       address="Out/Spec_")
  
  # Visualization 2: Heatmap (required more than one comparisons)2
  # Both FDR cutoff = 0.05 and fold change cutoff = 1.5
  groupComparisonPlots(data=GCF,
                       type="Heatmap",
                       FCcutoff=1.5,
                       which.Protein=as.vector(GCF$Protein),
                       address="Out/Spec_")
  
  # Visualization 3: Comparison plot
  groupComparisonPlots(data=GCF,
                       type="ComparisonPlot",
                       which.Protein=as.vector(GCF$Protein),
                       address="Out/Spec_")
  
  return(GCF)
}
