
#######################################################
################ Expression correlation ###############
#######################################################

dirName= "/Data/TCGA/cbioportal" #Put the location of cbioportal data of all tissue
dirNames= list.dirs(dirName,recursive=FALSE)
plotList= list()
result= data.frame()
for( i in 1:length(dirNames)){
  dataTCGA= read.table(paste(dirNames[i],"/tcga/data_RNA_Seq_v2_expression_median.txt",sep=""),stringsAsFactors = F,header = T)
  dataMYC_UBR5= data.frame(SampleID= colnames(dataTCGA[,-1:-2]), MYC_exp= as.numeric(dataTCGA[dataTCGA$Hugo_Symbol == "MYC",-1:-2]), UBR5_exp= as.numeric(dataTCGA[dataTCGA$Hugo_Symbol == "UBR5",-1:-2]))
  corrVal= cor.test(log2(as.numeric(dataMYC_UBR5$MYC_exp)),log2(as.numeric(dataMYC_UBR5[,3])))
  result= rbind(result,data.frame(Tissue= substr(dirNames[i],25,nchar(dirNames[i])),sampleSize= length(dataMYC_UBR5$SampleID),corr= corrVal$estimate,pVal= corrVal$p.value))
}


#######################################################
############### Co-amplification analysis #############
#######################################################

UBR5startLoc=102252273;
UBR5endLoc= 102412841;

MYCstartLoc=127735434;
MYCendLoc=127741434;
THRES= 0.2


dirName= "Data/TCGA/" #Put the location of TCGA tissue data
tissueNames= c("Breast","Pancreas","Ovary","Lymph")

for (tissueName in 1:length(tissueNames)){
  dirNames= list.dirs(paste(dirName,tissueNames[tissueName],"/",sep=""),recursive=FALSE);
  rm(dataBoth,dataUBR5,dataMYC)
  for (i in 1:length(dirNames)){
    fileNames= list.files(dirNames[i],pattern="*seg.txt")
    for (j in 1:length(fileNames)){
      dataSample= read.table(paste(dirNames[i],"/",fileNames[j],sep=""),header=TRUE,sep="\t",stringsAsFactors = FALSE)
      resultUBR5= subset(dataSample, Chromosome == 8 & Start <=  UBR5startLoc & End >= UBR5endLoc)
      resultMYC= subset(dataSample, Chromosome == 8 & Start <=  MYCstartLoc & End >= MYCendLoc)
      resultboth= subset(dataSample, Chromosome == 8 & Start <=  UBR5startLoc & End >= MYCendLoc)
      if (exists("dataBoth")){
        dataBoth= rbind(dataBoth, resultboth);
      }
      else{
        dataBoth= resultboth
      }
      if (exists("dataUBR5")){
        dataUBR5= rbind(dataUBR5, resultUBR5);
      }
      else{
        dataUBR5= resultUBR5
      }
      if (exists("dataMYC")){
        dataMYC= rbind(dataMYC, resultMYC);
      }
      else{
        dataMYC= resultMYC
      }
    }
  }
  mergedData= merge(dataUBR5,dataMYC, by= "Sample")
  mergedData= mergedData[-c(5,7,10)]
  colnames(mergedData)= c("Sample","Chromosome","Start_UBR5","End_UBR5","Segment_Mean_UBR5","Start_MYC","End_MYC","Segment_Mean_MYC")
  
  write.csv(mergedData,paste("Amplification_",tissueNames[tissueName],"_results.csv",sep=""))
  
  amplData= subset(mergedData,Segment_Mean_UBR5 > THRES | Segment_Mean_MYC > THRES)
  
  message(paste("***********",tissueNames[tissueName],"***********"))
  message(paste("% patients: ",100*(length(amplData$Sample)/(length(mergedData$Sample))),sep=""))
  message(paste("% Both: ",100*(sum(amplData$Segment_Mean_UBR5 > THRES & amplData$Segment_Mean_MYC > THRES))/(length(amplData$Sample)),sep=""))
  message(paste("% UBR5: ",100*(sum(amplData$Segment_Mean_UBR5 > THRES & amplData$Segment_Mean_MYC < THRES))/(length(amplData$Sample)),sep=""))
  message(paste("% MYC: ",100*(sum(amplData$Segment_Mean_UBR5 < THRES & amplData$Segment_Mean_MYC > THRES))/(length(amplData$Sample)),sep=""))
  
}


#######################################################
################ Essentiality analysis ################
#######################################################

zGarpscores <- read.table("Data/breast_zgarp.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t",fill=TRUE)
cnaData<- read.table("Data/breast_cna_lrr.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t",fill=TRUE)

essUBR5= subset(zGarpscores, symbol=="UBR5" | symbol=="MYC", select=c(-1:-3))
essUBR5= rbind(colnames(zGarpscores[,-1:-3]),essUBR5)
rownames(essUBR5)= c("Cell_line","Essentiality_MYC","Essentiality_UBR5")
essUBR5= t(essUBR5)
cna_gene_Breast <- subset(cnaData, symbol== "MYC" | symbol== "UBR5",select=c(-1:-3))
cna_gene_Breast= rbind(colnames(cnaData[,-1:-3]),cna_gene_Breast)
rownames(cna_gene_Breast)= c("Cell_line","Segment_Mean_MYC","Segment_Mean_UBR5")
cna_gene_Breast= t(cna_gene_Breast)
results <-merge(essUBR5,cna_gene_Breast,by="Cell_line")
write.csv(results,"Essentiality_Amplification_celllines_Marcotte.csv")


resultsData <- read.table("Essentiality_Amplification_celllines_Marcotte.csv",header=TRUE,stringsAsFactors = FALSE,sep=",",fill=TRUE)
minTest= cor.test(resultsData$Essentiality_UBR5[resultsData$Essentiality_UBR5 < -1.5],apply(rbind(resultsData$Segment_Mean_MYC[resultsData$Essentiality_UBR5 < -1.5],resultsData$Segment_Mean_UBR5[resultsData$Essentiality_UBR5 < -1.5]),2,min))
qplot(resultsData$Essentiality_UBR5[resultsData$Essentiality_UBR5 < -1.5],apply(rbind(resultsData$Segment_Mean_MYC[resultsData$Essentiality_UBR5 < -1.5],resultsData$Segment_Mean_UBR5[resultsData$Essentiality_UBR5 < -1.5]),2,min),xlab="Essentiality of UBR5",ylab="Min (MYC Amplifi.,UBR5 Amplfi.)",main=paste("correlation=",round(as.numeric(minTest$estimate),2),"(p-value=",round(as.numeric(minTest$p.value),2),")")) + geom_smooth(method=lm,se=FALSE)

