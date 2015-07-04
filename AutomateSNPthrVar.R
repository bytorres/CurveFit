#Load packages
#rSNPS
library("rsnps", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
#library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")


# Variables
t<-numeric(0)
nksnp_data<-ctrlsnp_all#sel_snp_gtype#ctrlsnp_all   # is the unique QTL hits
choosen_idx<-ANOVAList#c(1:ncol(sel_snp_gtype)) # idx of significant  ## for the selected analysis there is no specific idx 

symbol<-colnames(ctrlsnp_all)#colnames(ctrlsnp_all)
n=0
i=1

#Run a specific List
##list<-c('rs2843159','rs2842933','rs3094315','rs3855951')
#idx<-which(symbol %in% list)

#choosen_idx<-idx#AllSNPsel
#list<-c('rs10904915','rs7747253','rs11150882','rs7143764')

for(n in 1:length(choosen_idx))  #ncol(sel_snp_gtype)
{
  n
  a<-as.data.frame(table(nksnp_data[,choosen_idx[n]])) # Tally on # of 0,1,2
  #NKG7_nm<-NKG7[which(as.numeric(nksnp_data[,choosen_idx[n]])==2)]
  #RBC_nm<-RBC[which(as.numeric(nksnp_data[,choosen_idx[n]])==2)]
  #variance_2<-sqrt(var(NKG7_nm)+var(RBC_nm))
  
  
  
  if (a[3,2] > 5 && !is.na(a[3,2]) ) ##  #variance_2 < 2 && ##!is.na(variance_2)
  {
    
    #variance_2
    quartz(width=3,height=6.5)
    t[i]<-symbol[choosen_idx[n]]
    t0<-qplot(NKG7,RBC) + geom_point(color=ifelse(nksnp_data[,choosen_idx[n]]=='0',"#86B875",'gray'),size=3)+labs(title=symbol[choosen_idx[n]])
    t1<-qplot(NKG7,RBC) + geom_point(color=ifelse(nksnp_data[,choosen_idx[n]]=='1','#4CB9CC','gray'),size=3) 
    t2<-qplot(NKG7,RBC) + geom_point(color=ifelse(nksnp_data[,choosen_idx[n]]=='2','#CD99D8','gray'),size=3)
    multiplot(t0,t1,t2)
    remove(t0,t1,t2)
    i=i+1
  }
  
  
}


#snp_query_topresid<-NCBI_snp_query(t)
#tableau_heterosmall<-nksnp_data[,choosen_idx]
#View(tableau_heterosmall)
#tableau_heterosmall<-cbind(tableau_heterosmall,NKG7,RBC)
#write.table(tableau_heterosmall,'07_04_15_SigSNPSDynamicRange.txt',sep='\t')

#write.table(snp_query_topresid,'07_03_15_ALLSNPSNCBI.txt',sep='\t')


## Check how many of the hits have more than ten samples equal to 2
#b<-numeric(0)
#c<-numeric(0)


#for(n in 1:length(choosen_idx))
 # {

  #    a<-as.data.frame(table(nksnp_data[,choosen_idx[n]]))
   #   b[n]<-a[3,2]
    #  t[n]<-n

  #}
#NKG7_nm<-NKG7[which(as.numeric(nksnp_data[,choosen_idx[n]])==2)]
#RBC_nm<-RBC[which(as.numeric(nksnp_data[,choosen_idx[n]])==2)]
#variance_2<-sqrt(var(NKG7_nm)+var(RBC_nm))
#c[n]<-variance_2
#remove(RBC_nm)
#remove(NKG7_nm)




#table(b>10)


