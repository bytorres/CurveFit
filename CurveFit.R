## Better Curve Fit Analysis Part 1 CVR FIXED Changing it up

#Variables
residuals_nk_loessas=matrix(nrow=nrow(y),ncol=ncol(y))
fit_nk_loessas=matrix(nrow=nrow(y),ncol=ncol(y))
span_nk_loessas=matrix(nrow=nrow(y),ncol=1)

# Curve Fit Every Gene
for(n in 1:nrow(y))
{
  curvefit_nk_loessas=loess.as(as.numeric(x_nk),as.numeric(y[n,]),plot=TRUE)
  #Test by checking if the fit works
  fit_nk_loessas[n,]=as.numeric(curvefit_nk_loessas$fitted)
  #plot(as.numeric(x),fit_nk[34030,])
  residuals_nk_loessas[n,]=as.numeric(curvefit_nk_loessas$residuals)
  span_nk_loessas[n]=as.numeric(curvefit_nk_loessas$pars$span)
}

write.table(residuals_nk_loessas,'nk_resid_bfit.txt',sep='\t')
nk_resid_bfit_file_name=paste('nk_resid_bfit.txt')

# Find the Top %10 percent variation

nk_bfit_range=matrix(nrow=nrow(fit_nk_loessas),ncol=1)

#Find the Range of the Curve Fit
for(n in 1:nrow(fit_nk_loessas))
{
  nk_max_bfit<-max(fit_nk_loessas[n,])
  nk_min_bfit<-min(fit_nk_loessas[n,])
  nk_bfit_range[n]<-nk_max_bfit-nk_min_bfit
}


nk_bfit_top<-nk_bfit_range >.35
nk_bfit_top_index<-which(nk_bfit_range >.35)

#
nk_bfit_symbol<-Symbol[as.numeric(nk_bfit_top_index)]
nk_bfit_top_exp<-Exp1[nk_bfit_top_index,]
nk_bfit_top_fit<-fit_nk_loessas[nk_bfit_top_index,]
nk_bfit_top_resid<-residuals_nk_loessas[nk_bfit_top_index,]


write.table(nk_bfit_top_resid,'nk_bfit_top_resid.txt',sep="\t")
write.table(nk_bfit_top_exp,'nk_bfit_top_exp.txt',sep="\t")
nk_resid_bfit_file_name=paste('nk_bfit_top_resid.txt')

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste("cvrfix.txt");
covariates_file_name = character();

cvrtfix1 = SlicedData$new();
cvrtfix1$fileDelimiter = "\t"; # the TAB character
cvrtfix1$fileOmitCharacters = "NA"; # denote missing values;
cvrtfix1$fileSkipRows = 1; # one row of column labels
cvrtfix1$fileSkipColumns = 1; # one column of row labels
if(length(covariates_file_name)>0) {
  cvrtfix$LoadFile(covariates_file_name);
}

nkbfit = SlicedData$new();
#residual$fileDelimiter = "\t";# the TAB character
nkbfit$fileOmitCharacters = "NA"; # denote missing values;
nkbfit$fileSkipRows = 1;          # one row of column labels
nkbfit$fileSkipColumns = 1;       # one column of row labels
nkbfit$fileSliceSize = 2000;      # read file in slices of 2,000 rows
nkbfit$LoadFile(nk_resid_bfit_file_name);

expression_file_name = paste("nk_bfit_top_exp.txt");

gene_nk = SlicedData$new();
#gene$fileDelimiter = " ";      # the TAB character
gene_nk$fileOmitCharacters = "NA"; # denote missing values;
gene_nk$fileSkipRows = 1;          # one row of column labels
gene_nk$fileSkipColumns = 1;       # one column of row labels
gene_nk$fileSliceSize = 2000;      # read file in slices of 2,000 rowsB0002  B0003  B0005  B0006  B0007	B0008	B0009	B0010	B0011	B0012	B0013	B0014	B0015	B0016	B0017	B0018	B0019	B0020	B0021	B0022	B0023	B0024	B0026	B0025	B0027	B0028	B0029	B0030	B0031	B0032	B0033	B0034	B0035	B0036	B0038	B0040	B0041	B0043	B0044	B0045	B0046	B0047	B0048	B0049	B0051	B0052	B0053	B0054	B0055	C0001	C0002	C0003	C0004	C0005	C0006	C0007	C0008	C0009	C0010	C0011	C0012	C0013	C0014	C0016	C0017	C0018	C0019	C0020	C0021	C0022	C0023	C0024	C0025	C0026	C0027	C0028	C0029	C0030	C0031	C0032	C0033	C0034	C0035	C0036	C0037	C0038	C0039	HFS101	HFS105	HFS106	HFS109	HFS110	HFS111	HFS112	HFS113	HFS114	HFS115	HFS116	HFS117	HFS118	HFS119	HFS120	HFS122	HFS124	HFS126	HFS127	HFS129	HFS130	HFS132	HFS134	HFS135	HFS136	HFS137	HFS138	HFS139	HFS140	HFS141	HFS143	HFS144	HFS145	HFS146	HFS147	HFS148	HFS149	HFS150	HFS151	HFS152	HFS153	HFS154	HFS157	HFS158	HFS159	HFS160	HFS161	HFS162	HFS163	HFS164	HFS165	HFS166	HFS167	HFS168	HFS169	HFS170	HFS171	HFS172	HFS173
gene_nk$LoadFile(expression_file_name);

meq_nkbfit_gene_fix= Matrix_eQTL_engine(
  snps = snps2, # contains all samples (+ ctrl) 
  gene = gene_nk, # Contains all samples (+ctrl)
  cvrt = cvrtfix, 
  output_file_name = tempfile(),
  pvOutputThreshold = 1e-10, 
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = "qqplot");



View(meq_nkbfit$all$eqtls)
nk_bfit_hits<-meq_nkbfit$all$eqtls
write.table(nk_bfit_hits,'nk_bfit_hits1.txt',sep="\t")
write.table(nk_bfit_hits_fix,'nk_bfit_hits_fix.txt',sep="\t")
plot(meq_nkbfit, pch = 16, cex = 0.7)




nk_bfit_snps_fix<-as.character(meq_nkbfit_fix$all$eqtls$snps)
View(nk_bfit_snps_fix)
nk_bfit_resid_snps_fix<-as.numeric(nk_bfit_snps_fix)
View(nk_bfit_resid_snps_fix)


# SNP:Check to make sure you get the same snp hits
data2_add@gtdata@snpnames[193590]#[20836]
nk_bfit_gene_snp_fix<-as.numeric(data2_add[,nk_bfit_resid_snps_fix])
View(nk_bfit_resid_genotype_num_fix)
write.table(nk_bfit_resid_genotype_num_fix,'nk_bfit_resid_genotype_fix.txt',sep='\t')


# GENE SYMBOL: had to download the nk10_expression list []1
nk_bfit_resid_gene_fix<-meq_nkbfit_fix$all$eqtls$gene
nk_bfit_resid_gene_fix<-as.character(nk_bfit_resid_gene_fix)
nk_bfit_resid_gene_fix<-as.numeric(nk_bfit_resid_gene_fix)
nk_bfit_resid_symbol_fix<-nk_bfit_symbol[nk_bfit_resid_gene_fix]
write.table(nk_bfit_resid_symbol_fix,'nk_bfit_symbol_hit_fix.txt',sep="\t")



#GENE EXPRESSION
nk_bfit_resid_exp_fix<-nk_bfit_top_exp[nk_bfit_resid_gene_fix,]
write.table(nk_bfit_resid_exp_fix,'nk_bfit_resid_exp_fix.txt',sep="\t")
#GENE FIT
nk_bfit_fit_hits_fix<-nk_bfit_top_fit[nk_bfit_resid_gene_fix,]
write.table(nk_bfit_fit_hits_fix,'nk_bfit_fit_hits_fix.txt',sep="\t")



quartz()
par(mfrow=c(5,10))
for(n in 1:nrow(nk_bfit_top_exp))
{
  plot(as.numeric(x_nk),nk_bfit_top_exp[n,],pch=16)
  fitank1<-nk_bfit_top_exp[n,]
  t<-order(x_nk)
  lines(x_nk[t],fitank1[t],col='red',lwd=3)
  
}


