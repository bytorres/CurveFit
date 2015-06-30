## Better Curve Fit Analysis Using Loessas Function

#Variables
x_sel<-c[,2]
y_sel<-Sel_Exp

#Empty Matrixes
residuals_sel=matrix(nrow=nrow(y_sel),ncol=ncol(y_sel))
fit_sel=matrix(nrow=nrow(y_sel),ncol=ncol(y_sel))
span_sel=matrix(nrow=nrow(y_sel),ncol=1)

# Curve Fit Every Gene
for(n in 1:nrow(y_sel))
{
  curvefit_sel=loess.as(as.numeric(x_sel),as.numeric(y_sel[n,]),plot=TRUE)
  #Test by checking if the fit works
  fit_sel[n,]=as.numeric(curvefit_sel$fitted)
  #plot(as.numeric(x),fit_sel[34030,])
  residuals_sel[n,]=as.numeric(curvefit_sel$residuals)
  span_sel[n]=as.numeric(curvefit_sel$pars$span)
}

write.table(residuals_sel,'sel_resid.txt',sep='\t')

## Find the Top %10 percent variation
# Variable
#sel_range=matrix(nrow=nrow(fit_sel),ncol=1)

#Find the Range of the Curve Fit

#for(n in 1:nrow(fit_sel))
#{
 # sel_max_<-max(fit_sel[n,])
#  sel_min_<-min(fit_sel[n,])
 # sel_range[n]<-sel_max__sel_min_
#}


#sel_top<-sel_range >.35
#sel_top_index<-which(sel_range >.35)

#
#sel_symbol<-Symbol[as.numeric(sel_top_index)]
#sel_top_exp<-Exp1[sel_top_index,]
#sel_top_fit<-fit_sel[sel_top_index,]
#sel_top_resid<-residuals_sel[sel_top_index,]


#write.table(sel_top_resid,'sel_top_resid.txt',sep="\t")
#write.table(sel_top_exp,'sel_top_exp.txt',sep="\t")
#sel_resid_file_name=paste('sel_top_resid.txt') # Only top 10% residuals

# Covariates file name

# Set to character() for no covariates
#covariates_file_name = paste("cvrfix.txt");
covariates_file_name = character();

ncvrt = SlicedData$new();
ncvrt$fileDelimiter = "\t"; # the TAB character
ncvrt$fileOmitCharacters = "NA"; # denote missing values;
ncvrt$fileSkipRows = 1; # one row of column labels
ncvrt$fileSkipColumns = 1; # one column of row labels
if(length(covariates_file_name)>0) {
  ncvrt$LoadFile(covariates_file_name);
}



#selsnp_idx<-which(data2_add@gtdata@idnames %in% sel_patients)
#sel_snp_correct<-as.numeric(data2_add@gtdata[selsnp_idx,])
#dorder<-match(sel_patients,row.names(sel_snp_correct))
#Test<-sel_snp_correct[dorder,]
#write.table(Test,'sel_snp_ordered.txt',sep='\t')
sel_resid_file_name=paste('sel_snp_ordered.txt') 

sel = SlicedData$new();
#residual$fileDelimiter = "\t";# the TAB character
sel$fileOmitCharacters = "NA"; # denote missing values;
sel$fileSkipRows = 1;          # one row of column labels
sel$fileSkipColumns = 1;       # one column of row labels
sel$fileSliceSize = 2000;      # read file in slices of 2,000 rows
sel$LoadFile(sel_resid_file_name);




expression_file_name = paste("sel_resid.txt") # ALL residuals

gene_sel = SlicedData$new();
#gene$fileDelimiter = " ";      # the TAB character
gene_sel$fileOmitCharacters = "NA"; # denote missing values;
gene_sel$fileSkipRows = 1;          # one row of column labels
gene_sel$fileSkipColumns = 1;       # one column of row labels
gene_sel$fileSliceSize = 2000;      # read file in slices of 2,000 rowsB0002  B0003  B0005  B0006  B0007	B0008	B0009	B0010	B0011	B0012	B0013	B0014	B0015	B0016	B0017	B0018	B0019	B0020	B0021	B0022	B0023	B0024	B0026	B0025	B0027	B0028	B0029	B0030	B0031	B0032	B0033	B0034	B0035	B0036	B0038	B0040	B0041	B0043	B0044	B0045	B0046	B0047	B0048	B0049	B0051	B0052	B0053	B0054	B0055	C0001	C0002	C0003	C0004	C0005	C0006	C0007	C0008	C0009	C0010	C0011	C0012	C0013	C0014	C0016	C0017	C0018	C0019	C0020	C0021	C0022	C0023	C0024	C0025	C0026	C0027	C0028	C0029	C0030	C0031	C0032	C0033	C0034	C0035	C0036	C0037	C0038	C0039	HFS101	HFS105	HFS106	HFS109	HFS110	HFS111	HFS112	HFS113	HFS114	HFS115	HFS116	HFS117	HFS118	HFS119	HFS120	HFS122	HFS124	HFS126	HFS127	HFS129	HFS130	HFS132	HFS134	HFS135	HFS136	HFS137	HFS138	HFS139	HFS140	HFS141	HFS143	HFS144	HFS145	HFS146	HFS147	HFS148	HFS149	HFS150	HFS151	HFS152	HFS153	HFS154	HFS157	HFS158	HFS159	HFS160	HFS161	HFS162	HFS163	HFS164	HFS165	HFS166	HFS167	HFS168	HFS169	HFS170	HFS171	HFS172	HFS173
gene_sel$LoadFile(expression_file_name);

meq_sel= Matrix_eQTL_engine(
  snps = sel, # contains all samples (+ ctrl) 
  gene = gene_sel, # Contains all samples (+ctrl)
  cvrt = ncvrt, 
  output_file_name = tempfile(),
  pvOutputThreshold = 1e_10, 
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = "qqplot");



#View(meq_sel$all$eqtls)
#sel_hits<-meq_sel$all$eqtls
#write.table(sel_hits,'sel_hits.txt',sep="\t")


#sel_snps_fix<-as.character(meq_sel_fix$all$eqtls$snps)
#View(sel_snps_fix)
#sel_resid_snps_fix<-as.numeric(sel_snps_fix)
#View(sel_resid_snps_fix)


# SNP:Check to make sure you get the same snp hits
#data2_add@gtdata@snpnames[193590]#[20836]
#sel_gene_snp_fix<-as.numeric(data2_add[,sel_resid_snps_fix])
#View(sel_resid_genotype_num_fix)
#write.table(sel_resid_genotype_num_fix,'sel_resid_genotype_fix.txt',sep='\t')


# GENE SYMBOL: had to download the sel10_expression list []1
#sel_resid_gene_fix<-meq_sel_fix$all$eqtls$gene
#sel_resid_gene_fix<-as.character(sel_resid_gene_fix)
#sel_resid_gene_fix<-as.numeric(sel_resid_gene_fix)
#sel_resid_symbol_fix<-sel_symbol[sel_resid_gene_fix]
#write.table(sel_resid_symbol_fix,'sel_symbol_hit_fix.txt',sep="\t")



#GENE EXPRESSION
#sel_resid_exp_fix<-sel_top_exp[sel_resid_gene_fix,]
#write.table(sel_resid_exp_fix,'sel_resid_exp_fix.txt',sep="\t")
#GENE FIT
#sel_fit_hits_fix<-sel_top_fit[sel_resid_gene_fix,]
#write.table(sel_fit_hits_fix,'sel_fit_hits_fix.txt',sep="\t")




