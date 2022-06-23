#!/usr/bin/env Rscript

print("please make sure to install specL R library")

library(specL)

args = commandArgs(trailingOnly=TRUE)

o_dir = args[1]

input.file=file.path(o_dir, "best.msms.for.RT.txt")
output.file=file.path(o_dir, "RT.txt")


#########################################################################################################################################################
# This program uses the retention times and caluclates for each raw file a general linear model based on retention time versus hydrophobicy index (HI). #
# HI is calculated based on the peptide sequence using specL R library. outliers are identified in each raw file, based on he lower and upper quantiles #
#########################################################################################################################################################


#Input file example (can be generates also by Arrange.msms.inputs.py):

#Raw file        RID     Sequence        Retention time
#Seq57470_QE3    Seq57470_QE3@0  AAAAAAAAAAA     35.647
#QEP2_13046_YSA_921_1_301220     QEP2_13046_YSA_921_1_301220@1   AAAAAAAAATN     56.869
#QEP2_13046_YSA_917_1_301220     QEP2_13046_YSA_917_1_301220@2   AAAAAAAAAW      60.989
#QEP2_13046_YSA_917_1_301220     QEP2_13046_YSA_917_1_301220@3   AAAAAAAAAW      61.342
#QEP2_13046_YSA_917_1_301220     QEP2_13046_YSA_917_1_301220@4   AAAAAAAAAW      61.251
#QEP2_13046_YSA_918_1_301220     QEP2_13046_YSA_918_1_301220@5   AAAAAAAAAW      60.831


elements = read.csv(input.file, stringsAsFactors=FALSE, check.names = FALSE,sep="\t")

#####################
# this function caculates outliers based on linear model residuals. values greater than 3rd quartile+IQR are outliers
find_outliers = function (data){
	lowerq = quantile(data)[2] # <25%
	upperq = quantile(data)[4] # >75%
	iqr = upperq - lowerq #Or use IQR(data)
	#identify extreme outliers
	extreme.threshold.upper = (iqr * 1.5) + upperq
	(data > extreme.threshold.upper)
	}
###################

print("RT: Cacl HI")
elements$HI = sapply(elements$Sequence, function(seq) ssrc(seq))  # calculate HI for each peptide sequence method_col = 'ssrc_2004'

rt_col = grep('^Retention*', names(elements))  # find RT column index
# for each raw file separately, since RT is sample-dependent
uniq.raw=unique(elements[,1])
for (iraw in 1:length(uniq.raw) ) {
	print(uniq.raw[iraw])
	group.rows = (elements[,1]==uniq.raw[iraw])
	df=elements[group.rows,]

	rt.col.name = names(df)[rt_col] = 'current_RT' #rename retention time column
	model = if (nrow(df)>20) { #if more than 20 data points for model
			abs(lm(HI~current_RT, data = df)$residuals) #linear regression model, and save only the residuals
			}
		else { rep(0,nrow(df))
			} #assign zeros if not enough peptides for analysis

	elements[group.rows,"RT_residuals"] = model 
	elements[group.rows,"RT_outliers"] = find_outliers(model)
	model=NULL
	}

write.table(elements, output.file,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE) #save
