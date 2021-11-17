rm(list = ls())
source("FDPbound.R")


#================Data download and preprocessing=============
#The following data download and preprocessing code is partly adapted 
#from the example R script at https://hivdb.stanford.edu/pages/genopheno.dataset.html

#helper function: create the design matrix X with the input mutations/positions
buildX <- function(dat, mut, ps){
  X <- matrix(NA, nrow=nrow(dat), ncol=length(mut))
  # loop through all positions
  for(p in unique(ps)){
    p1 <- substr(dat[,p],1,1)  # first mutation at this position
    p2 <- substr(dat[,p],2,2)
    for(ind in which(ps==p)){
      X[,ind] <-  as.numeric(p1==as.character(mut[ind]) | 
                               p2==as.character(mut[ind]))  
    }
  }  
  colnames(X) <- paste0(ps,mut)  
  return(X)
}

#download the dataset for PI-type drugs from the website
dataset='PI'
drug='ATV'
#min.muts is the minimum number of sequences that a mutation must appear in.
min.muts=3 
dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt",
                  header=TRUE, sep="\t",stringsAsFactors=FALSE)
dat[dat=="."] <- NA
posu <- dat[,10:108]

#download PI complete mutations
PIcomplete=read.table("https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/MUTATIONLISTS/COMPLETE/PI",
                      header = F,sep="\t",stringsAsFactors=FALSE)
allPImut=c('test')
for (i in 1:nrow(PIcomplete)) {
  allPImut=append(allPImut,
                  paste0(PIcomplete[i,1],unique(toupper(unlist(strsplit(PIcomplete[i,2],split = ' ')))))
  )
}
allPImut=allPImut[-1]

muts.in=allPImut
# get the amino acids and positions for the mutations to be included in the model
mut <- ifelse(nchar(muts.in)==3,toupper(substr(muts.in,3,3)),
              toupper(substr(muts.in,2,2)))
ps <- suppressWarnings(ifelse(nchar(muts.in)==3,as.numeric(substr(muts.in,1,2)),
                              as.numeric(substr(muts.in,1,1))))

# construct design matrix for OLS
X <- buildX(posu, mut, ps)#1958 rows, 224 cols

# construct dependent variable
drugcol <- which(names(dat)==drug)    
Y <- as.numeric(dat[,drugcol])  # absolute measure
Ylog10 <- log10(Y)
df.log <- data.frame(Y=Ylog10, X=X)

# remove all rows with missing values
rem.rows <- unique(which(is.na(df.log),arr.ind=TRUE)[,1])
df.log.cc <- df.log[-rem.rows,]  # complete case
if(sum(is.infinite(df.log.cc$Y))>0){
  df.log.cc=df.log.cc[-which(is.infinite(df.log.cc$Y)),]
}

# remove mutations that are rare
rare.muts <- which(colSums(df.log.cc[,-1])<min.muts)
if(length(rare.muts)>0){
  message(paste0(muts.in[rare.muts],
                 " excluded from the model because it appears in fewer than ",
                 min.muts," sequences.\n"))
  df.log.cc <- df.log.cc[,-(rare.muts+1)]  
}
print(paste(drug,'dataset shape:',dim(df.log.cc)))#1083 rows, 211 cols

# check duplicated columns from X to allow for identifiability
X=as.matrix(df.log.cc[,-1])
if(length(which(duplicated(t(X))))>0){
  print('Warning! Duplicate columns!')
  return(NA)
}

# download TSM (the approximated ground truth) for PI-type of drugs
NPTSM_PI=read.table('https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/MUTATIONLISTS/NP_TSM/PI',
                    header = F,sep="\t",stringsAsFactors=FALSE)
allNPTSM_PImut=c('test')
for (i in 1:nrow(NPTSM_PI)) {
  allNPTSM_PImut=append(allNPTSM_PImut,
                        paste0(NPTSM_PI[i,1],unique(toupper(unlist(strsplit(NPTSM_PI[i,2],split = ' ')))))
  )
}
allNPTSM_PImut=allNPTSM_PImut[-1]
tp=paste0('X.',allNPTSM_PImut)#all NPTSM PI mutations



#=============Implement the FDX control methods==============
fdr.target=0.2# FDR<=0.2
fdp.target=0.8# Pr(PDP<=0.2)>=0.8

# number of candidate mutations in the dataset
p = ncol(df.log.cc)-1

#---run fdpc method
set.seed(123)
lambdas=qnorm(p=1-fdr.target/(2*p)*(1:p))#the lamgba_{bh} seqeuence
# index of the selected predictors
fdpc.selected=FDPcontrolFit(x=as.matrix(df.log.cc[,-1]), y=df.log.cc[,1],
                       lambdas=lambdas,fdr.target,FDP_target=0,sigma=NA,
                       version='fdpc',isOrthogonal=F,standardize=T)


print('selection results of fdpc method:')
c('total.num.of.discoveries'=length(fdpc.selected),
  'num.of.discoveries.in.TSM'=length(intersect(names(df.log.cc)[1+fdpc.selected],tp)),
  'num.of.discoveries.not.in.TSM'=length(setdiff(names(df.log.cc)[1+fdpc.selected],tp))
  )

# output:
# total.num.of.discoveries   num.of.discoveries.in.TSM     num.of.discoveries.not.in.TSM 
#        32                            28                             4 


#---run fdpc+
set.seed(123)
lambdas=qnorm(p=1-fdr.target/(2*p)*(1:p))#the lamgba_{bh} seqeuence
# index of the selected predictors
fdpcPlus.selected=FDPcontrolFit(x=as.matrix(df.log.cc[,-1]), y=df.log.cc[,1],
                           lambdas,fdr.target,FDP_target=fdp.target,sigma=NA,
                           version='fdpc+',isOrthogonal=F,standardize=T)

print('selection results of fdpc+ method:')
c('total.num.of.discoveries'=length(fdpcPlus.selected),
  'num.of.discoveries.in.TSM'=length(intersect(names(df.log.cc)[1+fdpcPlus.selected],tp)),
  'num.of.discoveries.not.in.TSM'=length(setdiff(names(df.log.cc)[1+fdpcPlus.selected],tp))
)

# output
# total.num.of.discoveries   num.of.discoveries.in.TSM     num.of.discoveries.not.in.TSM 
#         41                            33                             8 

