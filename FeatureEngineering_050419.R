#feature engineering
#28.03.19
#05.04.19 - set functions to iterate over files
#J Baxter
setwd('~/Desktop/fastas_for_raxml')

library(ggpubr)
library(bio3d)
library(seqinr)
library(ape)
library(protr)
library(ggpubr)
library(FactoMineR)


#Define functions (including Sam's functions)

cleanRownames <- function(tbl) {
  rn <- rownames(tbl)
  rn <- gsub("__","_",rn)
  rownames(tbl) <- rn
  return( tbl )
}

readProteinSeqs <- function( fname ) {
  protSeqs <- seqinr::read.fasta(fname, seqtype="AA")
  AATbl    <- t( matrix( unlist(protSeqs), length(protSeqs[[1]]), length(protSeqs) ) )
  rownames(AATbl) <- attributes(protSeqs)$names
  AATbl    <- toupper(AATbl)
  AATbl    <- cleanRownames(AATbl)
  return( AATbl )
}

majorityAA <- function(AA, useInds=1:length(AA), exAA=c("-","X","?"), return.max=1) {
  incAA     <- AA[useInds]
  exInds    <- c()
  for (j in 1:length(exAA)) {
    exInds    <- c(exInds,which(incAA==exAA[j]))
  }
  incAA     <- incAA[ setdiff(1:length(incAA),exInds)]
  
  tbl       <- sort(table(incAA), decreasing=TRUE)
  if (length(tbl)==1) {
    orderedAA <- attributes(tbl)$names
  } else {
    orderedAA <- rownames(tbl)
  }
  
  
  if (return.max==1) {
    return( orderedAA[1] )
  } else {
    if ( return.max < 0) {
      return( orderedAA )
    } else {
      if ( return.max <= length(orderedAA) ) {
        return( orderedAA[1:return.max] )
      } else {
        return (orderedAA )
      }
    }
  }
  
}

calcEntropy <- function( AA, exAA=c("X","-","?"), useInds=1:length(AA) ) {
  incAA     <- AA[useInds]
  exInds    <- c()
  for (j in 1:length(exAA)) {
    exInds    <- c(exInds,which(incAA==exAA[j]))
  }
  siteData<- incAA[ setdiff(1:length(incAA),exInds)]
  tempTbl <- table(siteData)
  prob_i  <- tempTbl/length(siteData)
  plogp   <- prob_i*log(prob_i)/log(2)
  w       <- -sum(plogp)
  return( w )
}



upperquart <- function(x){
  Threshold = quantile(x)[2]
  scatter.df = cbind.data.frame(1:length(x) , x)
  colnames(scatter.df) = c('Position' , 'Entropy')
  lowerq = cbind.data.frame(subset(scatter.df , scatter.df$Entropy<=Threshold), 'Bottom Quartile')
  colnames(lowerq) = c('Position' , 'Entropy' , 'Threshold')
  remainder = cbind.data.frame(subset(scatter.df , scatter.df$Entropy>Threshold), 'Remainder')
  colnames(remainder) = c('Position' , 'Entropy' , 'Threshold')
  scatter.groups.df = rbind.data.frame(lowerq , remainder)
  scatter.groups.df = scatter.groups.df[order(scatter.groups.df$Position , decreasing = TRUE),]
  return(list(scatter.groups.df,Threshold))
}

entropyPlot <- function(x , name){
  p = ggscatter(x[[1]], 'Position' , 'Entropy' , labels = x$Position, color = 'Threshold', palette = 'npg',
                ylab = 'Entropy' , xlab = 'AA Position' , size = 1.5)+ geom_hline(yintercept = as.numeric(x[[2]]) , color = 'black' , linetype = 2 )
  
  #save plot to file (by name and sysdate)
  pdf(file = paste0(name, '_entropyplot_', Sys.Date(), '.pdf'))
  print(p)
  dev.off()
  
  return(p)
}

chisq.sigsites <- function(x,y,z){
  sites.df= x[[1]]
  split = split.data.frame(sites.df, sites.df$Threshold)
  remainder = split$`Remainder`
  sites = remainder$Position
  site.matrix = lapply(sites , function(i) y[,i])
  labled.mat = lapply(site.matrix , function(i) cbind(i , z))
  contin.table = lapply(labled.mat , function(i) ftable(xtabs(~. , i)))
  chi = lapply(contin.table  , function(i) chisq.test(i , simulate.p.value = TRUE))
  names(chi) = sites
  sig = unlist(lapply(chi , function(i) i$p.value))
  sig.ordered = sort(sig , decreasing = TRUE)
  significant.sites = as.numeric(names(sig[sig<=0.05]))
  significant.mat= y[,significant.sites]
  colnames(significant.mat) = significant.sites
  return(significant.mat)
}


mapVectors <- function(x){
  require(protr)
  require(FactoMineR)
  data(AAMolProp)
  pca = PCA(AAMolProp)
  pchem.vec = pca$ind
  vectors = c(pchem.vec$dist , 0 , 0 )
  letter = c(names(vectors[1:20]), '*', '-')
  
  for (i in 1:length(vectors)){
    x[x ==letter[i]] <- vectors[i]
  }
  return(x)
}

label_features <- function(featureset, Host){
  sigsites.df = cbind.data.frame(Host , featureset)
  colnames(sigsites.df) = c('Host' , colnames(featureset))
  return(sigsites.df)
}


##Start##

##import fasta
FastaFiles = list.files(path = 'restricted_host_fasta/protein_fasta/')
FastaFiles = FastaFiles[grep('.fas' , FastaFiles)]
files = unlist(lapply(FastaFiles , function(x) paste0('restricted_host_fasta/protein_fasta/', x)))
seqmat_list = lapply(files, readProteinSeqs)

names = unlist(lapply(FastaFiles , function(x) {
  split = unlist(strsplit(x , '[/,_]+'))
  name = split[1]
  return(name)
}))
names(FastaFiles) = names

##import traits
TraitFiles = list.files(path = '/home/james/Desktop/fastas_for_raxml/traits/')
TraitFiles = TraitFiles[grep('.csv' , TraitFiles)]
tfiles= unlist(lapply(TraitFiles  , function(x) paste0('/home/james/Desktop/fastas_for_raxml/traits/', x)))
traits = lapply(tfiles , function(x) read.csv(file = x , stringsAsFactors = FALSE))
names(traits) = names
Host = lapply(traits , function(x) x$x.Host)

#calculate entropy
entropytest = lapply(seqmat_list  , function(x) apply(x, 2, calcEntropy))

upperq = lapply(entropytest , upperquart)

mapply(entropyPlot , upperq , names)

#test for assocation between residue state and host trait
chisite = mapply(chisq.sigsites , upperq, seqmat_list, Host , SIMPLIFY = FALSE)


#map vectors to sequence matrix
sigsite.pchem = lapply(chisite , mapVectors)

#label feature sets with host trait
featureset = mapply(label_features, sigsite.pchem, Host, SIMPLIFY = FALSE)

#export feature set
csv.names = lapply(names , function(x) paste0(x,'_featureset_', Sys.Date()))

names(featureset) <- csv.names

lapply(1:length(featureset), function(i) write.csv(featureset[[i]], 
                                                   file = paste0(names(featureset[i]), ".csv"), row.names = FALSE))


##END##

#Plotting Principle Components of amino acid vectors (append to mapVectors if required)
groups = c(
  'Aliphatic' , 'Sulfur Containing', 'Acidic',
  'Acidic', 'Aromatic' ,'Aliphatic',
  'Basic' , 'Aliphatic' ,'Basic' ,
  'Aliphatic', 'Sulfur Containing','Acidic Amide',
  'Cyclic', 'Acidic Amide' , 'Basic', 
  'Hydroxyl Containing','Hydroxyl Containing','Aliphatic',
  'Aromatic', 'Aromatic')

aa.num = cbind.data.frame(pchem.vec$coord[,1:2] , groups)

colnames(aa.num) = c('PC1' , 'PC2' ,'Groupings')


ggscatter(aa.num  , x = 'PC1' , y = 'PC2' , xlab = 'Principal Component 1 (66.95%)' , 
          ylab = 'Principal Component 2 (20.68%)', color = 'Groupings' , palette = 'npg' ,shape = 4 , label = rownames(aa.num) )

