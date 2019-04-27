#tree correlations
#J Baxter

library(ape)
library(seqinr)
library(ggpubr)
library(parallel)
library(doSNOW)

#Set current working directory
setwd('~/Desktop/fastas_for_raxml/test/')

#Sam's Functions
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


traitToFactor <- function(trait, other=c("Other","Z")) {
  tbl <- sort(table(trait),decreasing=TRUE)
  ll  <- rownames(tbl)
  inds <- match(ll,other)
  jj   <- which(!is.na(inds))
  ll   <- c(setdiff(ll, other[inds[jj]]), other[inds[jj]])
  trait<- factor(trait, levels=ll)
  return( trait )
}

cleanRownames <- function(tbl) {
  rn <- rownames(tbl)
  rn <- gsub("__","_",rn)
  rownames(tbl) <- rn
  return( tbl )
}

readTree <- function( fname ) {
  tr <- read.tree( fname )
  tr <- multi2di(tr)
  tr <- ladderize(tr)
  tr <- phytools::midpoint.root(tr)
  return( tr )
}

readProteinSeqs <- function( fname ) {
  protSeqs <- seqinr::read.fasta(fname, seqtype="AA")
  AATbl    <- t( matrix( unlist(protSeqs), length(protSeqs[[1]]), length(protSeqs) ) )
  rownames(AATbl) <- attributes(protSeqs)$names
  AATbl    <- toupper(AATbl)
  AATbl    <- cleanRownames(AATbl)
  return( AATbl )
}

ace_discrete <- function(trait=trait, tr=tr, model="ER", doTraitToFactor=TRUE) {
  if (doTraitToFactor) {
    trait <- traitToFactor(trait)
  }
  ace_model  <- ape::ace(trait, tr, type="discrete", model=model)
  max_anc    <- apply(ace_model$lik.anc, 1, which.max)
  anc_state  <- colnames(ace_model$lik.anc)[max_anc]
  anc_prob   <- ace_model$lik.anc[,max_anc]
  states     <- c(paste(trait), anc_state)
  probs      <- c(array(1,length(trait)), anc_prob)
  LL         <- logLik(ace_model)
  AIC        <- AIC(ace_model)
  rates      <- ace_model$rates
  rates_se   <- ace_model$rates_se
  return( list(states=states, probs=probs, LL=LL, AIC=AIC, rates=rates, rates_se=rates_se, ace_model=ace_model))
}

mappedStates <- function(ace_m=ace_m, tr=tr) {
  likanc <- ace_m$ace_model$lik.anc
  cn     <- colnames(likanc)
  ntips  <- length(tr$tip.label)
  liktip <- matrix(0, ntips, length(cn))
  colnames(liktip) <- cn
  for (j in 1:length(cn)) {
    tinds  <- which(ace_m$states[1:ntips]==cn[j])
    if (length(tinds)>0) {
      liktip[tinds,j] <- 1
    }
  }
  likanc <- rbind(liktip,likanc)
  return( likanc )
}

probAnyChange <- function(likanc=likanc, tr=tr) {
  fromNode <- tr$edge[,1]
  toNode   <- tr$edge[,2]
  psame    <- apply(likanc[fromNode,]*likanc[toNode,], 1, sum)
  pchange  <- 1-psame
  return(pchange)
}

mapTraitsToTree <- function(traits=traits, selinds=1:length(traits[1,]), tr=tr, model="ARD") {
  ace_res    <- apply(traits, 2, ace_discrete, tr=tr, model="ARD", doTraitToFactor=TRUE) 
  likanc_res <- lapply(ace_res, mappedStates, tr=tr)
  pchange_res<- lapply(likanc_res, probAnyChange, tr=tr)
  return(pchange_res)
}

plot.correlations = function(x,y){
  upper =  cbind.data.frame(x[which(x[,2] >= -log10(0.05)),], 'Upper')
  colnames(upper) = c('Site' , 'p.value' , 'Threshold')
  remainder = cbind.data.frame(x[which(x[,2] < -log10(0.05)),], 'Lower')
  colnames(remainder) = c('Site' , 'p.value', 'Threshold')
  signif.groups.df = rbind.data.frame(upper , remainder)
  signif.groups.df = signif.groups.df[order(signif.groups.df$Site , decreasing = FALSE),]
  labels = list(as.character(upper$Site))
  #plot primary sequence variablility
  library(ggpubr)
  p = ggscatter(signif.groups.df, 'Site' , 'p.value' , color = 'Threshold', palette = 'npg',
                ylab = '-log10(P.Value)' , xlab = 'Site' , label = 'Site', size = 1.5, ylim =c(0,17), label.select = list(criteria = "`y` > 1.30103") ,repel = TRUE)  
  p2 = p + geom_hline(yintercept = as.numeric(Threshold) , color = 'black' , linetype = 2 )
  
  return(p2)
} 

as.binarychange = function(x , threshold){
  new = list()
  for (i in 1:length(x)){
    if (x[i] >= threshold){new[i] = 1}
    else{new[i] = 0}}
  new = unlist(new)
  return(new)
}

bothchange  = function(x) {
  y = list()
  for (i in 1:nrow(x)){
    t = x[i,]
    if(t[1] & t[2] == 1){
      y[i] = 1
    }else{
      y[i] = 0
    }
  }
  return(unlist(y))
}

getchanges = function(x){
  y = list()
  for (i in 1:length(x)){
    if(x[i] == 1){
      y[[i]] = i
    }else{
      y[[i]] = NULL
    }
  }
  return(unlist(y))
}

##START##

#Import AA MSA
FastaFiles.pre = list.files(path = '/home/james/Desktop/fastas_for_raxml/restricted_host_fasta/protein_fasta/')
FastaFiles.pre = FastaFiles.pre[grep('.fas' , FastaFiles.pre)]
FastaFiles = unlist(lapply(FastaFiles.pre , function(x) paste0('/home/james/Desktop/fastas_for_raxml/restricted_host_fasta/protein_fasta/' , x)))
seqmat_list = lapply(FastaFiles , readProteinSeqs)

names = unlist(lapply(FastaFiles.pre , function(x) {
  split = unlist(strsplit(x , '[/,_]+'))
  name = split[1]
  return(name)
}))

names(seqmat_list) = names
seqmat_list = seqmat_list[-c(5,8,9,14,15,16,17,20,21)]

names = names[-c(5,8,9,14,15,16,17,20,21)]
#Import traits
TraitFiles.pre = list.files(path = '/home/james/Desktop/fastas_for_raxml/traits/')
TraitFiles.pre = TraitFiles.pre[grep('.csv' , TraitFiles.pre)]
TraitFiles = unlist(lapply(TraitFiles.pre , function(x) paste0('/home/james/Desktop/fastas_for_raxml/traits/' , x)))
traits = lapply(TraitFiles , function(x) read.csv(file = x , stringsAsFactors = FALSE))
names(traits) = names
traits = traits[-c(5,8,9,14,15,16,17,20,21)]

Host = lapply(traits , function(x) x$x.Host)

#Import raxml trees
TreeFiles.pre = list.files(path = '~/Desktop/fastas_for_raxml/restricted_host_fasta/raxml_tree')
TreeFiles.pre = TreeFiles.pre[grep('bestTree' , TreeFiles.pre)]
TreeFiles = unlist(lapply(TreeFiles.pre , function(x) paste0('/home/james/Desktop/fastas_for_raxml/restricted_host_fasta/raxml_tree/' , x)))
raxml_tree = lapply(TreeFiles, readTree)
names(raxml_tree) = names
raxml_tree = raxml_tree[-c(5,8,9,14,15,16,17,20,21)]

#Import Feature Importance
FIFiles.pre = list.files(path = '.')
FIFiles = FIFiles.pre[grep('_FI_' , FIFiles.pre)]
FI= lapply(FIFiles, function(x) read.csv(file = x , stringsAsFactors = FALSE))

SOI = lapply(FI , function(x) x$Residue)

#Truncated AA MSAs
SOI.mat = mapply(function(x,y) y[,x] , SOI , seqmat_list)
mapply(function(x ,y) rownames(x) = y$tip.label[order(y$tip.label, decreasing = FALSE)] , x = SOI.mat, y = raxml_tree)

#Ancestral Trait Reconstruction: Host
stopifnot(length(Host)==length(raxml_tree))
ace.host = mapply(function(x,y) ace_discrete(x , y, model="SYM", doTraitToFactor=TRUE), x = Host , y= raxml_tree , SIMPLIFY = FALSE)
likanc_res.host= mapply(function(x,y) mappedStates(x , y), x =ace.host , y= raxml_tree , SIMPLIFY = FALSE)
pchange_res.host= mapply(function(x,y) probAnyChange(x , y), x = likanc_res.host , y= raxml_tree , SIMPLIFY = FALSE)

#Ancestral Trait Reconstruction: AA sites
start = Sys.time()
print(start)
cl = makeCluster(6)
clusterExport(cl = cl, c('ace_discrete' , 'SOI.mat' , 'raxml_tree' , 'traitToFactor'), envir = .GlobalEnv)

ace.residues =  parallel::clusterMap(cl = cl ,function(x,y) apply(x, 2, ace_discrete, y, model="SYM", doTraitToFactor=TRUE) , x = SOI.mat, y= raxml_tree , SIMPLIFY = FALSE)

stopCluster(cl)
remove(cl)
duration = Sys.time() - start
print(duration)

likanc_res.residues = mapply(function(x,y) lapply(x , function(i) mappedStates(i , y)), x =ace.residues, y= raxml_tree , SIMPLIFY = FALSE)
pchange_res.residues = mapply(function(x,y) lapply(x , function(i) probAnyChange(i , y)), x = likanc_res.residues , y= raxml_tree , SIMPLIFY = FALSE)


#correlate AA mutations and host switch events
binary.host = lapply(pchange_res.host , as.binarychange ,0.7)
binary.site = lapply(pchange_res.residues , sapply, as.binarychange ,0.7)#NB cols = sites , rows = Nodes

compchanges = function(x,y){
  host = x
  sites = y
  chisq.assoc =apply(sites , 2 , function(i) -log10(chisq.test( host, i )$p.value))
  return(chisq.assoc)
}

chisq.results = mapply(compchanges , x = binary.host , y = binary.site , SIMPLIFY = FALSE)


signif.df = mapply(cbind.data.frame, SOI , chisq.results , SIMPLIFY = FALSE)
                
sigtest = lapply(signif.df , function(x)nrow(x[which(x[,2] >= -log10(0.05)),]))

signifsets = list()
for (i in 1:length(sigtest)){
  if(sigtest[[i]] == 0){
    signifsets[[i]] =NULL
  }else{
    signifsets[[i]] = signif.df[[i]]
  }
}

#filter out datasets with no significant assoication
names(signifsets) = names
signifsets[sapply(signifsets, is.null)] <- NULL

#remove Enterovirus & Lyssavirus
signifsets[[9]] = NULL

#plot logvalues
logplots = lapply(signifsets, plot.correlations)

ggarrange(logplots[[1]], logplots[[2]] , 
          logplots[[3]],logplots[[4]],
          logplots[[5]],logplots[[6]],
          logplots[[7]],logplots[[8]],
          logplots[[9]],ncol = 3 , nrow = 3,labels='AUTO' , common.legend = TRUE)

#Save sites and log scores to file (for functional domain ivestigation)

