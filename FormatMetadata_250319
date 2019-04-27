#metadata csv from fasta defline + host-order barplot
#250319
#J baxter

#set wd
setwd('~/Desktop/fastas_for_raxml/')

#import fasta
library(seqinr)
FastaFiles = list.files(path = '.')
FastaFiles = FastaFiles[grep('19.fas' , FastaFiles)]
seq.list = lapply(FastaFiles , function(x) read.fasta(x , seqtype = 'AA' , as.string = TRUE))
seq.names = lapply(seq.list , sapply, getAnnot)
names= unlist(lapply(FastaFiles , function(x) {
  split = unlist(strsplit(x , '[/,_]+'))
  name = split[1]
  return(name)
}))

#Create trait df
fasta.annotations = lapply(seq.list , getAnnot)
edit.annotations = lapply(fasta.annotations , sapply , function(x){
  annot1 = sub('>', '' , x)
  annot = gsub('_' , '|' , annot1)
  return(annot)
})

split.annot = lapply(edit.annotations , sapply, function(x){strsplit( x , '|' , fixed = TRUE)})
split.annot.tr = lapply(split.annot, purrr::transpose)

#dates to yyyy-MM-dd
dates.list = list()
for (i in 1:length(split.annot.tr)){ 
  df =  split.annot.tr[[i]] 
  dates.list[[i]] = df[[9]]
  df[[9]] = NULL
}

formatted_dates.list = lapply(dates.list ,  sapply , function(x){
  date = gsub('_', '-' , x)
})

edit.annot.tr = list()
for (i in 1:length(split.annot.tr)){ 
  edit.annot.tr[[i]] = c(split.annot.tr[[i]] , as.list(formatted_dates.list[i] ))
}
#rejoin dates and trait lists


annot.as.df = function(x){
  e = do.call(rbind, Map(data.frame, GenBank.Accession= x[[1]], Virus.Group = x[[2]],
                         Virus.Genus = x[[3]] , Virus.Species = x[[4]],
                         Host= x[[7]], Country = x[[8]], Collection.Date= x[[9]]))
              return(e)
              }

annotation_df = lapply(edit.annot.tr  , annot.as.df)

lapply(annotation_df, function(x){
  colnames(x) = c('GenBank.Accession' , 'Virus.Group' , 
                  'Virus.Genus' , 'Virus.Species' ,
                  'Host' ,  'Country' ,
                  'Collection.Date')
  return(x)
}
)

#rownames = annotation
trifle = function(x){
  df_args <- c(x, sep='|')
  rn = do.call(paste, df_args)
  rownames(x) = rn
  return(x)
}

named_annotation_df = lapply(annotation_df, trifle)
#Host Order

NA_chars = c('NA' , 'Unknown')

clean = function(x){
  a = x[! x$Host %in% NA_chars, ]
  b = a[! a$Country %in% NA_chars, ]
  c = b[! b$Collection.Date %in% NA_chars,]
  return(c)
}


na_removed.list = lapply(annotation_df, clean)
names(na_removed.list) <- names


##group phyla prior to stratification##

#form host table
host.df = lapply(na_removed.list , function(x) as.data.frame(x$Host , stringsAsFactors = FALSE))

host.tables = lapply(na_removed.list , function(x) as.data.frame(table(x$Host) , stringsAsFactors = FALSE))


library(stringr)
library(stringr)
group = function(i){
  x2 = lapply(i , function(x)str_replace_all(x ,regex('Mosquito|Midge|Tick|Insect|Fly(?!catcher)|Sandfly|culex.*|Sand_Fly|Hard_Bodied_Tick|Arachnid|Bug', ignore_case = TRUE) , 'Arthropoda'))
  x3 = lapply(x2 , function(x)str_replace_all(x ,regex('Mouse|Vole|Rat(?!.*)|Rhabdomys.*|Rice_Rat|Peromyscus.*|Microtus.*|Mesocricetus-auratus|Mesocricetus|Sciurus-griseus|Lophuromys.*|woodchuck|Hylomyscus|Cricetulus.*|Hamster|Mus-musculus|murine|Myodes.*|Marmota-monax|Rattus.*|Sigmodon.*|Sciurus-variegatoides|Guinea_Pig|Gerbil|Squirrel|Rodent|Hamster|Jerboa', ignore_case = TRUE) , 'Rodentia'))
  x4 = lapply(x3 , function(x)str_replace_all(x ,regex('domestic-pig|Pig(?!eon)|Swine.*|wild-boar|domesticated-swine|backyard-pig|Boar|porcine|Sus-scrofa-domesticus|Sus-scrofa.*', ignore_case = TRUE) , 'Artiodactyla:Suinamorpha'))
  x5 = lapply(x4 , function(x)str_replace_all(x , regex('Cattle|Buffalo|Yak|Bovine|Cow(?!bird)|Buffalo|bovine|Bos-taurus|calf' , ignore_case = TRUE), 'Artiodactyla:Ruminantiamorpha'))
  x6 = lapply(x5 , function(x)str_replace_all(x ,regex('Goat_Antelope', ignore_case = TRUE) , 'Artiodactyla:Ruminantiamorpha'))
  x7 = lapply(x6 , function(x)str_replace_all(x ,regex('Sheep|.*Goat(?!_)|Antelope|Ovis-aries|Dama-dama|capra-hircus|Rupicapra.*' , ignore_case = TRUE), 'Artiodactyla:Ruminantiamorpha'))
  x8 = lapply(x7 , function(x)str_replace_all(x ,regex('Horse.*|Donkey.*|Mule|Wild_Horse|Equus.*|equine.*|equidae' , ignore_case = TRUE), 'Perissodactyla'))
  x9 = lapply(x8 , function(x)str_replace_all(x ,regex('.*Dog.*|.*Wolf|Fox|Jackal.*|Coyote.*|.*Bear.*|Canis-lupus-familiaris|canis-lupus|gray-wolf', ignore_case = TRUE) , 'Carnivora:Caniformia'))
  x10 = lapply(x9 , function(x)str_replace_all(x ,regex('Cat(?!.*)|.Wild_Cat|Civet.*|Tiger|Cheetah|feline|felis.*|cats|(?!)cat(?!)|Mustela-vison', ignore_case = TRUE) , 'Carnivora:Feliforma'))
  x11 = lapply(x10 , function(x)str_replace_all(x ,regex('.*Skunk|.*Raccoon|Badger|Melogale-moschata|Mink|Ferret|Weasel|Polecat|Procyon-lotor|Nasua-narica', ignore_case = TRUE) , 'Carnivora:Caniformia'))
  x12 = lapply(x11 , function(x)str_replace_all(x ,regex('Deer|.*Gazelle', ignore_case = TRUE) , 'Artiodactyla:Ruminantiamorpha'))
  x13 = lapply(x12 , function(x)str_replace_all(x ,regex('Mongoose|Hyena|Felis-Carnivora:Feliformaus|Paradoxurus-zeylonensis' , ignore_case = TRUE), 'Carnivora:Feliforma'))
  x14 = lapply(x13 , function(x)str_replace_all(x ,regex('Snail', ignore_case = TRUE) , 'Gastropoda'))
  x15 = lapply(x14 , function(x)str_replace_all(x ,regex('Blue_Jay|Crow|Corvus.*|Blackbird|Red_Winged_Blackbird|Zonotrichia-leucophrys|Grackle|Sparrow(?!hawk)|Jay|Cyanocitta-cristata|Goldfinch|Finch|Kite|Brown_Headed_Cowbird|Blue_Tit|Flycatcher|Estrilda-troglodytes|Indigo_Bunting|Pipit|Chickadee|Vireo|Starling|Titmouse|Cardinal|Raven|Mockingbird|Magpie|Weaver(?!_bug)|Warbler|Swallow|Myna|Thrush|Rook', ignore_case = TRUE),
                                                'Passeriformes' ))
  x16 = lapply(x15 , function(x)str_replace_all(x ,regex('Hawk|Sparrowhawk|Northern_Goshawk|Eagle|Goshawk|Vulture|Accipitridae|Buteo-jamaicensis' , ignore_case = TRUE), 'Accipitriformes'))
  x17 = lapply(x16 , function(x)str_replace_all(x ,regex('Parrot|Parakeet|Psittacus-erithacus|Nymphicus-hollandicus|Primolius-auricollis|Canindae-Macaw|Cacatua.*|Ara-macao|Aratinga.*|Amazona-ventralis|Ara-ararauna|Nestor.*|Diopsittaca-nobilis|budgerigar', ignore_case = TRUE) , 'Psittaciformes'))
  x18 = lapply(x17 , function(x)str_replace_all(x ,regex('Woodpecker' , ignore_case = TRUE), 'Piciformes'))
  x19 = lapply(x18 , function(x)str_replace_all(x ,regex('Owl|Strix-varia' , ignore_case = TRUE) ,'Strigiformes'))
  x20 = lapply(x19 , function(x)str_replace_all(x ,regex('Flamingo|American_Flamingo', ignore_case = TRUE) , 'Phoenicopteriformes'))
  x21 = lapply(x20 , function(x)str_replace_all(x ,regex('Gull|Tern|Shorebird|Auk|Common_Snipe|Stint|Seabird|Calidris-ruficollis|Larus', ignore_case = TRUE) , 'Charadriiformes'))
  x22 = lapply(x21 , function(x)str_replace_all(x ,regex('Swan|.*Goose|.*Duck|Egyptian_Goose|Shelduck|Dabbling_Duck|Teal|Green_Winged_Teal|Anatidae|Canada-Goose|Anas-.*|mallard.*|Serinus.*|wigeon|pintail.*|Anser-sp.*|Alopochen.*|.*Pochard', ignore_case = TRUE) , 'Anseriformes'))
  x23 = lapply(x22 , function(x)str_replace_all(x ,regex('Merlin|Falcon|Common_Kestrel' , ignore_case = TRUE), 'Falconiformes'))
  x24 = lapply(x23 , function(x)str_replace_all(x ,regex('Pelican|Heron|Ibis(?!es)|Ibises|Egret', ignore_case = TRUE) , 'Pelicaniformes'))
  x25 = lapply(x24 , function(x)str_replace_all(x ,regex('Dove|Pigeon|Columba', ignore_case = TRUE) , 'Columbiformes'))
  x26 = lapply(x25 , function(x)str_replace_all(x ,regex('Llama|Camel.*|Alpaca|Vicugna.*|dromedary' , ignore_case = TRUE), 'Artiodactyla:Camelidamorpha'))
  x27 = lapply(x26 , function(x)str_replace_all(x ,regex('Pronghorn(?!_)|Tragelaphus-strepsiceros' , ignore_case = TRUE), 'Artiodactyla:Ruminantiamorpha'))
  x28 = lapply(x27 , function(x)str_replace_all(x ,regex('Lizard|Rattlesnake', ignore_case = TRUE) , 'Squamata'))
  x29 = lapply(x28 , function(x)str_replace_all(x ,regex('Cormorant', ignore_case = TRUE) , 'Suliformes'))
  x30 = lapply(x29 , function(x)str_replace_all(x ,regex('.*Panda|canine|Canis-familiaris', ignore_case = TRUE) , 'Carnivora:Caniformia'))
  x31 = lapply(x30 , function(x)str_replace_all(x ,regex('Rabbit' , ignore_case = TRUE), 'Lagomorphia'))
  x32 = lapply(x31 , function(x)str_replace_all(x ,regex('Sloth|Anteater', ignore_case = TRUE) , 'Pilosa'))
  x33 = lapply(x32 , function(x)str_replace_all(x ,regex('Monkey|Macaque.*|Baboon|Howler_Monkey|Cebus-apella|Mandrillus-sphinx|Primate.*|Pan-troglodytes|marmoset|Macaca.*|.*colobus|Gorilla|Chimpanzee|Black_Howler|Callithrix.*|Lemur|Papio-ursinus|Cercopithecus-mitis', ignore_case = TRUE) , 'Primates'))
  x34 = lapply(x33 , function(x)str_replace_all(x ,regex('Opossum' , ignore_case = TRUE), 'Didelphimorphia'))
  x35 = lapply(x34 , function(x)str_replace_all(x ,regex('Chicken|Grouse|Turkey|Quail|Pheasant|Peafowl|Partridge|Poultry|gallus.*' , ignore_case = TRUE), 'Galliformes'))
  x36 = lapply(x35 , function(x)str_replace_all(x ,regex('Whale|Dolphin|Porpoise', ignore_case = TRUE) , 'Artiodactyla:Cetancodontamorpha'))
  x37 = lapply(x36 , function(x)str_replace_all(x ,regex('Mississippi-Sandhill-Crane|Crane|Coot|Rail' , ignore_case = TRUE), 'Gruiformes'))
  x38 = lapply(x37 , function(x)str_replace_all(x ,regex('Salmon|Fish' , ignore_case = TRUE), 'Fish'))
  x38 = lapply(x37 , function(x)str_replace_all(x ,regex('California_Sea_Lion|sea-lion|walrus|Seal|Phoca-vitulina' , ignore_case = TRUE), 'Carinvora:Pinnipedia'))
  x39 = lapply(x38 , function(x)str_replace_all(x ,regex('Vero_Cell|Parasitic_Oomycote|Nematode|Crfk_Cell|Environmental' , ignore_case = TRUE), 'Experimental'))
  x40 = lapply(x39 , function(x)str_replace_all(x ,regex('Chiroptera.*|Bat|Eidolon.*|Rhinolophus*|Coleura.*|fruit-bat|Hipposideros.*|Myotis.*|Miniopterus.*|Lasiurus.*|Mops.*|Chaerephon.*|Eptesicus.*|Artibeus.*|ptero.*|Murina.*|Tadarida.*|Nyctinomops.*|Epomophorus.*|Sturnira.*|Rousettus.*|Pipistrellus.*' , ignore_case = TRUE), 'Chiroptera'))
  x41 = lapply(x40 , function(x)str_replace_all(x ,regex('Homo-Sapiens|Human' , ignore_case = TRUE), 'Homo sapiens'))
  x42 = lapply(x41 , function(x)str_replace_all(x ,regex('Penguin' , ignore_case = TRUE), 'Sphenisciformes'))
  x43 = lapply(x42 , function(x)str_replace_all(x ,regex('Crocidura.*|Asian_White_Toothed_Shrew|Erinaceus.*|Tupaia-belangeri-chinensis|Shrew|Mole|Sorex-araneus' , ignore_case = TRUE), 'Eulipotyphla'))
  x44 = lapply(x43 , function(x)str_replace_all(x ,regex('Macropus-agilis' , ignore_case = TRUE), 'Marsupialia'))
  return(x44)
}


grouped_t1 = lapply(host.df , group)

corrections = function(i){
 y2 = lapply(i , function(x)str_replace_all(x ,regex('cat|south-china-Carnivora:Feliforma|palm-Carnivora:Feliforma', ignore_case = TRUE) , 'Carnivora:Feliforma'))
 y3 = lapply(y2 , function(x)str_replace_all(x ,regex('domesticated-Artiodactyla:Suinamorpha|backyard-Artiodactyla:Suinamorpha|backyard-Artiodactyla:Suinamorpha|domestic-Artiodactyla:Suinamorpha|wild-Artiodactyla:Suinamorpha|Artiodactyla:Suinamorpha-domestica|domestiCarnivora:Feliformaed-Artiodactyla:Suinamorpha', ignore_case = TRUE) , 'Artiodactyla:Ruminantiamorpha'))
 y4 = lapply(y3 , function(x)str_replace_all(x ,regex('OtoChiroptera|Chiroptera-BF|Chiroptera-Chiroptera|Chiroptera-blasii|Chiroptera-landeri|Chiroptera-Ferrumequinum|Chiroptera-sinicus', ignore_case = TRUE) , 'Chiroptera'))
 y5 = lapply(y4 , function(x)str_replace_all(x ,regex('Passeriformes-olivaceus', ignore_case = TRUE) , 'Passeriformes'))
 y6 = lapply(y5 , function(x)str_replace_all(x ,regex('Rodentia-migratorius', ignore_case = TRUE) , 'Rodentia'))
 y7 = lapply(y6 , function(x)str_replace_all(x ,regex('gray-Carnivora:Caniformia', ignore_case = TRUE) , 'Carnivora:Caniformia'))
 y8 = lapply(y7 , function(x)str_replace_all(x ,regex('Charadriiformes-species', ignore_case = TRUE) , 'Charadriiformes'))
 return(y8)
}


grouped_t = lapply(grouped_t1 , corrections)
groups.df = lapply(grouped_t , function(x) as.data.frame(x , stringsAsFactors = FALSE))

bind_groups = function(x,y){
  frame = cbind.data.frame(x, y)
  frame$Host = NULL
  return(frame)
  
}

groupedhost.df = mapply(bind_groups , x = na_removed.list , y = groups.df , SIMPLIFY = FALSE)
t = lapply(groupedhost.df , function(x) table(x[,7]))
t
lapply(groupedhost.df,function(x) colnames(x) = c('GenBank.Accession' , 'Virus.Group', 'Virus.Group','Virus.Species', 'Country' , 'Collection.Date' , 'Host'))
groupedhost.tables = lapply(groupedhost.df , function(i) as.data.frame(prop.table(table(i$x.Host)) , stringsAsFactors = FALSE))


minor.hosts = purrr::flatten(lapply(groupedhost.tables , function(x) subset(x , Freq<0.075 , select = Var1)))

major_host_traits = mapply(function(x,y) x[!x[,7] %in% y, ] , x = groupedhost.df , y = minor.hosts , SIMPLIFY = FALSE)
lapply(major_host_traits , function(x) table(x[,7]))


plot.list = lapply(major_host_traits ,function(x)table(x[,7]))
plot.df = reshape2::melt(plot.list)
plot.df = plot.df[order(plot.df$Virus.Group , decreasing = TRUE),]
colnames(plot.df) = c('Host.Order', 'Frequency' , 'Virus.Group')
library(ggpubr)

###stacked bar plot###
library(viridis)
p = ggbarplot(plot.df, x = "Virus.Group", y = "Frequency", 
          fill = "Host.Order" , color = 'Host.Order', palette = 'simpsons' , xlab = 'Virus Group' , ylab = 'Frequency', orientation= 'horiz')
p

major_host_traits = major_host_traits[order(names(major_host_traits))]
major_host_traits = lapply(major_host_traits , function(x){ data <- x[c(1,2,3,4,5,7,6)]
                           return(data)})

#####Fix and save

annotations = lapply(major_host_traits , rownames)
#write new fastas (same as before, but minus incidental hosts)


### write stratified fasta files

#generate subsamples using original fastadeflines (rownames of strat_traits)

fastas = purrr::flatten(seq.list)

subset.fasta= lapply(annotations , function(x) fastas[names(fastas) %in% x])

names(subset.fasta) = names

#write New Deflines

new.def = function(x){
  df_args <- c(x, sep ='|')
  rn = do.call(paste, df_args)
  return(rn)
}

new.deflines = lapply(major_host_traits, new.def)
names(new.deflines) = names
new.deflines$Lyssavirus

a = lapply(new.deflines , length) 
b = lapply(subset.fasta , length)

for (i in 1:length(a)){stopifnot(a[[i]] == b[[i]])}

#write files

Strat_FastaFiles = lapply(names , function(x) paste0('restricted_host_fasta/',x,'_nuc_aligned_', Sys.Date()))


lapply(1:length(subset.fasta), function(i){new.fasta = seqinr::write.fasta(getSequence(subset.fasta[[i]]) ,names = new.deflines[[i]] ,file.out = paste0(Strat_FastaFiles[i], ".fasta"))})

#write phylip
#to matrix
prot = lapply(subset.fasta, sapply, function(x) seqinr::translate(getSequence(x)))
mat =lapply(subset.fasta , function(x) t(matrix(unlist(x), length(x[[1]]), length(x))))
accesion = lapply(major_host_traits , function(x) as.character(x[,1]))
mat2 = mapply(function(x,y){ rownames(x) = y
       return(x)}, x = mat , y = accesion )


lapply(1:length(subset.fasta), function(i){phylip= write.dna(mat2[[i]] , format = 'interleaved', file = paste0(Strat_FastaFiles[i], ".phy"))})

#Append Epidemic Potential Level

#Export .Csv

csv.names= lapply(names , function(x) paste0('traits/',x,'_defline_stratifiedtraits_', Sys.Date()))

names(major_host_traits) <- csv.names

for (i in 1:length(major_host_traits)){
  rownames(major_host_traits[[i]]) = new.deflines[[i]]
}
#write .csv to file

lapply(1:length(major_host_traits), function(i) write.csv(major_host_traits[[i]], ,####
                                                     file = paste0(names(major_host_traits[i]), ".csv"),
                                                     row.names = TRUE))
