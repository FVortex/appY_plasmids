load('/home/jane/Документы/Misha/Transcription_origin/df_7_phages_with_phiols.rda')

phiols_seqs <- c()
phiors_seqs <- c()
for (i in 1:nrow(df_7_phages_with_phiols)) {
  l <- df_7_phages_with_phiols$`PhiOL coordinate`[[i]]
  r <- df_7_phages_with_phiols$`PhiOR coordinate`[[i]]
  
  phiols_seqs <- c(phiols_seqs, paste0(df_7_phages_with_phiols$Sequence[[i]][(l-15):(l+7)], collapse = ''))
  phiors_seqs <- c(phiors_seqs, paste0(df_7_phages_with_phiols$Sequence[[i]][(r-15):(r+7)], collapse = ''))
}

library(ggseqlogo)
library(ggplot2)


svg('/home/jane/Рабочий стол/logos_7_phiols_phiors.svg', width = 8, height = 4)
pdf('/home/jane/Рабочий стол/logos_7_phiols_phiors.pdf', pointsize = 14, width = 12, height = 8)

L <- list(phiols_seqs, phiors_seqs, c(phiols_seqs, phiors_seqs))
#names(L) <- c('Промоторы phiOL', 'Промоторы phiOR' , 'Промоторы phiOL и phiOR')
names(L) <- c('A', 'Б' , 'В')
ggseqlogo(L,  ncol = 1)+ scale_y_continuous(name="Информационнное содержание (бит)")+
scale_x_discrete(name ="Последовательность относительно ТСТ (п.н.)", limits=1:23, labels = -17:5)
dev.off()

seqs <- sapply(dataset_pro, function(x) {x$seq})
which.forward <- which(lapply(dataset_pro, function(x) {x$strand})=='forward')
forward_seqs <- seqs[which.forward]

library(stringdist)

dists <- c("osa", "lv", "dl", "hamming", "lcs", "qgram",
  "cosine", "jaccard", "jw", "soundex")
ranges_dists <- c()
for (i in dists) {
  ranges_dists <- c(ranges_dists, range(stringdistmatrix(forward_seqs, method = i)))
}

forward_seqs_200nts <- sapply(forward_seqs, function(x) {substr(x, 100, 300)})
summary(stringdistmatrix(forward_seqs_200nts, method = 'lv'))

stringdistmatrix()


#whole genome

library(ape)
U00096.2 <- read.GenBank(access.nb = 'U00096.2', as.character = T)
NC_001604 <- read.GenBank(access.nb = 'NC_001604', as.character = T)
NC_001604 <- 
  load('/home/jane/Документы/Misha/NC_001604.1.rda')
t7_NC_001604.1 <- as.character(NC_001604.1)
library(stringdist)
library(stringi)
l_dist <- stringdistmatrix(phiols_seqs, method = 'dl')
r_dist <-stringdistmatrix(phiors_seqs, method = 'dl')
lr_dist <-stringdistmatrix(L[[3]], method = 'dl')

library(mefa)
image(l_dist)
image(r_dist)
image(lr_dist)


int <- 50
slidists <- c()
for (i in - int:(nchar(t7_NC_001604.1)-2*int)) {
  tmp <- substr(t7_NC_001604.1, i+int, i+2*int)
  slidists <- c(slidists, stringdist(tmp, stri_reverse(tmp), method = 'dl'))
}

summary(slidists)
tail(order(slidists),  100)

plot(1:nchar(t7_NC_001604.1),ty='n')
abline(v = tail(order(slidists),  which(sort(slidists)==15)[1]))

#better

promoters_t7 <- read.table('/home/jane/Документы/Misha/Lab/mol_a-16/t7/t7_rmd/promoters.txt', sep = '\t', header = T)
terminators <- c(7588, 24210)
plot(1:nchar(t7_NC_001604.1),ty='n', ylim =range(slidists))
abline(v = c(promoters_t7$TSS), lwd =3, lty = 3, col ='red')
abline(v = c(terminators), lwd =3, lty = 3, col ='blue')

#
library(zoo)
lines(rollapply(slidists, 2000, mean), type = 'l')
splitten <- split(slidists, f = gl(n = ceiling(length(slidists)/5000), k = length(slidists)/ceiling(length(slidists)/5000) ))

barplot(sapply(splitten, mean))
#
abline(v = which(slidists<18), col = 'royalblue')
abline(v = which(slidists<17), col = 'darkgreen')

abline(v = which(slidists<16), col = 'orange')
abline(v = which(slidists<15), col = 'tomato')

plot(sort(slidists),ty='l')
abline(h = 7)
abline(v = which(sort(slidists)==5),lwd = .1)
plot(slidists, type ='l')
plot(slidists[5500:6500], type ='l')

library(clusterProfiler)

load('/home/jane/Документы/Misha/mydf_genes_100_promoters_each_group_for_go_analysis.rda')

mydf_genes1 <- c()
mydf_genes1$Entrez <- as.character(mydf_genes$Entrez)
mydf_genes1$group <- as.character(mydf_genes$group)

genes_formula <- compareCluster(Entrez~group, data=mydf_genes1, fun='enrichGO', OrgDb='org.EcK12.eg.db')

genes_formula <- compareCluster(as.character(mydf_genesEntrez[1:100])~c(rep(1,50), rep(2,50)),mydf_genes, fun='enrichGO', OrgDb='org.EcK12.eg.db')

library("org.EcK12.eg.db")

gene=  as.character(mydf_genes$Entrez[1:15])

unlist(mget(x=gene))


#gene names converted to entrezids using David online tool
entrezids <- read.csv('https://david.ncifcrf.gov/data/download/conv_B7A961D52FA51502872079731.txt', sep ='\t')

intersect(unique(mydf_genes1$Entrez), unique(entrezids$From))
length(unique(mydf_genes1$Entrez))
length(unique(entrezids$From))
sort(mydf_genes1$Entrez)
sort(entrezids$From)

#to substutute names with ids
replace(seq(20), c(4,3,5), c(2,4,6))


vec <- as.character(mydf_genes$Entrez)
bitr(vec, fromType = 'GENENAME', toType = 'ENTREZID',  OrgDb='org.EcK12.eg.db')
library(plyr)

new_vec <- mapvalues(vec, from=entrezids$From, to=entrezids$To)
##new_vec <- replace(vec, as.character(unique((entrezids$From))), as.character(unique((entrezids$To))),)

mydf_genes$Entrezids <- new_vec

genes_formula_orgdb <- compareCluster(mydf_genes$Entrezids~as.character(mydf_genes$group), data=mydf_genes, fun='enrichGO', OrgDb='org.EcK12.eg.db')
summary(genes_formula)
genes_formula_default_ecolik12 <- compareCluster(mydf_genes$Entrezids~as.character(mydf_genes$group), data=mydf_genes, fun='enrichGO', organism = 'ecolik12')


genes_formula_orgdb <- compareCluster(mydf_genes$Entrezids~as.character(mydf_genes$group), data=mydf_genes, fun='enrichKEGG')
ggo <- groupGO(gene     = new_vec[1:58],
               OrgDb    = org.EcK12.eg.db,
               ont      = "MF",
               level    = 3,
               readable = TRUE)
View(ggo)
#other function

enrichGO(mydf_genes$Entrezids, OrgDb = 'org.EcK12.eg.db', keytype = "ENTREZID", ont = "MF")
ego <- enrichGO(gene          = new_vec[1:58],
               # universe      = names(geneList),
                OrgDb         = org.EcK12.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

#trying to repcroduce 
genes_old1 <- as.character(read.table('/home/jane/Загрузки/genes_to_promoters_from_1_cluster_on_aeos1.txt', sep = '\t')$V1)
genes_old2 <- as.character(read.table('/home/jane/Загрузки/genes_to_promoters_from_2_cluster_on_aeos1.txt', sep = '\t')$V1)

dfold <- c()
dfold$geneClusters = c(genes_old1, genes_old2)
dfold$group = c(rep('a', length(genes_old1)), rep('b', length(genes_old2)))
compareCluster(dfold$geneClusters, dfold$group,data = dfold, fun = 'enrichGO',  OrgDb = 'org.EcK12.eg.db')
               

#another translation of gene names

ids <- bitr(as.character(mydf_genes$Entrez), fromType='GENENAME', toType='ENTREZID',OrgDb = 'org.EcK12.eg.db')
                 

library(GOSemSim)
cluster1 <- mydf_genes$Entrezids[which(mydf_genes$group=='genes_004_005')]
cluster2 <- mydf_genes$Entrezids[which(mydf_genes$group=='genes_005_006')]
cluster3 <- mydf_genes$Entrezids[which(mydf_genes$group=='genes_006_007')]
cluster4 <- mydf_genes$Entrezids[which(mydf_genes$group=='genes_007_008')]



clusts1_2 <- clusterSim(cluster1, cluster2, ont = 'BP', organism = 'ecolik12')
##clusterSim(genes_old1, genes_old2, ont = 'BP', organism = 'ecolik12')
sims <- mclusterSim(list(cluster1, cluster2, cluster3, cluster4), organism = 'ecolik12', ont = 'BP')
image(sims)


load('/home/jane/Загрузки/mydf_genes_keeping300proms.rda')

all_entrezids <- read.csv('https://david.ncifcrf.gov/data/download/conv_B7A961D52FA51502872079731.txt', sep ='\t')

#from david

david300pros <- read.table('https://david.ncifcrf.gov/data/download/conv_20FD27233EBA1503049038597.txt', header = T, sep ='\t')
#to substutute names with ids
save(david300pros, file = '/home/jane/Документы/Misha/david300pros.rda')

vec <- as.character(mydf_genes$Entrez)
#bitr(vec, fromType = 'GENENAME', toType = 'ENTREZID',  OrgDb='org.EcK12.eg.db')

intersect(vec, david300pros$From)
which(vec%in% david300pros$From)

vec_fix <- vec
vec_fix[which(!vec_fix%in% david300pros$From)] <- '' #not presented are removed
library(plyr)

new_vec <- mapvalues(vec_fix, from=david300pros$From, to=david300pros$To)
##new_vec <- replace(vec, as.character(unique((entrezids$From))), as.character(unique((entrezids$To))),)

mydf_genes$Entrezids <- new_vec
genes_formula_orgdb <- compareCluster(as.character(mydf_genes$Entrezids)~as.character(mydf_genes$group), data=mydf_genes, fun='enrichGO', OrgDb='org.EcK12.eg.db')

summary(genes_formula)




library(ecoli2.db)
xx <- as.list(ecoli2ALIAS2PROBE)
x <- ecoli2ENTREZID


#activation energy data
which_004_005 <- which(mydf_genes$group=='genes_004_005')
which_005_006 <- which(mydf_genes$group=='genes_005_006')
which_006_007 <- which(mydf_genes$group=='genes_006_007')
which_007_008 <- which(mydf_genes$group=='genes_007_008')

entrezids1cl<-c()
for (i in vec[which_004_005]) {
  entrezids1cl<-c(entrezids1cl, xx[[i]])
  entrezids1cl<-na.omit(entrezids1cl)
}
entrezids1cl_done<- as.character(as.list(x[entrezids1cl]))


entrezids2cl<-c()
for (i in vec[which_005_006]) {
  entrezids2cl<-c(entrezids2cl, xx[[i]])
  entrezids2cl<-na.omit(entrezids2cl)
}
entrezids2cl_done<- as.character(as.list(x[entrezids2cl]))

entrezids3cl<-c()
for (i in vec[which_006_007]) {
  entrezids3cl<-c(entrezids3cl, xx[[i]])
  entrezids3cl<-na.omit(entrezids3cl)
}
entrezids3cl_done<- as.character(as.list(x[entrezids3cl]))


entrezids4cl<-c()
for (i in vec[which_007_008]) {
  entrezids4cl<-c(entrezids4cl, xx[[i]])
  entrezids4cl<-na.omit(entrezids4cl)
}
entrezids4cl_done<- as.character(as.list(x[entrezids4cl]))

diff <- setdiff(entrezids1cl_done, entrezids4cl_done)




mydf_aeos1 <- data.frame(Entrez=c(entrezids1cl_done, entrezids2cl_done, entrezids3cl_done, entrezids4cl_done),group = c(rep('1cl', length(entrezids1cl_done)),rep('2cl', length(entrezids2cl_done)),rep('3cl', length(entrezids3cl_done)),rep('4cl', length(entrezids4cl_done))) )
formula <- compareCluster(Entrez~group, data=mydf_aeos1, fun='enrichGO', OrgDb='org.EcK12.eg.db')


aeos1.formula <- compareCluster(Entrez~group, data=mydf_aeos1, fun='enrichGO', OrgDb='org.EcK12.eg.db')

aeos1.formula <- compareCluster(list(entrezids1cl_done, entrezids4cl_done), fun='enrichGO', OrgDb='org.EcK12.eg.db')

