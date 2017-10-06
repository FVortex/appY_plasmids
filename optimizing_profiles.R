rm(list = ls())
#whole genome

library(ape)
library(seqinr)
library(reldna)
U00096.2 <- as.character(read.fasta('/home/jane/Документы/Misha/Genome_to_dataset_pro_U00096.2.txt', as.string = T)[[1]])
mpots_U00096.2 <- reldna::lseqspline1D.BP(s = U00096.2, width = 1, bound = c(1, nchar(U00096.2)), ref = 1)


library(zoo)

rollmean <- rollapply(mpots_U00096.2$mpot, 50, mean)
int <- (which.min(rollmean)-100):(which.min(rollmean)+100)

plot(mpots_U00096.2$mpot[int], type = 'l')
lines(rollmean[int], lty = 2 )

int_seq <- substr(U00096.2, min(int), max(int))

#all possible snps
int_seq_copy <- int_seq

snps <- c()
for (i in 1:nchar(int_seq)) {
  nuc <- substr(int_seq, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    snps <- c(snps, int_seq_copy)
    int_seq_copy <- int_seq
    
  }
}


snps_mpots <- c()
for (i in seq_along(snps)) {
  res <- lseqspline1D(snps[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  snps_mpots <- cbind(snps_mpots, res$mpot[inds])
}


plot(apply(snps_mpots, 2, sum))
plot(apply(snps_mpots, 2, min))
minimorum1 <- which.min(apply(snps_mpots, 2, min))
min_rollmean60 <- apply(snps_mpots, 2, function(x){min(rollapply(x, 60, mean))})
plot(min_rollmean60)
matplot((snps_mpots), type = 'l', lty = 1, col = topo.colors(ncol(snps_mpots)))

par(mfrow = c(2,4))
plot(dend, xlab = NULL)
matplot(snps_mpots[,which(cut ==1)], type = 'l')
matplot(snps_mpots[,which(cut ==2)], type = 'l')
matplot(snps_mpots[,which(cut ==3)], type = 'l')
matplot(snps_mpots[,which(cut ==4)], type = 'l')

res <- lseqspline1D(int_seq, bound = c(1,201), width = 1, ref = 100)
inds <- which(res$x%in%-300:300)
mpot <- res$mpot[inds]
x <- -300:300
dir.create('/home/jane/Документы/Misha/lowest_ep_ecoli_all_snps')

# #unlink('/home/jane/Документы/Misha/lowest_ep_ecoli_all_snps',recursive = T)
for (i in 1:ncol(snps_mpots)) {
  png(paste0('/home/jane/Документы/Misha/lowest_ep_ecoli_all_snps/',  i, '.png'),width = 1000, height = 750, res =200)
  plot(x, snps_mpots[,i], type = 'l', xlab = 'Sequence (angstrom)', ylab ='EP vaue', ylim = range(snps_mpots), main = i)
  lines(x, mpot, col = 'orange')
  
  dev.off()
}
#and bendability
snps_bends <- c()
for (i in seq_along(snps)) {
  res <- bendability(snps[i], bound = c(1,201), width = 1)
  snps_bends <- cbind(snps_bends, res)
}


plot(apply(snps_bends, 2, sum))
plot(apply(snps_bends, 2, min))
minimorum1 <- which.min(apply(snps_bends, 2, min))
min_rollmean60 <- apply(snps_bends, 2, function(x){min(rollapply(x, 60, mean))})
plot(min_rollmean60)
matplot(snps_bends, type = 'l', lty = 1, col = topo.colors(201))

init_bend <- bendability(substr(mpots_U00096.2$seq, min(int), max(int)),  bound = c(1,201), width = 1)
x <- -100:100
dir.create('/home/jane/Документы/Misha/bends_lowest_ep_ecoli_all_snps')
# #unlink('/home/jane/Документы/Misha/bends_lowest_ep_ecoli_all_snps', recursive = T)
for (i in 1:ncol(snps_mpots)) {
  png(paste0('/home/jane/Документы/Misha/bends_lowest_ep_ecoli_all_snps/', i, '.png'),width = 1000, height = 750, res =200)
  plot(x, snps_bends[,i], type = 'l', xlab = 'Sequence (angstrom)', ylab ='Bendability', ylim = range(snps_bends, na.rm = T), main = i)
  lines(x, init_bend, col = 'orange')
  
  dev.off()
}
library(fastcluster)
library(dendextend)

colfunc <- colorRampPalette(c('tomato', 'royalblue'))
colnames(snps_mpots) <- 1:603
hclusted <- hclust.vector(snps_mpots, method = 'ward')
hclusted %>% as.dendrogram %>% color_branches(k = length(hclusted$order), col = colfunc(length(hclusted$order))[hclusted$order])%>%raise.dendrogram(2)-> dend
plot(dend, xlab = NULL)


cut <- cutree(hclusted, k=4)




#system('cd /home/jane/Документы/Misha/lowest_ep_ecoli_all_snps/
#convert -delay 10 *.png ani.gif')
system('cd /home/jane/Документы/Misha/lowest_ep_ecoli_all_snps/

convert -delay 10 $(for i in $(seq 0 1 603); do echo to_gif${i}.png; done) -loop 0 animated.gif')


#next iteration with minimorum 375
snps1[minimorum1]
#all possible snps1
int_seq_copy <- snps1[minimorum1]

snps1 <- c()
for (i in 1:nchar(int_seq_copy)) {
  nuc <- substr(int_seq_copy, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    snps1 <- c(snps1, int_seq_copy)
    
  }
}


snps1_mpots <- c()
for (i in seq_along(snps1)) {
  res <- lseqspline1D(snps1[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  snps1_mpots <- cbind(snps1_mpots, res$mpot[inds])
}


plot(apply(snps_mpots, 2, min))
points(apply(snps1_mpots, 2, min), col =2)

plot(apply(snps_mpots, 2, sum))
points(apply(snps1_mpots, 2, sum), col =2)

minimorum2 <- which.min(apply(snps1_mpots, 2, min))

matplot(snps_mpots, type = 'l', col =1, lty =1, ylim = range(snps1_mpots))
matlines(snps1_mpots, type = 'l', col =2, lty =1, ylim = range(snps1_mpots))

apply(snps_mpots, 2, min)[minimorum]
apply(snps1_mpots, 2, min)[minimorum1]


apply(snps1_mpots, 2, min)[minimorum1]

#next iteration with minimorum 221
snps1[minimorum2]
#all possible snps
int_seq_copy <- snps1[minimorum2]

snps2 <- c()
for (i in 1:nchar(int_seq_copy)) {
  nuc <- substr(int_seq_copy, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    snps2 <- c(snps2, int_seq_copy)
    
  }
}


snps2_mpots <- c()
for (i in seq_along(snps2)) {
  res <- lseqspline1D(snps2[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  snps2_mpots <- cbind(snps2_mpots, res$mpot[inds])
}


plot(apply(snps_mpots, 2, min))
points(apply(snps1_mpots, 2, min), col =2)
points(apply(snps2_mpots, 2, min), col =3)

plot(apply(snps_mpots, 2, sum))
points(apply(snps1_mpots, 2, sum), col =2)
points(apply(snps2_mpots, 2, sum), col =3)

minimorum2 <- which.min(apply(snps2_mpots, 2, min))

matplot(snps_mpots, type = 'l', col =1, lty =1, ylim = range(snps2_mpots))
matlines(snps1_mpots, type = 'l', col =2, lty =1, ylim = range(snps2_mpots))
matlines(snps2_mpots, type = 'l', col =3, lty =1, ylim = range(snps2_mpots))

apply(snps_mpots, 2, min)[minimorum]
apply(snps2_mpots, 2, min)[minimorum1]


apply(snps2_mpots, 2, min)[minimorum1]



#bulding a giant loop

for (iter in 2:7) {

#next iteration with minimorum 375
assign(int_seq_copy, 
       get(paste0('snps', iter-1, '_mpots'))[get(paste0('minimorum', iter-1))])
#all possible snps1
#int_seq_copy <- get(paste0('snps', i+1))

tmp <- c()
for (i in 1:nchar(int_seq_copy)) {
  nuc <- substr(int_seq_copy, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    tmp<- c(tmp, int_seq_copy)
    
  }
}

assign(paste0('snps', iter), tmp)

tmp1 <- c()
for (i in seq_along(tmp)) {
  res <- lseqspline1D(tmp[i], bound = c(1,201), width = 1, ref = 100)
  inds <- which(res$x%in%-300:300)
  tmp1<- cbind(tmp1, res$mpot[inds])
}

assign(paste0('snps', iter, '_mpots'), tmp1)

assign(paste0('minimorum', iter), which.min(apply(tmp1, 2, min)))
}##
plot(apply(snps_mpots, 2, min))
points(apply(snps1_mpots, 2, min), col =2)

plot(apply(snps_mpots, 2, sum))
points(apply(snps1_mpots, 2, sum), col =2)

minimorum2 <- which.min(apply(snps1_mpots, 2, min))

matplot(snps_mpots, type = 'n', col =1, lty =1, ylim = range(snps1_mpots))
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col =2)
for (i in 1:7) {
matlines(get(paste0('snps', i, '_mpots')), type = 'l', col =i+1, lty =1, ylim = range(snps1_mpots))
print(mean(apply(get(paste0('snps', i, '_mpots')), 2, min)))
}

#matplot(snps_mpots, type = 'l', col =1, lty =1, ylim = range(snps1_mpots))


library(RColorBrewer)
colfunc <- brewer.pal(n = 7, name = 'Dark2')
svg('/home/jane/Документы/Misha/search_for_lowest_ep_snp.svg', width = 8, height = 10)
par(mfrow = c (4,2),
    mar = rep(0.5, 4))
matplot(get(paste0('snps', '_mpots')), type = 'l', col = colfunc[i], lty =1, ylim = range(snps1_mpots), main = 'E. coli, lowest EP region')
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col =2)

for (i in 1:7) {
  matplot(get(paste0('snps', i, '_mpots')), type = 'l', col = colfunc[i], lty =1, ylim = range(snps1_mpots), main = paste0('Lowest EP SNP ', i,' iteration'))
  abline(h = mean(mpots_U00096.2$mpot), lty = 3, col =2)
}
dev.off()

apply(snps_mpots, 2, min)[minimorum]
apply(snps1_mpots, 2, min)[minimorum1]


apply(snps1_mpots, 2, min)[minimorum1]

#to look at initial lowest EP
library(ape)
library(seqinr)
as.DNAbin(strsplit(int_seq, ''))
GC(strsplit(int_seq, '')[[1]])


library(zoo)

rollmean <- rollapply(mpots_U00096.2$mpot, 50, mean)


rollmean200 <- rollapply(mpots_U00096.2$mpot, 200, mean)

minN <- function(x, N=2){
  x <- -x
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

mins <- -minN(rollmean200, N = 1:200)

mins_pos <- sapply(mins, function(x) {which(rollmean200 == x)})

#EP profile itself
mins <- -minN(mpots_U00096.2$mpot, N = 1:200)

mins_pos <- sapply(mins, function(x) {which(mpots_U00096.2$mpot == x)})

plot(rollmean200, type ='l')
abline(v = mins_pos, col =2)

plot(mpots_U00096.2$mpot[(mins_pos[1]-50):(mins_pos[1]+50)], type = 'n', ylim = range(mpots_U00096.2$mpot))
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col = 2)

nos <- which(diff(sort(mins_pos))!=1)
sapply(mins_pos[nos], function(x) {
  lines(mpots_U00096.2$mpot[(x-50):(x+50)], type = 'l')
  }
)

#plot(mpots_U00096.2$mpot[(mins_pos-50):(mins_pos+50)], type = 'l')
#substr(U00096.2, start = mins_pos-50, stop = mins_pos+50)

substrings_lowest_ep <- sapply(mins_pos, function(x) {substr(U00096.2, start = x-50, stop = x+50)})
as.DNAbin(strsplit(substrings_lowest_ep, ''))

library(dplyr)
library(janeaustenr)
library(tidyr)
library(ngram)
d <- data_frame(txt = substrings_lowest_ep)
d
#  unnest_tokens(substrings_lowest_ep, text, token = "ngrams", n = 2) 
d %>%
  unnest_tokens(output = bigram, input = txt, token = "ngrams", n = 2) -> d_bigrams

d_bigrams %>%
  separate(bigram, c("word1", "word2"), sep = " ")
d_bigrams %>%
  separate(bigram, c("word1", "word2"), sep = " ") ->
  tidy_chekhov_bigrams

sa <- unlist(strsplit(substrings_lowest_ep, split = ''))
sa1 <- paste(sa, collapse = ' ')
ng2 <- ngram(sa1, n=2)
ng3 <- ngram(sa1, n=3)
ng4 <- ngram(sa1, n=4)
ng5 <- ngram(sa1, n=5)
ng6 <- ngram(sa1, n=6)

print(ng, output=  "full")
rbind(get.phrasetable(ng2), get.phrasetable(ng3), get.phrasetable(ng4),get.phrasetable(ng5),get.phrasetable(ng6))  

phr_ng6 <- get.phrasetable(ng6)
#artificial_seq <- paste0(unlist(strsplit(paste0(phr_ng6$ngrams[1:16], collapse = ''), split = ' ')), collapse = '')
artificial_seq <- paste0(unlist(strsplit(babble(ng6 , 201, seed =10), split = ' ')), collapse = '')
#artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)
artificial_output <- lseqspline1D(artificial_seq, bound = c(1,201), width = 1, ref = 100)
inds <- which(artificial_output$x%in%-300:300)
artificial_mpots <- artificial_output$mpot[inds]

plot(artificial_mpots, type = 'l', ylim = range(mpots_U00096.2$mpot))
#artificial against natural EP-lowest
sapply(mins_pos[nos], function(x) {
  lines(mpots_U00096.2$mpot[(x-300):(x+300)], type = 'l')
}
)

abline(h = mean(mpots_U00096.2$mpot), lty = 3, col = 'seagreen')
abline(h = min(mpots_U00096.2$mpot), lty = 3, col = 'royalblue')
lines()

#searching through generated seqs

artificial_mpots <- c()
for (i in 1:15000) {
  
artificial_seq <- paste0(unlist(strsplit(babble(ng6 , 201, seed =i), split = ' ')), collapse = '')
#artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)
artificial_output <- lseqspline1D(artificial_seq, bound = c(1,201), width = 1, ref = 100)
inds <- which(artificial_output$x%in%-300:300)
tmp <- artificial_output$mpot[inds]
artificial_mpots <- cbind(artificial_mpots, tmp)
}

boxplot(apply(artificial_mpots, 2, min))
boxplot(apply(artificial_mpots, 2, mean))

which_min_min <- which.min(apply(artificial_mpots, 2, min))
which_min_mean <- which.min(apply(artificial_mpots, 2, mean))

plot(artificial_mpots[,which_min_min], type ='l')
lines(artificial_mpots[,which_min_mean], type ='l', col =2)
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col = 'seagreen')
abline(h = min(mpots_U00096.2$mpot), lty = 3, col = 'royalblue')


#dynamical properties
dynchars<-function(seq, interval_size) {
  if (missing(seq))
    stop("Need to specify sequence (as a vector of chars)")
  
  if (missing(interval_size))
    stop("Need to specify interval size")
  
  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  
  seq<-toupper(seq)
  seq<-c(seq, seq[2:(interval_size)])
  
  a<-3.4*10^(-10)
  I<-c(7.6, 4.8, 8.2, 4.1)*10^(-44)
  K<-c(227, 155, 220, 149)*10^(-20)
  V<-c(2.09, 1.43, 3.12, 2.12)*10^(-20)
  tau<-c(127, 99, 140, 84)
  
  csA<-cumsum(seq=='A') 
  csT<-cumsum(seq=='T')
  csG<-cumsum(seq=='G')
  csC<-cumsum(seq=='C')
  
  countA = csA[interval_size:length(csA)]-c(0, csA[1:(length(csA)-interval_size)])
  countT = csT[interval_size:length(csT)]-c(0, csT[1:(length(csT)-interval_size)])
  countG = csG[interval_size:length(csG)]-c(0, csG[1:(length(csG)-interval_size)])
  countC = csC[interval_size:length(csC)]-c(0, csC[1:(length(csC)-interval_size)])
  
  M<-cbind(countA, countT, countG, countC)/interval_size
  M_comp<-cbind(countT, countA, countC, countG)/interval_size
  M_comp<-apply(t(M_comp),1,rev) 
  Is<-as.numeric(M%*%I)#! numeric conversion
  Ks<-as.numeric(M%*%K)
  Vs<-as.numeric(M%*%V)
  
  E01<-(8*(Ks*Vs)^0.5)* 6E23 / 4184
  d1<-((Ks*a^2)/Vs)^(0.5)/a;
  c1<-(Ks*a^2/Is)^0.5
  m1<-E01/c1/6.011E-26
  taus1<-as.numeric(M%*%tau) #!as.numeric conversion
  gc1 = M[,3] + M[,4]
  
  Is<-as.numeric(M%*%I)#! numeric conversion
  Ks<-as.numeric(M%*%K)
  Vs<-as.numeric(M%*%V)
  
  
  E02<- 8*(Ks*Vs)^0.5  * 6E23 / 4184;
  d2<-((Ks*a^2)/Vs)^(0.5)/a;
  c2<-(Ks*a^2/Is)^0.5;
  m2<-E02/c2/6.011E-26;
  taus2<-as.numeric(M_comp%*%tau)
  gc2 = M_comp[,3] + M_comp[,4]
  
  dynchars_return<-list(E01=E01, d1=d1, c1=c1, m1=m1, taus1=taus1, gc1=gc1, E02=E02, d2=d2, c2=c2, m2=m2, taus2=taus2, gc2=gc2)
  
  return(dynchars_return)
  
}

snps_e0 <- c()
for (i in seq_along(snps)) {
  res <- dynchars(unlist(strsplit(snps[i], split = '')), interval_size = 50)$E01
  snps_e0 <- cbind(snps_e0, res)
}

matplot((snps_e0), type = 'l', lty = 1, col = topo.colors(201))

init_e0 <- dynchars(unlist(strsplit(substr(mpots_U00096.2$seq, min(int), max(int)), split = '')), interval_size = 50)$E01
x <- -100:100
dir.create('/home/jane/Документы/Misha/e0_lowest_ep_ecoli_all_snps')
# #unlink('/home/jane/Документы/Misha/bends_lowest_ep_ecoli_all_snps', recursive = T)
for (i in 1:ncol(snps_e0)) {
  png(paste0('/home/jane/Документы/Misha/e0_lowest_ep_ecoli_all_snps/', i, '.png'),width = 1000, height = 750, res =200)
  plot(x, snps_e0[,i], type = 'l', xlab = 'Sequence (angstrom)', ylab ='E01', ylim = range(snps_e0, na.rm = T), main = i)
  lines(x, init_e0, col = 'orange')
  
  dev.off()
}

#aproteins

library(ngram)
library(seqinr)
prot <- read.fasta('/home/jane/Документы/Misha/5xwy.fasta.txt', as.string = T, seqonly = T)[[1]]
prot_spaces <- paste0(unlist(strsplit(prot, split = '')), collapse = ' ')
ng10 <- ngram(prot_spaces, n=10)

print(ng10, output=  "full")

phr_ng10 <- get.phrasetable(ng10)
(phr_ng10[1:100,])$ngrams

library(ggseqlogo)
p <- ggseqlogo((phr_ng10[1:100,])$ngrams,)

p+geom_logo(data =(phr_ng10[1:100,])$ngrams, method = 'p')


ng2 <- ngram(prot_spaces, n=2)

print(ng2, output=  "full")

phr_ng2 <- get.phrasetable(ng2)
(phr_ng2[1:100,])$ngrams

library(ggseqlogo)
p <- ggseqlogo((phr_ng2[1:100,])$ngrams,)

p+geom_logo(data =(phr_ng2[1:100,])$ngrams, method = 'p')



#dynamic properties



dynchars<-function(seq, interval_size) {
  if (missing(seq))
    stop("Need to specify sequence (as a vector of chars)")
  
  if (missing(interval_size))
    stop("Need to specify interval size")
  
  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  
  seq<-toupper(seq)
  seq<-c(seq, seq[2:(interval_size)])
  
  a<-3.4*10^(-10)
  I<-c(7.6, 4.8, 8.2, 4.1)*10^(-44)
  K<-c(227, 155, 220, 149)*10^(-20)
  V<-c(2.09, 1.43, 3.12, 2.12)*10^(-20)
  tau<-c(127, 99, 140, 84)
  
  csA<-cumsum(seq=='A') 
  csT<-cumsum(seq=='T')
  csG<-cumsum(seq=='G')
  csC<-cumsum(seq=='C')
  
  countA = csA[interval_size:length(csA)]-c(0, csA[1:(length(csA)-interval_size)])
  countT = csT[interval_size:length(csT)]-c(0, csT[1:(length(csT)-interval_size)])
  countG = csG[interval_size:length(csG)]-c(0, csG[1:(length(csG)-interval_size)])
  countC = csC[interval_size:length(csC)]-c(0, csC[1:(length(csC)-interval_size)])
  
  M<-cbind(countA, countT, countG, countC)/interval_size
  M_comp<-cbind(countT, countA, countC, countG)/interval_size
  M_comp<-apply(t(M_comp),1,rev) 
  Is<-as.numeric(M%*%I)#! numeric conversion
  Ks<-as.numeric(M%*%K)
  Vs<-as.numeric(M%*%V)
  
  E01<-(8*(Ks*Vs)^0.5)* 6E23 / 4184
  d1<-((Ks*a^2)/Vs)^(0.5)/a;
  c1<-(Ks*a^2/Is)^0.5
  m1<-E01/c1/6.011E-26
  taus1<-as.numeric(M%*%tau) #!as.numeric conversion
  gc1 = M[,3] + M[,4]
  
  Is<-as.numeric(M%*%I)#! numeric conversion
  Ks<-as.numeric(M%*%K)
  Vs<-as.numeric(M%*%V)
  
  
  E02<- 8*(Ks*Vs)^0.5  * 6E23 / 4184;
  d2<-((Ks*a^2)/Vs)^(0.5)/a;
  c2<-(Ks*a^2/Is)^0.5;
  m2<-E02/c2/6.011E-26;
  taus2<-as.numeric(M_comp%*%tau)
  gc2 = M_comp[,3] + M_comp[,4]
  
  dynchars_return<-list(E01=E01, d1=d1, c1=c1, m1=m1, taus1=taus1, gc1=gc1, E02=E02, d2=d2, c2=c2, m2=m2, taus2=taus2, gc2=gc2)
  
  return(dynchars_return)
  
}
load('/home/jane/Документы/Misha/Lab/spline_dataset_pro.Rdata')

e01_U00096.2 <- dynchars(unlist(strsplit(U00096.2, split = '')), interval_size = 50)$E01
x<- which.min(e01_U00096.2)
substrings_accap <-substr(U00096.2, start = x-200, stop = x+200)

int_seq_copy <- substrings_accap
int_seq <- substrings_accap

snps <- c()
for (i in 1:nchar(int_seq)) {
  nuc <- substr(int_seq, i, i)
  
  for (j in setdiff(c('a', 'c', 'g', 't'), nuc)) {
    substr(int_seq_copy, i, i) <- j
    snps <- c(snps, int_seq_copy)
    int_seq_copy <- int_seq
    
  }
}


snps_e0 <- c()
for (i in seq_along(snps)) {
  tmp <- unlist(strsplit(snps[i], split = ''))
  tmp1 <- dynchars(seq = tmp, interval_size = 50 )$E01
  snps_e0 <- cbind(snps_e0, tmp1)
}


matplot((snps_e0), type = 'l', lty = 1, col = 1)

#artificail for e0


library(dplyr)

library(tidyr)
library(ngram)
d <- data_frame(txt = substrings_accap)
d


sa <- unlist(strsplit(substrings_accap, split = ''))
sa1 <- paste(sa, collapse = ' ')
ng2 <- ngram(sa1, n=2)
ng3 <- ngram(sa1, n=3)
ng4 <- ngram(sa1, n=4)
ng5 <- ngram(sa1, n=5)
ng6 <- ngram(sa1, n=6)

print(ng6, output=  "full")


phr_ng6 <- get.phrasetable(ng6)
#artificial_seq <- paste0(unlist(strsplit(paste0(phr_ng6$ngrams[1:16], collapse = ''), split = ' ')), collapse = '')
artificial_seq <- paste0(unlist(strsplit(babble(ng6 , 201, seed =10), split = ' ')), collapse = '')
#artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)

#searching through generated seqs

artificial_e0 <- c()
artificial_seqs <- c()
for (i in 1:30000) {
  tmp0 <- paste0(unlist(strsplit(babble(ng6 , 201, seed =i+3), split = ' ')), collapse = '')
  artificial_seqs <- c(artificial_seqs, tmp0)
  
  tmp <- unlist(strsplit(tmp0, split = ''))
 tmp1 <- dynchars(seq = tmp, interval_size = 50 )$E01
 #artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)
 
 artificial_e0 <- cbind(artificial_e0, tmp1)
}

boxplot(apply(artificial_e0, 2, min))
boxplot(apply(artificial_e0, 2, mean))
which_min_min <- which.min(apply(artificial_e0, 2, min))
which_min_mean <- which.min(apply(artificial_e0, 2, mean))

plot(artificial_e0[,which_min_min], type ='l')
lines(artificial_e0[,which_min_mean], type ='l', col =2)


abline(h = mean(e01_U00096.2), lty = 3, col = 'tomato')
abline(h = min(e01_U00096.2), lty = 3, col = 'royalblue')

artificial_lowest_e0 <- artificial_seqs[which_min_mean]
GC(unlist(strsplit(artificial_lowest_e0, split = '')))

as.DNAbin(unlist(strsplit(artificial_lowest_e0, split = '')))

#artificial_lowest_e0 ep and bendability

artificial_lowest_e0_ep <- lseqspline1D(artificial_lowest_e0, bound = c(1,201), width = 1, ref = 100)$mpot
plot(artificial_lowest_e0_ep[50:680], type = 'l')
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col ='tomato')
abline(h = min(mpots_U00096.2$mpot), lty = 3, col ='royalblue')

bend_U00096.2 <- (bendability(U00096.2, bound = c(1, nchar(U00096.2))))
artificial_lowest_e0_bend <- bendability(artificial_lowest_e0, bound = c(1,201), width = 1)
plot(artificial_lowest_e0_bend, type = 'l')
abline(h = mean(bend_U00096.2, na.rm = T), lty = 3, col ='tomato')
abline(h = min(bend_U00096.2, na.rm = T), lty = 3, col ='royalblue')


#artifical lowest nbendability
library(dplyr)

library(tidyr)
library(ngram)

x <- which.min(bend_U00096.2)

substring_lowest_bend <- substr(U00096.2, start = x-150, stop = x+50 )


sa <- unlist(strsplit(substring_lowest_bend, split = ''))

dnabin_lowest_bend <- as.DNAbin(sa)
s

sa1 <- paste(sa, collapse = ' ')
ng2 <- ngram(sa1, n=2)
ng3 <- ngram(sa1, n=3)
ng4 <- ngram(sa1, n=4)
ng5 <- ngram(sa1, n=5)
ng6 <- ngram(sa1, n=6)

print(ng6, output=  "full")


phr_ng6 <- get.phrasetable(ng6)
#artificial_seq <- paste0(unlist(strsplit(paste0(phr_ng6$ngrams[1:16], collapse = ''), split = ' ')), collapse = '')
artificial_seq <- paste0(unlist(strsplit(babble(ng6 , 201, seed =10), split = ' ')), collapse = '')
#artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)

#searching through generated seqs

artificial_bend <- c()
artificial_seqs <- c()
for (i in 1:10000) {
  tmp0 <- paste0(unlist(strsplit(babble(ng6 , 201, seed =i+3), split = ' ')), collapse = '')
  artificial_seqs <- c(artificial_seqs, tmp0)
  
  #tmp <- unlist(strsplit(tmp0, split = ''))
  tmp1 <- bendability(s = tmp0, bound = range(seq_along(tmp)), width = 1)
  #artificial_mpots <- (lseqspline1D(artificial_seq, bound = c(50, nchar(artificial_seq)-50), ref = 50)$mpot)
  
  artificial_bend <- cbind(artificial_bend, tmp1)
}

boxplot(apply(artificial_bend, 2, function(x) {min(x, na.rm = T)}))
boxplot(apply(artificial_bend, 2, function(x) {mean(x, na.rm = T)}))
which_min_min <- which.min(apply(artificial_bend, 2, function(x) {min(x, na.rm = T)}))
which_min_mean <- which.min(apply(artificial_bend, 2, function(x) {mean(x, na.rm = T)}))

plot(artificial_bend[,which_min_min], type ='l')
lines(artificial_bend[,which_min_mean], type ='l', col =2)


abline(h = mean(bend_U00096.2, na.rm = T), lty = 3, col = 'tomato')
abline(h = min(bend_U00096.2, na.rm = T), lty = 3, col = 'royalblue')

seq_bend_min_min <- artificial_seqs[which_min_min]
seq_bend_min_mean <- artificial_seqs[which_min_mean]

as.DNAbin(unlist(strsplit(seq_bend_min_min, split = '')))
as.DNAbin(unlist(strsplit(seq_bend_min_mean, split = '')))

#artificial_lowest_e0 ep and bendability

artificial_lowest_e0_ep <- lseqspline1D(artificial_lowest_e0, bound = c(1,201), width = 1, ref = 100)$mpot
plot(artificial_lowest_e0_ep[50:680], type = 'l')
abline(h = mean(mpots_U00096.2$mpot), lty = 3, col ='tomato')
abline(h = min(mpots_U00096.2$mpot), lty = 3, col ='royalblue')

bend_U00096.2 <- (bendability(U00096.2, bound = c(1, nchar(U00096.2))))
artificial_lowest_e0_bend <- bendability(artificial_lowest_e0, bound = c(1,201), width = 1)

#mol_a 2017

appY_green <-  readLines('/home/jane/Документы/Misha/appY_green.txt')

nchar(appY_green)

appY_green_output<- lseqspline1D(appY_green, bound = c(50, nchar(appY_green)-50), width = 1, ref = 1)
appY_green_mpot <- appY_green_output$mpot 

appY_green_bend <- bendability(appY_green, bound = c(50, nchar(appY_green)-50), width = 1)

appY_green_e0 <- dynchars(unlist(strsplit(appY_green, split = '')), interval_size = 50)$E01

appY_green_gc <- dynchars(unlist(strsplit(appY_green, split = '')), interval_size = 50)$gc1
#plasmid regions
#rough positions

coeff <- length(appY_green_mpot)/nchar(appY_green)
bg1 <- nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTA')
egfp <- nchar('CTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCAT')
bg2_with_primers <- nchar('AGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTAATTTTATGATAGGTGCAAGATGGATTATGTTTGCTCCGTAGTTTTCATCTGTCAATCATTTGATTTAATTATAAACAGGAGAGTTATCTCGTTCAAAAAAAATTCATTGTTTATTGTAAGCGACAAAATTAGAAGGGAGTTACCAGTATGCCCCTCTAAACTAAGATCTAACTAATTAAGTAGCTAGAGTTCTATCAAGAGGTAGCT')
mCherry <- nchar('ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG')
bg3 <- nchar('TAACTCGAGCACCACCACCACCACCACTGAGATCCGGCTGCTAACAAAGCCCGAAAGGAAGCTGAGTTGGCTGCTGCCACCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCTGAAAGGAGGAACTATATCCGGATTGGCGAATGGGACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGTTTACAATTTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAATTAATTCTTA')
can_resist <- nchar("GAAAAACTCATCGAGCATCAAATGAAACTGCAATTTATTCATATCAGGATTATCAATACCATATTTTTGAAAAAGCCGTTTCTGTAATGAAGGAGAAAACTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGGTATCGGTCTGCGATTCCGACTCGTCCAACATCAATACAACCTATTAATTTCCCCTCGTCAAAAATAAGGTTATCAAGTGAGAAATCACCATGAGTGACGACTGAATCCGGTGAGAATGGCAAAAGTTTATGCATTTCTTTCCAGACTTGTTCAACAGGCCAGCCATTACGCTCGTCATCAAAATCACTCGCATCAACCAAACCGTTATTCATTCGTGATTGCGCCTGAGCGAGACGAAATACGCGATCGCTGTTAAAAGGACAATTACAAACAGGAATCGAATGCAACCGGCGCAGGAACACTGCCAGCGCATCAACAATATTTTCACCTGAATCAGGATATTCTTCTAATACCTGGAATGCTGTTTTCCCGGGGATCGCAGTGGTGAGTAACCATGCATCATCAGGAGTACGGATAAAATGCTTGATGGTCGGAAGAGGCATAAATTCCGTCAGCCAGTTTAGTCTGACCATCTCATCTGTAACATCATTGGCAACGCTACCTTTGCCATGTTTCAGAAACAACTCTGGCGCATCGGGCTTCCCATACAATCGATAGATTGTCGCACCTGATTGCCCGACATTATCGCGAGCCCATTTATACCCATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTAGAGCAAGACGTTTCCCGTTGAATATGGCTCAT")
bg4 <- nchar('AACACCCCTTGTATTACTGTTTATGTAAGCAGACAGTTTTATTGTTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGGTCATGGCTGCGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGGCAGCTGCGGTAAAGCTCATCAGCGTGGTCGTGAAGCGATTCACAGATGTCTGCCTGTTCATCCGCGTCCAGCTCGTTGAGTTTCTCCAGAAGCGTTAATGTCTGGCTTCTGATAAAGCGGGCCATGTTAAGGGCGGTTTTTTCCTGTTTGGTCACTGATGCCTCCGTGTAAGGGGGATTTCTGTTCATGGGGGTAATGATACCGATGAAACGAGAGAGGATGCTCACGATACGGGTTACTGATGATGAACATGCCCGGTTACTGGAACGTTGTGAGGGTAAACAACTGGCGGTATGGATGCGGCGGGACCAGAGAAAAATCACTCAGGGTCAATGCCAGCGCTTCGTTAATACAGATGTAGGTGTTCCACAGGGTAGCCAGCAGCATCCTGCGATGCAGATCCGGAACATAATGGTGCAGGGCGCTGACTTCCGCGTTTCCAGACTTTACGAAACACGGAAACCGAAGACCATTCATGTTGTTGCTCAGGTCGCAGACGTTTTGCAGCAGCAGTCGCTTCACGTTCGCTCGCGTATCGGTGATTCATTCTGCTAACCAGTAAGGCAACCCCGCCAGCCTAGCCGGGTCCTCAACGACAGGAGCACGATCATGCGCACCCGTGGGGCCGCCATGCCGGCGATAATGGCCTGCTTCTCGCCGAAACGTTTGGTGGCGGGACCAGTGACGAAGGCTTGAGCGAGGGCGTGCAAGATTCCGAATACCGCAAGCGACAGGCCGATCATCGTCGCGCTCCAGCGAAAGCGGTCCTCGCCGAAAATGACCCAGAGCGCTGCCGGCACCTGTCCTACGAGTTGCATGATAAAGAAGACAGTCATAAGTGCGGCGACGATAGTCATGCCCCGCGCCCACCGGAAGGAGCTGACTGGGTTGAAGGCTTCAAGGGCATCGGTCGAGATCCCGGTGCCTAATGAGTGAGCTAACTTACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCA')


##egfp_plot <- appY_green_output$x[which(appY_green_output$x%in%(bg1:(bg1+egfp)))]
##mCherry_plot <- appY_green_output$x[which(appY_green_output$x%in%((bg1+egfp+bg2_with_primers):(bg1+egfp+bg2_with_primers+mCherry)))]

appY_green_sidd_online <- (read.table('/home/jane/Документы/Misha/appy_green_sidd_online.txt', sep = '\t', skip = 16, header = T))$P.x.

starts_red <- c(  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATA'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCAC'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTA'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTAATTTTATGAT'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTAATTTTATGATAGGTGCAAGATGG'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTAATTTTATGATAGGTGCAAGATGGATTATGT'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTAATTTTATGATAGGTGCAAGATGGATTATGT'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTAATTTTATGATAGGTGCAAGATGGATTATGTTTGCTCCGTAGTTTTCATCTGTC'),
  nchar('GCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCAGCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATCTCAGTGGTGGTGGTGGTGGTGCTCGAGTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCTGAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGTGCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTTGTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTTGCTCACCATAGCTACCTCTTGATAGAACTCTAGCTACTTAATTAGTTAAGATCTTTCAGGTGCGTTGTAGTGAGTTTATGTTAATAAAAAGCATAGTAAGCGTTGAAAAATGTAACTTTGAAATAAGTTAGAATAAAAAACAACATACATATAATAATTTAATCTTAAATGAAATTTATTAAAATTTGCAAACTATAATTTTGTGTATAAAAATATAAATGCACATCATCCTGATTATGATTGTGTATTTAATTGGTTGTTATTTGACTACTATCAACTTGTTTTAATTTTATGATAGGTGCAAGATGGATTATGTTTGCTCCGTAGTTTTCATCTGTCAATCATTTGATTTAATTATAAACAGGAGAGTTATCTCGTTCAAAAAAAATTCATTGTTTATTGTAAGCGACAAAATTAGAAGGGAGTTACCAGTATGCC'))
                
svg('/home/jane/Документы/Misha/appy_green_physics.svg', width = 6, height = 6)
par(mfrow = c (3,2))
#sliding window size
wind <- 100
plot(appY_green_mpots, type = 'l', col = 'lightgrey', main = 'Electrostatic potential', xlab = 'Sequence (angstrom)', ylab = 'EP value')
lines(rollapply(appY_green_mpots, wind*coeff, mean), col = 1, lwd =2)
segments(x0 =coeff*(bg1+egfp), x1=coeff*(bg1+egfp+bg2_with_primers), y0 = -.068, y1 = -0.068, col ='darkgreen', lwd =3)
segments(x0 =coeff*(bg1+egfp+bg2_with_primers), x1=coeff*(bg1+egfp+bg2_with_primers+mCherry), y0 = -.068, y1 = -0.068, col ='tomato', lwd =3)
segments(x0 =coeff*(bg1+egfp+bg2_with_primers+mCherry+bg3), x1=coeff*(bg1+egfp+bg2_with_primers+mCherry+bg3+can_resist), y0 = -.068, y1 = -0.068, col ='royalblue', lwd =3)
abline(v = starts_red*coeff, col = 'darkred', lty = 3)

plot(appY_green_bend, type = 'l', col = 'lightgrey', main = 'Bendability', xlab = 'Sequence (nts)', ylab = 'Value')
lines(rollapply(appY_green_bend, wind, mean), col = 1, lwd =2)
segments(x0 =(bg1+egfp), x1=(bg1+egfp+bg2_with_primers), y0 = -1.1, y1 = -1.1, col ='darkgreen', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers), x1=(bg1+egfp+bg2_with_primers+mCherry), y0 = -1.1, y1 = -1.1, col ='tomato', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers+mCherry+bg3), x1=(bg1+egfp+bg2_with_primers+mCherry+bg3+can_resist), y0 = -1.1, y1 = -1.1, col ='royalblue', lwd =3)
abline(v = starts_red, col = 'darkred', lty = 3)


plot(appY_green_e0, type = 'l', col = 'lightgrey', main = 'Open state activation energy', xlab = 'Sequence (nts)', ylab = 'kcal/mol')
lines(rollapply(appY_green_e0, wind, mean), col = 1, lwd =2)
y=min(appY_green_e0)+3
segments(x0 =(bg1+egfp), x1=(bg1+egfp+bg2_with_primers), y0 = y, y1 = y, col ='darkgreen', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers), x1=(bg1+egfp+bg2_with_primers+mCherry), y0 = y, y1 = y, col ='tomato', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers+mCherry+bg3), x1=(bg1+egfp+bg2_with_primers+mCherry+bg3+can_resist), y0 = y, y1 = y, col ='royalblue', lwd =3)
abline(v = starts_red, col = 'darkred', lty = 3)


plot(appY_green_gc, type = 'l', col = 'lightgrey', main = 'GC-content', xlab = 'Sequence (nts)', ylab = 'GC-content')
lines(rollapply(appY_green_gc, wind, mean), col = 1, lwd =2)
y=min(appY_green_gc)+0.02
segments(x0 =(bg1+egfp), x1=(bg1+egfp+bg2_with_primers), y0 = y, y1 = y, col ='darkgreen', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers), x1=(bg1+egfp+bg2_with_primers+mCherry), y0 = y, y1 = y, col ='tomato', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers+mCherry+bg3), x1=(bg1+egfp+bg2_with_primers+mCherry+bg3+can_resist), y0 = y, y1 = y, col ='royalblue', lwd =3)
abline(v = starts_red, col = 'darkred', lty = 3)

plot(appY_green_sidd_online, type = 'l', col = 'lightgrey', main = 'SIDD', xlab = 'Sequence (nts)', ylab = 'Opening probability')
lines(rollapply(appY_green_sidd_online, wind, mean), col = 1, lwd =2)
y=min(appY_green_gc)+0.02
segments(x0 =(bg1+egfp), x1=(bg1+egfp+bg2_with_primers), y0 = y, y1 = y, col ='darkgreen', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers), x1=(bg1+egfp+bg2_with_primers+mCherry), y0 = y, y1 = y, col ='tomato', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers+mCherry+bg3), x1=(bg1+egfp+bg2_with_primers+mCherry+bg3+can_resist), y0 = y, y1 = y, col ='royalblue', lwd =3)
abline(v = starts_red, col = 'darkred', lty = 3)

plot(1, type="n", axes=F, xlab="", ylab="")
legend('center', legend = c('egfp', 'mCherry', 'canamycine resistence'), bty ='n', lty = 1, lwd = 2,col =  c('darkgreen', 'tomato', 'royal blue'), cex= 1.25)
dev.off()

#egfp region appears to be GC-poor 
(as.DNAbin(strsplit(substr(appY_green, start = bg1, stop = bg1+egfp), split = '')))
GC(unlist(strsplit(substr(appY_green, start = bg1, stop = bg1+egfp), split = '')))
GC(unlist(strsplit(appY_green, split = '')))

abline(h = mean(bend_U00096.2, na.rm = T), lty = 3, col ='tomato')
abline(h = min(bend_U00096.2, na.rm = T), lty = 3, col ='royalblue')

#palyndromic rate

palyndrate <- function(x, method = 'lv') {
  library(stringdist)
  x_string <- paste0(x, collapse = '')
  x_rev_string <- paste0(rev(x), collapse = '')
  tmp <- stringdist(x_string, x_rev_string, method = method)
  return(tmp)
}


par(mfrow = c(4,4),
    mar = rep(0.9,4))
for (wind in seq(10, 75, by = 5)) {
#x <- unlist(strsplit(x, split = ''))
appY_green_palyndrate25 <- rollapply(unlist(strsplit(appY_green, split = '')), wind, function(x) {palyndrate(x, method = 'hamming')})
#summary(palyndrate_appy_green)

plot(appY_green_palyndrate25, type = 'l', col = 'lightgrey', main = wind, xlab = 'Sequence (nts)', ylab = 'Palyndromic rate')
lines(rollapply(appY_green_palyndrate25, wind, mean), col = 1, lwd =2)
y=min(appY_green_palyndrate25)+0.02
segments(x0 =(bg1+egfp), x1=(bg1+egfp+bg2_with_primers), y0 = y, y1 = y, col ='darkgreen', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers), x1=(bg1+egfp+bg2_with_primers+mCherry), y0 = y, y1 = y, col ='tomato', lwd =3)
segments(x0 =(bg1+egfp+bg2_with_primers+mCherry+bg3), x1=(bg1+egfp+bg2_with_primers+mCherry+bg3+can_resist), y0 = y, y1 = y, col ='royalblue', lwd =3)
abline(v = starts_red, col = 'darkred', lty = 3)
}
           
which.forward <- which(lapply(dataset_pro, FUN = function(x){x$strand})=='forward')
which.experimental <- which(lapply(dataset_pro, FUN = function(x){x$evidence})=='experimental')
nos <- intersect(which.forward, which.experimental)
forward_experimental_seqs <- sapply(dataset_pro, FUN = function(x){x$seq})[nos]


forward_experimental_palyndrates <- sapply(forward_experimental_seqs, function(x){palyndrate(unlist(strsplit(x, split = '')), method = 'lv')})
summary(forward_experimental_palyndrates)

load('/home/jane/Загрузки/spline_dataset_gen.Rdata')
gen_seqs <- sapply(dataset_gen, FUN = function(x){x$seq})
gen_palyndrates <- sapply(gen_seqs, function(x){palyndrate(unlist(strsplit(x, split = '')), method = 'lv')})

par(mfrow = c(1,2))
boxplot(forward_experimental_palyndrates, ylim = c(180,250))
boxplot(gen_palyndrates, ylim = c(180,250))

forward_experimental_roll_palyndrate25 <- sapply(forward_experimental_seqs, FUN = function(x) {rollapply(unlist(strsplit(x, split = '')), 25, palyndrate)})
#matplot(forward_experimental_roll_palydrate25, type = 'l', col =1)
gen_roll_palyndrate25 <- sapply(gen_seqs, FUN = function(x) {rollapply(unlist(strsplit(x, split = '')), 25, palyndrate)})
boxplot(apply(forward_experimental_roll_palyndrate25, 2, min), ylim = c(0,15))
boxplot(sapply(gen_roll_palyndrate25,  min), ylim = c(0,15))

#palyndromes are asessed for sliding window 25 nts
par(mfrow = c(1,2))
boxplot(as.numeric(forward_experimental_roll_palyndrate25))
boxplot(as.numeric(gen_roll_palyndrate25))


#loose palyndromes and T7
library(zoo)
load('/home/jane/Документы/Misha/Lab/mol_a-16/t7/t7_rmd/NC_001604.1.rda')
promoters <- read.table('/home/jane/Документы/Misha/Lab/mol_a-16/t7/t7_rmd/promoters.txt', sep = '\t', header = T)
t7_seq <- as.character(NC_001604.1)

t7_seq_char <- unlist(strsplit(t7_seq, split = ''))

t7_seq_roll_palyndrate25 <- rollapply(t7_seq_char, 15, palyndrate)
                    
summary(t7_seq_roll_palyndrate25)  

plot(t7_seq_roll_palyndrate25[1:1000], type = 'l')
abline(v = promoters$TSS, lty = 3, col = 'red')


#palyndromic rate

palyndrate <- function(x, method = 'lv') {
  library(stringdist)
  x_string <- paste0(x, collapse = '')
  x_rev_string <- paste0(rev(x), collapse = '')
  tmp <- stringdist(x_string, x_rev_string, method = method)
  return(tmp)
}