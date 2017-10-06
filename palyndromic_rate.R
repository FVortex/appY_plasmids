
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

t7_seq_roll_palyndrate25 <- rollapply(t7_seq_char, 10, palyndrate)

summary(t7_seq_roll_palyndrate25)  

plot(t7_seq_roll_palyndrate25[1:1000], type = 'l')
abline(v = promoters$TSS, lty = 3, col = 'red')

X <- -100:100
plot(X, t7_seq_roll_palyndrate25[X+500],type = 'n', ylim = c(0,25))
lapply(promoters$TSS, FUN = function(x) {lines(X,t7_seq_roll_palyndrate25[(x-100):(x+100)])})

#lowset values
which(t7_seq_roll_palyndrate25==4)
which(t7_seq_roll_palyndrate25==6)

plot(t7_seq_roll_palyndrate25, type = 'l')
abline(v = promoters$TSS, lty = 3, col = 'red')
abline(v = c(which(t7_seq_roll_palyndrate25==0), which(t7_seq_roll_palyndrate25==0)), lty = 3, col = 'blue')

#comparing differnet string mertics
Methods <- c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")
dif_metrics <- c()
for (i in Methods){
tmp <- rollapply(t7_seq_char, 15, function(x) {palyndrate(x, method = i)})
dif_metrics <- cbind(dif_metrics, tmp)
}

colnames(dif_metrics) <- Methods
#assign(paste0('t7_seq_roll_palyndrate25_', i), value = tmp)
matplot(dif_metrics, type = 'l')

N <- 25
dif_metrics_mins <- apply(dif_metrics, 2, function(x){order(x)[1:N]})

dif_metrics_mins_no_nulls <- dif_metrics_mins[,-c(6:8, 10)]
#plot(t7_seq_roll_palyndrate25, type = 'n')

par(mfrow =c (3,3))
X <- 1:1600
for (i in 1:ncol(dif_metrics_mins_no_nulls)){
  plot(t7_seq_roll_palyndrate25[X], type = 'n')
  abline(v =  dif_metrics_mins_no_nulls[,i], col = i)
}

plot(t7_seq_roll_palyndrate25, type = 'l', col =' lightgrey')
abline(v = promoters$TSS, lty = 3, col = 'red')
