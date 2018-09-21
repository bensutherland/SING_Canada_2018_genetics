# VCF to hierfstat
setwd("/Users/wayne/Documents/my_space/05_school-awards-talks/03_conferences_and_talks/2018-SING_genomics_SFU")

# source("https://bioconductor.org/biocLite.R")
# biocLite("VariantAnnotation")

# install.packages("vcfR")
library('vcfR')
# install.packages("shiny")
library("shiny")
# install.packages("adegenet")
library("adegenet")
# install.packages("hierfstat")
library("hierfstat")
# install.packages("ape")
library("ape")

vcf <- read.vcfR("butter_clam_out.recode.vcf")
vcf

x <- vcfR2genlight(vcf)
x
indNames(x)


pop(x) 
# currently the population is null, replace with the ind name component
pops <- gsub(pattern = ".*\\|", replacement = "", x =  indNames(x))
pops <- gsub(pattern = "\\_.*", replacement = "", x =  pops) 
pop(x) <- pops

pop(x)

# x <- x[indNames(x) != "BUTTERCLAM_01|RUSSELL"]
# x <- x[indNames(x) != "BUTTERCLAM_02|RUSSELL"]
# x <- x[indNames(x) != "BUTTERCLAM_31|PORTLAND_SHELLBE"]

x

# Run PCA
pca1 <- glPca(x, nf = 3)
scatter(x = pca1, posi = "bottomright", xax = 1, yax = 2)
title("PC1 (x) vs PC2 (y)", adj = 1)

# Run dapc
dapc1 <- dapc(x, n.pca = 10, n.da = 1)
scatter(dapc1, scree.da = F, bg = "white", legend = T
        , txt.leg=rownames(dapc1$means)
        , posi.leg = "topleft"
)

compoplot(dapc1
          #, lab="" # comment out if you want sample labels
          , txt.leg = rownames(dapc1$means)
          , posi = "topright"
          #          , cex = 0.7
          , cex.names = 0.6 
)

myFreq <- glMean(x)
hist(myFreq, proba=T, col="gold", xlab = "Allele frequencies"
     , main = "Distribution of (second) allele frequencies"
     , ylim = c(0,50)
)
text(x = 0.4, y = 7, labels = paste(nLoc(my.data), " loci", sep = "" ))
temp <- density(myFreq)
lines(temp$x, temp$y, lwd=3)


y <- seploc(x, n.block=10, parallel=F)

## Use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(y, function(e) dist(as.matrix(e)))
class(lD)
names(lD)
class(lD[[1]]) # lD is a list of distance matrices bw pairs of indiv

# Obtain general distance matrix by summing these
D <- Reduce("+", lD)

# Plot a neighbor-joining tree using ape's nj function on the distance matrix
par(mar=c(4,4,4,4))
plot(nj(D), type="fan", cex=0.7)


# Convert to genind format
my.data <- x

my.data.mat <- as.matrix(my.data)
my.data.mat[1:5,1:5]

# Second, translate the number of minor allele to genind format
my.data.mat[my.data.mat == 0] <- "1/1" #homozygote reference
my.data.mat[my.data.mat == 0] <- "1/2" #heterozygote
my.data.mat[my.data.mat == 1] <- "2/2" #homozygote alternate
#my.data.mat[my.data.mat == "NA"] <- "NA" # NA
my.data.mat[1:5,1:5]

# Convert matrix to genind
my.data.gid <- df2genind(my.data.mat, sep = "/", ploidy = 2) # convert df to genind

# Generate populations for genind
pop(my.data.gid) # currently null
indNames(my.data.gid) # The individual names actually contain pop IDs
pops <- gsub(pattern = ".*\\|", replacement = "", x =  indNames(x)) 
pop(my.data.gid) <- gsub(pattern = "\\_.*", replacement = "", x =  pops) 
pop(my.data.gid)

# Show sample size per population
par(mar=c(7,4,4,4))
table(pop(my.data.gid))
barplot(table(pop(my.data.gid)), col=funky(17)
        #, las=3, las = 1
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,25))
abline(h = c(10,20,30), lty=2)


# How many alleles? (should all be 2)
table(nAll(my.data.gid))

# Change from genind to hierfstat
all.data.hf <- genind2hierfstat(my.data.gid)
rownames(all.data.hf) <- indNames(my.data.gid)

# Remove the LITTLENECKCLAM from RUSSELL POP
all.data.hf <- all.data.hf[-64,]
rownames(all.data.hf)

# Pairwise Fst
pairwise.wc.fst <- pairwise.WCfst(all.data.hf)
pairwise.wc.fst
#write.csv(pairwise.wc.fst, file = "11-other_stats/all_data_wcfst.csv")

# Bootstrapping
# requires latest hierfstat (v0.04-29) otherwise get error
# library(devtools)
# install_github("jgx65/hierfstat")
# library("hierfstat")
boot.fst.all <- boot.ppfst(dat = all.data.hf, nboot = 100000, quant = c(0.025,0.975))
boot.fst.all
