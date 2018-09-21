# VCF to hierfstat
setwd("~/Documents/other_research_associates/breden_felix_2018/SING_Canada_2018_genetics")

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

#### 1. Import and set up genlight ####
# Import the recoded vcf
vcf <- read.vcfR("03_vcf/aln_maf.recode.vcf")
vcf

# Make genlight file (but note only loci with two alleles will be kept, not haplotypes. 1 locus is dropped)
my.data.gl <- vcfR2genlight(vcf)
my.data.gl

indNames(my.data.gl) # provides the individual names
pop(my.data.gl)  # provides the pop names (currently NULL)

# replace the pop with the relevant section of the individual names
pops <- gsub(pattern = ".*\\|", replacement = "", x =  indNames(my.data.gl))
# remove subsite within each location
pops <- gsub(pattern = "\\_.*", replacement = "", x =  pops) 
# replace null w/ new pops variable
pop(my.data.gl) <- pops
pop(my.data.gl)

# Drop the individual with an 'incomplete' designation
my.data.gl <- my.data.gl[pop(my.data.gl) != "INCOMPLETE"] # remove the incomplete record
pop(my.data.gl)
my.data.gl

#### 2. Run data exploration/multidimensional scaling stats ####
# MAF plot
myFreq <- glMean(my.data.gl)
hist(myFreq, proba=T, col="gold", xlab = "Allele frequencies"
     , main = "Distribution of (second) allele frequencies"
     , ylim = c(0,50)
)
text(x = 0.4, y = 7, labels = paste(nLoc(my.data), " loci", sep = "" ))
temp <- density(myFreq)
lines(temp$x, temp$y, lwd=3)

## Multidimensional scaling stats are not really working, probably because so few loci
# pca1 <- glPca(my.data.gl, nf = 2)
# scatter(x = pca1, posi = "bottomright", xax = 1, yax = 2)
# title("PC1 (x) vs PC2 (y)", adj = 1)

# Run dapc
## Also not really working, due to so few loci
# dapc1 <- dapc(my.data.gl, n.pca = 10, n.da = 1)
# scatter(dapc1, scree.da = F, bg = "white", legend = T
#         , txt.leg=rownames(dapc1$means)
#         , posi.leg = "topleft"
# )

# Compoplot
# compoplot(dapc1
#           #, lab="" # comment out if you want sample labels
#           , txt.leg = rownames(dapc1$means)
#           , posi = "topright"
#           #          , cex = 0.7
#           , cex.names = 0.6 
# )

#### 3. Neighbour joining tree ####
# Separate loci
y <- seploc(my.data.gl, n.block=10, parallel=F)

## Use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(y, function(e) dist(as.matrix(e)))
class(lD)
names(lD)
class(lD[[1]]) # lD is a list of distance matrices bw pairs of indiv

# Obtain general distance matrix by summing these
D <- Reduce("+", lD)

# Plot a neighbor-joining tree using ape's nj function on the distance matrix
pdf(file = "04_results/njt.pdf", width = 7, height = 5)
par(mar=c(4,3,5,3))
plot(nj(D), type="fan", cex=0.3)
dev.off()

#### 4. Set up for hierfstat ####
# Convert to genind format
# First by converting to matrix
my.data.mat <- as.matrix(my.data.gl)
my.data.mat[1:8,]

# Second, translate the number of minor allele to genind format
my.data.mat[my.data.mat == 0] <- "1/1" #homozygote reference
# my.data.mat[my.data.mat == 0] <- "1/2" #heterozygote are not possible here
my.data.mat[my.data.mat == 2] <- "2/2" #homozygote alternate
#my.data.mat[my.data.mat == "NA"] <- "NA" # NA
my.data.mat[1:8,]

# Convert matrix to genind
my.data.gid <- df2genind(my.data.mat, sep = "/", ploidy = 2) # convert df to genind (NOTE: ploidy = 2 here, but actually = 1)

# Generate population IDs for genind file
pop(my.data.gid) # currently null
indNames(my.data.gid) # The individual names actually contain pop IDs
pops <- gsub(pattern = ".*\\|", replacement = "", x =  indNames(my.data.gid)) 
pop(my.data.gid) <- gsub(pattern = "\\_.*", replacement = "", x =  pops) 
pop(my.data.gid)

# Barplot sample size per population
pdf(file = "04_results/sample_size_per_pop.pdf", width = 7, height = 5)
par(mar=c(7,4,4,4))
table(pop(my.data.gid))
barplot(table(pop(my.data.gid)), col=funky(17)
        #, las=3, las = 1
        , las=2
        , xlab=""
        , ylab="Sample size"
        , ylim = c(0,25))
abline(h = c(10,20,30), lty=2)
dev.off()

# How many alleles? (should all be 2)
table(nAll(my.data.gid))

# Change from genind to hierfstat
all.data.hf <- genind2hierfstat(my.data.gid)
rownames(all.data.hf) <- indNames(my.data.gid)

# # Remove the LITTLENECKCLAM from RUSSELL POP
# all.data.hf <- all.data.hf[-64,]
# rownames(all.data.hf)

# Pairwise Fst
pairwise.wc.fst <- pairwise.WCfst(all.data.hf)
pairwise.wc.fst
write.csv(pairwise.wc.fst, file = "04_results/all_data_wcfst.csv")

# Bootstrapping
# requires latest hierfstat (v0.04-29) otherwise get error
# library(devtools)
# install_github("jgx65/hierfstat")
# library("hierfstat")
options(scipen = 999999)
nboots <- 100000
boot.fst.all <- boot.ppfst(dat = all.data.hf, nboot = nboots, quant = c(0.025,0.975))
boot.fst.all

# Collect output
lower.limit <- t(boot.fst.all$ll)
upper.limit <- boot.fst.all$ul
upper.limit[is.na(upper.limit)] <- 0
lower.limit[is.na(lower.limit)] <- 0
boot.fst.all.output <- upper.limit + lower.limit
boot.fst.all.output

filename <- paste0("04_results/fst_nboot_", nboots, ".csv")
write.csv(x = boot.fst.all.output, file = filename)
