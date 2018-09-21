# SING_Canada_2018_genetics
Associated code for analysis of SING 2018 genetic data.    
More information about the Summer Internship for INdigenous Peoples in Genomics (SING) Canada program can be found here: http://indigenoussts.com/sing-canada/sing-canada-2018/

#### Requirements:    
`vcftools`    https://vcftools.github.io/index.html    
`adegenet`    https://cran.r-project.org/web/packages/adegenet/index.html    
`vcfR`    https://cran.r-project.org/web/packages/vcfR/index.html    
`hierfstat`     https://cran.r-project.org/web/packages/hierfstat/index.html    
`ape`     https://cran.r-project.org/web/packages/ape/index.html    
`jvarkit`    http://lindenb.github.io/jvarkit/MsaToVcf.html    

`java jdk 8.0`    http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
To use jvarkit, you will need java jdk 8.0. Check your java: 
Check java:
`/usr/bin/java -version`     

#### Inputs:
multiple sequence alignment file (msa)    

### 1. Prepare alignment file into a VCF
Put the .aln file into `02_raw_data`.        

Remove the individuals of the other species:    
`grep -vE 'LITTLENECK' 02_raw_data/butter__littleneck_clam_up.aln > 02_raw_data/butterclam.aln`

Move from the sequence alignment file into a vcf file (provide path to msa2vcf)    
`cat 02_raw_data/butterclam.aln | /Users/wayne/programs/jvarkit/dist/msa2vcf > 03_vcf/aln.vcf`

### 2. Filter the VCF
How many SNPs were identified before filtering?
`grep -vE '^#' 03_vcf/aln.vcf | wc -l`
111 sites are viewed in the vcf.   

And how many samples?    
`grep '|' 02_raw_data/butterclam.aln  | awk '{ print $1}' - | sort -n | uniq  | wc -l`    
85 samples are present in the vcf after removing the other species.   

With this dataset, there are 85 individual samples. If we want to make sure that at least 2 individuals have the SNP, then require 4 copies (because this is a diploid VCF with a haploid dataset, each individual is represented by 2 copies).  

Use VCFtools to filter the VCF and only keep sites that are present at least four times (2 individuals). For our data, this would equate to a minor allele frequency filter (MAF) of 4 / (85 * 2) = 0.0235 (or MAF > 0.02).   
`vcftools --vcf 03_vcf/aln.vcf --recode --maf 0.023 --out 03_vcf/aln_maf`     
Keeps 4 of 111 sites.   

(#todo: using the -m flag in msa2vcf allows haploid data, but will this be compatible to downstream?)

### 3. Generate stats on data (Fst)
Use the Rscript `01_scripts/vcf_to_hfstat.R` to perform some stats on the data as well as generating an Fst value. 
