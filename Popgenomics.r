#######################################################################
############ Population Genomics of South Center Lake##################
#######################################################################
#set the working directory:
setwd("~/SC_RAD_analysis/Pop_GenomicsR")
#######################################################################
#Read in data: 
library(vcfR)
#Read in the vcf Data:
SouthCenter <- read.vcfR("~/SouthCenter.snps.vcf.gz")
#Convert to file formats for Adegenet and Popr
#Genlight formant
SouthCenter.gl <- vcfR2genlight(SouthCenter)
#GenID format
SouthCenter.gid <- vcfR2genind(SouthCenter)
#Add the population information and plody levels: 
pop(SouthCenter.gl) <- as.factor(c("12-16", "12-16",	"12-16", "12-16", "12-16", "12-16",	"12-16", "12-16",	"12-16",	"12-16",	"20-24",	"20-24", "20-24", "20-24",	"20-24",	"20-24", "20-24",	"20-24",	"20-24",	"20-24",	"4-8",	"4-8",	"4-8",	"4-8",	"4-8", "4-8",	"4-8",	"4-8",	"4-8",	"4-8",	"52-56",	"52-56", "60-64",	"lake", "lake",	"lake",	"lake", "lake", "lake",	"lake",	"lake",	"lake",	"lake"))
popNames(SouthCenter.gl)
ploidy(SouthCenter.gl) <- 2
#Now same for the GenID:
pop(SouthCenter.gid) <- as.factor(c("12-16", "12-16",	"12-16", "12-16", "12-16", "12-16",	"12-16", "12-16",	"12-16",	"12-16",	"20-24",	"20-24", "20-24", "20-24",	"20-24",	"20-24", "20-24",	"20-24",	"20-24",	"20-24",	"4-8",	"4-8",	"4-8",	"4-8",	"4-8", "4-8",	"4-8",	"4-8",	"4-8",	"4-8",	"52-56",	"52-56", "60-64",	"lake", "lake",	"lake",	"lake", "lake", "lake",	"lake",	"lake",	"lake",	"lake"))
popNames(SouthCenter.gid)
ploidy(SouthCenter.gid) <- 2
#Sanity check
SouthCenter.gl
SouthCenter.gid
##################################################################
############ Visalizing Population Structure #####################
##################################################################

tree.sc <- aboot(SouthCenter.gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
cols <- brewer.pal(n = nPop(SouthCenter.gl), name = "Dark2")
library(RColorBrewer)
cols <- brewer.pal(n = nPop(SouthCenter.gl), name = "Dark2")
plot.phylo(tree.sc, cex = 0.5, font = 2, adj = 0, tip.color =  cols[pop(SouthCenter.gl)])
nodelabels(tree.sc$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.5,font = 3, xpd = TRUE)
#legend('topleft', legend = c( "12-16", "20-24", "4-8", "52-56", "60-64", "lake"), fill = cols, border = FALSE, bty = "n", cex = 0.8)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
###############################################################################
#Principal Components Analysis:
SC.pca <- glPca(SouthCenter.gl)
#scatter plotting the PCA
SC.pca.scores <- as.data.frame(SC.pca$scores)
SC.pca.scores$pop <- pop(SouthCenter.gl)
set.seed(6)
scplot <- ggplot(SC.pca.scores, aes(x= PC1, y= PC2, colour= pop))
scplot <- scplot + geom_point(size = 2)
scplot <- scplot + stat_ellipse(level = 0.95, size = 1)
scplot <- scplot + geom_hline(yintercept = 0)
scplot <- scplot + geom_vline(xintercept = 0)
scplot <- scplot + theme_bw()
scplot
#Done
#########################################################################################################
#Discriminant Analysis of PCA
#Use cross validation to get correct number of PCs to retain: 
set.seed(999)
SCx <- xvalDapc(tab(SouthCenter.gl, NA.method = "mean"), pop(SouthCenter.gl))
#Converge to smaller number of axes:
set.seed(999)
system.time(SCxval <- xvalDapc(tab(SouthCenter.gl, NA.method + "mean"), pop(SouthCenter.gl), n.pca = 5:15, n.rep = 1000, parallel = "multicore", ncpus = 4L))

#Select the lowest core based analysis:
SCxval[-1]

#########################################################################################
#Run the DAPC and pot the scatter:
SC.dapc <- dapc(SouthCenter.gl, var.contrib = T, scale = F, n.pca = 11, n.da = 5)
scatter(SC.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
posi.pca = "topright", cleg = 0.75)
compoplot(SC.dapc, col = cols, posi = 'top')
dapc.results <- as.data.frame(SC.dapc$posterior)
dapc.results$pop <- pop(SouthCenter.gl)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop", "Sample", "Assigned_Pop", "Posterior_Membership_Prob")
bp <- ggplot(dapc.results, aes(x=Sample, y=Posterior_Membership_Prob, fill=Assigned_Pop))
bp <- bp + geom_bar(stat='identity')
bp <- bp + scale_fill_manual(values = cols)
bp <- bp + facet_grid(~Original_Pop, scales = "free")
bp <- bp + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
bp
#############################################################################
#Structure like Bar-plots of Postior Probs of assignment
library(reshape2)
compoplot(SC.dapc, col = cols, posi = 'top')
dapc.results <- as.data.frame(SC.dapc$posterior)
dapc.results$pop <- pop(SouthCenter.gl)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop", "Sample", "Assigned_Pop", "Posterior_Membership_Prob")
bp <- ggplot(dapc.results, aes(x=Sample, y=Posterior_Membership_Prob, fill=Assigned_Pop))
bp <- bp + geom_bar(stat='identity')
bp <- bp + scale_fill_manual(values = cols)
bp <- bp + facet_grid(~Original_Pop, scales = "free")
bp <- bp + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
bp
#################################################################################
####################Tests for selection##########################################
#################################################################################
#visualizing fst in a simple outlier test Fst calculated with the Vcftools on the command line:
library(tidyverse)
fst <- read_tsv("Allpops.weir.fst", col_name = T)
names(fst)
head(fst)
#take a look at the data
ggplot(fst, aes(POS, WEIR_AND_COCKERHAM_FST)) + geom_point()
quantile(fst$WEIR_AND_COCKERHAM_FST, c(0.975, 0.995), na.rm = T)
#Identify the values
threshold <- quantile(fst$WEIR_AND_COCKERHAM_FST, 0.975, na.rm = T)
#id outliers and add to the tibble data frame
fst <- fst %>% mutate(outlier = if_else(fst$WEIR_AND_COCKERHAM_FST > threshold, "outlier", "background"))
#Plot the results
fst %>% group_by(outlier) %>% tally()
ggplot(fst, aes(POS, WEIR_AND_COCKERHAM_FST, colour = outlier)) + geom_point()

#########################################################################################################
#Slidinw window analysis of Fst/ Pi
#Generated from Pyton script on command line
SC.window <- read.csv("~/SC_RAD_analysis/Genome_Scan/SouthCenter.window.csv", header =T)
names(SC.window)
length(SC.window$scaffold)
# Add the means FST for the window across populations: # remove the 60-64 cm pop bc of sample size
SC.window2 = SC.window %>% mutate(Mean_fst = rowMeans(select(., starts_with("Fst"), contains(c("52.56", "lake", "4.8", "12.16", "20.24")))))

#Calc the outlier threshold
threshold.meanfst <- quantile(SC.window2$Mean_fst, 0.975, na.rm = T)
#id outliers and add to the tibble data frame
SC.window2<- SC.window2 %>% mutate(outlier_fst = if_else(SC.window2$Mean_fst > threshold.meanfst, T, F))
#Plot the results
SC.window2 %>% group_by(outlier_fst) %>% tally()
grep("T", SC.window2$outlier_fst)
Man.sc_fst <- ggplot(SC.window2, aes(SC.window2$scaffold, SC.window2$Mean_fst, ymin = 0, ymax = 1.0, colour = SC.window2$outlier_fst)) + geom_point()
Man.sc_fst <- Man.sc_fst + geom_abline(slope = 0, intercept = threshold.meanfst)
Man.sc_fst
######################################################################################################
#Manhattan Plots pi
SC.window2 = SC.window2 %>% mutate(Mean_Pi = rowMeans(select(., starts_with("pi")), na.rm = T))
threshold.meanpi <- quantile(SC.window2$Mean_Pi, 0.975, na.rm = T)
#id outliers and add to the tibble data frame
SC.window2<- SC.window2 %>% mutate(outlier_pi = if_else(SC.window2$Mean_Pi > threshold.meanpi, "outlier", "background"))
#Plot the results
SC.window2 %>% group_by(outlier_pi) %>% tally()

Man.sc_pi <- ggplot(SC.window2, aes(SC.window2$scaffold, SC.window2$Mean_Pi, ymin = 0, ymax = 1.0, colour = SC.window2$outlier_pi)) + geom_point()
Man.sc_pi <- Man.sc_pi + geom_abline(slope = 0, intercept = threshold.meanpi)
Man.sc_pi
#############################################################################################
#########################Selection with PCAdapt and Outflank#################################
#############################################################################################
library("devtools")
library("OutFLANK")
library("pcadapt")
library("qvalue")
library("vcfR")
###########################################################################################
#Read in the data:
SC.adapt <- read.pcadapt("~/VCF_files/ALL_SNPs/SouthCenter.snps.vcf.gz", type = "vcf")
#Visualize:
Scpcadapt <- pcadapt(input = SC.adapt, K = 20)
plot(Scpcadapt, option = "screeplot")
#Add some Pop info
poplist.names <- c("12-16", "12-16",	"12-16", "12-16", "12-16", "12-16",	"12-16", "12-16",	"12-16",	"12-16",	"20-24",	"20-24", "20-24", "20-24",	"20-24",	"20-24", "20-24",	"20-24",	"20-24",	"20-24",	"4-8",	"4-8",	"4-8",	"4-8",	"4-8", "4-8",	"4-8",	"4-8",	"4-8",	"4-8",	"52-56",	"52-56", "60-64",	"lake", "lake",	"lake",	"lake", "lake", "lake",	"lake",	"lake",	"lake",	"lake")
plot(Scpcadapt, option = "scores", pop = poplist.names)
# Looks like K =2 makes the most sense here.
Scpcadapt <- pcadapt(input = SC.adapt, K = 2)
#inspect and plot the analysis:
summary(Scpcadapt)
plot(Scpcadapt, option ="manhattan")
plot(Scpcadapt, option = "qqplot")
hist(Scpcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(Scpcadapt, option = "stat.distribution")

#Correct for false discovery using FDR
qval <- qvalue(Scpcadapt$pvalues)$qvalues
alpha <- 0.1
outliers_FDR <- which(qval < alpha)
length(outliers_FDR)
print(outliers_FDR)
#Outliers with BH procedure
padj <- p.adjust(Scpcadapt$pvalues,method="BH")
alpha <- 0.1
outliers_BH <- which(padj < alpha)
length(outliers_BH)
print(outliers_BH)
#Bonferroni- Most Conservative
padj <- p.adjust(Scpcadapt$pvalues,method="bonferroni")
alpha <- 0.1
outliers_BF <- which(padj < alpha)
length(outliers_BF)
print(outliers_BF)

#Plot loadings to see if anything looks off:

par(mfrow = c(1,1))
for (i in 1:2)
plot(Scpcadapt$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
########################################################################################
#########################################Outflank here 
########################################################################################
#read in the data and set up the dataframe from the documentation: All snps
SouthCenter <- read.vcfR("/home/matt/VCF_files/ALL_SNPs/SouthCenter.all_rm6064.vcf.gz")
geno.all <- extract.gt(SouthCenter)
position.all <- getPOS(SouthCenter)
chromosome.all <- getCHROM(SouthCenter)
Geno.SCall <- matrix(NA, nrow = nrow(geno.all), ncol = ncol(geno.all))
Geno.SCall[geno.all %in% c("0/0", "0|0")] <- 0
Geno.SCall[geno.all %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
Geno.SCall[geno.all %in% c("1/1", "1|1")] <- 2
Geno.SCall[is.na(Geno.SCall)] <- 9
table(as.vector(Geno.SCall))
#####################################################################
#VCF that has been LD pruned by writing 1 snp per tag 
SouthCenter.LD <- read.vcfR("/home/matt/VCF_files/Single_SNP/SouthCenter.single_rm6064.vcf.gz")
geno.LD <- extract.gt(SouthCenter.LD)
position.LD <- getPOS(SouthCenter.LD)
chromosome.LD <- getCHROM(SouthCenter.LD)
Geno.SCLD <- matrix(NA, nrow = nrow(geno.LD), ncol = ncol(geno.LD))
Geno.SCLD[geno.LD %in% c("0/0", "0|0")] <- 0
Geno.SCLD[geno.LD %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
Geno.SCLD[geno.LD %in% c("1/1", "1|1")] <- 2
Geno.SCLD[is.na(Geno.SCLD)] <- 9
table(as.vector(Geno.SCLD))

#####################################################################
#Calculate Fst for All SNPs

colnames(Geno.SCall) <- c("12-16", "12-16",	"12-16", "12-16", "12-16", "12-16",	"12-16", "12-16",	"12-16",	"12-16",	"20-24",	"20-24", "20-24", "20-24",	"20-24",	"20-24", "20-24",	"20-24",	"20-24",	"20-24",	"4-8",	"4-8",	"4-8",	"4-8",	"4-8", "4-8",	"4-8",	"4-8",	"4-8",	"4-8",	"52-56",	"52-56",	"lake", "lake",	"lake",	"lake", "lake", "lake",	"lake",	"lake",	"lake",	"lake")

colnames(Geno.SCall)
SC_fst <- MakeDiploidFSTMat(t(Geno.SCall), locusNames = position.all, popNames = colnames(Geno.SCall))
head(SC_fst)
#QC plots for for sanity check:
plot(SC_fst$He, SC_fst$FST)
plot(SC_fst$FST, SC_fst$FSTNoCorr)
abline(0,1)
#Fst for LD SNPs
colnames(Geno.SCLD) <- c("12-16", "12-16",	"12-16", "12-16", "12-16", "12-16",	"12-16", "12-16",	"12-16",	"12-16",	"20-24",	"20-24", "20-24", "20-24",	"20-24",	"20-24", "20-24",	"20-24",	"20-24",	"20-24",	"4-8",	"4-8",	"4-8",	"4-8",	"4-8", "4-8",	"4-8",	"4-8",	"4-8",	"4-8",	"52-56",	"52-56",	"lake", "lake",	"lake",	"lake", "lake", "lake",	"lake",	"lake",	"lake",	"lake")

colnames(Geno.SCLD)
SC_fst_LD <- MakeDiploidFSTMat(t(Geno.SCLD), locusNames = position.LD, popNames = colnames(Geno.SCLD))
head(SC_fst_LD)
plot(SC_fst_LD$He, SC_fst_LD$FST)
plot(SC_fst_LD$FST, SC_fst_LD$FSTNoCorr)
abline(0,1)
#######################################################################
##Null Distribution with LD pruned SNPs

Pruned_fst<- OutFLANK(SC_fst_LD, NumberOfSamples = 42, qthreshold = 0.05, Hmin = 0.1)
str(Pruned_fst)
head(Pruned_fst$results)

OutFLANKResultsPlotter(Pruned_fst, withOutliers = T, NoCorr = T, Hmin = 0.1, binwidth = 0.001, Zoom = F, titletext = NULL)

OutFLANKResultsPlotter(Pruned_fst, withOutliers = T, NoCorr = T, Hmin = 0.1, binwidth = 0.001, Zoom = T, RightZoomFraction = 0.15, titletext = NULL)

hist(Pruned_fst$results$pvaluesRightTail)


All_loci <- pOutlierFinderChiSqNoCorr(SC_fst, Fstbar = Pruned_fst$FSTNoCorrbar, dfInferred = Pruned_fst$dfInferred, qthreshold = 0.05, Hmin = 0.1)

head(All_loci)

All_My_Outliers <- All_loci$OutlierFlag == T

plot(All_loci$He, All_loci$FST, pch = 19, col = rgb(0,0,0,0.1))
points(All_loci$He[All_My_Outliers], All_loci$FST[All_My_Outliers], col = "blue")

hist(All_loci$pvaluesRightTail)

plot(All_loci$LocusName[All_loci$He>0.1], All_loci$FST[All_loci$He>0.1], 
     xlab = "Position", ylab = "FST", col = rgb(0,0,0,0.2))
points(All_loci$LocusName[All_My_Outliers], All_loci$FST[All_My_Outliers], col = "red", pch = 20)
#To do: 

# other popgen simulations
###################################

