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
SC.window <- read.csv("SouthCenter.popstats.slidewind.csv", header =T)
names(SC.window)
length(SC.window$scaffold)
# Add the means FST for the window across populations:
SC.window2 = SC.window %>% mutate(Mean_fst = rowMeans(select(., starts_with("Fst"))
View(SC.window2)
#Calc the outlier threshold
threshold.meanfst <- quantile(SC.window2$Mean_fst, 0.975, na.rm = T)
#id outliers and add to the tibble data frame
SC.window2<- SC.window2 %>% mutate(outlier = if_else(SC.window2$Mean_fst > threshold.meanfst, "outlier", "background"))
#Plot the results
SC.window2 %>% group_by(outlier) %>% tally()

Man.sc_fst <- ggplot(SC.window2, aes(SC.window2$scaffold, SC.window2$Mean_fst, ymin = 0, ymax = 1.0, colour = SC.window2$outlier)) + geom_point()
Man.sc_fst <- Man.sc_fst + geom_abline(slope = 0, intercept = threshold.meanfst)
Man.sc_fst
######################################################################################################
#Manhattan Plots pi
SC.window2 = SC.window2 %>% mutate(Mean_Pi = rowMeans(select(., starts_with("pi")), na.rm = T))
threshold.meanpi <- quantile(SC.window2$Mean_Pi, 0.975, na.rm = T)
#id outliers and add to the tibble data frame
SC.window2<- SC.window2 %>% mutate(outlier = if_else(SC.window2$Mean_Pi > threshold.meanpi, "outlier", "background"))
#Plot the results
SC.window2 %>% group_by(outlier) %>% tally()

Man.sc_pi <- ggplot(SC.window2, aes(SC.window2$scaffold, SC.window2$Mean_Pi, ymin = 0, ymax = 1.0, colour = SC.window2$outlier)) + geom_point()
Man.sc_pi <- Man.sc_fst + geom_abline(slope = 0, intercept = threshold.meanpi)
Man.sc_pi
#############################################################################################
#########################Selection with PCAdapt and Outflank#################################
#############################################################################################
library("devtools")
library("OutFLANK")
library("pcadapt")
###########################################################################################
#Read in the data:
SC.adapt <- read.pcadapt("~/SC_RAD_analysis/VCF/SouthCenter.snps.vcf", type = "vcf")
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
outliers <- which(qval < alpha)
length(outliers)
#Outliers with BH procedure
padj <- p.adjust(Scpcadapt$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
#Bonferroni- Most Conservative
padj <- p.adjust(Scpcadapt$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
print(outliers)

#Plot loadings to see if anything looks off:

par(mfrow = c(1,1))
for (i in 1:2)
plot(Scpcadapt$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
########################################################################################
#Outflank here 
########################################################################################
#read in the data and set up the dataframe from the documentation:
SouthCenter <- read.vcfR("~/SC_RAD_analysis/VCF/SouthCenter.snps.vcf")
geno <- extractchromosome <- getCHROM(SouthCenter)
Geno.SC <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
Geno.SC[geno %in% c("0/0", "0|0")] <- 0
Geno.SC[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
Geno.SC[geno %in% c("1/1", "1|1")] <- 2
Geno.SC[is.na(Geno.SC)] <- 9
table(as.vector(Geno.SC))
#####################################################################
#Calculate Fst
poplist.names <- c("12-16", "12-16",	"12-16", "12-16", "12-16", "12-16",	"12-16", "12-16",	"12-16",	"12-16",	"20-24",	"20-24", "20-24", "20-24",	"20-24",	"20-24", "20-24",	"20-24",	"20-24",	"20-24",	"4-8",	"4-8",	"4-8",	"4-8",	"4-8", "4-8",	"4-8",	"4-8",	"4-8",	"4-8",	"52-56",	"52-56", "60-64",	"lake", "lake",	"lake",	"lake", "lake", "lake",	"lake",	"lake",	"lake",	"lake")
colnames(Geno.SC) <- c("12-16", "12-16",	"12-16", "12-16", "12-16", "12-16",	"12-16", "12-16",	"12-16",	"12-16",	"20-24",	"20-24", "20-24", "20-24",	"20-24",	"20-24", "20-24",	"20-24",	"20-24",	"20-24",	"4-8",	"4-8",	"4-8",	"4-8",	"4-8", "4-8",	"4-8",	"4-8",	"4-8",	"4-8",	"52-56",	"52-56", "60-64",	"lake", "lake",	"lake",	"lake", "lake", "lake",	"lake",	"lake",	"lake",	"lake")
SC_fst <- MakeDiploidFSTMat(Geno.SC, locusNames = position, popNames = poplist.names)
SC_fst <- MakeDiploidFSTMat(t(Geno.SC), locusNames = position, popNames = colnames(Geno.SC))
head(SC_fst)
#QC plots for for sanity check:
plot(SC_fst$He, SC_fst$FST)
plot(Sc_fst$FST, SC_fst$FSTNoCorr)
abline(0,1)
plot(SC_fst$FST, SC_fst$FSTNoCorr)
abline(0,1)

#To do: 
# insert and LD pruned dataset for null distro 
# do selection tests for Fst ouliers
# other popgen simulations
###################################

