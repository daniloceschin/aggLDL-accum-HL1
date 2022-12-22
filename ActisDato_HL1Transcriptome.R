rm(list=ls())


# ipak function for intalling and calling several packages at once
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# package's list for using

#packages <- c("DESeq2", "ggplot2", "reshape2", "biomaRt", "RColorBrewer",
#              "factoextra", "gplots", "openxlsx", "edgeR", "corrplot", "limma", 
#              "org.Mm.eg.db", "GSEABase", "EnhancedVolcano", "heatmaply", "rWikiPathways",
#              "heatmap3")

packages <- c("openxlsx","biomaRt","edgeR", "ggplot2", "ggfortify", "corrplot",  "RColorBrewer", 
              "EnhancedVolcano")
# "reshape2", "biomaRt", "RColorBrewer",
# "factoextra", "gplots", "openxlsx", "edgeR", "corrplot", "limma", 
# "org.Mm.eg.db", "GSEABase", "EnhancedVolcano", "heatmaply", "rWikiPathways",
# "heatmap3")

# installation and/or calling packages
ipak(packages)

# set working directory
setwd("~/Documents/OneDrive/01.CIC/01.Proyectos/01.Investigadores/VickyAD/Transcriptome_HL1/Blind_Analysis/ScriptsPapiro")

DT1 <- read.table(file = "DT1_pub.csv", header = T, sep = ",", quote = "\"")
rownames(DT1) <- DT1$Geneid

##Blind analysis in Salmona

###----------------------------------------------###
### Building the experimental design information ###
###----------------------------------------------###

### count matrix
mycounts <- DT1[,c(9:24)] # count matrix
myannot <- DT1[,c(-2:-5,-9:-24)] # annotations matrix

## rearrange table following heatmap analysis
## 1. delete B3 sample
## 2. B1 = sampleA = SA
## 3. C1 + A1 = sampleB = SB
## 4. D1 + C2 = sampleC = SC
## 5. B2 + A2 = sampleD = SD

## Delete B3
#mycounts <- mycounts[,c(-17:-20)]

### Factors
### read metadata
### first time to build the metadata file
##write.csv(colnames(mycounts), file = "MetaHL1.csv")

met <- read.csv("MetaHL1_pub.csv")
met2 <- as.data.frame(colnames(mycounts))
colnames(met2)[1] <- "Seq_id"
met3 <- merge(met2, met, by="Seq_id", all.x = TRUE, sort=FALSE)
met <- met3
myfactors <- met
rm("met2", "met3")

myfactors$SelectCond <- c("HIL", "HIL", "HIL", "HIL",
                          "HI", "HI", "HI", "HI",
                          "HL", "HL", "HL", "HL",
                          "HM", "HM", "HM", "HM")

#### Differential Expression gene

## DGEList preparation and annotation

x <- DGEList(counts=mycounts)
colnames(x)

## rename cols using the blind annotations
mynameid <- c("HIL_1", "HIL_2", "HIL_3", "HIL_4",
              "HI_1", "HI_2", "HI_3", "HI_4",
              "HL_1", "HL_2", "HL_3", "HL_4",
              "HM_1", "HM_2", "HM_3", "HM_4")


colnames(x)<- mynameid
colnames(x)

##---------------------------------------------------
## Building the experimental design information  ----
##---------------------------------------------------
## blind
##group <- as.factor(as.character(myfactors[,"TreatBlind"]))

## final condition selected
group <- as.factor(as.character(myfactors[,"SelectCond"]))

x$samples$group <- group

lane <- met[,6]
x$samples$lane <- lane

id <- myfactors[, "SelectCond"]
x$samples$id <- id

x$samples

geneid <- rownames(x)
head(x$counts)


##-------------------------------------
## filtering low expressed genes
##----------------------------------

keep <- filterByExpr(x)
x <- x[keep, , keep.lib.sizes=FALSE]

### Add annotations

xannot <- as.data.frame(rownames(x))
colnames(xannot)[1] <- "Geneid"
xannot <- merge(xannot, myannot, by="Geneid", all.x = TRUE, sort = FALSE)
xannot <- xannot[!duplicated(xannot$Geneid),]

x$genes <- xannot

###-----------------------------------
### Normalization
### ---------------------------

x <- calcNormFactors(x, method = "TMM")

###-------------------
### building the design
###--------------------

#design <- model.matrix(~group)

design <- model.matrix(~0+group, data=x$samples)
colnames(design) <- levels(x$samples$group)
design

###----------------------------------------
### estimating dispersion
###----------------------------------------

x <- estimateDisp(x, design)


##-----------------------------------
### Testing for DE genes
###-----------------------------------

fit <- glmQLFit(x, design)

x$samples

### contrast with conditions

matrix.contr <- makeContrasts(HL.vs.HM = HL - HM,
                              HI.vs.HM = HI - HM,
                              HIL.vs.HM = HIL - HM,
                              levels = design)

#### comparison with conditions

qlf.HL.vs.HM <- glmQLFTest(fit, contrast = matrix.contr[,"HL.vs.HM"])
qlf.HI.vs.HM <- glmQLFTest(fit, contrast = matrix.contr[,"HI.vs.HM"])
qlf.HIL.vs.HM <- glmQLFTest(fit, contrast = matrix.contr[,"HIL.vs.HM"])

summary(decideTests(qlf.HL.vs.HM))
summary(decideTests(qlf.HI.vs.HM))
summary(decideTests(qlf.HIL.vs.HM))

#kable(summary(decideTests(qlf.HL.vs.HM))[c(3,1,2), 1:4], caption = "Summary of Tests", align = "c")

###----------------------------------------
### CPMs
###--------------------------------------

x.cpm <- cpm(x)
x.lcpm <- cpm(x, log = TRUE)

x$cpms <- x.cpm
x$lcpms <- x.lcpm

## table sfor conditions
HL.vs.HM <- topTags(qlf.HL.vs.HM, n = Inf)
HI.vs.HM <- topTags(qlf.HI.vs.HM, n = Inf)
HIL.vs.HM <- topTags(qlf.HIL.vs.HM, n = Inf)

### log cpms Average by conditions
HL.avg <- as.data.frame(rowMeans(x$lcpms[, c(9:12)]))
HIL.avg <- as.data.frame(rowMeans(x$lcpms[, c(1:4)]))
HM.avg <- as.data.frame(rowMeans(x$lcpms[, c(13:16)]))
HI.avg <- as.data.frame(rowMeans(x$lcpms[, c(5:8)]))

## Tables with cpms average by selected condition

HL.vs.HM$table <- cbind(HL.vs.HM$table[,1:10], 
                        Log2.CPM.HL = HL.avg[HL.vs.HM$table$Geneid,], 
                        Log2.CPM.HM = HM.avg[HL.vs.HM$table$Geneid,],
                        HL.vs.HM$table[,c(11, 14:15)])

HI.vs.HM$table <- cbind(HI.vs.HM$table[,1:10], 
                        Log2.CPM.HI = HI.avg[HI.vs.HM$table$Geneid,], 
                        Log2.CPM.HM = HM.avg[HI.vs.HM$table$Geneid,],
                        HI.vs.HM$table[,c(11, 14:15)])

HIL.vs.HM$table <- cbind(HIL.vs.HM$table[,1:10], 
                         Log2.CPM.HIL = HIL.avg[HIL.vs.HM$table$Geneid,], 
                         Log2.CPM.HM = HM.avg[HIL.vs.HM$table$Geneid,],
                         HIL.vs.HM$table[,c(11, 14:15)])


#### Figures Paper

### Data

mat.fig <- as.data.frame(x$counts[, c(5:16)]) ## matrix to generate papres' figures


### PCA

### PCA for paper
col.4paper <- as.factor(c("HI", "HI", "HI", "HI",
                          "HL", "HL", "HL", "HL",
                          "HM", "HM", "HM", "HM"))

lab.paper <- c("Insulin", "aggLDL", "Control")
col.4paper <- c(rep("slateblue3", 4), rep("firebrick2", 4), rep("olivedrab4", 4))


mat.fig$HI_1[ mat.fig$HI_1 == 0] <- 1 ## avoid variance = 0 for PCA
mat.fig$HL_1[ mat.fig$HL_1 == 0] <- 1
mat.fig$HM_1[ mat.fig$HM_1 == 0] <- 1


pca.fig <- prcomp(t(as.matrix(mat.fig)), scale. = TRUE)

# autoplot(pca.fig, data = pca.fig, colour = col.4paper, alpha = 1/5, variance_percentage = F, size =5, label =T) +
#   coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) +
#   theme_bw()

autoplot(pca.fig, data = pca.fig, colour = col.4paper, alpha = 1/5, variance_percentage = F, size =5, title = "PCA") +
  annotate("text", x= -0.26, y = 0.18, label= "Insulin", col= "slateblue3") +
  annotate("text", x= 0.38, y = -0.02, label= "aggLDL", col= "firebrick2") +
  annotate("text", x= -0.13, y = -0.27, label= "Control", col= "olivedrab4") +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) +
  theme_bw()


## PCA Complementary

mat.fig.comp <- as.data.frame(x$counts)
col.4comp <- c(rep("coral2", 4), rep("slateblue3", 4), rep("firebrick2", 4), rep("olivedrab4", 4))

lab.comp <- c("aggLDL + Insulin", "Insulin", "aggLDL", "Control")


mat.fig.comp$HIL_1[ mat.fig.comp$HIL_1 == 0] <- 1 ## avoid variance = 0 for PCA
mat.fig.comp$HI_1[ mat.fig.comp$HI_1 == 0] <- 1 
mat.fig.comp$HL_1[ mat.fig.comp$HL_1 == 0] <- 1
mat.fig.comp$HM_1[ mat.fig.comp$HM_1 == 0] <- 1

lab.pos <- c(-1, 1, 2, -2)
pca.fig.comp <- prcomp(t(as.matrix(mat.fig.comp)), scale. = TRUE)

# autoplot(pca.fig.comp, data = pca.fig.comp, colour = col.4comp, alpha = 1/5, variance_percentage = F, 
#          size =5, label =T, label.size = 3, label.hjust = -1) +
#   coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) +
#   theme_bw()

autoplot(pca.fig.comp, data = pca.fig.comp, colour = col.4comp, alpha = 1/5, variance_percentage = F, size =5) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) +
  annotate("text", x= -0.29, y = 0.18, label= "aggLDL + Insulin", col= "coral2") + 
  annotate("text", x= 0.31, y = -0.15, label= "Insulin", col= "slateblue3") +
  annotate("text", x= -0.18, y = -0.27, label= "aggLDL", col= "firebrick2") +
  annotate("text", x= 0.15, y = 0.1, label= "Control", col= "olivedrab4") +
  theme_bw()


### Heatmaps 


# select data for the 100 most highly expressed genes
sel100 <- order(rowMeans(x$cpms), decreasing=TRUE)[1:100]
highexprgenes_counts <- as.matrix(x$counts[sel100,])
colnames(highexprgenes_counts) <- colnames(x$counts)
higherorlog <- log(highexprgenes_counts,10)

## Heatmap manuscript
heatmap( 
  higherorlog[,c(5:16)],
  col=colorRampPalette(brewer.pal(11, "RdYlBu"))(100), 
  margin=c(8,5),
  labRow =  F
)
heat1.ggp <- ggplot2::ggplot_gtable(heat1)
heat1.grob <- ggplotGrob(heat1.ggp)

## Heatmap complementary
heatmap( 
  higherorlog,
  col=colorRampPalette(brewer.pal(11, "RdYlBu"))(100), 
  margin=c(8,5),
  labRow =  F
)

##### Volcano plots

## HI vs HM
HI.vs.HM.vol <- HI.vs.HM$table[, c(3, 13, 14)]
head(HI.vs.HM.vol)

keyvals.HI.vs.HM.vol <- ifelse(
  HI.vs.HM.vol$logFC < -2 & HI.vs.HM.vol$PValue < 0.01, 'blue',
  ifelse( HI.vs.HM.vol$logFC > 2 & HI.vs.HM.vol$PValue < 0.01, 'red',
  'grey'))

keyvals.HI.vs.HM.vol[is.na(keyvals.HI.vs.HM.vol)] <- 'grey'
names(keyvals.HI.vs.HM.vol)[keyvals.HI.vs.HM.vol == 'red'] <- 'up-regulated: FC > 2; p-val < 0.01'
names(keyvals.HI.vs.HM.vol)[keyvals.HI.vs.HM.vol == 'grey'] <- 'no sig: -2 < FC < 2; p-val > 0.01'
names(keyvals.HI.vs.HM.vol)[keyvals.HI.vs.HM.vol == 'blue'] <- 'down-regulated: FC < -2; p-val < 0.01'


##### Volcano plots
EnhancedVolcano(HI.vs.HM.vol, lab = "",#HI.vs.HM.vol[,1],
                x = "logFC", y = "PValue",
                pCutoff = 0.01, FCcutoff = 2,
                title = "HI vs HM",
                subtitle = "Actis Dato, V et al.",
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                colCustom = keyvals.HI.vs.HM.vol,
                ylim= c(0, 40),
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')


## HL vs HM
HL.vs.HM.vol <- HL.vs.HM$table[, c(3, 13, 14)]
head(HL.vs.HM.vol)

keyvals.HL.vs.HM.vol <- ifelse(
  HL.vs.HM.vol$logFC < -2 & HL.vs.HM.vol$PValue < 0.01, 'blue',
  ifelse( HL.vs.HM.vol$logFC > 2 & HL.vs.HM.vol$PValue < 0.01, 'red',
          'grey'))

keyvals.HL.vs.HM.vol[is.na(keyvals.HL.vs.HM.vol)] <- 'grey'
names(keyvals.HL.vs.HM.vol)[keyvals.HL.vs.HM.vol == 'red'] <- 'up-regulated: FC > 2; p-val < 0.01'
names(keyvals.HL.vs.HM.vol)[keyvals.HL.vs.HM.vol == 'grey'] <- 'no sig: -2 < FC < 2; p-val > 0.01'
names(keyvals.HL.vs.HM.vol)[keyvals.HL.vs.HM.vol == 'blue'] <- 'down-regulated: FC < -2; p-val < 0.01'

##### Volcano plots
EnhancedVolcano(HL.vs.HM.vol, lab = "",# HL.vs.HM.vol[,1],
                x = "logFC", y = "PValue",
                pCutoff = 0.01, FCcutoff = 2,
                title = "HL vs HM",
                subtitle = "Actis Dato, V et al.",
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                colCustom = keyvals.HL.vs.HM.vol,
                ylim= c(0, 40),
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')


## HIL vs HM
HIL.vs.HM.vol <- HIL.vs.HM$table[, c(3, 13, 14)]
head(HIL.vs.HM.vol)

keyvals.HIL.vs.HM.vol <- ifelse(
  HIL.vs.HM.vol$logFC < -2 & HIL.vs.HM.vol$PValue < 0.01, 'blue',
  ifelse( HIL.vs.HM.vol$logFC > 2 & HIL.vs.HM.vol$PValue < 0.01, 'red',
          'grey'))

keyvals.HIL.vs.HM.vol[is.na(keyvals.HIL.vs.HM.vol)] <- 'grey'
names(keyvals.HIL.vs.HM.vol)[keyvals.HIL.vs.HM.vol == 'red'] <- 'up-regulated: FC > 2; p-val < 0.01'
names(keyvals.HIL.vs.HM.vol)[keyvals.HIL.vs.HM.vol == 'grey'] <- 'no sig: -2 < FC < 2; p-val > 0.01'
names(keyvals.HIL.vs.HM.vol)[keyvals.HIL.vs.HM.vol == 'blue'] <- 'down-regulated: FC < -2; p-val < 0.01'

##### Volcano plots
EnhancedVolcano(HIL.vs.HM.vol, lab = "", # = HIL.vs.HM.vol[,1],
                x = "logFC", y = "PValue",
                pCutoff = 0.01, FCcutoff = 2,
                title = "HIL vs HM",
                subtitle = "Actis Dato, V et al.",
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                colCustom = keyvals.HIL.vs.HM.vol,
                ylim= c(0, 40),
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

