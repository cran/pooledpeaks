## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(pooledpeaks)

## ----message=FALSE------------------------------------------------------------
library(Fragman)
library(ape)
library(magrittr)
library(tibble)
library(dplyr)

## -----------------------------------------------------------------------------
file_path <- system.file("extdata", package = "pooledpeaks")
eggcount <- data.frame(
    ID = c("X104.1",  "X1084.1", "X1084.3", "X1086.3", "X1087.3", "X1205.3",
           "X121.3",  "X1222.3", "X1354.3", "X1453.3", "X1531.3", "X1540.1",
           "X1550.3", "X1796.1", "X1809.1", "X1968.1", "X1968.3", "X2100.1",
           "X2462.1", "X2463.1", "X473.1",  "X620.1",  "X620.3", "X679.1",
           "X910.1",  "X910.3"),
    n = c(192, 126, 185, 171, 140, 20, 46, 80, 156, 154, 122, 19, 45, 117, 75,
          22, 175, 100, 97, 183, 67, 90, 157, 104, 195, 145)
  )
Shae10 <- c(161,164,167,170,173,176,179,182,185,188,191,194,197,200,203,206,209,
            212,215,218)
mic_SMMS2 <- c(211, 215, 219, 223, 227, 231, 235, 239)
GS600LIZ <- c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220,
              240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414,
              420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600)

## -----------------------------------------------------------------------------
fsa_data <- fsa_batch_imp(file_path, channels = 5, rawPlot = FALSE,
                              fourier = TRUE, saturated = TRUE,
                              lets.pullup = FALSE)
fsa_data <- associate_dyes(fsa_data, file_path)

## ----message=FALSE,eval=FALSE-------------------------------------------------
# ladder.info.attach(stored = fsa_data,ladder = GS600LIZ,
#                    ladd.init.thresh = 200, prog = FALSE, draw = FALSE)
# corro <- unlist(sapply(list.data.covarrubias, function(x){x$corr}))
# bad <- which(corro < .999)

## ----message=FALSE------------------------------------------------------------
scores_SMMS2 <- score_markers_rev3(my.inds = fsa_data,
                                   channel = 1,
                                   channel.ladder = 5,
                                   panel = "mic_SMMS2",
                                   ladder = GS600LIZ,
                                   init.thresh = 100,
                                   ploidy = length(mic_SMMS2),
                                   shift = 1,
                                   windowL = 1,
                                   windowR= 1,
                                   left.cond = c(0, 2.5),
                                   right.cond = 0,
                                   pref = 1,
                                   plotting = FALSE
                                   )

scores_Shae10 <- score_markers_rev3(my.inds = fsa_data,
                                   channel = 1,
                                   channel.ladder = 5,
                                   panel = "Shae10",
                                   ladder = GS600LIZ,
                                   init.thresh = 100,
                                   ploidy = length(Shae10),
                                   shift = 1,
                                   windowL = 1,
                                   windowR= 1,
                                   left.cond = c(0, 2.5),
                                   right.cond = 0,
                                   pref = 1,
                                   plotting = FALSE
                                   )

## -----------------------------------------------------------------------------
scores_SMMS2_lf<-clean_scores(scores_SMMS2, pattern1 = "_FA.*",replacement1 = "",
                              pattern2 = "_Sample.*", replacement2 = "")

scores_Shae10_lf<-clean_scores(scores_Shae10, pattern1 = "_FA.*",replacement1 = "",
                              pattern2 = "_Sample.*", replacement2 = "")

## -----------------------------------------------------------------------------
scores_SMMS2_tdf <- lf_to_tdf(scores_SMMS2_lf)

scores_Shae10_tdf <- lf_to_tdf(scores_Shae10_lf)

## ----eval=FALSE---------------------------------------------------------------
# write.table(scores_SMMS2_lf, file = "scores_SMMS2_lfex.txt", col.names = NA,
#             quote = FALSE, row.names = TRUE, sep = "\t")
# write.table(scores_SMMS2_tdf, file = "scores_SMMS2_tdfex.txt", col.names = NA,
#             quote = FALSE, row.names = TRUE, sep = "\t")
# 
# write.table(scores_Shae10_lf, file = "scores_Shae10_lfex.txt", col.names = NA,
#             quote = FALSE, row.names = TRUE, sep = "\t")
# write.table(scores_Shae10_tdf, file = "scores_Shae10_tdfex.txt", col.names = NA,
#             quote = FALSE, row.names = TRUE, sep = "\t")
# 

## -----------------------------------------------------------------------------
SMMS2<- read.delim("./scores_SMMS2_tdfex.txt")

## ----include=FALSE------------------------------------------------------------
SMMS2<- SMMS2%>%
  column_to_rownames(var = "X")%>%
  select(-contains(".fsa"))

## -----------------------------------------------------------------------------
head(SMMS2[, 1:9])

## -----------------------------------------------------------------------------
SMMS2_IDM <- data_manipulation(SMMS2, threshold = 200)
head(SMMS2_IDM[, 1:9])

## -----------------------------------------------------------------------------
SMMS2_repcheck <- Rep_check(SMMS2_IDM)
head(SMMS2_repcheck[, 3:11])

## -----------------------------------------------------------------------------
SMMS2_PCM<-PCDM(SMMS2_repcheck,eggcount,'SMMS2')
head(SMMS2_PCM[,1:6])

## ----eval=FALSE---------------------------------------------------------------
# combined3<-rbind.fill(SMMS2_PCM, SMMS13_PCM, SMMS16_PCM)
# 
# write.table(combined3, file = "combined3.txt", col.names = NA,
#             quote = FALSE, row.names = TRUE, sep = "\t")

## -----------------------------------------------------------------------------
gends <- LoadData("./combined3.txt")
head(gends[1:8])

## -----------------------------------------------------------------------------
N <- TypedLoci(gends)
head(N[,1:5])


## -----------------------------------------------------------------------------
J <- GeneIdentityMatrix(gends,N)
head(J[,1:5])


## -----------------------------------------------------------------------------
D <- GeneticDistanceMatrix(J)
head(D[,1:5])

## -----------------------------------------------------------------------------
print(head(GST(J)[,1:5]))

print(head(JostD(J)[,1:5]))


## ----fig.width=6, fig.height=4------------------------------------------------
M <- MDSplot(D,pcs=c(1,2))

## ----fig.width=6, fig.height=4------------------------------------------------
Tr <- nj(D)
Tr <- ladderize(Tr)
plot(Tr,cex=0.5,no.margin = TRUE,type='phylogram')

