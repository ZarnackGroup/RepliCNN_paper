options(scipen = 999)

OKseqHMM <- function(
    chrsizes,
    fileOut,
    thresh,
    winS,
    binSize,
    fwdFiles,
    revFiles,
    hwinS=winS/2,
    st=c("D", "L", "H", "U"),
    sym=c("V", "W", "X", "Y", "Z"),
    pstart=rep(1/4, 4),
    pem=t(matrix(c(0.383886256, 0.255924171, 0.170616114, 0.113744076, 0.075829384,
                   .10,.20,.40,.20,.10,
                   .10,.20,.40,.20,.10,
                   0.022222222, 0.033333333, 0.066666667, 0.211111111, 0.666666667),
                   ncol=4)),
    ptrans=t(matrix(c(0.9999,0.000020,0,0.000080,
                      0,0.999,0,0.001,
                      0.001,0,0.999,0,
                      0.000080,0,0.000020,0.9999),
                      ncol=4)),
   quant=c(-1, -0.0082058939609862, -0.00141890249101162, 0.00103088286465956, 0.00800467305420799, 1))

{
  require(HMM)
  
  chromNames <- read.table(
    chrsizes,
    header=FALSE,
    sep="\t",
    comment.char="#",
    stringsAsFactors=FALSE
  )
  
  chr.sizes <- data.frame(
    chr=chromNames[,1],
    size=chromNames[,2]
  )
  
  # only use chromosomes for which files exist
  chrom <- intersect(
    chromNames[,1],
    names(fwdFiles)
  )
  
  for (i in c(1:length(chrom))){
    chr.name <- chrom[i]
    print(chr.name)
    chr.length <- chr.sizes[chr.sizes$chr == chr.name,2]
    print(chr.length)
    if (chr.length > winS*binSize*100){
      
      fileInF <- fwdFiles[[chr.name]]
      tags <- scan(fileInF, what=integer(), quiet=TRUE)
      tags[tags<=0] <- 1
      breaks <- seq(0, chr.length+binSize, by=binSize)
      h <- hist(tags, breaks=breaks, plot=FALSE)
      c <- h$counts
      
      fileInR <- revFiles[[chr.name]]
      tags <- scan(fileInR, what=integer(), quiet=TRUE)
      tags[tags<=0] <- 1
      breaks <- seq(0, chr.length+binSize, by=binSize)
      h <- hist(tags, breaks=breaks, plot=FALSE)
      w <- h$counts
      
      # raw polarity for later
      polar <- c/(c+w)
      polar[c<thresh & w<thresh] <- NA

      # 1kb RFD:
      rfd <- (c-w)/(w+c)
      rfd[is.na(rfd)] <- 0
      rfd[w<thresh & c<thresh] <- 0
      rfd[rfd > 1] <- 1
      rfd[rfd < -1] <- -1
      
      start_pos <- as.integer(breaks[1:length(breaks)-1])
      end_pos <- as.integer(breaks[2 : length(breaks)])
      chrName <- rep(chr.name,length(rfd))
      df <-data.frame(chr=chrName,startPos = start_pos,endPos = end_pos,rd_nb=rfd)

      #restrict the last position is the chr.length, not exceed that.
      df$endPos[nrow(df)] <- chr.length
      write.table(df, file = paste0(fileOut,"_",chr.name,"_RFD_cutoff_",thresh,"_bs_",binSize/1000,"kb.bedgraph", sep=""), append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)

      # smoothing RFD
      print(paste("Smoothing window size:", winS, "kb"))
      sw <- cumsum(w)
      lg <- length(w)
      from <- (-hwinS+2):(lg-hwinS+1)
      to <- from+winS-1
      from[from<1] <- 1
      to[to>lg] <- lg
      
      print(paste("Number of bins:", length(w)))
      
      win <- matrix(c(from, to), ncol=2)
      ws <- apply(win ,1, function(x) { (sw[x[2]]-sw[x[1]])/winS } )
      
      sc <- cumsum(c)
      lg <- length(c)
      cs <- apply(win ,1, function(x) { (sc[x[2]]-sc[x[1]])/winS } )
      
      
      print(paste("Cutoff:",thresh))
      rfd <- (cs-ws)/(ws+cs)
      rfd[is.na(rfd)] <- 0
      rfd[ws<thresh & cs<thresh] <- 0
      rfd[rfd > 1] <- 1
      rfd[rfd < -1] <- -1
      
      start_pos <- as.integer(breaks[1:length(breaks)-1])
      end_pos <- as.integer(breaks[2 : length(breaks)])
      chrName <- rep(chr.name,length(rfd))
      df <-data.frame(chr=chrName,startPos = start_pos,endPos = end_pos,rd_nb=rfd)
      #restrict the last position is the chr.length, not exceed that.
      df$endPos[nrow(df)] <- chr.length
      write.table(df, file = paste0(fileOut,"_",chr.name,"_RFD_cutoff_",thresh,"_bs_",binSize/1000,"kb_sm_",winS,"kb.bedgraph", sep=""), append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
      
      # derive
      bias <- cs/(ws+cs)
      bias[is.na(bias)] <- 0.5
      bias[ws<thresh & cs<thresh] <- 0.5
      
      delta <- c(0,bias[-1]-bias[-length(bias)])
      delta[is.na(delta)] <- 0.5
      
      # affect symbols
      if (is.na(quant[1])) { quant <- quantile(delta, probs = seq(0, 1, 0.20)) }
      quant[1] <- -1
      quant[length(quant)] <- 1
      
      print("Quantile borders :")
      print(quant)
      dx  <- unlist(sapply(delta, function(x) { ix <- which(x>=quant); ix[length(ix)] }))
      dx[dx>5] <- 5
      
      # HMM ==============================
      hmm1 <- initHMM(States=st, Symbols=sym, startProbs=pstart, transProbs=ptrans, emissionProbs=pem)
      print(hmm1)
      print("Viterbi...")
      seg <- viterbi(hmm=hmm1, observation=dx)
      
      # for showing the states profiles
      prof <- rep(NA, times=length(seg))
      prof[seg=="U"] <- 1
      prof[seg=="H"] <- 0.2
      prof[seg=="L"] <- -0.2
      prof[seg=="D"] <- -1
      prof[is.na(w[-length(w)])] <- -1.1
      prof[is.na(c[-length(c)])] <- -1.1
      prof[is.na(prof)] <- -.1
      
      write.table(prof, file = paste(fileOut,"_",chr.name,"_HMM.txt", sep=""), append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
      
      # adding the probability curve =========================
      print("Probabilities…")
      post <- posterior(hmm1,dx)
      prb <- rep(NA, length(seg))
      for (i in 1:length(seg))
      {
        prb[i] <- post[seg[i],i]
      }
      prb[prb>1] <- 1
      
      write.table(as.integer(prb*1000), file = paste(fileOut,"_",chr.name,"_HMMproba.txt", sep=""), append = T,
                  quote = FALSE, sep = "\t", col.names=F, row.names=F)
      
      # calculate the region coordinates by profil viterbi ========================
      
      left <- seg[-length(seg)]
      right <- seg[-1]
      ix <- which(left!=right)
      from <- c(1,ix+1)
      to <- c(ix,length(seg))
      
      # for the table, adding states and calculating lengths, slopes and associated probabilities ====================
      # calculate the slope of the polarity of the states (between the left and right part of the state))
      
      states <- meanPol <- ymax <- ymin <- p <- inc <- napc <- corr <- rep(NA, length(from))
      states <- seg[from]
      
      for (i in 1:length(from))
      {
        pos <- from[i]:to[i]
        lgpos <- length(pos)
        
        m <- prb[pos[!is.na(prb[pos])]]
        if (length(m)>0)
        {
          p[i] <- mean(m, na.rm=T)
        }
        meanPol[i] <- mean(polar[pos], na.rm=T)
        
        pos2 <- pos[!is.na(polar[pos])]
        lgpos2 <- length(pos2)
        napc[i] <- 100-round(lgpos2/lgpos*100)
        if (lgpos2<2)
        {
          res <- data.frame(NA, NA)
          inc[i] <- ymin[i] <- ymax[i] <- NA
        } else {
          realPos <- pos2*binSize
          res <- lm(polar[pos2] ~ realPos)
          inc[i] <- res[[1]][2]
          ymin[i] <- 	(inc[i]*realPos[1])+res[[1]][1]
          ymax[i] <- 	(inc[i]*realPos[lgpos2])+res[[1]][1]
          corr[i] <- cor(x=realPos, y=polar[pos2])
        }
      }
      
      ymin[ymin<0] <- 0
      ymax[ymax<0] <- 0
      ymin[ymin>1] <- 1
      ymax[ymax>1] <- 1
      
      # for display ==================================
      
      inc1 <- round(inc*10^8)		# as a percentage of RFD per megabase
      p <- round(p*100)
      chr <- rep(chr.name,length(from))
      from1 <- (from-1)*binSize+1
      to1 <- to*binSize
      lg1 <- to1-from1+1
      meanPol1 <- round(meanPol*100)
      ymin1 <- round(ymin*100)
      ymax1 <- round(ymax*100)
      corr1 <- round(corr*100)
      
      # adjustment of the border polarities by the average and recalculation of the slope ======================
      
      polL <- ymin1
      polR <- ymax1
      
      polD <- polR[-length(polR)]
      polG <- polL[-1]
      
      polD[is.na(polD)] <- polG[is.na(polD)]
      polG[is.na(polG)] <- polD[is.na(polG)]
      
      polM <- round((polG+polD)/2)
      polL[-1] <- polR[-length(polR)] <- polM
      
      # adjusted slope: 10 ^ 6 for clarity of the table display (% by megabase)
      slope_adj <- round(10^6*(polR-polL)/(to1-from1))
      
      
      # writing =============================
      dataOut <- data.frame(chr, from=as.integer(from1), to=as.integer(to1), state=states, length=lg1, slope=inc1,
                            p, pol_mean=meanPol1, pol_left=ymin1, pol_right=ymax1,
                            na=napc, cor=corr1, slope_adj=slope_adj, pol_adj_left=polL, pol_adj_right=polR)
      dataOut_U <- dataOut[dataOut$state == "U",]
      dataOut_D <- dataOut[dataOut$state == "D",]
      dataOut_HF <- dataOut[dataOut$state == "H",]
      dataOut_LF <- dataOut[dataOut$state == "L",]
      
      
      write.table(dataOut_U, file = paste(fileOut,"_",chr.name,"_HMMsegments_IZ.txt", sep=""), append = T, quote = FALSE, sep = "\t", col.names=T, row.names=F)
      write.table(dataOut_U[,1:3], file = paste(fileOut,"_",chr.name,"_HMMsegments_IZ.bed", sep=""), append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
      write.table(dataOut_D, file = paste(fileOut,"_",chr.name,"_HMMsegments_TZ.txt", sep=""), append = T, quote = FALSE, sep = "\t", col.names=T, row.names=F)
      write.table(dataOut_D[,1:3], file = paste(fileOut,"_",chr.name,"_HMMsegments_TZ.bed", sep=""), append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
      write.table(dataOut_HF, file = paste(fileOut,"_",chr.name,"_HMMsegments_highFlatZone.txt", sep=""), append = T, quote = FALSE, sep = "\t", col.names=T, row.names=F)
      write.table(dataOut_HF[,1:3], file = paste(fileOut,"_",chr.name,"_HMMsegments_highFlatZone.bed", sep=""), append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
      write.table(dataOut_LF, file = paste(fileOut,"_",chr.name,"_HMMsegments_LowFlatZone.txt", sep=""), append = T, quote = FALSE, sep = "\t", col.names=T, row.names=F)
      write.table(dataOut_LF[,1:3], file = paste(fileOut,"_",chr.name,"_HMMsegments_LowFlatZone.bed", sep=""), append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
    }
    else
    {
      i= i+1
    }
  }
}


# end of the function =========================================================================
setwd("PATH/okseq_hmm")

samples <- c(
  "imb_ulrich_2024_04_01_traelseq_hct116_asy",
  "sfb_ulrich_2025_01_01_hct116_wt_S1",
  "sfb_ulrich_2025_01_02_hct116_rad21_wt_S2"
)

for (sample in samples) {
  
  print(paste0("Running OKseqHMM for: ", sample))
  
  # input files
  fwdFiles <- setNames(
    paste0(
      "./data/by_chrom/",
      sample,
      ".umi.umidedup_sorted.fwd_chr",
      1:22,
      ".txt"
    ),
    paste0("chr", 1:22)
  )
  
  revFiles <- setNames(
    paste0(
      "./data/by_chrom/",
      sample,
      ".umi.umidedup_sorted.rev_chr",
      1:22,
      ".txt"
    ),
    paste0("chr", 1:22)
  )
  
  
  # output directory
  outdir <- paste0("./logs/", sample)
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  
  OKseqHMM(
    chrsizes="./hg38.chrom.sizes",
    fileOut=paste0(outdir, "/hmm"),
    thresh=0,
    winS=15,
    binSize=1000,
    fwdFiles=fwdFiles,
    revFiles=revFiles
  )
  
  print(paste0("Finished: ", sample))
}
