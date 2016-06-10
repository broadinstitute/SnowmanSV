source("batch.functions.R")
library(skitools)
library(gUtils)
library(gTrack)
library(data.table)
library(skidb)

Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")
           
.load_snowman_vcf_datatable <- function(x) {
  
  suppressWarnings(ff <- fread(paste("grep -v '^#'", x),sep='\t'))
  setnames(ff, paste0("V",seq(10)), c("chr","pos","id","REF","ALT","QUAL","FILTER","INFO","GENOINFO","GENOTYPE"))
  ff[, c("REF","QUAL","FILTER","GENOINFO") := NULL]  
  ff[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
  ff[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
  ff[, HOMSEQ := ifelse(grepl("HOMSEQ", INFO), gsub(".*?HOMSEQ=([A-Z]+).*", "\\1", INFO), "")]
  ff[, REPSEQ := ifelse(grepl("REPSEQ", INFO), gsub(".*?REPSEQ=([A-Z]+).*", "\\1", INFO), "")]
  ff[, SCTG   := gsub(".*?SCTG=(.*?);.*", "\\1", INFO)]
  ff[, MAPQ   := as.numeric(gsub(".*?MAPQ=([0-9]+).*", "\\1", INFO))]
  ff[, INFO := NULL]
  ff[, RARID := gsub("([0-9]+):.*","\\1", id)]
  ff[, c("id") := NULL]
  
  ff[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')] 
  ff[, altstrand := rev(strand), by=c("RARID")]
  ff[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
  ff[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
  
  return(ff)
}

.load_bedpe2 <- function(x) {
  ff <- fread(x)
  setnames(ff, paste0("V", seq(10)), c("seqnames","start","end","altchr","altpos","altend","ID","d","strand","altstrand"))
  ff[ , seqnames := gsub("chr","", seqnames)]
  ff[ , altchr := gsub("chr","", altchr)]
  
  ff2 <- data.table::copy(ff) 
  setnames(ff2,c("seqnames","start","end","strand","altchr","altpos","altend","altstrand"),c("altchr","altpos","altend", "altstrand","seqnames", "start","end","strand"))
  ff2$strand <- "-"
  dels <- sort(gr.fix(dt2gr(ff), si))
  ff_double <- rbind(ff,ff2)
  dels_double <- sort(gr.fix(gr.nochr(dt2gr(ff_double)), si))
  grl.ff <- split(dels_double, dels_double$ID)
  
  setnames(ff, c("seqnames", "start"), c("chr","pos"))
  ff[, c("end","altend") := NULL]
  
  setkey(ff, chr, pos)
  
  return(list(grl=grl.ff, dt=ff))
}

.prepare_bps <- function(x) {
  
  ff <- fread(paste0("gunzip -c data/160413_",x,"x_sim.bps.txt.gz"))
  setnames(ff, c("n001_/xchip/gistic/Jeremiah/Projects/SnowmanPaper/Benchmark/160410/normfull_0.6_subsampled.bam",paste0("t000_/xchip/gistic/Jeremiah/Projects/SnowmanPaper/Benchmark/160410/sim_wnormal_",x,"x_d1.rh.bam")), c("normal","tumor"))
  ff[, N_ALT := as.integer(gsub("0/1:([0-9]+):.*","\\1", normal))]
  ff[, T_ALT := as.integer(gsub("0/1:([0-9]+):.*","\\1", tumor))]
  ff[, N_COV := as.integer(gsub("0/1:[0-9]+:([0-9]+):.*","\\1", normal))]
  ff[, T_COV := as.integer(gsub("0/1:[0-9]+:([0-9]+):.*","\\1", tumor))]
  ff[, DBSNP2 := nchar(DBSNP) > 1]
  ff[, c("normal","tumor", "quality", "DBSNP", "pon_samples","graylist", "secondary_alignment") := NULL]
  ff2 <- ff[N_ALT <= 3 & T_ALT >= 2 & confidence!="NOLOCAL"] # & pmax(mapq1, mapq2) >= 30 & T_COV < 400 & N_COV < 400 & confidence!="NOLOCAL"]

  setkey(ff2, chr1, pos1, strand1, chr2, pos2, strand2, N_ALT, T_ALT)
  setkey(ff2, chr1, pos1, strand1, chr2, pos2, strand2)
  ff2 <- unique(ff2)
  saveRDS(ff2, paste0("data/160413_",x,"x_sim.snowman.bps.rds"))
}

.load_bedpe <- function(x) {
  ff <- fread(x)
  setnames(ff, paste0("V",seq(11)), c("seqnames", "start1","end1", 'altchr','altpos','altend',"ID","SPAN","strand","altstrand","TYPE")) 
  ff[ , seqnames := gsub("chr","", seqnames)]
  ff[ , altchr := gsub("chr","", altchr)]
  ff2 <- data.table::copy(ff) 
  setnames(ff2,c("seqnames","start1","end1","strand","altchr","altpos","altend","altstrand"),c("altchr","altpos","altend", "altstrand","seqnames", "start1","end1","strand"))
  ff2$strand <- "-"
  #ff2[, altstart := mean(start1, end1), by=ID]
  #dels <- sort(gr.fix(dt2gr(ff), si))
  ff_double <- rbind(ff,ff2)
  ff_double[, start := (start1+end1)/2, by=ID]
  ff_double[, end := start]
  
  dels_double <- sort(gr.fix(gr.nochr(dt2gr(ff_double)), si))
  grl.ff <- split(dels_double, dels_double$ID)
  
  #setnames(ff, c("seqnames", "start"), c("chr","pos"))
  #ff[, c("end","altend","TYPE") := NULL]
  
  #setkey(ff, chr, pos)
  return(list(grl=grl.ff, dt=ff))
}

## load the truth NA12878 from lumpy paper
truth2.NA12878 <- readRDS("data/lumpy_na12878_set2.rds")
truth.NA12878 <- .load_bedpe("data/na12878.deletion.lumpyfile4.bedpe")

## funciton for analyzing na12878 deletions
.na12878_olap <- function(x) {
  
  # if ("indel" %in% names(x)) {
  #   
  #   x2 <- x$indel[x$indel$SPAN > 50 & x$indel$ETYPE == "del"]
  #   x3_1 <- GRanges(seqnames(x2), IRanges(start(x2), width=1), strand="+", id=seq_along(x2))
  #   x3_2 <- GRanges(seqnames(x2), IRanges(end(x2), width=1), strand="-", id=seq_along(x2))
  #   b <- c(x3_1, x3_2)
  #   b <- split(b, b$id)
  #   mcols(b)$deltype <- TRUE
  #   mcols(b)$SPAN <- x2$SPAN
  #   mcols(b)$pacbio <- FALSE
  #   x$grl <- grlbind(x$grl, b)
  #   x$dt$id = x$dt$SCTG
  #     
  # }
    
  if ("deltype" %in% colnames(mcols(x$grl))) {
    dd = x$grl[mcols(x$grl)$deltype]
  } else if ("SVTYPE" %in% colnames(x$grl)) { ## delly
    dd = x$grl[mcols(x$grl)$SVTYPE=="DEL"]
  } else {
    dd = x$grl
    x$dt$deltype <- TRUE
    mcols(dd)$deltype <- TRUE
  }
  
  if ("SVLEN" %in% colnames(mcols(dd))) {
      mcols(dd)$SPAN = mcols(dd)$SVLEN
      x$dt[x$dt$ALT >= 6]
      dd <- dd[!is.na(mcols(dd)$ALT) & mcols(dd)$ALT >= 6]
  }
  
  if (!"id" %in% colnames(x$dt))
      x$dt$id <- seq(nrow(x$dt))
    
  if ("pacbio" %in% colnames(mcols(dd)))
    pacbio <- mcols(dd)$pacbio
  else
    pacbio <- rep(FALSE, length(dd))
  print(paste("Pacbio hits", sum(pacbio)))
  
  #suppressWarnings(ro1 <- ra.overlaps(dd,  truth.NA12878$grl, pad=2e2, ignore.strand=TRUE))
  #TP1e <- dd[unique(c(ro1[,'ra1.ix'], which(pacbio)))]
  #TP1 <- length(unique(ro1[,'ra2.ix'])) ## TRUE POSITIVE
  suppressWarnings(ro2 <- ra.overlaps(dd, truth2.NA12878$grl, pad=2e2, ignore.strand=TRUE))
  TP2e <- dd[unique(c(ro2[,'ra1.ix'], which(pacbio)))]
  TP2 <- length(unique(ro2[,'ra2.ix'])) ## TRUE POSITIVE
  #FP1e <- dd[setdiff(seq_along(dd), c(ro1[,'ra1.ix'], which(pacbio)))]
  #FP1 = sum(mcols(FP1e)$SPAN >= 100 & mcols(FP1e)$SPAN < 1e5)
  FP2e <- dd[setdiff(seq_along(dd), c(ro2[,'ra1.ix'], which(pacbio)))]
  FP2 = sum(mcols(FP2e)$SPAN >= 100 & mcols(FP2e)$SPAN < 1e5)

  TP1=FP1=FP1e=TP1e=0

  #setkey(x$dt, seqnames, start)
  return(list(dels=x$dt[!duplicated(x$dt$id) & x$dt$deltype & x$dt$SPAN >= 50], dels.grl=x$grl[mcols(x$grl)$deltype], FP1=FP1, FP2=FP2, TP1=TP1, TP2=TP2, ro1=ro1, ro2=ro2, TP1e=TP1e,TP2e=TP2e, FP1e=FP1e,FP2e=FP2e))
}

# ## do the overlaps with truth set SNOWMAN NA12878
#saveRDS(snowman12878 <- .na12878_olap(readRDS("data/snowman_na12878.rds")), "data/snowman_processed_na12878.rds")
snowman12878 <- readRDS("data/snowman_processed_na12878.rds")
#saveRDS(lumpy12878 <- .na12878_olap(readRDS("data/lumpy_na12878.rds")), "data/lumpy_processed_na12878.rds")
lumpy12878 <- readRDS("data/lumpy_processed_na12878.rds")
#saveRDS(pindel12878 <- .na12878_olap(readRDS("data/pindel_na12878.rds")), "data/pindel_processed_na12878.rds")
pindel12878 <- readRDS("data/pindel_processed_na12878.rds")

## NA12878 ROC
lum <- readRDS("data/lumpy_na12878.rds")
lum.out <- rbindlist(lapply(2:20, function(x) {
  print(x)
  
  llo <- lum
  llo$grl <- lum$grl[ix <- mcols(lum$grl)$SR + mcols(lum$grl)$PE >= x]
  ll2 <- .na12878_olap(llo)
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="LUMPY"))
  return(dt)
  
}))
saveRDS(lum.out, "data/lumpy_roc_na12878.rds")
lumpy_roc_na12878 <- readRDS("data/lumpy_roc_na12878.rds")

## SNOWMAN
ss <- readRDS("data/snowman_na12878.rds")
ss$grl <- ss$grl[-unique(c(which(pmax(mcols(ss$grl)$DISC_MAPQ, mcols(ss$grl)$MAPQ) <= 10),which(pmax(mcols(ss$grl)$NM, mcols(ss$grl)$MATENM) >= 10)))]
ss.out <- rbindlist(lapply(seq(2,20,1), function(x) {
  print(x)
  
  llo <- ss
  llo$grl <- ss$grl[ix <- mcols(ss$grl)$NORMALT >= x & mcols(ss$grl)$EVDNC %in% c("ASSMB", "ASDIS")]
  #llo$indel <- ss$indel[ss$indel$DP >= x]
  ll2 <- .na12878_olap(llo)

  llo <- ss
  llo$grl <- ss$grl[ix <- mcols(ss$grl)$NORMALT >= x]
  #llo$indel <- ss$indel[ss$indel$DP >= x]
  ll3 <- .na12878_olap(llo)
  
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="Snowman (AS)"))
  dt <- rbind(dt, with(ll3, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="Snowman (RP + AS)")))
  return(dt)
}))
saveRDS(ss.out, "data/snowman_roc_na12878.rds")
snowman_roc_na12878 <- readRDS("data/snowman_roc_na12878.rds")

## PINDEL
ss <- readRDS("data/pindel_na12878.rds")
ss.out <- rbindlist(lapply(2:20, function(x) {
  print(x)
  
  llo <- ss
  llo$grl <- ss$grl[ix <- !is.na(mcols(ss$grl)$ALT) & mcols(ss$grl)$ALT >= x]
  ll2 <- .na12878_olap(llo)
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="Pindel"))
}))
saveRDS(ss.out, "data/pindel_roc_na12878.rds")
pindel_roc_na12878 <- readRDS("data/pindel_roc_na12878.rds")

## DELLY
ss <- readRDS("data/delly_na12878.rds")
ss$grl <- ss$grl[mcols(ss$grl)$SVTYPE =="DEL" & mcols(ss$grl)$SPAN > 50]
ss.out <- rbindlist(lapply(2:20, function(x) {
  print(x)
  llo <- ss
  llo$grl <- ss$grl[ix <- !is.na(mcols(ss$grl)$DISC) & (mcols(ss$grl)$DISC + mcols(ss$grl)$SPLIT) >= x]
  ll2 <- .na12878_olap(llo)
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="DELLY (RP + SR)"))
  
  llo <- ss
  llo$grl <- ss$grl[ix <- !is.na(mcols(ss$grl)$DISC) & mcols(ss$grl)$DISC >= x]
  ll2 <- .na12878_olap(llo)
  dt <- rbind(dt, with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="DELLY (RP)")))
  
}))
saveRDS(ss.out, "data/delly_roc_na12878.rds")
delly_roc_na12878 <- readRDS("data/delly_roc_na12878.rds")

ggplot(data=rbindlist(list(delly_roc_na12878, snowman_roc_na12878, lumpy_roc_na12878, pindel_roc_na12878)), aes(x=FP1, y=TP1, color=caller)) + geom_point() + geom_line() + theme_bw() + scale_x_continuous(breaks=seq(0, 3000, by=500)) + xlab("Not in validation set") + ylab("True Positive") #+ scale_color_manual(values=("snowman"="blue","Pindel"="red","LUMPY"="darkgreen"), name="")
ggplot(data=rbindlist(list(delly_roc_na12878, snowman_roc_na12878, lumpy_roc_na12878, pindel_roc_na12878)), aes(x=FP2, y=TP2, color=caller)) + geom_point() + geom_line() + theme_bw() + xlab("Not in validation set") + scale_x_continuous(breaks=seq(0, 3000, by=250)) + ylab("True Positive") + scale_y_continuous(breaks=seq(1200,2900,by=100)) + coord_cartesian(xlim=c(0,3000), ylim=c(1500,2900))

###########################
## SIMLULATED
###########################

print("...loading simulated data")
events.d1 <- fread("data/events.d1.txt")
setnames(events.d1, paste0("V",seq(11)), c("chr","pos","altchr","altpos","strand","altstrand","dummy","span","class","ins_seq","ID"))
events.d1[class == "del"]$ins_seq <- ""
events.d1[class == "ins"]$ins_seq <- events.d1[class == "ins"]$dummy
events.d1[, dummy := NULL]
setkey(events.d1, chr, pos)
gr.events <- sort(gr.fix(with(events.d1, GRanges(c(chr,altchr), IRanges(c(pos,altpos), width=1), strand=c(strand,altstrand), id=c(ID,ID), span=c(span,span), type=c(class,class), ins_seq=c(ins_seq, ins_seq))),si))

#.prepare_bps(2)
#.prepare_bps(5)
#.prepare_bps(10)

## load pindel
#pindel <- readRDS("data/pindel_sim.rds")

snow_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr=dat
  datr$sv = datr$sv[mcols(datr$sv)$TUMALT >= x]
  datr$indel = datr$indel[mcols(datr$indel)$AD >= x]
  dt <- flag.plot(xsv=datr$sv, xindel=datr$indel, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Snowman"
  dt$ALT = x
  return(dt)
})
rb <- rbindlist(snow_sim_results)

ff <- readRDS("data/pindel_sim.rds")[[2]]
pindel_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[ff$TUMALT >= x & ff$NORMAL <= 0]
  dt <- flag.plot(xindel=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Pindel"
  dt$ALT = x
  return(dt)
})
pindel_sim_results <- rbindlist(pindel_sim_results)

ff <- readRDS("data/platypus_sim.rds")[[2]]
platypus_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[ff$NV_t >= x]
  dt <- flag.plot(xindel=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Platypus"
  dt$ALT = x
  return(dt)
})
platypus_sim_results <- rbindlist(platypus_sim_results)

ff <- readRDS("data/lumpy_sim.rds")[[2]]
mcols(ff)$TUMALT <- mcols(ff)$SR + mcols(ff)$PE
mcols(ff)$SPAN <- mcols(ff)$span
lumpy_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[mcols(ff)$TUMALT >= x]
  dt <- flag.plot(xsv=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Lumpy"
  dt$ALT = x
  return(dt)
})
lumpy_sim_results <- rbindlist(lumpy_sim_results)

ff <- readRDS("data/strelka_sim.rds")[[2]]
ff$SPAN <- ff$CALC_SPAN
strelka_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[mcols(ff)$TALT >= x]
  dt <- flag.plot(xindel=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Strelka"
  dt$ALT = x
  return(dt)
})
strelka_sim_results <- rbindlist(strelka_sim_results)

ggplot(data=rbindlist(list(rb[Group=="SV"],lumpy_sim_results[Group=="SV"], pindel_sim_results[Group=="SV"])), aes(x=FP, y=TP, color=caller))+ geom_point() + geom_line() + theme_bw() + xlab("False Positive") + scale_x_continuous(breaks=seq(0, 15, by=1)) + coord_cartesian(xlim=c(0,15))
ggplot(data=rbindlist(list(rb[Group=="Indel"],strelka_sim_results[Group=="Indel"],pindel_sim_results[Group=="Indel" & FP < 1000],platypus_sim_results[Group=="Indel"])), aes(x=FP, y=TP, color=caller))+ geom_point() + geom_line() + theme_bw() + xlab("False Positive") + ylab("True Positive") 
ggplot(data=rbindlist(list(rb[Group=="Medium"],lumpy_sim_results[Group=="Medium"], pindel_sim_results[Group=="Medium" & FP < 500],platypus_sim_results[Group=="Medium"])), aes(x=FP, y=TP, color=caller))+ geom_point() + geom_line() + theme_bw() + xlab("False Positive") + ylab("True Positive") # + scale_y_continuous(breaks=seq(0,3000,500)) + scale_x_continuous(breaks=seq(0,5,20))

print("done loading data")

