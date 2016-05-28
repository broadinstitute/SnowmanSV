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
  ff2 <- ff[N_ALT <= 3 & T_ALT >= 2 & pmax(mapq1, mapq2) >= 30 & T_COV < 400 & N_COV < 400 & confidence!="NOLOCAL"]
  setkey(ff2, chr1, pos1, strand1, chr2, pos2, strand2, N_ALT, T_ALT)
  setkey(ff2, chr1, pos1, strand1, chr2, pos2, strand2)
  ff2 <- unique(ff2)
  saveRDS(ff2, paste0("data/160413_",x,"x_sim.snowman.bps.rds"))
}

.load_bedpe <- function(x) {
  ff <- fread(x)
  setnames(ff, paste0("V",seq(11)), c("seqnames", "start","end", 'altchr','altpos','altend',"ID","SPAN","strand","altstrand","TYPE")) 
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
  ff[, c("end","altend","TYPE") := NULL]
  
  setkey(ff, chr, pos)
  return(list(grl=grl.ff, dt=ff))
}

## load other data
#map <- readRDS("data/wgEncodeCrgMapabilityAlign100mer.gr.rds")

### load the NA12878 data
print("...loading snowman NA12878")
snow.NA12878 <- .load_snowman_vcf_datatable("data/NA12878_160525.snowman.sv.vcf")
snow.NA12878[, deltype := strand[1] == '+' && altstrand[1] == "-" && altchr[1]==chr[1] && SPAN[1] < 1e6, by=RARID]
snow.dels.NA12878 <- snow.NA12878[!duplicated(RARID) & deltype]

## annotate lumpy data
truth.NA12878 <- readRDS("data/lumpy_na12878.rds")
#truth.NA12878 <- .load_bedpe("data/na12878.deletion.lumpyfile4.bedpe")
#gr <- with(truth.NA12878$dt, GRanges(chr, IRanges(pos,pos)))
#fo <- gr.val(gr, map, 'score')
#truth.NA12878$dt$mean_mappability <- fo$score
#gr <- with(truth.NA12878$dt, GRanges(altchr, IRanges(altpos,altpos)))
#fo <- gr.val(gr, map, 'score')
#truth.NA12878$dt$mean_mappability_alt <- fo$score
#saveRDS(truth.NA12878, "data/lumpy_na12878.rds")

## annotate lumpy data
truth2.NA12878 <- readRDS("data/lumpy_na12878_set2.rds")
truth2.NA12878 <- .load_bedpe2("data/lumpy_na12878_truthset4.bedpe")
#gr <- with(truth.NA12878$dt, GRanges(chr, IRanges(pos,pos)))
#fo <- gr.val(gr, map, 'score')
#truth.NA12878$dt$mean_mappability <- fo$score
#gr <- with(truth.NA12878$dt, GRanges(altchr, IRanges(altpos,altpos)))
#fo <- gr.val(gr, map, 'score')
#truth.NA12878$dt$mean_mappability_alt <- fo$score
#saveRDS(truth2.NA12878, "data/lumpy_na12878_set2.rds")

## prepare as GRangesLists
tmp <- snow.NA12878[snow.NA12878$deltype]
gr <- with(tmp, GRanges(chr, IRanges(pos, pos), strand=strand, RARID=RARID))
grl.snow.dels.na12878 <- split(gr, gr$RARID)
mcols(grl.snow.dels.na12878) <- tmp[!duplicated(RARID),.(SPAN)]
stopifnot(as.numeric(names(table(elementLengths(grl.snow.dels.na12878))))==2)

## do the overlaps with truth set
ro_snowman_set1 <- ra.overlaps(grl.snow.dels.na12878,  truth.NA12878$grl, pad=2e3)
TP_snowman1 <- length(unique(ro_snowman_set1[,'ra1.ix'])) ## TRUE POSITIVE
ro_snowman_set2 <- ra.overlaps(grl.snow.dels.na12878, truth2.NA12878$grl, pad=2e3)
TP_snowman2 <- length(unique(ro_snowman_set2[,'ra1.ix'])) ## TRUE POSITIVE
FP_snowman1 <- sum(ix <- mcols(grl.snow.dels.na12878)$SPAN > 100) - TP_snowman1
FP_snowman2 <- sum(ix <- mcols(grl.snow.dels.na12878)$SPAN > 100) - TP_snowman2

### NA12878 Pindel
ff <- readRDS("data/pindel_NA12878.rds")
pindel.dels <- ff[SVLEN > 30 & SVTYPE=="DEL"]
gr <- with(pindel.dels, GRanges(c(V1, V1), IRanges(c(V2, V2), width=1), strand=rep(c("+","-"), each=nrow(pindel.dels)), id=rep(seq(nrow(pindel.dels)),2)))
grl.pindel <- split(gr, gr$id) 
ro_pindel_set1 <- ra.overlaps(grl.pindel,  truth.NA12878$grl, pad=2e3)
TP_pindel1 <- length(unique(ro_pindel_set1[,'ra1.ix'])) ## TRUE POSITIVE
ro_pindel_set2 <- ra.overlaps(grl.pindel, truth2.NA12878$grl, pad=2e3)
TP_pindel2 <- length(unique(ro_pindel_set2[,'ra1.ix'])) ## TRUE POSITIVE
FP_pindel1 <- nrow(pindel.dels[SVLEN > 100]) - TP_pindel1
FP_pindel2 <- nrow(pindel.dels[SVLEN > 100]) - TP_pindel2

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

setkey(snow.NA12878, RARID)
print("done loading data")

