# -----
# arg
# -----

args <- commandArgs(TRUE)
in.file <- paste0(args[which(args == "--infile")+1])

# -----
# lib
# -----

suppressMessages(library("DNAcopy"))
suppressMessages(library("jsonlite"))

# -----
# main
# -----

# process data

input <- read_json(in.file)
ratio <- as.numeric(unlist(input$results_r))
weights <- as.numeric(unlist(input$results_w))
gender <- input$ref_gender
alpha <- as.numeric(input$alpha)
out.file <- as.character(input$outfile)

if (gender == "M"){
    chrs = 1:24
} else {
    chrs = 1:23
}

bins.per.chr <- sapply(chrs, FUN = function(x) length(unlist(input$results_r[x])))
chr.end.pos <- c(0)
for (chr in chrs){
  l = bins.per.chr[chr]
  chr.end.pos <- c(chr.end.pos, l + chr.end.pos[length(chr.end.pos)])
}

# prepare for CBS

ratio[ratio == 0] = NA
weights[weights == 0] = 0.00001 # omit DNAcopy weirdness -- weight cannot be 0
# Only segments will be further processed with this data.
# Software will not account for 0's (these bins had no reference)

for.cbs <- as.data.frame(ratio)
chr.rep <- c()
chr.rep.2 <- c()
for (chr in chrs){
  chr.rep <- c(chr.rep, rep(chr, chr.end.pos[chr + 1] - chr.end.pos[chr]))
  chr.rep.2 <- c(chr.rep.2, 1:(chr.end.pos[chr + 1] - chr.end.pos[chr]))
}
for.cbs$chromosome <- chr.rep; for.cbs$x <- chr.rep.2
for.cbs <- for.cbs[, c(2,3,1)] ; colnames(for.cbs)[3] <- "y"

# Check for complete NA/chr

cbs.mask <- c()
for (chr in chrs){
    check <- which(for.cbs$chromosome == chr)
    if(!(all(is.na(for.cbs$y[check])))){
        cbs.mask <- c(cbs.mask, check)
    }
}
for.cbs <- for.cbs[cbs.mask,]

# CBS

CNA.object <- CNA(for.cbs$y, for.cbs$chromosome, for.cbs$x, data.type = "logratio", sampleid = "X")
f = file()
sink(file=f) ## silence output
CNA.object <- invisible(segment(CNA.object, alpha = as.numeric(alpha), verbose=1, weights=weights[cbs.mask])$output)
sink() ## undo silencing
close(f)

# Write output
write_json(t(CNA.object), out.file)