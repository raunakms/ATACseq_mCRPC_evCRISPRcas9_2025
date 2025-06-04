########################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                      ##
## RAUNAK SHRESTHA, PhD (https://github.com/raunakms/ATACseq_mCRPC_evCRISPRcas9_2025) ##
########################################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.wrk <- getwd()
dir.data <- file.path(dir.wrk, "data")

### DEFINE FILES ---
# NOTE GENCODEv28 REFERENCE ANNOTAITON DOWNLOADED FROM: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
# DOWNLOAD GENCODEv28 GTF FILE AND PLACE IT IN THE 'data' DIRECTORY ---
file.annot <- file.path("data/gencode.v28.annotation.gtf.gz") 

###################################################################################
### FUNCTION: parseGTF() ---
parseGTF <- function(file.gtf, feature.type){
    # CHROMOSOMES ---
    chromosome <- paste("chr", c(1:22,"X","Y"), sep="")

    # LOAD GTF FILE ---
    annot <- data.table::fread(file=file.gtf, sep="\t", header=FALSE, nThread=50, data.table=TRUE, stringsAsFactors=FALSE, skip=5, verbose=FALSE)
    annot <- subset(annot, annot$V1 %in% chromosome)

    # PARSE GTF COLUMNS ---
    if(feature.type == "gene"){    
        annot <- subset(annot, annot$V3 == "gene")
    }else if(feature.type == "exon"){
        annot <- subset(annot, annot$V3 == "exon")
    }else if(feature.type == "transcript"){
        annot <- subset(annot, annot$V3 == "transcript")
    }

    list.annot <- noquote(stringr::str_split(annot$V9, "; "))
    annot$EnsemblID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_id")]))), " "), function(x) x[2]))
    annot$Gene <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_name")]))), " "), function(x) x[2]))
    annot$GeneType <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_type")]))), " "), function(x) x[2]))

    # TRIM DATA ---
    if(feature.type == "gene"){ 
        annot <- subset(annot, annot$GeneType == "protein_coding")   
        items <- c("V1","V4","V5","Gene","EnsemblID","V7")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","EnsemblID","Strand")
    }else if(feature.type == "exon"){
        annot <- subset(annot, annot$GeneType == "protein_coding")
        list.annot <- stringr::str_split(annot$V9, "; ")
        annot$TranscriptID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "transcript_id")]))), " "), function(x) x[2]))
        annot$ExonID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "exon_id")]))), " "), function(x) x[2]))  
        items <- c("V1","V4","V5","Gene","EnsemblID","TranscriptID","ExonID","V7")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","GeneID","TranscriptID","ExonID","Strand")
    }else if(feature.type == "transcript"){
        annot$TranscriptID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "transcript_id")]))), " "), function(x) x[2]))
        items <- c("V1","V4","V5","Gene","EnsemblID","V7","TranscriptID","GeneType")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","EnsemblGeneID","Strand","EnsemblTranscriptID","GeneType")
    }        

    return(df)
}

###################################################################################
### LOAD ANNOTATION DATA ---
annot <- parseGTF(file.gtf=file.annot, feature.type="exon")

### EXTRACT EXONS FOR AR GENE ---
annot_exons <- annot[which(annot$Gene == "AR"),]

### WRITE OUTPUT ---
file.annot_exons <- file.path(dir.data, "exons_ar.tsv")
write.table(annot_exons, file.annot_exons, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
