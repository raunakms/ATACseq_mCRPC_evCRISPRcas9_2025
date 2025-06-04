########################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                      ##
## RAUNAK SHRESTHA, PhD (https://github.com/raunakms/ATACseq_mCRPC_evCRISPRcas9_2025) ##
########################################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("Gviz")
library("GenomicRanges")
library("IRanges")

### DEFINE PATH ---
dir.wrk <- getwd()
dir.data <- file.path(dir.wrk, "data")
dir.plots <- file.path(dir.wrk, "plots")

### DEFINE FILES ---
file.bed_sgRNA <- file.path(dir.data, "sgRNA_ar.bed")
file.annot_exons <- file.path(dir.data, "exons_ar.tsv")
file.atacseq_peaks <- file.path(dir.data, "mcrpc_atacseq_peaks_ar.rds")

###############################################################################################################
sampleids_cell_lines <- c("22Rv1_1","22Rv1_2",
							"C4-2_1","C4-2_2",
							"LNCaP_1","LNCaP_2",
							"VCaP_1","VCaP_2")


sampleids_pdx <- c("MSK-PCa2_1","MSK-PCa2_2",
					"MSK-PCa19_1","MSK-PCa19_2",
					"MSK-PCa22",
					"PDX09a","PDX10a",
					"PDX16a","PDX16b",
					"PDX34a")

sampleids_mcrpc <- c("DTB-090-PRO","DTB-156-BL","DTB-128-BL","DTB-176-BL","DTB-069-BL",
						"DTB-149-BL","DTB-101-BL","DTB-119-PRO",
						"DTB-009-BL","DTB-019-BL","DTB-019-PRO","DTB-022-BL","DTB-035-BL",
						"DTB-077-PRO","DTB-085-BL","DTB-111-PRO","DTB-127-BL","DTB-127-PRO",
						"DTB-129-BL","DTB-143-BL","DTB-146-BL","DTB-148-BL","DTB-167-BL",
						"DTB-172-BL","DTB-206-BL","DTB-222-BL")



########################################################################################
### LOAD AR EXONS ---
annot_exons <- data.table::fread(file=file.annot_exons, sep="\t", header=TRUE, nThread=1, data.table=FALSE, stringsAsFactors=FALSE, verbose=FALSE)
annot_exons <- subset(annot_exons, select=c("Chrom","Start","End","Strand","GeneID","ExonID","TranscriptID","Gene"))
colnames(annot_exons) <- c("chromosome","start","end","strand","gene","exon","transcript","symbol")


### LOAD sgRNA BED ---
dat_bed_sgrna <- data.table::fread(file=file.bed_sgRNA, sep="\t", header=FALSE, nThread=1, data.table=FALSE, stringsAsFactors=FALSE, verbose=FALSE)
colnames(dat_bed_sgrna) <- c("chr","start","end","name","strand","seq")
dat_bed_sgrna$strand <- "*"
dat_bed_sgrna$id <- c("sgRNA1-LBD","sgRNA2-NTD","sgRNA3-NTD")
gr_bed_sgrna <- GenomicRanges::makeGRangesFromDataFrame(dat_bed_sgrna, keep.extra.columns=TRUE)



### LOAD ATAC-SEQ PEAKS ---
list.atacseq_peaks <- readRDS(file=file.atacseq_peaks)
list.gr_tang <- list.atacseq_peaks$tang
list.gr_wcdt <- list.atacseq_peaks$wcdt




########################################################################################
############### FIGURE 3A-1: sgRNA1-LBD ###############
# PLOT PARAMS ---
plot_chr <- "chrX"
plot_start <- 67716900
plot_end <- 67718300

# CONSTRUCT PLOTTING TRACK ---
track_ideogram <- Gviz::IdeogramTrack(genome="hg38", chromosome="chrX", sizes=1.5)
track_genome_axis <- Gviz::GenomeAxisTrack(col="#000000", from=67717000, to=67718200, ticksAt=seq(67717000,67718200, by=300), lwd=0.8, sizes=1)
track_gene_region <- Gviz::GeneRegionTrack(annot_exons, genome="hg38", chromosome=plot_chr, stacking="dense", name="AR gene", lwd=0.5, col="#000000", fill="#ffff99", rot.title=0, sizes=1)
track_peaks_sgrna <- Gviz::AnnotationTrack(range=gr_bed_sgrna[1], stacking="full", name="sgRNA1-LBD", group="sgRNA1-LBD", groupAnnotation="group", just.group="below", fontsize.group=9, col="#000000", fill="#000000", rot.title=0)
track_peaks_1 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[1] ]], stacking="squish", name=sampleids_cell_lines[1], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_2 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[2] ]], stacking="squish", name=sampleids_cell_lines[2], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_3 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[3] ]], stacking="squish", name=sampleids_cell_lines[3], col="#1f78b4", fill="#1f78b4", rot.title=0, sizes=1)
track_peaks_4 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[4] ]], stacking="squish", name=sampleids_cell_lines[4], col="#1f78b4", fill="#1f78b4", rot.title=0, sizes=1)
track_peaks_5 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[5] ]], stacking="squish", name=sampleids_cell_lines[5], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_6 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[6] ]], stacking="squish", name=sampleids_cell_lines[6], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_7 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[7] ]], stacking="squish", name=sampleids_cell_lines[7], col="#b15928", fill="#b15928", rot.title=0, sizes=1)
track_peaks_8 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[8] ]], stacking="squish", name=sampleids_cell_lines[8], col="#b15928", fill="#b15928", rot.title=0, sizes=1)

track_highlight <- Gviz::HighlightTrack(trackList = list(track_gene_region, 
                        track_peaks_1, track_peaks_2, track_peaks_3, track_peaks_4,
                        track_peaks_5, track_peaks_6, track_peaks_7, track_peaks_8,
                        track_peaks_sgrna),
                        start=67717473, end=67717522, chromosome = plot_chr, col="#d9d9d9", fill="#d9d9d9", alpha=0.8, lwd=0.5)



# PLOT ---
file_plot <- file.path(dir.plots, "fig_3a_1.pdf")
pdf(file_plot, height=2.5, width=4)
    Gviz::plotTracks(list(track_ideogram, track_genome_axis, track_highlight), 
                        from = plot_start, to = plot_end,
                        #extend.left = 100, extend.right = 100,
                        cex=0.5, cex.axis=0.3, cex.main=0.3, cex.title=0.5,
                        margin=8, innerMargin=1, title.width=2, 
                        col.axis="#000000", col.title="#000000", background.title="#FFFFFF")
dev.off()





########################################################################################
############### FIGURE 3A-2: sgRNA2-NTD ###############
# PLOT PARAMS ---
plot_chr <- "chrX"
plot_start <- 67544150
plot_end <- 67545950

# CONSTRUCT PLOTTING TRACK ---
track_ideogram <- Gviz::IdeogramTrack(genome="hg38", chromosome="chrX", sizes=1.5)
track_genome_axis <- Gviz::GenomeAxisTrack(col="#000000", lwd=0.8, sizes=1)
track_gene_region <- Gviz::GeneRegionTrack(annot_exons, genome="hg38", chromosome=plot_chr, stacking="dense", name="AR gene", lwd=0.5, col="#000000", fill="#ffff99", rot.title=0, sizes=1)
track_peaks_sgrna <- Gviz::AnnotationTrack(range=gr_bed_sgrna[2], stacking="full", name="sgRNA2-NTD", group="sgRNA2-NTD", groupAnnotation="group", just.group="below", fontsize.group=9, col="#000000", fill="#000000", rot.title=0)
track_peaks_1 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[1] ]], stacking="squish", name=sampleids_cell_lines[1], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_2 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[2] ]], stacking="squish", name=sampleids_cell_lines[2], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_3 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[3] ]], stacking="squish", name=sampleids_cell_lines[3], col="#1f78b4", fill="#1f78b4", rot.title=0, sizes=1)
track_peaks_4 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[4] ]], stacking="squish", name=sampleids_cell_lines[4], col="#1f78b4", fill="#1f78b4", rot.title=0, sizes=1)
track_peaks_5 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[5] ]], stacking="squish", name=sampleids_cell_lines[5], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_6 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[6] ]], stacking="squish", name=sampleids_cell_lines[6], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_7 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[7] ]], stacking="squish", name=sampleids_cell_lines[7], col="#b15928", fill="#b15928", rot.title=0, sizes=1)
track_peaks_8 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_cell_lines[8] ]], stacking="squish", name=sampleids_cell_lines[8], col="#b15928", fill="#b15928", rot.title=0, sizes=1)

track_highlight <- Gviz::HighlightTrack(trackList = list(track_gene_region, 
                        track_peaks_1, track_peaks_2, track_peaks_3, track_peaks_4,
                        track_peaks_5, track_peaks_6, track_peaks_7, track_peaks_8,
                        track_peaks_sgrna),
                        start=67545671, end=67545720, chromosome = plot_chr, col="#d9d9d9", fill="#d9d9d9", alpha=0.8, lwd=0.5)

# PLOT ---
file_plot <- file.path(dir.plots, "fig_3a_2.pdf")
pdf(file_plot, height=2.5, width=4)
    Gviz::plotTracks(list(track_ideogram, track_genome_axis, track_highlight), 
                        from = plot_start, to = plot_end,
                        #extend.left = 100, extend.right = 100,
                        cex=0.5, cex.axis=0.3, cex.main=0.3, cex.title=0.5,
                        margin=8, innerMargin=1, title.width=2,
                        col.axis="#000000", col.title="#000000", background.title="#FFFFFF")
dev.off()

########################################################################################





########################################################################################
############### FIGURE S4A: sgRNA2-NTD ###############
# PLOT PARAMS ---
plot_chr <- "chrX"
plot_start <- 67544150
plot_end <- 67545950

# CONSTRUCT PLOTTING TRACK ---
track_ideogram <- Gviz::IdeogramTrack(genome="hg38", chromosome="chrX", sizes=1.5)
track_genome_axis <- Gviz::GenomeAxisTrack(col="#000000", lwd=0.8, sizes=1)
track_gene_region <- Gviz::GeneRegionTrack(annot_exons, genome="hg38", chromosome=plot_chr, stacking="dense", name="AR gene", lwd=0.5, col="#000000", fill="#ffff99", rot.title=0, sizes=1)
track_peaks_sgrna <- Gviz::AnnotationTrack(range=gr_bed_sgrna[2], stacking="full", name="sgRNA2-NTD", group="sgRNA2-NTD", groupAnnotation="group", just.group="below", fontsize.group=9, col="#000000", fill="#000000", rot.title=0)

track_peaks_1 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[1] ]], stacking="squish", name=sampleids_pdx[1], col="#33a02c", fill="#33a02c", rot.title=0, sizes=1)
track_peaks_2 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[2] ]], stacking="squish", name=sampleids_pdx[2], col="#33a02c", fill="#33a02c", rot.title=0, sizes=1)
track_peaks_3 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[3] ]], stacking="squish", name=sampleids_pdx[3], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_4 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[4] ]], stacking="squish", name=sampleids_pdx[4], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_5 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[5] ]], stacking="squish", name=sampleids_pdx[5], col="#cab2d6", fill="#cab2d6", rot.title=0, sizes=1)
track_peaks_6 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[6] ]], stacking="squish", name=sampleids_pdx[6], col="#fb9a99", fill="#fb9a99", rot.title=0, sizes=1)
track_peaks_7 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[7] ]], stacking="squish", name=sampleids_pdx[7], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_8 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[8] ]], stacking="squish", name=sampleids_pdx[8], col="#b15928", fill="#b15928", rot.title=0, sizes=1)
track_peaks_9 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[9] ]], stacking="squish", name=sampleids_pdx[9], col="#b15928", fill="#b15928", rot.title=0, sizes=1)
track_peaks_10 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[10] ]], stacking="squish", name=sampleids_pdx[10], col="#b2df8a", fill="#b2df8a", rot.title=0, sizes=1)

track_highlight <- Gviz::HighlightTrack(trackList = list(track_gene_region, 
                        track_peaks_1, track_peaks_2, track_peaks_3, track_peaks_4, track_peaks_5, 
                        track_peaks_6, track_peaks_7, track_peaks_8, track_peaks_9, track_peaks_10, 
                        track_peaks_sgrna),
                        start=67545671, end=67545720, chromosome = plot_chr, col="#d9d9d9", fill="#d9d9d9", alpha=0.8, lwd=0.5)


# PLOT ---
file_plot <- file.path(dir.plots, "fig_S4a.pdf")
pdf(file_plot, height=2.5, width=4)
    Gviz::plotTracks(list(track_ideogram, track_genome_axis, track_highlight), 
                        from = plot_start, to = plot_end,
                        #extend.left = 100, extend.right = 100,
                        cex=0.5, cex.axis=0.3, cex.main=0.3, cex.title=0.5,
                        margin=8, innerMargin=1, title.width=2,
                        col.axis="#000000", col.title="#000000", background.title="#FFFFFF")
dev.off()







########################################################################################
############### FIGURE S4B: sgRNA1-LBD ###############
# PLOT PARAMS ---
plot_chr <- "chrX"
plot_start <- 67716900
plot_end <- 67718300

# CONSTRUCT PLOTTING TRACK ---
track_ideogram <- Gviz::IdeogramTrack(genome="hg38", chromosome="chrX", sizes=1.5)
track_genome_axis <- Gviz::GenomeAxisTrack(col="#000000", from=67717000, to=67718200, ticksAt=seq(67717000,67718200, by=300), lwd=0.8, sizes=1)
track_gene_region <- Gviz::GeneRegionTrack(annot_exons, genome="hg38", chromosome=plot_chr, stacking="dense", name="AR gene", lwd=0.5, col="#000000", fill="#ffff99", rot.title=0, sizes=1)
track_peaks_sgrna <- Gviz::AnnotationTrack(range=gr_bed_sgrna[1], stacking="full", name="sgRNA1-LBD", group="sgRNA1-LBD", groupAnnotation="group", just.group="below", fontsize.group=9, col="#000000", fill="#000000", rot.title=0)

track_peaks_1 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[1] ]], stacking="squish", name=sampleids_pdx[1], col="#33a02c", fill="#33a02c", rot.title=0, sizes=1)
track_peaks_2 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[2] ]], stacking="squish", name=sampleids_pdx[2], col="#33a02c", fill="#33a02c", rot.title=0, sizes=1)
track_peaks_3 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[3] ]], stacking="squish", name=sampleids_pdx[3], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_4 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[4] ]], stacking="squish", name=sampleids_pdx[4], col="#e31a1c", fill="#e31a1c", rot.title=0, sizes=1)
track_peaks_5 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[5] ]], stacking="squish", name=sampleids_pdx[5], col="#cab2d6", fill="#cab2d6", rot.title=0, sizes=1)
track_peaks_6 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[6] ]], stacking="squish", name=sampleids_pdx[6], col="#fb9a99", fill="#fb9a99", rot.title=0, sizes=1)
track_peaks_7 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[7] ]], stacking="squish", name=sampleids_pdx[7], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_8 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[8] ]], stacking="squish", name=sampleids_pdx[8], col="#b15928", fill="#b15928", rot.title=0, sizes=1)
track_peaks_9 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[9] ]], stacking="squish", name=sampleids_pdx[9], col="#b15928", fill="#b15928", rot.title=0, sizes=1)
track_peaks_10 <- Gviz::AnnotationTrack(range=list.gr_tang[[ sampleids_pdx[10] ]], stacking="squish", name=sampleids_pdx[10], col="#b2df8a", fill="#b2df8a", rot.title=0, sizes=1)


track_highlight <- Gviz::HighlightTrack(trackList = list(track_gene_region, 
                        track_peaks_1, track_peaks_2, track_peaks_3, track_peaks_4, track_peaks_5, 
                        track_peaks_6, track_peaks_7, track_peaks_8, track_peaks_9, track_peaks_10, 
                        track_peaks_sgrna),
                        start=67717473, end=67717522, chromosome = plot_chr, col="#d9d9d9", fill="#d9d9d9", alpha=0.8, lwd=0.5)


# PLOT ---
file_plot <- file.path(dir.plots, "fig_S4b.pdf")
pdf(file_plot, height=2.5, width=4)
    Gviz::plotTracks(list(track_ideogram, track_genome_axis, track_highlight), 
                        from = plot_start, to = plot_end,
                        #extend.left = 100, extend.right = 100,
                        cex=0.5, cex.axis=0.3, cex.main=0.3, cex.title=0.5,
                        margin=8, innerMargin=1, title.width=2,
                        col.axis="#000000", col.title="#000000", background.title="#FFFFFF")
dev.off()





########################################################################################
############### FIGURE S4C: sgRNA2-NTD ###############
# PLOT PARAMS ---
plot_chr <- "chrX"
plot_start <- 67544150
plot_end <- 67545950

# CONSTRUCT PLOTTING TRACK ---
track_ideogram <- Gviz::IdeogramTrack(genome="hg38", chromosome="chrX", sizes=1.5)
track_genome_axis <- Gviz::GenomeAxisTrack(col="#000000", lwd=0.8, sizes=1)
track_gene_region <- Gviz::GeneRegionTrack(annot_exons, genome="hg38", chromosome=plot_chr, stacking="dense", name="AR gene", lwd=0.5, col="#000000", fill="#ffff99", rot.title=0, sizes=1)
track_peaks_sgrna <- Gviz::AnnotationTrack(range=gr_bed_sgrna[2], stacking="full", name="sgRNA2-NTD", group="sgRNA2-NTD", groupAnnotation="group", just.group="below", fontsize.group=9, col="#000000", fill="#000000", rot.title=0)

track_peaks_1 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[1] ]], stacking="squish", name=sampleids_mcrpc[1], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_2 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[2] ]], stacking="squish", name=sampleids_mcrpc[2], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_3 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[3] ]], stacking="squish", name=sampleids_mcrpc[3], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_4 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[4] ]], stacking="squish", name=sampleids_mcrpc[4], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_5 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[5] ]], stacking="squish", name=sampleids_mcrpc[5], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_6 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[6] ]], stacking="squish", name=sampleids_mcrpc[6], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_7 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[7] ]], stacking="squish", name=sampleids_mcrpc[7], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_8 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[8] ]], stacking="squish", name=sampleids_mcrpc[8], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_9 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[9] ]], stacking="squish", name=sampleids_mcrpc[9], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_10 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[10] ]], stacking="squish", name=sampleids_mcrpc[10], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_11 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[11] ]], stacking="squish", name=sampleids_mcrpc[11], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_12 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[12] ]], stacking="squish", name=sampleids_mcrpc[12], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_13 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[13] ]], stacking="squish", name=sampleids_mcrpc[13], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_14 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[14] ]], stacking="squish", name=sampleids_mcrpc[14], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_15 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[15] ]], stacking="squish", name=sampleids_mcrpc[15], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_16 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[16] ]], stacking="squish", name=sampleids_mcrpc[16], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_17 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[17] ]], stacking="squish", name=sampleids_mcrpc[17], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_18 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[18] ]], stacking="squish", name=sampleids_mcrpc[18], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_19 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[19] ]], stacking="squish", name=sampleids_mcrpc[19], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_20 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[20] ]], stacking="squish", name=sampleids_mcrpc[20], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_21 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[21] ]], stacking="squish", name=sampleids_mcrpc[21], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_22 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[22] ]], stacking="squish", name=sampleids_mcrpc[22], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_23 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[23] ]], stacking="squish", name=sampleids_mcrpc[23], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_24 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[24] ]], stacking="squish", name=sampleids_mcrpc[24], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_25 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[25] ]], stacking="squish", name=sampleids_mcrpc[25], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_26 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[26] ]], stacking="squish", name=sampleids_mcrpc[26], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)

track_highlight <- Gviz::HighlightTrack(trackList = list(track_gene_region, 
                        track_peaks_1, track_peaks_2, track_peaks_3, track_peaks_4, track_peaks_5, 
                        track_peaks_6, track_peaks_7, track_peaks_8, track_peaks_9, track_peaks_10, 
                        track_peaks_11, track_peaks_12, track_peaks_13, track_peaks_14, track_peaks_15, 
                        track_peaks_16, track_peaks_17, track_peaks_18, track_peaks_19, track_peaks_20,
                        track_peaks_21, track_peaks_22, track_peaks_23, track_peaks_24, track_peaks_25, track_peaks_26,
                        track_peaks_sgrna),
                        start=67545671, end=67545720, chromosome = plot_chr, col="#d9d9d9", fill="#d9d9d9", alpha=0.8, lwd=0.5)


# PLOT ---
file_plot <- file.path(dir.plots, "fig_S4c.pdf")
pdf(file_plot, height=4, width=4)
    Gviz::plotTracks(list(track_ideogram, track_genome_axis, track_highlight), 
                        from = plot_start, to = plot_end,
                        #extend.left = 100, extend.right = 100,
                        cex=0.5, cex.axis=0.3, cex.main=0.3, cex.title=0.5,
                        margin=8, innerMargin=1, title.width=2,
                        col.axis="#000000", col.title="#000000", background.title="#FFFFFF")
dev.off()



########################################################################################
############### FIGURE S4D: sgRNA1-LBD ###############
# PLOT PARAMS ---
plot_chr <- "chrX"
plot_start <- 67716900
plot_end <- 67718300

# CONSTRUCT PLOTTING TRACK ---
track_ideogram <- Gviz::IdeogramTrack(genome="hg38", chromosome="chrX", sizes=1.5)
track_genome_axis <- Gviz::GenomeAxisTrack(col="#000000", from=67717000, to=67718200, ticksAt=seq(67717000,67718200, by=300), lwd=0.8, sizes=1)
track_gene_region <- Gviz::GeneRegionTrack(annot_exons, genome="hg38", chromosome=plot_chr, stacking="dense", name="AR gene", lwd=0.5, col="#000000", fill="#ffff99", rot.title=0, sizes=1)
track_peaks_sgrna <- Gviz::AnnotationTrack(range=gr_bed_sgrna[1], stacking="full", name="sgRNA1-LBD", group="sgRNA1-LBD", groupAnnotation="group", just.group="below", fontsize.group=9, col="#000000", fill="#000000", rot.title=0)

track_peaks_1 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[1] ]], stacking="squish", name=sampleids_mcrpc[1], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_2 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[2] ]], stacking="squish", name=sampleids_mcrpc[2], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_3 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[3] ]], stacking="squish", name=sampleids_mcrpc[3], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_4 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[4] ]], stacking="squish", name=sampleids_mcrpc[4], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_5 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[5] ]], stacking="squish", name=sampleids_mcrpc[5], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_6 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[6] ]], stacking="squish", name=sampleids_mcrpc[6], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_7 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[7] ]], stacking="squish", name=sampleids_mcrpc[7], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_8 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[8] ]], stacking="squish", name=sampleids_mcrpc[8], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_9 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[9] ]], stacking="squish", name=sampleids_mcrpc[9], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_10 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[10] ]], stacking="squish", name=sampleids_mcrpc[10], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_11 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[11] ]], stacking="squish", name=sampleids_mcrpc[11], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_12 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[12] ]], stacking="squish", name=sampleids_mcrpc[12], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_13 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[13] ]], stacking="squish", name=sampleids_mcrpc[13], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_14 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[14] ]], stacking="squish", name=sampleids_mcrpc[14], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_15 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[15] ]], stacking="squish", name=sampleids_mcrpc[15], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_16 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[16] ]], stacking="squish", name=sampleids_mcrpc[16], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_17 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[17] ]], stacking="squish", name=sampleids_mcrpc[17], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_18 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[18] ]], stacking="squish", name=sampleids_mcrpc[18], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_19 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[19] ]], stacking="squish", name=sampleids_mcrpc[19], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_20 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[20] ]], stacking="squish", name=sampleids_mcrpc[20], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_21 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[21] ]], stacking="squish", name=sampleids_mcrpc[21], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_22 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[22] ]], stacking="squish", name=sampleids_mcrpc[22], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_23 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[23] ]], stacking="squish", name=sampleids_mcrpc[23], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_24 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[24] ]], stacking="squish", name=sampleids_mcrpc[24], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_25 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[25] ]], stacking="squish", name=sampleids_mcrpc[25], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)
track_peaks_26 <- Gviz::AnnotationTrack(range=list.gr_wcdt[[ sampleids_mcrpc[26] ]], stacking="squish", name=sampleids_mcrpc[26], col="#a6cee3", fill="#a6cee3", rot.title=0, sizes=1)

track_highlight <- Gviz::HighlightTrack(trackList = list(track_gene_region, 
                        track_peaks_1, track_peaks_2, track_peaks_3, track_peaks_4, track_peaks_5, 
                        track_peaks_6, track_peaks_7, track_peaks_8, track_peaks_9, track_peaks_10, 
                        track_peaks_11, track_peaks_12, track_peaks_13, track_peaks_14, track_peaks_15, 
                        track_peaks_16, track_peaks_17, track_peaks_18, track_peaks_19, track_peaks_20,
                        track_peaks_21, track_peaks_22, track_peaks_23, track_peaks_24, track_peaks_25, track_peaks_26,
                        track_peaks_sgrna),
                        start=67717473, end=67717522, chromosome = plot_chr, col="#d9d9d9", fill="#d9d9d9", alpha=0.8, lwd=0.5)

# PLOT ---
file_plot <- file.path(dir.plots, "fig_S4d.pdf")
pdf(file_plot, height=4, width=4)
    Gviz::plotTracks(list(track_ideogram, track_genome_axis, track_highlight), 
                        from = plot_start, to = plot_end,
                        #extend.left = 100, extend.right = 100,
                        cex=0.5, cex.axis=0.3, cex.main=0.3, cex.title=0.5,
                        margin=8, innerMargin=1, title.width=2,
                        col.axis="#000000", col.title="#000000", background.title="#FFFFFF")
dev.off()



