#!/usr/bin/env Rscript
# Zone 1 / Zone 4 cell-cell interaction CHORD diagrams — FULL pipeline, r=200,
# Robust_95>=0.5. All directed sender->receiver channels as ribbons around the ring
# (width = interaction count), coloured by SENDER, directional arrows, self-loops
# on the diagonal. No center crossings.
suppressPackageStartupMessages({library(data.table); library(circlize)})

vf <- fread("/tmp/ef/validated_full.csv")[Scale==200 & Robust_95>=0.5]
cell_colors <- c(MSN="#E74C3C", Astrocyte="#3498DB", Microglia="#27AE60",
                 OPC="#F1C40F", Oligo="#9B59B6", Vascular="#E67E22")
cells <- names(cell_colors)
vf <- vf[Sender %in% cells & Receiver %in% cells]
z1 <- vf[Mean_Real_Delta > 0]; z4 <- vf[Mean_Real_Delta < 0]

mkmat <- function(d){
  m <- matrix(0, length(cells), length(cells), dimnames=list(cells,cells))
  cc <- d[,.(n=.N),by=.(Sender,Receiver)]
  for(i in seq_len(nrow(cc))) m[cc$Sender[i], cc$Receiver[i]] <- cc$n[i]
  m
}
m1 <- mkmat(z1); m4 <- mkmat(z4)
col_by_sender <- function(m){
  cm <- matrix(adjustcolor(cell_colors[rownames(m)], alpha.f=0.65), nrow(m), ncol(m))
  dimnames(cm) <- dimnames(m); cm
}

plot_chord <- function(m, title){
  circos.clear()
  circos.par(gap.degree=8, start.degree=90, points.overflow.warning=FALSE,
             canvas.xlim=c(-1,1), canvas.ylim=c(-1,1))
  chordDiagram(m, order=cells, grid.col=cell_colors, col=col_by_sender(m),
               directional=1, direction.type=c("diffHeight","arrows"),
               link.arr.type="big.arrow", diffHeight=convert_height(2,"mm"),
               self.link=1, annotationTrack=c("name","grid"),
               annotationTrackHeight=convert_height(c(3.5,4),"mm"),
               link.sort=TRUE, link.largest.ontop=TRUE)
  title(main=title, cex.main=1.7, line=0.2)
}

W <- 16; H <- 9
render <- function(){
  par(mfrow=c(1,2), mar=c(2,1,3.5,1), oma=c(2,0,0,0))
  plot_chord(m1, "Zone 1 Interaction Network")
  plot_chord(m4, "Zone 4 Interaction Network")
  mtext("ribbon width = interaction count · coloured by sender · arrows = sender → receiver",
        side=1, outer=TRUE, line=-1.5, cex=0.95, col="grey40")
  circos.clear()
}
pdf("/tmp/ef/Fig_chord_maps.pdf", width=W, height=H, family="Helvetica"); render(); dev.off()
svglite::svglite("/tmp/ef/Fig_chord_maps.svg", width=W, height=H); render(); dev.off()
png("/tmp/ef/Fig_chord_maps.png", width=W, height=H, units="in", res=200); render(); dev.off()
cat("wrote Fig_chord_maps.{pdf,svg,png}\n")
