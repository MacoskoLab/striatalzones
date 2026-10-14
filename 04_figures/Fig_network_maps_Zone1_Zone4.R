#!/usr/bin/env Rscript
# Zone 1 / Zone 4 cell-cell interaction network maps — FULL pipeline, r=200,
# Robust_95>=0.5. Nodes = cell types (size = total outgoing robust interactions,
# incl. autocrine self-loops); directed edges = sender->receiver robust-interaction
# count; split by zone-enrichment sign of Mean_Real_Delta. Recreates the original
# igraph circular-layout style (Force_Directed_Jan_19.R) at the new criteria.
suppressPackageStartupMessages({library(data.table); library(igraph)})

vf <- fread("/tmp/ef/validated_full.csv")[Scale==200 & Robust_95>=0.5]
cell_colors <- c(MSN="#E74C3C", Astrocyte="#3498DB", Microglia="#27AE60",
                 OPC="#F1C40F", Oligo="#9B59B6", Vascular="#E67E22")
cells <- names(cell_colors)
vf <- vf[Sender %in% cells & Receiver %in% cells]
z1 <- vf[Mean_Real_Delta > 0]        # Zone 1 (Patch) enriched
z4 <- vf[Mean_Real_Delta < 0]        # Zone 4 (Matrix) enriched   (self-loops kept)

# global maxima for shared scaling across the two panels
gmax_node <- max(c(z1[,.N,by=Sender]$N, z4[,.N,by=Sender]$N))
pc <- function(d) d[,.(n=.N),by=.(Sender,Receiver)]
gmax_pair <- max(c(pc(z1)$n, pc(z4)$n, 1))
node_size_func <- function(x,g) 16 + (x/g)*66
edge_width_func <- function(x) 0.5 + sqrt(x)*0.8

plot_zonal <- function(d, title, g_node, g_pair, leg=c("node","edge"), minedge=1){
  leg <- match.arg(leg)
  node_df <- d[,.(total=.N),by=Sender]                 # node size from ALL outgoing
  edge_df <- pc(d)[Sender %in% cells & Receiver %in% cells][n >= minedge]  # draw only >= minedge
  net <- graph_from_data_frame(edge_df, vertices=data.frame(name=cells), directed=TRUE)
  v_tot <- node_df$total[match(V(net)$name, node_df$Sender)]; v_tot[is.na(v_tot)] <- 0
  V(net)$size  <- node_size_func(v_tot, g_node)
  V(net)$color <- cell_colors[V(net)$name]
  E(net)$width      <- edge_width_func(edge_df$n)
  E(net)$color      <- adjustcolor(cell_colors[edge_df$Sender], alpha.f=0.72)  # colour by SENDER
  E(net)$arrow.size <- 0.5; E(net)$arrow.width <- 1.2; E(net)$curved <- 0.45
  # self-loop directions: radial-outward per node, but Oligo points LEFT
  node_ang <- setNames(2*pi*(seq_along(cells)-1)/length(cells), cells)
  node_ang["Oligo"] <- pi
  el <- ends(net, E(net), names=TRUE); isloop <- el[,1]==el[,2]
  la <- rep(0, ecount(net)); if(any(isloop)) la[isloop] <- node_ang[el[isloop,1]]
  E(net)$loop.angle <- la
  ld <- c(MSN=5.2, Astrocyte=3.0, Oligo=2.6, OPC=3.4, Microglia=2.6, Vascular=2.6)
  vdist <- ld[V(net)$name]; vdist[is.na(vdist)] <- 3.0

  par(mar=c(6,1,4,1))
  l <- layout_in_circle(net); vang <- -atan2(l[,2], l[,1]); l <- l*0.72
  plot(net, layout=l, rescale=FALSE, xlim=c(-1.55,1.55), ylim=c(-2.05,1.45),
       vertex.label.degree=vang, vertex.label.dist=vdist, vertex.label.cex=1.45,
       vertex.label.color="black", vertex.label.family="Helvetica", vertex.label.font=2,
       vertex.frame.color="grey30", main=title, cex.main=1.7)
  mtext("all edges, coloured by sender", side=3, line=-0.2, cex=0.9, col="grey45")
  if(leg=="node"){                                   # node-size legend (under Zone 1)
    nv <- c(100,300,500); nv <- nv[nv<=g_node]; ns <- node_size_func(nv, g_node)
    text(x=0, y=-1.40, labels="Outgoing Signals", cex=1.25, font=2)
    x0 <- -0.62
    for(i in seq_along(nv)){ xp <- x0 + (i-1)*0.62
      symbols(x=xp, y=-1.80, circles=ns[i]/300, inches=FALSE, add=TRUE, bg="grey90", fg="grey55")
      text(x=xp, y=-2.10, labels=nv[i], cex=1.05, font=2) }
  } else {                                           # edge-width legend (under Zone 4)
    ev <- c(1,100,200,300)
    text(x=0, y=-1.42, labels="Interactions", cex=1.25, font=2)
    x0 <- -0.95
    for(i in seq_along(ev)){ xp <- x0 + (i-1)*0.65
      segments(xp-0.18, -1.75, xp+0.18, -1.75, lwd=edge_width_func(ev[i]), col="grey55")
      text(x=xp, y=-2.02, labels=ev[i], cex=1.05, font=2) }
  }
}

# side-by-side (K/L) — larger canvas; PDF + SVG are vector, PNG raster preview
draw_both <- function(){
  par(mfrow=c(1,2))
  plot_zonal(z1,"Zone 1 Interaction Network",gmax_node,gmax_pair,leg="node")
  plot_zonal(z4,"Zone 4 Interaction Network",gmax_node,gmax_pair,leg="edge")
}
W <- 18; H <- 9.5
pdf("/tmp/ef/Fig_network_maps.pdf", width=W, height=H, family="Helvetica"); draw_both(); dev.off()
svglite::svglite("/tmp/ef/Fig_network_maps.svg", width=W, height=H); draw_both(); dev.off()
png("/tmp/ef/Fig_network_maps.png", width=W, height=H, units="in", res=200); draw_both(); dev.off()
cat("wrote Fig_network_maps.{pdf,svg,png}\n")
cat(sprintf("gmax_node=%d gmax_pair=%d\n", gmax_node, gmax_pair))
