#!/usr/bin/env Rscript
# ==============================================================================
# fig_4panel_combo.R  —  FINAL figure (self-contained).
#
# Four thin panels (full-page-width strip), left to right:
#   a  Robustness vs radius : robust-interaction count vs neighbourhood radius,
#                             each cell type as % of its own maximum.
#   b  200 vs 400 µm        : per-sender turnover (same/new/lost) of the robust
#                             set, 200 µm vs 400 µm, as % of the larger (max) set.
#   c  Neighbor distances   : mean distance to the k-th nearest partner (k=1..5),
#                             per cell type.
#   d  K=2 vs K=5 radius     : per-sender turnover between each sender's K=2 and
#                             K=5 density-matched radii (nearest 100 µm), % of max.
#
# Run from THIS folder:  Rscript fig_4panel_combo.R
# Reads:  data/robustness_by_radius.csv , data/kdist.csv
# Writes: Fig_4panel_combo.pdf , Fig_4panel_combo.png
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork); library(scales) })
OUTDIR <- "."
INK <- "#14171C"; MUT <- "#6A7178"; GRY <- "#B9BEB8"
pal <- c(MSN="#2a78d6", Astrocyte="#eb6834", Oligo="#1baf7a", OPC="#eda100", Vascular="#4a3aa7", Microglia="#e34948")
ord <- c("MSN","Astrocyte","Oligo","OPC","Vascular","Microglia")
abbr <- c(MSN="MSN", Astrocyte="Astro", Oligo="Oligo", OPC="OPC", Vascular="Vasc", Microglia="Micro")
catcol <- c("same"="#9AA0A6","new"="#2E6E6A","lost"="#c94a4a")

# ---- inputs -----------------------------------------------------------------
# robustness_by_radius.csv: one row per (Sender,Receiver,Interaction,Scale) with
#   Robust_95 = fraction of 18 donors whose |Zone1-Zone4 delta| beats its own
#   permutation null at that neighbourhood radius (µm). "Robust" = Robust_95>=0.5.
rbr <- fread("data/robustness_by_radius.csv")
rbr[, key := paste(Sender,Receiver,Interaction,sep="||")]
# kdist.csv: mean distance (µm) from a cell to its k-th nearest neighbour of each type
kd  <- fread("data/kdist.csv")

# set turnover between two radii for one sender, normalised to the larger (max) set
mkchurn <- function(dt, s, sc_a, sc_b){
  a <- dt[Sender==s & Scale==sc_a & Robust_95>=0.5, key]; b <- dt[Sender==s & Scale==sc_b & Robust_95>=0.5, key]
  data.table(same=length(intersect(a,b)), new=length(setdiff(b,a)), lost=length(setdiff(a,b)),
             n0=max(length(a),length(b))) }
churn_long <- function(sc_a, sc_b_fun){
  rbindlist(lapply(ord, function(s){ d <- mkchurn(rbr, s, if(is.function(sc_a)) sc_a(s) else sc_a,
                                                          if(is.function(sc_b_fun)) sc_b_fun(s) else sc_b_fun)
    data.table(Sender=s, cat=c("same","new","lost"), pct=100*c(d$same,d$new,d$lost)/d$n0) }))[
    , Sender:=factor(Sender,ord)][, cat:=factor(cat,c("same","new","lost"))] }

# ---- panel data -------------------------------------------------------------
ra <- rbr[Scale %in% c(50,100,200,300,400,500,700,1000) & Robust_95>=0.5, .N, by=.(Sender,Scale)][Sender %in% ord]
ra[, pctmax := 100*N/max(N), by=Sender]; ra[, cell := factor(Sender, levels=ord)]

Lb <- churn_long(200, 400)                                          # panel b

dk <- kd[k<=5, .(um=mean(med_dist)), by=.(receiver,k)][receiver %in% ord]; dk[, cell := factor(receiver, levels=ord)]

dkr <- kd[, .(um=mean(med_dist)), by=.(receiver,k)]                 # per-cell-type radii for panel d
R2 <- sapply(ord, function(x) round(dkr[receiver==x & k==2, um]/100)*100)
R5 <- sapply(ord, function(x) round(dkr[receiver==x & k==5, um]/100)*100)
Ld <- churn_long(function(s) R2[[s]], function(s) R5[[s]])          # panel d

# ---- theme ------------------------------------------------------------------
thin <- theme_classic(base_size=8) +
  theme(axis.text=element_text(size=6.3, colour=MUT), axis.text.x=element_text(angle=40, hjust=1),
        axis.title=element_text(size=6.8, colour=MUT),
        axis.line=element_line(linewidth=0.3, colour=GRY), axis.ticks=element_line(linewidth=0.3, colour=GRY),
        plot.title=element_text(size=8.2, face="bold", colour=INK),
        plot.subtitle=element_text(size=5.9, colour=MUT, margin=margin(b=3)),
        legend.position="bottom", legend.text=element_text(size=6.4, colour=INK),
        legend.title=element_blank(), legend.key.size=unit(7,"pt"),
        legend.margin=margin(0,0,0,0), plot.margin=margin(5,5,4,5))

pa <- ggplot(ra, aes(Scale, pctmax, colour=cell)) +
  geom_line(linewidth=0.6) + geom_point(size=0.9, stroke=0) +
  scale_colour_manual(values=pal, name=NULL) +
  scale_x_log10(breaks=c(50,100,200,400,1000)) +
  scale_y_continuous(limits=c(0,103), breaks=c(0,50,100), labels=function(x)paste0(x,"%")) +
  labs(title="a  Robustness vs radius", subtitle="% of each cell's own max", x="radius (µm, log)", y="% of max") + thin

pb <- ggplot(Lb, aes(Sender, pct, fill=cat)) +
  geom_col(position=position_dodge(width=0.8), width=0.74) +
  scale_fill_manual(values=catcol, name=NULL) + scale_x_discrete(labels=abbr) +
  scale_y_continuous(limits=c(0,104), breaks=c(0,50,100), labels=function(x)paste0(x,"%"), expand=expansion(mult=c(0,0.04))) +
  labs(title="b  200 vs 400 µm", subtitle="turnover, % of larger set", x=NULL, y="% of max") + thin

pc <- ggplot(dk, aes(k, um, colour=cell)) +
  geom_line(linewidth=0.6) + geom_point(size=1.0, stroke=0) +
  scale_colour_manual(values=pal, name=NULL, guide="none") +
  scale_x_continuous(breaks=1:5) + scale_y_continuous(limits=c(0,NA), labels=comma) +
  labs(title="c  Neighbor distances", subtitle="dist to k-th partner", x="k-th partner", y="distance (µm)") + thin

pd <- ggplot(Ld, aes(Sender, pct, fill=cat)) +
  geom_col(position=position_dodge(width=0.8), width=0.74) +
  scale_fill_manual(values=catcol, name=NULL) + scale_x_discrete(labels=abbr) +
  scale_y_continuous(limits=c(0,104), breaks=c(0,50,100), labels=function(x)paste0(x,"%"), expand=expansion(mult=c(0,0.04))) +
  labs(title="d  K=2 vs K=5 radius", subtitle="per-sender radii, % of larger", x=NULL, y="% of max") + thin

fig <- (pa | pb | pc | pd) + plot_layout(guides="collect") & theme(legend.position="bottom", legend.box="horizontal")
ggsave(file.path(OUTDIR,"Fig_4panel_combo.pdf"), fig, width=7.2, height=3.0)                 # vector
ggsave(file.path(OUTDIR,"Fig_4panel_combo.svg"), fig, width=7.2, height=3.0)                 # vector (editable)
ggsave(file.path(OUTDIR,"Fig_4panel_combo.png"), fig, width=7.2, height=3.0, dpi=300)        # raster
cat("wrote Fig_4panel_combo.{pdf,svg,png}\n")
