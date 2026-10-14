#!/usr/bin/env Rscript
# Panel H — Notch / WNT / BMP ligands FROM MSN, full pipeline r=200, Robust_95>=0.5.
# Per ligand: number of robust interactions, diverging by zone (Zone 4 left / Zone 1
# right), stacked & coloured by RECEIVER cell type.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
vf <- fread("/tmp/ef/validated_full.csv")[Scale==200 & Robust_95>=0.5 & Sender=="MSN"]
vf[, zone := fifelse(Mean_Real_Delta>0,"Zone 1","Zone 4")]
fam <- function(lg) fifelse(grepl("^JAG|^DLL|^DLK",lg),"Notch",
                     fifelse(grepl("^WNT",lg),"WNT", fifelse(grepl("^BMP",lg),"BMP",NA_character_)))
vf[, family := fam(Ligand)]
cnt <- vf[!is.na(family), .(n=.N), by=.(family,Ligand,zone,Receiver)]
cnt[, plot_val := fifelse(zone=="Zone 1", n, -n)]

rec_col <- c(Astrocyte="#3498DB", MSN="#E74C3C", OPC="#F1C40F",
             Microglia="#27AE60", Oligo="#9B59B6", Vascular="#E67E22")
cnt[, Receiver := factor(Receiver, levels=names(rec_col))]

INK<-"#14171C"; MUT<-"#6A7178"
mk <- function(fm){
  d <- cnt[family==fm]
  ord <- d[, .(net=sum(plot_val)), by=Ligand][order(net), Ligand]   # Zone1-dominant on top (after coord_flip)
  d[, Ligand := factor(Ligand, levels=ord)]
  rng <- max(abs(d[, .(s=sum(plot_val)), by=.(Ligand,zone)]$s))
  ggplot(d, aes(Ligand, plot_val, fill=Receiver)) +
    geom_col(width=0.72, colour="white", linewidth=0.15) +
    geom_hline(yintercept=0, colour="black", linewidth=0.4) +
    coord_flip() +
    scale_fill_manual(values=rec_col, drop=FALSE, name="Receiver:") +
    scale_y_continuous(labels=abs, limits=c(-rng*1.08, rng*1.08)) +
    annotate("text", x=Inf, y=-rng*0.55, label="Zone 4", colour="#2bb3c4", fontface=2, size=3.4, vjust=1.6) +
    annotate("text", x=Inf, y= rng*0.55, label="Zone 1", colour="#e9776d", fontface=2, size=3.4, vjust=1.6) +
    labs(title=paste(fm,"Ligands from MSN"), x=NULL, y="Number of Interactions") +
    theme_classic(base_size=11) +
    theme(axis.text.y=element_text(face="italic", size=9.5, colour=INK),
          axis.text.x=element_text(colour=MUT),
          plot.title=element_text(face="bold", size=11, hjust=0.5),
          legend.position="bottom",
          plot.margin=margin(10,8,4,4))
}
pN <- mk("Notch"); pW <- mk("WNT"); pB <- mk("BMP")
fig <- (pN / pB | pW) + plot_layout(guides="collect", widths=c(1,1)) &
       theme(legend.position="bottom")
ggsave("/tmp/ef/Fig_panelH_ligands.pdf", fig, width=10.5, height=7.2)
ggsave("/tmp/ef/Fig_panelH_ligands.png", fig, width=10.5, height=7.2, dpi=200)
ggsave("/tmp/ef/Fig_panelH_ligands.svg", fig, width=10.5, height=7.2)
cat("wrote Fig_panelH_ligands.{pdf,png,svg}\n")
