#!/usr/bin/env Rscript
# Panel G — enrichment (Delta Z1-Z4) of the custom HH / TGFb interactions, full
# pipeline r=200, Robust_95>=0.5. Bar = mean donor enrichment; dots = per-donor.
# TGFB2_TGFBR3 (Astro->Astro, Zone 1) replaces the TGFB2_TGFBR1/2_TGFBR3 complex;
# SHH_PTCH2 / SHH_PTCH1 (MSN->Astro, Zone 4). Two panels = two enrichment scales.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
shh <- fread("/tmp/ef/pg_shh.csv")[, .(Donor, Interaction, Delta=Delta_ZA_ZB)]
tg  <- fread("/tmp/ef/pg_tgfb2.csv")[, .(Donor, Interaction, Delta=Delta_ZA_ZB)]
zcol <- c("Zone 1"="#E9776D","Zone 4"="#2BB3C4")

panel <- function(d, order_levels){
  bars <- d[,.(mean=mean(Delta,na.rm=TRUE)),by=Interaction]
  bars[, Zone := fifelse(mean>0,"Zone 1","Zone 4")]
  d[, Interaction := factor(Interaction, levels=order_levels)]
  bars[, Interaction := factor(Interaction, levels=order_levels)]
  ggplot() +
    geom_col(data=bars, aes(Interaction, mean, fill=Zone), width=0.6, alpha=0.9) +
    geom_jitter(data=d, aes(Interaction, Delta), width=0.14, height=0, size=2,
                colour="#2a2a2a", alpha=0.7, stroke=0) +
    geom_hline(yintercept=0, colour="black", linewidth=0.35) +
    coord_flip() +
    scale_fill_manual(values=zcol, drop=FALSE, name=NULL) +
    labs(x=NULL, y="Enrichment (Delta Z1-Z4)") +
    theme_classic(base_size=11) +
    theme(axis.text.y=element_text(face="italic", size=10.5, colour="#14171C"),
          axis.text.x=element_text(colour="#6A7178"),
          legend.position="bottom",
          plot.margin=margin(6,10,4,4))
}
pT <- panel(tg,  "TGFB2_TGFBR3")
pS <- panel(shh, c("SHH_PTCH1","SHH_PTCH2"))

fig <- (pT / pS) + plot_layout(heights=c(1,2), guides="collect") &
       theme(legend.position="bottom")
ggsave("/tmp/ef/Fig_panelG_custom.pdf", fig, width=4.6, height=3.8)
ggsave("/tmp/ef/Fig_panelG_custom.png", fig, width=4.6, height=3.8, dpi=200)
ggsave("/tmp/ef/Fig_panelG_custom.svg", fig, width=4.6, height=3.8)
cat("wrote Fig_panelG_custom.{pdf,png,svg}\n")
