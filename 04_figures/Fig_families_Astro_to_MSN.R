#!/usr/bin/env Rscript
# Astrocyte -> MSN pathway-family zonal landscape, FULL pipeline, r=200, Rob95>=0.5.
# Top 4 Zone-1 families + top 4 Zone-4 families. Diverging bars: Zone 1 Patch (red,
# right) / Zone 4 Matrix (teal, left). Matches Fig_msn_astro_top4 style (ticks at 5).
suppressPackageStartupMessages({library(data.table); library(stringr); library(ggplot2)})
ROB<-0.5
assign_family <- function(p){
  if(is.na(p)) return("Unknown")
  if(str_detect(p,"WNT")) return("Wnt"); if(str_detect(p,"NOTCH|JAG|DLL")) return("Notch")
  if(str_detect(p,"HH|SHH|IHH|DHH")) return("Hedgehog"); if(str_detect(p,"FGF")) return("FGF")
  if(str_detect(p,"EGF|TGFA|AREG|HBEGF|ERBB")) return("EGF"); if(str_detect(p,"PDGF")) return("PDGF")
  if(str_detect(p,"IGF|INS")) return("IGF/Insulin"); if(str_detect(p,"VEGF|PGF")) return("VEGF")
  if(str_detect(p,"BMP|GDF")) return("BMP"); if(str_detect(p,"TGF")) return("TGF-beta")
  if(str_detect(p,"NRG")) return("Neuregulins"); if(str_detect(p,"HGF|MET")) return("HGF")
  if(str_detect(p,"PTN|MDK")) return("PTN"); if(str_detect(p,"SEMA")) return("Semaphorins")
  if(str_detect(p,"EPH|EFN")) return("Ephrins"); if(str_detect(p,"SLIT|ROBO|SLITRK")) return("Slit/Robo")
  if(str_detect(p,"CNTN")) return("Contactins"); if(str_detect(p,"NTN|Netrin|UNC5|DCC|FLRT")) return("Netrin/UNC5")
  if(str_detect(p,"ADGR|C1QL|LPHN")) return("aGPCRs"); if(str_detect(p,"PTPR|LAR")) return("PTP Phosphatases")
  if(str_detect(p,"NECTIN|PVR|NECL")) return("Nectins"); if(str_detect(p,"CADM")) return("SynCAM")
  if(str_detect(p,"JAM|F11R")) return("JAM"); if(str_detect(p,"NCAM|L1CAM|NFASC")) return("CAMs (NCAM/L1)")
  if(str_detect(p,"CDH|PCDH")) return("Cadherins"); if(str_detect(p,"LAM|COL|FN1|THBS|TNC|SPP1|ITG")) return("ECM & Integrins")
  if(str_detect(p,"MIF|ACKR")) return("MIF"); if(str_detect(p,"APP|CD47|SIRP")) return("APP/CD47")
  if(str_detect(p,"GAS|PROS")) return("GAS6/PROS1"); if(str_detect(p,"MHC|HLA")) return("MHC")
  return(paste("Other:",p)) }
fam<-function(v) vapply(v, assign_family, character(1))

vf<-fread("/tmp/ef/validated_full.csv")[Scale==200 & Sender=="Astrocyte" & Receiver=="MSN" & Robust_95>=ROB]
vf[,Family:=fam(Pathway)][,Bias:=fifelse(Mean_Real_Delta>0,"Zone 1 (Patch)","Zone 4 (Matrix)")]
cnt<-vf[!grepl("^Other:|Unknown",Family),.(Z1=sum(Bias=="Zone 1 (Patch)"),Z4=sum(Bias=="Zone 4 (Matrix)")),by=Family]

z1<-cnt[Z1>Z4][order(-Z1)][1:4, Family]
z4<-cnt[Z4>Z1][order(-Z4)][1:4, Family]
sel<-cnt[Family %in% c(z1,z4)]
cat("Zone-1 top4:",paste(z1,collapse=", "),"\nZone-4 top4:",paste(z4,collapse=", "),"\n\n"); print(sel)

sel[,SortMetric:=fifelse(Z1>=Z4, Z1, -Z4)]
ord<-sel[order(SortMetric),Family]
L<-rbind(sel[,.(Family,Bias="Zone 1 (Patch)",Count=Z1)], sel[,.(Family,Bias="Zone 4 (Matrix)",Count=-Z4)])
L[,Family:=factor(Family,levels=ord)]

INK<-"#14171C"; MUT<-"#6A7178"
p<-ggplot(L,aes(Family,Count,fill=Bias))+
  geom_col(width=.68,colour="white",linewidth=.25)+coord_flip()+
  geom_hline(yintercept=0,colour="black",linewidth=.35)+
  scale_fill_manual(values=c("Zone 1 (Patch)"="#e9776d","Zone 4 (Matrix)"="#2bb3c4"))+
  scale_y_continuous(labels=abs, breaks=seq(-20,20,5))+
  labs(title="Astrocyte → MSN signaling families",
       subtitle="full pipeline, 200 µm · robust (Rob95≥0.5) · top 4 families per zone",
       x=NULL,y="robust interactions",fill=NULL)+
  theme_minimal(base_size=12)+
  theme(legend.position="bottom",panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_text(face="bold",size=11,colour=INK),
        axis.text.x=element_text(colour=MUT),
        plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=8.5,colour=MUT))
ggsave("/tmp/ef/Fig_astro_msn_top4.pdf",p,width=4.95,height=4.6)
ggsave("/tmp/ef/Fig_astro_msn_top4.png",p,width=4.95,height=4.6,dpi=200)
cat("\nwrote Fig_astro_msn_top4.{pdf,png}\n")
