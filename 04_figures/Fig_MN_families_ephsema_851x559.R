#!/usr/bin/env Rscript
# Combined 2x2: M = MSN<->Astro signalling families (left col); N = Ephrin/Semaphorin
# dot plots (right col). Full pipeline r=200, Robust_95>=0.5. Arial. Output 851x559 px.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg); library(stringr)})
FONT<-"Arial"; TS<-15; AX<-13
vf <- fread("/tmp/ef/validated_full.csv")[Scale==200 & Robust_95>=0.5]
raw<- fread("/tmp/ef/raw_ephsema_full.csv")
ZC <- c("Zone 1"="#e9776d","Zone 4"="#2bb3c4")

assign_family<-function(p){ if(is.na(p))return("Unknown")
  if(str_detect(p,"WNT"))return("Wnt"); if(str_detect(p,"NOTCH|JAG|DLL"))return("Notch")
  if(str_detect(p,"HH|SHH|IHH|DHH"))return("Hedgehog"); if(str_detect(p,"FGF"))return("FGF")
  if(str_detect(p,"EGF|TGFA|AREG|HBEGF|ERBB"))return("EGF"); if(str_detect(p,"PDGF"))return("PDGF")
  if(str_detect(p,"IGF|INS"))return("IGF/Insulin"); if(str_detect(p,"VEGF|PGF"))return("VEGF")
  if(str_detect(p,"BMP|GDF"))return("BMP"); if(str_detect(p,"TGF"))return("TGF-beta")
  if(str_detect(p,"NRG"))return("Neuregulins"); if(str_detect(p,"HGF|MET"))return("HGF")
  if(str_detect(p,"PTN|MDK"))return("PTN"); if(str_detect(p,"SEMA"))return("Semaphorins")
  if(str_detect(p,"EPH|EFN"))return("Ephrins"); if(str_detect(p,"SLIT|ROBO|SLITRK"))return("Slit/Robo")
  if(str_detect(p,"CNTN"))return("Contactins"); if(str_detect(p,"NTN|Netrin|UNC5|DCC|FLRT"))return("Netrin/UNC5")
  if(str_detect(p,"ADGR|C1QL|LPHN"))return("aGPCRs"); if(str_detect(p,"PTPR|LAR"))return("PTP Phosphatases")
  if(str_detect(p,"NECTIN|PVR|NECL"))return("Nectins"); if(str_detect(p,"CADM"))return("SynCAM")
  if(str_detect(p,"JAM|F11R"))return("JAM"); if(str_detect(p,"NCAM|L1CAM|NFASC"))return("CAMs (NCAM/L1)")
  if(str_detect(p,"CDH|PCDH"))return("Cadherins"); if(str_detect(p,"LAM|COL|FN1|THBS|TNC|SPP1|ITG"))return("ECM & Integrins")
  if(str_detect(p,"MIF|ACKR"))return("MIF"); if(str_detect(p,"APP|CD47|SIRP"))return("APP/CD47")
  if(str_detect(p,"GAS|PROS"))return("GAS6/PROS1"); if(str_detect(p,"MHC|HLA"))return("MHC")
  return(paste("Other:",p)) }
fam<-function(v) vapply(v,assign_family,character(1))
basetheme<-function() theme_classic(base_size=AX,base_family=FONT)+
  theme(axis.text.x=element_text(size=AX-1),axis.title.x=element_text(size=AX),
        plot.title=element_text(size=TS,hjust=.5),legend.position="bottom",
        legend.text=element_text(size=AX),legend.key.size=unit(10,"pt"),plot.margin=margin(3,5,2,2))

fampanel<-function(snd,rcv,title,drop_z4=character(0)){
  d<-vf[Sender==snd&Receiver==rcv]; d[,Family:=fam(Pathway)][,Bias:=fifelse(Mean_Real_Delta>0,"Zone 1","Zone 4")]
  cnt<-d[!grepl("^Other:|Unknown",Family),.(Z1=sum(Bias=="Zone 1"),Z4=sum(Bias=="Zone 4")),by=Family]
  z1<-cnt[Z1>Z4][order(-Z1)][1:4,Family]; z4<-cnt[Z4>Z1 & !(Family%in%drop_z4)][order(-Z4)][1:4,Family]
  sel<-cnt[Family%in%c(z1,z4)]; sel[,SM:=fifelse(Z1>=Z4,Z1,-Z4)]; ord<-sel[order(SM),Family]
  L<-rbind(sel[,.(Family,Bias="Zone 1",Count=Z1)],sel[,.(Family,Bias="Zone 4",Count=-Z4)]); L[,Family:=factor(Family,levels=ord)]
  ggplot(L,aes(Family,Count,fill=Bias))+geom_col(width=.66,colour="white",linewidth=.2)+
    geom_hline(yintercept=0,colour="black",linewidth=.3)+coord_flip()+
    scale_fill_manual(values=ZC,name=NULL)+scale_y_continuous(labels=abs,breaks=seq(-20,20,5))+
    labs(title=title,x=NULL,y="Interaction Count")+basetheme()+
    theme(axis.text.y=element_text(size=AX)) }

ve<-vf[Pathway %like% "EPH|SEMA"]
ephpanel<-function(snd,rcv,title){
  rob<-ve[Sender==snd&Receiver==rcv,Interaction]
  d<-raw[Sender==snd&Receiver==rcv&Interaction%in%rob,.(Donor,Interaction,Delta=Delta_ZA_ZB)]
  bars<-d[,.(mean=mean(Delta,na.rm=TRUE)),by=Interaction]; bars[,Bias:=fifelse(mean>0,"Zone 1","Zone 4")]
  bars[,fm:=fifelse(grepl("^SEMA",Interaction),"sema","ephrin")][,rk:=frank(-abs(mean),ties.method="first"),by=fm]
  keep<-bars[rk<=3,Interaction]; bars<-bars[Interaction%in%keep]; d<-d[Interaction%in%keep]
  ord<-bars[order(mean),Interaction]; d[,Interaction:=factor(Interaction,levels=ord)]; bars[,Interaction:=factor(Interaction,levels=ord)]
  ggplot()+geom_col(data=bars,aes(Interaction,mean,fill=Bias),width=.66,alpha=.9)+
    geom_jitter(data=d,aes(Interaction,Delta),width=.14,height=0,size=1.4,colour="#0d0d0d",alpha=.7,stroke=0)+
    geom_hline(yintercept=0,colour="black",linewidth=.3)+coord_flip()+
    scale_fill_manual(values=ZC,name=NULL)+
    labs(title=title,y="Enrichment (Delta Z1-Z4)",x=NULL)+basetheme()+
    theme(axis.text.y=element_text(size=AX-1,face="italic")) }

fMA<-fampanel("MSN","Astrocyte","MSN to Astrocyte Interaction Pathways",drop_z4="Netrin/UNC5")
fAM<-fampanel("Astrocyte","MSN","Astrocyte to MSN Interaction Pathways")
eMA<-ephpanel("MSN","Astrocyte","MSN to Astrocyte Ephrin/Semaphorin Signaling")
eAM<-ephpanel("Astrocyte","MSN","Astrocyte to MSN Ephrin/Semaphorin Signaling")

noleg <- theme(legend.position="none")
combined <- (fMA+noleg) + (eMA+noleg) + (fAM+noleg) + (eAM+noleg) +
  plot_layout(design="AB\nCD")
agg_png("/tmp/ef/Fig_MN_combined.png", width=851, height=559, units="px", res=72, background="white")
print(combined); dev.off()
svglite::svglite("/tmp/ef/Fig_MN_combined.svg", width=851/72, height=559/72); print(combined); dev.off()
# vector PDF with REAL editable Arial text (extrafont-registered ArialMT)
suppressMessages(library(extrafont)); loadfonts(device="pdf", quiet=TRUE)
pdf("/tmp/ef/Fig_MN_combined.pdf", width=851/72, height=559/72, useDingbats=FALSE)
print(combined)
dev.off()
cat("wrote Fig_MN_combined.{png,svg,pdf}; PDF = real editable Arial text (ArialMT)\n")
