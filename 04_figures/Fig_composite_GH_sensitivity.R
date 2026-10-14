#!/usr/bin/env Rscript
# Composite single row: G (custom HH/TGFb) | H (ligands from MSN) | sensitivity a-d.
# Arial, +50% fonts, NO bold, black axes+labels, 0.5pt strokes. base pdf() = 0 clip masks.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg)})
BLK<-"#000000"; INK<-BLK; MUT<-BLK; GRY<-BLK; FONT<-"Arial"
BS<-34; TT<-33; AXT<-30; AXL<-32; SUB<-27; LEG<-29; TAG<-42   # font sizes (pt) = prior x1.5, no bold
ZANN<-10.5                                                     # zone in-panel label size (7 x1.5)
LW<-0.5/2.1339567   # ggplot linewidth that renders as exactly 0.5 pt stroke on base pdf() (~0.234)
at <- function() theme_classic(base_size=BS, base_family=FONT, base_line_size=LW, base_rect_size=LW)
MICRO<-"#0E6B3C"    # microglia: dark forest green (old #27AE60 blended into white bar outlines)

## ===== PANEL G =====
shh <- fread("/tmp/ef/pg_shh.csv")[, .(Donor, Interaction, Delta=Delta_ZA_ZB)]
tg  <- fread("/tmp/ef/pg_tgfb2.csv")[, .(Donor, Interaction, Delta=Delta_ZA_ZB)]
zcol <- c("Zone 1"="#E9776D","Zone 4"="#2BB3C4")
gpanel <- function(d, lv){
  bars <- d[,.(mean=mean(Delta,na.rm=TRUE)),by=Interaction]; bars[, Zone := fifelse(mean>0,"Zone 1","Zone 4")]
  d[, Interaction:=factor(Interaction,levels=lv)]; bars[, Interaction:=factor(Interaction,levels=lv)]
  ggplot()+geom_col(data=bars,aes(Interaction,mean,fill=Zone),width=0.6,alpha=0.9)+
    geom_jitter(data=d,aes(Interaction,Delta),width=0.13,height=0,size=4.5,colour="#2a2a2a",alpha=0.7,stroke=0)+
    geom_hline(yintercept=0,colour="black",linewidth=LW)+coord_flip()+
    scale_fill_manual(values=zcol,drop=FALSE,name=NULL)+labs(x=NULL,y="Enrichment (Delta Z1-Z4)")+
    at()+theme(axis.text.y=element_text(face="italic",size=AXT,colour=BLK),axis.text.x=element_text(size=AXT,colour=BLK),
      axis.title.x=element_text(size=AXL,colour=BLK),legend.position="bottom",legend.text=element_text(size=LEG,colour=BLK),
      legend.key.size=unit(26,"pt"),plot.margin=margin(6,10,4,4)) }
pG <- gpanel(tg,"TGFB2_TGFBR3") / gpanel(shh,c("SHH_PTCH1","SHH_PTCH2")) +
      plot_layout(heights=c(1,2), guides="collect") & theme(legend.position="bottom")

## ===== PANEL H =====
vf <- fread("/tmp/ef/validated_full.csv")[Scale==200 & Robust_95>=0.5 & Sender=="MSN"]
vf[, zone := fifelse(Mean_Real_Delta>0,"Zone 1","Zone 4")]
famf <- function(lg) fifelse(grepl("^JAG|^DLL|^DLK",lg),"Notch",fifelse(grepl("^WNT",lg),"WNT",fifelse(grepl("^BMP",lg),"BMP",NA_character_)))
vf[, family:=famf(Ligand)]
cnt <- vf[!is.na(family),.(n=.N),by=.(family,Ligand,zone,Receiver)]; cnt[, plot_val:=fifelse(zone=="Zone 1",n,-n)]
rec_col <- c(Astrocyte="#3498DB",MSN="#E74C3C",OPC="#F1C40F",Microglia=MICRO,Oligo="#9B59B6",Vascular="#E67E22")
cnt[, Receiver:=factor(Receiver,levels=names(rec_col))]
hmk <- function(fm){
  d<-cnt[family==fm]; ord<-d[,.(net=sum(plot_val)),by=Ligand][order(net),Ligand]; d[,Ligand:=factor(Ligand,levels=ord)]
  Ntop<-length(ord); rng<-max(abs(d[,.(s=sum(plot_val)),by=.(Ligand,zone)]$s))
  dummy<-data.table(Ligand=factor(ord[1],levels=ord),plot_val=0,Receiver=factor(names(rec_col),levels=names(rec_col)))
  ggplot(d,aes(Ligand,plot_val,fill=Receiver))+geom_col(data=dummy,width=0)+geom_col(width=0.72,colour="white",linewidth=0.2)+
    annotate("segment",x=0.6,xend=0.6,y=-rng*1.10,yend=rng*1.10,colour="black",linewidth=LW)+          # bottom count axis (full)
    annotate("segment",x=0.6,xend=Ntop+0.45,y=-rng*1.10,yend=-rng*1.10,colour="black",linewidth=LW)+   # left axis stops at top bar
    annotate("segment",x=0.6,xend=Ntop+0.45,y=0,yend=0,colour="black",linewidth=LW)+                    # x=0 divider stops at top bar
    coord_flip(clip="off")+
    scale_fill_manual(values=rec_col,drop=FALSE,name="Receiver:",guide=guide_legend(nrow=3))+
    scale_y_continuous(labels=abs,limits=c(-rng*1.10,rng*1.10))+
    scale_x_discrete(expand=expansion(add=c(0.4,2.1)))+          # headroom above top bar for zone labels
    annotate("text",x=Ntop+1.35,y=-rng*0.5,label="Zone 4",colour="#2bb3c4",fontface=1,size=ZANN)+
    annotate("text",x=Ntop+1.35,y=rng*0.5,label="Zone 1",colour="#e9776d",fontface=1,size=ZANN)+
    labs(title=paste(fm,"Ligands from MSN"),x=NULL,y="Number of Interactions")+
    at()+theme(axis.text.y=element_text(face="italic",size=AXT,colour=BLK),axis.text.x=element_text(size=AXT,colour=BLK),
      axis.title.x=element_text(size=AXL,colour=BLK),plot.title=element_text(size=TT,hjust=0.5,colour=BLK),
      axis.line=element_blank(),                              # all frame lines removed; redrawn truncated as segments above
      legend.position="bottom",legend.text=element_text(size=LEG,colour=BLK),legend.title=element_text(size=LEG,colour=BLK),
      legend.key.size=unit(26,"pt"),plot.margin=margin(6,8,4,4)) }
pH <- (hmk("Notch") + theme(legend.position="none")) +
      (hmk("BMP") + theme(legend.position="none")) + hmk("WNT") + plot_layout(design="AC\nBC")

## ===== SENSITIVITY a/b/c/d =====
pal <- c(MSN="#E74C3C",Astrocyte="#3498DB",Oligo="#9B59B6",OPC="#F1C40F",Vascular="#E67E22",Microglia=MICRO)
ord <- c("MSN","Astrocyte","Oligo","OPC","Vascular","Microglia")
abbr <- c(MSN="MSN",Astrocyte="Astro",Oligo="Oligo",OPC="OPC",Vascular="Vasc",Microglia="Micro")
catcol <- c("same"="#9AA0A6","new"="#2E6E6A","lost"="#c94a4a")
rbr <- fread("/Users/akraft/PD_multiome/LR_sensitivity/FINAL/data/robustness_by_radius.csv"); rbr[, key:=paste(Sender,Receiver,Interaction,sep="||")]
kd <- fread("/Users/akraft/PD_multiome/LR_sensitivity/FINAL/data/kdist.csv")
mkchurn<-function(dt,s,a,b){aa<-dt[Sender==s&Scale==a&Robust_95>=0.5,key];bb<-dt[Sender==s&Scale==b&Robust_95>=0.5,key]
  data.table(same=length(intersect(aa,bb)),new=length(setdiff(bb,aa)),lost=length(setdiff(aa,bb)),n0=max(length(aa),length(bb)))}
churn_long<-function(a,bf){rbindlist(lapply(ord,function(s){d<-mkchurn(rbr,s,if(is.function(a))a(s) else a,if(is.function(bf))bf(s) else bf)
  data.table(Sender=s,cat=c("same","new","lost"),pct=100*c(d$same,d$new,d$lost)/d$n0)}))[,Sender:=factor(Sender,ord)][,cat:=factor(cat,c("same","new","lost"))]}
ra<-rbr[Scale%in%c(50,100,200,300,400,500,700,1000)&Robust_95>=0.5,.N,by=.(Sender,Scale)][Sender%in%ord]
ra[,pctmax:=100*N/max(N),by=Sender][,cell:=factor(Sender,levels=ord)]
Lb<-churn_long(200,400)
dk<-kd[k<=5,.(um=mean(med_dist)),by=.(receiver,k)][receiver%in%ord][,cell:=factor(receiver,levels=ord)]
dkr<-kd[,.(um=mean(med_dist)),by=.(receiver,k)]
R2<-sapply(ord,function(x)round(dkr[receiver==x&k==2,um]/100)*100);R5<-sapply(ord,function(x)round(dkr[receiver==x&k==5,um]/100)*100)
Ld<-churn_long(function(s)R2[[s]],function(s)R5[[s]])
thin<-at()+theme(axis.text=element_text(size=AXT,colour=BLK),axis.text.x=element_text(angle=40,hjust=1),
  axis.title=element_text(size=AXL,colour=BLK),axis.line=element_line(linewidth=LW,colour=BLK),axis.ticks=element_line(linewidth=LW,colour=BLK),
  plot.title=element_text(size=TT,colour=BLK),
  legend.position="bottom",legend.text=element_text(size=LEG,colour=BLK),legend.title=element_blank(),legend.key.size=unit(26,"pt"),plot.margin=margin(6,8,5,8))
pa<-ggplot(ra,aes(Scale,pctmax,colour=cell))+geom_line(linewidth=LW)+geom_point(size=2,stroke=0)+scale_colour_manual(values=pal,name=NULL)+
  scale_x_log10(breaks=c(50,100,200,400,1000))+scale_y_continuous(limits=c(0,103),breaks=c(0,50,100),labels=function(x)paste0(x,"%"))+
  labs(title="Robustness vs radius",x="radius (µm, log)",y="% of max")+thin
pb<-ggplot(Lb,aes(Sender,pct,fill=cat))+geom_col(position=position_dodge(width=0.8),width=0.74)+scale_fill_manual(values=catcol,name=NULL)+scale_x_discrete(labels=abbr)+
  scale_y_continuous(limits=c(0,104),breaks=c(0,50,100),labels=function(x)paste0(x,"%"),expand=expansion(mult=c(0,0.04)))+
  labs(title="200 vs 400 µm",x=NULL,y="% of max")+thin
pc<-ggplot(dk,aes(k,um,colour=cell))+geom_line(linewidth=LW)+geom_point(size=2.1,stroke=0)+scale_colour_manual(values=pal,name=NULL,guide="none")+
  scale_x_continuous(breaks=1:5)+scale_y_continuous(limits=c(0,NA),labels=comma)+
  labs(title="Neighbor distances",x="k-th partner",y="distance (µm)")+thin
pd<-ggplot(Ld,aes(Sender,pct,fill=cat))+geom_col(position=position_dodge(width=0.8),width=0.74)+scale_fill_manual(values=catcol,name=NULL)+scale_x_discrete(labels=abbr)+
  scale_y_continuous(limits=c(0,104),breaks=c(0,50,100),labels=function(x)paste0(x,"%"),expand=expansion(mult=c(0,0.04)))+
  labs(title="K=2 vs K=5 radius",x=NULL,y="% of max")+thin
sens <- pa + pb + pc + pd + plot_layout(design="AC\nBD", guides="collect") & theme(legend.position="bottom")

## ===== COMPOSE =====
composite <- (pG | sens | pH) + plot_layout(widths=c(1, 2.2, 2.6))
W<-48; H<-15
suppressMessages(library(extrafont)); loadfonts(device="pdf",quiet=TRUE)
agg_png("/tmp/ef/Fig_composite_arial.png", width=W, height=H, units="in", res=150, background="white")
print(composite); dev.off()
# PRIMARY editable PDF: base pdf() = 0 clip masks + editable Arial (referenced ArialMT)
pdf("/tmp/ef/Fig_composite_arial.pdf", width=W, height=H, family="Arial", useDingbats=FALSE)
print(composite); dev.off()
# SUBMISSION copy: quartz embeds Arial (adds ~150 clip masks -- do not edit)
grDevices::quartz(type="pdf", file="/tmp/ef/Fig_composite_arial_embedded.pdf", width=W, height=H)
print(composite); dev.off()
cat("wrote PNG + clean editable PDF (0 clips) + _embedded submission PDF\n")
