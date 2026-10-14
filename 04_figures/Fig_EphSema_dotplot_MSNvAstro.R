#!/usr/bin/env Rscript
# Ephrin/Semaphorin guidance-cue dot plot, FULL pipeline, r=200, Rob95>=0.5.
# Bar = mean per-donor enrichment (Delta Z1-Z4); dots = individual donors.
# Zone 1 (Patch) red, Zone 4 (Matrix) teal. Reproduces the original targeted plot.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROB<-0.5
raw<-fread("/tmp/ef/raw_ephsema_full.csv")                       # per-donor Delta_ZA_ZB
vf <-fread("/tmp/ef/validated_full.csv")[Scale==200 & Robust_95>=ROB & Pathway %like% "EPH|SEMA"]

chan<-function(snd,rcv){
  rob<-vf[Sender==snd & Receiver==rcv, .(Interaction, Robust_95)]
  d<-raw[Sender==snd & Receiver==rcv & Interaction %in% rob$Interaction,
         .(Donor,Interaction,Delta=Delta_ZA_ZB)]
  d[,Sender:=snd][,Receiver:=rcv]
  d
}
D<-rbind(chan("MSN","Astrocyte"), chan("Astrocyte","MSN"))
D[,chan:=factor(fifelse(Sender=="MSN","MSN → Astrocyte","Astrocyte → MSN"),
                levels=c("MSN → Astrocyte","Astrocyte → MSN"))]
bars<-D[,.(mean=mean(Delta,na.rm=TRUE)),by=.(chan,Interaction)]
bars[,Bias:=fifelse(mean>0,"Zone 1 (Patch)","Zone 4 (Matrix)")]

# keep the top-3 (by |mean| enrichment) EPHRIN and top-3 SEMAPHORIN interactions per channel
bars[,fam:=fifelse(grepl("^SEMA",Interaction),"sema","ephrin")]
bars[,rk:=frank(-abs(mean),ties.method="first"),by=.(chan,fam)]
keep <- bars[rk<=3, paste(chan,Interaction)]
bars <- bars[paste(chan,Interaction) %in% keep][,c("fam","rk"):=NULL]
D    <- D[paste(chan,Interaction) %in% keep]
# order interactions within each panel by mean (ascending -> Zone1 on top after coord_flip)
bars<-bars[order(chan,mean)]
bars[,ikey:=paste(chan,Interaction,sep="__")]        # unique key (an interaction can be robust in both channels)
D<-merge(D,bars[,.(chan,Interaction,Bias,mean,ikey)],by=c("chan","Interaction"))
lev<-bars$ikey
D[,ikey:=factor(ikey,levels=lev)]; bars[,ikey:=factor(ikey,levels=lev)]
ilab<-setNames(as.character(bars$Interaction),as.character(bars$ikey))

p<-ggplot()+
  geom_col(data=bars,aes(ikey,mean,fill=Bias),width=.72,alpha=.9)+
  geom_jitter(data=D,aes(ikey,Delta),width=.16,height=0,size=2.1,
              colour="#0d0d0d",alpha=.7,stroke=0)+
  geom_hline(yintercept=0,colour="black",linewidth=.3)+
  coord_flip()+
  facet_grid(chan~.,scales="free_y",space="free_y",switch="y")+
  scale_x_discrete(labels=ilab)+
  scale_fill_manual(values=c("Zone 1 (Patch)"="#f8766d","Zone 4 (Matrix)"="#00bfc4"))+
  labs(title="Ephrin / Semaphorin guidance cues — full pipeline, 200 µm",
       subtitle="bar = mean donor enrichment; dots = individual donors; robust (Rob95≥0.5)",
       y="Enrichment (Delta Z1-Z4)",x=NULL,fill=NULL)+
  theme_minimal(base_size=10)+
  theme(legend.position="bottom",panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=7.5,face="italic"),
        strip.text.y.left=element_text(angle=0,face="bold",size=9.5),
        strip.placement="outside",
        plot.title=element_text(face="bold",size=11.5))
ggsave("/tmp/ef/Fig_ephsema_dotplot.pdf",p,width=8.2,height=8.8)
ggsave("/tmp/ef/Fig_ephsema_dotplot.png",p,width=8.2,height=8.8,dpi=200)
cat("wrote Fig_ephsema_dotplot.{pdf,png}\n")
cat("\n=== robust eph/sema, mean enrichment (full pipeline, r=200) ===\n")
print(bars[order(chan,-mean),.(chan,Interaction,mean=round(mean,3),Bias)])
