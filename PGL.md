---
Title: Visualisation of differential gene expression
Author: Stanislav Vosolsobě
---

# VISUALISATION OF DIFFERENTIAL GENE EXPRESSION

Published in: **POLYGALACTURONASES REGULATED BY AUXIN facilitate root cell elongation in Arabidopsis thaliana via pectin remodelling**

Monika Kubalová, Anna Kampová, Stanislav Vosolsobě, Karel Raabe, Brigita Simonaviciene, Yoselin Benitez-Alfonso, Karel Müller, Eva Medvecká, Matyáš Fendrych

doi: https://doi.org/10.1101/2025.05.07.652666 

Requirements

- *R* with packages **scales** and **shades**
- *Bash* console
- *Inkscape*

```r
require(scales)
require(shades)
```

## Loading of input data


```r
setwd("/media/Home/home/standa/Plocha/Mates/Monika/PGL/")

pgl <- read.table("data_AXR3-1_genelevels_simple.csv", header = T, sep="\t",quote = "",comment.char = "")
pll <- read.table("pll.csv", header = T, sep="\t",quote = "",comment.char = "")
PLL7 <- read.table("PLL7.csv", header = T, sep="\t",quote = "",comment.char = "")
```

## Initilisation of variables

sifni - list of significant genes
qval  - respective q-values
fold  - respective fold-changes

```r
signi <- list()
qval <- list()
fold <- list()
```

## Acquisition of DEGs from original dataset

###Selection of the list of DEGs

```r
signi[["Col_noS_mock"]] <- pgl$gene_ID[pgl$Col_noS_mock.sig!="-"]
signi[["PLL10_noS_mock"]] <- pgl$gene_ID[pgl$PLL10_noS_mock.sig!="-"]
signi[["mock_PLL10_Col"]] <- pgl$gene_ID[pgl$mock_PLL10_Col.sig!="-"]
signi[["mock_PLL7_Col"]] <- PLL7$gene_ID[PLL7$mock.PLL7.vs.Col.0.II!="-"]
signi[["pll_Col"]] <- pll$gene_ID[pll$pII_Col.sign!="-"]
```

###Selection of the statistical parameters for DEGs

q-values

```r
qval[["Col_noS_mock"]] <- pgl$Col_noS_mock.qval[pgl$Col_noS_mock.sig!="-"]
qval[["PLL10_noS_mock"]] <- pgl$PLL10_noS_mock.qval[pgl$PLL10_noS_mock.sig!="-"]
qval[["mock_PLL10_Col"]] <- pgl$mock_PLL10_Col.qval[pgl$mock_PLL10_Col.sig!="-"]
qval[["mock_PLL7_Col"]] <- PLL7$qval[PLL7$mock.PLL7.vs.Col.0.II!="-"]
qval[["pll_Col"]] <- pll$pII_Col.qval[pll$pII_Col.sign!="-"]

names(qval[["Col_noS_mock"]]) <- signi[["Col_noS_mock"]]
names(qval[["PLL10_noS_mock"]]) <- signi[["PLL10_noS_mock"]]
names(qval[["mock_PLL10_Col"]]) <- signi[["mock_PLL10_Col"]]
names(qval[["mock_PLL7_Col"]]) <- signi[["mock_PLL7_Col"]]
names(qval[["pll_Col"]]) <- signi[["pll_Col"]]
```

Fold-change

```r
fold[["Col_noS_mock"]] <- pgl$Col_noS_mock.b[pgl$Col_noS_mock.sig!="-"]
fold[["PLL10_noS_mock"]] <- pgl$PLL10_noS_mock.b[pgl$PLL10_noS_mock.sig!="-"]
fold[["mock_PLL10_Col"]] <- pgl$mock_PLL10_Col.b[pgl$mock_PLL10_Col.sig!="-"]
fold[["mock_PLL7_Col"]] <- PLL7$b[PLL7$mock.PLL7.vs.Col.0.II!="-"]
fold[["pll_Col"]] <- pll$pII_Col.b[pll$pII_Col.sign!="-"]

names(fold[["Col_noS_mock"]]) <- signi[["Col_noS_mock"]]
names(fold[["PLL10_noS_mock"]]) <- signi[["PLL10_noS_mock"]]
names(fold[["mock_PLL10_Col"]]) <- signi[["mock_PLL10_Col"]]
names(fold[["mock_PLL7_Col"]]) <- signi[["mock_PLL7_Col"]]
names(fold[["pll_Col"]]) <- signi[["pll_Col"]]
```

Merging the data into a final dataset

```r
finframe_col_noS <- data.frame(gene=signi[["Col_noS_mock"]],qval=qval[["Col_noS_mock"]],fold=fold[["Col_noS_mock"]])
finframe_PLL10_noS <- data.frame(gene=signi[["PLL10_noS_mock"]],qval=qval[["PLL10_noS_mock"]],fold=fold[["PLL10_noS_mock"]])
finframe_PLL10_col <- data.frame(gene=signi[["mock_PLL10_Col"]],qval=qval[["mock_PLL10_Col"]],fold=fold[["mock_PLL10_Col"]])
finframe_PLL7_col <- data.frame(gene=signi[["mock_PLL7_Col"]],qval=qval[["mock_PLL7_Col"]],fold=fold[["mock_PLL7_Col"]])
finframe_pll_col <- data.frame(gene=signi[["pll_Col"]],qval=qval[["pll_Col"]],fold=fold[["pll_Col"]])


finframe_PLLall <- data.frame(gene=unique(unlist(signi[c("mock_PLL7_Col","mock_PLL10_Col")])),qval=NA,PLL7=NA,PLL10=NA)


for(j in 1:nrow(finframe_PLLall)){
  g <- finframe_PLLall$gene[j]
  finframe_PLLall$qval[j] <- min(unlist(lapply(qval[c("mock_PLL7_Col","mock_PLL10_Col")],`[`,g)),na.rm = T)
  ff <- unlist(lapply(fold[c("mock_PLL7_Col")],`[`,g))
  fm <- ff[which.max(abs(ff))]
  if(length(fm > 0)){
    finframe_PLLall$PLL7[j] <- fm
  }
  ff <- unlist(lapply(fold[c("mock_PLL10_Col")],`[`,g))
  fm <- ff[which.max(abs(ff))]
  if(length(fm > 0)){
    finframe_PLLall$PLL10[j] <- ff[which.max(abs(ff))]
  }
}
finframe_PLLall$PLLm <- rowMeans(finframe_PLLall[,c("PLL7","PLL10")],na.rm = T)

plot(finframe_PLLall$PLL10~finframe_PLLall$PLL7)
abline(b=1,a = 0)

finframe_PLL10_col$qval[finframe_PLL10_col$qval<1e-100] <- 1e-100
finframe_PLL7_col$qval[finframe_PLL7_col$qval<1e-100] <- 1e-100
finframe_PLLall$qval[finframe_PLLall$qval<1e-100] <- 1e-100



write.table(finframe_col_noS,"finframe_col_noS.tab",quote = F,row.names = F,sep = "\t")
write.table(finframe_PLL10_noS,"finframe_PLL10_noS.tab",quote = F,row.names = F,sep = "\t")
write.table(finframe_PLL10_col,"finframe_PLL10_col.tab",quote = F,row.names = F,sep = "\t")
write.table(finframe_PLL7_col,"finframe_PLL7_col.tab",quote = F,row.names = F,sep = "\t")
write.table(finframe_pll_col,"finframe_pll_col.tab",quote = F,row.names = F,sep = "\t")
write.table(finframe_PLLall,"finframe_PLLall.tab",quote = F,row.names = F,sep = "\t")
```

Now upload the table of genes to STRING, download the map file and SVG with the network

### Formating table for network

```r
mapp2fmt <- function(map,frame,fmt_out,mult=10){
  mapp <- read.delim(map,comment.char = "")  # There are ' chars in names and annotations!
  rng <- ceiling(max(abs(c(range(frame$fold, na.rm = T),range(frame$fold, na.rm = T)))))
  #rampcols <- colorRampPalette(colors = c("green", "white", "magenta"), space = "rgb")(rng*20+1)
  #rampcols <- colorRampPalette(colors = c("green", "magenta"), space = "rgb")(rng*20+1)
  rampcols <- colorRampPalette(colors = c("green","wheat", "magenta"), space = "rgb")(rng*20+1)
  #rc1 <- colorRampPalette(colors = c("green", "wheat"), space = "rgb")(20)
  #rc2 <- colorRampPalette(colors = c("wheat", "magenta"), space = "rgb")(20)
  #rampcols <- c(rc1, rc2)
  require(plotfunctions)
  pdf(file = paste(fmt_out,"pdf",sep = "."),width=1,height = 3)
  par(mar=c(0,0,0,0))
  plot(1,type="n",xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  gradientLegend(valRange = c(-rng,rng),rampcols, dec = 1,pos = c(-1,-1,0,1),coords = T,n.seg = ((-rng/2):(rng/2))*2)
  dev.off()
  swatch(rampcols)
  frame_col <- round(frame$fold*10+rng*10+1)
  frame_col[is.na(frame_col)] <- rng*10+1
  mapname <- NULL
  agi <- NULL
  mixcol <- NULL
    for(c in 1:nrow(frame)){
    mixcol[c] <-rampcols[frame_col[c]]
    mn <- mapp$preferredName[mapp$queryItem==frame$gene[c]]
    if(length(mn>0)){
      mapname[c] <- mn
      agi[c] <- mapp$queryItem[mapp$queryItem==frame$gene[c]] 
    }
  }
  diam <- mult*sqrt(-log10(frame$qval))
  fmt <- data.frame(mapname,agi,mixcol,diam)
  write.table(fmt,fmt_out,quote = F,sep = "\t",row.names = F,col.names = F)
}
```

### Two colours only NEW FUNCTION
```r
mapp2fmt <- function(map,frame,fmt_out,mult){
  require(scales)
  require(shades)
  require(interp)
  mapp <- read.delim(map,comment.char = "")  # There are ' chars in names and annotations!
  
  
  na="gray"
  up="magenta"
  do="#07ff07"
  cols <- c(up,na,do)
  colm <- matrix(cols,3,1)
  
  
  
  pdf(file = paste(fmt_out,"pdf",sep = "."),width=4,height = 12)
  par(mar=c(4,4,2,1),mfcol=c(3,1))
  plot(1,type="n",xlim=c(0,3),ylim=c(0,1),xlab="Expression",ylab="",xaxt="n",yaxt="n",bty="n") #
  for(a in 1:3){
    m=1
    rect(a-1,m-1,a,m,col=colm[a,m],border = NA)
  }
  axis(side = 1, at = (1:3)-0.5, labels = c("up","NA","down"),pos=0,lwd=0)
  
  plot(1,type="n",xlim=c(0,3),ylim=c(0,3),xlab="Log2FC",ylab="",xaxt="n",yaxt="n",bty="n") #
  circles(x=rep(1.5,6),y=(1:6)/5,r=(1:6)/5)
  text(x = rep(1.5,6),y=2*(1:6)/5-1/5,labels = 1:6)
  
  plot(1,type="n",xlim=c(1,100),ylim=c(0,1),xlab="-log10(FDR)",ylab="Colour density")
  for(i in 1:100){
    pcol <- submix("white","purple",amount = sqrt(-log10(10^(-i)))/10)
    points(i,sqrt(-log10(10^(-i)))/10,pch=16, col=pcol )
  }
  
  dev.off()
  
  
  diam <- NULL
  mapname <- NULL
  agi <- NULL
  mixcol <- NULL
  
  for(c in 1:nrow(frame)){
    if(is.na(frame$fold[c])){
      mixcol[c]<-na
    }else{
      if(frame$fold[c]>0){  
        mixcol[c]<-up
      }else{
        mixcol[c]<-do
      }
    }
    mixcol[c] <- submix("white",mixcol[c],amount = sqrt(-log10(frame$qval[c]))/10)
    
    mn <- mapp$preferredName[mapp$queryItem==frame$gene[c]]
    if(length(mn>0)){
      mapname[c] <- mn
      agi[c] <- mapp$queryItem[mapp$queryItem==frame$gene[c]] 
    }
    diam[c] <- mult*abs(frame$fold[c])
  }
  
  outtable <- data.frame(mapname,agi,mixcol,diam)
  
  write.table(outtable,fmt_out,quote = F,sep = "\t",row.names = F,col.names = F)
  
}
```



### Double scale
```r
mapp3fmt <- function(map,frame,fmt_out,mult){
  require(scales)
  require(shades)
  require(interp)
  mapp <- read.delim(map,comment.char = "")  # There are ' chars in names and annotations!
  
  
  PLL7naPLL10na="gray"
  PLL7naPLL10do="#faf317"
  PLL7naPLL10up="#00c5ff"
  
  PLL7upPLL10na="magenta"
  PLL7upPLL10do="#ff8d00"
  PLL7upPLL10up="purple"
  
  PLL7doPLL10na="#008000"
  PLL7doPLL10do="#07ff07"
  PLL7doPLL10up="#00a09f"
  
  cols <- c(PLL7upPLL10up,PLL7naPLL10up,PLL7doPLL10up,
            PLL7upPLL10na,PLL7naPLL10na,PLL7doPLL10na,
            PLL7upPLL10do,PLL7naPLL10do,PLL7doPLL10do)
  colm <- matrix(cols,3,3)
  
  
  
  pdf(file = paste(fmt_out,"pdf",sep = "."),width=4,height = 12)
  par(mar=c(4,4,2,1),mfcol=c(3,1))
  plot(1,type="n",xlim=c(0,3),ylim=c(0,3),xlab="PLL7",ylab="PLL10",xaxt="n",yaxt="n",bty="n") #
  for(a in 1:3){
    for(m in 1:3){
      rect(a-1,m-1,a,m,col=colm[a,m],border = NA)
    }
  }
  axis(side = 1, at = (1:3)-0.5, labels = c("up","NA","down"),pos=0,lwd=0)
  axis(side = 2, at = (1:3)-0.5, labels = c("up","NA","down"),pos=0,lwd=0)
  plot(1,type="n",xlim=c(0,3),ylim=c(0,3),xlab="Log2FC",ylab="",xaxt="n",yaxt="n",bty="n") #
  circles(x=rep(1.5,6),y=(1:6)/5,r=(1:6)/5)
  text(x = rep(1.5,6),y=2*(1:6)/5-1/5,labels = 1:6)
  
  plot(1,type="n",xlim=c(1,100),ylim=c(0,1),xlab="-log10(FDR)",ylab="Colour density")
  for(i in 1:100){
    pcol <- submix("white","purple",amount = sqrt(-log10(10^(-i)))/10)
    points(i,sqrt(-log10(10^(-i)))/10,pch=16, col=pcol )
  }
  
  dev.off()
  
  
  diam <- NULL
  mapname <- NULL
  agi <- NULL
  mixcol <- NULL
  
  for(c in 1:nrow(frame)){
    if(is.na(frame$PLL7[c])){
      if(is.na(frame$PLL10[c])){mixcol[c]<-PLL7naPLL10na}else{
        if(frame$PLL10[c]>0){mixcol[c]<-PLL7naPLL10up}else{mixcol[c]<-PLL7naPLL10do}
      }
    }else{
      if(frame$PLL7[c]>0){
        if(is.na(frame$PLL10[c])){mixcol[c]<-PLL7upPLL10na}else{
          if(frame$PLL10[c]>0){mixcol[c]<-PLL7upPLL10up}else{mixcol[c]<-PLL7upPLL10do}
        }
        
      }else{
        if(is.na(frame$PLL10[c])){mixcol[c]<-PLL7doPLL10na}else{
          if(frame$PLL10[c]>0){mixcol[c]<-PLL7doPLL10up}else{mixcol[c]<-PLL7doPLL10do}
        }
      }
    }
    mixcol[c] <- submix("white",mixcol[c],amount = sqrt(-log10(frame$qval[c]))/10)
    
    mn <- mapp$preferredName[mapp$queryItem==frame$gene[c]]
    if(length(mn>0)){
      mapname[c] <- mn
      agi[c] <- mapp$queryItem[mapp$queryItem==frame$gene[c]] 
    }
    diam[c] <- mult*abs(frame$PLLm[c])
  }
  
  outtable <- data.frame(mapname,agi,mixcol,diam)
  
  write.table(outtable,fmt_out,quote = F,sep = "\t",row.names = F,col.names = F)
  
}
```




TSV file - downloaded mapping file from string
Frame - gene table
FMT file - output file

```r
mapp2fmt("Col_noS_mock.tsv",finframe_col_noS,"fmt_col_noS.fmt")
mapp2fmt("PLL10_noS_mock.tsv",finframe_PLL10_noS,"fmt_PLL10_noS.fmt")
mapp2fmt("pll_Col.tsv",finframe_pll_col,"fmt_pll_col.fmt")

mapp2fmt("mock_PLL10_Col.tsv",frame = finframe_PLL10_col,"fmt_PLL10_col.fmt",mult = 20)
mapp2fmt("mock_PLL7_Col.tsv",frame = finframe_PLL7_col,"fmt_PLL7_col.fmt",mult = 20)


mapp3fmt(map = "PLLall.tsv",frame = finframe_PLLall,fmt_out = "PLLall.fmt",mult=35)
```


##


```sh
cd /media/Home/home/standa/Plocha/Mates/Monika/PGL

  # File 'Anthers_005-STRING-100' was manually edited to move distant nodes inside main frame >> Anthers_005-STRING-100u
filename="Col_noS_mock" # without .svg, 
filename="PLL10_noS_mock" # without .svg,
filename="mock_PLL10_Col_core" # without .svg,
filename="mock_PLL10_Col"
filename="mock_PLL7_Col_core" # without .svg,
filename="mock_PLL7_Col"
filename="pll_Col"

#filename="PLL_all_notext_core"
filename="PLL_all_notext_full"
filename="PLL7_notext_full"
filename="PLL10_notext_full"

inkscape --export-filename=${filename}u.svg --export-type="svg" ${filename}.svg

cat ${filename}u.svg | tr -d $'\n' | sed -e "s/>/>\n/g" | sed -e "s/font-size: 12px/font-size: 20px/" -e "s/url(#filter_bg_textFlat)/url(#filter_bg_text)/" -e "s/style=\"fill:url(#radialGradient.*\"//" -e "s/opacity: 0.4/opacity: 1/" > ${filename}u_formated.svg
while read l
do
echo $l
name=`echo $l | cut -f 1 -d " "`
colfos=`echo $l | cut -f 3 -d " "`
radfos=`echo $l | cut -f 4 -d " "`
if [ "$colfos" != "NA" ]
then
sed -i "/data-safe_div_label=\"${name}\"/,/<\/g>/ s/r=\"20\"/r=\"${radfos}\"/g" ${filename}u_formated.svg
# CHECK NUMBER OF SPACES IN RESPECTIVE FILES
sed -i "/data-safe_div_label=\"${name}\"/,/${name}<\/text>/ s/nwbubblecoloredcircle\(.*\)fill=\".*\"         r=/nwbubblecoloredcircle\1fill=\"${colfos}\"         r=/" ${filename}u_formated.svg
fi
done < fmt_PLL10_col.fmt  # HERE UPDATE

inkscape ${filename}u_formated.svg
```


