#!/usr/bin/env Rscript
### CABANA

options(warn=1)
options(stringsAsFactors=F)

parameters = list(case=list('casefile', 'string', '.bam', T), control=list('controlfile', 'string', '.bam', T), bed=list('bedfile', 'string', '.bed', T))
args = commandArgs(trailingOnly=T)
for (arg in args) {
       argName = sub('=.*', '', sub('^--', '', arg, perl=T), perl=T)
       argValue = sub('^[^=]*=', '', arg, perl=T)
       assign(parameters[[argName]][[1]], argValue)
}


system (paste0('samtools depth -a -b ', bedfile, ' ', casefile, ' > case.depth.txt') )
controls=strsplit(controlfile, ',', fixed=T)[[1]]

n.control=length(controls)
n.w=n.control+1

for (k in 1:n.control){
 system (paste0 ('samtools depth -a -b ', bedfile, ' ', controls[k], ' > control', k, '.depth.txt') )
}

t.bed=read.delim(bedfile, header=F)
colnames (t.bed)[1:4]=c('chrom', 'start', 'end', 'gene')
genes=sort(levels(as.factor(as.character(t.bed$gene))))


cs = paste0('control',1:n.control)
ws = c('case', cs)
depths=paste0 (ws, '.depth.txt')


dp=NULL; 
for(v in 1: n.w ) {
 if (v==1) {dp [[v]]=read.delim (paste0 (depths[v]), header=F)}
 if (v>1) {dp [[v]]=read.delim (paste0 (depths[v]), header=F)}
 colnames (dp [[v]]) [1:3]=c('chrom', 'pos', 'depth')
}

dp.total=NULL
for(v in 1:n.w ) {
 dp.total [[v]]=sum (as.numeric (dp [[v]]$depth) , na.rm=T)
}
dp.mean=mean (dp.total)


dp.s=NULL; p.s=NULL; p.s.0 =NULL; nd=NULL;

for (k in 1:length (genes)) {
tg.s=subset(t.bed, gene==genes[k])
nchrom=tg.s$chrom[1]
tg.s.s=cbind (tg.s, paste0 (nchrom, ':', tg.s$start, '-', tg.s$end) )
colnames (tg.s.s) [ ncol (tg.s.s) ]='chrom.pos'
ir.min=min(tg.s$start,tg.s$end,na.rm=T)
ir.max=max(tg.s$start,tg.s$end,na.rm=T)

plot.r=NULL
for (h in 1 : nrow (tg.s) ) {
e.s=as.numeric (as.character (tg.s$start [h]) ) : as.numeric (as.character (tg.s$end [h]) )
ess=cbind(e.s, h )
if (h==1) esa=ess else { esa=rbind (esa, ess) }
if (h==1) plotrs=cbind(h, nrow (ess) ) else plotrs=cbind(h, plot.r [h-1,2] + nrow (ess) )
plot.r=rbind (plot.r,  plotrs)
}
colnames (plot.r) [2]='plot.pos'
rownames (plot.r)=plot.r [,1]
plot.r=as.data.frame (plot.r)
plot.mid=( plot.r$plot.pos +c(0, plot.r$plot.pos [0: (nrow (plot.r) -1) ] ) ) / 2
esa=as.data.frame (esa)
esa=cbind (esa, 1:nrow (esa) )
colnames (esa)=c('pos', 'region.n', 'plot.pos' )

for(v in 1:n.w ) {
 dp.s [[v]]=subset (dp [[v]], chrom==nchrom & pos > (ir.min -100) & pos < (ir.max +100) ) 
 p.s.0 [[v]] =merge (esa, dp.s [[v]] , by='pos', sort=F, all.x=T, all.y=F) 
 p.s [[v]]=subset (p.s.0 [[v]], ! is.na (chrom) ) 
 p.s [[v]]=p.s [[v]]  [ order (p.s [[v]]$plot.pos), ]
 p.s [[v]]$depth [ which (is.na (p.s [[v]]$depth), arr.ind=T ) ]=0
 p.s [[v]]$n.depth=(p.s [[v]]$depth * dp.mean) / dp.total [[v]]
 if (v==1) {
   nd=p.s [[v]]
   colnames (nd) [ncol(nd)]=ws [v] 
   } else {
   nd.s=p.s [[v]] [ , c('pos', 'n.depth')]
   colnames (nd.s) [2]=ws [v]
   if(identical(nd$pos, nd.s$pos)) {
    nd=cbind (nd, nd.s [,2])
    colnames (nd) [ncol(nd)]=ws [v]
    } else {
    nd=merge (nd, nd.s, by='pos', all=T, sort=F)
   }
 }
}

m.nd=rowMeans (nd [ cs ] , na.rm=T)
nd.2=nd
for(v in 1:n.w ) {
 nd.2 [ ws [v] ]= round ((nd [ ws [v] ]-m.nd) / m.nd, 3)
}

w.col='#1E4B85'
w.bg.col='#1E4B8522'
l.col='#D6D6D6'
line.col='#96760B77'
n.order=c(1, rep(2,n.control))
plot.col=c(w.col, line.col) [n.order]


cairo_pdf( paste0(genes [k], '.CNV.pdf'), width=11.7, height=8.3)
######## layered depth ###########
y=NULL
x.start=min (p.s [[1]]$plot.pos); x.end=max (p.s [[1]]$plot.pos); x.width=x.end-x.start;
x.lim=c (x.start, x.end )
for(v in 1:n.w ) { y [[v]]=p.s [[v]]$depth }
y.max=max(unlist (y), na.rm=T )
y.lim=c (0, y.max )

par(fig=c(0,1,0.74,1), mar=c (0.5,5,3,2) )
plot (x.lim, y.lim, type='n', xlim=x.lim, ylim=y.lim, xlab=NA, ylab=NA,axes=F, xpd=F, xaxs='i', yaxs='i' )
abline(v=plot.r$plot.pos, col=l.col )
for(v in 2:n.w ) {
lines(p.s[[v]]$plot.pos, p.s[[v]]$depth,  xlim=x.lim,  ylim =y.lim,  col=line.col,  lwd=1.2, xlab='',ylab='', xaxs='i', yaxs='i' )
}
polygon(c(p.s[[1]]$plot.pos, rev (p.s[[1]]$plot.pos) ), c(p.s[[1]]$depth, rep(0, nrow (p.s[[1]])) ),  xlim=x.lim,  ylim =y.lim,  col=w.bg.col,  border=NA,  xlab='',ylab='', xaxs='i', yaxs='i' )
axis (2, las=2, mgp=c (0, 0.7, 0), cex.axis=0.6 )
mtext ('Coverage' , side=2, line= 2.5, cex=0.8, font=2)
mtext (genes [k] , side=3, line= 0.4, cex=1.6, font=4)
box()


######## normalized.depth plot ###########
nd.all=NULL

for (v in 1:n.w) {nd.all=c (nd.all, as.numeric (nd.2 [ws[v]] [,1])) }
nd.all.2=nd.all[which( ! is.infinite(nd.all) & ! is.na(nd.all) )]
if(length(which(nd.all.2 > 1)) > 50) {nd.all.3=nd.all.2} else {nd.all.3 = 1 }
y.n.max=max (c(1, nd.all.3), na.rm=T)
y.n.lim=c (-1, y.n.max)
par(fig=c(0,1, 0.32, 0.73 ), mar=c (8,5,0,2), new=T )
plot (x.lim, y.n.lim, type='n', xlim=x.lim, ylim=y.n.lim, xlab=NA, ylab=NA,axes=F, xpd=F, xaxs='i', yaxs='i' )
abline(v=plot.r$plot.pos, col=l.col )

for(v in 1:n.w ) {
lines(nd.2$plot.pos, as.numeric (nd.2 [ws[v]] [,1]),  xlim=x.lim,  ylim =y.n.lim,  col=plot.col [v],  lwd=1.2, xlab='',ylab='', xaxs='i', yaxs='i' )
}
axis (2, las=2, mgp=c (0, 0.7, 0), cex.axis=0.6 )
mtext ('Normalized' , side=2, line= 2.7, cex=0.8, font=2)
mtext ('depth' , side=2, line= 1.9, cex=0.8, font=2)
box()

par(xpd=T)
for (kq in 1:length(plot.mid)) text(plot.mid [kq], -1.1, labels=tg.s.s$chrom.pos [kq], adj=1, srt=90, cex=0.5)
legend('bottomright', c('case', 'control'),  col=plot.col[1:2], lty=1, cex=0.6, bty='n')

dev.off()


}

