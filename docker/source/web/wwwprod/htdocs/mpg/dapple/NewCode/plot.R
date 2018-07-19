##Script to plot networks from BuildNetwork_copy.py

args <- commandArgs(trailingOnly = TRUE)
#1 = COLLAPSEHERE
#2 = KEYWORD
#3 = manualColoring

source("/home/unix/rossin/DALY/PPI/NewCode/3Dcircle.R")
source("/home/unix/rossin/DALY/PPI/NewCode/collapseCI.R")
collapsePlot = args[1]
manual.coloring=args[3]
if (manual.coloring == TRUE) {
	coloring.scheme=as.matrix(read.table(args[4]))
}
test = scan("./tmpDirect",what="character")
if (length(test) > 0) {in.set.ppi = as.matrix(read.table("./tmpDirect",sep="\t"))}
if (length(test) == 0) {in.set.ppi = matrix(nrow=0,ncol=2)}
SNP.gene = as.matrix(read.table("./tmpSnpGene"))
SNP.gene[,1] = unlist(lapply(SNP.gene[,1],function(x){if(unlist(strsplit(x,"_"))[1]==""){return(unlist(strsplit(x,"_"))[2])}else{return(unlist(strsplit(x,"_"))[1])}}))
SNP.gene = cbind(SNP.gene,as.numeric(as.factor(SNP.gene[,1])))
disease.genes = SNP.gene[,2]
library(igraph)
if (file.exists("zoomedGeneFile")) {
	zoomedGenes = scan("zoomedGeneFile",what="character")
} else {
	zoomedGenes = "none"
}

if (nrow(in.set.ppi)>0) {
	unique.list<-unique(c(t(in.set.ppi)))
	net<-graph.edgelist(in.set.ppi,directed=F)

	if (nrow(in.set.ppi)==1) {
		in.set.ppi = rbind(in.set.ppi,c(in.set.ppi[,2],in.set.ppi[,1]))
	}

	if (nrow(in.set.ppi)==2) {
#		tmp<-matrix(0,nrow=length(unique(c(in.set.ppi))),ncol=length(unique(c(in.set.ppi))))
#		rownames(tmp)=sort(unique(c(in.set.ppi)))
#		colnames(tmp)=sort(unique(c(in.set.ppi)))
#		tmp[in.set.ppi[1,1],in.set.ppi[1,2]]=tmp[in.set.ppi[1,2],in.set.ppi[1,1]]=1
#		tmp[in.set.ppi[2,1],in.set.ppi[2,2]]=tmp[in.set.ppi[2,2],in.set.ppi[2,1]]=1
		tmp = rbind(in.set.ppi,cbind(in.set.ppi[,2],in.set.ppi[,1]))
		net=graph.edgelist(tmp)
	}	
	
	#Figure out coloring of vertices
	if(manual.coloring==FALSE) {
		vertex.col=as.numeric(unlist(lapply(unique.list,function(x){return(SNP.gene[which(SNP.gene[,2]==x),3])})))
		num.colors<-max(SNP.gene[,3],na.rm=T)
		colorpalette = seq(0,1,1/max(vertex.col)); colorpalette = sample(colorpalette,length(colorpalette))
		label.col<-rep("black",length(unique.list))
		legendcols = hsv(colorpalette[as.numeric(unique(vertex.col))],.4,1); legendlabels = unlist(lapply(as.numeric(unique(vertex.col)),function(x){return(SNP.gene[which(SNP.gene[,3]==x),1][1])}))
	}
		
	if(manual.coloring==TRUE) {
		groups<-coloring.scheme
		vertex.col = groups[match(unique.list,groups[,1]),2]
		vertex.col[is.na(vertex.col)]="grey"
		vertex.border=vertex.col
		}
	#Figure out plot dimensions
	expand = 1
	for (label in legendlabels) {
		n = length(unlist(strsplit(label,"_")))
		if (n>expand) { expand = n } }

	pdf(paste(args[2],"_plot.pdf",sep=""))
	coords=layout.fruchterman.reingold(net,xmin=-1,xmax=3)
	coords2=layout.norm(coords,xmin=-1,xmax=1,ymin=-1,ymax=1)
	rownames(coords2) = unique.list
	plot.new()
	plot.window(xlim=c(-1.1,1.3+expand/3),ylim=c(-1.1,1.1))

	#draw edges
	for (i in 1:nrow(in.set.ppi)) {
		p1 = in.set.ppi[i,1]
		p2 = in.set.ppi[i,2]
		segments(coords2[p1,1],coords2[p1,2],coords2[p2,1],coords2[p2,2],col=hsv(.6,.05,.7))
	}
	
	#draw circles
	r = 0.04*(2^log10(30/(length(vertex.col))))
	text.cex=.55*(2^log10(30/(length(vertex.col))))
	num=floor(20*(2^log10(30/(length(vertex.col)))))
	for (i in 1:length(vertex.col)) {
		draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r*1.2,hue="grey",num=0)		
		draw3Dcircle(coords2[i,1],coords2[i,2],r=r,hue=colorpalette[vertex.col[i]],num=num)
	}
	text(coords2,label=unique.list,cex=text.cex,font=2,pos=4,offset=0.4)
	legend(1.2,0,yjust=0.5,fill=legendcols,legend=legendlabels,cex=.5,bty='n')
	if(manual.coloring==TRUE){
		legend(1.2,0,yjust=0.5,fill=unique(vertex.col),legend=rep("",length(unique(vertex.col))),cex=.75,bty='n')
	}	
}

if (nrow(in.set.ppi)==0) {pdf(paste(args[2],"_plot.pdf",sep=""))}

pairs = scan("./tmpIndirect",what="character")
if (length(pairs)>0){
	pairs = as.matrix(read.table("./tmpIndirect",sep="\t"))
} else {
	pairs = matrix(0,nrow=0,ncol=2)
}
pairs = rbind(pairs,in.set.ppi)
net=graph.edgelist(pairs,directed=F)
pairsOriginal = pairs
coords=layout.fruchterman.reingold(net,area = (vcount(net))^2.3)
plotSize=1.8
coordsOriginal=layout.norm(coords,xmin=-plotSize,xmax=plotSize,ymin=-plotSize,ymax=plotSize)
coords = layout.graphopt(net)
coordsToWrite = coords
coordsToWrite[,1]=coordsToWrite[,1]+abs(min(coordsToWrite[,1]));coordsToWrite[,1]=coordsToWrite[,1]*600/max(coordsToWrite)
coordsToWrite[,2]=coordsToWrite[,2]+abs(min(coordsToWrite[,2]));coordsToWrite[,2]=coordsToWrite[,2]*400/max(coordsToWrite)
cat("",file="PlotCoordinates")
unique.list<-unique(c(t(pairs)))
unique.listOriginal = unique.list
rownames(coordsOriginal)=unique.list
for (i in 1:length(unique.list)) {
	cat(c(unique.list[i],coordsToWrite[i,1],coordsToWrite[i,2]),sep="\t",append=T,file="PlotCoordinates")
	cat("\n",append=T,file="PlotCoordinates")
}
n=1
if (collapsePlot) {
	new=collapseCI(pairs,SNP.gene)
	pairs=new$pairs
	n = new$count
	proteins = unique(c(pairs))
	labelcount = 0 
	plot.new()
	par(mar=c(0,0,0,0)+.1)
	plot.window(xlim=c(0,1),ylim=c(0,1))
	text(0,.99,labels="You chose to simplify the indirect plot my merging common interactors.",pos=4,font=2,cex=.95)
	text(0,.95,labels="This merges any common interactors that bind the same disease proteins.",pos=4,cex=.95)
	text(0,.91,labels=paste(n," substitutions have been made. The details were placed in ",args[2],"_CIsimplifyDetails.",sep=""),pos=4,cex=.95)
	cat("You chose to simplify the indirect plot my merging common interactors.\n",file=paste(args[2],"_CIsimplifyDetails",sep=""))
	cat("This merges any common interactors that bind the same disease proteins.\n",file=paste(args[2],"_CIsimplifyDetails",sep=""),append=T)
	cat("The following substitutions have been made:\n",file=paste(args[2],"_CIsimplifyDetails",sep=""),append=T)
	cat("New Label : Proteins that were merged\n",file=paste(args[2],"_CIsimplifyDetails",sep=""),append=T)
#text(0,.85-.85*count/(n-1),labels="New Label : Proteins that were merged",cex=.75,pos=4)
	for (i in 1:length(proteins)) {
		x=unlist(strsplit(proteins[i],","))
		if (length(x)>1) {
			labelcount = labelcount+1
			cat(paste(labelcount," : ",paste(x,collapse=","),"\n",sep=""),file=paste(args[2],"_CIsimplifyDetails",sep=""),append=T)
			pairs[pairs[,1]==proteins[i],1]=as.character(labelcount)
			pairs[pairs[,2]==proteins[i],2]=as.character(labelcount)
		}
	}
	if (labelcount==0){
		text(0,.5,labels="No substitutions could be made.",pos=4)
		cat("No substitutions could be made.\n",file=paste(args[2],"_CIsimplifyDetails",sep=""),append=T)
	}
}	

net=graph.edgelist(pairs,directed=F)

if (nrow(pairs)>10000) {
	plot.new()
	plot.window(xlim=c(-1.1,1.1),ylim=c(-1.1,1.1))
	text(0,0,labels="The indirect network contains more than 10,000 edges\nand was not plotted.\nContact dapple@broadinstitute.org for assistance.")	
	dev.off();quit()}

edge.col=c(rep("grey",nrow(pairs)),rep("black",nrow(in.set.ppi)))

unique.list<-unique(c(t(pairs)))

if(nrow(pairs)==2) {pairs = rbind(pairs,cbind(pairs[,2],pairs[,1]))}
net=graph.edgelist(pairs,directed=F)
vertex.radius=unlist(lapply(unique.list,function(x){if(x%in%disease.genes){return(.04)}else{return(.02)}}))
label.col="black"

	
if(manual.coloring==TRUE) {
	groups = coloring.scheme
	vertex.col = groups[match(unique.list,groups[,1]),2]
	vertex.col[is.na(vertex.col)]="grey"
	vertex.border=vertex.col
	}

if(manual.coloring==FALSE) {
	vertex.col=as.numeric(unlist(lapply(unique.list,function(x){if(x%in%SNP.gene[,2]){return(SNP.gene[which(SNP.gene[,2]==x),3])} else {return(0)}})))
	num.colors<-max(SNP.gene[,3],na.rm=T)
	colorpalette = seq(0,1,1/max(vertex.col)); colorpalette = sample(colorpalette,length(colorpalette))
	legendcols = hsv(colorpalette[as.numeric(unique(vertex.col[vertex.col!=0]))],.4,1); legendlabels = unlist(lapply(as.numeric(unique(vertex.col[vertex.col!=0])),function(x){return(SNP.gene[which(SNP.gene[,3]==x),1][1])}))
	legendcols = c(legendcols,"grey50"); legendlabels=c(legendlabels,"Common Interactor")
}


# To plot both labeled and unlabeled

plotSize = 1.8
coords=layout.fruchterman.reingold(net,area = (vcount(net))^2.3)
coords2=layout.norm(coords,xmin=-plotSize,xmax=plotSize,ymin=-plotSize,ymax=plotSize)
rownames(coords2) = unique.list

plot.new()
plot.window(xlim=c(-plotSize,plotSize+.5),ylim=c(-plotSize,plotSize))

#draw edges
for (i in 1:nrow(pairs)) {
	p1 = pairs[i,1]
	p2 = pairs[i,2]
	segments(coords2[p1,1],coords2[p1,2],coords2[p2,1],coords2[p2,2],col=hsv(.6,.05,.7))
}

r = 0.08*(2^log10(30/(length(vertex.col))))
rCI = 0.04*(2^log10(30/(length(vertex.col))))
text.cex=.75*(2^log10(30/(length(vertex.col))))
text.cexCI=.55*(2^log10(30/(length(vertex.col))))
num=floor(20*(2^log10(30/(length(vertex.col)))))

#draw circles
for (i in 1:length(vertex.col)) {
	if (vertex.col[i]==0) {
		if (unique.list[i] %in% as.character(1:(n+1))) {
			draw3Dcircle(coords2[i,1],coords2[i,2],r=0.06*(2^log10(30/(length(vertex.col)))),hue="grey",num=0)
		} else {
			draw3Dcircle(coords2[i,1],coords2[i,2],r=rCI,"grey",num=0)
		}
	}
}
for (i in 1:length(vertex.col)) {
	if (vertex.col[i]!=0) {
		draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r*1.2,hue="grey",num=0)		
		draw3Dcircle(coords2[i,1],coords2[i,2],r=r,colorpalette[vertex.col[i]],num=0)
	}
}
textsize=rep(text.cex,length(vertex.col)); textsize[vertex.col==0] = text.cexCI
text(coords2,label=unique.list,cex=textsize,font=1)
legend(2,0,yjust=0.5,fill=legendcols,legend=legendlabels,cex=.5,bty='n')

plot.new()
plot.window(xlim=c(-plotSize,plotSize+.5),ylim=c(-plotSize,plotSize))

#draw edges
for (i in 1:nrow(pairs)) {
	p1 = pairs[i,1]
	p2 = pairs[i,2]
	segments(coords2[p1,1],coords2[p1,2],coords2[p2,1],coords2[p2,2],col=hsv(.6,.05,.7),lwd=.5)
}
	
#draw circles
for (i in 1:length(vertex.col)) {
	if (vertex.col[i]==0) {
		draw3Dcircle(coords2[i,1],coords2[i,2],r=rCI,"grey",num=0)
	}
}
for (i in 1:length(vertex.col)) {
	if (vertex.col[i]!=0) {
		draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r*1.2,hue="grey",num=0)
		draw3Dcircle(coords2[i,1],coords2[i,2],r=r,colorpalette[vertex.col[i]],num=0)
	}
}
new.labels=unique.list; new.labels[vertex.col==0]=""
text(coords2,label=new.labels,cex=text.cex,font=2)
legend(2,0,yjust=0.5,fill=legendcols,legend=legendlabels,cex=.5,bty='n')


if(manual.coloring==TRUE){
	legend(1.2,0,yjust=0.5,fill=unique(vertex.col),legend=rep("",length(unique(vertex.col))),cex=.75,bty='n')
}	

colors = tapply(vertex.col,as.factor(vertex.col),length)
if (manual.coloring==TRUE) {for (i in 1:length(colors)) {text(-1,-1.5+.1*i,labels=paste(names(colors)[i]," = ",colors[i]))}}




###Plot zoomed plot
if (zoomedGenes[1]=="none" | length(zoomedGenes)==0) {
	dev.off()
	quit()
}

plot.new()
plot.window(xlim=c(-plotSize,plotSize+.5),ylim=c(-plotSize,plotSize))

#draw edges
unique.list.reduced = c()
pairs = pairsOriginal
coords2 = coordsOriginal
unique.list=unique.listOriginal
for (i in 1:nrow(pairs)) {
	p1 = pairs[i,1]
	p2 = pairs[i,2]
	if (p1 %in% zoomedGenes | p2 %in% zoomedGenes) {
		segments(coords2[p1,1],coords2[p1,2],coords2[p2,1],coords2[p2,2],col=hsv(.6,.05,.7),lwd=.5)
		unique.list.reduced = c(unique.list.reduced,p1,p2)
	}
}
	
#draw circles

for (i in 1:length(vertex.col)) {
	if (unique.list[i]%in%unique.list.reduced){
		if (vertex.col[i]==0) {
			draw3Dcircle(coords2[i,1],coords2[i,2],r=rCI,"grey",num=0)
		}	
	}	
}	
for (i in 1:length(vertex.col)) {
	if (unique.list[i]%in%unique.list.reduced) {
		if (vertex.col[i]!=0) {
			draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r*1.2,hue="grey",num=0)
			draw3Dcircle(coords2[i,1],coords2[i,2],r=r,colorpalette[vertex.col[i]],num=0)
		}
	}
}
new.labels=c()
legendcols2=c()
coords3 = matrix(0,nrow=0,ncol=2)
print(unique.list)
print(unique.list.reduced)
for (i in 1:nrow(coords2)) {
	if (unique.list[i]%in%unique.list.reduced) {
		coords3 = rbind(coords3,coords2[i,])
		new.labels = c(new.labels, unique.list[i])
		legendcols2 = c(legendcols2,legendcols[i])
	}
}
text(coords3,label=new.labels,cex=text.cex,font=2)
legend(2,0,yjust=0.5,fill=legendcols2,legend=legendlabels,cex=.5,bty='n')


if(manual.coloring==TRUE){
	legend(1.2,0,yjust=0.5,fill=unique(vertex.col),legend=rep("",length(unique(vertex.col))),cex=.75,bty='n')
}	

colors = tapply(vertex.col,as.factor(vertex.col),length)
if (manual.coloring==TRUE) {for (i in 1:length(colors)) {text(-1,-1.5+.1*i,labels=paste(names(colors)[i]," = ",colors[i]))}}

dev.off();quit()

