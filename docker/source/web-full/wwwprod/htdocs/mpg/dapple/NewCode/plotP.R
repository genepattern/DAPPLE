##Script to plot networks from BuildNetwork_copy.py
source("/home/unix/rossin/DALY/PPI/NewCode/3Dcircle.R")
source("/home/unix/rossin/DALY/PPI/NewCode/collapseCI.R")
collapsePlot = COLLAPSEHERE

manual.coloring=FALSE
if (manual.coloring == TRUE) {
	coloring.scheme=as.matrix(read.table("coloringScheme"))
}
test = scan("./tmpDirect",what="character")
if (length(test) > 0) {in.set.ppi = as.matrix(read.table("./tmpDirect",sep="\t"))}
if (length(test) == 0) {in.set.ppi = matrix(nrow=0,ncol=2)}
SNP.gene = as.matrix(read.table("./tmpSnpGene"))
SNP.gene = cbind(SNP.gene,as.numeric(as.factor(SNP.gene[,1])))
disease.genes = SNP.gene[,2]
library(igraph)
divBy=2
scores = na.omit(as.matrix(read.table("KEYWORDHERE_seedScores",header=T)))

grade = 0.01
#colorMat = matrix(nrow=1/grade,ncol=3)
#colorMat[,1] = c(seq(0,1-grade,grade))
#colorMat[,2] = c(seq(grade,1,grade))
#colorMat[,3] = c(seq(grade/2,1-grade/2,grade))/2

colorMat = matrix(nrow=17,ncol=3)
colorMat[,1] = c(1/10001,0.00025,0.0005,0.00075,0.001,0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,.1,.2,.3,.4,.5)
colorMat[,2] = c(0.00025,0.0005,0.00075,0.001,0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,.1,.2,.3,.4,.5,1)
colorMat[,3] = seq(0,.2,.2/16)

colorByScore = unlist(lapply(scores[,4],function(x){return(colorMat[colorMat[,1]<as.numeric(x)&colorMat[,2]>=as.numeric(x),3])}))
names(colorByScore)=scores[,1]

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
		vertex.col=unlist(lapply(unique.list,function(x){return(SNP.gene[which(SNP.gene[,2]==x),3])}))
		legendcols = as.numeric(as.factor(unique(vertex.col))); legendlabels = unlist(lapply(legendcols,function(x){return(SNP.gene[which(SNP.gene[,3]==x),1][1])}))
		vertex.border<-vertex.col
		num.colors<-max(SNP.gene[,3],na.rm=T)
		palette(rainbow(num.colors,s=.4))
		label.col<-rep("black",length(unique.list))
	}
		
	if(manual.coloring==TRUE) {
		groups<-coloring.scheme
		vertex.col = groups[match(unique.list,groups[,1]),2]
		vertex.col[is.na(vertex.col)]="grey"
		vertex.border=vertex.col
		}

	pdf("KEYWORDHERE_plot.pdf")
	coords=layout.fruchterman.reingold(net,xmin=-1,xmax=3)
	coords2=layout.norm(coords,xmin=-1,xmax=1,ymin=-1,ymax=1)
	coloring = colorByScore[unique.list]

	rownames(coords2) = unique.list
	plot.new()
	plot.window(xlim=c(-1.1,1.5),ylim=c(-1.1,1.1))
	
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
	if (num < 5) { num=0 }
	for (i in 1:length(vertex.col)) {
		draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r*1.2,hue="grey",num=0)
		if (unique.list[i]%in%names(colorByScore)) {
			draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r,hue=colorByScore[unique.list[i]],sat=.55,num=num)
		} else {
			draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r,hue="grey",sat=.55,num=num)			
		}
	}
	text(coords2,label=unique.list,cex=text.cex,font=2,pos=4,offset=0.5)
	
	
#	colorByScoreLegend = hsv(colorMat[,3],.55,1)
	colorByScoreLegend = hsv(colorMat[,3],.55,1)
#	y=(length(colorByScoreLegend)*grade)/2
	y=17*.05/2
	for (i in 1:length(colorByScoreLegend)) {
#		rect(1.25,y-i*(grade),1.35,y-i*grade+grade,col=colorByScoreLegend[i],border=NA)
		rect(1.25,y-i*(0.05),1.35,y-i*0.05+0.05,col=colorByScoreLegend[i],border=NA)
		if (i %in% seq(1,17,2)){text(1.35,y-i*0.05+.025,labels=signif(colorMat[i,2],1),cex=.7,offset=0,pos=4)}
	}
#	text(1.35,y-grade/2,label=paste("p=",signif(1/npermute,1)),cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/20)*grade,label="p=0.05",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/4)*grade,label="p=0.25",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/2)*grade,label="p=0.5",cex=.7,offset=0,pos=4)
#	text(1.35,y-(3*length(colorByScoreLegend)/4)*grade,label="p=0.75",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend))*grade+grade/2,label="p=1",cex=.7,offset=0,pos=4)
	
	if(manual.coloring==TRUE){
		legend(1.2,0,yjust=0.5,fill=unique(vertex.col),legend=rep("",length(unique(vertex.col))),cex=.75,bty='n')
	}	
}
#dev.off();quit();		
if (nrow(in.set.ppi)==0) {pdf("KEYWORDHERE_plot.pdf")}
pairs = as.matrix(read.table("./tmpIndirect",sep="\t"))
pairs = rbind(pairs,in.set.ppi)
pairsOriginal = pairs
net=graph.edgelist(pairs,directed=F)
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
for (i in 1:length(unique.list)) {
	cat(c(unique.list[i],coordsToWrite[i,1],coordsToWrite[i,2]),sep="\t",append=T,file="PlotCoordinates")
	cat("\n",append=T,file="PlotCoordinates")
}	
if (collapsePlot) {
	new=collapseCI(pairs,SNP.gene)
	pairs=new$pairs
	n = new$count
	proteins = unique(c(pairs))
	count = 0
	labelcount = 0 
	plot.new()
	par(mar=c(0,0,0,0)+.1)
	plot.window(xlim=c(0,1),ylim=c(0,1))
	text(0,.99,labels="You chose to simplify the indirect plot my merging common interactors.",pos=4,font=2,cex=.95)
	text(0,.95,labels="This merges any common interactors that bind the same disease proteins.",pos=4,cex=.95)
	text(0,.91,labels=paste(n," substitutions have been made. The details were placed in KEYWORDHERE_CIsimplifyDetails.",sep=""),pos=4,cex=.95)
	cat("You chose to simplify the indirect plot my merging common interactors.\n",file="KEYWORDHERE_CIsimplifyDetails")
	cat("This merges any common interactors that bind the same disease proteins.\n",file="KEYWORDHERE_CIsimplifyDetails",append=T)
	cat("The following substitutions have been made:\n",file="KEYWORDHERE_CIsimplifyDetails",append=T)
	cat("New Label : Proteins that were merged\n",file="KEYWORDHERE_CIsimplifyDetails",append=T)
	#text(0,.85-.85*count/(n-1),labels="New Label : Proteins that were merged",cex=.75,pos=4)
	count=count+1
	for (i in 1:length(proteins)) {
		x=unlist(strsplit(proteins[i],","))
		if (length(x)>1) {
			count = count+1
			labelcount = labelcount+1
#			if(length(x)>15) {
				#text(0,.85-.85*count/(n-1),labels=paste(labelcount,paste(x[1:15],collapse=","),sep=" : "),cex=.7,pos=4)
				cat(paste(labelcount," : ",paste(x,collapse=","),"\n",sep=""),file="KEYWORDHERE_CIsimplifyDetails",append=T)
				count=count+1
#if(length(x)>30) {
#					text(.05,.85-.85*count/(n-1),labels=paste(x[16:30],collapse=","),cex=.7,pos=4)
#					count=count+1
#					text(.05,.85-.85*count/(n-1),labels=paste(x[31:length(x)],collapse=","),cex=.7,pos=4)
#				} else {
#					text(.05,.85-.85*count/(n-1),labels=paste(x[16:length(x)],collapse=","),cex=.7,pos=4)
#				}
#			} else {
#				text(0,.85-.85*count/(n-1),labels=paste(labelcount,proteins[i],sep=" : "),cex=.7,pos=4)
#			}
			pairs[pairs[,1]==proteins[i],1]=as.character(count)
			pairs[pairs[,2]==proteins[i],2]=as.character(count)
		}
	}
	if (count==0){
		text(0,.5,labels="No substitutions could be made.",pos=4)
		cat("No substitutions could be made.\n",file="KEYWORDHERE_CIsimplifyDetails",append=T)
	}
}	
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
params=list()
params$area = (2*vcount(net))^2
#coords=layout.fruchterman.reingold(net,params)
coords = layout.graphopt(net)
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
		draw3Dcircle(coords2[i,1],coords2[i,2],r=rCI,"grey",num=0)
	}
}
for (i in 1:length(vertex.col)) {
	if (vertex.col[i]!=0) {
		draw3Dcircle(x=coords2[i,1],y=coords2[i,2],r=r*1.2,hue="grey",num=0)		
		if (unique.list[i] %in% names(colorByScore)) {
			draw3Dcircle(coords2[i,1],coords2[i,2],r=r,hue=colorByScore[unique.list[i]],num=0)
		} else {
			draw3Dcircle(coords2[i,1],coords2[i,2],r=r,hue="grey",num=0)
		}		
	}
}
textsize=rep(text.cex,length(vertex.col)); textsize[vertex.col==0] = text.cexCI
text(coords2,label=unique.list,cex=textsize,font=1,pos=1,offset=0)
#	colorByScoreLegend = hsv(colorMat[,3],.55,1)
colorByScoreLegend = hsv(colorMat[,3],.55,1)
#	y=(length(colorByScoreLegend)*grade)/2
y=17*.05/2
for (i in 1:length(colorByScoreLegend)) {
#		rect(1.25,y-i*(grade),1.35,y-i*grade+grade,col=colorByScoreLegend[i],border=NA)
	rect(2,y-i*(0.05),2.1,y-i*0.05+0.05,col=colorByScoreLegend[i],border=NA)
	if (i %in% seq(1,17,2)){text(2.1,y-i*0.05+.025,labels=signif(colorMat[i,2],1),cex=.7,offset=0,pos=4)}
}
#	text(1.35,y-grade/2,label=paste("p=",signif(1/npermute,1)),cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/20)*grade,label="p=0.05",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/4)*grade,label="p=0.25",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/2)*grade,label="p=0.5",cex=.7,offset=0,pos=4)
#	text(1.35,y-(3*length(colorByScoreLegend)/4)*grade,label="p=0.75",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend))*grade+grade/2,label="p=1",cex=.7,offset=0,pos=4)

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
		if (unique.list[i] %in% names(colorByScore)) {
			draw3Dcircle(coords2[i,1],coords2[i,2],r=r,hue=colorByScore[unique.list[i]],num=0)
		} else {
			draw3Dcircle(coords2[i,1],coords2[i,2],r=r,hue="grey",num=0)
		}		
	}
}
new.labels=unique.list; new.labels[vertex.col==0]=""
text(coords2,label=new.labels,cex=text.cex,font=2,pos=1,offset=0)
#	colorByScoreLegend = hsv(colorMat[,3],.55,1)
colorByScoreLegend = hsv(colorMat[,3],.55,1)
#	y=(length(colorByScoreLegend)*grade)/2
y=17*.05/2
for (i in 1:length(colorByScoreLegend)) {
#		rect(1.25,y-i*(grade),1.35,y-i*grade+grade,col=colorByScoreLegend[i],border=NA)
	rect(2,y-i*(0.05),2.1,y-i*0.05+0.05,col=colorByScoreLegend[i],border=NA)
	if (i %in% seq(1,17,2)){text(2.1,y-i*0.05+.025,labels=signif(colorMat[i,2],1),cex=.7,offset=0,pos=4)}
}
#	text(1.35,y-grade/2,label=paste("p=",signif(1/npermute,1)),cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/20)*grade,label="p=0.05",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/4)*grade,label="p=0.25",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend)/2)*grade,label="p=0.5",cex=.7,offset=0,pos=4)
#	text(1.35,y-(3*length(colorByScoreLegend)/4)*grade,label="p=0.75",cex=.7,offset=0,pos=4)
#	text(1.35,y-(length(colorByScoreLegend))*grade+grade/2,label="p=1",cex=.7,offset=0,pos=4)


###Plot zoomed plot
zoomedGenes = scan("zoomedGeneFile",what="character")
if (zoomedGenes[1]=="none" | length(zoomedGenes) == 0) {
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
rownames(coords2)=unique.list
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


