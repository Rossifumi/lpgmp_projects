normalize.rows.f=function (Mat) t(apply(Mat, 1, function(y) (y-mean(y))/sd(y)))

sum_series<-function(x)			#//here x is 'size'
{
  sum_se<-c(x[1])
  sx<-x[1]
  for(i in 2:length(x))
  {
    sx<-x[i]+sx
    sum_se<-c(sum_se,sx)
  }
  return(sum_se)					#//here sum_se is 'sum of size(length(sum_se)==length(size))'
}

k_means<-function(mx,cln,itermax,filename,age)	#//(a,?,10,Name,c(1:12))
{
  library(cluster)
  km<-kmeans(mx,centers=cln,iter.max=itermax)
  index<-km$cluster
  means<-km$centers
  size<-km$size
  sort.index<-sort(index, index.return = T)$ix	#//index的序号(from 1:length)
  sizes<-sum_series(size)
  
  
  #//write cluster and means
  write.table(index,file=paste(filename,"_cluster.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\n")
  write.table(means,file=paste(filename,"_means.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
  
  #//get ids for each cluster
  clusters_id<-c(paste("clusters_id_",c(1:cln),sep=""))		#//filenames of ens of each cluster
  assign(clusters_id[1],Mat_Ori[,1][sort.index[1:sizes[1]]])	
  write.table(get(clusters_id[1]),"clusters_id_1.txt",quote=F,col.names=F,row.names=F,sep="\n")
  for (i in 2:cln)
  {
    assign(clusters_id[i],Mat_Ori[,1][sort.index[(sizes[i-1]+1):sizes[i]]])
    write.table(get(clusters_id[i]),paste("clusters_id_",i,".txt",sep=""),quote=F,col.names=F,row.names=F,sep="\n")
  }
  
  clusters<-c(paste("cluster",c(1:cln),sep=""))		#//filenames
  assign(clusters[1],mx[sort.index[1:sizes[1]],])		#//get expression data of each cluster(names in 'clusters')exp,dim 17 12
  
  #//write "cluster1" in the final txt file
  
  write.table(get(clusters[1]),file=paste(filename,".txt",sep=""),quote=F,col.names=F,row.names=T,append=T,sep="\t")
  
  for(i in 2:cln)		#//appending results 
  {
    assign(clusters[i],mx[sort.index[(sizes[i-1]+1):sizes[i]],])
    #write.table(paste("cluster",i,sep=""),file=paste(filename,".txt",sep=""),quote=F,col.names=F,row.names=F,append=T,sep="\t")
    write.table(get(clusters[i]),file=paste(filename,".txt",sep=""),quote=F,col.names=F,row.names=T,append=T,sep="\t")
  }
  
  #plot 
  pdf(paste(filename,".pdf",sep=""))
  par(mfrow=c(2,2),las=1)		#//one page four subplot by row
  yrange<-c((round(min(mx))-0.1),(round(max(mx))+0.1))
  
  for(i in  1: cln)
  {
    plot(age,get(clusters[i])[1,],type="o",cex=.4,col="yellow",ylim=yrange,xlab="condition")	#//draw first row of expression data of each cluster
    
    for(j in i:dim(get(clusters[i]))[1])		#//dim[1]==rows number,for those rows,plot
    {
      lines(age,get(clusters[i])[j,],type="o",cex=.4,col="grey",ylim=yrange)			#//target补上去的线
    }
    
    lines(age,means[i,],type="l",cex=.9,ylim=yrange)	#//draw the mean plot
    mtext(paste("Size",size[i]," No",i,sep=""))		#//text on each subplot
  }
  
  dev.off()
  Ages<-age
  
  for(i in 1:cln)
  {
    pdf(paste("cluster",i,".pdf",sep=""))	#//for each cluster,plot each gene's expression data
    for(j in 1:dim(get(clusters[i]))[1])
    {
      plot(Ages,get(clusters[i])[j,],xlab="condition",ylab="")
      lines(smooth.spline(Ages,get(clusters[i])[j,]),col="blue")
      lines(smooth.spline(Ages,get(clusters[i])[j,],df= 10),lty=2,col="red")
    }
    dev.off()
  }
  
  
  #//mainn show plot
  pdf(paste(filename,".all.pdf",sep=""),width=12,height=10)
  par(mfrow=c(5,6))
  par(mai=c(0,0,0,0))
  yrange<-c((round(min(mx))-0.1),(round(max(mx))+0.1))
  
  
  for(i in  1: cln)
  {
    plot(age,get(clusters[i])[1,],type="o",cex=.4,col="grey",ylim=yrange,xaxt="n",yaxt="n")	#//draw first row of expression data of each cluster
    
    for(j in i:dim(get(clusters[i]))[1])		#//dim[1]==rows number,for those rows,plot
    {
      lines(age,get(clusters[i])[j,],type="o",cex=.4,col="grey",ylim=yrange,xaxt="n",yaxt="n")			#//target补上去的线
    }
    
    lines(age,means[i,],type="l",cex=.9,col="red",ylim=yrange,xaxt="n",yaxt="n")	#//draw the mean plot
    mtext(side=3,line=-2,paste("Size:",size[i]," No",i,sep=""),cex=.9)		#//text on each subplot
  }
  
  dev.off()
  
  #//same plot in png format
  png(paste(filename,".all.png",sep=""),width = 1200, height = 1000, units = "px", pointsize = 12,bg = "white")
  par(mfrow=c(5,6))
  par(mai=c(0,0,0,0))
  yrange<-c((round(min(mx))-0.1),(round(max(mx))+0.1))
  
  
  for(i in  1: cln)
  {
    plot(age,get(clusters[i])[1,],type="o",cex=.4,col="grey",ylim=yrange,xaxt="n",yaxt="n")	#//draw first row of expression data of each cluster
    
    for(j in i:dim(get(clusters[i]))[1])		#//dim[1]==rows number,for those rows,plot
    {
      lines(age,get(clusters[i])[j,],type="o",cex=.4,col="grey",ylim=yrange,xaxt="n",yaxt="n")			#//target补上去的线
    }
    
    lines(age,means[i,],type="l",cex=.9,col="red",ylim=yrange,xaxt="n",yaxt="n")	#//draw the mean plot
    mtext(side=3,line=-2,paste("Size:",size[i]," No",i,sep=""),cex=.9)		#//text on each subplot
  }
  
  dev.off()
  
  

  #//plot means for each cluster
  pdf(paste(filename,"_means.pdf",sep=""))
  for (i in 1:cln)
  {
    plot(means[i,],type="p",xlab="condition",ylab=paste(filename,"_means_",i,sep=""))
    lines(smooth.spline(means[i,]),col="red")
  }
  dev.off()
}

