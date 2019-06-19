# title: Functions
# Author: Yiwen Wang
# Date: Jun. 16th


# ---------------------------------------
#  TSS (Total Sum Scaling) Normalisation
# ---------------------------------------
TSS.divide = function(x){
  (x)/sum(x)
}


#---------------------------------------------------------------------
# Principal component analysis (PCA) with density plots per component
#---------------------------------------------------------------------
Scatter_Density <- function(data = data,batch = batch, trt = trt,expl.var = expl.var,
                            xlim=xlim,ylim=ylim, batch.legend.title = 'Batch', 
                            trt.legend.title = 'Treatment', density.lwd = 0.2,
                            title = NULL, title.cex = 1.5, legend.cex = 0.7, legend.title.cex =0.75){
  data = as.data.frame(data)
  pMain <- ggplot(data = data, aes(x=data[,1], y=data[,2], colour = batch,shape = trt)) + 
    geom_point() + xlab(paste0('PC1: ',round(as.numeric(expl.var[1])*100),'% expl.var')) + 
    ylab(paste0('PC2: ',round(as.numeric(expl.var[2])*100),'% expl.var')) + 
    scale_color_manual(values=color.mixo(1:10)) + theme_bw() + xlim(xlim[1],xlim[2]) + 
    ylim(ylim[1],ylim[2]) + labs(colour=batch.legend.title,shape = trt.legend.title) 
  
  pTop <- ggplot(data,aes(x=data[,1], fill=batch,linetype = trt)) + 
    geom_density(size = density.lwd,alpha=0.5) + ylab('Density') + 
    theme(axis.title.x  = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
          plot.title = element_text(hjust = 0.5,size = rel(title.cex)), 
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.background = element_blank()) + scale_fill_manual(values=color.mixo(1:10)) +
    xlim(xlim[1],xlim[2]) + labs(title = title)
  
  pRight <- ggplot(data,aes(x=data[,2], fill=batch,linetype = trt)) + 
    geom_density(size=density.lwd,alpha=0.5) +  coord_flip() + ylab('Density') +
    theme(axis.title.x = element_text(size = rel(0.8)), 
          axis.title.y  = element_blank(), axis.line = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.background = element_blank()) + scale_fill_manual(values=color.mixo(1:10)) +
    xlim(ylim[1],ylim[2])
  
  g <- ggplotGrob(pMain + theme(legend.position="right",legend.box='horizontal',
                                legend.direction = 'vertical', 
                                legend.key.height = unit(0.2, 'cm'),
                                legend.key.width = unit(0.1, 'cm'),
                                legend.title = element_text(size = rel(legend.title.cex)),
                                legend.spacing.x = unit(0.1, 'cm'),
                                legend.spacing.y = unit(0.1, 'cm'),
                                legend.text = element_text(size = rel(legend.cex))))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  grid.arrange(pTop+ theme(legend.position="none"), legend, pMain + 
                 theme(legend.position="none"), pRight+ theme(legend.position="none"),
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  
}

#----------
# Box plot
#----------

box_plot_fun = function(data = data,x = x,y=y,title=NULL, batch.legend.title = 'Batch',x.angle = 0,
                        x.hjust = 0.5){
  ggplot(data = data, aes(x=x, y=y, fill = x)) + stat_boxplot(geom = "errorbar", width = 0.4) + 
    geom_boxplot() + scale_fill_manual(values = color.mixo(1:10)) + theme_bw() + 
    theme(axis.text.x = element_text(angle = x.angle, hjust = x.hjust), panel.grid = element_blank(),
          axis.title.x = element_blank(), axis.text = element_text(size=10),
          axis.title = element_text(size=12),
          plot.title = element_text(hjust = 0.5,size = rel(1))) + 
    labs(fill= batch.legend.title,y='value',title = title) 
} 

#-------------
# RLE plot
#------------

RleMicroRna2 <- function (object, maintitle = NULL, batch = batch, xlab = NA,
                          legend = TRUE,cex.lab = 1.2, cex.xaxis = 1, 
                          cex.yaxis = 1, abline.lwd=0.5,legend.cex = 0.8,
                          xaxis.dist.ratio = 0.1, outcex = 1, title.cex = 1.3) 
{
  colorfill = color.mixo(batch)
  nARR = dim(object)[2]
  nGEN = dim(object)[1]
  y = apply(object, 1, median)
  mva = matrix(nrow = nGEN, ncol = nARR)
  for (i in 1:nARR) {
    x = object[, i]
    mva[, i] = (x - y)
  }
  med = apply(mva, 2, median)
  MIN = min(mva, na.rm = TRUE)
  MAX = max(mva, na.rm = TRUE)
  par(las = 3)
  plot(med, xlim = c(0, nARR + 1), ylim = c(MIN, MAX), axes = FALSE, 
       xlab = xlab, ylab = "Deviations",cex.lab = cex.lab)
  colnames(mva) = colnames(object)
  res = boxplot(data.frame(mva), outline = TRUE, add = TRUE, col = colorfill,
                xaxt = 'n', outcex = outcex, cex.axis = cex.yaxis) #outcex for outlier
  axis(1,cex.axis=cex.xaxis,at = 1:ncol(object), labels = NA)
  points(med, type = "p", col = "blue",cex = outcex)
  lines(med, type = "l", col = "blue", lty = "dotted")
  title(main = maintitle, cex.main = title.cex)
  abline(0, 0, col = "red",lwd = abline.lwd)
  par(las = 0)
  end_point = 0.5 + ncol(object)  # add degrees to the x axis
  box.max = max(max(res$stats),max(res$out))
  box.min = min(min(res$stats),min(res$out))
  box.range = box.max - box.min
  text(seq(1.2,end_point,by=1), par("usr")[3]-xaxis.dist.ratio*box.range, 
       srt = 60, adj= 1, xpd = TRUE,
       labels = paste(colnames(object)), cex=cex.xaxis)
  if(legend == TRUE){
    legend('topright',legend = unique(batch), pch=15,col = unique(colorfill),cex=legend.cex)
  }
}


#--------------------------
# Percentile normalisation
#--------------------------
percentileofscore = function(df,control.index){
  df.percentile = df
  df.percentile[1:nrow(df),1:ncol(df)] = NA
  for(i in 1:ncol(df)){
    control = sort(df[control.index,i])
    for(j in 1:nrow(df)){
      percentile.strick = sum(control < df[j,i])/length(control)
      percentile.weak = (length(control) - sum(control > df[j,i]))/length(control)
      percentile = (percentile.strick + percentile.weak)/2
      df.percentile[j,i] = percentile
      
    }
  }
  return(df.percentile)
}

################

percentile_norm = function(data = data, batch = batch, trt = trt){
  batch = as.factor(batch)
  trt = as.factor(trt)
  trt.list = list()
  data.pn.df = data.frame()
  for(i in 1:nlevels(batch)){
    trt.each.b = trt[batch == levels(batch)[i]]
    trt.list[[i]] = trt.each.b
    data.each.b.pn = percentileofscore(data[batch == levels(batch)[i],], which(trt.each.b == levels(trt.each.b)[1]))
    data.pn.df = rbind(data.pn.df,data.each.b.pn)
  }
  names(trt.list) = levels(batch)
  data.pn.df.reorder = data.pn.df[rownames(data),]
  return(data.pn.df.reorder)
}


#--------------------------
# Silhouette coefficient
#--------------------------
# function that calculates the silhouette coefficient based on a known cluster (i.e. batch or treatment)
# calculates silhouette width average
calc.sil = function(
  x, # the PC variates
  y1, y2 = NULL, # factor of interest, e.g. known batch info or known treatment info
  name.y1, name.y2 = NULL # character of the factor of interest
){
  library(cluster)
  # calculate the distance, here euclidean is appropriate for PCA, NOT for t-SNE
  dist.res = daisy(x, metric = 'euclidean')
  # for factor 1
  sil.batch.res1 = silhouette(x = as.numeric(y1), dist = dist.res)
  # if factor 2 is provided
  if(!is.null(y2))  sil.batch.res2 = silhouette(x = as.numeric(y2), dist = dist.res)
  
  # extract average width silhouette per level
  res1 = c(summary(sil.batch.res1)["clus.avg.widths"]$clus.avg.widths)
  names(res1) = levels(y1)
  
  
  if(!is.null(y2)){
    res2 = c(summary(sil.batch.res2)["clus.avg.widths"]$clus.avg.widths)
    names(res2) = levels(y2)
  }
  
  # output data for plotting
  if(!is.null(y2)){
    silh.coeff = c(res1, res2)
    Cluster = c(levels(y1), levels (y2))
    Type = c(rep(name.y1, nlevels(y1)), rep(name.y2, nlevels(y2)))
    data.plot = data.frame(silh.coeff, Cluster, Type)
    
  }else{
    silh.coeff = c(res1)
    Cluster = c(levels(y1))
    Type = rep(name.y1, nlevels(y1))
    data.plot = data.frame(silh.coeff, Cluster, Type)
  }
  
  return(invisible(data.plot))
}
