all <- read.tree("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Trees/Amphibolurinae_MCMCTree.tre")
#all$tip.label <- sapply(all$tip.label, function(x) strex::str_after_first(x,"_"))
all <- drop.tip(all, c("Moloch_sp.","Diporiphora_linga2","Pogona_minor_minor2",
                       "Sphenodon_punctatus","Phrynocephalus_putjatai"))  
all$edge.length <- all$edge.length*100
  
agam.tree2 <- ladderize(all)

arc_height<-0.4
h<-max(nodeHeights(agam.tree2))
plotTree2(agam.tree2,type="arc",lwd=1,fsize=0.45,ftype="i",
         arc_height=arc_height,ylim=c(-0.1*h,1.1*(1+arc_height)*h))
labs<-seq(0,h,by=5)
a1<-axis(1,pos=-0.02*h,at=h-labs+arc_height*h,
         labels=labs,cex.axis=0.8,lwd=2,lend=2)
text(mean(a1),-0.23*h,"million years bp",font=3)
a2<-axis(1,pos=-0.02*h,at=-h+labs-arc_height*h,
         labels=labs,cex.axis=0.8,lwd=2,lend=2)
text(mean(a2),-0.23*h,"million years bp",font=3)
draw.arc(0,0,radius=h-labs[2:length(labs)]+arc_height*h,
         angle1=,angle2=pi,col=make.transparent("blue",0.4),
         lty="dotted")



# UPDATED FUNCTIONS FOR PLOTTING TREE


# function to plot simmap tree in type "fan"
# written by Liam J. Revell 2013-2017
plotFan2<-function(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,ylim,tips,maxY,lend,plot,offset){
  if(!plot) cat("plot=FALSE option is not permitted for type=\"fan\". Tree will be plotted.\n")
  if(is.null(offset)) offset<-1
  # reorder
  cw<-reorder(tree)
  pw<-reorder(tree,"pruningwise")
  # count nodes and tips
  n<-Ntip(cw)
  m<-cw$Nnode 
  # get Y coordinates on uncurved space
  Y<-vector(length=m+n)
  if(is.null(tips)) tips<-1:n
  if(part<1.0) Y[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)
  else Y[cw$edge[cw$edge[,2]<=n,2]]<-tips
  nodes<-unique(pw$edge[,1])
  for(i in 1:m){
    desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
    Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
  }
  if(is.null(maxY)) maxY<-max(Y)
  Y<-setNames(Y/maxY*2*pi,1:(n+m))
  Y<-part*cbind(Y[as.character(cw$edge[,2])],Y[as.character(cw$edge[,2])])
  R<-nodeHeights(cw)
  # now put into a circular coordinate system
  x<-R*cos(Y)
  y<-R*sin(Y)
  # optimize x & y limits
  par(mar=mar)
  offsetFudge<-1.37 # empirically determined
  OFFSET<-0
  pp<-par("pin")[1]
  sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
    offsetFudge*OFFSET*fsize*strwidth("W",units="inches") 
  alp<-optimize(function(a,H,sw,pp) (2*a*1.04*max(H)+2*sw-pp)^2,H=R,sw=sw,pp=pp,
                interval=c(0,1e6))$minimum
  if(part<=0.25) x.lim<-y.lim<-c(0,max(R)+sw/alp)
  else if(part>0.25&&part<=0.5){ 
    x.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
    y.lim<-c(0,max(R)+sw/alp)
  } else x.lim<-y.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
  if(is.null(xlim)) xlim<-x.lim
  if(is.null(ylim)) ylim<-y.lim
  # plot tree
  if(!add) plot.new()
  plot.window(xlim=xlim,ylim=ylim,asp=1)
  # plot radial lines (edges)
  ## first, the lines emerging from the root (if there are only two):
  jj<-which(cw$edge[,1]==(Ntip(cw)+1))
  if(length(jj)==2){
    m.left<-cumsum(cw$maps[[jj[1]]])/sum(cw$maps[[jj[1]]])
    xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
    yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
    m.right<-cumsum(cw$maps[[jj[2]]])/sum(cw$maps[[jj[2]]])
    xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
    yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
    xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
    yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
    col<-colors[c(names(m.left)[length(m.left):1],names(m.right))]
    segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
             col=col,lwd=lwd,lend=lend)
  } else jj<-NULL
  for(i in 1:nrow(cw$edge)){
    if(i%in%jj==FALSE){
      maps<-cumsum(cw$maps[[i]])/sum(cw$maps[[i]])
      xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
      yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
      for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colors[names(maps)[i]],
                                       lwd=lwd,lend=lend)
    }
  }
  # plot circular lines
  for(i in 1:m+n){
    r<-R[match(i,cw$edge)]
    a1<-min(Y[which(cw$edge==i)])
    a2<-max(Y[which(cw$edge==i)])
    plotrix::draw.arc(0,0,r,a1,a2,lwd=lwd,col=colors[names(cw$maps[[match(i,cw$edge[,1])]])[1]])
  }
  # plot labels
  for(i in 1:n){
    ii<-which(cw$edge[,2]==i)
    aa<-Y[ii,2]/(2*pi)*360
    adj<-if(aa>-1&&aa<181) c(1,0.25) else c(0,0.25)
    tt<-if(aa>-1&&aa<181) paste(paste(rep(" ",offset),collapse=""),cw$tip.label[i],sep="") else paste(cw$tip.label[i],paste(rep(" ",offset),collapse=""),sep="")
    aa<-if(aa>-1&&aa<181) aa+180 else aa
    if(ftype) text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
  }
  if(setEnv){
    PP<-list(type="fan",use.edge.length=TRUE,node.pos=1,
             show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
             font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
             x.lim=xlim,y.lim=ylim,direction="rightwards",tip.color="black",
             Ntip=Ntip(cw),Nnode=cw$Nnode,edge=tree$edge,
             xx=c(x[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],x[1,1],
                  if(m>1) x[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()),
             yy=c(y[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],y[1,1],
                  if(m>1) y[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()))
    assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
  }
}


## functions plot stochastic character mapped trees
## written by Liam Revell 2011-2023

plotSimmap2<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,
                      pts=FALSE,node.numbers=FALSE,mar=NULL,add=FALSE,offset=NULL,direction="rightwards",
                      type="phylogram",setEnv=TRUE,part=if(type=="arc") 0.5 else 1.0,xlim=NULL,ylim=NULL,
                      nodes="intermediate",tips=NULL,maxY=NULL,hold=TRUE,split.vertical=FALSE,lend=2,asp=NA,
                      outline=FALSE,plot=TRUE,underscore=FALSE,arc_height=2){
  if(inherits(tree,"multiPhylo")){
    par(ask=TRUE)
    for(i in 1:length(tree)) plotSimmap2(tree[[i]],colors=colors,fsize=fsize,
                                         ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers,mar,add,offset,
                                         direction,type,setEnv,part,xlim,ylim,nodes,tips,maxY,hold,split.vertical,
                                         lend,asp,outline,plot,underscore)
  } else {
    # check tree
    if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\"")
    if(is.null(tree$maps)) stop("tree should contain mapped states on edges.")
    # check font
    ftype<-which(c("off","reg","b","i","bi")==ftype)-1
    if(!ftype) fsize=0 
    # check colors
    if(is.null(colors)){
      st<-sort(unique(unlist(sapply(tree$maps,names))))
      colors<-palette()[1:length(st)]
      names(colors)<-st
      if(length(st)>1){
        cat("no colors provided. using the following legend:\n")
        print(colors)
      }
    }
    # swap out "_" character for spaces (assumes _ is a place holder)
    if(!underscore) tree$tip.label<-gsub("_"," ",tree$tip.label)
    # get margin
    if(is.null(mar)) mar=rep(0.1,4)
    if(hold) null<-dev.hold()
    if(type=="phylogram"){
      if(direction%in%c("upwards","downwards")){
        if(outline){
          fg<-par()$fg
          par(fg="transparent")
          black<-colors
          black[]<-fg
          updownPhylogram(tree,colors=black,fsize,ftype,lwd=lwd+2,pts,
                          node.numbers,mar,add,offset,direction,setEnv,xlim,ylim,nodes,
                          tips,split.vertical,lend,asp,plot,underscore)
          par(fg=fg)
        }
        updownPhylogram(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
                        add=if(outline) TRUE else add,offset,direction,setEnv,xlim,ylim,nodes,
                        tips,split.vertical,lend,asp,plot,underscore)
      } else {
        if(outline){
          fg<-par()$fg
          par(fg="transparent")
          black<-colors
          black[]<-fg
          plotPhylogram(tree,colors=black,fsize,ftype,lwd=lwd+2,pts,
                        node.numbers,mar,add,offset,direction,setEnv,xlim,ylim,nodes,
                        tips,split.vertical,lend,asp,plot,underscore)
          par(fg=fg)
        }
        plotPhylogram(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
                      add=if(outline) TRUE else add,offset,direction,setEnv,xlim,ylim,nodes,
                      tips,split.vertical,lend,asp,plot,underscore)
      }
    } else if(type=="fan"){
      if(outline){
        fg<-par()$fg
        par(fg="transparent")
        black<-colors
        black[]<-fg
        plotFan2(tree,colors=black,fsize,ftype,lwd=lwd+2,mar,add,part,setEnv,
                 xlim,ylim,tips,maxY,lend,plot,offset)
        par(fg=fg)
      }
      plotFan2(tree,colors,fsize,ftype,lwd,mar,add=if(outline) TRUE else add,part,
               setEnv,xlim,ylim,tips,maxY,lend,plot,offset)
    } else if(type=="arc"){
      if(outline){
        fg<-par()$fg
        par(fg="transparent")
        black<-colors
        black[]<-fg
        arcPhylogram(tree,colors=black,fsize,ftype,lwd=lwd+2,mar,add,part,setEnv,
                     xlim,ylim,tips,maxY,lend,plot,offset,arc_height)
        par(fg=fg)
      }
      arcPhylogram(tree,colors,fsize,ftype,lwd,mar,add=if(outline) TRUE else add,part,
                   setEnv,xlim,ylim,tips,maxY,lend,plot,offset,arc_height)
    } else if(type=="cladogram"){
      if(outline){
        fg<-par()$fg
        par(fg="transparent")
        black<-colors
        black[]<-fg
        plotCladogram(tree,colors=black,fsize,ftype,lwd=lwd+2,mar,add,offset,
                      direction,xlim,ylim,nodes,tips,lend,asp,plot)
        par(fg=fg)
      }
      plotCladogram(tree,colors,fsize,ftype,lwd,mar,add=if(outline) TRUE else add,
                    offset,direction,xlim,ylim,nodes,tips,lend,asp,plot)
    }
    if(hold) null<-dev.flush()
  }
}

# function plots a tree; in the new version this is just a wrapper for plotSimmap
# written by Liam Revell 2012-2017, 2023
plotTree2<-function(tree,...){
  if(hasArg(color)) color<-list(...)$color
  else color<-NULL
  if(hasArg(fsize)) fsize<-list(...)$fsize
  else fsize<-1.0
  if(hasArg(ftype)) ftype<-list(...)$ftype
  else ftype<-"reg"
  if(hasArg(lwd)) lwd<-list(...)$lwd
  else lwd<-2
  if(hasArg(pts)) pts<-list(...)$pts
  else pts<-FALSE
  if(hasArg(node.numbers)) node.numbers<-list(...)$node.numbers
  else node.numbers<-FALSE
  if(hasArg(mar)) mar<-list(...)$mar
  else mar<-NULL
  if(hasArg(add)) add<-list(...)$add
  else add<-FALSE
  if(hasArg(offset)) offset<-list(...)$offset
  else offset<-NULL
  if(hasArg(type)) type<-list(...)$type
  else type<-"phylogram"
  if(hasArg(direction)) direction<-list(...)$direction
  else direction<-"rightwards"
  if(hasArg(setEnv)) setEnv<-list(...)$setEnv
  else setEnv<-TRUE
  if(hasArg(part)) part<-list(...)$part
  else part<-if(type=="arc") 0.5 else 1
  if(hasArg(xlim)) xlim<-list(...)$xlim
  else xlim<-NULL
  if(hasArg(ylim)) ylim<-list(...)$ylim
  else ylim<-NULL
  if(hasArg(nodes)) nodes<-list(...)$nodes
  else nodes<-"intermediate"
  if(hasArg(tips)) tips<-list(...)$tips
  else tips<-NULL
  if(hasArg(maxY)) maxY<-list(...)$maxY
  else maxY<-NULL
  if(hasArg(hold)) hold<-list(...)$hold
  else hold<-TRUE
  if(hasArg(lend)) lend<-list(...)$lend
  else lend<-2
  if(hasArg(asp)) asp<-list(...)$asp
  else asp<-NA
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(underscore)) underscore<-list(...)$underscore
  else underscore<-FALSE
  if(hasArg(arc_height)) arc_height<-list(...)$arc_height
  else arc_height<-2
  if(inherits(tree,"multiPhylo")){
    par(ask=TRUE)
    if(!is.null(color)) names(color)<-"1"
    for(i in 1:length(tree)) plotTree(tree[[i]],color=color,fsize=fsize,ftype=ftype,
                                      lwd=lwd,pts=pts,node.numbers=node.numbers,mar=mar,add=add,offset=offset,
                                      direction=direction,type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,
                                      nodes=nodes,tips=tips,maxY=maxY,hold=hold,lend=lend,asp=asp,plot=plot,
                                      underscore=underscore,arc_height=arc_height)
  } else {
    if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
    tree$maps<-as.list(tree$edge.length)
    for(i in 1:length(tree$maps)) names(tree$maps[[i]])<-c("1")
    if(!is.null(color)) names(color)<-"1"
    plotSimmap2(tree,colors=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,
                node.numbers=node.numbers,mar=mar,add=add,offset=offset,direction=direction,
                type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,nodes=nodes,tips=tips,maxY=maxY,
                hold=hold,lend=lend,asp=asp,plot=plot,underscore=underscore,arc_height=arc_height)
  }
}

## this is a wrapper of plotFan
## written by Liam J. Revell 2023
arcPhylogram<-function(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,
                       xlim,ylim,tips,maxY,lend,plot,offset,arc_height){
  tree<-reorder(tree,"cladewise")
  edge<-tree$edge
  edge[edge>Ntip(tree)]<-edge[edge>Ntip(tree)]+1
  edge<-rbind(Ntip(tree)+c(1,2),edge)
  Nnode<-Nnode(tree)+1
  edge.length<-c(arc_height*max(nodeHeights(tree)),tree$edge.length)
  maps<-c(vector(length=1,mode="numeric"),tree$maps)
  maps[[1]]<-setNames(edge.length[1],"NULO")
  colors<-c(setNames("transparent","NULO"),colors)
  object<-list(edge=edge,Nnode=Nnode,tip.label=tree$tip.label,
               edge.length=edge.length,maps=maps)
  class(object)<-class(tree)
  attr(object,"map.order")<-attr(object,"map.order")
  plotFan2(object,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,
           ylim,tips,maxY,lend,plot,offset)
  if(setEnv){
    PP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    ROOT<-PP$Ntip+1
    PP$Nnode<-PP$Nnode-1
    PP$edge<-PP$edge[2:nrow(PP$edge),]
    PP$edge[PP$edge>PP$Ntip]<-PP$edge[PP$edge>PP$Ntip]-1
    PP$xx<-PP$xx[-ROOT]
    PP$yy<-PP$yy[-ROOT]
    assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
  }
}



plotTree2(agam.tree2,type="arc",lwd=1,fsize=0.45,ftype="i",
          arc_height=arc_height,ylim=c(-0.1*h,1.1*(1+arc_height)*h))

