tree <- agam.tree

## create animated plot:
t=35
h<-max(nodeHeights(tree))
png(file="Amphibolurinae-%03d.png",width=1000,height=600,res=144)
for(i in 1:100){
  dev.hold()
  if(i<100){ 
    tt<-make.era.map(tree,c(0,i*h/100))
    plot(tt,colors=setNames(c('blue','transparent'),1:2),
         ftype="off",lwd=2,direction="leftwards",
         xlim=c(t,0),
         mar=c(4.1,1.1,1.1,1.1))
  } else { 
    tt<-paintSubTree(tree,Ntip(tree)+1,"1","2")
    tt$mapped.edge<-cbind(tt$mapped.edge,
                          rep(0,nrow(tt$mapped.edge)))
    plot(tt,colors=setNames('blue',1),
         ftype="off",lwd=2,direction="leftwards",
         xlim=c(t,0),
         mar=c(4.1,1.1,1.1,1.1))
  }
  axis(1)
  title(xlab="time before present")
  dev.flush()
}
dev.off()

system("magick convert -delay 10 -loop 0 *.png Amphibolurinae-anim.gif")
file.remove(list.files(pattern=".png"))
