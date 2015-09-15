library(readxl)
library(reader)
library(reshape)

main.dir <- c("/chiswick/data/ncooper/AEs","~/Documents/AEs")[2]

Ksheet.names <- c("CRP","sCD25","sIL6R","IL2","IL6","IFN","Treg","Temp","sSiglec","TregDIL","CD25")
sheets <- length(Ksheet.names)
nxt <- vector("list",sheets)

for (cc in 1:sheets) {
  nxt[[cc]] <- read_excel(cat.path(main.dir,"AEDATA2.xls"),sheet=cc)
}

# now go into each sheet and extract the proper data
# sheet names are


names(nxt) <- Ksheet.names

Kday <- c(0.001,.25,1,2,3,4,7,9,14,21,60)
Knms <- c("V0_Pre","V0_Post",paste0("V",1:9))

dose.key <- reader(cat.path(main.dir,"ID_DOSE.tab")); 
colnames(dose.key) <- c("ID","DOSE")[2]

# adverse event data
safety <- reader(cat.path(main.dir,"AEs/safety.csv"))
#colnames(safety)-->
#"ID" "StartVisit" "EndVisit" "TimeSinceMed" "Severity"  "Relatedness" "Category" "Description" 

#input raw from SS
#Knew.dat <- matrix(numeric(),ncol=3+length(Ksheet.names)*2,nrow=length(Kday)*nrow(dose.key))
#colnames(Knew.dat) <- c("ID","DAY","DOSE",paste(rep(Ksheet.names,each=2),c("raw","norm"),sep="_"))

Knew.dat <- matrix(numeric(),ncol=3+length(Ksheet.names),nrow=length(Kday)*nrow(dose.key))
colnames(Knew.dat) <- c("ID","DAY","DOSE",Ksheet.names)

rownames(Knew.dat) <- paste(rep(rownames(dose.key),each=length(Kday)),rep(Knms,length(nrow(dose.key))),sep="_")
Knew.dat <- as.data.frame(Knew.dat)
Knew.dat$ID <- rep(rownames(dose.key),each=length(Kday))
Knew.dat$DAY <- rep(Kday,nrow(dose.key))
Knew.dat$DOSE <- rep(dose.key[,1,drop=T],each=length(Kday))


####### DEFINE AE FUNCTIONS #########


# function to detect cyto-events

Diff <- function(x) { return(c(0,diff(x))) }

Z.matrix <- function(X) {
  # standardize by col means, sds
  if(!is(X)[1] %in% c("matrix","data.frame")) { stop("X should be a matrix or data.frame") }
  mns <- apply(X,2,mean,na.rm=T)
  sds <- apply(X,2,sd,na.rm=T)
  Zs <- t((t(X)-mns)/sds)
  return(Zs)
}

summary2 <- function(x) {
  ii <- summary(x)
  jj <- sd(x,na.rm=T)
  n <- length(ii)
  ii[n] <- jj
  names(ii)[n] <- "StDev"
  mn <- ii["Mean"]
  plus2 <- mn + (2*jj)
  minus2 <- mn - (2*jj)
  ii <- c(ii,minus2,plus2)
  names(ii)[(n+1):(n+2)] <- c("-2SD","+2SD")
  return(round(ii,2))
}

detect.peak <- function(x,or.dip=TRUE,thr=2) {
  if(length(x)==0) { warning("x was empty") ; return(NULL) }
  if(!length(Dim(x))==1) { stop("x must be a vector") }
  dif <-Diff(x)
  SD <- sd(dif,na.rm=T)
  if(or.dip) { pks <- which(abs(dif)>(thr*SD)) } else { pks <- which(dif>(thr*SD)) }
  if(length(pks)>0) {
    pks <- unique(sort(c(pks,pks-1)))
  }
  return(pks)
}

detect.peaks <- function(X,or.dip=TRUE,dif=TRUE,thr=2,start.end=FALSE) {
  # X is matrix, rows are samples, cols are timepoints
  if(length(Dim(X))!=2) { stop("x must be a matrix/data.frame ") }
  if(dif) {
    X <- t(apply(X,1,Diff))
  }
  Zs <- Z.matrix(X)
  pks <- which( (if(or.dip) { (abs(Zs)) } else { Zs })>thr, arr.ind=T)
  pks <- pks[order(pks[,1]),]
  if(start.end) {
    idz <- rownames(pks); Pks <- pks; rownames(pks) <- NULL
    pks <- as.data.frame(pks); pks[["id"]] <- idz
    pks[["end"]] <- pks[["start"]] <- pks$col
    pks[["dif"]] <- rep(0,nrow(pks))
    for (cc in 2:nrow(pks)) {
      if(pks$id[cc]==pks$id[cc-1]) { if(pks$col[cc]==pks$col[cc-1]+1) { pks$dif[cc] <- 1 } }
    }
    for(cc in nrow(pks):2) {
      if(pks$dif[cc]==1) {
        pks$end[cc-1] <- pks$end[cc]
        pks$start[cc] <- pks$end[cc] <- NA
      }
    }
    pks <- pks[-which(is.na(pks$start)),-6]
  }
  return(pks)
}



# function to plot events, AE or cyto

std.vec <- function(X,norm=c("raw","low3","base","base3","mean","log")[1]) {
  # standardize X according to a number of possible methods #
  norm <- tolower(norm[1])
  if(!norm %in% c("raw","low3","base","base3","mean","log","pcb")) { stop("invalid norm value") }
  if(norm=="low3") { X <- X/(mean(sort(narm(X[-2:-4]))[1:3],na.rm=T)) } # -2:-4 is to get rid of dose visits
  if(norm=="base") { X <- X/(narm(X)[1]) }
  if(norm=="base3") {  X <- X/(mean(narm(X[c(1,10,11)]),na.rm=T)) }
  if(norm=="mean") {  X <- X/(mean(narm(X),na.rm=T)) }
  if(norm=="pcb") { X <- (X-(narm(X)[1]))/(narm(X)[1]) }
  if(norm=="log") {  X <- log(X+1); X[!is.finite(X)] <- NA }
  return(X)  
}


plots <- function(id="DILT1D040",days=Kday,measure="IL2",
                  norm=c("raw","low3","base","base3","mean")[1],
                  lab.above=NULL,events=FALSE,thr=1.5,dif=TRUE,or.dip=FALSE) {
  # plot a set of ids, across a set of days for a given measure, using a set normalization scheme (or raw)
  X <- Knew.dat
  measure <- measure[1]
  #prv(X)
  vec <- std.vec(X[(X$DAY %in% days) & (X$ID %in% id[1]),measure],norm=norm)
  #print(vec)
  dat <- (vec(id=id,days=days,measure=measure[1],norm=norm))
  distr <- as.numeric(dat)
  yl <- summary(distr)[c(1,6)] # min, max
  #prv(yl)
  x.labs <- (Knms[match(days,Kday)])
  #print(x.labs)
  sub <- if(norm[1]=="raw") { "" } else { paste0("[normalised using '",norm,"']") }
  top <- paste(id,collapse=",")
  if(nchar(top)>100) { top <- "DILT1D Samples" }
  top <- paste(top,"for",measure)
  plot(days,vec,main=rbind(top,sub),xlab="day",ylab=measure,type="l",
       ylim=yl,log=if(norm=="raw") { "xy" } else { "x" },xaxt="n")
  axis(side = 1, at=days, labels=x.labs)
  if(length(id)>1) {
    for(cc in 2:length(id)) {
      vec <-  std.vec(X[(X$DAY %in% days) & (X$ID %in% id[cc]),measure],norm=norm)
      #print(vec)
      lines(days,vec,col=cc)
    }
    legend("topleft",legend=id,col=1:length(id),lwd=1,bty="n",cex=.75)
  }
  if(!is.null(lab.above)) {
    Dat <- t(apply(dat,1,function(x) { x[x!=max(x,na.rm=T)] <- 0 ; return(x) } ))
    ii <- which(Dat>lab.above,arr.ind=T)
    if(length(ii)>0) {
      for(cc in 1:nrow(ii)) {
        #cat(cc,"\n")
        text(x=days[ii[cc,2]],y=dat[ii[cc,1],ii[cc,2]],labels = rownames(dat)[ii[cc,1]],cex=0.75)
      }
    }
  }
  if(events) {
    pks <- detect.peaks(dat,dif=dif,or.dip=or.dip,thr=thr) 
    #prv(pks)
    #return(list(dat,pks))
    plot.events(dat,pks)
  }
}


plot.events <- function(X,row.col,col="blue",pch="+",cex=2) {
  if(length(Dim(row.col))!=2) { stop("row.col should have at least 1 row, and should have exactly 2 columns") }
  #prv(X,row.col)
  #print(X[row.col[,1],row.col[,2]])
  ex <- Kday[row.col[,2]]
  #prv(as.matrix(X),row.col)
  wy <- as.matrix(X)[row.col]
  #prv(ex,wy)
  points(ex,wy,pch=pch,col=col,cex=cex)
}

vec <- function(id="DILT1D040",days=Kday,measure="IL2",norm=c("raw","low3","base","base3","mean")[1]) {
  X <- Knew.dat
  measure <- measure[1]
  #prv(X)
  vec <- std.vec(X[(X$DAY %in% days) & (X$ID %in% id[1]),measure],norm=norm)
  #print(vec)
  vec <- t(as.matrix(vec))
  if(length(id)>1) {
    for(cc in 2:length(id)) {
      vec <- rbind(vec,std.vec(t(as.matrix(X[(X$DAY %in% days) & (X$ID %in% id[cc]),measure])),norm=norm))
      #print(vec)
    }
  }
  #prv(vec)
  rownames(vec) <- id
  colnames(vec) <- make.names(paste(days))
  return(vec)
}

plot.all <-  function(measures=Ksheet.names,...) {
  ii <- list(...)
  if("measure" %in% names(ii)) { stop("argument measure not allowed for plot.all") }
  
  fn <- cat.path(main.dir,"all.cyto.plots.pdf")
  pdf(fn)
  # iterate over each sheet
  for (jj in 1:length(measures)) {
    plots(...,measure=measures[jj])
  }
  dev.off()
}

code.intervals <- function(ids,starts,ends=starts,id.list=unique(ids),start.code=2) {
  ign <- is.na(starts) | is.na(ends)
  if(any(starts[ign]>ends[ign])) { stop("starts can not be greater than ends") }
  obz <- sort(unique(c(starts,ends)))
  #valz <- c(head(obz,1):tail(obz,1))
  newmat <- matrix(0,nrow=length(id.list),ncol=length(obz))
  rownames(newmat) <- id.list
  colnames(newmat) <- paste0("t",obz)
  #prv(newmat)
  for (jj in 1:length(starts)) {
    newmat[id.list==ids[jj],c(which(obz==starts[jj]):which(obz==ends[jj]))] <- 1
    newmat[id.list==ids[jj],c(which(obz==starts[jj]):which(obz==starts[jj]))] <- start.code
  }
  return(newmat)
}


event.mix <- function(pis, cis, per.plot=8, max1=TRUE) {
  if(any(dim(pis)!=dim(cis))) { stop("pis and cis must be same dimension")}
  if(any(rownames(pis)!=rownames(cis))) { stop("pis and cis didn't have the same rownames names")}
  if(any(colnames(pis)!=colnames(cis))) { warnings("pis and cis didn't have the same column names")}
  if(max1) { pis[pis>1] <- 1; cis[cis>1] <- 1 }
  x.ax <- colnames(pis)
  xs <- length(x.ax)
  n.ids <- nrow(pis)
  n.plots <- (n.ids %/% per.plot) + if((n.ids %% per.plot)>0) { 1 } else { 0 }
  prv(n.ids,n.plots)
  par(mfrow=c(n.plots,1),mar=0.1+c(5,6,4,2))
  for (pp in 1:n.plots){
    ids <- ((pp-1)*per.plot)+c(1:per.plot)
    idz <- range(ids)
    if(any(idz>n.ids)) { idz[2] <- n.ids ; ids <- ids[-which(ids>n.ids)]  }
    plot(1:xs,rep(1,xs),ylim=c(0,per.plot+2),main=paste("Adverse Events and Cytokine Events: IDs",idz[1],"to",idz[2]),
         xlab="timepoint",ylab="",col="white",bty="l",xaxt="n",yaxt="n")
    legend("top",ncol=2,legend=c("Adverse events","Cytokine events"),col=c("red","blue"),lwd=2)
    axis(1,at=1:xs,labels=x.ax)
    n <- 0
    print(ids)
    axis(2,at=1:length(ids),labels=rownames(cis)[ids],las=2)
    for(cc in ids) {
      n <- n+1
      ev.vec <- cis[cc,]
      cy.vec <- pis[cc,]
      prv(ev.vec,cy.vec)
      lines(1:xs,rep(n,xs),col="grey")
      if(any(ev.vec!=0)) {
        #lines(1:xs,n+(ev.vec*.33),col="grey")
        ind <- unique(c(which(ev.vec>0))) #,which(ev.vec[-1]>0))) # get bits to make red (events)
        XX <- 1:xs; YY <- n+(ev.vec*.44)
        if(length(ind)>0) { XX[-ind] <- NA; YY[-ind] <- NA }
        lines(XX,YY,col="red")
      }
      lines(1:xs,rep(n,xs),col="grey",lty="dotted")
      if(any(cy.vec!=0)) { 
        #lines(1:xs,n+(cy.vec*.33),col="grey",lty="dotted")
        ind <- unique(c(which(cy.vec>0))) #,which(cy.vec[-1]>0))) # get bits to make red (events)
        XX <- 1:xs; YY <- n+(cy.vec*.33)
        if(length(ind)>0) { XX[-ind] <- NA; YY[-ind] <- NA }
        lines(XX,YY,col="blue")
      }
      #  text(x=1.5,y=n,rownames(pis)[cc])
      print(n)
    }
    
  }
}


colfunc <- colorRampPalette(c("blue", "red"))


plots2 <- function(X,norm=c("lowest 3 values","raw")[1],add.mean=FALSE,xlab="day since AE",log=NULL,leg.by.dose=FALSE) {
  # plot a set of ids, across a set of days for a given measure, using a set normalization scheme (or raw)
  cn <- colnames(X)
  datz <- cn[!cn %in% c("ID","DOSE","events","visit")]
  dat <- as.matrix(X[,datz])
  distr <- as.numeric(dat)
  yl <- summary(distr)[c(1,6)] # min, max
  #prv(yl)
  x.labs <- lapply(strsplit(datz,".",fixed=T),"[",-1)
  x.labs <- sapply(x.labs,paste,collapse=".")
  meas <- unique(unlist(sapply(strsplit(datz,".",fixed=T),"[",1)))[1]
  #print(x.labs)
  id <- X$ID
  dose <- X$DOSE
  doses <- sort(unique(dose))
  if(leg.by.dose) {
    coloz <- colfunc(length(doses))
    colz <- coloz[match(dose,doses)]
  } else {
    colz <- get.distinct.cols(length(id))
  }
  top <- paste(id,collapse=",")
  if(nchar(top)>100) { top <- "DILT1D Samples" }
  top <- paste(top,"for",meas)
  #yyy <- as.numeric(x.labs); xxx <- dat[1,]; prv(xxx,yyy)
  logo <- if(is.null(log)) { if(norm=="raw") { "y" } else { "" } } else { log }
  plot(as.numeric(x.labs),dat[1,],main=rbind(top,paste("normalized using",norm)),xlab=xlab,ylab=meas,type="l",
       ylim=yl,log=logo,xaxt="n",bty="l")
  axis(side = 1, at=as.numeric(x.labs), labels=x.labs)
  if(length(id)>0) {
    for(cc in 1:length(id)) {
      lines(as.numeric(x.labs),dat[cc,],col=colz[cc])
    }
    if(leg.by.dose) {
      legend("top",legend=doses,col=coloz,lwd=1,bty="n",cex=.75,ncol=4)
    } else {
      legend("top",legend=id,col=colz,lwd=1,bty="n",cex=.75,ncol=4)
    }
  }
  if(add.mean) {
    lines(as.numeric(x.labs),apply(dat[,],2,mean,na.rm=T),col="black",lwd=3)
    legend("bottom",legend="mean",col="black",lwd=3,bty="n",cex=.75)
  }
}


std.dat <- function(X,norm="mean") {
  idz <- unique(X$ID)
  #will fail if dataset not ordered by visit within ids
  X <- X[order(X$DAY),]
  X <- X[order(X$ID),]
  cols <- 1:ncol(X); cols <- cols[-which(colnames(X) %in% c("ID","DAY","DOSE","visit","event"))]
  for(cc in 1:length(idz)) {
    XX <- X[X$ID==idz[cc],cols]
    #XX <- XX[order(XX$visit),]
    #prv(XX)
    chunk <- apply(XX,2,std.vec,norm=norm)
    X[X$ID==idz[cc],cols] <- chunk
  }
  return(X)
}


interpolate.days <- function(X,days=Kday,not.custom=FALSE) {
  # simple linear interpolation of days not measured to fill out series #
  # specific to DILT1D only!! #
  if(length(which(!is.na(X)))<2) { return(X) }
  must.use.package("zoo")
  x <- zoo(X,Sys.Date()+days)
  if(not.custom) {
    x.new <- zoo(NA, c(Sys.Date()+days[1],seq(Sys.Date()+1, to =  Sys.Date()+tail(days,1), by = "day")) )
  } else {
    x.new <- zoo(NA, c(Sys.Date()+0.001,seq(Sys.Date()+1, to =  Sys.Date()+60, by = "day")) )
  }
  x.new[time(x)] <- x 
  #  prv(x.new)
  # prv(x)
  y <- na.approx(x.new,na.rm=F)
  # prv(y)
  y <- aggregate(y, identity, head, 1)
  #  prv(y)
  out <- as.numeric(y)
  # prv(y)
  nams <- round(time(na.approx(x.new))-Sys.Date(),3)
  nams <- gsub(".000","",paste(nams),fixed=T)
  names(out) <- nams
  #prv(out)
  return(out)
}

calc.base <- function(X,rowcol=1) {
  # calculate baseline (rowcol tells whether to use columns or rows [see apply()])
  .internalF <- function(x) {
    # prv(x)
    if(length(x)!=11) { warning("x was not length 11, first obs used as baseline"); return(x[1]) }
    return(min(x[c(1,6:11)],na.rm=T))
  }
  return( if(length(Dim(X))==1) { .internalF(X) } else { apply(X,rowcol,.internalF) } )
}

calc.peak <- function(X,rng=1:length(X),wh=FALSE) {
  mns <- apply(X,2,mean,na.rm=T)
  mx <- max(mns[rng],na.rm=T)
  ii <- which(mns[rng]==mx); oo <- (1:length(X))[rng][ii]
  cat("peak found at timepoint:",ii,"\n")
  if(wh) { return(oo) } else { return(X[,oo]) } 
}

visit.to.day <- function(visit) {
  Kday[visit+1] # e.g, baselne = 0, last visit = 10 = day60 
}

# conduct a test of proportions
t.o.p <- function(p1,p2,n1,n2) {
  p <- (p1 * n1 + p2 * n2) / (n1 + n2)
  SE <- sqrt((p*(1 - p)) * ((1/n1) + (1/n2)))
  z <- (p1 - p2) / SE
  return(list(Z=z,p=Z.to.p(z)))
}   

# coefficient of variation
c.o.v <- function(x) { sd(x,na.rm=T)/abs(mean(x,na.rm=T)) }

###### END FUNCTIONS #######


## READ EXCEL DATA INTO A LONG FORMAT DATAFILE ##

for (cc in 1:length(Ksheet.names)) {
 # Header(Ksheet.names[cc])
  X <- nxt[[cc]]

  XX <- apply(X,2,function(x) { substr(x,1,6) })

  ii <- which(XX=="DILT1D",arr.ind=T)
  #prv(ii)
  for (dd in 1:nrow(ii)){
    x <- ii[dd,1] ; y <- ii[dd,2]
    nxt.id <- X[x,y]
    P1 <- as.numeric(X[(x+2):(x+(length(Kday)+1)),y])
    if(length(P1)==length(Kday)) { Knew.dat[Knew.dat$ID==nxt.id,3+(cc)] <- P1 }
   # P1 <- as.numeric(X[(x+2):(x+(length(Kday)+1)),y])
   # if(length(P1)==length(Kday)) { Knew.dat[Knew.dat$ID==nxt.id,2+(cc*2)] <- P1 }
   # P2 <- as.numeric(X[(x+2):(x+(length(Kday)+1)),y+1])
   # if(length(P2)==length(Kday)) { Knew.dat[Knew.dat$ID==nxt.id,3+(cc*2)] <- P2 }
  }
}

#now Knew.dat is ready, is what we want.


### define some vars ###

all.ids <- unique(Knew.dat$ID) # paste0("DILT1D0",c(paste0("0",1:9),10:45))

AE.ints <- with(safety,code.intervals(ID,StartVisit,EndVisit,all.ids))
Knew.dat[["visit"]] <- as.numeric(as.factor(Knew.dat$DAY))-1


## Define adverse events data.frame ##
safe <- safety[which(safety$Category %in% c("Viral","Infection") & safety$StartVisit>2),]
safe <- safe[!safe$ID %in% c("DILT1D024","DILT1D036"),]

## Add 'event' relative visit column variable to main data 
Knew.dat[["event"]] <- rep(NA,nrow(Knew.dat))
Knew.dat[["event"]] <- round(Knew.dat$visit-as.numeric(paste(safe$StartVisit[match(Knew.dat$ID,safe$ID)])))


# normal control data for IL2
normal.il2 <- c(26.67604,24.01834,37.59469,70.07785,18.09988,21.51989,22.887,14.83851,
                19.30518,8.452769,16.00315,29.41683,24.31079,28.88869,20.69033,38.0976,
                37.04048,21.23852,50.76191,27.87039,36.42799,40.21718,48.02901,15.67291,
                30.79769,18.43563,29.13794,23.923,34.32385,41.56363,27.09747,14.12415)

## calculate which samples inside outside thresholds for IL2, etc ##



il2 <- (Knew.dat$IL2)
not.dose <- !Knew.dat$DAY %in% c(.25,1,2,3,4)
il2 <- il2[not.dose]
il2.sd <- sd(normal.il2,na.rm=T)
il2.mn <- mean(normal.il2,na.rm=T)
base.il2 <- Knew.dat$IL2[Knew.dat$DAY==0.001]
t1d.il2.sd <- sd(base.il2,na.rm=T)
t1d.il2.mn <- mean(base.il2,na.rm=T)
ttt <- (base.il2); hhh <- (normal.il2)
pdf(cat.path(main.dir,"T1DvsCTRL_IL2_dists.pdf"))
plot(density(c(ttt,hhh)),col="white",xlab="IL2 raw values",main="")
lines(density(ttt),col="red");lines(density(hhh),col="blue")
legend("topright",legend=c("T1D","Healthy Controls"),col=c("red","blue"),bty="n",lwd=2)
dev.off()
cat("T1D mean: ",t1d.il2.mn,"versus healthy mean: ",il2.mn,"\n")
cat("T1D SD: ",t1d.il2.sd,"versus healthy SD: ",il2.sd,"\n")
cat("t.test of the mean difference: p=",t.test(normal.il2,base.il2)$p.value,"\n")

il2.log <- log(il2)
il2.log.mn <- mean(log(normal.il2),na.rm=T)
il2.log.sd <- sd(log(normal.il2),na.rm=T)
thresh <- 1.96

# get data for IL2 for samples/timepoints exceeding the mean+2SD (events!)
cat("\nComparison of T1D non-drug timepoints to healthy normal IL2 range\n")
cat("Log-version\n")
il2_plussers <- Knew.dat[not.dose,][(which(il2.log>(il2.log.mn+(thresh*il2.log.sd)))),]
n.up <- nrow(il2_plussers)
n.all <- nrow(Knew.dat[not.dose,])
cat(out.of(n.up,n.all),"p=",round(t.o.p(n.up/n.all,Z.to.p(thresh)/2,n.all,10^8)$p,6),"\n")

cat("Raw-version\n")
il2_plussers <- Knew.dat[not.dose,][(which(il2>(il2.mn+(thresh*il2.sd)))),]
n.up <- nrow(il2_plussers)
n.all <- nrow(Knew.dat[not.dose,])
cat(out.of(n.up,n.all),"p=",round(t.o.p(n.up/n.all,Z.to.p(thresh)/2,n.all,10^8)$p,6),"\n")




## create 'event' variable as number of days since/before adverse event ##
V <- Knew.dat
V[["event"]] <- rep(NA,nrow(V))
V[["event"]] <- round(V$DAY-as.numeric(paste(il2_plussers$DAY[match(V$ID,il2_plussers$ID)])))

X2 <- std.dat(V,norm="pcb")
X3 <- std.dat(V,norm="raw")
ii <- "Treg"  # TregDIL


### Main interesting analyses ##

samples.to.exclude <- c("DILT1D024","DILT1D040")  # add nonsense character to each to avoid excluding

for (ii in Ksheet.names) {
  Header(ii)
  Vv2 <- reshape(X2[!X2$ID %in% samples.to.exclude,],idvar="ID",v.names=ii,
               timevar="DAY",direction="wide",drop=c("event",Ksheet.names[!Ksheet.names %in% ii]))
  Vv3 <- reshape(X3[!X3$ID %in% samples.to.exclude,],idvar="ID",v.names=ii,
                 timevar="DAY",direction="wide",drop=c("event",Ksheet.names[!Ksheet.names %in% ii]))
  cln <- if(ii=="Treg") { -1 } else { 0 }
  bs <- Vv2[,4]; pk <- calc.peak(Vv2[,-1:-3],3:6) ; 
  cat(ii,"\n")
  Dose <- scale(Vv2$DOSE,scale=F); first <- Vv2[,6+cln] # 6 = day 1, 5 = 90mins
  # t-test results for IL2-effect at timepoint of detected peak (not specifying dose)
  lll <- (t.test(bs,pk,paired=T)[c(1,2,5,3)]); # peak versus baseline ignoring dose
  cat("t:",lll[[1]],"; df:",lll[[2]],"; dif:",lll[[3]],"; p:",lll[[4]],"\n")
  ## Linear regression results against 'first' [e.g. day1 or 90mins] ##
  print(round(summary(lm(first ~ Dose))[[4]],5)) # test effect of dose on day 'n' levels
  if(ii %in% c("Treg","Temp")) {
    # no 90mins data for these measures, so shift dataset appropriately
     Vv2[,6:14] <- Vv2[,5:13]; Vv2[,5]  <- apply(Vv2[,4:5],1,mean,na.rm=T)
     Vv3[,6:14] <- Vv3[,5:13]; Vv3[,5]  <- apply(Vv3[,4:5],1,mean,na.rm=T)
  }
  CoV1 <- apply(Vv3[,c(-1:-3,-5:-9)],1,c.o.v) # baseline timepoints
  CoV2 <- apply(Vv3[,c(-1:-3,-10:-14)],1,c.o.v) # dose related timepoints
  cat("median inactive coefficient of variation=",round(median(CoV1,na.rm=T),3)," versus dose:",round(median(CoV2,na.rm=T),3),"\n")
  ### Select which timepoints to add to the correlation data.frame ###
  #select1 <- c(-1:-9)[-4][-2:-4]  # base, 21, 60 - events not due to il2
  select1 <- c(-1:-9)[-4][-2:-3]  # base, 14, 21, 60 - events not due to il2
  #select <- c(-1:-14)[-4]  # base
  select2 <- 4:7 # baseline and active doses
  if(ii==Ksheet.names[1]) { cor.dat1 <- Vv2[,select1] } else { cor.dat1 <- cbind(cor.dat1,Vv2[,select1]) }
  if(ii==Ksheet.names[1]) { cor.dat2 <- Vv2[,select2] } else { cor.dat2 <- cbind(cor.dat2,Vv2[,select2]) }
  ###
}

Kday <- c(0.001,.25,1,2,3,4,7,9,14,21,60)

###### CORRELATIONS ######

# Correlation Analysis for Non-active timepoints (0,14,21,60)
tmz <- c(0.001,14,21,60)  # baseliney a
#tmz <- c(0.001,21,60)     # baseliney b
long.way1 <- reshape(cor.dat1,varying=list(paste("CRP",tmz,sep="."), paste("sCD25",tmz,sep="."),
 paste("sIL6R",tmz,sep="."),paste("IL2",tmz,sep="."), paste("IL6",tmz,sep="."),
 paste("IFN",tmz,sep="."), paste("Treg",tmz,sep="."), paste("Temp",tmz,sep="."),
 paste("sSiglec",tmz,sep="."), paste("TregDIL",tmz,sep="."), paste("CD25",tmz,sep=".")),
 direction="long")
 
corz1 <- round(cor(long.way1[,c(-1:-45,-ncol(long.way1))],use="pairwise.complete"),1)
#corz1 <- round(cor(long.way1[,c(-1:-45,-79)],use="pairwise.complete"),1)
colnames(corz1) <- gsub(".0.001","",colnames(corz1))
rownames(corz1) <- gsub(".0.001","",rownames(corz1))

corz1[!lower.tri(corz1)] <- ""
cat("\nCorrelation Analysis for Non-active timepoints (0,14,21,60)\n")
print(corz1[-1,-ncol(corz1)],quote=F)

##########

# Correlation Analysis for Active timepoints (0,90m,1,2)
tmz <- c(0.001,0.25,1,2)  # active drug
long.way2 <- reshape(cor.dat2,varying=list(paste("CRP",tmz,sep="."), paste("sCD25",tmz,sep="."),
                 paste("sIL6R",tmz,sep="."),paste("IL2",tmz,sep="."), paste("IL6",tmz,sep="."),
                 paste("IFN",tmz,sep="."), paste("Treg",tmz,sep="."), paste("Temp",tmz,sep="."),
                 paste("sSiglec",tmz,sep="."), paste("TregDIL",tmz,sep="."), paste("CD25",tmz,sep=".")),
                    direction="long")

corz2 <- round(cor(long.way2[,c(-1,-ncol(long.way))],use="pairwise.complete"),1)
colnames(corz2) <- gsub(".0.001","",colnames(corz2))
rownames(corz2) <- gsub(".0.001","",rownames(corz2))

corz2[!lower.tri(corz2)] <- ""
cat("\nCorrelation Analysis for Active timepoints (0,90m,1,2)\n")
print(corz2[-c(1,nrow(corz2)),c(-ncol(corz2),-(ncol(corz2)-1))],quote=F)

##########

if(F) {
big.cor <- cor(cor.dat,use="pairwise.complete")
diag(big.cor[grep("CD25",rownames(big.cor)),grep("CRP",colnames(big.cor))])
diag(big.cor[grep("IL2",rownames(big.cor)),grep("sCD25",colnames(big.cor))])

with(long.way,cor.test(sSiglec.0.001,CRP.0.001)) # show that correlation .3 is .0002 (passes bonf for 60 tests)
#with(cor.dat,cor.test(Temp,CRP)) # show that correlation .3 is .0002 (passes bonf for 60 tests)
}


if(F) {
###### Hi/Low counts ######

cat("\nHigh and low counts for exceeding threshold\n\n")
X2 <- std.dat(V,norm="pcb")

for (ii in Ksheet.names) {
  #Header(ii)  
  DAT <- X2 # Knew.dat
  cat(ii,"\n")
  il2 <- (DAT[,ii])
  il2n <- il2[DAT$DAY %in% c(0.001,60)]
  il2 <- il2[!DAT$DAY %in% c(0.001,.25,1,2,3)]
  if(T | ii %in% c("Treg","TregDIL")) { il2.log <- il2; il2.logn <- il2n  } else { il2.log <- log(il2); il2.logn <- log(il2n) }
  il2.log.mn <- mean(il2.logn,na.rm=T)
  il2.log.sd <- sd(il2.logn,na.rm=T)
# textogram((il2.log-il2.log.mn)/il2.log.sd)
  
  thresh <- 1.96
  cnt <- length(which(il2.log>(il2.log.mn+(thresh*il2.log.sd))))
  cat(out.of(cnt,length(il2.log)),"p=",round(t.o.p(cnt/length(il2.log),Z.to.p(thresh),length(il2.log),10^8)$p,6),"\n")
  cnt <- length(which(il2.log<(il2.log.mn-(thresh*il2.log.sd))))
  cat(out.of(cnt,length(il2.log)),"p=",round(t.o.p(cnt/length(il2.log),Z.to.p(thresh),length(il2.log),10^8)$p,6),"\n")

  #il2_plussers <- Knew.dat[!Knew.dat$DAY %in% c(0.001,.25,1,2,3),][(which(il2.log>(il2.log.mn+(1.5*il2.log.sd)))),]
}


}