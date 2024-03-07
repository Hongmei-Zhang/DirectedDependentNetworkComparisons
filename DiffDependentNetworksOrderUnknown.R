  basefn<-"10nodes10and8edges30sampCor061to10Dep" 

# Number of data sets
  starts<-1
  ends<-10

  library(pscl)
  library(psych)
  library(MCMCpack)
  library(invgamma)

  dataAll<-read.table("simuData10nodes10and8edges30samp06Corbeta1.txt",header=F)
  deltaAll<-read.table("gams10nodes10and8edges30samp06Corbeta1.txt",header=F)
  betaAll<-read.table("Beta10nodes10and8edges30samp06Corbeta1.txt",header=F)

  J<-nodes<-ncol(dataAll)
  numData<-100
  n<-sampsize<-nrow(dataAll)/numData
  
  possParents<-seq(0,(nodes-1))
  sampbetaBC<-function(betaBC,data,n,j,l,deltaBC,Sigj)
  {

    nx<-n/2
    loc<-sum(possParents[1:(j-1)])+1
    
    deltaTmp<-deltaBC
    deltaTmp[loc+l-1]<-0
    betaTmp<-betaBC*deltaTmp
    betaTmp<-betaTmp[loc:(loc+possParents[j]-1)]
    betaTmp<-betaTmp[-l]
    dataTmp<-as.matrix(data[,1:(j-1)])
    dataTmpX<-as.matrix(dataTmp[1:nx,])
    dataTmpY<-as.matrix(dataTmp[(nx+1):n,])
    
    XYj<-cbind(data[1:nx,j],data[(nx+1):n,j])
    XYl<-cbind(dataTmpX[,l],dataTmpY[,l])
    
    part1<-XYl%*%solve(Sigj)
    part2<-t(XYl)
    tmp<-0
    for (ii in 1:nx)
    {
      tmp<-tmp+part1[ii,]%*%part2[,ii]
    }
    sig2B<-1/(tmp+sig2bc^(-1))
    
    if (j>2)
    {
      dataTmpX1<-as.matrix(dataTmpX[,-l])
      dataTmpY1<-as.matrix(dataTmpY[,-l])

      betaMat<-matrix(rep(as.vector(betaTmp),nx),byrow=TRUE,ncol=(j-2))
      betaX1<-rowSums(matrix(as.numeric(betaMat),ncol=(j-2),byrow=FALSE)*dataTmpX1)
      betaY1<-rowSums(matrix(as.numeric(betaMat),ncol=(j-2),byrow=FALSE)*dataTmpY1)
      
      C1<--cbind(betaX1,betaY1)
    }
    else C1=0
    
    tmp<-0
    part2<-t(XYj+C1)
    for (ii in 1:nx)
    {
      tmp<-tmp+part1[ii,]%*%part2[,ii]
    }
    mubeta<-sig2B*tmp
    temp<-rnorm(1,mubeta,(sqrt(sig2B)))
    betaBC[loc+l-1]<-temp
  
    return(betaBC)
  }
 
  sampdeltaBC<-function(data,betaBC,deltaBC,priorSig,SigAllJ)
  {
    for (kk in 2:nodes)
    {
      loc<-sum(possParents[1:(kk-1)])+1
      dataTmp<-as.matrix(data[,1:(kk-1)])
      dataTmpX<-as.matrix(dataTmp[1:nx,])
      dataTmpY<-as.matrix(dataTmp[(nx+1):n,])
          
      for (l in 1:(kk-1))
      {
        deltaTmp<-deltaBC
        
        betaTmp1<-betaBC*deltaTmp
        betaTmp11<-betaTmp1[loc:(loc+possParents[kk]-1)]  

        betaMat<-matrix(rep(as.vector(betaTmp11),nx),byrow=TRUE,ncol=(kk-1))
        betaX<-rowSums(matrix(as.numeric(betaMat),ncol=(kk-1),byrow=FALSE)*dataTmpX)
        betaY<-rowSums(matrix(as.numeric(betaMat),ncol=(kk-1),byrow=FALSE)*dataTmpY)
        
        Rj1<-cbind(betaX,betaY)    
        XYj<-cbind(data[1:nx,kk],data[(nx+1):n,kk])
        Qj<-0
        for (i in 1:nx)
        {
          Qj<-Qj+(as.matrix(XYj[i,]-Rj1[i,]))%*%t(XYj[i,]-Rj1[i,])
        }
        Sigj<-SigAllJ[(2*(kk-1)-1):(2*(kk-1)),]
       
        invSigj<-solve(Sigj)
        part1<--tr(Qj%*%invSigj)/2
        
        deltaTmp[loc+l-1]<-1-deltaBC[loc+l-1]
        betaTmp2<-betaBC*deltaTmp
        betaTmp22<-betaTmp2[loc:(loc+possParents[kk]-1)]  
        
        betaMat<-matrix(rep(as.vector(betaTmp22),nx),byrow=TRUE,ncol=(kk-1))
        betaX<-rowSums(matrix(as.numeric(betaMat),ncol=(kk-1),byrow=FALSE)*dataTmpX)
        betaY<-rowSums(matrix(as.numeric(betaMat),ncol=(kk-1),byrow=FALSE)*dataTmpY)
        
        Rj2<-cbind(betaX,betaY)
        
        Qj<-0
        for (i in 1:nx)
        {
          Qj<-Qj+(as.matrix(XYj[i,]-Rj2[i,]))%*%t(XYj[i,]-Rj2[i,])
        }
        
        part2<--tr(Qj%*%invSigj)/2
        prop<-0.1
        if (deltaBC[loc+l-1]==1)
        {
           ratio<-1/(exp(part2-part1+betaTmp1[loc+l-1]^2/(2*priorSig))*(1-prop)/prop+1)
        }
         else 
         {
           ratio<-1/(exp(part1-part2+betaTmp2[loc+l-1]^2/(2*priorSig))*(1-prop)/prop+1)
         }
        judge<-runif(1)
        if (ratio>judge)
        {
          deltaBC[loc+l-1]<-1
        }
        else deltaBC[loc+l-1]<-0
      }
    }
    return(deltaBC)
  }
  
  sampgam<-function(gam,betaGam,data,n,j,l,deltaG,deltaBG,Sigj)
  {
    nx<-n/2
    loc<-sum(possParents[1:(j-1)])+1
    deltaTmp<-deltaG
    deltaTmp[loc+l-1]<-0
    gamTmp<-gam*deltaTmp
    gamTmp<-gamTmp[loc:(loc+possParents[j]-1)]
    gamTmp<-gamTmp[-l]
    
    deltaTmp<-deltaBG
    betaGamTmp<-as.numeric(betaGam*deltaTmp)
    betaGamTmp<-betaGamTmp[loc:(loc+possParents[j]-1)]
    
    dataTmp<-as.matrix(data[,1:(j-1)])
    dataTmpX<-as.matrix(dataTmp[1:nx,])
    dataTmpY<-as.matrix(dataTmp[(nx+1):n,])
    
    XYj<-cbind(data[1:nx,j],data[(nx+1):n,j])
    XYl<-cbind(rep(0,nx), dataTmpY[,l])
    
    part1<-XYl%*%solve(Sigj)
    part2<-t(XYl)
    tmp<-0
    for (ii in 1:nx)
    {
      tmp<-tmp+part1[ii,]%*%part2[,ii]
    }
    sig2B<-1/(tmp+sig2g^(-1))
    
    betaGamMat<-matrix(rep(as.vector(betaGamTmp),nrow(dataTmpX)),byrow=TRUE,ncol=ncol(dataTmpX))
    
    betaX11<-rowSums(betaGamMat*dataTmpX)
    C21<--cbind(betaX11,rep(0,nx))
    
    if (j>2)
    {
      dataTmpY1<-as.matrix(dataTmpY[,-l])
      
      gamMat<-matrix(rep(as.vector(gamTmp),(n-nx)),byrow=TRUE,ncol=(j-2))
      
      gamY1<-rowSums(matrix(as.numeric(gamMat),ncol=(j-2),byrow=FALSE)*dataTmpY1)
      C22<--cbind(rep(0,nx),gamY1)
    }
    else C22=0
    C2<-C21+C22
    
    tmp<-0
    part2<-t(XYj+C2)
    for (ii in 1:nx)
    {
      tmp<-tmp+part1[ii,]%*%part2[,ii]
    }
    
    mugam<- sig2B*tmp
    temp<-rnorm(1,mugam,(sqrt(sig2B)))
    gam[loc+l-1]<-temp
    return(gam)
  }
  
  sampdeltaG<-function(data,betaGam,deltaBG,gam,deltaG,priorSig,SigAllJ)
  {
    for (kk in 2:nodes)
    {
      loc<-sum(possParents[1:(kk-1)])+1
      
      for (l in 1:(kk-1))
      {	
        deltaBGTmp<-deltaBG
        deltaGTmp<-deltaG
        
        dataTmp<-as.matrix(data[,1:(kk-1)])
        dataTmpX<-as.matrix(dataTmp[1:nx,])
        dataTmpY<-as.matrix(dataTmp[(nx+1):n,])
        
        betaGamTmp<-betaGam*deltaBGTmp
        betaGamTmp<-betaGamTmp[loc:(loc+possParents[kk]-1)]  
        
        gamTmp1<-gam*deltaGTmp
        gamTmp11<-gamTmp1[loc:(loc+possParents[kk]-1)]  
        
        betaGamMat<-matrix(rep(as.vector(betaGamTmp),nx),byrow=TRUE,ncol=(kk-1))
        betaX<-rowSums(matrix(as.numeric(betaGamMat),ncol=(kk-1),byrow=FALSE)*dataTmpX)
        
        gamMat<-matrix(rep(as.vector(gamTmp11),(n-nx)),byrow=TRUE,ncol=(kk-1))
        betaY<-rowSums(matrix(as.numeric(gamMat),ncol=(kk-1),byrow=FALSE)*dataTmpY)
        
        Rj1<-cbind(betaX,betaY)    
        XYj<-cbind(data[1:nx,kk],data[(nx+1):n,kk])
        Qj<-0
        for (i in 1:nx)
        {
          Qj<-Qj+(as.matrix(XYj[i,]-Rj1[i,]))%*%t(XYj[i,]-Rj1[i,])
        }
        Sigj<-SigAllJ[(2*(kk-1)-1):(2*(kk-1)),]
        invSigj<-solve(Sigj)
        
        part1<--tr(Qj%*%invSigj)/2
        
        deltaGTmp[loc+l-1]<-1-deltaG[loc+l-1]
        gamTmp2<-gam*deltaGTmp
        gamTmp22<-gamTmp2[loc:(loc+possParents[kk]-1)]  

        gamMat<-matrix(rep(as.vector(gamTmp22),(n-nx)),byrow=TRUE,ncol=(kk-1))
        betaY<-rowSums(matrix(as.numeric(gamMat),ncol=(kk-1),byrow=FALSE)*dataTmpY)
        
        Rj2<-cbind(betaX,betaY)
        
        Qj<-0
        for (i in 1:nx)
        {
          Qj<-Qj+(as.matrix(XYj[i,]-Rj2[i,]))%*%t(XYj[i,]-Rj2[i,])
        }
        
        part2<--tr(Qj%*%invSigj)/2
        prop<-0.1
        if (deltaG[loc+l-1]==1)
        {
          ratio<-1/(exp(part2-part1+gamTmp1[loc+l-1]^2/(2*priorSig))*(1-prop)/prop+1)
        }
        else
        {
          ratio<-1/(exp(part1-part2+gamTmp2[loc+l-1]^2/(2*priorSig))*(1-prop)/prop+1)
        }
        judge<-runif(1)
        if (ratio>judge)
        {
          deltaG[loc+l-1]<-1
        }
        else deltaG[loc+l-1]<-0
      }
    }
    return(deltaG)
  }
  
  sampbetaGam<-function(gam,betaGam,data,n,j,l,deltaG,deltaBG,Sigj)
  {
    nx<-n/2
    loc<-sum(possParents[1:(j-1)])+1
    deltaTmp<-deltaBG
    deltaTmp[loc+l-1]<-0
    betaGamTmp<-betaGam*deltaTmp
    betaGamTmp<-betaGamTmp[loc:(loc+possParents[j]-1)]
    betaGamTmp<-betaGamTmp[-l]
    
    deltaTmp<-deltaG
    gamTmp<-as.numeric(gam*deltaTmp)
    gamTmp<-gamTmp[loc:(loc+possParents[j]-1)]
    
    dataTmp<-as.matrix(data[,1:(j-1)])
    dataTmpX<-as.matrix(dataTmp[1:nx,])
    dataTmpY<-as.matrix(dataTmp[(nx+1):n,])
    
    XYj<-cbind(data[1:nx,j],data[(nx+1):n,j])
    XYl<-cbind(dataTmpX[,l],rep(0,nx))
    
    part1<-XYl%*%solve(Sigj)
    part2<-t(XYl)
    tmp<-0
    for (ii in 1:nx)
    {
      tmp<-tmp+part1[ii,]%*%part2[,ii]
    }
    sig2B<-1/(tmp+sig2bg^(-1))
    
    gamMat<-matrix(rep(as.vector(gamTmp),(n-nx)),byrow=TRUE,ncol=(j-1))
    betaY11<-rowSums(gamMat*dataTmpY)
    C21<--cbind(rep(0,nx),betaY11)
    
    if (j>2)
    {
      dataTmpX1<-as.matrix(dataTmpX[,-l])
      
      betaGamMat<-matrix(rep(as.vector(betaGamTmp),nx),byrow=TRUE,ncol=(j-2))
      
      betaGamX1<-rowSums(matrix(as.numeric(betaGamMat),ncol=(j-2),byrow=FALSE)*dataTmpX1)
      C22<--cbind(betaGamX1,rep(0,nx))
    }
    else C22=0
    C2<-C21+C22
    
    tmp<-0
    part2<-t(XYj+C2)
    for (ii in 1:nx)
    {
      tmp<-tmp+part1[ii,]%*%part2[,ii]
    }
    
    mugam<- sig2B*tmp

    temp<-rnorm(1,mugam,(sqrt(sig2B)))
    betaGam[loc+l-1]<-temp
    return(betaGam)
  }
  
  sampdeltaBG<-function(data,betaGam,deltaBG,gam,deltaG,priorSig,SigAllJ)
  {
    for (kk in 2:nodes)
    {
      loc<-sum(possParents[1:(kk-1)])+1
      
      for (l in 1:(kk-1))
      {	
        deltaBGTmp<-deltaBG
        deltaGTmp<-deltaG
        
        dataTmp<-as.matrix(data[,1:(kk-1)])
        dataTmpX<-as.matrix(dataTmp[1:nx,])
        dataTmpY<-as.matrix(dataTmp[(nx+1):n,])
        
        betaGamTmp1<-betaGam*deltaBGTmp
        betaGamTmp11<-betaGamTmp1[loc:(loc+possParents[kk]-1)]  
        gamTmp<-gam*deltaGTmp
        gamTmp<-gamTmp[loc:(loc+possParents[kk]-1)]  
        
        betaGamMat<-matrix(rep(as.vector(betaGamTmp11),nx),byrow=TRUE,ncol=(kk-1))
        betaX<-rowSums(matrix(as.numeric(betaGamMat),ncol=(kk-1),byrow=FALSE)*dataTmpX)
        
        gamMat<-matrix(rep(as.vector(gamTmp),(n-nx)),byrow=TRUE,ncol=(kk-1))
        betaY<-rowSums(matrix(as.numeric(gamMat),ncol=(kk-1),byrow=FALSE)*dataTmpY)
        
        Rj1<-cbind(betaX,betaY)    
        XYj<-cbind(data[1:nx,kk],data[(nx+1):n,kk])
        Qj<-0
        for (i in 1:nx)
        {
          Qj<-Qj+(as.matrix(XYj[i,]-Rj1[i,]))%*%t(XYj[i,]-Rj1[i,])
        }
        Sigj<-SigAllJ[(2*(kk-1)-1):(2*(kk-1)),]
        invSigj<-solve(Sigj)
        
        part1<--tr(Qj%*%invSigj)/2
  
        deltaBGTmp[loc+l-1]<-1-deltaBG[loc+l-1]
        betaGamTmp2<-betaGam*deltaBGTmp
        betaGamTmp22<-betaGamTmp2[loc:(loc+possParents[kk]-1)]  
        
        betaGamMat<-matrix(rep(as.vector(betaGamTmp22),nx),byrow=TRUE,ncol=(kk-1))
        betaX<-rowSums(matrix(as.numeric(betaGamMat),ncol=(kk-1),byrow=FALSE)*dataTmpX)
        Rj2<-cbind(betaX,betaY)
        
        Qj<-0
        for (i in 1:nx)
        {
          Qj<-Qj+(as.matrix(XYj[i,]-Rj2[i,]))%*%t(XYj[i,]-Rj2[i,])
        }
        part2<--tr(Qj%*%invSigj)/2
        prop<-0.1
        if (deltaBG[loc+l-1]==1)
        {
          ratio<-1/(exp(part2-part1+betaGamTmp1[loc+l-1]^2/(2*priorSig))*(1-prop)/prop+1)
        }
        else
        {
          ratio<-1/(exp(part1-part2+betaGamTmp2[loc+l-1]^2/(2*priorSig))*(1-prop)/prop+1)
        }

        judge<-runif(1)
        if (ratio>judge)
        {
          deltaBG[loc+l-1]<-1
        }
        else deltaBG[loc+l-1]<-0
      }
    }
    return(deltaBG)
  }
  
  sampEta<-function(data,eta,betaBC,deltaBC,gam,deltaG,betaGam,deltaBG,a1,a2,n,J)
  {
    nx<-n/2
    sumCom<-0
    sumDiff<-0
    SigAllJ<-NULL
    betaDiff<-betaGam*deltaBG
    gamDiff<-gam*deltaG
    
    for (j in 2:J)
    {
      loc<-sum(possParents[1:(j-1)])+1
      XYj<-cbind(data[1:nx,j],data[(nx+1):n,j])
      datajX<-as.matrix(data[1:nx,1:(j-1)])
      datajY<-as.matrix(data[(nx+1):n,1:(j-1)])
      
      betaTmp<-betaDiff[loc:(loc+possParents[j]-1)]
      gamTmp<-gamDiff[loc:(loc+possParents[j]-1)]
      
      betaMat<-matrix(rep(as.vector(betaTmp),nx),byrow=TRUE,ncol=(j-1))
      betaGamX<-rowSums(matrix(as.numeric(betaMat),ncol=(j-1),byrow=FALSE)*datajX)
      C20<-cbind(betaGamX,rep(0,nx))
      
      gamMat<-matrix(rep(as.vector(gamTmp),(n-nx)),byrow=TRUE,ncol=(j-1))
      gamY11<-rowSums(matrix(as.numeric(gamMat),ncol=(j-1),byrow=FALSE)*datajY)
      C21<-cbind(rep(0,nx),gamY11)
      
      Rj<-C20+C21
      Qj<-0
  
      for (i in 1:nx)
      {
        Qj<-Qj+as.matrix(XYj[i,]-Rj[i,])%*%t(XYj[i,]-Rj[i,])
      }
      S<-2*nu*diag(c(1/a1[j],1/a2[j]))+Qj
  
      degre<-nu+nx+1
      Sigj<-riwish(degre, S)
      SigAllJ<-rbind(SigAllJ,Sigj)
      
      sumj<-cbind(rowSums(matrix(as.numeric(betaMat),ncol=(j-1),byrow=FALSE)*datajX),
                    rowSums(matrix(as.numeric(gamMat),ncol=(j-1),byrow=FALSE)*datajY))
      
      diff<-(XYj-sumj)%*%solve(Sigj)
      tmpsumVec<-NULL
      for (i in 1:nx)
      {
        tmp<-t(diff[i,])%*%((XYj-sumj)[i,])
        tmpsumVec<-c(tmpsumVec,tmp)
      }
      
      sumDiff<-sumDiff-sum(tmpsumVec)/2-nx/2*log(det(Sigj))
  
      betaCom<-betaBC*deltaBC
      betaTmp<-betaCom[loc:(loc+possParents[j]-1)]
  
      betaMat<-matrix(rep(as.vector(betaTmp),nx),byrow=TRUE,ncol=(j-1))
      sumj<-cbind(rowSums(matrix(as.numeric(betaMat),ncol=(j-1),byrow=FALSE)*datajX),
                    rowSums(matrix(as.numeric(betaMat),ncol=(j-1),byrow=FALSE)*datajY))
      
      diff<-(XYj-sumj)%*%solve(Sigj)
      tmpsumVec<-NULL
      for (i in 1:nx)
      {
        tmp<-t(diff[i,])%*%((XYj-sumj)[i,])
        tmpsumVec<-c(tmpsumVec,tmp)
      }
      sumCom<-sumCom-sum(tmpsumVec)/2-nx/2*log(det(Sigj))
    }
    sumDeltaCom<-sum(deltaBC)
    sumDeltaY<-sum(deltaG)
    sumDeltaX<-sum(deltaBG)
    
    r<-1+exp(sumDiff-sumCom+0.5*log(nx)*(sumDeltaCom-sumDeltaX-sumDeltaY))
    r<-r^(-1)
    eta<-1-sum(round(r,digits=1)>=0.5)
    return(list(eta=eta,SigAllJ=SigAllJ,A=sumCom,B=sumDiff))
  }
  
  MHStep<-function(order,k,iii,nodes,etaSampled,betaGamSampled,betaGamSampledUnique,deltaBGSampled,deltaBGSampledUnique,
                   gamSampled,gamSampledUnique,deltaGSampled,deltaGSampledUnique,
                   betaBCSampled,betaBCSampledique,deltaBCSampled,deltaBCSampledUnique,
                   SigAllJUnique,nx,data,probTmp)
  {
    if (iii==((K-k)*(B+N)+1))
    {
      orderTmp<-order[k,1,]
      deltaBG<-deltaBGSampled[k,1,2:ncol(deltaBGSampled[k,,])]
      betaGam<-betaGamSampled[k,1,2:ncol(betaGamSampled[k,,])]
      
      deltaG<-deltaGSampled[k,1,2:ncol(deltaGSampled[k,,])]
      gam<-gamSampled[k,1,2:ncol(gamSampled[k,,])]
      
      betaBC<-betaBCSampled[k,1,2:ncol(betaBCSampled[k,,])]
      
      SigAllJ<-SigAllJUnique[1:(2*(J-1)),]
    }
    else 
    {
      orderTmp<-order[k,(iii-1),]
      deltaBG<-deltaBGSampled[k,(iii-1),2:ncol(deltaBGSampled[k,,])]
      betaGam<-betaGamSampled[k,(iii-1),2:ncol(betaGamSampled[k,,])]
      
      deltaG<-deltaGSampled[k,(iii-1),2:ncol(deltaGSampled[k,,])]
      gam<-gamSampled[k,(iii-1),2:ncol(gamSampled[k,,])]
      
      betaBC<-betaBCSampled[k,(iii-1),2:ncol(betaBCSampled[k,,])]
      
      SigLoc<-order[k,(iii-1),1]
  
      SigAllJ<-SigAllJUnique[(2*(J-1)*(SigLoc-1)+1):(2*(J-1)*SigLoc),]
    }
    dataTmp<-data[,orderTmp[2:length(orderTmp)]]		
    
    probO<-OrderProb(deltaBG,betaGam,deltaG,gam,betaBC,n,SigAllJ,dataTmp)
    probOld<-probO$prob
  
  	deltaBGOld<-probO$avgdeltaBG
    betaGamOld<-probO$avgbetaGam
    
    deltaGOld<-probO$avgdeltaG
    gamOld<-probO$avggam
    
    deltaBCOld<-probO$avgdeltaBC
    betaBCOld<-probO$avgbetaBC
    
    SigAllJOld<-probO$SigAllJ
    etaOld<-probO$eta

     loc<-sort(sample(seq(1,nodes), 2, replace = FALSE, prob = NULL))

     orderTmp1<-orderTmp
     tmp<-orderTmp[(loc[1]+1)]
     for (mm in 2:length(loc))
     {
       orderTmp1[(loc[mm-1]+1)]<-orderTmp[(loc[mm]+1)]
     }
     orderTmp1[(loc[mm]+1)]<-tmp

    if (iii==((K-k)*(B+N)+20))
    {
      orderTmp1[2:length(orderTmp1)]<-seq(1,nodes)
    }

     dataTmp<-data[,orderTmp1[2:length(orderTmp1)]]

  for (kk in 1:K)
  {
    if (kk==1)
    {
      uniqueOrder<-unique(order[kk,1:(iii-1),1])
    }
    else
    {
      uniqueOrder<-unique(c(uniqueOrder,unique(order[kk,1:(iii-1),1])))
    }
  }
     NAloc<-which(is.na(uniqueOrder)==1)
     if (length(NAloc!=0))	uniqueOrder<-uniqueOrder[-NAloc]

     for (kk in 1:length(uniqueOrder))
     {
       check<-0
       for (t in 1:K)
       {
         loc<-max(which(order[t,,1]==uniqueOrder[kk]))
         if (loc>0)
         {
           diff<-orderTmp1[2:length(orderTmp1)]-order[t,loc,2:(nodes+1)]
           if ((sum(abs(diff)))==0)
           {
             check<-1
             locPick<-loc
             kChose<-t
             break
           }
         }
         if (check==1) break
       }
       if (check==1) break
     }# end for (kk in 1:length(uniqueOrder))
     if (check==1)
     {
       betaGam<-betaGamSampledUnique[order[kChose,loc,1],]
       deltaBG<-deltaBGSampledUnique[order[kChose,loc,1],]

       gam<-gamSampledUnique[order[kChose,loc,1],]
       deltaG<-deltaGSampledUnique[order[kChose,loc,1],]

       betaBC<-betaBCSampledUnique[order[kChose,loc,1],]
       deltaBC<-deltaBCSampledUnique[order[kChose,loc,1],]

       SigLoc<-order[kChose,loc,1]

       SigAllJ<-SigAllJUnique[(2*(J-1)*(SigLoc-1)+1):(2*(J-1)*SigLoc),]
     }
     else
     {
       # Time 1 network
       deltaBG<-NULL
       for (i in 2:J)
       {
         if (i==2)
         {
           deltaBG<-rbinom((i-1),1,prob=0.5)
         }
         else
         {
           deltaBG<-c(deltaBG,rbinom((i-1),1,prob=0.5))
         }
       } # end for (i in 2:nodes)
       betaGam<-deltaBG

       # Time 2 network
       deltaG<-NULL
       for (i in 2:J)
       {
         if (i==2)
         {
           deltaG<-rbinom((i-1),1,prob=0.5)
         }
         else
         {
           deltaG<-c(deltaG,rbinom((i-1),1,prob=0.5))
         }
       } # end for (i in 2:nodes)
       gam<-deltaG

       # identical network
       deltaBC<-as.numeric(deltaBG|deltaG)
       betaBC<-deltaBC

       SigAllJ<-matrix(rep(SigjInit,(J-1)),ncol=2,byrow=T)
     } # end else

     probO<-OrderProb(deltaBG,betaGam,deltaG,gam,betaBC,n,SigAllJ,dataTmp)

     lpprobk<- -max(-probO$prob,H[k])/T[k]
     lpprobOldk<- -max(-probOld,H[k])/T[k]
     ratio<-min(1,exp(lpprobk-lpprobOldk))
  if (ratio>runif(1))
  {
    # assign the new order to the iteration iii.
    order[k,iii,]<-orderTmp1
    if (check==1)
    {
      order[k,iii,1]<-order[kChose,locPick,1]
      betaGamSampledUnique[order[kChose,locPick,1],]<-probO$avgbetaGam
      deltaBGSampledUnique[order[kChose,locPick,1],]<-probO$avgdeltaBG

      gamSampledUnique[order[kChose,locPick,1],]<-probO$avggam
      deltaGSampledUnique[order[kChose,locPick,1],]<-probO$avgdeltaG

      betaBCSampledUnique[order[kChose,locPick,1],]<-probO$avgbetaBC
      deltaBCSampledUnique[order[kChose,locPick,1],]<-probO$avgdeltaBC

      SigLoc<-order[kChose,locPick,1]

      SigAllJUnique[(2*(J-1)*(SigLoc-1)+1):(2*(J-1)*SigLoc),]<-probO$SigAllJ

      betaGamSampled[k,iii,1]<-order[k,iii,1]
      betaGamSampled[k,iii,2:ncol(betaGamSampled[k,,])]<-probO$avgbetaGam
      deltaBGSampled[k,iii,1]<-order[k,iii,1]
      deltaBGSampled[k,iii,2:ncol(deltaBGSampled[k,,])]<-probO$avgdeltaBG

      gamSampled[k,iii,1]<-order[k,iii,1]
      gamSampled[k,iii,2:ncol(gamSampled[k,,])]<-probO$avggam
      deltaGSampled[k,iii,1]<-order[k,iii,1]
      deltaGSampled[k,iii,2:ncol(deltaGSampled[k,,])]<-probO$avgdeltaG

      betaBCSampled[k,iii,1]<-order[k,iii,1]
      betaBCSampled[k,iii,2:ncol(betaBCSampled[k,,])]<-probO$avgbetaBC
      deltaBCSampled[k,iii,1]<-order[k,iii,1]
      deltaBCSampled[k,iii,2:ncol(deltaBCSampled[k,,])]<-probO$avgdeltaBC

      etaSampled[k,iii,1]<-order[k,iii,1]
      etaSampled[k,iii,2]<-probO$eta

      probTmp<-probO$prob
    }
    else
    {
      order[k,iii,1]<-max(uniqueOrder)+1
      betaGamSampledUnique<-rbind(betaGamSampledUnique,probO$avgbetaGam)
      deltaBGSampledUnique<-rbind(deltaBGSampledUnique,probO$avgdeltaBG)

      gamSampledUnique<-rbind(gamSampledUnique,probO$avggam)
      deltaGSampledUnique<-rbind(deltaGSampledUnique,probO$avgdeltaG)

      betaBCSampledUnique<-rbind(betaBCSampledUnique,probO$avgbetaBC)
      deltaBCSampledUnique<-rbind(deltaBCSampledUnique,probO$avgdeltaBC)

      SigAllJUnique<-rbind(SigAllJUnique,probO$SigAllJ)

      betaGamSampled[k,iii,1]<-order[k,iii,1]
      betaGamSampled[k,iii,2:ncol(betaGamSampled[k,,])]<-probO$avgbetaGam
      deltaBGSampled[k,iii,1]<-order[k,iii,1]
      deltaBGSampled[k,iii,2:ncol(deltaBGSampled[k,,])]<-probO$avgdeltaBG

      gamSampled[k,iii,1]<-order[k,iii,1]
      gamSampled[k,iii,2:ncol(gamSampled[k,,])]<-probO$avggam
      deltaGSampled[k,iii,1]<-order[k,iii,1]
      deltaGSampled[k,iii,2:ncol(deltaGSampled[k,,])]<-probO$avgdeltaG

      betaBCSampled[k,iii,1]<-order[k,iii,1]
      betaBCSampled[k,iii,2:ncol(betaBCSampled[k,,])]<-probO$avgbetaBC
      deltaBCSampled[k,iii,1]<-order[k,iii,1]
      deltaBCSampled[k,iii,2:ncol(deltaBCSampled[k,,])]<-probO$avgdeltaBC

      etaSampled[k,iii,1]<-order[k,iii,1]
      etaSampled[k,iii,2]<-probO$eta

      probTmp<-probO$prob
    }
  } # end if (ratio>runif(1))
  else
  {
      if (iii==((K-k)*(B+N)+1))
      {
        order[k,iii,]<-order[k,1,]
      }
      else	
      {
        order[k,iii,]<-order[k,(iii-1),]
      }
      betaGamSampled[k,iii,1]<-order[k,iii,1]
      betaGamSampled[k,iii,2:ncol(betaGamSampled[k,,])]<-betaGamOld
      deltaBGSampled[k,iii,1]<-order[k,iii,1]
      deltaBGSampled[k,iii,2:ncol(deltaBGSampled[k,,])]<-deltaBGOld
      
      gamSampled[k,iii,1]<-order[k,iii,1]
      gamSampled[k,iii,2:ncol(gamSampled[k,,])]<-gamOld
      deltaGSampled[k,iii,1]<-order[k,iii,1]
      deltaGSampled[k,iii,2:ncol(deltaGSampled[k,,])]<-deltaGOld
      
      betaBCSampled[k,iii,1]<-order[k,iii,1]
      betaBCSampled[k,iii,2:ncol(betaBCSampled[k,,])]<-betaBCOld
      deltaBCSampled[k,iii,1]<-order[k,iii,1]
      deltaBCSampled[k,iii,2:ncol(deltaBCSampled[k,,])]<-deltaBCOld
      
      etaSampled[k,iii,1]<-order[k,iii,1]
      etaSampled[k,iii,2]<-etaOld
      
      probTmp<-probOld
     }
    return(list(order=order,betaGamSampled=betaGamSampled,betaGamSampledUnique=betaGamSampledUnique,deltaBGSampled=deltaBGSampled,
                deltaBGSampledUnique=deltaBGSampledUnique,
                gamSampled=gamSampled,gamSampledUnique=gamSampledUnique,deltaGSampled=deltaGSampled,
                deltaGSampledUnique=deltaGSampledUnique,
                betaBCSampled=betaBCSampled,betaBCSampledUnique=betaBCSampledUnique,deltaBCSampled=deltaBCSampled,
                deltaBCSampledUnique=deltaBCSampledUnique,
                etaSampled=etaSampled,SigAllJUnique=SigAllJUnique,
                probTmp=probTmp))
  }
  
  OrderProb<-function(deltaBG,betaGam,deltaG,gam,betaBC,n,SigAllJ,data)
  {
    deltaBGSum<-0
    sumbetaGam<-0
    
    deltaGSum<-0
    sumgam<-0
    
    deltaBCSum<-0
    sumbetaBC<-0
    
    a1sum<-a2sum<-0
    
    eta<-0
    
    dataX<-data[1:nx,]
    dataY<-data[(nx+1):n,]
    nx<-n/2
    a1<-a2<-rep(0,J)
    
    for (ii in 1:loops)
    {
      deltaBG<-sampdeltaBG(data,betaGam,deltaBG,gam,deltaG,sig2bg,SigAllJ)

      for (i in 2:J)
      {
        loc<-sum(possParents[1:(i-1)])+1
        nonZero<-which(deltaBG[loc:(loc+possParents[i]-1)]!=0)
        for (k1 in nonZero)
        {
          betaGam<-sampbetaGam(gam,betaGam,data,n,i,k1,deltaG,deltaBG,SigAllJ[(2*(i-1)-1):(2*(i-1)),])
        }
      }
      deltaG<-sampdeltaG(data,betaGam,deltaBG,gam,deltaG,sig2g,SigAllJ)
      for (i in 2:J)
      {
        loc<-sum(possParents[1:(i-1)])+1
        nonZero<-which(deltaG[loc:(loc+possParents[i]-1)]!=0)
        for (k1 in nonZero)
        {
          gam<-sampgam(gam,betaGam,data,n,i,k1,deltaG,deltaBG,SigAllJ[(2*(i-1)-1):(2*(i-1)),])
        }
      }
      deltaBC<-as.numeric(deltaBG|deltaG)
      deltaBC<-sampdeltaBC(data,betaBC,deltaBC,sig2bc,SigAllJ)
      for (i in 2:J)
      {
        loc<-sum(possParents[1:(i-1)])+1
        nonZero<-which(deltaBC[loc:(loc+possParents[i]-1)]!=0)
        for (k1 in nonZero)
        {
          betaBC<-sampbetaBC(betaBC,data,n,i,k1,deltaBC,SigAllJ[(2*(i-1)-1):(2*(i-1)),])
        }
      }
      for (i in 2:J)
      {
        shape1<-shape2<-(nu+2)/2
        
        scale1<-nu*SigAllJ[(2*(i-1)-1),1]+alp1
        a1[i]<-rigamma(1,shape1,scale1)
        
        scale2<-nu*SigAllJ[(2*(i-1)),2]+alp2
        a2[i]<-rigamma(1,shape2,scale2)
      }
      
      if (ii>(loops/2)) 
      {
        deltaBGSum<-deltaBGSum+deltaBG
        sumbetaGam<-sumbetaGam+betaGam
        
        deltaGSum<-deltaGSum+deltaG
        sumgam<-sumgam+gam
        
        deltaBCSum<-deltaBCSum+deltaBC
        sumbetaBC<-sumbetaBC+betaBC
        
        a1sum<-a1sum+a1
        a2sum<-a2sum+a2
      }
    } # end number of loops
    avgdeltaBG<-round(deltaBGSum/(loops/2))
    avgbetaGam<-sumbetaGam/(loops/2)*avgdeltaBG

    avgdeltaG<-round(deltaGSum/(loops/2))
    avggam<-sumgam/(loops/2)*avgdeltaG

    avgdeltaBC<-round(deltaBCSum/(loops/2))
    avgbetaBC<-sumbetaBC/(loops/2)*avgdeltaBC

    avga1<-a1sum/(loops/2)
    avga2<-a2sum/(loops/2)
    etaSig<-sampEta(data,eta,avgbetaBC,avgdeltaBC,avggam,avgdeltaG,avgbetaGam,avgdeltaBG,avga1,avga2,n,J)
    eta<-etaSig$eta
    SigAllJ<-etaSig$SigAllJ
    # B is for differential networks and A is for identical networks
    # eta=0 means identical networks
    if (eta==0)
    {
      prob<-etaSig$A
    }
    else
    {
      prob<-etaSig$B
    }
    return(list(eta=eta,avgdeltaBG=avgdeltaBG,avgbetaGam=avgbetaGam,avgdeltaG=avgdeltaG,avggam=avggam,
                avgdeltaBC=avgdeltaBC,avgbetaBC=avgbetaBC,avga1=avga1,avga2=avga2,
                SigAllJ=SigAllJ,prob=prob))
  }
  
  # This function is only used to get a rough estimate of graph probability, used to set H values. 
  probGraphBC<-function(deltaBC,betaBC,SigAllJ,nx,data)
  {
    prob<-0
    
    for (kk in 2:J)
    {
      datajX<-data[1:nx,1:(kk-1)]
      datajY<-data[(nx+1):n,1:(kk-1)]
      XYj<-cbind(data[1:nx,kk],data[(nx+1):n,kk])
      
      loc<-sum(possParents[1:(kk-1)])+1
      deltaBCTmp<-deltaBC
      betaBCTmp<-betaBC[loc:(loc+possParents[kk]-1)]*deltaBC[loc:(loc+possParents[kk]-1)]
  
      betaBCX<-rowSums(as.matrix(betaBCTmp*datajX))
      betaBCY<-rowSums(as.matrix(betaBCTmp*datajY))
      Rj<-cbind(betaBCX,betaBCY)
      
      Qj<-0
      for (i in 1:nx)
      {
        Qj<-Qj+as.matrix(XYj[i,]-Rj[i,])%*%t(XYj[i,]-Rj[i,])
      }
      Sigj<-SigAllJ[(2*kk-3):(2*kk-2),]
      invSigj<-solve(Sigj)
      
      prob<-prob-nx/2*log(det(Sigj))-tr(Qj%*%invSigj)/2
    }
    return(prob)
  }
  
  
  start_time <- Sys.time()
  
  # number of MCMC iterations
  loops<-10
  
  # K is the number of chains
  K<-5
  # N is the number of iterations to be used to collect results
  B<-20
  N<-20
  orderloops<-(B+N)*K+B+N
  
  choice<-rep(0,numData)
  pee<-0.1
  
  accurdeltaBG<-accurdeltaG<-accurdeltaBC<-rep(0,numData)
  fpdeltaBG<-fpdeltaG<-fpdeltaBC<-rep(0,numData)
  tpdeltaBG<-tpdeltaG<-tpdeltaBC<-rep(0,numData)
  
  avgEtaEst<-rep(0,numData)
  
  for (jj in starts:ends)
  {
    cat("data ",jj,"\n")
    set.seed(12345*jj)
    
    data<-dataAll[(1+(jj-1)*sampsize):(jj*sampsize),]
    n<-nrow(data)
    nx<-n/2
    
    #  this is for separated networks
    stoppoint<-nodes*(nodes-1)/2
    deltaBGTrue<-deltaAll[jj,1:stoppoint]
    deltaGTrue<-deltaAll[jj,(stoppoint+1):(2*stoppoint)]
    
    truebetaGam<-betaAll[jj,1:stoppoint]
    truegam<-betaAll[jj,(stoppoint+1):(2*stoppoint)]
    
    #  this will be used for combined networks
    deltaBCTrue<-deltaAll[jj,1:stoppoint]
    truebetaBC<-betaAll[jj,1:stoppoint]

    #  Initial values
    
    # pre-specified large variances in the mixed distribution of beta_ij
    sig2bc<-sig2bg<-sig2g<-100
    
    #initial values
    deltaBC<-NULL
    for (i in 2:nodes)
    {
      if (i==2)
      {
        deltaBC<-rbinom((i-1),1,prob=0.5)
      }
      else
      {
        deltaBC<-c(deltaBC,rbinom((i-1),1,prob=0.5))
      }
    }

    # X network
    deltaBG<-NULL
    for (i in 2:nodes)
    {
      if (i==2)
      {
        deltaBG<-rbinom((i-1),1,prob=0.5)
      }
      else
      {
        deltaBG<-c(deltaBG,rbinom((i-1),1,prob=0.5))
      }
    }
    
    # Y network
    deltaG<-NULL
    for (i in 2:nodes)
    {
      if (i==2)
      {
        deltaG<-rbinom((i-1),1,prob=0.5)
      }
      else
      {
        deltaG<-c(deltaG,rbinom((i-1),1,prob=0.5))
      }
    }
    
    # initinal values for beta
    betaBC<-rep(1,length(deltaBC))
    betaGam<-rep(1,length(deltaBG))
    gam<-rep(1,length(deltaG))
    
    # parameter in the prior of Sigma (the covariance matrix in the linear mixed model)
    nu<-2
    # parameters in the hyper-prior of a_j1 and a_j2
    alp1<-alp2<-0.0001
    
    deltaBCSum<-rep(0,length(deltaBC))
    deltaBGSum<-rep(0,length(deltaBG))
    deltaGSum<-rep(0,length(deltaG))
    sumbetaBC<-sumbetaGam<-sumgam<-0
    
    betaGamTruePos<-gamTruePos<-betaBCTruePos<-0
    
    #  set H0, H1, ..., HK
    
    # D is the rings containing samples of orders. In total, K rings. 
    # the first element is for sample index, the second is for order of nodes (nodes), the order index (1), and the chain index (1), 
    # and the third is for ring index
    
    D<-array(NA,dim=c((K*orderloops),(nodes+2),K))
    Dindex<-rep(0,K)
    
    # order is for the ordering collected in each chain from each iteration. The first element is the chain index, the second is for 
    # the iteration index, and the third is the ordering (orderIndex + ordering). 
    order<-array(NA,dim=c(K, K*orderloops,(nodes+1)))
    
    # a beta matrix corresponding to each order with the first column an index for order
    # this index should corresponding to the order sampled. 
    OrderIndex<-1
    betaGamSampled<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    betaGamSampledUnique<-matrix(betaGam,nrow=1)
    deltaBGSampled<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    deltaBGSampledUnique<-matrix(deltaBG,nrow=1)
    
    gamSampled<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    gamSampledUnique<-matrix(gam,nrow=1)
    deltaGSampled<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    deltaGSampledUnique<-matrix(deltaG,nrow=1)
    
    betaBCSampled<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    betaBCSampledUnique<-matrix(betaBC,nrow=1)
    deltaBCSampled<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    deltaBCSampledUnique<-matrix(deltaBC,nrow=1)
    
    SigjInit<-matrix(c(1,0,0,1),byrow=TRUE,nrow=2)
    SigAllJUnique<-SigAllJUnique0<-matrix(rep(SigjInit,(J-1)),ncol=2,byrow=T)

    etaSampled<-array(NA,dim=c(K, K*orderloops,2))
    
    # The following value of H0 is based on initial runs and check the values of prob. 
    
    H0<-10000/8000
    
    geo<-2
    H1<-H0*geo
    H2<-H1*geo
    H3<-H2*geo
    H4<-H3*geo
    H<-c(H0,H1,H2,H3,H4)
    
    #  set T0 to TK
    T0<-1
    T1<-T0*geo
    T2<-T1*geo
    T3<-T2*geo
    T4<-T3*geo
    T<-c(T0,T1,T2,T3,T4)
    
    for (iii in 1:K)
    {
      order[iii,1,]<-c(OrderIndex,seq(1,nodes))

      betaGamSampled[iii,1,]<-c(OrderIndex,betaGam)
      deltaBGSampled[iii,1,]<-c(OrderIndex,deltaBG)
      
      gamSampled[iii,1,]<-c(OrderIndex,gam)
      deltaGSampled[iii,1,]<-c(OrderIndex,deltaG)
      
      betaBCSampled[iii,1,]<-c(OrderIndex,betaBC)
      deltaBCSampled[iii,1,]<-c(OrderIndex,deltaBC)
      
      etaSampled[iii,1,]<-c(OrderIndex,1)
    }

    probTmp<-0
    
    for (iii in 2:orderloops)
    {
      for (k in K:1)
      {
        if (iii>((K-k)*(B+N)))
        {
          if (k==K)
          {
            out<-MHStep(order,k,iii,nodes,etaSampled,betaGamSampled,betaGamSampledUnique,deltaBGSampled,deltaBGSampledUnique,
                        gamSampled,gamSampledUnique,deltaGSampled,deltaGSampledUnique,
                        betaBCSampled,betaBCSampledUnique,deltaBCSampled,deltaBCSampledUnique,
                        SigAllJUnique,nx,data,probTmp)		
            order<-out$order
            betaGamSampled<-out$betaGamSampled
            betaGamSampledUnique<-out$betaGamSampledUnique
            deltaBGSampled<-out$deltaBGSampled
            deltaBGSampledUnique<-out$deltaBGSampledUnique
            
            gamSampled<-out$gamSampled
            gamSampledUnique<-out$gamSampledUnique
            deltaGSampled<-out$deltaGSampled
            deltaGSampledUnique<-out$deltaGSampledUnique
            
            betaBCSampled<-out$betaBCSampled
            betaBCSampledUnique<-out$betaBCSampledUnique
            deltaBCSampled<-out$deltaBCSampled
            deltaBCSampledUnique<-out$deltaBCSampledUnique
            
            SigAllJUnique<-out$SigAllJUnique
            probTmp<-out$probTmp
            
            etaSampled<-out$etaSampled
          } # end if (k==K)
          else
          {
            judge<-runif(1)
            
            if (pee<judge)
            {
              out<-MHStep(order,k,iii,nodes,etaSampled,betaGamSampled,betaGamSampledUnique,deltaBGSampled,deltaBGSampledUnique,
                          gamSampled,gamSampledUnique,deltaGSampled,deltaGSampledUnique,
                          betaBCSampled,betaBCSampledUnique,deltaBCSampled,deltaBCSampledUnique,
                          SigAllJUnique,nx,data,probTmp)			
              
              order<-out$order
              betaGamSampled<-out$betaGamSampled
              betaGamSampledUnique<-out$betaGamSampledUnique
              deltaBGSampled<-out$deltaBGSampled
              deltaBGSampledUnique<-out$deltaBGSampledUnique
              
              gamSampled<-out$gamSampled
              gamSampledUnique<-out$gamSampledUnique
              deltaGSampled<-out$deltaGSampled
              deltaGSampledUnique<-out$deltaGSampledUnique
              
              betaBCSampled<-out$betaBCSampled
              betaBCSampledUnique<-out$betaBCSampledUnique
              deltaBCSampled<-out$deltaBCSampled
              deltaBCSampledUnique<-out$deltaBCSampledUnique
              
              SigAllJUnique<-out$SigAllJUnique
              probTmp<-out$probTmp
              
              etaSampled<-out$etaSampled				
            } # end if (pee<runif(1))
            else
            {
              if (iii==((K-k)*(B+N)+1))
              {
                orderTmp<-order[k,1,]
                deltaBG<-deltaBGSampledUnique[1,]
                betaGam<-betaGamSampledUnique[1,]
                
                deltaG<-deltaGSampledUnique[1,]
                gam<-gamSampledUnique[1,]
                
                deltaBC<-deltaBCSampledUnique[1,]
                betaBC<-betaBCSampledUnique[1,]
                
                SigAllJ<-SigAllJUnique[1:(2*(J-1)),]
              }
              else 
              {
                orderTmp<-order[k,(iii-1),]
                deltaBG<-deltaBGSampled[k,(iii-1),2:ncol(deltaBGSampled[k,,])]
                betaGam<-betaGamSampled[k,(iii-1),2:ncol(betaGamSampled[k,,])]
                
                deltaG<-deltaGSampled[k,(iii-1),2:ncol(deltaGSampled[k,,])]
                gam<-gamSampled[k,(iii-1),2:ncol(gamSampled[k,,])]
                
                deltaBC<-deltaBCSampled[k,(iii-1),2:ncol(deltaBCSampled[k,,])]
                betaBC<-betaBCSampled[k,(iii-1),2:ncol(betaBCSampled[k,,])]
                
                SigLoc<-order[k,(iii-1),1]
                SigAllJ<-SigAllJUnique[(2*(J-1)*(SigLoc-1)+1):(2*(J-1)*SigLoc),]
              }
              dataTmp<-data[,orderTmp[2:length(orderTmp)]]		
              
              probO<-OrderProb(deltaBG,betaGam,deltaG,gam,betaBC,n,SigAllJ,dataTmp)
              
              for (j in 1:(K-1))
              {
                if (H[j]<(-probO$prob) && (-probO$prob)<H[j+1])
                {
                  ring<-j
                  break
                }
              }
              if ((-probO$prob)>H[K])
              {
                ring<-K
              }
              
              # sample from rings 
              selectRing<-sample(seq(ring,K),1)
              len<-length(D[,2,selectRing])-sum(is.na(D[,2,selectRing]))
              if (len==0)
              {
                out<-MHStep(order,k,iii,nodes,etaSampled,betaGamSampled,betaGamSampledUnique,deltaBGSampled,deltaBGSampledUnique,
                            gamSampled,gamSampledUnique,deltaGSampled,deltaGSampledUnique,
                            betaBCSampled,betaBCSampledUnique,deltaBCSampled,deltaBCSampledUnique,
                            SigAllJUnique,nx,data,probTmp)	
                
                order<-out$order
                betaGamSampled<-out$betaGamSampled
                betaGamSampledUnique<-out$betaGamSampledUnique
                deltaBGSampled<-out$deltaBGSampled
                deltaBGSampledUnique<-out$deltaBGSampledUnique
                
                gamSampled<-out$gamSampled
                gamSampledUnique<-out$gamSampledUnique
                deltaGSampled<-out$deltaGSampled
                deltaGSampledUnique<-out$deltaGSampledUnique
                
                betaBCSampled<-out$betaBCSampled
                betaBCSampledUnique<-out$betaBCSampledUnique
                deltaBCSampled<-out$deltaBCSampled
                deltaBCSampledUnique<-out$deltaBCSampledUnique
                
                etaSampled<-out$etaSampled
                
                SigAllJUnique<-out$SigAllJUnique
              }
              else
              {
                
                loc1<-which(is.na(D[,2,selectRing]))
                cand<-D[-loc1,2,selectRing]
                if (length(cand)==1)
                {
                  loc<-which((D[,2,selectRing])==cand[1])
                }
                else
                {
                  sampleIndex<-sample(x=cand,size=1)
                  loc<-which((D[,2,selectRing])==sampleIndex)
                }
                
                orderIndex<-D[loc[1],2,selectRing]
                
                SigAllJTmp<-SigAllJUnique[(2*(J-1)*(orderIndex-1)+1):(2*(J-1)*orderIndex),]

                deltaBGTmp<-deltaBGSampledUnique[orderIndex,]
                betaGamTmp<-betaGamSampledUnique[orderIndex,]
                
                deltaGTmp<-deltaGSampledUnique[orderIndex,]
                gamTmp<-gamSampledUnique[orderIndex,]
                
                deltaBCTmp<-deltaBCSampledUnique[orderIndex,]
                betaBCTmp<-betaBCSampledUnique[orderIndex,]
                
                orderTmp<-D[loc[1],2:(nodes+2),selectRing]
                dataTmp<-data[,orderTmp[2:length(orderTmp)]]		
                
                probTmp<-OrderProb(deltaBGTmp,betaGamTmp,deltaGTmp,gamTmp,betaBCTmp,n,SigAllJTmp,dataTmp)
                
                lpyCurrent<--max(-probTmp$prob,H[(k)])/T[(k)]
                lpyHigh<--max(-probTmp$prob,H[(k+1)])/T[(k+1)]
                
                lpxCurrent<--max(-probO$prob,H[(k)])/T[(k)]
                
                lpxHigh<--max(-probO$prob,H[(k+1)])/T[(k+1)]
                lratio<-lpyCurrent+lpxHigh-lpxCurrent-lpyHigh
                if (exp(lratio)>runif(1))
                {
                  order[k,iii,]<-D[loc[1],2:(nodes+2),selectRing]
                  betaGamSampledUnique[order[k,iii,1],]<-probTmp$avgbetaGam
                  deltaBGSampledUnique[order[k,iii,1],]<-probTmp$avgdeltaBG
                  
                  gamSampledUnique[order[k,iii,1],]<-probTmp$avggam
                  deltaGSampledUnique[order[k,iii,1],]<-probTmp$avgdeltaG
                  
                  betaBCSampledUnique[order[k,iii,1],]<-probTmp$avgbetaBC
                  deltaBCSampledUnique[order[k,iii,1],]<-probTmp$avgdeltaBC
                  
                  SigLoc<-order[k,iii,1]
                  
                  SigAllJUnique[(2*(J-1)*(SigLoc-1)+1):(2*(J-1)*SigLoc),]<-probTmp$SigAllJ
                  
                  betaGamSampled[k,iii,2:ncol(betaGamSampled[k,,])]<-probTmp$avgbetaGam
                  betaGamSampled[k,iii,1]<-order[k,iii,1]
                  deltaBGSampled[k,iii,2:ncol(deltaBGSampled[k,,])]<-probTmp$avgdeltaBG
                  deltaBGSampled[k,iii,1]<-order[k,iii,1]
                  
                  gamSampled[k,iii,2:ncol(gamSampled[k,,])]<-probTmp$avggam
                  gamSampled[k,iii,1]<-order[k,iii,1]
                  deltaGSampled[k,iii,2:ncol(deltaGSampled[k,,])]<-probTmp$avgdeltaG
                  deltaGSampled[k,iii,1]<-order[k,iii,1]
                  
                  betaBCSampled[k,iii,2:ncol(betaBCSampled[k,,])]<-probTmp$avgbetaBC
                  betaBCSampled[k,iii,1]<-order[k,iii,1]
                  deltaBCSampled[k,iii,2:ncol(deltaBCSampled[k,,])]<-probTmp$avgdeltaBC
                  deltaBCSampled[k,iii,1]<-order[k,iii,1]
                  
                  etaSampled[k,iii,2]<-probTmp$eta
                  etaSampled[k,iii,1]<-order[k,iii,1]
                }
                else
                {
                  if (iii==((K-k)*(B+N)+1))
                  {
                    order[k,iii,]<-order[k,1,]
                  }
                  else
                  {
                    order[k,iii,]<-order[k,iii-1,]
                  }
                  betaGamSampledUnique[order[k,iii,1],]<-probO$avgbetaGam
                  deltaBGSampledUnique[order[k,iii,1],]<-probO$avgdeltaBG
                  
                  gamSampledUnique[order[k,iii,1],]<-probO$avggam
                  deltaGSampledUnique[order[k,iii,1],]<-probO$avgdeltaG
                  
                  betaBCSampledUnique[order[k,iii,1],]<-probO$avgbetaBC
                  deltaBCSampledUnique[order[k,iii,1],]<-probO$avgdeltaBC
                  
                  SigLoc<-order[k,iii,1]
                  
                  SigAllJUnique[(2*(J-1)*(SigLoc-1)+1):(2*(J-1)*SigLoc),]<-probO$SigAllJ
                  
                  betaGamSampled[k,iii,2:ncol(betaGamSampled[k,,])]<-probO$avgbetaGam
                  betaGamSampled[k,iii,1]<-order[k,iii,1]
                  deltaBGSampled[k,iii,2:ncol(deltaBGSampled[k,,])]<-probO$avgdeltaBG
                  deltaBGSampled[k,iii,1]<-order[k,iii,1]
                  
                  gamSampled[k,iii,2:ncol(gamSampled[k,,])]<-probO$avggam
                  gamSampled[k,iii,1]<-order[k,iii,1]
                  deltaGSampled[k,iii,2:ncol(deltaGSampled[k,,])]<-probO$avgdeltaG
                  deltaGSampled[k,iii,1]<-order[k,iii,1]
                  
                  betaBCSampled[k,iii,2:ncol(betaBCSampled[k,,])]<-probO$avgbetaBC
                  betaBCSampled[k,iii,1]<-order[k,iii,1]
                  deltaBCSampled[k,iii,2:ncol(deltaBCSampled[k,,])]<-probO$avgdeltaBC
                  deltaBCSampled[k,iii,1]<-order[k,iii,1]
                  
                  etaSampled[k,iii,2]<-probO$eta
                  etaSampled[k,iii,1]<-order[k,iii,1]
                } # end else (exp(lratio)>runif(1))
              } # end else	(len==0)				
            } # end else (pee<runif(1))
            #		cat("k ",k," iii ",iii," etaSampled ",etaSampled[k,iii,2],"\n")
          } # end else corresponding to if (k==K)
        } # end if (iii>((K-k)*(B+N))
        if (iii>((K-k)*(B+N)+B))
        {
          orderIndex<-order[k,iii,1]
          SigAllJ<-SigAllJUnique[(2*(J-1)*(orderIndex-1)+1):(2*(J-1)*orderIndex),]
          
          deltaBG<-deltaBGSampled[k,iii,2:ncol(deltaBGSampled[k,,])]
          betaGam<-betaGamSampled[k,iii,2:ncol(betaGamSampled[k,,])]
          
          deltaG<-deltaGSampled[k,iii,2:ncol(deltaGSampled[k,,])]
          gam<-gamSampled[k,iii,2:ncol(gamSampled[k,,])]
          
          deltaBC<-deltaBCSampled[k,iii,2:ncol(deltaBCSampled[k,,])]
          betaBC<-betaBCSampled[k,iii,2:ncol(betaBCSampled[k,,])]
          
          orderTmp<-order[k,iii,]
          dataTmp<-data[,orderTmp[2:length(orderTmp)]]		
          
          probO<-OrderProb(deltaBG,betaGam,deltaG,gam,betaBC,n,SigAllJ,dataTmp)
          
          for (j in 1:(K-1))
          {
            if (H[j]<(-probO$prob)&& (-probO$prob)<H[(j+1)])
            {
              D[Dindex[j],1,j]<-k
              D[Dindex[j],2,j]<-order[k,iii,1]
              D[Dindex[j],3:(nodes+2),j]<-order[k,iii,2:(nodes+1)]
              Dindex[j]<-Dindex[j]+1
              break
            }
          }
          if ((-probO$prob)>H[K])
          {
            D[Dindex[K],1,K]<-k
            D[Dindex[K],2,K]<-order[k,iii,1]
            D[Dindex[K],3:(nodes+2),K]<-order[k,iii,2:(nodes+1)]
            Dindex[K]<-Dindex[K]+1
          }
        } # end if (iii>((K-k)*(B+N)+B)
      } # end for (k in K:1)
    } # end for (iii in 2:orderloops)

    parentsTruedeltaBG<-parentsTruedeltaG<-parentsTruedeltaBC<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)
    parentsEstdeltaBG<-parentsEstdeltaG<-parentsEstdeltaBC<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)
    parentsTruebetaGam<-parentsTruegam<-parentsTruebetaBC<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)
    parentsEstbetaGam<-parentsEstgam<-parentsEstbetaBC<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)
    
    orderEst<-rep(0,nodes)
    
    for (kk in 0:(nodes-2))
    {
      start<-sum(0:kk)+1
      end<-start+kk
      t<-1
      for (kkk in start:end)
      {
        parentsTruedeltaBG[(kk+2),t]<-deltaBGTrue[1,kkk]
        parentsTruedeltaG[(kk+2),t]<-deltaGTrue[1,kkk]
        parentsTruedeltaBC[(kk+2),t]<-deltaBCTrue[1,kkk]
        
        parentsTruebetaGam[(kk+2),t]<-truebetaGam[1,kkk]
        parentsTruegam[(kk+2),t]<-truegam[1,kkk]
        parentsTruebetaBC[(kk+2),t]<-truebetaBC[1,kkk]
        
        t<-t+1
      }
    }
    
    Sumeta<-sumEtaEst<-0
    counts<-0
    for (iii in ((K-1)*(B+N)+B):((K-1)*(B+N)+B+B+2*N))
    {
      for (t in 1:K)
      {
        if (counts==0)
        {
          sumEtaEst<-etaSampled[t,iii,2]
          counts<-counts+1
        }
        else 
        {
          sumEtaEst<-etaSampled[t,iii,2]+sumEtaEst
          counts<-counts+1
        }
      }
    }
    avgEtaEst[jj]<-sumEtaEst/counts
    cat("avgEtaEst jj ",jj," ",avgEtaEst[jj],"\n")
    
    counts<-0
    
    for (iii in ((K-1)*(B+N)+B):((K-1)*(B+N)+B+B+2*N))
    {
      for (t in 1:K)
      {
        for (kk in 0:(nodes-2))
        {
          orderEst[1:(kk+2)]<-order[t,iii, 2:(kk+3)]
          start<-sum(0:kk)+1+1
          end<-start+kk
          tt<-1
          for (kkk in start:end)
          {
            if (avgEtaEst[jj]>0.5)
            {
              if (etaSampled[t,iii,2]==1)
              {
                parentsEstdeltaBG[orderEst[kk+2],orderEst[tt]]<-deltaBGSampled[t,iii,kkk]
                parentsEstdeltaG[orderEst[kk+2],orderEst[tt]]<-deltaGSampled[t,iii,kkk]
                parentsEstbetaGam[orderEst[kk+2],orderEst[tt]]<-betaGamSampled[t,iii,kkk]
                parentsEstgam[orderEst[kk+2],orderEst[tt]]<-gamSampled[t,iii,kkk]
                tt<-tt+1
              }
            }
            else
            {
              if (etaSampled[t,iii,2]==0)
              {
                parentsEstdeltaBC[orderEst[kk+2],orderEst[tt]]<-deltaBCSampled[t,iii,kkk]
                parentsEstbetaBC[orderEst[kk+2],orderEst[tt]]<-betaBCSampled[t,iii,kkk]
                tt<-tt+1
              }
            }	
          }
        }
        # avgEtaEst[jj]>0.5 means differential; <=0.5 identical. 
        if (counts==0)
        {
          if (avgEtaEst[jj]>0.5 && etaSampled[t,iii,2]==1)
          {
            sumParentsEstdeltaBG<-parentsEstdeltaBG
            sumParentsEstdeltaG<-parentsEstdeltaG
            sumParentsEstbetaGam<-parentsEstbetaGam
            sumParentsEstgam<-parentsEstgam
            counts<-counts+1
          }
          if (avgEtaEst[jj]<=0.5 && etaSampled[t,iii,2]==0)
          {
            sumParentsEstdeltaBC<-parentsEstdeltaBC
            sumParentsEstbetaBC<-parentsEstbetaBC
            counts<-counts+1
          }
        }
        else 
        {
          if (avgEtaEst[jj]>0.5 && etaSampled[t,iii,2]==1)
          {			
            sumParentsEstdeltaBG<-parentsEstdeltaBG+sumParentsEstdeltaBG
            sumParentsEstdeltaG<-parentsEstdeltaG+sumParentsEstdeltaG
            sumParentsEstbetaGam<-parentsEstbetaGam+sumParentsEstbetaGam
            sumParentsEstgam<-parentsEstgam+sumParentsEstgam
            counts<-counts+1
          }
          if (avgEtaEst[jj]<=0.5 && etaSampled[t,iii,2]==0)
          {
            sumParentsEstdeltaBC<-parentsEstdeltaBC+sumParentsEstdeltaBC
            sumParentsEstbetaBC<-parentsEstbetaBC+sumParentsEstbetaBC
            counts<-counts+1
          }
        }
      }
    }
    
    if (avgEtaEst[jj]>0.5)
    {
      avgParentsEstdeltaBG<-round(sumParentsEstdeltaBG/counts,0)
      diffdeltaBG<-avgParentsEstdeltaBG-parentsTruedeltaBG
      accurdeltaBG[jj]<-1-sum(abs(diffdeltaBG))/(nodes*(nodes-1)/2)
      fpdeltaBG[jj]<-sum(diffdeltaBG>0)/(nodes*(nodes-1)/2-sum(parentsTruedeltaBG))
      tpdeltaBG[jj]<-1-abs(sum(diffdeltaBG<0))/sum(parentsTruedeltaBG)
      
      avgParentsEstdeltaG<-round(sumParentsEstdeltaG/counts,0)
      diffdeltaG<-avgParentsEstdeltaG-parentsTruedeltaG
      accurdeltaG[jj]<-1-sum(abs(diffdeltaG))/(nodes*(nodes-1)/2)
      fpdeltaG[jj]<-sum(diffdeltaG>0)/(nodes*(nodes-1)/2-sum(parentsTruedeltaG))
      tpdeltaG[jj]<-1-abs(sum(diffdeltaG<0))/sum(parentsTruedeltaG)
    }
    else
    {
      avgParentsEstdeltaBC<-round(sumParentsEstdeltaBC/counts,0)
      diffdeltaBC<-avgParentsEstdeltaBC-parentsTruedeltaBC
      accurdeltaBC[jj]<-1-sum(abs(diffdeltaBC))/(nodes*(nodes-1)/2)
      fpdeltaBC[jj]<-sum(diffdeltaBC>0)/(nodes*(nodes-1)/2-sum(parentsTruedeltaBC))
      tpdeltaBC[jj]<-1-abs(sum(diffdeltaBC<0))/sum(parentsTruedeltaBC)
    }
  } # end data=jj

end_time <- Sys.time()
end_time - start_time

save(avgEtaEst,accurdeltaBC,fpdeltaBC,tpdeltaBC,accurdeltaBG,fpdeltaBG,tpdeltaBG,accurdeltaG,fpdeltaG,tpdeltaG,file=paste(basefn,".Rdata",sep=""))

