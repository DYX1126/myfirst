library(dplyr)
library(rootSolve)
library(nleqslv)
library(statmod)
library(stats4)
re=10
s40=1/(273.15+40)
s65=1/(273.15+65)

s100=1/(273.15+100)
w65=(s65-s40)/(s100-s40)

w100=1

G.a2 <- -2.3532

G.b2 <- 1.9201

G.alpha <-  1.3696

G.d <- 0.4980


G.u40 <- exp(G.a2)
G.u65 <- exp(G.a2+G.b2*w65)
G.u100 <- exp(G.a2+G.b2)

w=30
q<-0.01
GAqTRUE <- uniroot(function(t)  pgamma(w, shape=G.alpha*(t^G.d),scale =G.u80/G.alpha,lower = FALSE)-q, lower = 0, upper =25*10^20, tol = 0.001)$root

##############--------------------
n=6
n65=n
n100=n
mis_ratio<-matrix(0,30,30)
##############-----------------------------------------

data=read.csv("C:/Users/user/Desktop/misspecification/yang.csv",header = T)
da=as.matrix(data)


#65«×
t65=da[1:12,2]
y65=da[1:12,3:8]

#85«×
t85=da[13:23,2]
y85=da[13:23,3:8]
#100«×
t100=da[24:34,2]
y100=da[24:34,3:8]

###
z65=matrix(0,11,6)
for(j in 1:6){
  for(i in 1:11){
    z65[i,j]=y65[i+1,j]-y65[i,j]
  }
}

z85=matrix(0,10,6)
for(j in 1:6){
  for(i in 1:10){
    z85[i,j]=y85[i+1,j]-y85[i,j]
  }
  
}

z100=matrix(0,10,6)
for(j in 1:6){
  for(i in 1:10){
    z100[i,j]=y100[i+1,j]-y100[i,j]
  }
}


n65=ncol(y65)
m65=nrow(y65)-1

n85=ncol(y85)
m85=nrow(y85)-1

n100=ncol(y100)
m100=nrow(y100)-1

#####-----------------------------------------------
library(statmod)
    G.d65 <- (t65[-1]^G.d)-(t65[-length(t65)]^G.d)
    G.d100 <- (t100[-1]^G.d)-(t100[-length(t100)]^G.d)
    seed<-NA
    dis<-NA
    for(iii in 1:1000){
      set.seed(iii)
      
      z65.simRE=matrix(0,(length(t65)-1),n)
      for(i in 1:(length(t65)-1)){
        z65.simRE[i,]= round( rgamma(n,shape = G.alpha *G.d65[i],scale = G.u65/G.alpha),re)
      }
      
      z100.simRE=matrix(0,(length(t100)-1),n)
      for(i in 1:(length(t100)-1)){
        z100.simRE[i,]= round( rgamma(n,shape = G.alpha *G.d100[i],scale = G.u100/G.alpha),re)
      }
      
      
      ###############################---------------------------------------------------
      MLEADT =function(a,lb,llam,ldt) {
        b=exp(lb)
        mu65=exp(a+b*w65)
        mu100=exp(a+b*w100)
        lam=exp(llam)
        dt=exp(ldt)
        k65=sapply(1:(length(t65)-1),
                    function(i){
                      (t65[i+1]^dt)-(t65[i]^dt)})
        k100=sapply(1:(length(t100)-1),
                    function(i){
                      (t100[i+1]^dt)-(t100[i]^dt)})
        return(-(
          0.5*n65*(m65*log(lam)+sum(log(k65^2)))+
            sum((-lam*(z65.simRE-mu65*(k65))^2)/(2*(mu65^2)*z65.simRE))-
            0.5*sum(log(2*pi*z65.simRE^3))+
            0.5*n100*(m100*log(lam)+sum(log(k100^2)))+
            sum((-lam*(z100.simRE-mu100*(k100))^2)/(2*(mu100^2)*z100.simRE))-
            0.5*sum(log(2*pi*z100.simRE^3))
        ))
      }
      
      parsADT=mle(MLEADT, start = list(a=G.a2,
                                       lb=log(G.b2),
                                       llam=log(0.36),
                                       ldt=log(G.d)),control=list(maxit=1000000)) 
      
      #MLEADT(IG.a2,log(IG.b2),log(IG.lam),log(IG.d))
      AIC_IG_ADT<-8+2*parsADT@details$value
      BIC_IG_ADT<-log(n)*4+2*parsADT@details$value
      
      igADT=exp(parsADT@coef)
      igADT[1]=parsADT@coef[1]
      #     (igADT=round(igADT,digits = 6))
      ########-------------------------
      library(stats4)
      qmleADT =function(aa,lb,la,ld) {
        b=exp(lb)
        lmu65=(aa+b*w65)
        mu65=exp(aa+b*w65)
        lmu100=(aa+b*w100)
        mu100=exp(aa+b*w100)
        a=exp(la) 
        d=exp(ld)
        k65=sapply(1:(length(t65)-1),
                    function(i){
                      (t65[i+1]^d)-(t65[i]^d)})
        k100=sapply(1:(length(t100)-1),
                    function(i){
                      (t100[i+1]^d)-(t100[i]^d)})
        return(-( 
          (n65*a*(-(lmu65-la))*(sum(k65)))-
            n65*sum(lgamma(a*k65))+a*sum(k65*log(z65.simRE))-
            exp(-(lmu65-la))*sum(z65.simRE)-sum(log(z65.simRE))+
            (n100*a*(-(lmu100-la))*(sum(k100)))-
            n100*sum(lgamma(a*k100))+a*sum(k100*log(z100.simRE))-
            exp(-(lmu100-la))*sum(z100.simRE)-sum(log(z100.simRE)))
        )
      }
      GAparsADT=mle(qmleADT, start = list(aa=G.a2,
                                          lb=log(G.b2),
                                          la=log(G.alpha),
                                          ld=log(G.d)),
                    control=list(maxit=10000)) 
      # str(GAparsADT)
      AIC_GA_ADT<- 8+2*GAparsADT@details$value
      BIC_GA_ADT<-log(n)*4+2*GAparsADT@details$value
      
      gammaADT=exp(GAparsADT@coef)
      gammaADT[1]=GAparsADT@coef[1]
      #     (gammaADT=round(gammaADT,digits = 6))
      seed[iii] <- AIC_IG_ADT<AIC_GA_ADT
      dis[iii]<- AIC_IG_ADT-AIC_GA_ADT
      print(iii)
    }
    
  mis_ratio<-length(which(seed==TRUE))/1000 
  #0.015
  
  
#     write.csv(seed,file=paste0("C:/Users/user/Desktop/misspecified/MisspecificationWithAIC/GAMMA_DATA/n_",n65,"/mL_",m65,"_mH_",m100,"_seed.csv"))
#     write.csv(dis,file=paste0("C:/Users/user/Desktop/misspecified/MisspecificationWithAIC/GAMMA_DATA/n_",n65,"/mL_",m65,"_mH_",m100,"_dis.csv"))
#     print(paste0("m100=",m100,"; m65=",m65))    
# 
# write.csv(mis_ratio,file=paste0("C:/Users/user/Desktop/misspecified/MisspecificationWithAIC/GAMMA_DATA/n_",n65,"/mis_ratio.csv"))


