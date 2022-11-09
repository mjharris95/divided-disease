library(deSolve)
library(tidyverse)
library(reshape2)

full_aware<-function(h=NA, epsilon=NA, mu=NA, kappa=NA, beta=NA, rho=NA, 
                     theta=5, phi=.1, kappa_a=.1,  kappa_b=.1, rho_a=.1,  
                     mu_a=.05, mu_b=.05, beta_a=.5, beta_b=.5, rho_b=.1,
                     ha=.99, hb=.99, epsilon_a=.5, epsilon_b=.5, ell=1, I0_a=.005/2, I0_b=.005/2, v_val=0, v_start=0, time=500, theta_a=.2, theta_b=.2,
                     get_params=FALSE){
  if(!is.na(beta)){
    beta_a<-beta
    beta_b<-beta
  }
  
  if(!is.na(kappa)){
    kappa_a<-kappa
    kappa_b<-kappa
  }
  
  if(!is.na(mu)){
    mu_a<-mu
    mu_b<-mu
  }
  
  if(!is.na(rho)){
    rho_a<-rho
    rho_b<-rho
  }
  
  if(!is.na(h)){
    ha<-h
    hb<-h
  }
  
  if(!is.na(epsilon)){
    epsilon_a<-epsilon
    epsilon_b<-epsilon
  }
  
  if(theta != 5){
    theta_a<-theta
    theta_b<-theta
  }
  
  params<-c("kappa_a"=kappa_a, "kappa_b"=kappa_b,  "theta"=theta, 
            "beta_a"=beta_a, "beta_b"=beta_b, "rho"=rho, "phi"=phi, 
            "rho_a"=rho_a, "rho_b"=rho_b,
            "mu_a"=mu_a, "mu_b"=mu_b, "epsilon_a"=epsilon_a, "epsilon_b"=epsilon_b,  "ha"=ha, "hb"=hb, "ell"=ell, "v_start"=v_start, "v_val"=v_val, "time"=time, "I0_a"=I0_a, "I0_b"=I0_b, "theta_a"=theta_a, "theta_b"=theta_b) 
  state<-c( 
    
    #group a
    SUa<-1-I0_a,
    SPa<-0,
    IUa<-I0_a,
    IPa<-0,
    RUa<-0,
    RPa<-0,
    DUa<-0,
    DPa<-0,
    
    
    #group b
    SUb<-1-I0_b,
    SPb<-0,
    IUb<-I0_b,
    IPb<-0,
    RUb<-0,
    RPb<-0,
    DUb<-0,
    DPb<-0

   
  )
  
  names(state)<- c("SUa", "SPa", "IUa", "IPa", "RUa", "RPa", "DUa", "DPa", 
                   "SUb", "SPb", "IUb", "IPb", "RUb", "RPb", "DUb", "DPb")
  
  sir_upu<-function(t, state, parameter){
    
    with(as.list(c(parameter, state)), {
      
      
      if(t<ell){
        lag<-rep(0, 16)
      }
      else{
        lag<-lagvalue(t-ell)
      }
      
      #transitions for group a
      if(t>v_start){
        v <- v_val
      }
      else{
        v<-0
      }
      
      dSUa <- -sqrt(beta_a)*SUa*(sqrt(beta_a)*ha*(IUa+kappa_a*IPa)+sqrt(beta_b)*(1-ha)*(IUb+kappa_b*IPb)) - 
        theta_a * SUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) + phi * SPa  
      
      dSPa <- -sqrt(beta_a)*kappa_a*SPa*(sqrt(beta_a)*ha*(IUa+kappa_a*IPa)+sqrt(beta_b)*(1-ha)*(IUb+kappa_b*IPb)) + 
        theta_a * SUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) - phi * SPa  -
        v*SPa
      
      dIUa <- sqrt(beta_a)*SUa*(sqrt(beta_a)*ha*(IUa+kappa_a*IPa)+sqrt(beta_b)*(1-ha)*(IUb+kappa_b*IPb)) - 
        theta_a * IUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) + phi * IPa -
        (rho_a) * IUa
      
      dIPa <- sqrt(beta_a)*kappa_a*SPa*(sqrt(beta_a)*ha*(IUa+kappa_a*IPa)+sqrt(beta_b)*(1-ha)*(IUb+kappa_b*IPb)) + 
        theta_a * IUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16])  ) - phi * IPa -
        (rho_a) * IPa
      
      dRUa <- rho_a * (1-mu_a) * IUa - 
        theta_a * RUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) + phi * RPa  
      
      dRPa <- rho_a * (1-mu_a) * IPa + v * SPa +
        theta_a * RUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16])  ) - phi * RPa  
      
      dDUa <- rho_a*mu_a*IUa
      
      dDPa <- rho_a*mu_a*IPa
  
      #transitions for group b
  
      dSUb <- -sqrt(beta_b)*SUb*(sqrt(beta_a)*(1-hb)*(IUa+kappa_a*IPa)+sqrt(beta_b)*(hb)*(IUb+kappa_b*IPb)) - 
        theta_b * SUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16]) ) + phi * SPb
      
      dSPb <- -sqrt(beta_b)*kappa_b*SPb*(sqrt(beta_a)*(1-hb)*(IUa+kappa_a*IPa)+sqrt(beta_b)*(hb)*(IUb+kappa_b*IPb)) + 
        theta_b * SUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16])  ) - phi * SPb -
        v * SPb
      
      dIUb <- sqrt(beta_b)*SUb*(sqrt(beta_a)*(1-hb)*(IUa+kappa_a*IPa)+sqrt(beta_b)*(hb)*(IUb+kappa_b*IPb)) - 
        theta_b * IUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16])  ) + phi * IPb -
       (rho_b) * IUb
      
      dIPb <- sqrt(beta_b)*kappa_b*SPb*(sqrt(beta_a)*(1-hb)*(IUa+kappa_a*IPa)+sqrt(beta_b)*(hb)*(IUb+kappa_b*IPb)) +
        theta_b * IUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16]) ) - phi * IPb -
        (rho_b) * IPb
      
      dRUb <- (1-mu_b) * rho_b * IUb - 
        theta_b * RUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16]) ) + phi * RPb
      
      dRPb <- (1-mu_b) * rho_b * IPb + v * SPb +
        theta_b * RUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16])  ) - phi * RPb
      
      dDUb <- rho_b*mu_b*IUb
      
      dDPb <- rho_b*mu_b*IPb
      

      
      return(list(c(dSUa, dSPa, dIUa, dIPa, dRUa, dRPa, dDUa, dDPa,
                    dSUb, dSPb, dIUb, dIPb, dRUb, dRPb, dDUb, dDPb
      )))
    })
  }
  
  times<-seq(from=0, to=time, by=1)

  as.data.frame(dede(state, times, sir_upu, params))->sim
  
  if(get_params){
    return(list(sim=sim, params=params))
  }
  
  else{
    return(sim)
  }

}

# ls<-full_aware(h=.99, epsilon=.5, mu=.01, kappa=.3, beta=.2, rho=.1, theta=100, ell=1, phi=0, I0_a=.001, I0_b=0, get_params=TRUE)
# sim<-ls$sim
# params<-ls$params
#
# 
# Smax<-max(c(sim$SUa, sim$SPa, sim$SUb, sim$SPb))
# Imax<-max(c(sim$IUa, sim$IPa, sim$IUb, sim$IPb))
# Rmax<-max(c(sim$RUa, sim$RPa, sim$RUb, sim$RPb))
# 
# par(mfcol=c(2, 5), lwd=2)
# 
# 
# plot(SUa~time, data=sim, type="l", col="#3D85BD", ylim=c(0, Smax), ylab="Susceptible", main="a")
# lines(SPa~time, data=sim, col="#3D85BD", lty=2)
# plot(SUb~time, data=sim, col="#3D85BD", type="l", ylim=c(0, Smax), main="b", ylab="Susceptible")
# lines(SPb~time, data=sim, col="#3D85BD", lty=2)
# 
# plot(IUa~time, data=sim, type="l", col="#7CAA2D", ylim=c(0, Imax), ylab="Infected")
# lines(IPa~time, data=sim, col="#7CAA2D", lty=2)
# plot(IUb~time, data=sim, col="#7CAA2D", type="l", ylim=c(0, Imax), ylab="Infected")
# lines(IPb~time, data=sim, col="#7CAA2D", lty=2)
# 
# plot(RUa~time, data=sim, type="l", col="#F5E356", ylim=c(0, Rmax), ylab="Recovered")
# lines(RPa~time, data=sim, col="#F5E356", lty=2)
# plot(RUb~time, data=sim, col="#F5E356", type="l", ylim=c(0, Rmax), ylab="Recovered")
# lines(RPb~time, data=sim, col="#F5E356", lty=2)
# 
# plot(SUa+IUa+RUa~time, data=sim, type="l", col="purple", ylim=c(0, Amax), ylab="Attitude")
# plot(SUb+IUb+RUb~time, data=sim, col="purple", type="l", ylim=c(0, Amax), ylab="Attitude")
# 
# sim %>% melt(id="time") %>% filter(variable %in% c("DUa", "DPa")) %>% group_by(time) %>% summarise(Da=sum(value)) -> df
# sim %>% melt(id="time") %>% filter(variable %in% c("DUb", "DPb")) %>% group_by(time) %>% summarise(Db=sum(value)) %>% right_join(df) -> df
# 
# lapply(unique(df$time), function(x) c(ifelse(x-params[["ell"]]<0, 0, x-params[["ell"]]),x)) %>%
#   lapply(function(times) filter(df, time %in% times) %>% summarise(alpha_a=params[["epsilon_a"]]*diff(Da)+(1-params[["epsilon_a"]])*diff(Db),
#                                                                    alpha_b=(1-params[["epsilon_b"]])*diff(Da)+(params[["epsilon_b"]])*diff(Db),
#                                                                    time=tail(times, 1))) %>%
#   do.call(rbind, .) -> awareness
# 
# a_lim<-range(awareness$alpha_a, awareness$alpha_b)
# plot(alpha_a~time, data=awareness, type="l", col="#CB6318", ylim=a_lim, ylab="Awareness")
# plot(alpha_b~time, data=awareness, col="#CB6318", type="l", ylim=a_lim, ylab="Awareness")

