library(deSolve)
library(tidyverse)
library(reshape2)


full_aware_vax<-function(h=NA, epsilon=NA, mu=NA, kappa=.1, beta=NA, rho=.1, 
                     theta=NA, phi=.01,
                     mu_a=.05, mu_b=.05, beta_a=.5, beta_b=.5,
                     ha=.99, hb=.99, epsilon_a=.5, epsilon_b=.5, ell=1, 
                     I0_a=.005/2, I0_b=.005/2, v_start=0, time=500, 
                     theta_a=50, theta_b=50, zeta=0.1, omega=.01,
                     get_params=FALSE){
  if(!is.na(beta)){
    beta_a<-beta
    beta_b<-beta
  }
  
  if(!is.na(mu)){
    mu_a<-mu
    mu_b<-mu
  }
  
  if(!is.na(h)){
    ha<-h
    hb<-h
  }
  
  if(!is.na(epsilon)){
    epsilon_a<-epsilon
    epsilon_b<-epsilon
  }
  
  if(!is.na(theta)){
    theta_a<-theta
    theta_b<-theta
  }
  
  params<-c("kappa"=kappa, "omega"=omega,
            "beta_a"=beta_a, "beta_b"=beta_b, "rho"=rho, "phi"=phi, 
            "mu_a"=mu_a, "mu_b"=mu_b, "zeta"=zeta, "epsilon_a"=epsilon_a, 
            "epsilon_b"=epsilon_b,  "ha"=ha, "hb"=hb, "ell"=ell, 
            "v_start"=v_start, "time"=time, "I0_a"=I0_a, "I0_b"=I0_b, 
            "theta_a"=theta_a, "theta_b"=theta_b) 
  state<-c( 
    
    #group a
    SUa<-1-I0_a,
    STa<-0,
    SMa<-0,
    IUa<-I0_a,
    ITa<-0,
    IMa<-0,
    RTa<-0,
    DUa<-0,
    DTa<-0,
    DMa<-0,
    
    
    #group b
    SUb<-1-I0_b,
    STb<-0,
    SMb<-0,
    IUb<-I0_b,
    ITb<-0,
    IMb<-0,
    RTb<-0,
    DUb<-0,
    DTb<-0,
    DMb<-0,
    
    COa<-0,
    COb<-0,
    OVa<-0,
    OVb<-0, 
    OBa<-0,
    OBb<-0
    
  )
  
  names(state)<- c("SUa", "STa", "SMa", "IUa", "ITa", "IMa", "RTa", "DUa", "DTa", "DMa",
                   "SUb", "STb", "SMb", "IUb", "ITb", "IMb", "RTb", "DUb", "DTb", "DMb",
                   "COa", "COb", "OVa", "OVb", "OBa", "OBb")
  
  sir_upu<-function(t, state, parameter){
    
    with(as.list(c(parameter, state)), {
      
      
      if(t<ell){
        lag<-rep(0, 20)
      }
      else{
        lag<-lagvalue(t-ell)
      }
      
      
      if(t<v_start){
        theta_a_v <- 0
        theta_b_v <- 0
      }
      else{
        theta_a_v <-theta_a
        theta_b_v <- theta_b
      }
      
      dSUa <- -sqrt(beta_a)*SUa*(sqrt(beta_a)*ha*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(1-ha)*(IUb+IMb+kappa*ITb)) - 
        theta_a_v * SUa * (epsilon_a * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (1-epsilon_a) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) 
      
      dSTa <- -sqrt(beta_a)*kappa*STa*(sqrt(beta_a)*ha*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(1-ha)*(IUb+IMb+kappa*ITb)) - 
        theta_a_v * STa * (epsilon_a * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (1-epsilon_a) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) -
        phi * STa + (omega)*(RTa) 
        
      
      dSMa <- -sqrt(beta_a)*SMa*(sqrt(beta_a)*ha*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(1-ha)*(IUb+IMb+kappa*ITb)) - 
        theta_a_v * SMa * (epsilon_a * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (1-epsilon_a) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) +
        phi * STa 
      
      
      dIUa <- sqrt(beta_a)*SUa*(sqrt(beta_a)*ha*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(1-ha)*(IUb+IMb+kappa*ITb)) - 
        (rho) * IUa 
      
      dITa <- sqrt(beta_a)*kappa*STa*(sqrt(beta_a)*ha*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(1-ha)*(IUb+IMb+kappa*ITb)) - 
        (rho) * ITa 
      
      dIMa <- sqrt(beta_a)*SMa*(sqrt(beta_a)*ha*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(1-ha)*(IUb+IMb+kappa*ITb)) -
        (rho) * IMa 
      
      
      dRTa <- rho*(1-mu_a)*(IUa) + rho*(1-zeta*mu_a)* (ITa+IMa) - omega * (RTa) +
        theta_a_v * (SMa+SUa+STa) * (epsilon_a * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (1-epsilon_a) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) )
      
      dDUa <- mu_a * rho * IUa
      
      dDTa <- mu_a * zeta * rho * ITa
      
      dDMa <- mu_a * zeta* rho * IMa
      
      #transitions for group b
      
      dSUb <- -sqrt(beta_b)*SUb*(sqrt(beta_a)*(1-hb)*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*hb*(IUb+IMb+kappa*ITb)) - 
        theta_b_v * SUb * ((1-epsilon_b) * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (epsilon_b) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) 
      
      dSTb <- -sqrt(beta_b)*kappa*STb*(sqrt(beta_a)*(1-hb)*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(hb)*(IUb+IMb+kappa*ITb))  -
        theta_b_v * STb * ((1-epsilon_b) * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (epsilon_b) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) -
        phi * STb + (omega)*(RTb)
      
      dSMb <- -sqrt(beta_b)*SMb*(sqrt(beta_a)*(1-hb)*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(hb)*(IUb+IMb+kappa*ITb)) - 
        theta_b_v * SMb * ((1-epsilon_b) * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (epsilon_b) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) +
        phi * STb 
      
      
      dIUb <- sqrt(beta_b)*SUb*(sqrt(beta_a)*(1-hb)*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(hb)*(IUb+IMb+kappa*ITb)) - 
        (rho) * IUb 
      
      dITb <- sqrt(beta_b)*kappa*STb*(sqrt(beta_a)*(1-hb)*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(hb)*(IUb+IMb+kappa*ITb)) - 
        (rho) * ITb
      
      dIMb <- sqrt(beta_b)*SMb*(sqrt(beta_a)*(1-hb)*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(hb)*(IUb+IMb+kappa*ITb)) -
        (rho) * IMb 
      
      
      dRTb <- rho*(1-mu_b)*(IUb) + rho*(1-zeta*mu_b)* (ITb+IMb) - omega * (RTb) +
        theta_b_v * (SMb+SUb+STb) * ((1-epsilon_b) * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (epsilon_b) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) )
      
      dDUb <- mu_b * rho * IUb
      
      dDTb <- mu_b * zeta * rho * ITb
      
      dDMb <- mu_b * zeta* rho * IMb
      
      dCOa<-sqrt(beta_a)*(SUa+kappa*STa+SMa)*(sqrt(beta_a)*ha*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*(1-ha)*(IUb+IMb+kappa*ITb))
      dCOb<-sqrt(beta_b)*(SUb+kappa*STb+SMb)*(sqrt(beta_a)*(1-hb)*(IUa+IMa+kappa*ITa)+sqrt(beta_b)*hb*(IUb+IMb+kappa*ITb))
      
      dOVa <- theta_a_v * (SUa) * (epsilon_a * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (1-epsilon_a) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) 
      dOVb <- theta_b_v * (SUb) * ((1-epsilon_b) * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + epsilon_b * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) )
      
      dOBa <- theta_a_v * (SUa+SMa+STa) * (epsilon_a * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + (1-epsilon_a) * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) ) 
      dOBb <- theta_b_v * (SUb+SMb+STb) * ((1-epsilon_b) * (DUa+DTa+DMa-lag[8]-lag[9]-lag[10]) + epsilon_b * (DUb+DTb+DMb-lag[18]-lag[19]-lag[20]) )
      
      
      
      return(list(c(dSUa, dSTa, dSMa, dIUa, dITa, dIMa, dRTa, dDUa, dDTa, dDMa,
                    dSUb, dSTb, dSMb, dIUb, dITb, dIMb, dRTb, dDUb, dDTb, dDMb,
                    dCOa, dCOb, dOVa, dOVb, dOBa, dOBb
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

# ls<-full_aware_vax(h=.99, epsilon=.5, mu_a=.02, mu_b=.01, kappa=.05, zeta=0.05, beta=.2, rho=.1, omega=0.01, theta=20, ell=30, phi=.1, I0_a=.001, I0_b=.001, get_params=TRUE)
# sim<-ls$sim
# params<-ls$params
# 
# 
# Smax<-max(c(sim$SUa, sim$STa, sim$SMa, sim$SUb, sim$STb, sim$SMb))
# Imax<-max(c(sim$IUa, sim$ITa, sim$IMa, sim$IUb, sim$ITb, sim$IMb))
# Rmax<-max(c(sim$RTa, sim$RTb))
# Bmax<-max(c(sim$OBa, sim$OBb))
# 
# par(mfcol=c(2, 5), lwd=2)
# 
# 
# plot(SUa~time, data=sim, type="l", col="#3D85BD", ylim=c(0, Smax), ylab="Susceptible", main="a")
# lines(STa~time, data=sim, col="#3D85BD", lty=2)
# lines(SMa~time, data=sim, col="#3D85BD", lty=3)
# plot(SUb~time, data=sim, col="#3D85BD", type="l", ylim=c(0, Smax), main="b", ylab="Susceptible")
# lines(STb~time, data=sim, col="#3D85BD", lty=2)
# lines(SMb~time, data=sim, col="#3D85BD", lty=3)
# 
# plot(IUa~time, data=sim, type="l", col="#7CAA2D", ylim=c(0, Imax), ylab="Infected")
# lines(ITa~time, data=sim, col="#7CAA2D", lty=2)
# lines(IMa~time, data=sim, col="#7CAA2D", lty=3)
# plot(IUb~time, data=sim, col="#7CAA2D", type="l", ylim=c(0, Imax), ylab="Infected")
# lines(ITb~time, data=sim, col="#7CAA2D", lty=2)
# lines(IMb~time, data=sim, col="#7CAA2D", lty=3)
# 
# plot(RTa~time, data=sim, type="l", col="#F5E356", ylim=c(0, Rmax), ylab="Recovered", lty=2)
# plot(RTb~time, data=sim, type="l", col="#F5E356", ylim=c(0, Rmax), ylab="Recovered", lty=2)
# 
# 
# plot(OBa~time, data=sim, type="l", col="purple", ylim=c(0, Bmax), ylab="Vaccinations")
# plot(OBb~time, data=sim, col="purple", type="l", ylim=c(0, Bmax), ylab="Vaccinations")
# 
# sim %>% melt(id="time") %>% filter(variable %in% c("DUa", "DTa", "DMa")) %>% group_by(time) %>% summarise(Da=sum(value)) -> df
# sim %>% melt(id="time") %>% filter(variable %in% c("DUb", "DTb", "DMb")) %>% group_by(time) %>% summarise(Db=sum(value)) %>% right_join(df) -> df
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
