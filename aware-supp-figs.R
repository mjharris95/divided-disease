#### Supplementary Figure 3 ####
library(cowplot)
library(grid)
library(magrittr)
source("aware-eqs.R")

theme_set(
  theme_classic(base_size = 16)
)

#can change the values of h and epsilon. Both vectors can be any length.
h<-c(.5, .99)
epsilon<-c(.5, .99)

expand.grid(h=h, epsilon=epsilon)  %>%
  data.frame() %>%
  apply(1, function(par) full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=par[["h"]],time=200, kappa=.3, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=par[["epsilon"]]) %>% 
          cbind(h=par[["h"]], epsilon=par[["epsilon"]])) %>% 
  do.call(rbind, .) %>%
  mutate(eps_descript = ifelse(epsilon==.5, "Uniform Awareness", 
                               ifelse(epsilon==.99,"Separated Awareness",
                                      " "))) %>%
  mutate(h_descript = ifelse(h==.5, "Uniform \nMixing", 
                             ifelse(h==.99,"Separated \nMixing",
                                    " "))) -> df

png(file="../../Figs/aware-supp3.png",width=7, height=6, res=300, units="in")

ggplot(df)+geom_line(aes(x=time, y=DUa+DPa), size=line_sz, color="#f4a896")+
  geom_line(aes(x=time, y=DUb+DPb, size=as.factor(h)), color="#358597")+
  geom_line(aes(x=time, y=(DUa+DUb+DPa+DPb)/2), color="black")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  #geom_hline(yintercept=c(0,.04,.08, .12), linetype="dotted", color="gray")+
  facet_grid(h+h_descript~epsilon+eps_descript, 
             labeller=label_bquote(
               cols=atop(atop(phantom(), .(eps_descript)), atop("("~epsilon==.(epsilon)~")", phantom())),
               rows=atop(atop(.(h_descript), ""), atop("("~h==.(h)~")", phantom()))))+
  ylab("Cumulative Deaths")+
  theme(strip.text.y = element_text(angle = 0),
        panel.grid.major = element_line(size=.6, color="grey90"),
        panel.grid.minor = element_line(size=.3, color="grey90"))+
  scale_size_manual(values=c("0.5"=.3*line_sz, "0.99"=line_sz), guide="none", na.value = line_sz)+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b)), "full population"),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597", "full population"="black"))+
  theme(strip.text=element_text(size=20))

dev.off()

#### Supplementary Figure 4 ####
library(cowplot)
library(grid)
library(magrittr)
source("aware-eqs.R")

theme_set(
  theme_classic(base_size = 12)
)

#can change the values of h and epsilon. Both vectors can be any length.
h<-c(.5, .9, .97, .99)
epsilon<-c(.5, .75, .9, .99)

expand.grid(h=h, epsilon=epsilon)  %>%
  data.frame() %>%
  apply(1, function(par) full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=par[["h"]],time=200, kappa=.3, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=par[["epsilon"]]) %>% 
          cbind(h=par[["h"]], epsilon=par[["epsilon"]])) %>% 
  do.call(rbind, .) %>%
  mutate(eps_descript = ifelse(epsilon==.5, "Uniform Awareness", 
                               ifelse(epsilon==.99,"Separated Awareness",
                                     " "))) %>%
  mutate(h_descript = ifelse(h==.5, "Uniform \nMixing", 
                             ifelse(h==.99,"Separated \nMixing",
                                    " "))) -> df


pdf(file="../../Figs/aware-supp4.pdf",width=8, height=8)

ggplot(df)+geom_line(aes(x=time, y=IUa+IPa), size=line_sz, color="#f4a896")+geom_line(aes(x=time, y=IUb+IPb, size=as.factor(h)), color="#358597")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  geom_hline(yintercept=c(0,.04,.08, .12), linetype="dotted", color="gray")+
  facet_grid(h+h_descript~epsilon+eps_descript, 
             labeller=label_bquote(
               cols=atop(atop(phantom(), .(eps_descript)), atop("("~epsilon==.(epsilon)~")", phantom())),
               rows=atop(atop(.(h_descript), ""), atop("("~h==.(h)~")", phantom()))))+
  ylab("Infections")+
  theme(strip.text.y = element_text(angle = 0))+scale_y_continuous(breaks=c(0, .04, .08, .12))+
  scale_size_manual(values=c("0.5"=.3*line_sz, "0.75"=.3*line_sz, "0.99"=line_sz), guide="none", na.value = line_sz)+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(strip.text=element_text(size=12))

dev.off()

#### Supplementary Figure 5 ####
source("aware-eqs.R")
library(ggnewscale)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(scales)
library(dplyr)

theme_set(
  theme_classic(base_size = 16)
)

all_eps<-lapply(seq(.5, .99, .01), function(x)  full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=.99,time=200, kappa=.3, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=x) %>% cbind(epsilon=x)) %>% do.call(rbind, .)

all_eps %>% group_by(epsilon) %>% dplyr::summarize(Itota=tail(1-(SUa+SPa),1), Itotb=tail(1-(SUb+SPb),1), Itot=tail(2-(SUa+SPa+SUb+SPb),1)*.5, Imaxa=max(IUa+IPa), Imaxb=max(IUb+IPb), Imax=max(IUa+IPa+IUb+IPb)*.5, Iwhenmaxa=which.max(IUa+IPa), Iwhenmaxb=which.max(IUb+IPb), Iwhenmax=which.max(IUa+IPa+IUb+IPb), Pfina=tail(SPa+IPa+RPa,1), Pfinb=tail(SPb+IPb+RPb,1), Pfin=.5*tail(SPa+IPa+RPa+SPb+IPb+RPb,1), sanitize=FALSE) %>% ungroup() -> sum_stats_eps

all_h<-lapply(seq(.5, .99, .01), function(x)  full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=x,time=200, kappa=.3, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=.5) %>% cbind(h=x)) %>% do.call(rbind, .)

all_h %>% group_by(h) %>% dplyr::summarize(Itota=tail(1-(SUa+SPa),1), Itotb=tail(1-(SUb+SPb),1), Itot=tail(2-(SUa+SPa+SUb+SPb),1)*.5, Imaxa=max(IUa+IPa), Imaxb=max(IUb+IPb), Imax=max(IUa+IPa+IUb+IPb)*.5, Iwhenmaxa=which.max(IUa+IPa), Iwhenmaxb=which.max(IUb+IPb), Iwhenmax=which.max(IUa+IPa+IUb+IPb), Pfina=tail(SPa+IPa+RPa,1), Pfinb=tail(SPb+IPb+RPb,1), Pfin=.5*tail(SPa+IPa+RPa+SPb+IPb+RPb,1), sanitize=FALSE) %>% ungroup() -> sum_stats_h

Itot_lims<-range(sum_stats_eps$Itot,sum_stats_eps$Itota,sum_stats_eps$Itotb,
                 sum_stats_h$Itot,sum_stats_h$Itota,sum_stats_h$Itotb)
Imax_lims<-range(sum_stats_eps$Imax,sum_stats_eps$Imaxa,sum_stats_eps$Imaxb,
                 sum_stats_h$Imax,sum_stats_h$Imaxa,sum_stats_h$Imaxb)
Iwhenmax_lims<-range(sum_stats_eps$Iwhenmax,sum_stats_eps$Iwhenmaxa,sum_stats_eps$Iwhenmaxb,
                     sum_stats_h$Iwhenmax,sum_stats_h$Iwhenmaxa,sum_stats_h$Iwhenmaxb)


Itot_eps<-ggplot(sum_stats_eps)+
  geom_line(aes(x=epsilon, y=Itota, color="group *a*"))+
  geom_line(aes(x=epsilon, y=Itotb, color="group *b*"))+
  geom_line(aes(x=epsilon, y=Itot, color="full population"))+
  ylab("Total Infections")+
  ylim(Itot_lims)+
  xlab(expression(paste("Awareness Separation (",epsilon, ")")))+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b)), "full population"),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597", "full population"="black"))+
  theme(legend.position="bottom")

leg<-get_legend(Itot_eps)
Itot_eps <- Itot_eps+theme(legend.position="none")

Imax_eps<-ggplot(sum_stats_eps)+
  ylim(Imax_lims)+
  geom_line(aes(x=epsilon, y=Imaxa), color="#f4a896")+
  geom_line(aes(x=epsilon, y=Imaxb), color="#358597")+
  geom_line(aes(x=epsilon, y=Imax), color="black")+
  ylab("Peak Infections")+
  xlab(expression(paste("Awareness Separation (",epsilon, ")")))

Iwhenmax_eps<-ggplot(sum_stats_eps)+
  ylim(Iwhenmax_lims)+
  geom_line(aes(x=epsilon, y=Iwhenmaxa), color="#f4a896")+
  geom_line(aes(x=epsilon, y=Iwhenmaxb), color="#358597")+
  geom_line(aes(x=epsilon, y=Iwhenmax), color="black")+
  ylab("Peak Infection Date")+
  xlab(expression(paste("Awareness Separation (",epsilon, ")")))

Itot_h<-ggplot(sum_stats_h)+
  ylim(Itot_lims)+
  geom_line(aes(x=h, y=Itot), color="#f4a896", size=1.5)+
  geom_line(aes(x=h, y=Itotb), color="#358597")+
  geom_line(aes(x=h, y=Itot), color="black")+
  ylab("Total Infections")+
  xlab("Mixing Separation (h)")

Imax_h<-ggplot(sum_stats_h)+
  ylim(Imax_lims)+
  geom_line(aes(x=h, y=Imaxa), color="#f4a896")+
  geom_line(aes(x=h, y=Imaxb), color="#358597")+
  geom_line(aes(x=h, y=Imax), color="black")+
  ylab("Peak Infections")+
  xlab("Mixing Separation (h)")

Iwhenmax_h<-ggplot(sum_stats_h)+
  ylim(Iwhenmax_lims)+
  geom_line(aes(x=h, y=Iwhenmaxa), color="#f4a896")+
  geom_line(aes(x=h, y=Iwhenmaxb), color="#358597")+
  geom_line(aes(x=h, y=Iwhenmax), color="black")+
  ylab("Peak Infection Date")+
  xlab("Mixing Separation (h)")

plot_all<-plot_grid(Itot_eps, Itot_h, Imax_eps, Imax_h, 
                    Iwhenmax_eps, Iwhenmax_h, 
                    nrow=3, labels="AUTO",
                    scale=.9,
                    hjust=0, label_size=18)

pdf(file="../../Figs/aware-supp5.pdf",width=7, height=10)
plot_grid(plot_all, leg, nrow=2, rel_heights=c(9,1))
dev.off()

#### Supplementary Figure 6-8 Function ####
library(cowplot)
library(grid)
library(magrittr)
source("aware-eqs.R")

theme_set(
  theme_classic(base_size = 12)
)

fig1_diff_init<-function(beta_a=.2, beta_b=.2, rho_a=.1, rho_b=.1, mu_a=.01, mu_b=.01,
                         line_scale=FALSE, plot_type="I"){
  
  expand.grid(h=c(.5, .99), epsilon=c(.5,.99))  %>%
    data.frame() %>%
    apply(1, function(par) full_aware(beta_a=beta_a, beta_b=beta_b, 
                                      rho_a=rho_a, rho_b=rho_b,
                                      mu_a=mu_a, mu_b=mu_b,
                                      theta=100, time=200, kappa=.3, 
                                      I0_a=.001/2, I0_b=.001/2, ell=1, phi=0,
                                      epsilon=par[["epsilon"]], h=par[["h"]]) %>%
            mutate(Da=DUa+DPa-lag(DUa)-lag(DPa),
                   Db=DUb+DPb-lag(DUb)-lag(DPb)) %>%
            cbind(h= par[["h"]], epsilon=par[["epsilon"]])) %>%
    do.call(rbind, .) %>%
    mutate(eps_descript = ifelse(epsilon==.5, "Uniform", "Separated")) %>%
    mutate(h_descript = ifelse(h==.5, "Uniform", "Separated")) -> df
  
  scale_sz<-ifelse(line_scale, .3, 1)
  
  #Infections plot
  if(plot_type=="I"){
    this_plot<-ggplot(df)+geom_line(aes(x=time, y=IUa+IPa), size=.8, color="#f4a896")+geom_line(aes(x=time, y=IUb+IPb, size=as.factor(epsilon)), color="#358597")+
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
      geom_hline(yintercept=c(0,.04,.08, .12), linetype="dotted", color="gray")+
      facet_grid(h+h_descript~epsilon+eps_descript, 
                 labeller=label_bquote(
                   cols=atop(atop(phantom(), .(eps_descript)~"Awareness"), atop("("~epsilon==.(epsilon)~")", phantom())),
                   rows=atop(atop(.(h_descript), "Mixing"), atop("("~h==.(h)~")", phantom()))))+
      ylab("Infections")+
      theme(strip.text.y = element_text(angle = 0))+scale_y_continuous(breaks=c(0, .04, .08, .12))+
      scale_size_manual(values=c("0.5"=scale_sz*.8, "0.99"=.8), guide="none")+
      scale_colour_manual(name="",
                          labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                          values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
      theme(strip.text=element_text(size=14))
  }
  
  #Deaths plot
  else{
    
    this_plot<-ggplot(df)+geom_line(aes(x=time, y=Da), size=.8, color="#f4a896")+geom_line(aes(x=time, y=Db, size=as.factor(epsilon)), color="#358597")+
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
      geom_hline(yintercept=c(0,5e-5,1e-4), linetype="dotted", color="gray")+
      facet_grid(h+h_descript~epsilon+eps_descript, 
                 labeller=label_bquote(
                   cols=atop(atop(phantom(), .(eps_descript)~"Awareness"), atop("("~epsilon==.(epsilon)~")", phantom())),
                   rows=atop(atop(.(h_descript), "Mixing"), atop("("~h==.(h)~")", phantom()))))+
      ylab("Deaths")+
      theme(strip.text.y = element_text(angle = 0))+ scale_y_continuous(breaks=c(0,5e-5,1e-4))+
      scale_size_manual(values=c("0.5"=scale_sz*.8, "0.99"=.8), guide="none")+
      scale_colour_manual(name="", 
                          labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                          values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
      theme(strip.text=element_text(size=14))
  }
  
  this_plot
}

#### Supplementary Figure 6 ####
betaI<-fig1_diff_init(beta_a=.21, beta_b=.19,plot_type="I")+ggtitle(expression(paste("Transmission Coefficient (", beta[a], " = 0.21; ", beta[b], " = 0.19)")))
betaD<-fig1_diff_init(beta_a=.21, beta_b=.19,plot_type="D")+ggtitle(" ")

pdf(file="../../Figs/aware-supp6.pdf", height=6.5, width=5)
plot_grid(betaI, betaD, nrow=2, labels="AUTO", label_y=.8)
dev.off()

#### Supplementary Figure 7 ####
rhoI<-fig1_diff_init(rho_a=.09, rho_b=.11, plot_type = "I")+ggtitle(expression(paste("Infectious Period (1/", rho[a], " = 11.11; 1/", rho[b], " = 9.09)")))
rhoD<-fig1_diff_init(rho_a=.09, rho_b=.11, plot_type = "D")+ggtitle(" ")

pdf(file="../../Figs/aware-supp7.pdf", height=6.5, width=5)
plot_grid(rhoI, rhoD, nrow=2, labels="AUTO", label_y=.8, align="h")
dev.off()

#### Supplementary Figure 8 ####
muI<-fig1_diff_init(mu_a=.015, mu_b=.0067,plot_type="I", line_scale=TRUE)+ggtitle(expression(paste("Infection Fatality Rate (", mu[a], " = 0.015; ", mu[b], " = 0.0067)")))
muD<-fig1_diff_init(mu_a=.015, mu_b=.0067,plot_type="D")+ggtitle(" ")

pdf(file="../../Figs/aware-supp8.pdf", height=6.5, width=5)
plot_grid(muI, muD, nrow=2, labels="AUTO", label_y=.8)
dev.off()

#### Supplementary Figure 9 ####
library(cowplot)
library(grid)
library(magrittr)
library(ggtext)
source("aware-eqs.R")

theme_set(
  theme_classic(base_size = 12)
)


beta<-seq(from=.22, to=.55, by=.11/2)
kappa<-.2*(0:5)
theta<-c(0, 100, 500, 1000)

beta_ts<-c(.22, .55)
theta_ts<-c(0,100, 500, 1000)

expand.grid(beta=beta, kappa=kappa, theta=theta)  %>%
  data.frame() %>%
  apply(1, function(par) full_aware(theta=par[["theta"]], beta=par[["beta"]], rho=.1, mu=.01, h=.99, time=200, kappa=par[["kappa"]], I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=.5) %>%
          filter(time==200) %>%
          mutate(rat=(1-(SUa+SPa))/(1-(SUb+SPb))) %>%
          select(rat) %>%
          cbind(beta=par[["beta"]], kappa=par[["kappa"]], theta=par[["theta"]])) %>% 
  do.call(rbind, .) %>%
  mutate(theta_text=plyr::mapvalues(theta, 
                                    from=c(0,100,500,1000),
                                    to=c("No Resp.", "Low Resp.", "Med Resp.", "High Resp.")))-> df

ratios_plot<-ggplot()+
  geom_line(data=df, aes(y=rat, x=beta, group=kappa, colour=as.factor(kappa), size=ifelse(kappa==.2, ".2", "false")))+
  geom_point(data=df %>% filter(beta %in% beta_ts & theta %in% theta_ts & kappa==.2), 
             aes(y=rat, x=beta))+
  geom_text(data=df %>% filter(beta %in% beta_ts & theta %in% theta_ts & kappa==.2), 
            aes(y=rat, x=beta, label=paste('beta', "==", beta)), 
            nudge_y=.3, parse=TRUE)+
  facet_wrap(~theta+theta_text, nrow=4,
             labeller=label_bquote(
               rows=atop(.(theta_text), "("~theta==.(theta)~")"))
  )+
  scale_size_manual(values=c(.9, .6))+
  theme(legend.position="bottom", legend.text = element_markdown())+
  ylab("Ratio of Cumulative Infections (a/b)")+
  xlab(expression(paste("Transmission Coefficient (", beta, ")")))+
  scale_color_grey(expression(paste("Protection Efficacy (", kappa, ")")),
                   labels=c("0", "**0.2**", "0.4", "0.6", "0.8", "1"))+
  guides(color=guide_legend(title.position="top", override.aes=list(size=5)),
         size="none")+
  xlim(c(.17, .6))

#fix kappa = .2, grid of beta =.11, .22, 1.1 and theta = 200, 600, 1000
#fix all labels



expand.grid(beta=beta_ts, theta=theta_ts)  %>%
  data.frame() %>%
  apply(1, function(par) full_aware(theta=par[["theta"]], beta=par[["beta"]], rho=.1, mu=.01, h=.99, time=200, kappa=.2, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=.5) %>%
          cbind(beta=par[["beta"]], theta=par[["theta"]])) %>% 
  do.call(rbind, .) %>%
  mutate(theta_text=plyr::mapvalues(theta, 
                                    from=c(0,100,500,1000),
                                    to=c("No Resp.", "Low Resp.", "Med Resp.", "High Resp.")))-> df_ts


ts_plot<-ggplot(df_ts)+geom_line(aes(x=time, y=IUa+IPa, color="group *a*"))+geom_line(aes(x=time, y=IUb+IPb, color="group *b*"))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  ylab("Infections")+
  facet_grid(theta+theta_text~beta, 
             labeller=label_bquote(
               rows=atop(.(theta_text), "("~theta==.(theta)~")"),
               col=beta==.(beta)))+
  theme(strip.text.y = element_text(angle = 0))+
  geom_hline(yintercept=c(0, .1, .2, .3, .4, .5), linetype="dotted", color="gray")+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(strip.text=element_text(size=12), legend.position="bottom")

plots<-plot_grid(ratios_plot+theme(legend.position="none"),
                 ts_plot+theme(legend.position="none"),
                 ncol=2, labels="AUTO", rel_widths=c(1,2))

ratios_legend<-get_legend(ratios_plot)
ts_legend<-get_legend(ts_plot)

legends<-plot_grid(ratios_legend,
                   ts_legend,
                   ncol=2, rel_widths=c(1,2))

pdf(file="../../Figs/aware-supp9.pdf",width=8, height=11)
plot_grid(plots, legends, nrow=2, rel_heights=c(10,1))
dev.off()

#### Supplementary Figure 10 ####
source("aware-eqs.R")
library(cowplot)

theme_set(
  theme_classic(base_size = 16)
)

#can change the values of epsilon here. Vectors must have length 2.
epsilon_vals<-c(.75, .9)
epsilon_labs<-ifelse(epsilon_vals==.5, "Uniform", 
                     ifelse(epsilon_vals==.99,"Separated",
                            ""))

long_uni<-full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=0.99, time=1000, I0_a=.001, I0_b=0, kappa=.3, ell=30, phi=.02, epsilon=epsilon_vals[1]) 

long_sep<-full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=0.99, time=1000, I0_a=.001, I0_b=0, kappa=.3, ell=30, phi=.02, epsilon=epsilon_vals[2]) 

long_sep %>% 
  ggplot()+
  geom_line(aes(x=time, y=IUa+IPa, color="group *a*"), size=.8)+
  geom_line(aes(x=time, y=IUb+IPb, color="group *b*"), size=.8)+
  ylab("Infections")+
  ggtitle(expr(paste("B. ", !!epsilon_labs[2], " (",epsilon," = ", !!epsilon_vals[2], ")")))+  
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(legend.position="bottom",
        title=element_text(size=12))+
  scale_y_continuous(limits=range(-.002, long_uni$IUa+long_uni$IPa),
                     breaks=c(0, .02, .04)) -> long_sep_plot

long_leg<-get_legend(long_sep_plot)

long_sep_plot<-long_sep_plot+theme(legend.position="none")

pdf(file="../../Figs/aware-supp10.pdf",width=6, height=5.5)

long_uni %>% 
  ggplot()+
  geom_line(aes(x=time, y=IUa+IPa, color="group *a*"), size=.8, linetype="solid")+
  geom_line(aes(x=time, y=IUb+IPb, color="group *b*"), size=.8, linetype="solid")+
  ylab("Infections")+
  ggtitle(expr(paste("A. ", !!epsilon_labs[1], " (",epsilon," = ", !!epsilon_vals[1], ")")))+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(legend.position="none",
        title=element_text(size=12))+
  scale_y_continuous(limits=range(-.002, long_uni$IUa+long_uni$IPa),
                     breaks=c(0, .02, .04), expand=c(0,0)) -> long_uni_plot

long_plots<-plot_grid(long_uni_plot, long_sep_plot, nrow=2)
plot_grid(long_plots, long_leg, nrow=2, rel_heights=c(10, 1))

dev.off()

#### Supplementary Figure 11 ####
epsilon_vals<-c(.5, .99)
epsilon_labs<-ifelse(epsilon_vals==.5, "Uniform Awareness", 
                     ifelse(epsilon_vals==.99,"Separated Awareness",
                            ""))

long_uni<-full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=0.99, time=1000, I0_a=.001, I0_b=0, kappa=.3, ell=30, phi=.02, epsilon=epsilon_vals[1]) 

long_sep<-full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=0.99, time=1000, I0_a=.001, I0_b=0, kappa=.3, ell=30, phi=.02, epsilon=epsilon_vals[2]) 

theme_set(
  theme_classic(base_size = 12)
)


ggplot()+
  geom_line(data=long_uni, aes(x=time, y=(IUa+IPa+IUb+IPb)/2), linetype="dashed")+
  geom_line(data=long_sep, aes(x=time, y=(IUa+IPa+IUb+IPb)/2), linetype="solid")+
  ylab("Infections") +
  ggtitle("A.") -> long_I

ggplot()+
  geom_line(data=long_uni, aes(x=time, y=(DUa+DPa+DUb+DPb)/2, linetype="0.5"))+
  geom_line(data=long_sep, aes(x=time, y=(DUa+DPa+DUb+DPb)/2, linetype="0.99"))+
  ylab("Cumulative Deaths") +
  ggtitle("B.")+
  scale_linetype_manual("",
                        labels=c(expression(paste("Uniform Awareness (",epsilon, "=0.5)")),
                                 expression(paste("Separated Awareness (", epsilon, "=0.99)"  ))),
                        values=c("0.5"="dashed", "0.99"="solid"))+
  theme(legend.position="bottom") -> long_D


pdf(file="../../Figs/aware-supp11.pdf",width=6, height=5.5)

plot_grid(long_I, long_D, nrow=2)

dev.off()

#### Supplementary Figure 12 ####
library(cowplot)
source("split-vax-eqs.R")
source("aware-vax-process.R")

source("aware-eqs.R")
theme_set(
  theme_classic(base_size = 16)
)


long_uni_vax<-full_aware_vax(theta=40, beta=.2, rho=.1, mu_a=.02, mu_b=.01, zeta=.05, kappa=.05, h=.99, time=2000, v_start=200, I0_a=.001/2, I0_b=0.001/2, phi=0.01, omega=0.01, ell=30, epsilon=.5)

vax_df<-lapply(c(.5, .99), function(x) full_aware_vax(theta=40, beta=.2, rho=.1, mu_a=.02, mu_b=.01, zeta=.05, kappa=.05, h=.99, time=2000, v_start=200, I0_a=.001/2, I0_b=0.001/2, phi=0.01, omega=0.01, ell=30, epsilon=x) %>% 
                 mutate(Ia=IUa+ITa+IMa, 
                        Ib=IUb+ITb+IMb,
                        Da=DUa+DTa+DMa,
                        Db=DUb+DTb+DMb) %>%
                 mutate(Da=Da-lag(Da),
                        Db=Db-lag(Db)) %>%
                 mutate(epsilon=x)) %>% do.call(rbind, .) %>% 
  select(time, Ia, Ib, Da, Db, epsilon) %>%
  pivot_longer(!c(time, epsilon), values_to="y", names_to="cat") %>%
  mutate(group = substring(cat, 2),
         epi=substring(cat, 1,1)) %>%
  mutate(eps_lab=ifelse(epsilon==.5, "Uniform", "Separated"))


Ivax_dyn<-ggplot()+
  geom_line(data=vax_df %>% filter(group=="a" & epi== "I"), aes(x=time, y=y, color="group *a*"), size=.8)+
  geom_line(data=vax_df %>% filter(group=="b" & epi== "I"), aes(x=time, y=y, color="group *b*"), size=.8)+
  ylab("Infections")+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(legend.position="bottom")+
  annotate("segment", x = 200, xend = 200, y =.1, yend = .01,
           arrow = arrow(type="closed", length=unit(0.30,"cm"), angle=20),
           size=.8)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  annotate("text", x=200, y=.11, label="Vaccination", hjust=0)+
  facet_grid(~eps_lab+epsilon, 
             labeller=label_bquote(
               cols=.(eps_lab)~"Awareness ("*epsilon==.(epsilon)*")",
               rows=NULL))+theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ylim(c(0, .2))

Dvax_dyn<-ggplot()+
  geom_line(data=vax_df %>% filter(group=="a" & epi== "D"), aes(x=time, y=y, color="group *a*"), size=.8)+
  geom_line(data=vax_df %>% filter(group=="b" & epi== "D"), aes(x=time, y=y, color="group *b*"), size=.8)+
  ylab(expression(paste("Deaths (", log[10], ")")))+
  scale_y_log10(limits=c(4e-7,1e-3))+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(legend.position="bottom",
        plot.title=element_text(size=16))+
  annotate("segment", x = 200, xend = 200, y=9e-5, yend = 5e-6,
           arrow = arrow(type="closed", length=unit(0.30,"cm"), angle=20),
           size=.8)+
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=0, yend=Inf, size=1)+
  annotate("text", x=200, y=2e-4, label="Vaccination", hjust=0)+
  facet_grid(~eps_lab+epsilon)+
  theme(strip.text.y = element_blank(),
        strip.text.x = element_blank())

pdf(file="../../Figs/aware-supp12.pdf",width=6.7, height=6.7)
plot_grid(Ivax_dyn+theme(legend.position="none"), Dvax_dyn, get_legend(Ivax_dyn), nrow=3, align="hv", rel_heights = c(6,8,1))
dev.off()

#### Supplementary Figures 13 and 14 ####
source("split-vax-eqs.R")
source("aware-vax-process.R")
library(magrittr)

#can change the values of epsilon here. Vectors must have length 2.
epsilon_vals<-c(.5, .99)

options(dplyr.summarise.inform = FALSE)

full_df<-lapply(c(0, 50, 100, 200), function(vtime) lapply(epsilon_vals, function(epsilon) lapply(c(.05, seq(from=0, to=1, by=.1)), function(x) full_aware_vax(theta=40, beta=.2, rho=.1, mu_a=.02, mu_b=.01, zeta=x, kappa=x, h=.99, time=2200, v_start=vtime, I0_a=.001/2, I0_b=0.001/2, phi=0.01, omega=0.01, ell=30, epsilon=epsilon) %>% 
                                                                                                    vaxdf_process() %>%
                                                                                                    cbind(veff=x, epsilon=epsilon, vtime=vtime)) %>%
                                                             do.call(rbind, .)) %>%
                  do.call(rbind, .)) %>%
  do.call(rbind, .)

full_df %>%
  group_by(dis, time, group, veff, epsilon, vtime) %>%
  summarize(y=sum(y)) -> dis_df

full_df %>%
  group_by(imm, time, group, veff, epsilon, vtime) %>%
  summarize(y=sum(y)) -> imm_df

library(ggforce)

theme_set(
  theme_classic(base_size = 10)
)

vtime_sum_plots<-list()
vtime_dyn_plots<-list()

v_start_vals<-c(0, 50, 100, 200)

for(i in 1:4){
  
  v_start_val<-v_start_vals[i]
  
  Dvax_lims<-range(filter(dis_df, veff==0.05 & dis=="D" & epsilon %in% c(.5, .99))$y)
  
  Dvax_dyn<-ggplot()+
    geom_line(data=dis_df %>% filter(group=="a" & vtime==v_start_val & veff==0.05 & dis=="D" & epsilon %in% c(.5, .99)), 
              aes(x=time, y=y, color="group *a*", linetype=as.factor(epsilon)), size=.8)+
    geom_line(data=dis_df %>% filter(group=="b" & vtime==v_start_val & veff==0.05 & dis=="D" & epsilon %in% c(.5, .99)), 
              aes(x=time, y=y, color="group *b*", linetype=as.factor(epsilon)), size=.8)+
    ylab("Cumulative Deaths")+
    xlab("Time")+
    scale_colour_manual(name="", 
                        labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                        values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
    scale_linetype_manual(name="", 
                          labels=c(expression(paste("Uniform Awareness (", epsilon, " = 0.5)")),
                                   expression(paste("Separated Awareness (", epsilon, " = 0.99)"))),
                          values=c("0.5"="dashed", "0.99"="solid"))+
    guides(linetype = guide_legend(override.aes = list(color="black", size=.5)))+
    theme(legend.position="bottom")+    
    ylim(Dvax_lims)+
    geom_vline(xintercept=v_start_val, linetype="dashed", color="grey40")+
    annotate("text", x=v_start_val, label=paste0("vaccination at t = ", v_start_val[i]), 
             y=0, hjust=0, size=3)
  
  legend_dyn<-get_legend(Dvax_dyn)
  Dvax_dyn<-Dvax_dyn+theme(legend.position="none")
  
  inf_tot<-filter(dis_df,  time %in% c(0,2200) & dis=="C") %>%
    group_by(group, veff, epsilon, vtime) %>%
    summarize(y=diff(range(y))) 
  
  I_plot<-veff_plot(filter(inf_tot, vtime==v_start_val), legend=TRUE)+
    ylim(range(inf_tot$y))+
    ylab("Cumulative Infections")
  
  legend_sum<-get_legend(I_plot)
  
  I_plot<-I_plot+theme(legend.position="none")#+ylim(c(1, 3.6))
  
  death_tot<-filter(dis_df, time %in% c(0,2000) & dis=="D") %>%
    group_by(group, veff, epsilon, vtime) %>%
    summarize(y=diff(range(y))) 
  
  D_plot<-veff_plot(death_tot %>% filter(vtime==v_start_val))+
    #scale_y_continuous(trans='log2', 
    #                   limits = c(.001, .018),
    #                   labels = scales::label_number(accuracy = 0.001))+
    theme(legend.position="none",
          panel.grid.major = element_line(size = .2, colour = NA),
          panel.grid.minor = element_line(size = .2, colour = NA))+
    ylim(range(death_tot$y))+
    ylab("Cumulative Deaths")+
    geom_text(aes(x=.55, y=.039), label=paste("Vaccination\n at t =", v_start_vals[i]))
  
  B_tot<-filter(imm_df, time==2200 & imm=="B") %>%
    select(group, y, veff, epsilon, vtime)
  
  B_plot<-veff_plot(B_tot %>% filter(vtime==v_start_val), is_B_plot=TRUE)+
    ylim(range(B_tot$y))+
    ylab("Cumulative Vaccinations")
  
  vtime_dyn_plots[[i]]<-Dvax_dyn
  vtime_sum_plots[[(i-1)*3+1]]<-D_plot
  vtime_sum_plots[[(i-1)*3+2]]<-B_plot
  vtime_sum_plots[[(i-1)*3+3]]<-I_plot
}


sum_grid<-plot_grid(plotlist=vtime_sum_plots,
                    nrow=3,
                    # greedy=TRUE,
                    hjust=0,
                    align="v",
                    scale=1,
                    label_x=.3,
                    label_size=14,
                    byrow=FALSE,
                    labels="AUTO")

dyn_grid<-plot_grid(plotlist=vtime_dyn_plots, 
                    hjust=0, align="v", scale=1, ncol=4,
                    byrow=FALSE,
                    labels="AUTO")

pdf(file="../../Figs/aware-supp13.pdf",width=8, height=4)

plot_grid(dyn_grid, legend_dyn, nrow=2, rel_heights=c(7,1))

dev.off()

pdf(file="../../Figs/aware-supp14.pdf",width=8, height=10)

plot_grid(sum_grid, legend_sum, nrow=2, rel_heights=c(7,1))

dev.off()

