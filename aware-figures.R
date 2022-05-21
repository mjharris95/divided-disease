#### Figure 1 #### 

library(cowplot)
library(grid)
library(magrittr)
source("aware-eqs.R")

theme_set(
  theme_classic(base_size = 16)
)

h<-c(.5, .99)
epsilon<-c(.5, .99)

expand.grid(h=h, epsilon=epsilon)  %>%
  data.frame() %>%
  apply(1, function(par) full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=par[["h"]],time=200, kappa=.3, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=par[["epsilon"]]) %>% 
          cbind(h=par[["h"]], epsilon=par[["epsilon"]])) %>% 
  do.call(rbind, .) %>%
  mutate(eps_descript = ifelse(epsilon==.5, "Uniform", "Separated")) %>%
  mutate(h_descript = ifelse(h==.5, "Uniform", "Separated")) -> df


line_sz<-.8

ggplot(df)+geom_line(aes(x=time, y=IUa+IPa), size=line_sz, color="maroon")+geom_line(aes(x=time, y=IUb+IPb, size=as.factor(h)), color="blue")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  geom_hline(yintercept=c(0,.04,.08, .12), linetype="dotted", color="gray")+
  facet_grid(h+h_descript~epsilon+eps_descript, 
             labeller=label_bquote(
               cols=atop(atop(phantom(), .(eps_descript)~"Awareness"), atop("("~epsilon==.(epsilon)~")", phantom())),
               rows=atop(atop(.(h_descript), "Mixing"), atop("("~h==.(h)~")", phantom()))))+
  ylab("Infections")+
  theme(strip.text.y = element_text(angle = 0))+scale_y_continuous(breaks=c(0, .04, .08, .12))+
  scale_size_manual(values=c("0.5"=.3*line_sz, "0.99"=line_sz), guide="none")+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="maroon", "group *b*"="blue"))+
  theme(strip.text=element_text(size=20))

#### Figure 2 ####
source("aware-eqs.R")
library(cowplot)

theme_set(
  theme_classic(base_size = 16)
)

long_uni<-full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=0.99, time=1000, I0_a=.001, I0_b=0, kappa=.3, ell=30, phi=.02, epsilon=.5) 

long_sep<-full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=0.99, time=1000, I0_a=.001, I0_b=0, kappa=.3, ell=30, phi=.02, epsilon=.99) 

long_sep %>% 
  ggplot()+
  geom_line(aes(x=time, y=IUa+IPa, color="group *a*"), size=.8)+
  geom_line(aes(x=time, y=IUb+IPb, color="group *b*"), size=.8)+
  ylab("Infections")+
  ggtitle(expression(paste("B. Separated Awareness (",epsilon," = 0.99)")))+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="maroon", "group *b*"="blue"))+
  theme(legend.position="bottom",
        title=element_text(size=12))+
  scale_y_continuous(limits=range(-.002, long_uni$IUa+long_uni$IPa),
                     breaks=c(0, .02, .04)) -> long_sep_plot

long_leg<-get_legend(long_sep_plot)

long_sep_plot<-long_sep_plot+theme(legend.position="none")

long_uni %>% 
  ggplot()+
  geom_line(aes(x=time, y=IUa+IPa, color="group *a*"), size=.8, linetype="solid")+
  geom_line(aes(x=time, y=IUb+IPb, color="group *b*"), size=.8, linetype="solid")+
  ylab("Infections")+
  ggtitle(expression(paste("A. Uniform Awareness (",epsilon," = 0.5)")))+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="maroon", "group *b*"="blue"))+
  theme(legend.position="none",
        title=element_text(size=12))+
  scale_y_continuous(limits=range(-.002, long_uni$IUa+long_uni$IPa),
                     breaks=c(0, .02, .04), expand=c(0,0)) -> long_uni_plot

long_plots<-plot_grid(long_uni_plot, long_sep_plot, nrow=2)
plot_grid(long_plots, long_leg, nrow=2, rel_heights=c(10, 1))

#### Figure 3 ####
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
                      values = c("group *a*"="maroon", "group *b*"="blue"))+
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
                      values = c("group *a*"="maroon", "group *b*"="blue"))+
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


plot_grid(Ivax_dyn+theme(legend.position="none"), Dvax_dyn, get_legend(Ivax_dyn), nrow=3, align="hv", rel_heights = c(6,8,1))

#### Figure 4 ####
source("split-vax-eqs.R")
source("aware-vax-process.R")

options(dplyr.summarise.inform = FALSE)

full_df<-lapply(c(.5, .99), function(epsilon) lapply(c(1/(2^c((0:14)*.5)),.6, .7, .8), function(x) full_aware_vax(theta=40, beta=.2, rho=.1, mu_a=.02, mu_b=.01, zeta=x, kappa=x, h=.99, time=2000, v_start=200, I0_a=.001/2, I0_b=0.001/2, phi=0.01, omega=0.01, ell=30, epsilon=epsilon) %>% 
                                                       vaxdf_process() %>%
                                                       cbind(veff=x, epsilon=epsilon)) %>%
                  do.call(rbind, .)) %>%
  do.call(rbind, .)


full_df %>%
  group_by(dis, time, group, veff, epsilon) %>%
  summarize(y=sum(y)) -> dis_df

full_df %>%
  group_by(imm, time, group, veff, epsilon) %>%
  summarize(y=sum(y)) -> imm_df

library(ggforce)

theme_set(
  theme_classic(base_size = 10)
)

inf_tot<-filter(dis_df,  time %in% c(200,2000) & dis=="C") %>%
  group_by(group, veff, epsilon) %>%
  summarize(y=diff(range(y))) 

I_plot<-veff_plot(inf_tot, legend=TRUE)+ylim(c(.4, 3.3))
legend_epi<-get_legend(I_plot)

I_plot<-I_plot+theme(legend.position="none")#+ylim(c(1, 3.6))

death_tot<-filter(dis_df, time %in% c(200,2000) & dis=="D") %>%
  group_by(group, veff, epsilon) %>%
  summarize(y=diff(range(y))) 

D_plot<-veff_plot(death_tot)+
  #scale_y_continuous(trans='log2', 
  #                   limits = c(.001, .018),
  #                   labels = scales::label_number(accuracy = 0.001))+
  theme(legend.position="none",
        panel.grid.major = element_line(size = .2, colour = NA),
        panel.grid.minor = element_line(size = .2, colour = NA))

B_tot<-filter(imm_df, time==2000 & imm=="B") %>%
  select(group, y, veff, epsilon)

B_plot<-veff_plot(B_tot, is_B_plot=TRUE)#+ylim(c(.66, 8))


ADplot<-plot_grid(D_plot, B_plot, I_plot, nrow=1,
                  labels=c("A. Deaths",
                           "B. Vaccinations ", "C. Infections" 
                  ),
                  # greedy=TRUE,
                  hjust=0,
                  align="v",
                  scale=1,
                  label_x=.3,
                  label_size=14
)

plot_grid(ADplot, legend_epi, nrow=2, rel_heights=c(7,2))

#### Supplementary Figure 1 ####
li_size1<-.8
li_size2<-1.5


h<-c(.5, .99)
epsilon<-c(.5, .99)



expand.grid(h=h, epsilon=epsilon)  %>%
  data.frame() %>%
  apply(1, function(par) full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=par[["h"]],time=200, kappa=.3, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=par[["epsilon"]]) %>% 
          cbind(h=par[["h"]], epsilon=par[["epsilon"]])) %>% 
  do.call(rbind, .) %>%
  mutate(eps_descript = ifelse(epsilon==.5, "Uniform", "Separated")) %>%
  mutate(h_descript = ifelse(h==.5, "Uniform", "Separated")) -> df


df %>% filter(h==.99, time<=80) %>%
  mutate(Pb=SPb+IPb+RPb,
         Pa=SPa+IPa+RPa,
         Ib=1-(SPb+SUb),
         Ia=1-(SPa+SUa)) -> t80

Pmax<-max(c(t80$Pa, t80$Pb))
Imax<-max(c(t80$Ia, t80$Ib))


Pb_plot<-ggplot(t80) + 
  geom_line(aes(x=time, y=Pb, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="lightblue")+
  scale_linetype_manual("Awareness",
                        labels=c(expression(paste("Uniform (",epsilon, "=0.5)")),
                                 expression(paste("Separated (", epsilon, "=0.99)"  ))),
                        values=c("0.5"="dashed", "0.99"="solid"))+
  scale_size_manual(expression(paste("Awareness (", epsilon, ")")),
                    labels=c("Uniform", "Separated"),
                    values=c("0.5"=li_size2, "0.99"=li_size1),
                    guide="none")+
  ylab("Protective")+
  theme(legend.position="bottom")+
  ylim(c(0, Pmax))+
  guides(linetype = guide_legend(override.aes = list(color="black", size=.5)))

legend<-get_legend(Pb_plot)

Pb_plot<-Pb_plot+theme(legend.position="none")

Ib_plot<-ggplot(t80) + 
  geom_line(aes(x=time, y=Ib, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="blue")+
  scale_linetype_manual(expression(paste("Awareness (", epsilon, ")")),
                        labels=c("Uniform", "Separated"),
                        values=c("0.5"="dashed", "0.99"="solid"))+
  scale_size_manual(expression(paste("Awareness (", epsilon, ")")),
                    labels=c("Uniform", "Separated"),
                    values=c("0.5"=li_size2, "0.99"=li_size1))+
  ylab("Infections")+
  theme(legend.position="none")+
  ylim(c(0, Imax))



Pa_plot<-ggplot(t80) + 
  geom_line(aes(x=time, y=Pa, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="pink")+
  scale_linetype_manual(expression(paste("Awareness (", epsilon, ")")),
                        labels=c("Uniform", "Separated"),
                        values=c("0.5"="dashed", "0.99"="solid"))+
  scale_size_manual(expression(paste("Awareness (", epsilon, ")")),
                    labels=c("Uniform", "Separated"),
                    values=c("0.5"=li_size2, "0.99"=li_size1))+
  ylab("Protective")+
  #scale_x_continuous(limits=c(0, 150), breaks=c(0, 100, 200))+
  ylim(c(0, .5))+
  theme(legend.position="none")+
  ylim(c(0, Pmax))

Ia_plot<-ggplot(t80) + 
  geom_line(aes(x=time, y=Ia, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="maroon")+
  scale_linetype_manual("Awareness",
                        values=c("0.5"="dashed", "0.99"="solid"))+
  scale_size_manual(expression(paste("Awareness (", epsilon, ")")),
                    labels=c("Uniform", "Separated"),
                    values=c("0.5"=li_size2, "0.99"=li_size1))+
  ylab("Infections")+
  theme(legend.position="none")+
  ylim(c(0, Imax))

ts_grid<-plot_grid(Pa_plot, Pb_plot, 
                   Ia_plot, Ib_plot, nrow=2,
                   labels=c("A. Protective \n    \t(group a)", "B. Protective \n    \t(group b)",
                            "C. Infections \n    \t(group a)", "D. Infections \n    \t(group b)"),
                   label_x=.1, label_size=16)

a_vals<-filter(t80, time==80) %>%
  group_by(epsilon) %>%
  summarize(Pa=SPa+IPa+RPa,
            Ia=1-(SPa+SUa)) %>%
  select(Pa, Ia, epsilon) 


b_vals<-filter(t80, time==80) %>%
  group_by(epsilon) %>%
  summarize(Pb=SPb+IPb+RPb,
            Ib=1-(SPb+SUb)) %>%
  select(Pb, Ib, epsilon)

phase<-ggplot(data=filter(df, h==.99))+
  geom_line(aes(y=SPa+IPa+RPa, x=1-(SPa+SUa), linetype=as.factor(epsilon), color="group *a*"), size=li_size1)+
  geom_line(aes(y=SPb+IPb+RPb, linetype=as.factor(epsilon), x=1-(SPb+SUb), size=as.factor(epsilon), color="group *b*"))+
  xlab("Total Infections")+ylab("Protective Attitude")+
  scale_linetype_manual(expression(paste("Awareness (", epsilon, ")")),
                        labels=c("Uniform", "Separated"),
                        values=c("0.5"="dashed", "0.99"="solid"))+
  scale_size_manual(guide="none",
                    expression(paste("Awareness (", epsilon, ")")),
                    labels=c("Uniform", "Separated"),
                    values=c("0.5"=li_size1, "0.99"=li_size1*.2))+
  scale_colour_manual(name="Group", 
                      labels = c("a", "b"),
                      values = c("group *a*"="maroon", "group *b*"="blue"),
                      guide="legend")+
  theme(legend.position="none")+
  geom_point(data=filter(t80, time==80), aes(x=Ia, y=Pa), color="maroon", size=3)+
  geom_point(data=filter(t80, time==80), aes(x=Ib, y=Pb), color="blue", size=3)+
  annotate("segment", 
           xend=a_vals$Ia[1],
           x=a_vals$Ia[1],
           yend=a_vals$Pa[2],
           y=a_vals$Pa[1],
           arrow=arrow(angle=20, length = unit(0.2, "inches")),
           color="grey70")+
  annotate("segment", 
           xend=b_vals$Ib[1],
           x=b_vals$Ib[1],
           yend=b_vals$Pb[2],
           y=b_vals$Pb[1],
           arrow=arrow(angle=20, length = unit(0.2, "inches")),
           color="grey70")+
  annotate("segment", 
           xend=a_vals$Ia[2],
           x=a_vals$Ia[1],
           yend=a_vals$Pa[2],
           y=a_vals$Pa[2],
           arrow=arrow(angle=20, length = unit(0.2, "inches")),
           color="black")+
  annotate("segment", 
           xend=b_vals$Ib[2],
           x=b_vals$Ib[1],
           yend=b_vals$Pb[2],
           y=b_vals$Pb[2],
           arrow=arrow(angle=20, length = unit(0.2, "inches")),
           color="black")+
  geom_point(aes(x=.4, y=.01), size=3)+
  annotate("segment", 
           x=.37,
           xend=.43,
           yend=.04,
           y=.04,
           arrow=arrow(angle=20, length = unit(0.2, "inches")),
           color="grey70")+
  annotate("segment", 
           x=.37,
           xend=.43,
           yend=.07,
           y=.07,
           arrow=arrow(angle=20, length = unit(0.1, "inches")),
           color="black")+
  geom_text(x=.45, y=.01, hjust="left", label="t = 80")+
  geom_text(x=.45, y=.04, hjust="left", label="Change in Protective")+
  geom_text(x=.45, y=.07, hjust="left", label="Change in Infections")+
  geom_text(aes(label="B", x=b_vals$Ib[1], y=mean(b_vals$Pb)), nudge_x=-.015, color="grey70", size=6)+
  geom_text(aes(label="D", x=mean(b_vals$Ib), y=b_vals$Pb[2]), nudge_y=-.020, nudge_x=-.02, color="black", size=6)+
  geom_text(aes(label="A", x=a_vals$Ia[1], y=mean(a_vals$Pa)), nudge_x=.015, color="grey70", size=6)+
  geom_text(aes(label="C", x=mean(a_vals$Ia), y=a_vals$Pa[2]), nudge_y=.020, nudge_x=.01, color="black", size=6)

plot_grid(ts_grid, phase, legend, labels=c("", "E. Phase Portrait", ""), nrow=3, rel_heights=c(8, 8, 1), label_size=16, label_x=.05, label_y=.98)


#### Supplementary Figure 4-6 Function ####
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
    this_plot<-ggplot(df)+geom_line(aes(x=time, y=IUa+IPa), size=.8, color="maroon")+geom_line(aes(x=time, y=IUb+IPb, size=as.factor(epsilon)), color="blue")+
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
                          values = c("group *a*"="maroon", "group *b*"="blue"))+
      theme(strip.text=element_text(size=20))
  }
  
  #Deaths plot
  else{
    
    this_plot<-ggplot(df)+geom_line(aes(x=time, y=Da), size=.8, color="maroon")+geom_line(aes(x=time, y=Db, size=as.factor(epsilon)), color="blue")+
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
                          values = c("group *a*"="maroon", "group *b*"="blue"))+
      theme(strip.text=element_text(size=20))
  }
  
  this_plot
}

#### Supplementary Figure 4 ####
betaI<-fig1_diff_init(beta_a=.21, beta_b=.19,plot_type="I")+ggtitle(expression(paste("Transmission Coefficient (", beta[a], " = 0.21; ", beta[b], " = 0.19)")))
betaD<-fig1_diff_init(beta_a=.21, beta_b=.19,plot_type="D")+ggtitle(" ")
plot_grid(betaI, betaD, nrow=2, labels="AUTO", label_y=.8)


#### Supplementary Figure 5 ####
rhoI<-fig1_diff_init(rho_a=.09, rho_b=.11, plot_type = "I")+ggtitle(expression(paste("Infectious Period (1/", rho[a], " = 11.11; 1/", rho[b], " = 9.09)")))
rhoD<-fig1_diff_init(rho_a=.09, rho_b=.11, plot_type = "D")+ggtitle(" ")
plot_grid(rhoI, rhoD, nrow=2, labels="AUTO", label_y=.8, align="h")


#### Supplementary Figure 6 ####
muI<-fig1_diff_init(mu_a=.015, mu_b=.0067,plot_type="I", line_scale=TRUE)+ggtitle(expression(paste("Infection Fatality Rate (", mu[a], " = 0.015; ", mu[b], " = 0.0067)")))
muD<-fig1_diff_init(mu_a=.015, mu_b=.0067,plot_type="D")+ggtitle(" ")
plot_grid(muI, muD, nrow=2, labels="AUTO", label_y=.8)

                                                                                                                                                                                                                                                                                                                                                                                                                                 ```{r diff-def, eval=TRUE}
                                                                                                                                                                                                                                                                                                                                                                                                                                    
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
                                                                                                                                                                                                                                                                                                                                                                                                                                        this_plot<-ggplot(df)+geom_line(aes(x=time, y=IUa+IPa), size=.8, color="maroon")+geom_line(aes(x=time, y=IUb+IPb, size=as.factor(epsilon)), color="blue")+
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                              values = c("group *a*"="maroon", "group *b*"="blue"))+
                                                                                                                                                                                                                                                                                                                                                                                                                                          theme(strip.text=element_text(size=20))
                                                                                                                                                                                                                                                                                                                                                                                                                                      }
                                                                                                                                                                                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                      #Deaths plot
                                                                                                                                                                                                                                                                                                                                                                                                                                      else{
                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                        this_plot<-ggplot(df)+geom_line(aes(x=time, y=Da), size=.8, color="maroon")+geom_line(aes(x=time, y=Db, size=as.factor(epsilon)), color="blue")+
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                              values = c("group *a*"="maroon", "group *b*"="blue"))+
                                                                                                                                                                                                                                                                                                                                                                                                                                          theme(strip.text=element_text(size=20))
                                                                                                                                                                                                                                                                                                                                                                                                                                      }
                                                                                                                                                                                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                      this_plot
                                                                                                                                                                                                                                                                                                                                                                                                                                    }
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    ```
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    ```{r diff-beta, fig.height=8, fig.cap="**Separated awareness reduces differences in epidemic size between groups in epidemic size that arise from differences in transmission rates coupled with separated mixing.** Plots of [A] infections and [B] deaths over time in group *a* (maroon) and group *b* (blue). We consider different levels of awareness separation [left column: uniform awareness ($\\epsilon = 0.5$); right column: separated awareness ($\\epsilon = 0.99$)] and mixing separation [top row: uniform mixing ($h = 0.5$); bottom row: separated mixing ($h = 0.99$)]. The groups are initialized so that group *a* has a greater transmission coefficient than group *b* ($\\beta_a = 0.21$ and $\\beta_b = 0.19$). We assume the pathogen is introduced in both groups at prevalence $0.0005$. All other parameter values are the same as those used in Figure 1: infectious period ($\\frac{1}{\\rho} = 10$), infection fatality rate ($\\mu = 0.01$), protective measure efficacy ($\\kappa = 0.3$), responsiveness ($\\theta = 100$), memory ($\\ell = 1$), and fatigue ($\\phi = 0$)."}
                                                                                                                                                                                                                                                                                                                                                                                                                                    betaI<-fig1_diff_init(beta_a=.21, beta_b=.19,plot_type="I")+ggtitle(expression(paste("Transmission Coefficient (", beta[a], " = 0.21; ", beta[b], " = 0.19)")))
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    betaD<-fig1_diff_init(beta_a=.21, beta_b=.19,plot_type="D")+ggtitle(" ")
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    plot_grid(betaI, betaD, nrow=2, labels="AUTO", label_y=.8)
                                                                                                                                                                                                                                                                                                                                                                                                                                    ```
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    ```{r diff-rho, fig.height=8, fig.cap="**Separated awareness reduces differences in epidemic size between groups in epidemic size that arise from differences in infectious period coupled with separated mixing.** Plots of [A] infections and [B] deaths over time in group *a* (maroon) and group *b* (blue). We consider different levels of awareness separation [left column: uniform awareness ($\\epsilon = 0.5$); right column: separated awareness ($\\epsilon = 0.99$)] and mixing separation [top row: uniform mixing ($h = 0.5$); bottom row: separated mixing ($h = 0.99$)]. The groups are initialized so that group *a* has a longer infectious period than group *b* ($\\frac{1}{\\rho_a} = 11.11$ and $\\frac{1}{\\rho_b} = 9.09$). We assume the pathogen is introduced in both groups at prevalence $0.0005$. All other parameter values are the same as those used in Figure 1: transmission coefficient ($\\beta$), infection fatality rate ($\\mu = 0.01$), protective measure efficacy ($\\kappa = 0.3$), responsiveness ($\\theta = 100$), memory ($\\ell = 1$), and fatigue ($\\phi = 0$)."}
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    rhoI<-fig1_diff_init(rho_a=.09, rho_b=.11, plot_type = "I")+ggtitle(expression(paste("Infectious Period (1/", rho[a], " = 11.11; 1/", rho[b], " = 9.09)")))
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    rhoD<-fig1_diff_init(rho_a=.09, rho_b=.11, plot_type = "D")+ggtitle(" ")
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    plot_grid(rhoI, rhoD, nrow=2, labels="AUTO", label_y=.8, align="h")
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    ```
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    ```{r diff-mu, fig.height=8, fig.cap="**Separated awareness reduces differences in mortality between groups arising from differences in their infection fatality rates and causes differences in infections between the groups.** Plots of [A] infections and [B] deaths over time in group *a* (maroon) and group *b* (blue). We consider different levels of awareness separation [left column: uniform awareness ($\\epsilon = 0.5$); right column: separated awareness ($\\epsilon = 0.99$)] and mixing separation [top row: uniform mixing ($h = 0.5$); bottom row: separated mixing ($h = 0.99$)]. The groups are initialized so that group *a* has a higher infection fatality rate than group *b* ($\\mu_a = 0.015$ and $\\mu_b = 0.0067$). We assume the pathogen is introduced in both groups at prevalence $0.0005$. All other parameter values are the same as those used in Figure 1: transmission coefficient ($\\beta$),infectious period ($\\frac{1}{\\rho} = 10$), per-capita recovery rate ($\\rho = 0.1$), protective measure efficacy ($\\kappa = 0.3$), responsiveness ($\\theta = 100$), memory ($\\ell = 1$), and fatigue ($\\phi = 0$)."}
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    muI<-fig1_diff_init(mu_a=.015, mu_b=.0067,plot_type="I", line_scale=TRUE)+ggtitle(expression(paste("Infection Fatality Rate (", mu[a], " = 0.015; ", mu[b], " = 0.0067)")))
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    muD<-fig1_diff_init(mu_a=.015, mu_b=.0067,plot_type="D")+ggtitle(" ")
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    plot_grid(muI, muD, nrow=2, labels="AUTO", label_y=.8)
                                                                                                                                                                                                                                                                                                                                                                                                                                    ```
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    \newpage
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                    ## Works Cited
                                                                                                                                                                                                                                                                                                                                                                                                                                    