#### Figure 1 #### 

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


line_sz<-.8

png(file="../../Figs/fig1.png",width=7, height=6, res=300, units="in")

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
  scale_size_manual(values=c("0.5"=.3*line_sz, "0.99"=line_sz), guide="none", na.value = line_sz)+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(strip.text=element_text(size=20))

dev.off()

#### Figure 2 ####
theme_set(
  theme_classic(base_size = 12)
)


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
  geom_line(aes(x=time, y=Pb, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="#358597")+
  scale_linetype_manual("Awareness",
                        labels=c(expression(paste("Uniform (",epsilon, "=0.5)")),
                                 expression(paste("Separated (", epsilon, "=0.99)"  ))),
                        values=c("0.5"="dashed", "0.99"="solid"))+
  scale_size_manual(expression(paste("Awareness (", epsilon, ")")),
                    labels=c("Uniform", "Separated"),
                    values=c("0.5"=li_size2, "0.99"=li_size1),
                    guide="none")+
  ylab("Protective")+
  ylim(c(0, Pmax))+
  theme(legend.position="none")

Ib_plot<-ggplot(t80) + 
  geom_line(aes(x=time, y=Ib, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="#358597")+
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
  geom_line(aes(x=time, y=Pa, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="#f4a896")+
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
  geom_line(aes(x=time, y=Ia, linetype=as.factor(epsilon), size=as.factor(epsilon)), color="#f4a896")+
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
                   labels=c("A", "B",
                            "C", "D"),
                   label_size=12, hjust=0)

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
  scale_linetype_manual(name="",
                        labels=c(expression(paste("Uniform Awareness (", epsilon, " = 0.5)")), 
                                 expression(paste("Separated Awareness (", epsilon, " = 0.99)"))),
                        values=c("0.5"="dashed", "0.99"="solid"))+
  scale_size_manual(name="",
                    labels=c(expression(paste("Uniform Awareness (", epsilon, " = 0.5)")), 
                             expression(paste("Separated Awareness (", epsilon, " = 0.99)"))),
                    values=c("0.5"=li_size1, "0.99"=li_size1*.2))+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"),
                      guide="legend")+
  guides(linetype = guide_legend(nrow=2, override.aes = list(color="black", size=.5)),
         color = guide_legend(nrow = 2),
         size="none")+
  geom_point(data=filter(t80, time==80), aes(x=Ia, y=Pa), color="#f4a896", size=3)+
  geom_point(data=filter(t80, time==80), aes(x=Ib, y=Pb), color="#358597", size=3)+
  annotate("segment", 
           xend=a_vals$Ia[1],
           x=a_vals$Ia[1],
           yend=a_vals$Pa[2],
           y=a_vals$Pa[1],
           arrow=arrow(angle=20, length = unit(0.1, "inches")),
           color="grey70")+
  annotate("segment", 
           xend=b_vals$Ib[1],
           x=b_vals$Ib[1],
           yend=b_vals$Pb[2],
           y=b_vals$Pb[1],
           arrow=arrow(angle=20, length = unit(0.1, "inches")),
           color="grey70")+
  annotate("segment", 
           xend=a_vals$Ia[2],
           x=a_vals$Ia[1],
           yend=a_vals$Pa[2],
           y=a_vals$Pa[2],
           arrow=arrow(angle=20, length = unit(0.1, "inches")),
           color="black")+
  annotate("segment", 
           xend=b_vals$Ib[2],
           x=b_vals$Ib[1],
           yend=b_vals$Pb[2],
           y=b_vals$Pb[2],
           arrow=arrow(angle=20, length = unit(0.1, "inches")),
           color="black")+
  geom_point(aes(x=.4, y=.01), size=3)+
  annotate("segment", 
           x=.37,
           xend=.41,
           yend=.04,
           y=.04,
           arrow=arrow(angle=20, length = unit(0.1, "inches")),
           color="grey70")+
  annotate("segment", 
           x=.37,
           xend=.41,
           yend=.07,
           y=.07,
           arrow=arrow(angle=20, length = unit(0.1, "inches")),
           color="black")+
  annotate("text", x=.42, y=.01, hjust="left", label="t = 80")+
  annotate("text", x=.42, y=.04, hjust="left", label="Change in Protective")+
  annotate("text", x=.42, y=.07, hjust="left", label="Change in Infections")+
  annotate("text", label="B", x=b_vals$Ib[1]-.015, y=mean(b_vals$Pb), color="grey70", size=4)+
  annotate("text", label="D", x=mean(b_vals$Ib)-.02, y=b_vals$Pb[2]-.02, color="black", size=4)+
  annotate("text", label="A", x=a_vals$Ia[1]+.015, y=mean(a_vals$Pa), color="grey70", size=4)+
  annotate("text", label="C", x=mean(a_vals$Ia)+.01, y=a_vals$Pa[2]+.02, color="black", size=4)+
  theme(legend.position="bottom")


phase_leg<-get_legend(phase)

phase<-phase+theme(legend.position="none")

pdf(file="../../Figs/aware-fig2.pdf", height=7, width=5)

plot_grid(ts_grid, phase, phase_leg, labels=c("", "E", ""), nrow=3, 
          rel_heights=c(4,4,1), label_size=12, hjust=0)

dev.off()

#### Figure 3 ####
source("aware-eqs.R")
library(cowplot)

theme_set(
  theme_classic(base_size = 16)
)

#can change the values of epsilon here. Vectors must have length 2.
epsilon_vals<-c(.5, .99)
epsilon_labs<-ifelse(epsilon_vals==.5, "Uniform Awareness", 
                     ifelse(epsilon_vals==.99,"Separated Awareness",
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

pdf(file="../../Figs/aware-fig3.pdf",width=6, height=5.5)

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

#### Figure 4 ####
library(cowplot)
source("split-vax-eqs.R")
source("aware-vax-process.R")

source("aware-eqs.R")
theme_set(
  theme_classic(base_size = 16)
)

#can change vaccination start time here
v_start_val<-200

#can change the values of epsilon here. Vectors must have length 2.
epsilon_vals<-c(.5, .99)

vax_df<-lapply(epsilon_vals, function(x) full_aware_vax(theta=40, beta=.2, rho=.1, mu_a=.02, mu_b=.01, zeta=.05, kappa=.05, h=.99, time=2000+v_start_val, v_start=v_start_val, I0_a=.001/2, I0_b=0.001/2, phi=0.01, omega=0.01, ell=30, epsilon=x) %>% 
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
  mutate(eps_lab=ifelse(epsilon==.5, "Uniform", 
                        ifelse(epsilon==.99,"Separated",
                               ifelse(epsilon>.5 & epsilon <.99, "Intermediate", "")))) %>%
  mutate(t_postvax=time-v_start_val)


Ivax_dyn<-ggplot()+
  geom_line(data=vax_df %>% filter(group=="a" & epi== "I" & t_postvax>=0), aes(x=t_postvax, y=y, color="group *a*"), size=.8)+
  geom_line(data=vax_df %>% filter(group=="b" & epi== "I" & t_postvax>=0), aes(x=t_postvax, y=y, color="group *b*"), size=.8)+
  ylab("Infections")+
  xlab("Time (post-vaccine)")+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  # theme(legend.position="bottom")+
  #    annotate("segment", x = 200, xend = 200, y =.1, yend = .01,
  #          arrow = arrow(type="closed", length=unit(0.30,"cm"), angle=20),
  #         size=.8)+
  # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  # annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  # annotate("text", x=200, y=.11, label="Vaccination", hjust=0)+
  facet_grid(~epsilon+eps_lab, 
             labeller=label_bquote(
               cols=.(eps_lab)~"Awareness ("*epsilon==.(epsilon)*")",
               rows=NULL))+theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+#+ylim(c(0, .2))
  geom_hline(yintercept=seq(from=.005, to=.025, by=.005), linetype="dotted", color="gray")

Dvax_dyn<-ggplot()+
  geom_line(data=vax_df %>% filter(group=="a" & epi== "D" & t_postvax>=0), aes(x=t_postvax, y=y, color="group *a*"), size=.8)+
  geom_line(data=vax_df %>% filter(group=="b" & epi== "D" & t_postvax>=0), aes(x=t_postvax, y=y, color="group *b*"), size=.8)+
  ylab("Deaths")+
  xlab("Time (post-vaccine)")+
  #scale_y_log10(limits=c(4e-7,1e-3))+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(legend.position="bottom",
        plot.title=element_text(size=16))+
  # annotate("segment", x = 200, xend = 200, y=9e-5, yend = 5e-6,
  #          arrow = arrow(type="closed", length=unit(0.30,"cm"), angle=20),
  #         size=.8)+
  # annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, size=1)+
  # annotate("segment", x=-Inf, xend=-Inf, y=0, yend=Inf, size=1)+
  # annotate("text", x=200, y=2e-4, label="Vaccination", hjust=0)+
  facet_grid(~epsilon+eps_lab)+
  theme(strip.text.y = element_blank(),
        strip.text.x = element_blank())+
  geom_hline(yintercept=c(1e-6, 2e-6,3e-6, 4e-6, 5e-6, 6e-6), linetype="dotted", color="gray")

pdf(file="../../Figs/aware-fig4.pdf",width=6.7, height=6.7)

plot_grid(Ivax_dyn+theme(legend.position="none"), Dvax_dyn, get_legend(Ivax_dyn), nrow=3, align="hv", rel_heights = c(6,8,1))

dev.off()

#### Figure 5 ####

source("split-vax-eqs.R")
source("aware-vax-process.R")

#can change vaccination start time here
v_start_val<-200

#can change the values of epsilon here. Vectors must have length 2.
epsilon_vals<-c(.5, .99)

options(dplyr.summarise.inform = FALSE)

full_df<-lapply(epsilon_vals, function(epsilon) lapply(c(1/(2^c((0:14)*.5)),.6, .7, .8), function(x) full_aware_vax(theta=40, beta=.2, rho=.1, mu_a=.02, mu_b=.01, zeta=x, kappa=x, h=.99, time=2000+v_start_val, v_start=v_start_val, I0_a=.001/2, I0_b=0.001/2, phi=0.01, omega=0.01, ell=30, epsilon=epsilon) %>% 
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
  theme_classic(base_size = 12)
)

inf_tot<-filter(dis_df,  time %in% c(v_start_val,v_start_val+2000) & dis=="C") %>%
  group_by(group, veff, epsilon) %>%
  summarize(y=diff(range(y))) 

I_plot<-veff_plot(inf_tot, legend=TRUE)#+ylim(c(.4, 3.3))
legend_epi<-get_legend(I_plot)

I_plot<-I_plot+theme(legend.position="none")#+ylim(c(1, 3.6))

death_tot<-filter(dis_df, time %in% c(v_start_val,v_start_val+2000) & dis=="D") %>%
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

png(file="../../Figs/fig5.png",width=4000, height=2250, res=500)

plot_grid(ADplot, legend_epi, nrow=2, rel_heights=c(7,2))

dev.off()


