
vaxdf_process<-function(df){
  
  if("STa" %in% names(df)){ 
    pivot_longer(df, SUa:OBb, 
                        values_to="y",
                        names_sep=c(1, 2),
                        names_to=c("dis", "imm", "group"))
  }
  
  else{
    rename(df, "RR" = "RT") %>%
      pivot_longer(SU:OB, 
                   values_to="y",
                   names_sep=c(1),
                   names_to=c("dis", "imm"))
  }
}



vaxdf_summarise<-function(df, v_start=0){
  
  t_end<-max(df$time)
  
  return(list(
    "total_inf_a" = filter(df, group=="a" & dis=="C" & time==t_end) %>% select(x) %>% sum(),
    "total_inf_b" = filter(df, group=="a" & dis=="C" & time==t_end) %>% select(x) %>% sum(),
    "total_death_a" = filter(df, group=="a" & dis=="C" & time==t_end) %>% select(x) %>% sum(),
    "total_death_b" = filter(df, group=="b" & dis=="D" & time==t_end) %>% select(x) %>% sum(),
    "peak_inf_a" = filter(df, group=="a" & dis=="I") %>% group_by(time) %>% summarise(inf_all=sum(x)) %>% select(inf_all) %>% max(),
    "peak_inf_b" = filter(df, group=="b" & dis== "I") %>% group_by(time) %>% summarise(inf_all=sum(x)) %>% select(inf_all) %>% max(),
    "whenpeak_inf_a" = filter(df, group=="a" & dis=="I") %>% group_by(time) %>% summarise(inf_all=sum(x)) %>% select(inf_all) %>% 
      unlist() %>% which.max() %>% as.numeric(),
    "whenpeak_inf_b" = filter(df, group=="b" & dis=="I") %>% group_by(time) %>% summarise(inf_all=sum(x)) %>% select(inf_all) %>% 
      unlist() %>% which.max() %>% as.numeric(),
    
    
    "total_imm_a" = filter(df, group=="a" & imm %in% c("T", "M") &  dis != "D" & time==t_end) %>% select(x) %>% sum(),
    "total_imm_b" = filter(df, group=="b" & imm %in% c("T", "M") & dis != "D" & time==t_end) %>% select(x) %>% sum(),
    "postvax_deaths_a" = filter(df, group=="a" & imm == "U" & dis=="D" & time > v_start) %>% select(x) %>% range() %>% diff(),
    "postvax_deaths_b" = filter(df, group=="b" & imm == "U" & dis=="D" & time > v_start) %>% select(x) %>% range () %>% diff(),
    "total_first_vax_a" = filter(df, group=="a" & imm == "V" & time==t_end) %>% select(x) %>% sum(),
    "total_first_vax_b" = filter(df, group=="b" & imm == "V" & time==t_end) %>% select(x) %>% sum()
    
  ))
}




plot_single_vax<-function(df, name=name){
  S_plot<-ggplot()+
    geom_line(df, mapping=aes(x=time, y=SU+ST+SM, color=as.factor(var)), size=.8)+
    ylab("Susceptible")+
    scale_color_grey(bquote(.(name)), end=.9)+ylim(c(0,1))+theme(legend.position="bottom")
  
  legend<-get_legend(S_plot)
  
  S_plot<-S_plot+theme(legend.position="none")
  
  I_plot<-ggplot()+
    geom_line(df, mapping=aes(x=time, y=IU+IT+IM, color=as.factor(var)), size=.8)+
    ylab("Infected")+
    scale_color_grey(bquote(.(name)), end=.9)+theme(legend.position="none")
  
  T_plot<-ggplot()+
    geom_line(df, mapping=aes(x=time, y=ST+IT+RT, color=as.factor(var)), size=.8)+
    ylab("Transmissiong-Reducing \nImmunity")+
    scale_color_grey(bquote(.(name)), end=.9)+ylim(c(0,1))+theme(legend.position="none")
  
  M_plot<-ggplot()+
    geom_line(df, mapping=aes(x=time, y=SM+IM+ST+IT+RT, color=as.factor(var)), size=.8)+
    ylab("Mortality-Reducing \nImmunity")+
    scale_color_grey(bquote(.(name)), end=.9)+ylim(c(0,1))+theme(legend.position="none")
  
  
  V_plot<-ggplot()+
    geom_line(df, mapping=aes(x=time, y=OV, color=as.factor(var)), size=.8)+
    ylab("First-Time Vaccines")+
    scale_color_grey(bquote(.(name)), end=.9)+theme(legend.position="none")
  
  B_plot<-ggplot()+
    geom_line(df, mapping=aes(x=time, y=OB, color=as.factor(var)), size=.8)+
    ylab("Total Vaccines")+
    scale_color_grey(bquote(.(name)), end=.9)+theme(legend.position="none")
  
  
  all_plots<-plot_grid(S_plot, I_plot, T_plot, M_plot, V_plot, B_plot, nrow=3)
  
  
  plot_grid(all_plots, legend, nrow=2, rel_heights = c(12,1))
}

plot_single_vax_2var<-function(df, var1_name=var1_name, var2_name=var2_name){
  library(cowplot)
  
  I1_plot<-ggplot()+geom_line(df %>% filter(var1==unique(df$var1)[1]), mapping=aes(x=time, y=IU+IT+IM, color=as.factor(var2)), size=.8)+
    ylab("Infected")+
    scale_color_grey(bquote(.(var2_name)), end=.9)+theme(legend.position="bottom")+
    ggtitle(paste(var1_name, "=", unique(df$var1)[1]))
  
  legend<-get_legend(I1_plot)
  
  
  I1_plot<-I1_plot+theme(legend.position="none")
  
  I2_plot<-ggplot()+geom_line(df %>% filter(var1==unique(df$var1)[2]), mapping=aes(x=time, y=IU+IT+IM, color=as.factor(var2)), size=.8)+
    ylab("Infected")+
    scale_color_grey(end=.9)+theme(legend.position="none")+
    ggtitle(paste(var1_name, "=", unique(df$var1)[2]))
  
  T1_plot<-ggplot()+geom_line(df %>% filter(var1==unique(df$var1)[1]), mapping=aes(x=time, y=ST+IT+RT, color=as.factor(var2)), size=.8)+
    ylab("Transmission-Reducing \nImmunity")+
    scale_color_grey(end=.9)+ylim(c(0,1))+theme(legend.position="none")
  
  T2_plot<-ggplot()+geom_line(df %>% filter(var1==unique(df$var1)[2]), mapping=aes(x=time, y=ST+IT+RT, color=as.factor(var2)), size=.8)+
    ylab("Transmission-Reducing \nImmunity")+
    scale_color_grey(end=.9)+ylim(c(0,1))+theme(legend.position="none")
  
  M1_plot<-ggplot()+geom_line(df %>% filter(var1==unique(df$var1)[1]), mapping=aes(x=time, y=ST+IT+RT+SM+IM, color=as.factor(var2)), size=.8)+
    ylab("Mortality-Reducing \nImmunity")+
    scale_color_grey(end=.9)+ylim(c(0,1))+theme(legend.position="none")
  
  M2_plot<-ggplot()+geom_line(df %>% filter(var1==unique(df$var1)[2]), mapping=aes(x=time, y=ST+IT+RT+SM+IM, color=as.factor(var2)), size=.8)+
    ylab("Mortality-Reducing \nImmunity")+
    scale_color_grey(end=.9)+ylim(c(0,1))+theme(legend.position="none")
  
  all_plots<-plot_grid(I1_plot, I2_plot, T1_plot, T2_plot, M1_plot, M2_plot, nrow=3)
  
  
  plot_grid(all_plots, legend, nrow=2, rel_heights = c(12,1))
}


veff_plot<-function(df, is_B_plot=FALSE, legend=FALSE, eps_val=NA){
  my_epsilons<-unlist(sort(unique(df$epsilon)))
  lab_epsilons<-ifelse(my_epsilons==.5, "Uniform Awareness", 
                       ifelse(my_epsilons==.99,"Separated Awareness",
                              ifelse(my_epsilons>.5 & my_epsilons <.99, "Intermediate Awareness", "Other Awareness")))
  
  if(legend==TRUE){
    legend.pos="bottom"
  }
  else{
    legend.pos="none"
  }
  
  if(!is.na(eps_val)){
    df %<>% filter(epsilon==eps_val)  
  }
  
  #calculate totals
  #df %>% 
  #  group_by(veff, epsilon) %>%
  #  summarize(y=sum(y)*.5) %>%
    #mutate(group="total") %>%
  #  rbind(df) %>%
    
  
  if(is_B_plot){
    df %<>% mutate(size_dummy = ifelse(group=="b" & epsilon==.5, "small", "normal"))
  } else {
    df %<>% mutate(size_dummy = "normal")
  }
ggplot()+geom_line(data=filter(df, group=="a"), 
                   aes(x=veff, y=y, color="group *a*", 
                       linetype=as.factor(epsilon), 
                       size=size_dummy))+
    geom_line(data=filter(df, group=="b"), 
              aes(x=veff, y=y, color="group *b*",
             size=size_dummy,
              linetype=as.factor(epsilon)))+
    ylab("")+
    #xlab(expression(paste("Vaccine Efficacy ", "\n(1-", zeta,"; 1-",kappa, ")")))+
    xlab(expression(paste("Immune Protection")))+
    scale_color_manual("", 
                       labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                       values=c("group *a*"="#f4a896", "group *b*"="#358597"))+
    scale_linetype_manual("", labels=c(expr(paste(!!lab_epsilons[1], " (", epsilon, " = ", !!my_epsilons[1], ")")),
                                       expr(paste(!!lab_epsilons[2], " (", epsilon, " = ", !!my_epsilons[2], ")"))),
                          values=c("dashed", "solid"))+
    scale_size_manual(values=c("normal"=1.1, "small"=.8),
                      guide="none")+
    theme(legend.position=legend.pos,
          legend.box="horizontal", 
          panel.grid.major.x = element_line(size = .2, colour = NA),
          panel.grid.minor.x = element_line(size = .2, colour = NA),
          panel.grid.major.y = element_line(size = .2, colour = NA))+
    scale_x_continuous(trans = trans_reverser('identity'))+
    scale_y_continuous(trans="identity")+
    #scale_x_continuous(trans = trans_reverser('log2'), breaks=c(.5, .125, .031, .008))+
    #scale_y_continuous(trans = 'log2')+
    guides(linetype = guide_legend(override.aes = list(color="black", size=.5)))
        
}

veff_plot_present<-function(df, is_B_plot=FALSE, legend=FALSE, no.asrt=FALSE){
  if(legend==TRUE){
    legend.pos="bottom"
  }
  else{
    legend.pos="none"
  }
  
  if(no.asrt){
    df %<>% filter(epsilon==.5)  
  }
  
  #calculate totals
  #df %>% 
  #  group_by(veff, epsilon) %>%
  #  summarize(y=sum(y)*.5) %>%
  #mutate(group="total") %>%
  #  rbind(df) %>%
  
  
  if(is_B_plot){
    df %<>% mutate(size_dummy = ifelse(group=="b" & epsilon==.5, "small", "normal"))
  } else {
    df %<>% mutate(size_dummy = "normal")
  }
  ggplot()+geom_line(data=filter(df, group=="a"), 
                     aes(x=veff, y=y, color="group *a*", 
                         linetype=as.factor(epsilon), 
                         size=size_dummy))+
    geom_line(data=filter(df, group=="b"), 
              aes(x=veff, y=y, color="group *b*",
                  size=size_dummy,
                  linetype=as.factor(epsilon)))+
    ylab("")+
    #xlab(expression(paste("Vaccine Efficacy ", "\n(1-", zeta,"; 1-",kappa, ")")))+
    xlab(expression(paste("Immune Protection")))+
    scale_color_manual("", 
                       labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                       values=c("group *a*"="#F4A896", "group *b*"="#358597"),
                       guide="none")+
    scale_linetype_manual("", labels=c(expr(paste(!!lab_epsilons[1], " (", epsilon, " = ", !!my_epsilons[1], ")")),
                                       expr(paste(!!lab_epsilons[2], " (", epsilon, " = ", !!my_epsilons[2], ")"))),
                          values=c("dashed", "solid"))+
    scale_size_manual(values=c("normal"=.9, "small"=.7),
                      guide="none")+
    theme(legend.position=legend.pos,
          legend.box="horizontal", 
          panel.grid.major.x = element_line(size = .2, colour = NA),
          panel.grid.minor.x = element_line(size = .2, colour = NA),
          panel.grid.major.y = element_line(size = .2, colour = NA))+
    scale_x_continuous(trans = trans_reverser("identity"), breaks=c(0, .25, .5, .75, 1), labels=c("1.00", "0.75", "0.50", "0.25", "0.00"))+
    scale_y_continuous(trans="identity")+
    #scale_x_continuous(trans = trans_reverser('log2'), breaks=c(.5, .125, .031, .008))+
    #scale_y_continuous(trans = 'log2')+
    guides(linetype = guide_legend(override.aes = list(color="black", size=.5)))
  
}


veff_plot_poster<-function(df, is_B_plot=FALSE, legend=FALSE, no.asrt=FALSE){
  if(legend==TRUE){
    legend.pos="bottom"
  }
  else{
    legend.pos="none"
  }
  
  if(no.asrt){
    df %<>% filter(epsilon==.5)  
  }
  
  #calculate totals
  #df %>% 
  #  group_by(veff, epsilon) %>%
  #  summarize(y=sum(y)*.5) %>%
  #mutate(group="total") %>%
  #  rbind(df) %>%
  
  
  if(is_B_plot){
    df %<>% mutate(size_dummy = ifelse(group=="b" & epsilon==.5, "small", "normal"))
  } else {
    df %<>% mutate(size_dummy = "normal")
  }
  ggplot()+geom_line(data=filter(df, group=="a"), 
                     aes(x=veff, y=y, color="a", 
                         linetype=as.factor(epsilon), 
                         size=size_dummy))+
    geom_line(data=filter(df, group=="b"), 
              aes(x=veff, y=y, color="b",
                  size=size_dummy,
                  linetype=as.factor(epsilon)))+
    ylab("")+
    #xlab(expression(paste("Vaccine Efficacy ", "\n(1-", zeta,"; 1-",kappa, ")")))+
    xlab(expression(paste("Immune Protection")))+
    scale_color_manual("Group",
                       values=c("a"="#f4a896", "b"="#358597"),
                       guide="none")+
    scale_linetype_manual("Awareness",
                          labels=c("Uniform", "Separated"),
                          values=c("0.5"="dashed", "0.99"="solid"))+
    scale_size_manual(values=c("normal"=1.8, "small"=1),
                      guide="none")+
    theme(legend.position=legend.pos,
          legend.box="vertical", 
          panel.grid.major.x = element_line(size = .2, colour = NA),
          panel.grid.minor.x = element_line(size = .2, colour = NA),
          panel.grid.major.y = element_line(size = .2, colour = NA))+
    scale_x_continuous(trans = trans_reverser('identity'))+
    guides(linetype = guide_legend(override.aes = list(color="black", size=.5)))
  
  
}

vdiff_group_plot<-function(df, var, name, legend=FALSE){
  if(legend==TRUE){
    legend.pos="bottom"
  }
  else{
    legend.pos="none"
  }
  
  df$var<-df[[var]]
  
  ggplot(df)+geom_line(aes(y=var, x=theta_mult, color=group, linetype=as.factor(epsilon)))+
    scale_color_manual(" ",
                       labels=c("Group a", "Group b"),
                       values=c("a"="#f4a896", "b"="#358597"))+
    scale_linetype_manual("Awareness",
                          labels=c(expression(paste("Uniform (", epsilon, " = 0.5)")),
                                   expression(paste("Separated (", epsilon, " = 0.99)"))),
                          values=c("0.5"="dashed", "0.99"="solid"))+
    scale_x_continuous(trans="log2")+
    scale_y_continuous(trans="log2")+
    theme(legend.position=legend.pos,
          legend.box="vertical")+
    ylab(name)+
    xlab(expression(paste("Responsiveness Ratio (", theta[a],"/",theta[b], ")")))+
    guides(linetype = guide_legend(override.aes = list(color="black", size=.5)))
}

vdiff_ratio_plot<-function(df, var, name, square=TRUE){

  df$var<-df[[var]]
  
  ggplot(df)+geom_line(aes(y=var, x=theta_mult, linetype=as.factor(epsilon)))+
       scale_linetype_manual("Awareness",
                          labels=c(expression(paste("Uniform (", epsilon, " = 0.5)")),
                                   expression(paste("Separated (", epsilon, " = 0.99)"))),
                          values=c("0.5"="dashed", "0.99"="solid"))+
    scale_x_continuous(trans="log2",
                       limits=c(ifelse(square, 
                                       min(c(df$var, df$theta_mult), na.rm=TRUE),
                                       min(df$theta_mult)),
                                ifelse(square, 
                                       max(c(df$var, df$theta_mult), na.rm=TRUE),
                                       max(df$theta_mult)))
                       )+
    scale_y_continuous(trans="log2",
                       limits=c(ifelse(square, 
                                      min(c(df$var, df$theta_mult), na.rm=TRUE),
                                      min(df$var)),
                                ifelse(square, 
                                     max(c(df$var, df$theta_mult), na.rm=TRUE),
                                     max(df$var)))
                       )+
    theme(legend.position="none")+
    ylab(paste(name, "(a/b)"))+
    xlab(expression(paste("Responsiveness Ratio (", theta[a],"/",theta[b], ")")))+
    geom_hline(yintercept=1, linetype="dotted", color="grey70")#+
    #geom_abline(intercept=0, slope=1, linetype="dotted", color="grey70")
   
}

