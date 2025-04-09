##### e.vancouveri abiotic stress analysis #####
#reading in packaages and data
pacman::p_load(tidyverse, lme4, MuMIn, emmeans, patchwork, glmmTMB)
light_temp <- read_csv("alex_anais_tidy.csv") 
turb <- read_csv("tidy_dirtnworms.csv") %>%
  mutate(trial=as.factor(trial), 
         site=as.factor(site)) %>%
  rename(emerge_s=CR_CE_time) %>%
  filter(emerge_s != 0) #removing 0s...
ev<- light_temp %>%
  rename(tube_diameter_mm = tube_diameter_cm) %>%
  filter(species == "EV") %>%
  mutate(tube_diameter_cm = tube_diameter_mm/10,
         retract_s = replace_na(retract_s, 0), 
         CR_CE_time = emerge_s - retract_s, 
         retract_s = round(retract_s, digits=0),
         trt_comb =str_c(temp_trt, light_trt))

ev_retract <- ev %>% filter(retract_s > 0)
ev_emerge <- ev %>% filter(emerge_s> 0)
########## Dirt n' Worms data #########
#looking at reemergence time
hist(log(turb$emerge_s))
shapiro.test(log(turb$emerge_s)) #not normal, using log transformation
hist(turb$tube_diameter_cm)
shapiro.test(sqrt(turb$tube_diameter_cm)) #this improves it?
#running gamma because it's rate data
mg1 <- glmmTMB(emerge_s ~ trt + sqrt(tube_diameter_cm) + (1 | trial), data = turb, family="Gamma"(link = "log"))
summary(mg1)

#back transforming the estimates:)
backtransformed_turb<- emmeans(mg1, ~ trt, type="response") 
plot_df <-confint(backtransformed_turb, adjust="bonferroni", level=0.95)

#getting the tubediameter estimate
emmeans(mg1, ~ sqrt(tube_diameter_cm),type="response", adjust="bonferroni", level=0.95)


############# Alex's data ###########
#looking at reemergence time
hist(ev$emerge_s)
shapiro.test(ev$emerge_s) #not normal, using log transformation
hist(ev$tube_diameter_cm)
shapiro.test(sqrt(ev$tube_diameter_cm)) #this improves it?

mg2 <- glmmTMB(emerge_s ~ temp_trt*light_trt + sqrt(tube_diameter_cm) + (1 | trial), data = ev_emerge, family="Gamma"(link = "log"))
summary(mg2)

backtransformed_ev_abiotic<- emmeans(mg2, ~ light_trt|temp_trt, type="response", adjust="bonferroni", level=0.95) 
plot_df2 <-confint(backtransformed_ev_abiotic, adjust="bonferroni", level=0.95) %>%
  as.data.frame() %>%
  mutate(trt_comb =str_c(temp_trt, light_trt))

#need to double check whether this is transforms the sqrt() too
emmeans(mg2, ~sqrt(tube_diameter_cm), type="response", adjust="bonferroni", level=0.95)

#looking at retraction time
hist(ev$retract_s)
shapiro.test(ev$retract_s) #not normal, using log transformation

#need to think about this more, exluding these zeros makes this a very different sample size
mg3 <- glmmTMB(retract_s ~ temp_trt*light_trt + sqrt(tube_diameter_cm) + (1 | trial), data = ev_retract, family="Gamma"(link = "log"))
summary(mg3)

backtransformed_ev_ret<- emmeans(mg3, ~ temp_trt*light_trt, type = "response", adjust="bonferroni", level=0.95) 

plot_df3 <-confint(backtransformed_ev_ret, adjust="bonferroni", level=0.95) %>%
  as.data.frame() %>%
  mutate(trt_comb =str_c(temp_trt, light_trt))



emmeans(mg3, ~sqrt(tube_diameter_cm), type="response", adjust="bonferroni", level=0.95)


####### plots #######
p1 <- ggplot()+
  geom_point(data=turb, aes(x=trt, y=emerge_s, colour=trt), alpha=0.03)+
  #geom_violin(data=turb, aes(x=trt, y=CR_CE_time, colour=trt), alpha=0.03)+
  geom_pointrange(data=plot_df, aes(x=trt, y=response, ymin=asymp.LCL, ymax=asymp.UCL, colour=trt), size=0.75)+
  scale_colour_manual(values = c("black", "navajowhite3", "lightsalmon4"), breaks=c("control", "moderate  sediment", "high sediment"))+
  scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
  coord_flip()+
  theme_classic()+
  ylab("Re-emergence time (s)")+
  xlab("")
  
  
 p2 <-  ggplot()+
    geom_point(data=ev_emerge,aes(x=temp_trt, y=emerge_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
   #geom_violin(data=ev,aes(x=temp_trt, y=emerge_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
   geom_pointrange(data=plot_df2, aes(x=temp_trt, y=response, ymin=asymp.LCL, 
                                                     ymax=asymp.UCL, colour=trt_comb),
                    position = position_dodge(width = 0.75), size=0.75)+
    scale_colour_manual(values = c("black", "darkgrey","red4", "tomato3"))+
    #scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
    coord_flip()+
    theme_classic()+
    ylab("Re-emergence time (s)")+
   xlab("")+
   theme(legend.position='none')

 p3 <- ggplot()+
    geom_point(data=ev_retract,aes(x=temp_trt, y=retract_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
    geom_pointrange(data=plot_df3, aes(x=temp_trt, y=response, ymin=asymp.LCL, 
                                                     ymax=asymp.UCL, colour=trt_comb),
                    position = position_dodge(width = 0.75), size=0.75)+
    scale_colour_manual(values=c("black", "darkgrey","red4", "tomato3"))+
    #scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
    coord_flip()+
    theme_classic()+
    ylab("Retraction time (s)")+
    xlab("")


fig <- (p3 / p2 / p1)
ggsave(filename="figures/fig3.png", fig, dpi=800, width=7, height=5)

hist(sqrt(turb$tube_diameter_cm))
shapiro.test(sqrt(turb$tube_diameter_cm))


ggplot()+
  geom_point(data=ev, aes(x=sqrt(tube_diameter_cm), y=emerge_s), colour="grey", alpha=0.5)+
  geom_point(data=turb,aes(x=sqrt(tube_diameter_cm), y=emerge_s), colour="black", alpha=0.5)+
  geom_abline(slope=5.624, intercept=84.5, colour="grey")+
  geom_abline(slope=5.684, intercept=233, colour="black")+
  #scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
  #coord_flip()+
  theme_classic()+
  ylab("Retraction time (s)")+
  xlab("sqrt(tube diameter)")
