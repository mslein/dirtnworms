##### e.vancouveri abiotic stress analysis #####
#reading in packaages and data
pacman::p_load(tidyverse, lme4, MuMIn, emmeans, patchwork, glmmTMB, car, ggsignif)
light_temp <- read_csv("alex_anais_tidy.csv") 
turb <- read_csv("tidy_dirtnworms.csv") %>%
  mutate(trial=as.factor(trial), 
         site=as.factor(site)) %>%
  rename(emerge_s=CR_CE_time) 
ev<- light_temp %>%
  rename(tube_diameter_mm = tube_diameter_cm) %>%
  filter(species == "EV") %>%
  mutate(tube_diameter_cm = tube_diameter_mm/10,
         retract_s = replace_na(retract_s, 0),
         retract_s = round(retract_s, digits=0),
         adj_retract = retract_s+1,
         trt_comb =str_c(temp_trt, light_trt))

########## Dirt n' Worms data #########
#looking at reemergence time
hist(log(turb$emerge_s))
shapiro.test(log(turb$emerge_s)) #not normal, using log transformation
hist(sqrt(turb$tube_diameter_cm))
shapiro.test(sqrt(turb$tube_diameter_cm)) #this improves it?
#running gamma because it's rate data, lots of zeros so did zero inflated gamma
mg1 <- glmmTMB(emerge_s ~ trt + sqrt(tube_diameter_cm) + (1 | trial), ziformula = ~1,
               data = turb, family = ziGamma(link = "log"))
Anova(mg1, type=2)
df.residual(mg1) #calculating for F stats
#back transforming by hand because of the asymmetric CIs
emm_1 <- emmeans(mg1, specs = ~ trt, type = "link")
df1 <-summary(emm_1, type = "response") %>% as.data.frame()  # Back-transform manually if needed
pairs(emm_1, adjust = "bonferroni") # to get p-values between the different levels

#getting the back transformed tube diameter estimate
emmeans(mg1, ~sqrt(tube_diameter_cm), type="response", level=0.95)

emm_1b <- emmeans(mg1, specs = ~ tube_diameter_cm, type = "link")
df1b <-summary(emm_1b, type = "response")


############# Alex's data ###########
#looking at reemergence time
hist(ev$emerge_s)
shapiro.test(ev$emerge_s) #not normal, using log transformation
hist(sqrt(ev$tube_diameter_cm))
shapiro.test(sqrt(ev$tube_diameter_cm)) #this improves it?
#running gamma because it's rate data, lots of zeros so did zero inflated gamma
mg2 <- glmmTMB(emerge_s ~ temp_trt*light_trt + sqrt(tube_diameter_cm) + (1 | trial), ziformula = ~1,
               data = ev, family = ziGamma(link = "log"))
Anova(mg2, type=3)
df.residual(mg2)
#back transforming by hand because of the asymmetric CIs
emm_2 <- emmeans(mg2, specs = ~ temp_trt*light_trt, type = "link")
df2 <-summary(emm_2, type = "response") %>% as.data.frame() %>%
  mutate(trt_comb =str_c(temp_trt, light_trt))
pairs(emm_2, adjust = "bonferroni") # to get p-values between the different levels
#need to double check whether this is transforms the sqrt() too
emmeans(mg2, ~sqrt(tube_diameter_cm), type="response", level=0.95)

emm_2b <- emmeans(mg2, specs = ~ tube_diameter_cm, type = "link")
df2b <-summary(emm_2b, type = "response")



#### retraction time (SI)
#looking at retraction time
hist(ev$retract_s)
shapiro.test(ev$retract_s) #not normal, using log transformation

#back transforming by hand because of the asymmetric CIs
mg3 <- glmmTMB(retract_s ~ temp_trt*light_trt + sqrt(tube_diameter_cm) + (1 | trial), ziformula = ~1,
               data = ev, family = ziGamma(link = "log"))
Anova(mg3, type=3)
#back transforming by hand because of the asymmetric CIs
emm_3 <- emmeans(mg3, specs = ~ temp_trt*light_trt, type = "link")
pairs(emm_3) # to get p-values between the different levels
df3 <-summary(emm_3, type = "response") %>% as.data.frame() %>%
  mutate(trt_comb =str_c(temp_trt, light_trt))
#need to double check whether this is transforms the sqrt() too
emmeans(mg3, ~tube_diameter_cm, type="response", level=0.95)

emm_3b <- emmeans(mg3, specs = ~ tube_diameter_cm, type = "link")
df3b <-summary(emm_3b, type = "response")


ggplot(df3b, aes(x = tube_diameter_cm, y = response)) +
  geom_line() +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.2) +
  labs(x = "Tube Diameter (cm)", y = "Predicted Response") +
  theme_minimal()


####### plots #######
p1 <- ggplot(data=turb, aes(x=trt, y=emerge_s, colour=trt), )+
  geom_point(alpha=0.03)+
  #geom_violin(data=turb, aes(x=trt, y=CR_CE_time, colour=trt), alpha=0.03)+
  geom_pointrange(data=df1, aes(x=trt, y=response, ymin=asymp.LCL, ymax=asymp.UCL, colour=trt), size=0.75)+
  scale_colour_manual(values = c("black", "navajowhite3", "lightsalmon4"), breaks=c("control", "moderate  sediment", "high sediment"))+
  scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"), 
                   labels=c("low sediment", "moderate sediment", "high sediment"))+
  coord_flip()+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("Re-emergence time (s)")+
  xlab("")+
  geom_segment(aes(x = 1, xend = 3, y = 650, yend = 650), colour="black") +   # A vs B bracket
  geom_text(aes(x = 2, y = 670, label = "*"), colour="black")
  
  
 p2 <-  ggplot(data=ev,aes(x=trt_comb, y=emerge_s, colour=trt_comb), position = position_dodge(width = 0.75))+
    geom_point(alpha=0.15)+
   #geom_violin(data=ev,aes(x=temp_trt, y=emerge_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
   geom_pointrange(data=df2, aes(x=trt_comb, y=response, ymin=asymp.LCL, 
                                                     ymax=asymp.UCL, colour=trt_comb),
                    position = position_dodge(width = 0.75), size=0.75)+
    scale_colour_manual(values = c("black", "darkgrey","red4", "tomato3"))+
    scale_x_discrete(labels=c("cold/dark", "cold/light", "warm/dark", "warm/light"))+
    coord_flip()+
    theme_classic()+
    ylab("Re-emergence time (s)")+
   xlab("")+
   theme(legend.position='none')+
   geom_segment(aes(x = 1, xend = 3, y = 600, yend = 600), colour="black") +   # A vs B bracket
   geom_text(aes(x = 2, y = 620, label = "*"), colour="black") +          # A vs B label
   geom_segment(aes(x = 2, xend = 4, y = 550, yend = 550), colour="black") +   # A vs C bracket
   geom_text(aes(x = 3, y = 570, label = "*"), colour="black") 




fig <- ( p2 / p1) + plot_annotation(tag_levels = "a")
ggsave(filename="figures/fig3.png", fig, dpi=800, width=7, height=5)

hist(sqrt(turb$tube_diameter_cm))
shapiro.test(sqrt(turb$tube_diameter_cm))



#SI figures

ps2 <- ggplot()+
  geom_point(data=ev,aes(x=trt_comb, y=retract_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
  geom_pointrange(data=df3, aes(x=trt_comb, y=response, ymin=asymp.LCL, 
                                ymax=asymp.UCL, colour=trt_comb),
                  position = position_dodge(width = 0.75), size=0.75)+
  scale_colour_manual(values=c("black", "darkgrey","red4", "tomato3"))+
  scale_x_discrete(labels=c("cold/dark", "cold/light", "warm/dark", "warm/light"))+
  coord_flip()+
  theme_classic()+
  theme(legend.position="none")+
  ylab("Retraction time (s)")+
  xlab("")

ggsave(filename="figures/sfig2.png", ps2, dpi=800, width=5, height=3)



#figure for tube diameter
ps3 <- ggplot()+
  geom_point(data=ev, aes(x=sqrt(tube_diameter_cm), y=emerge_s), colour="grey", alpha=0.5)+
  geom_point(data=turb,aes(x=sqrt(tube_diameter_cm), y=emerge_s), colour="black", alpha=0.5)+
  geom_abline(slope=114, intercept=0.546, colour="grey")+
  geom_abline(slope=233, intercept=0.684, colour="black")+
  #scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
  #coord_flip()+
  theme_classic()+
  ylab("Re-emergence time (s)")+
  xlab("sqrt(tube diameter)")

ggsave(filename="figures/sfig3.png", ps3, dpi=800, width=5, height=3)


ps4 <- ggplot()+
  geom_point(data=ev, aes(x=sqrt(tube_diameter_cm), y=retract_s), colour="grey", alpha=0.5)+
  geom_abline(slope=2.63, intercept=0.546, colour="grey")+
  #scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
  #coord_flip()+
  theme_classic()+
  ylab("Retraction time (s)")+
  xlab("sqrt(tube diameter)")

ggsave(filename="figures/sfig4.png", ps4, dpi=800, width=5, height=3)
