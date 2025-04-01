##### e.vancouveri abiotic stress analysis #####
#reading in packaages and data
pacman::p_load(tidyverse, lme4, MuMIn, emmeans, patchwork)
light_temp <- read_csv("alex_anais_tidy.csv") 
turb <- read_csv("tidy_dirtnworms.csv") %>%
  mutate(trial=as.factor(trial), 
         site=as.factor(site))
ev<- light_temp %>%
  rename(tube_diameter_mm = tube_diameter_cm) %>%
  filter(species == "EV") %>%
  mutate(tube_diameter_cm = tube_diameter_mm/10,
         retract_s = replace_na(retract_s, 0), 
         CR_CE_time = emerge_s - retract_s, 
         retract_s = round(retract_s, digits=0),
         trt_comb =str_c(temp_trt, light_trt))

########## Dirt n' Worms data #########
#looking at reemergence time
hist(log(turb$CR_CE_time))
shapiro.test(turb$CR_CE_time) #not normal, using log transformation

#trying a glmm w/ a poisson distribution
mod1 <- glmer(CR_CE_time ~ trt + tube_diameter_cm + (1|trial), family=poisson(link="log"), data=turb)
summary(mod1)
confint(mod1)
#looking at overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(mod1) #grossly overdispersed
#trying a glmm w/ a negative binomial distribution
mnb1 <- glmer.nb(CR_CE_time ~ trt + tube_diameter_cm + (1 | trial), data = turb)
summary(mnb1)
confint.merMod(mnb1, method="Wald")
overdisp_fun(mnb1) #not overdispersed:)

#using the negative binomial and back transforming the estimates:)
backtransformed_turb<- emmeans(mnb1, ~ trt + tube_diameter_cm, type = "response", adjust="bonferroni", level=0.95) %>% as.data.frame()

############# Alex's data ###########
#looking at reemergence time
hist(ev$emerge_s)
shapiro.test(ev$emerge_s) #not normal, using log transformation
mod2 <- glmer(emerge_s ~ temp_trt*light_trt + tube_diameter_cm + (1|trial), family=poisson(link="log"), data=ev)
summary(mod2)
confint(mod2)
overdisp_fun(mod2) #grossly overdispersed
#trying a glmm w/ a negative binomial distribution
mnb2 <- glmer.nb(emerge_s ~ temp_trt*light_trt + tube_diameter_cm + (1 | trial), data = ev)
summary(mnb2)
confint.merMod(mnb2, method="Wald")
overdisp_fun(mnb2) #less overdispersed
backtransformed_ev_rem<- emmeans(mnb2, ~ temp_trt*light_trt + tube_diameter_cm , type = "response", adjust="bonferroni", level=0.95) %>%
  as.data.frame() %>%
  mutate(trt_comb =str_c(temp_trt, light_trt))

#looking at retraction time
hist(ev$retract_s)
shapiro.test(ev$retract_s) #not normal, using log transformation
#trying a glmm w/ a poison distribution
mod3 <- glmer(retract_s ~ temp_trt*light_trt + tube_diameter_cm + (1|trial), family=poisson(link="log"), data=ev)
summary(mod3)
confint(mod3)
overdisp_fun(mod3) #overdispersed
#trying a glmm w/ a negative binomial distribution
mnb3 <- glmer.nb(retract_s ~ temp_trt*light_trt + tube_diameter_cm + (1 | trial), data = ev, control=glmerControl(optimizer="bobyqa"))
summary(mnb3)
overdisp_fun(mnb3) # a little over dispersed
#using the negative binomial and back transforming the estimates:)
backtransformed_ev_ret<- emmeans(mnb3, ~ temp_trt*light_trt + tube_diameter_cm , type = "response", adjust="bonferroni", level=0.95) %>%
  as.data.frame() %>%
  mutate(trt_comb =str_c(temp_trt, light_trt))
confint(backtransformed_ev_ret, adjust = "bonferroni", level = 0.95) #angry??

## this works?? ###

backtransformed_ev_ret1 <- emmeans(mnb3, ~temp_trt*light_trt + tube_diameter_cm, type = "response", adjust="bonferroni", level = 0.95)

backtransformed_ev_ret1_df <- confint(backtransformed_ev_ret1, adjust = "bonferroni", level = 0.95) %>%
  as.data.frame() %>%
  mutate(trt_comb = str_c(temp_trt, light_trt))
####### visualizing data ###############

p1 <- ggplot()+
  geom_point(data=turb, aes(x=trt, y=CR_CE_time, colour=trt), alpha=0.03)+
  #geom_violin(data=turb, aes(x=trt, y=CR_CE_time, colour=trt), alpha=0.03)+
  geom_pointrange(data=backtransformed_turb, aes(x=trt, y=response, ymin=asymp.LCL, ymax=asymp.UCL, colour=trt), size=0.75)+
  scale_colour_manual(values = c("snow3", "snow4", "grey30"), breaks=c("control", "moderate  sediment", "high sediment"))+
  scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
  coord_flip()+
  theme_classic()+
  ylab("Re-emergence time (s)")+
  xlab("")
  
  
 p2 <-  ggplot()+
    geom_point(data=ev,aes(x=temp_trt, y=emerge_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
   #geom_violin(data=ev,aes(x=temp_trt, y=emerge_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
   geom_pointrange(data=backtransformed_ev_rem, aes(x=temp_trt, y=response, ymin=asymp.LCL, 
                                                     ymax=asymp.UCL, colour=trt_comb),
                    position = position_dodge(width = 0.75), size=0.75)+
    scale_colour_manual(values = c("black", "grey30","snow4", "snow3"))+
    #scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
    coord_flip()+
    theme_classic()+
    ylab("Re-emergence time (s)")+
   xlab("")+
   theme(legend.position='none')
 
 p1 / p2

  
 p3 <- ggplot()+
    geom_point(data=ev,aes(x=temp_trt, y=retract_s, colour=trt_comb), position = position_dodge(width = 0.75), alpha=0.15)+
    geom_pointrange(data=backtransformed_ev_ret, aes(x=temp_trt, y=response, ymin=asymp.LCL, 
                                                     ymax=asymp.UCL, colour=trt_comb),
                    position = position_dodge(width = 0.75))+
    scale_colour_manual(values = c("black", "grey30","snow4", "snow3"))+
    #scale_x_discrete(limits=c("control", "moderate  sediment", "high sediment"))+
    coord_flip()+
    theme_classic()+
    ylab("Retraction time (s)")+
    xlab("")


p3/p2/p1
