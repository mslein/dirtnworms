pacman::p_load(tidyverse, lme4, MuMIn)
light_temp <- read_csv("alex_anais_tidy.csv") 
turb <- read_csv("tidy_dirtnworms.csv")
ev<- light_temp %>%
  rename(tube_diameter_mm = tube_diameter_cm) %>%
  filter(species == "EV") %>%
  mutate(tube_diameter_cm = tube_diameter_mm/10,
         retract_s = replace_na(retract_s, 0), 
         CR_CE_time = emerge_s - retract_s, 
         retract_s = round(retract_s, digits=0))


########## Dirt n' Worms data #########

#looking at reemergence time
hist(log(turb$CE_Time_s))
shapiro.test(turb$CE_Time_s) #not normal, using log transformation

mod1 <- glmer(CE_Time_s ~ trt + tube_diameter_cm + (1|trial), family=poisson(link="log"), data=turb)
summary(mod1)
confint(mod1)

#looking at retraction time
hist(turb$CR_Time_s)
shapiro.test(turb$CR_Time_s) #not normal, using log transformation

mod2 <- glmer(CR_Time_s ~ trt + tube_diameter_cm + (1|trial), family=poisson(link="log"), data=turb)
summary(mod2)
confint(mod2)


############# Alex's data ###########
#looking at retraction time
hist(ev$retract_s)
shapiro.test(ev$retract_s) #not normal, using log transformation

mod3 <- glmer(retract_s ~ temp_trt*light_trt + tube_diameter_cm + (1|trial), family=poisson(link="log"), data=ev)
summary(mod3)
confint(mod3)

#looking at reemergence time
hist(ev$emerge_s)
shapiro.test(ev$emerge_s) #not normal, using log transformation
mod4 <- glmer(emerge_s ~ temp_trt*light_trt + tube_diameter_cm + (1|trial), family=poisson(link="log"), data=ev)
summary(mod4)
confint(mod4)



####### visualizing data 
turb %>%
  ggplot(aes(x=trt, y=CR_Time_s))+
  geom_point(alpha=0.5)+
  coord_flip()+
  theme_classic()+
  geom_boxplot()

ev %>%
  ggplot(aes(x=light_trt, y=retract_s, colour=temp_trt))+
  geom_point(position = position_jitterdodge(jitter.width = 0.05, jitter.height = 0.05, dodge.width = 0.8))+
  #geom_jitter(alpha=0.5, width = 0.007)+
  coord_flip()+
  theme_classic()




