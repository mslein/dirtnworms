pacman::p_load(tidyverse, lme4)
light_temp <- read_csv("alex_anais_tidy.csv") 
turb <- read_csv("tidy_dirtnworms.csv")
ev<- light_temp %>%
  rename(tube_diameter_mm = tube_diameter_cm) %>%
  filter(species == "EV") %>%
  mutate(tube_diameter_cm = tube_diameter_mm/10,
         retract_s = replace_na(retract_s, 0), 
         CR_CE_time = emerge_s - retract_s)



#i think this is the right random effects structure?
mod1 <- lmer(CR_CE_time ~ trt + tube_diameter_cm +(1|site/trial/indiv), data=turb)
summary(mod1)
confint(mod1)


#also singular
moda <- lmer(CR_Time_s ~ trt + tube_diameter_cm +(1|site/indiv), data=turb) #trying to match Alex's structure is futile


mod2 <- lmer(emerge_s ~ light_trt*temp_trt + tube_diameter_cm + (1|indiv), data=ev) #this is always singular
mod3 <- lmer(retract_s ~  light_trt*temp_trt + tube_diameter_cm + (1|indiv), data=ev) #idk about this random effect
mod4 <- lmer(CR_CE_time ~ light_trt*temp_trt + tube_diameter_cm + (1|indiv), data=ev) #this is also always singular

