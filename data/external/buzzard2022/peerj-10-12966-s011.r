# Figure 4. assessing the environmental drivers of H2 uptake.

library(readr)
library(tidyverse)
library(gvlma)

data<- read_csv("Buzzard-RWH-H2-insitu-Fluxes20211119.csv")
soil.h2<-as.data.frame(data)

### h2 All by soil moisture
#control
soil.h2.ctrl<- soil.h2 %>%
  filter(Treatment == "Control")
lm.ctrl<- lm(soil.h2.ctrl$Hydrogen.flux ~ soil.h2.ctrl$sm)
summary(lm.ctrl) # significant p-value: 0.001247 Adjusted R-squared: 0.5037 slope: -26.351 


soil.h2.bpc<- soil.h2 %>%
  filter(Treatment == "Passive")
lm.bpc<- lm(soil.h2.bpc$Hydrogen.flux ~ soil.h2.bpc$sm)
summary(lm.bpc) # significant p-value: 0.0007153 Adjusted R-squared: 0.5399 slope: -34.306 

soil.h2.bai<- soil.h2 %>%
  filter(Treatment == "Active")
lm.bai<- lm(soil.h2.bai$Hydrogen.flux ~ soil.h2.bai$sm)
summary(lm.bai) # not significant p-value: 0.09438 Adjusted R-squared: 0.1289 slope: -10.5853

soil.h2.gwi<- soil.h2 %>%
  filter(Treatment == "Greywater")
lm.gwi<- lm(soil.h2.gwi$Hydrogen.flux ~ soil.h2.gwi$sm)
summary(lm.gwi) # significant p-value: 0.02404 Adjusted R-squared: 0.2647 slope: -29.530


# Reorder data along treatment gradient new data frame
treatment_ordered.h2 <- soil.h2 %>%
  mutate(Treatment = fct_relevel(Treatment, 
                                 "Control", "Passive", "Active", 
                                 "Greywater"))

## generalized multiple mixed linear model (soil hydrogen ~ sm and treatment)
library(rsq)
glm <- glm(treatment_ordered.h2$Hydrogen.flux ~ treatment_ordered.h2$sm +treatment_ordered.h2$Treatment)
summary(glm) # significant p-value h2 sm < 0.0001 
rsq(glm) # R-squared:  0.3971
rsq(glm, adj = T) # R-squared:  0.3563

## generalized linear model (soil hydrogen ~ sm)
glm.sm <- glm(treatment_ordered.h2$Hydrogen.flux ~ treatment_ordered.h2$sm)
summary(glm.sm) # significant p-value h2 sm < 0.001 
rsq(glm.sm) # r-squared = 0.1944
rsq(glm.sm, adj = T) # R-squared:  0.1814

## generalized multiple mixed linear model (soil hydrogen ~ temp and treatment)
glm.t <- glm(treatment_ordered.h2$Hydrogen.flux ~ treatment_ordered.h2$temp +treatment_ordered.h2$Treatment)
summary(glm.t) # significant p-value < 0.0001 
rsq(glm.t) # R-squared:  0.2844 
rsq(glm.t, adj = T) # r-squared = 0.2359

## generalized linear model (soil hydrogen ~ temp)
glm.temp <- glm(treatment_ordered.h2$Hydrogen.flux ~ treatment_ordered.h2$temp)
summary(glm.temp) # significant p-value < 0.0001 
rsq(glm.temp) # R-squared:  0.2453 
rsq(glm.temp, adj = T) # r-squared = 0.2331

######## modeled using the new h2 fluxes that are calculated using the sampling event response factor

## generalized multiple mixed linear model (soil hydrogen ~ sm and treatment)
library(rsq)
glm <- glm(treatment_ordered.h2$Hydrogen.Uptake ~ treatment_ordered.h2$sm +treatment_ordered.h2$Treatment)
summary(glm) # significant p-value h2 sm < 0.001 
rsq(glm) # R-squared:  0.3346395
rsq(glm, adj = T) # R-squared:  0.2895303

## generalized linear model (soil hydrogen ~ sm)
glm.sm <- glm(treatment_ordered.h2$Hydrogen.Uptake ~ treatment_ordered.h2$sm)
summary(glm.sm) # significant p-value h2 sm < 0.01 
rsq(glm.sm) # r-squared = 0.1547107
rsq(glm.sm, adj = T) # R-squared:  0.141077

## generalized multiple mixed linear model (soil hydrogen ~ temp and treatment)
glm.t <- glm(treatment_ordered.h2$Hydrogen.Uptake ~ treatment_ordered.h2$temp +treatment_ordered.h2$Treatment)
summary(glm.t) # significant p-value < 0.0001 
rsq(glm.t) # R-squared:  0.1885558
rsq(glm.t, adj = T) # r-squared = 0.1335426

## generalized linear model (soil hydrogen ~ temp)
glm.temp <- glm(treatment_ordered.h2$Hydrogen.Uptake ~ treatment_ordered.h2$temp)
summary(glm.temp) # significant p-value < 0.01 
rsq(glm.temp) # R-squared:  0.1189634
rsq(glm.temp, adj = T) # r-squared = 0.1047532



library(ggpubr)
library(cowplot)
library(ggplot2)
library(wesanderson)


coef(glm(soil.h2$Hydrogen.flux~soil.h2$sm+soil.h2$Treatment))
coef(glm(soil.h2$Hydrogen.flux~soil.h2$sm))

#with boxplots
sm.hydro<-ggscatterhist(
  treatment_ordered.h2, x = "sm", y = "Hydrogen.flux", 
  color = "Treatment", size = 3,  shape="Season",
  palette = rev(wes_palette("Moonrise2")),
  margin.plot = "boxplot", xlab =expression("Soil Moisture"*" (m"^"3"*"m"^"-3"*")"),
  ggtheme = theme_bw(), 
  legend = "right", 
  ylab =expression("Hydrogen Fluxes (nmol H "['2']*" m "^"-2"*" s "^"-1"*")")
)

sm.hydro$sp <- sm.hydro$sp + 
  geom_abline(intercept = -0.6967, slope = -12.7466)
sm.hydro



#with boxplots for new response factor calculated for each sampling event
sm.hydro<-ggscatterhist(
  treatment_ordered.h2, x = "sm", y = "Hydrogen.Uptake", 
  color = "Treatment", size = 3,  shape="Season",
  palette = rev(wes_palette("Moonrise2")),
  margin.plot = "boxplot", xlab =expression("Soil Moisture"*" (m"^"3"*"m"^"-3"*")"),
  ggtheme = theme_bw(), 
  legend = "right", 
  ylab =expression("Hydrogen Fluxes (nmol H "['2']*" m "^"-2"*" s "^"-1"*")")
)

coef(glm(soil.h2$Hydrogen.Uptake~soil.h2$sm+soil.h2$Treatment))
coef(glm(soil.h2$Hydrogen.Uptake~soil.h2$sm))

sm.hydro$sp <- sm.hydro$sp + 
  geom_abline(intercept = -1.201065, slope = -9.481161)
sm.hydro

ggplot2::ggsave("Fig4-H2-sm-p-single-RF20211129.tiff",width = 6, height = 4, dpi = 600, units = "in")

coef(glm(soil.h2$Hydrogen.flux~soil.h2$temp+soil.h2$Treatment))
coef(glm(soil.h2$Hydrogen.flux~soil.h2$temp))

st.hydro<-ggscatterhist(
  treatment_ordered.h2, x = "temp", y = "Hydrogen.flux", 
  color = "Treatment", size = 3,  shape="Season",
  palette = rev(wes_palette("Moonrise2")),
  margin.plot = "boxplot", xlab="Soil Temperature °C",
  ggtheme = theme_bw(), legend = "right", 
  ylab =expression("Hydrogen Fluxes (nmol H "['2']*" m "^"-2"*" s "^"-1"*")")
)
st.hydro$sp <- st.hydro$sp + 
  geom_abline(intercept = -4.28, slope = 0.108)
st.hydro


coef(glm(soil.h2$Hydrogen.Uptake~soil.h2$temp+soil.h2$Treatment))
coef(glm(soil.h2$Hydrogen.Uptake~soil.h2$temp))

st.hydro<-ggscatterhist(
  treatment_ordered.h2, x = "temp", y = "Hydrogen.Uptake", 
  color = "Treatment", size = 3,  shape="Season",
  palette = rev(wes_palette("Moonrise2")),
  margin.plot = "boxplot", xlab="Soil Temperature °C",
  ggtheme = theme_bw(), legend = "right", 
  ylab =expression("New-Hydrogen Fluxes (nmol H "['2']*" m "^"-2"*" s "^"-1"*")")
)
st.hydro$sp <- st.hydro$sp + 
  geom_abline(intercept = -3.61682284, slope = 0.03237263)
st.hydro


ggplot2::ggsave("Fig4-H2-temp-single-RF20211129.tiff",width = 6, height = 4, dpi = 600, units = "in")

