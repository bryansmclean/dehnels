####################################################################################
## Dehnel's Phenomenon - Sorex cinereus analysis
## B.S McLean et al. "Seasonal body size plasticity and the generality of Dehnel’s Phenomenon in Sorex shrews."
## created 13 Dec 2023; last updated 20 jan 2025
####################################################################################

library(ggplot2)


#########################################
# read specimen data & metadata
#########################################

soci.df <- read.csv("SOCI_metadata_crania_femora_ageclass-v2.csv")
soci.df$season <- factor(soci.df$season, levels = c("spring", "summer", "fall", "winter"))
soci.df$ageclass <- factor(soci.df$ageclass, levels = c("subad_sp","subad_su","subad_fa","subad_wi","ad_sp","ad_su","ad_fa","ad_wi"))


#########################################
# average traits by season
#########################################

aggregate(soci.df[,c("weight_g")] ~ soci.df$season, FUN = mean) #mass
aggregate(soci.df[,c("BCH")] ~ soci.df$season, FUN = mean) #braincase height
aggregate(soci.df[,c("TL")] ~ soci.df$season, FUN = mean) # femur length
aggregate((soci.df$TL/soci.df$total_length_mm) ~ soci.df$season, FUN = mean)


#########################################
# ANOVAs by season
#########################################

weight.aov <- aov(weight_g ~ ageclass, data = soci.df)
TukeyHSD(weight.aov)
bodycond.aov <- aov(weight_g/total_length_mm ~ ageclass, data = soci.df)
TukeyHSD(bodycond.aov)
bch.aov <- aov(BCH ~ ageclass, data = soci.df)
TukeyHSD(bch.aov)
tl.aov <- aov(TL ~ ageclass, data = soci.df)
TukeyHSD(tl.aov)

tl.aov <- aov(TL/weight_g ~ ageclass, data = soci.df)
TukeyHSD(tl.aov)
summary(lm(soci.df$TL ~ soci.df$weight_g))


#########################################
# BOXPLOTS of seasonal change
#########################################
seas.cols <- c("#E2D200", "#E58601", "#B40F20", "#46ACC8")

## MASS BOXPLOTS ##

ggplot(
	soci.df, 
	aes(x = body_length_mm, y = TL, colour = ageclass, fill = ageclass)
	) + 
geom_point()

# body mass
ggplot(
	soci.df, 
	aes(x = ageclass, y = weight_g/total_length_mm, colour = season, fill = season)
	) + 
geom_boxplot(
	alpha = 0.1
	) + 
geom_jitter(
	width = .1, cex = 2, alpha = 0.4
	) +
scale_colour_manual(values = seas.cols) +
scale_fill_manual(values = seas.cols) +
theme_bw()
	
## BRAINCASE HEIGHT BOXPLOTS ##

ggplot(
	soci.df, 
	aes(x = ageclass, y = (BCH*2)/(SKL+BCW), colour = season, fill = season)
	) + 
geom_boxplot(
	alpha = 0.1
	) + 
geom_jitter(
	width = .1, cex = 2, alpha = 0.4
	) +
scale_colour_manual(values = seas.cols) +
scale_fill_manual(values = seas.cols) +
theme_bw()

## FEMUR LENGTH BOXPLOTS ##

ggplot(
	soci.df, 
	aes(x = ageclass, y = TL/total_length_mm, colour = season, fill = season)
	) + 
geom_boxplot(
	alpha = 0.1
	) + 
geom_jitter(
	width = .1, cex = 2, alpha = 0.4
	) +
scale_colour_manual(values = seas.cols) +
scale_fill_manual(values = seas.cols) +
theme_bw()
