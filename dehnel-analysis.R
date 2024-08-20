####################################################################################
## Dehnel's Phenomenon - Sorex cinereus analysis
## created 13 Dec 2023; last updated 13 Feb 2024
####################################################################################

library(ggplot2)

#########################################
# read specimen data & metadata
#########################################
soci.df <- read.csv("SOCI_metadata_crania_femora_ageclass.csv")
soci.df$season <- factor(soci.df$season, levels = c("spring", "summer", "fall", "winter"))
soci.df$ageclass <- factor(soci.df$ageclass, levels = c("subad_sp","subad_su","subad_fa","subad_wi","ad_sp","ad_su","ad_fa","ad_wi"))

#########################################
# averages by season
#########################################

aggregate(soci.df[,c("weight_g")] ~ soci.df$season, FUN = mean)
aggregate(soci.df[,c("BCH")] ~ soci.df$season, FUN = mean)
aggregate(soci.df[,c("TL")] ~ soci.df$season, FUN = mean)

#########################################
# ANOVAs by season
#########################################

weight.aov <- aov(weight_g ~ ageclass, data = soci.df)
TukeyHSD(weight.aov)
bch.aov <- aov(BCH ~ ageclass, data = soci.df)
TukeyHSD(bch.aov)
tl.aov <- aov(TL ~ ageclass, data = soci.df)
TukeyHSD(tl.aov)


#########################################
# BOXPLOTS of seasonal change
#########################################
seas.cols <- c("#E2D200", "#E58601", "#B40F20", "#46ACC8")

# body mass
ggplot(
	soci.df, 
	aes(x = ageclass, y = weight_g, colour = season, fill = season)
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
	
# braincase heigh ratio
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

# femur length
ggplot(
	soci.df, 
	aes(x = ageclass, y = TL, colour = season, fill = season)
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
