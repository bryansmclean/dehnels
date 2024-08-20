####################################################################################
## Dehnel's Phenomenon meta-analysis
## created 13 Dec 2023; last updated 13 Feb 2024
####################################################################################

library(ggplot2)
library(ggbiplot)
library(ggfortify)
library(viridis)
library(ape)
library(phytools)
library(lme4)
library(rstanarm)
library(loo)

#########################################
# read all-Sorex data set
#########################################

sorex <- read.csv("~/Library/CloudStorage/OneDrive-UNCG/Projects/Sorex-Dehnels/_data_/_data-extractions_and_metadata_bioclimate_20240328_.csv")

#########################################
# analysis of paired climate data
#########################################

# PCA of climate vars from all records (including cinereus); append to data file
clim.fields <- colnames(sorex)[grep("bio", colnames(sorex))]
sorex.climpca <- prcomp(sorex[,clim.fields], center = T, scale. = T) # scaling all vars to unit var, only analyzing unique localities
sorex <- data.frame(sorex, sorex.climpca$x[, 1:4], scale(sorex.climpca$x[, 1:4]))
colnames(sorex) <- gsub(".1", ".scaled", colnames(sorex), fixed = T)

#sorex.climpca$rotation
#abs(sorex.climpca$rotation[,1])[order(abs(sorex.climpca$rotation[,1]),decreasing = T)] #PC1 loadings
#abs(sorex.climpca$rotation[,2])[order(abs(sorex.climpca$rotation[,2]),decreasing = T)] #PC2 loadings
#summary(sorex.climpca) #screedata
#biplot(sorex.climpca, c(1,2)) #standard biplot

#########################################
# PLOTS of climate data
#########################################

# create a custom color and shape palette
cont <- c(sorex$continent)
cont.shp <- gsub("North America", 16, cont)
cont.shp <- gsub("Europe", 15, cont.shp)
cont.shp <- gsub("Asia", 17, cont.shp)
cont.col <- gsub("North America", "#3B9AB2", cont)
cont.col <- gsub("Europe", "#E1AF00", cont.col)
cont.col <- gsub("Asia", "#F21A00", cont.col)

# plot highlighting axis loadings
autoplot(
	object = sorex.climpca,
	# points arguments
	col = "grey60",
	alpha = .33,
	size = 3.75,
	shape = as.numeric(cont.shp),
	# loadings arguments
	loadings = TRUE,
	loadings.label = TRUE, 
	loadings.colour = 'deepskyblue3',
	loadings.label.size = 3.5,
	loadings.label.colour = "red2"
	) + 
	theme_classic()

# plot highlighting data localities
autoplot(
	object = sorex.climpca,
	# points arguments
	col = as.vector(cont.col),
	alpha = .66,
	size = 3.75,
	shape = as.numeric(cont.shp),
	# loadings arguments
	loadings = FALSE
	) + 
	theme_classic()

#########################################
# create datasets for different models
#########################################

# mass-specific and BCH-specific datasets
# MASS
sorex.mass <- sorex[which(is.na(sorex$massing_summer.winter_decrease_percent) == F),]
# set obs without a Dehnel's effect (ie positive change between summer and winter) to delta = 0
sorex.mass$massing_summer.winter_decrease_percent[c(which(sorex.mass$massing_summer.winter_decrease_percent > 0))] <- 0
# subset S. araneus and non-S. araneus for future comparison
sorex.mass.Soar <- sorex.mass[which(sorex.mass$acceptedScientificName == "Sorex araneus"),]
sorex.mass.noSoar <- sorex.mass[which(sorex.mass$acceptedScientificName != "Sorex araneus"),]

# BRAINCASE HEIGHT
sorex.bch <- sorex[which(is.na(sorex$BCH_summer.winter_decrease_percent) == F),]
# subset S. araneus and non-S. araneus for future comparison
sorex.bch.Soar <- sorex.bch[which(sorex.bch$acceptedScientificName == "Sorex araneus"),]
sorex.bch.noSoar <- sorex.bch[which(sorex.bch$acceptedScientificName != "Sorex araneus"),]

#########################################
# analysis of phylogenetic signal
#########################################

## read tree and match taxa
mammal.tre <- read.nexus("~/Library/CloudStorage/OneDrive-UNCG/Projects/Sorex-Dehnels/_data_/mammal-tree.tre")
sorex.spp <- gsub(" ", "_", unique(sorex$acceptedScientificName))
sorex.spp <- paste(sorex.spp, "SORICIDAE_EULIPOTYPHLA", sep = "_")

## check matches
# sorex.spp %in% mammal.tre$tip.label

## prune tree and match taxa
sorex.tre <- keep.tip(mammal.tre, sorex.spp)
sorex.tre$tip.label <- gsub("_SORICIDAE_EULIPOTYPHLA", "", sorex.tre$tip.label)
sorex.tre$tip.label <- gsub("_", " ", sorex.tre$tip.label)
plot(ladderize(sorex.tre, right = F))

# MASS - test of phylogenetic signal
mass.avg <- aggregate(massing_summer.winter_decrease_percent ~ acceptedScientificName, sorex.mass, FUN = mean)[,2]
names(mass.avg) <- aggregate(massing_summer.winter_decrease_percent ~ acceptedScientificName, sorex.mass, FUN = mean)[,1]
phylosig(
	tree = keep.tip(sorex.tre, intersect(sorex.tre$tip.label, names(mass.avg))),
	x = mass.avg,
	method = "K",
	test = T
	)
#contMap(tree = ladderize(keep.tip(sorex.tre, intersect(sorex.tre$tip.label, names(mass.avg))), F), x = mass.avg)

# BRAINCASEHEIGHT - test of phylogenetic signal
bch.avg <- aggregate(BCH_summer.winter_decrease_percent ~ acceptedScientificName, sorex.bch, FUN = mean)[,2]
names(bch.avg) <- aggregate(BCH_summer.winter_decrease_percent ~ acceptedScientificName, sorex.bch, FUN = mean)[,1]
phylosig(
	tree = keep.tip(sorex.tre, intersect(sorex.tre$tip.label, names(bch.avg))),
	x = bch.avg,
	method = "K",
	test = T
	)
#contMap(tree = ladderize(keep.tip(sorex.tre, intersect(sorex.tre$tip.label, names(bch.avg))), F), x = bch.avg)


#########################################
# Bayesian linear/mixed models
#########################################

## MASS ##
# full MASS model with all species/observations
mixmod.mass.1pc <- stan_lmer(
	massing_summer.winter_decrease_percent ~ 
	PC1.scaled + 
	(1 | acceptedScientificName), 
	data = sorex.mass)
mixmod.mass.2pc <- update(mixmod.mass.1pc, formula = . ~ PC1.scaled + PC2.scaled + (1 | acceptedScientificName))
mixmod.mass.3pc <- update(mixmod.mass.2pc, formula = . ~ PC1.scaled + PC2.scaled + PC3.scaled + (1 | acceptedScientificName))

# compare models containing PC1-PC3
loo1 <- loo(mixmod.mass.1pc, k_threshold = 0.7)
loo2 <- loo(mixmod.mass.2pc, k_threshold = 0.7)
loo3 <- loo(mixmod.mass.3pc, k_threshold = 0.7)
loo_compare(loo1, loo2, loo3) # model with only PC1 preferred, justifies use of only PC1 in further models

# extract posterior estimates of coefs from full model
mixmod.mass.1pc.coef <- as.data.frame(mixmod.mass.1pc)[,1:2]

# MASS model for S. araneus
mixmod.mass.soar <- stan_glm(
	massing_summer.winter_decrease_percent ~ 
	PC1.scaled, 
	data = sorex.mass.Soar
	)
mixmod.mass.soar.coef <- as.data.frame(mixmod.mass.soar)[,1:2]

# MASS model for non-S. araneus
mixmod.mass.nosoar <- stan_lmer(
	massing_summer.winter_decrease_percent ~ 
	PC1.scaled + 
	(1 | acceptedScientificName), 
	data = sorex.mass.noSoar
	)
mixmod.mass.nosoar.coef <- as.data.frame(mixmod.mass.nosoar)[,1:2]

colnames(mixmod.mass.1pc.coef) <- colnames(mixmod.mass.soar.coef) <- colnames(mixmod.mass.nosoar.coef) <- c("intercept", "PC1.scaled")

## BRAINCASE HEIGHT ##

# full BCH model with all species/observations
mixmod.bch.1pc <- stan_lmer(
	BCH_summer.winter_decrease_percent ~ 
	PC1.scaled + 
	(1 | acceptedScientificName), 
	data = sorex.bch
	)
mixmod.bch.2pc <- update(mixmod.bch.1pc, formula = . ~ PC1.scaled + PC2.scaled + (1 | acceptedScientificName))
mixmod.bch.3pc <- update(mixmod.bch.2pc, formula = . ~ PC1.scaled + PC2.scaled + PC3.scaled + (1 | acceptedScientificName))

# compare models with PC1-PC3
loo1 <- loo(mixmod.bch.1pc, k_threshold = 0.7)
loo2 <- loo(mixmod.bch.2pc, k_threshold = 0.7)
loo3 <- loo(mixmod.bch.3pc, k_threshold = 0.7)
loo_compare(loo1, loo2, loo3)

# extract posterior estimates of coefs from full model
mixmod.bch.1pc.coef <- as.data.frame(mixmod.bch.1pc)[,1:2]

# BCH model for S. araneus
mixmod.bch.soar <- stan_glm(
	BCH_summer.winter_decrease_percent ~ 
	PC1.scaled, 
	data = sorex.bch.Soar
	)
mixmod.bch.soar.coef <- as.data.frame(mixmod.bch.soar)[,1:2]

# BCH model for non-S. araneus (no Soci)
mixmod.bch.nosoar <- stan_lmer(
	BCH_summer.winter_decrease_percent ~ 
	PC1.scaled + 
	(1 | acceptedScientificName), 
	data = sorex.bch.noSoar
	)
mixmod.bch.nosoar.coef <- as.data.frame(mixmod.bch.nosoar)[,1:2]

colnames(mixmod.bch.1pc.coef) <- colnames(mixmod.bch.soar.coef) <- colnames(mixmod.bch.nosoar.coef) <- c("intercept", "PC1.scaled")

#########################################
# analyze posterior predictive intervals 
# for SOCI observations
#########################################

# create new MASS model but removing SOCI
mixmod.mass.1pc.nosoci <- stan_lmer(
	massing_summer.winter_decrease_percent ~ 
	scale(PC1) + 
	(1 | acceptedScientificName), 
	data = sorex.mass[-grep("mclean", sorex.mass$measurementID),],
	iter = 10000
	)

# predict MASS delta interval given new SOCI observations
predictive_interval(
  object = mixmod.mass.1pc.nosoci,
  newdata = sorex.mass[grep("mclean", sorex.mass$measurementID),],
  prob = 0.90,
  draws = NULL,
  re.form = NULL
)
#original values
#           5%        95%
#  -0.2070092 0.04102968

# predict MASS delta given new SOCI observations
mixmod.mass.1pc.nosoci.pp <- posterior_predict(
  object = mixmod.mass.1pc.nosoci,
  newdata = sorex.mass[grep("mclean", sorex.mass$measurementID),],
  draws = NULL,
  re.form = NULL
)
# plot posterior predictive interval w/ the new SOCI data point
hist(mixmod.mass.1pc.nosoci.pp, col = "cyan3", breaks = 50)
pt <- sorex.mass[grep("mclean", sorex.mass$measurementID),]
arrows(
	pt$massing_summer.winter_decrease_percent,
	10000000,
	pt$massing_summer.winter_decrease_percent,
	0,
	col = "red",
	lwd = 2,
	length = 0.2
	)
	
	
# create new BRAINCASE HEIGHT model but removing SOCI
mixmod.bch.1pc.nosoci <- stan_lmer(
	BCH_summer.winter_decrease_percent ~ 
	scale(PC1) + 
	(1 | acceptedScientificName), 
	data = sorex.bch[-grep("mclean", sorex.mass$measurementID),],
	iter = 10000
	)

# predict BRAINCASE HEIGHT delta interval given new SOCI observations
predictive_interval(
  object = mixmod.bch.1pc.nosoci,
  newdata = sorex.bch[grep("mclean", sorex.bch$measurementID),],
  prob = 0.90,
  draws = NULL,
  re.form = NULL
)
#original values
#           5%         95%
#   -0.150609 -0.03760405

# predict BRAINCASE HEIGHT delta given new SOCI observations
mixmod.bch.1pc.nosoci.pp <- posterior_predict(
  object = mixmod.bch.1pc.nosoci,
  newdata = sorex.bch[grep("mclean", sorex.bch$measurementID),],
  draws = NULL,
  re.form = NULL
)

# plot posterior predictive interval w/ the new SOCI data point
hist(mixmod.bch.1pc.nosoci.pp, col = "cyan3", breaks = 45)
pt <- sorex.bch[grep("mclean", sorex.bch$measurementID),]
arrows(
	pt$BCH_summer.winter_decrease_percent,
	10000000,
	pt$BCH_summer.winter_decrease_percent,
	0,
	col = "red",
	lwd = 2,
	length = 0.2
	)

#########################################
# analyze inter/intra predictive intervals 
#########################################

# do interspecific data predict intraspecific data (for S. araneus)?
# MASS
mass.nosoar.pred.soar <- predictive_interval(
  object = mixmod.mass.nosoar,
  newdata = sorex.mass.Soar[order(sorex.mass.Soar$decimallatitude),],
  prob = 0.90,
  draws = NULL,
  re.form = NULL
)
# MASS with stricter cutoff
mass.nosoar.pred.soar2 <- predictive_interval(
  object = mixmod.mass.nosoar,
  newdata = sorex.mass.Soar[order(sorex.mass.Soar$decimallatitude),],
  prob = 0.70,
  draws = NULL,
  re.form = NULL
)

plot(1:nrow(mass.nosoar.pred.soar),mass.nosoar.pred.soar[,1], ylim = c(-.4,.1), col = "white")
segments(
	1:nrow(mass.nosoar.pred.soar),
	mass.nosoar.pred.soar[,1],
	1:nrow(mass.nosoar.pred.soar),
	mass.nosoar.pred.soar[,2],
	lwd = 7,
	lend = 1,
	col = 'cyan3'
	)
segments(
	1:nrow(mass.nosoar.pred.soar2),
	mass.nosoar.pred.soar2[,1],
	1:nrow(mass.nosoar.pred.soar2),
	mass.nosoar.pred.soar2[,2],
	lwd = 0.5
	)
points(1:nrow(mass.nosoar.pred.soar),sorex.mass.Soar[order(sorex.mass.Soar$decimallatitude),]$massing_summer.winter_decrease_percent, col = "red", pch = 16)

# BRAINCASE HEIGHT
bch.nosoar.pred.soar <- predictive_interval(
  object = mixmod.bch.nosoar,
  newdata = sorex.bch.Soar[order(sorex.bch.Soar$decimallatitude),],
  prob = 0.90,
  draws = NULL,
  re.form = NULL
)
# BRAINCASE HEIGHT with stricter cutoff
bch.nosoar.pred.soar2 <- predictive_interval(
  object = mixmod.bch.nosoar,
  newdata = sorex.bch.Soar[order(sorex.bch.Soar$decimallatitude),],
  prob = 0.70,
  draws = NULL,
  re.form = NULL
)

plot(1:nrow(bch.nosoar.pred.soar),bch.nosoar.pred.soar[,1], ylim = c(-.25,0), col = "white")
segments(
	1:nrow(bch.nosoar.pred.soar),
	bch.nosoar.pred.soar[,1],
	1:nrow(bch.nosoar.pred.soar),
	bch.nosoar.pred.soar[,2],
	lwd = 12,
	lend = 1,
	col = 'cyan3'
	)
segments(
	1:nrow(bch.nosoar.pred.soar2),
	bch.nosoar.pred.soar2[,1],
	1:nrow(bch.nosoar.pred.soar2),
	bch.nosoar.pred.soar2[,2],
	lwd = 0.5
	)
points(1:nrow(bch.nosoar.pred.soar),sorex.bch.Soar[order(sorex.bch.Soar$decimallatitude),]$BCH_summer.winter_decrease_percent, col = "red", pch = 16)

#########################################
# PLOT of the full MASS model
#########################################

spp.col <- viridis(length(unique(sorex$acceptedScientificName)), option = "mako")
names(spp.col) <- unique(sorex$acceptedScientificName)[order(unique(sorex$acceptedScientificName))]
spp.col[length(spp.col)-1] <- "#00CC99"
spp.col[length(spp.col)] <- "#00FFCC"

ggplot(data = sorex.mass, 
	aes(x = PC1.scaled, y = massing_summer.winter_decrease_percent,
	colour = factor(acceptedScientificName),
	shape = factor(continent)
	)
	) +
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(
    data = mixmod.mass.1pc.coef[sample(nrow(mixmod.mass.1pc.coef), 1000),], 
	aes(intercept = intercept, slope = PC1.scaled), 
    alpha = 0.1,
    col = "grey"
	) + 
  # Plot median prediction
  geom_abline(
	aes(intercept = fixef(mixmod.mass.1pc)[1], slope = fixef(mixmod.mass.1pc)[2]), 
    col = "black"
	) + 
	geom_point(size = 3, alpha = 1) + 
	scale_colour_manual(values = spp.col) +
	theme_bw()	

#########################################
# PLOT of the full BRAINCASE HEIGHT model
#########################################

ggplot(data = sorex.bch, 
	aes(x = PC1.scaled, y = BCH_summer.winter_decrease_percent,
	colour = factor(acceptedScientificName),
	shape = factor(continent)
	)
	) +
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(
    data = mixmod.bch.1pc.coef[sample(nrow(mixmod.bch.1pc.coef), 1000),], 
	aes(intercept = intercept, slope = PC1.scaled), 
    alpha = 0.1,
    col = "grey"
	) + 
  # Plot median prediction
  geom_abline(
	aes(intercept = fixef(mixmod.bch.1pc)[1], slope = fixef(mixmod.bch.1pc)[2]), 
    col = "black"
	) + 
	geom_point(size = 3, alpha = 1) + 
	scale_colour_manual(values = spp.col) +
	theme_bw()	
	
#########################################
# get a full species legend
#########################################

ggplot(data = sorex, 
	aes(x = PC1, y = PC2,
	colour = factor(acceptedScientificName),
	shape = factor(continent))
	) +
	geom_point(alpha = .8, size = 3) + 
	scale_colour_manual(values =  spp.col) +
	theme_bw()


#########################################
# PLOTS of MASS~climate slopes inter/intra
#########################################

post.mass <- rbind(
	data.frame(
		post = mixmod.mass.1pc.coef[,2], 
		mod = rep("full", nrow(mixmod.bch.1pc.coef))
		),
	data.frame(
		post = mixmod.mass.soar.coef[,2], 
		mod = rep("soar", nrow(mixmod.mass.soar.coef))
		),
	data.frame(
		post = mixmod.mass.nosoar.coef[,2], 
		mod = rep("nosoar", nrow(mixmod.mass.nosoar.coef))
		)
	)

# quick visualize of all 3 posteriors
ggplot(post.mass, aes(post, fill = mod, color = mod)) +
	geom_histogram(position="identity", alpha = 0.3) + 
	theme_classic()

#########################################
# PLOTS of BCH~climate slopes inter/intra
#########################################

post.bch <- rbind(
	data.frame(
		post = mixmod.bch.1pc.coef[,2], 
		mod = rep("full", nrow(mixmod.bch.1pc.coef))
		),
	data.frame(
		post = mixmod.bch.soar.coef[,2], 
		mod = rep("soar", nrow(mixmod.bch.soar.coef))
		),
	data.frame(
		post = mixmod.bch.nosoar.coef[,2], 
		mod = rep("nosoar", nrow(mixmod.bch.nosoar.coef))
		)
	)

# quick visualize of all 3 posteriors
ggplot(post.bch, aes(post, fill = mod, color = mod)) +
	geom_histogram(position="identity", alpha = 0.3) + 
	theme_classic()

