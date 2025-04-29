
setwd("/Users/ThomasMac_1/Desktop/2_Projects/A_Other_Projects/Bowerbird_bowers")

library(ape)
library(phytools)

# Import the dated bowerbird consensus tree from Ericson et al. 2022. 
bowerbird.tree <- read.nexus("Bower.tree_Ericson.nex")
str(bowerbird.tree)

# Quick plot
plot(bowerbird.tree)

# Make sure the tree isultrametric (all tips the same distance from the root)
is.ultrametric(bowerbird.tree)

# List of species to drop
# (These are duplicate samples of the same species included in the Ericson et al. 2020 tree)
species_to_drop <- c(
  "CcervAM",
  "CcervF",
  "CgutcM",
  "CnucnuF",
  "CnucorF",
  "PnewF",
  "PviomiM",
  "SaureNRM",
  "SchryF" 
)

# Drop the undesired species
bowerbird.tree <- drop.tip(bowerbird.tree, species_to_drop)
str(bowerbird.tree)

# Plot the pruned tree
plot(bowerbird.tree)
write.nexus(bowerbird.tree, file = "bowerbird_tree_edited.nexus")


# Read the dataframe
ZZ<-read.csv("Bower.data.ASR.csv",row.names=1)
str(ZZ)

#######################
### Reconstructions ###
#######################

# Run either one of these for respective reconstructions: 
ZZ$Construction <- as.factor(ZZ$Construction) # Bower-construction
ZZ$Display.structure <- as.factor(ZZ$Display.structure) # Display using some kind of structure

bowers<-setNames(ZZ[,2],rownames(ZZ)) # 1 = bower construction and 2 = display around centralized structure 
bowers

# First, we test whether an ARD or ER model best fits the data #

# Estimate ancestral states under a ER model
fitER<-ace(bowers,bowerbird.tree,model="ER",type="discrete")
fitER

fitER$lik.anc

# Estimate ancestral states under a ARD model
fitARD<-ace(bowers,bowerbird.tree,model="ARD",type="discrete")
fitARD

fitARD$lik.anc

### Model comparisons 

# fit of models using AIC
AIC<-setNames(sapply(list(fitER,fitARD),AIC),c("ER", "ARD"))
AIC
aic.w(AIC)

# ER was the best fitting model for both traits, so we use that one for stochastic character mapping #

# Stochastic character mapping #

set.seed(123)

# Run MCMC
ERtrees <- make.simmap(bowerbird.tree, 
                      bowers, 
                      model="ER", 
                      nsim = 1000, 
                      Q = "mcmc",
                      vQ = 0.01, 
                      prior = list(use.empirical = TRUE), samplefreq = 10)

# set plot margins 
par(mar = c(5.1, 4.1, 2.1, 2.1))
# plot posterior density from stochastic mapping
plot(d <- density(sapply(ERtrees, function(x) x$Q[1, 2]), 
                  bw = 0.005), bty = "n", main = "", xlab = "q", xlim = c(0, 0.5), 
                  ylab = "Posterior density from MCMC", las = 1, 
                  cex.axis = 0.8)
polygon(d, col = make.transparent("black", 0.25))

# Summarize stochastic maps 
pd <- summary(ERtrees)
pd

# Plotting 
cols<-setNames(c("grey90", "black"),levels(bowers))
plot(pd, colors = cols, fsize = 0.6, ftype = "i", lwd = 1.2, 
     offset = 0.4, ylim = c(-1, Ntip(bowerbird.tree)), 
     cex = c(0.5, 0.3))
# Add legend
legend("bottomleft", legend = levels(bowers), pch = 22, pt.cex = 1.5, pt.bg = cols, bty = "n", cex = 1.2)

# Create density map 
bower.densityMap <- densityMap(ERtrees, 
                               states = levels(bowers) [2:1], plot = FALSE)
# Update colour gradient 
bower.densityMap <- setMap(bower.densityMap, cols[2:1])
# Plot it 
plot(bower.densityMap, fsize = c(0.5, 0.7), lwd = c(3, 4))





