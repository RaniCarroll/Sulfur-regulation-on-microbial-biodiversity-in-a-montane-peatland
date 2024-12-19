###########################################################################################
### Carroll et al. Sulfur regulation on microbial biodiversity in a montane peatland_RDA analysis script

#####Install packages
install.packages("tidyverse")
install.packages("psych")
install.packages("adespatial")

# Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)

# Import genetic data
allele_freqs = read.csv("MB_ASV_20240624.csv", row.names = 1, check.names = FALSE)
# Hellinger transform the community data (as per Johnson et al. 2023)
allele_freqs.spe.hel <- decostand(allele_freqs, method = "hellinger")

# Import spatial data
dbmem.raw = read.csv("MB_Cores_20240624.csv", row.names = 1)

# Import environmental data
env.raw = read.csv("MB_Soil_20241212.csv", row.names = 1)

# Set seed
set.seed(123)

#
# Multicollinearity checks
#
#--------------#

# Plot and run correlation test on environmental variables
pairs.panels(env.raw, scale = TRUE)

# Remove correlated variables
env.data = subset(env.raw, select = -c(NA., TAA))
pairs.panels(env.data, scale = TRUE)

############################################ 

#--------------#
#
# Redundancy analysis 
#
#--------------#

# Perform RDA with all variables
rda1 = rda(allele_freqs.spe.hel ~ ., data = env.raw, scale = TRUE)
rda1

# Model summaries
RsquareAdj(rda1) # adjusted Rsquared 
vif.cca(rda1) # variance inflation factor (<10 OK)
anova.cca(rda1, permutations = 1000) # full model
anova.cca(rda1, permutations = 1000, by="margin") # per variable 

# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "constrained"))
screeplot(rda1)

# Create a dataframe to correctly colour regions
col_dframe = data.frame("site" = rownames(allele_freqs))

# Function to add regional labels to dataframe (based on zone - oxidised, transition and reduced)
addregion = function(x){
  # If pop label is present function will output the region
  if(x=="H1-0-10"|x=="H2-0-10"|x=="H3-0-10"|x=="H4-0-10"|x=="H5-0-10"|x=="H6-0-10") y = "Oxidised"
  if(x=="H1-40-50"|x=="H2-40-50"|x=="H3-40-50-1"|x=="H3-40-50-2"|x=="H4-40-50"|x=="H5-40-50-1"|x=="H5-40-50-2"|x=="H6-40-50-1"|x=="H6-40-50-2") y = "Mid-depth"
  if(x=="H1-90-100"|x=="H2-70-80-1"|x=="H2-70-80-2"|x=="H3-60-70"|x=="H4-130-140"|x=="H5-70-80"|x=="H6-100-110") y = "Reduced"
  return(y)
}

# Add regional labels
col_dframe$region = sapply(col_dframe$site, addregion)

# Add factor levels
region_order = c("Oxidised","Mid-depth","Reduced")
col_dframe$region = factor(col_dframe$region, levels = region_order)

# Create colour scheme
# blue=#377EB8, green=#7FC97F, orange=#FDB462, red=#E31A1C
cols = c("red","blue","orange")

# Visualise results of RDA
##to save directly as png file

pdf("rda_11.pdf", width = 8, height = 7)
plot(rda1, type="n", scaling = 3, xlab="RDA1 (21% explained)", ylab="RDA2 (15% explained)")
title("Mount Banks redundancy analysis")
# SITES
points(rda1, display="sites", pch=21, scaling=3, cex=1.5, col="black",
       bg=cols[col_dframe$region]) # sites
# text(rda1, display="sites", scaling = 3, col="black", font=2, pos=4)
# PREDICTORS
text(rda1, display="bp", scaling=3, col="black", cex=1, lwd=2)
# SNPS
# text(rda1, display="species", scaling = 3, col="blue", cex=0.7, pos=4) # SNPs
# LEGEND
legend("bottomleft", legend=levels(col_dframe$region), bty="n", col="black",
       pch=21, cex=1.2, pt.bg=cols)
# OTHER LABELS
adj.R2 = round(RsquareAdj(rda1)$adj.r.squared, 3)
mtext(bquote(italic("R")^"2"~"= "~.(adj.R2)), side = 3, adj = 0.5)
dev.off()

#########################################################################################################

### RDA analysis with plot by sites
# Perform RDA with all variables
rda2 = rda(allele_freqs.spe.hel ~ ., data = env.raw, scale = TRUE)
rda2

# Model summaries
RsquareAdj(rda2) # adjusted Rsquared 
vif.cca(rda2) # variance inflation factor (<10 OK)
anova.cca(rda2, permutations = 1000) # full model
anova.cca(rda2, permutations = 1000, by="margin") # per variable 

# Variance explained by each canonical axis
summary(eigenvals(rda2, model = "constrained"))
screeplot(rda2)

# Create a dataframe to correctly colour regions
col_dframe = data.frame("site" = rownames(allele_freqs))

# Function to add regional labels to dataframe (by core/site)
addregion = function(x){
  # If pop label is present function will output the region
  if(x=="H1-0-10"|x=="H1-40-50"|x=="H1-90-100") y = "H1"
  if(x=="H2-0-10"|x=="H2-40-50"|x=="H2-70-80-1"|x=="H2-70-80-2") y = "H2"
  if(x=="H3-0-10"|x=="H3-40-50-1"|x=="H3-40-50-2"|x=="H3-60-70") y = "H3"
  if(x=="H4-0-10"|x=="H4-40-50"|x=="H4-130-140") y = "H4"
  if(x=="H5-0-10"|x=="H5-40-50-1"|x=="H5-40-50-2"|x=="H5-70-80") y = "H5"
  if(x=="H6-0-10"|x=="H6-40-50-1"|x=="H6-40-50-2"|x=="H6-100-110") y = "H6"
  return(y)
}

# Add regional labels
col_dframe$region = sapply(col_dframe$site, addregion)

# Add factor levels
region_order = c("H1","H2","H3","H4", "H5", "H6")
col_dframe$region = factor(col_dframe$region, levels = region_order)

# Create colour scheme
# blue=#377EB8, green=#7FC97F, orange=#FDB462, red=#E31A1C
cols = c("royalblue","maroon","palegreen3", "mediumpurple2", "skyblue1", "orange")

# Visualise results of RDA
##to save directly as png file

pdf("rda_14.pdf", width = 8, height = 7)
plot(rda2, type="n", scaling = 3, xlab="RDA1 (21% explained)", ylab="RDA2 (15% explained)")
title("Mount Banks redundancy analysis")
# SITES
points(rda2, display="sites", pch=21, scaling=3, cex=1.5, col="black",
       bg=cols[col_dframe$region]) # sites
# text(rda1, display="sites", scaling = 3, col="black", font=2, pos=4)
# PREDICTORS
text(rda2, display="bp", scaling=3, col="black", cex=1, lwd=2)
# SNPS
# text(rda1, display="species", scaling = 3, col="blue", cex=0.7, pos=4) # SNPs
# LEGEND
legend("bottomleft", legend=levels(col_dframe$region), bty="n", col="black",
       pch=21, cex=1.2, pt.bg=cols)
# OTHER LABELS
adj.R2 = round(RsquareAdj(rda2)$adj.r.squared, 3)
mtext(bquote(italic("R")^"2"~"= "~.(adj.R2)), side = 3, adj = 0.5)
dev.off()

#########################################################################################################
