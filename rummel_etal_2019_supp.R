library(ape)
library(stringr)
library(phytools)
library(dplyr)

belize <- read.csv("Belize_data_upload.csv", header= TRUE)

# add rectal-muscle temp differentials to dataframe

belize$Rec_for_diff <- belize$Rectal - belize$Forearm
belize$Rec_bic_diff <- belize$Rectal - belize$Biceps
belize$Rec_pec_diff <- belize$Rectal - belize$Pectoralis

species.counts <- table(belize$Species)
species <- species.counts[species.counts > 1]

# for each species with n > 1
# make a subset of just those species
names(species)

coeff.df <- data.frame(Species = names(species),
                       Forearm.coeff = rep(NA, length(species)),
                       Pect.coeff = rep(NA, length(species)),
                       Bic.coeff = rep(NA, length(species)))
muscle <- c("Rec_for_diff","Rec_bic_diff","Rec_pec_diff")
# for each species 
par(mfrow = c(4, 4))
for (i in 1:length(names(species))) {
  
  sub.bat <- filter(belize, Species == names(species)[i]) # filter by species for each iteration
  coeff.df$Forearm.coeff[i] <- lm(sub.bat$Rec_for_diff ~ sub.bat$Air)$coefficients[2] # regress forediff on air temp for subset of values for each species
  coeff.df$Bic.coeff[i] <- lm(sub.bat$Rec_bic_diff ~ sub.bat$Air)$coefficients[2]
  coeff.df$Pect.coeff[i]<- lm(sub.bat$Rec_pec_diff ~ sub.bat$Air)$coefficients[2] 
  # pick out slope from each species' regression, put in dataframe called coeff.df
  #  jpeg(file = "forearm.jpg")

  for (j in 1:3) {
    
    temp.bat <- sub.bat[[muscle[j]]]
    plot(sub.bat$Air, temp.bat, 
         xlab = "Air temp",
         ylab = "",
         main = paste(coeff.df$Species[i]))
    abline(lm(temp.bat ~ sub.bat$Air))
    
  }
  dev.off()
}


coeff.df$Species <- str_replace(coeff.df$Species," ","_")

coeff.df$Species[5] <- "Artibeus_phaeotis"
coeff.df$Species[6] <- "Dermanura_incomitata"
coeff.df$Species[10] <- "Pteronotus_parnellii"
coeff.df$Species[12] <- "Sturnira_lilium"

# negative relationship between distal muscle temperature differentials and air temp
# is a feature of things that thermoregulate and not a special adaptation of bats
# this is especially a thing in bats because of their weird morphology 

# coeff.df.long <- melt(coeff.df, id.vars = "Species", measure.vars = c("Pect.coeff","Bic.coeff","Forearm.coeff"))

batTree<-read.tree("ShiRabosky.tree") #ShiRaboskyTree
batTree<-drop.tip(batTree, batTree$tip.label
                  [-na.omit(match(coeff.df$Species, batTree$tip.label))]) #trims tree to my taxa
plotTree(batTree,ftype="i",fsize=0.6,lwd=1, color="black")

# coeff.comp <- comparative.data(data=coeff.df, phy=batTree, names.col="Species")
# NEED TO SORT COEFF.DF SO ORDER MATCHES TIP.LABELS
order.species <- batTree$tip.label

df.test <- data.frame(Species = order.species)
coeff.sorted <- coeff.df[match(df.test$Species,coeff.df$Species),]
coeff.sorted <- remove_rownames(coeff.sorted)

forearm.sig <- phylosig(batTree,coeff.sorted$Forearm.coeff,method="lambda",test=TRUE)
biceps.sig <- phylosig(batTree,coeff.sorted$Bic.coeff,method="lambda",test=TRUE)
pectoralis.sig <- phylosig(batTree,coeff.sorted$Pect.coeff,method="lambda",test=TRUE)


################# DROP NON-PHYLLOSTOMID SPECIES

coeff.phyll <- coeff.sorted[-c(13:14),]
phyllTree <- drop.tip(batTree, batTree$tip.label
                      [-na.omit(match(coeff.phyll$Species, batTree$tip.label))])

forearmphyll.sig <- phylosig(phyllTree,coeff.phyll$Forearm.coeff,method="lambda",test=TRUE)
bicepsphyll.sig <- phylosig(phyllTree,coeff.phyll$Bic.coeff,method="lambda",test=TRUE)
pectoralisphyll.sig <- phylosig(phyllTree,coeff.phyll$Pect.coeff,method="lambda",test=TRUE)

