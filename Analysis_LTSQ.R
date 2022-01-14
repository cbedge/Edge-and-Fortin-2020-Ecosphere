library(dplyr)
library(ggplot2)
library(lme4)
library(effects)
library(ggpubr)
library(tidyverse)

setwd("N:/Population Model Results") #Set directory to store outputs
setwd("Z:/cedge/Population Model Results") #Laptop

full_rand <- read.csv(file="Full_Random_LatinSQ.csv", header=T)
names(full_rand)

full_rand$nDisturbed <- (full_rand$nDisturbed/60) * 100

###Statistical analyses###

##FULL NETWORK##
hist(full_rand$metaR) #left skew
hist(full_rand$ExtinctTime) #right skew with a lot of values at 500
hist(full_rand$Attract)

full_rand_center <- full_rand # new dataset for scaling

#METAPOPULATION GROWTH RATE#

#Scale all the variables to SD and variance
full_rand_center$metaR <- scale(full_rand_center$metaR) #mean centre and scale to SD
full_rand_center$nDisturbed <- scale(full_rand_center$nDisturbed)
full_rand_center$Sev <- scale(full_rand_center$Sev)
full_rand_center$Attract <- scale(full_rand_center$Attract)
full_rand_center$Disperse <- scale(full_rand_center$Disperse)
#hist(full_rand_center$nDisturbed)
#hist(full_rand_center$Sev)

hist(full_rand_center$metaR) #left skew with negative values
min(full_rand_center$metaR) #value to shift the data to the right to get rid of negative values

#mod1.meta <- glm((metaR + 6.513615) ~ nDisturbed * Sev * Attract * Disperse, family = Gamma(link="inverse"), data = full_rand_center)
mod1.meta <- glm(metaR ~ nDisturbed * Sev * Attract * Disperse, data = full_rand_center)
summary(mod1.meta)
write.csv(summary(mod1.meta)$coefficients, file="full_rand_meta_coef.csv")

#mod2 and mod3 are analyses to look at model simplificaiton
mod2.meta <- glm((metaR+6.513615) ~ nDisturbed + Sev + Attract + Disperse, data = full_rand_center, family = Gamma(link="inverse"))
summary(mod2.meta)

mod3.meta <- glm((metaR+6.513615) ~ nDisturbed + Sev + Attract + Disperse + nDisturbed:Sev + nDisturbed:Disperse + nDisturbed:Sev:Disperse + nDisturbed:Attract:Disperse + nDisturbed:Sev:Attract:Disperse, data = full_rand_center, family = Gamma(link="inverse"))
summary(mod3.meta)

metacoef <- as.data.frame(cbind(summary(mod1.meta)$coefficients[,1], confint.default(mod1.meta)))
metacoef <- metacoef[-1,]
colnames(metacoef) <- c("mean_MetaR", "lowerCI", "upperCI")
metacoef$var <- rownames(metacoef)
metacoef <- metacoef[order(-metacoef$mean),]
metacoef$var <- factor(metacoef$var, levels=metacoef$var) 
metacoef$Network <- "Full"
metacoef$mSelect <- "Random"

min(metacoef$lowerCI)
max(metacoef$upperCI)

p.meta.full <- ggplot(metacoef, aes(x=var, y=mean_MetaR)) +
  geom_point() +
  geom_errorbar(metacoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Full - Metapopulation Growth Rate") +
  #scale_y_continuous(breaks=c(-0.0002, 0, 0.0002, 0.0004), labels = c("-0.002", "0", "0.002", "0.004"), name = "Effect Size") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )
  

jpeg('full_meta.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.meta.full
dev.off()




#EXTINCTION#
hist(full_rand_center$ExtinctTime) #right skew with lots of values at one location
max(full_rand$ExtinctTime)


#probability of extinction
full_rand_center$Extinct[full_rand_center$ExtinctTime == 501] <- 0 #survive = 1
full_rand_center$Extinct[full_rand_center$ExtinctTime < 501] <- 1 #extinct = 0

mod1.ex <- glm(Extinct ~ nDisturbed * Sev * Attract * Disperse, data = full_rand_center, family = binomial(link = "logit"))
summary(mod1.ex)
write.csv(summary(mod1.ex)$coefficients, file="full_rand_exprob_coef.csv")


excoef <- as.data.frame(cbind(summary(mod1.ex)$coefficients[,1], confint.default(mod1.ex)))
excoef <- excoef[-1,]
colnames(excoef) <- c("mean_Exprob", "lowerCI", "upperCI")
excoef$var <- rownames(excoef)
excoef <- excoef[order(-excoef$mean),]
excoef$var <- factor(excoef$var, levels=excoef$var) 
excoef$Network <- "Full"
excoef$mSelect <- "Random"

min(excoef$lowerCI)
max(excoef$upperCI)

p.exprob.full <- ggplot(excoef, aes(x=var, y=mean_Exprob)) +
  geom_point() +
  geom_errorbar(excoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Full - Probability of Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('full_exprob.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.exprob.full
dev.off()




#Time to extinction
mod2.ex <- glm(ExtinctTime ~ nDisturbed * Sev * Attract * Disperse, data = full_rand_center)
summary(mod2.ex)
write.csv(summary(mod2.ex)$coefficients, file="full_rand_extime_coef.csv")

excoef.mod2 <- as.data.frame(cbind(summary(mod2.ex)$coefficients[,1], confint.default(mod2.ex)))
excoef.mod2 <- excoef.mod2[-1,]
colnames(excoef.mod2) <- c("mean_ExTime", "lowerCI", "upperCI")
excoef.mod2$var <- rownames(excoef.mod2)
excoef.mod2 <- excoef.mod2[order(-excoef.mod2$mean),]
excoef.mod2$var <- factor(excoef.mod2$var, levels=excoef.mod2$var) 
excoef.mod2$Network <- "Full"
excoef.mod2$mSelect <- "Random"

p.extime.full <- ggplot(excoef.mod2, aes(x=var, y=mean_ExTime)) +
  geom_point() +
  geom_errorbar(excoef.mod2, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Full - Time to Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )

jpeg('full_extime.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.extime.full
dev.off()


####ERDOS RENYI NETWORK
erdos_rand <- read.csv(file="Erdos_Random_LatinSQ.csv", header=T)
names(erdos_rand)

erdos_rand$nDisturbed <- (erdos_rand$nDisturbed/60) * 100

hist(erdos_rand$metaR) #left skew
hist(erdos_rand$ExtinctTime) #right skew with a lot of values at 500
hist(erdos_rand$Attract)

erdos_rand_center <- erdos_rand # new dataset for scaling

#METAPOPULATION GROWTH RATE#

#Scale all the variables to SD and variance
erdos_rand_center$metaR <- scale(erdos_rand_center$metaR) #mean centre and scale to SD
erdos_rand_center$nDisturbed <- scale(erdos_rand_center$nDisturbed)
erdos_rand_center$Sev <- scale(erdos_rand_center$Sev)
erdos_rand_center$Attract <- scale(erdos_rand_center$Attract)
erdos_rand_center$Disperse <- scale(erdos_rand_center$Disperse)
#hist(full_rand_center$nDisturbed)
#hist(full_rand_center$Sev)

hist(erdos_rand_center$metaR) #left skew with negative values
min(erdos_rand_center$metaR) #value to shift the data to the right to get rid of negative values

#mod1.meta <- glm((metaR + 6.513615) ~ nDisturbed * Sev * Attract * Disperse, family = Gamma(link="inverse"), data = full_rand_center)
erdos.meta.mod1 <- glm(metaR ~ nDisturbed * Sev * Attract * Disperse, data = erdos_rand_center)
summary(erdos.meta.mod1)
write.csv(summary(erdos.meta.mod1)$coefficients, file="erdos_rand_meta_coef.csv")

erdos.metacoef <- as.data.frame(cbind(summary(erdos.meta.mod1)$coefficients[,1], confint.default(erdos.meta.mod1)))
erdos.metacoef <- erdos.metacoef[-1,]
colnames(erdos.metacoef) <- c("mean_MetaR", "lowerCI", "upperCI")
erdos.metacoef$var <- rownames(erdos.metacoef)
erdos.metacoef <- erdos.metacoef[order(-erdos.metacoef$mean),]
erdos.metacoef$var <- factor(erdos.metacoef$var, levels=erdos.metacoef$var) 
erdos.metacoef$Network <- "Erdos"
erdos.metacoef$mSelect <- "Random"
  
min(erdos.metacoef$lowerCI)
max(erdos.metacoef$upperCI)

p.erdos.meta.rand <- ggplot(erdos.metacoef, aes(x=var, y=mean_MetaR)) +
  geom_point() +
  geom_errorbar(erdos.metacoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos - Metapopulation Growth Rate") +
  #scale_y_continuous(breaks=c(-0.0002, 0, 0.0002, 0.0004), labels = c("-0.002", "0", "0.002", "0.004"), name = "Effect Size") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('erdos_rand_meta.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.erdos.meta.rand
dev.off()

#EXTINCTION#
hist(erdos_rand_center$ExtinctTime) #right skew with lots of values at one location
max(erdos_rand$ExtinctTime)


#probability of extinction
erdos_rand_center$Extinct[erdos_rand_center$ExtinctTime == 501] <- 0 #survive = 0
erdos_rand_center$Extinct[erdos_rand_center$ExtinctTime < 501] <- 1 #extinct = 1

erdos.ex.mod1 <- glm(Extinct ~ nDisturbed * Sev * Attract * Disperse, data = erdos_rand_center, family = binomial(link = "logit"))
summary(erdos.ex.mod1)
write.csv(summary(erdos.ex.mod1)$coefficients, file="erdos_rand_exprob_coef.csv")

erdos.excoef <- as.data.frame(cbind(summary(erdos.ex.mod1)$coefficients[,1], confint.default(erdos.ex.mod1)))
erdos.excoef <- erdos.excoef[-1,]
colnames(erdos.excoef) <- c("mean_Exprob", "lowerCI", "upperCI")
erdos.excoef$var <- rownames(erdos.excoef)
erdos.excoef <- erdos.excoef[order(-erdos.excoef$mean),]
erdos.excoef$var <- factor(erdos.excoef$var, levels=erdos.excoef$var) 
erdos.excoef$Network <- "Erdos"
erdos.excoef$mSelect <- "Random"

min(erdos.excoef$lowerCI)
max(erdos.excoef$upperCI)

p.erdos.exprob <- ggplot(erdos.excoef, aes(x=var, y=mean_Exprob)) +
  geom_point() +
  geom_errorbar(erdos.excoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos - Probability of Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('erdos_exprob.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.erdos.exprob
dev.off()


#Time to extinction
erdos.ex.mod2 <- glm(ExtinctTime ~ nDisturbed * Sev * Attract * Disperse, data = erdos_rand_center)
summary(erdos.ex.mod2)
write.csv(summary(erdos.ex.mod2)$coefficients, file="erdos_rand_extime_coef.csv")

erdos.excoef.mod2 <- as.data.frame(cbind(summary(erdos.ex.mod2)$coefficients[,1], confint.default(erdos.ex.mod2)))
erdos.excoef.mod2 <- erdos.excoef.mod2[-1,]
colnames(erdos.excoef.mod2) <- c("mean_ExTime", "lowerCI", "upperCI")
erdos.excoef.mod2$var <- rownames(erdos.excoef.mod2)
erdos.excoef.mod2 <- erdos.excoef.mod2[order(-erdos.excoef.mod2$mean),]
erdos.excoef.mod2$var <- factor(erdos.excoef.mod2$var, levels=erdos.excoef.mod2$var) 
erdos.excoef.mod2$Network <- "Erdos"
erdos.excoef.mod2$mSelect <- "Random"

p.extime.erdos <- ggplot(erdos.excoef.mod2, aes(x=var, y=mean_ExTime)) +
  geom_point() +
  geom_errorbar(erdos.excoef.mod2, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos - Time to Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )

jpeg('erdos_extime.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.extime.erdos
dev.off()


####TREE HABITAT NETWORK

tree_rand <- read.csv(file="Tree_Random_LatinSQ.csv", header=T)
names(tree_rand)

tree_rand$nDisturbed <- (tree_rand$nDisturbed/60) * 100

hist(tree_rand$metaR) #left skew
hist(tree_rand$ExtinctTime) #right skew with a lot of values at 500
hist(tree_rand$Attract)

tree_rand_center <- tree_rand # new dataset for scaling

#METAPOPULATION GROWTH RATE#

#Scale all the variables to SD and variance
tree_rand_center$metaR <- scale(tree_rand_center$metaR) #mean centre and scale to SD
tree_rand_center$nDisturbed <- scale(tree_rand_center$nDisturbed)
tree_rand_center$Sev <- scale(tree_rand_center$Sev)
tree_rand_center$Attract <- scale(tree_rand_center$Attract)
tree_rand_center$Disperse <- scale(tree_rand_center$Disperse)
#hist(tree_rand_center$nDisturbed)
#hist(tree_rand_center$Sev)

hist(tree_rand_center$metaR) #left skew with negative values
min(tree_rand_center$metaR) #value to shift the data to the right to get rid of negative values

#mod1.meta <- glm((metaR + 6.513615) ~ nDisturbed * Sev * Attract * Disperse, family = Gamma(link="inverse"), data = full_rand_center)
tree.meta.mod1 <- glm(metaR ~ nDisturbed * Sev * Attract * Disperse, data = tree_rand_center)
summary(tree.meta.mod1)
write.csv(summary(tree.meta.mod1)$coefficients, file="tree_rand_meta_coef.csv")

tree.metacoef <- as.data.frame(cbind(summary(tree.meta.mod1)$coefficients[,1], confint.default(tree.meta.mod1)))
tree.metacoef <- tree.metacoef[-1,]
colnames(tree.metacoef) <- c("mean_MetaR", "lowerCI", "upperCI")
tree.metacoef$var <- rownames(tree.metacoef)
tree.metacoef <- tree.metacoef[order(-tree.metacoef$mean),]
tree.metacoef$var <- factor(tree.metacoef$var, levels=tree.metacoef$var) 
tree.metacoef$Network <- "Tree"
tree.metacoef$mSelect <- "Random"

min(tree.metacoef$lowerCI)
max(tree.metacoef$upperCI)

p.tree.meta.rand <- ggplot(tree.metacoef, aes(x=var, y=mean_MetaR)) +
  geom_point() +
  geom_errorbar(tree.metacoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree - Metapopulation Growth Rate") +
  #scale_y_continuous(breaks=c(-0.0002, 0, 0.0002, 0.0004), labels = c("-0.002", "0", "0.002", "0.004"), name = "Effect Size") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('tree_rand_meta.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.tree.meta.rand
dev.off()




#EXTINCTION#
hist(tree_rand_center$ExtinctTime) #right skew with lots of values at one location
max(tree_rand$ExtinctTime)


#probability of extinction
tree_rand_center$Extinct[tree_rand_center$ExtinctTime == 501] <- 0 #survive = 0
tree_rand_center$Extinct[tree_rand_center$ExtinctTime < 501] <- 1 #extinct = 1

tree.ex.mod1 <- glm(Extinct ~ nDisturbed * Sev * Attract * Disperse, data = tree_rand_center, family = binomial(link = "logit"))
summary(tree.ex.mod1)
write.csv(summary(tree.ex.mod1)$coefficients, file="tree_rand_exprob_coef.csv")


tree.excoef <- as.data.frame(cbind(summary(tree.ex.mod1)$coefficients[,1], confint.default(tree.ex.mod1)))
tree.excoef <- tree.excoef[-1,]
colnames(tree.excoef) <- c("mean_Exprob", "lowerCI", "upperCI")
tree.excoef$var <- rownames(tree.excoef)
tree.excoef <- tree.excoef[order(-tree.excoef$mean),]
tree.excoef$var <- factor(tree.excoef$var, levels=tree.excoef$var) 
tree.excoef$Network <- "Tree"
tree.excoef$mSelect <- "Random"

min(tree.excoef$lowerCI)
max(tree.excoef$upperCI)

p.tree.exprob <- ggplot(tree.excoef, aes(x=var, y=mean_Exprob)) +
  geom_point() +
  geom_errorbar(tree.excoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree - Probability of Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('tree_exprob.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.tree.exprob
dev.off()


#Time to extinction
tree.ex.mod2 <- glm(ExtinctTime ~ nDisturbed * Sev * Attract * Disperse, data = tree_rand_center)
summary(tree.ex.mod2)
write.csv(summary(tree.ex.mod2)$coefficients, file="tree_rand_extime_coef.csv")

tree.excoef.mod2 <- as.data.frame(cbind(summary(tree.ex.mod2)$coefficients[,1], confint.default(tree.ex.mod2)))
tree.excoef.mod2 <- tree.excoef.mod2[-1,]
colnames(tree.excoef.mod2) <- c("mean_ExTime", "lowerCI", "upperCI")
tree.excoef.mod2$var <- rownames(tree.excoef.mod2)
tree.excoef.mod2 <- tree.excoef.mod2[order(-tree.excoef.mod2$mean),]
tree.excoef.mod2$var <- factor(tree.excoef.mod2$var, levels=tree.excoef.mod2$var) 
tree.excoef.mod2$Network <- "Tree"
tree.excoef.mod2$mSelect <- "Random"

p.extime.tree <- ggplot(tree.excoef.mod2, aes(x=var, y=mean_ExTime)) +
  geom_point() +
  geom_errorbar(tree.excoef.mod2, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree - Time to Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )

jpeg('tree_extime.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.extime.tree
dev.off()



####COMBINING FIGURES FROM ALL NETWORK TYPES

#mean metapopulation growth rate, random selection
metar.rand <- rbind(metacoef, erdos.metacoef, tree.metacoef)
metar.rand$Network <- factor(metar.rand$Network, levels = c("Full", "Erdos", "Tree"))

metar.rand$var <- gsub("nDisturbed", "Disturbed", metar.rand$var)
metar.rand$var <- gsub("Sev", "Severity", metar.rand$var)

metar.plot.rand <- ggplot(metar.rand, aes(x=reorder(var, -mean_MetaR), y=mean_MetaR, colour=Network, order=Network)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(metar.rand, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Mean Metapopulation Growth Rate") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Network Type") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('MetaR_Rand.jpeg', width=10, height=8, units="cm", res=300)
metar.plot.rand
dev.off()

#Extinction probability random selection
exprob.rand <- rbind(excoef, erdos.excoef, tree.excoef)
exprob.rand$Network <- factor(exprob.rand$Network, levels = c("Full", "Erdos", "Tree"))

exprob.rand$var <- gsub("nDisturbed", "Disturbed", exprob.rand$var)
exprob.rand$var <- gsub("Sev", "Severity", exprob.rand$var)

exprob.plot.rand <- ggplot(exprob.rand, aes(x=reorder(var, mean_Exprob), y=mean_Exprob, colour=Network)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(exprob.rand, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Extinction Probability") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Network Type") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.8, 0.75),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('EXProb_Rand.jpeg', width=10, height=8, units="cm", res=300)
exprob.plot.rand
dev.off()

#Time to extinction random selection

extime.rand <- rbind(excoef.mod2, erdos.excoef.mod2, tree.excoef.mod2)
extime.rand$Network <- factor(extime.rand$Network, levels = c("Full", "Erdos", "Tree"))
extime.rand$var <- gsub("nDisturbed", "Disturbed", extime.rand$var)
extime.rand$var <- gsub("Sev", "Severity", extime.rand$var)

extime.plot.rand <- ggplot(extime.rand, aes(x=reorder(var, -mean_ExTime), y=mean_ExTime, colour=Network)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(extime.rand, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Time to Extinction") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Network Type") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('ExTime_Rand.jpeg', width=10, height=8, units="cm", res=300)
extime.plot.rand
dev.off()

#Writing all three figures together
tiff('Rand_Selection.tif', width=10, height=20, units="cm", res=300)
ggarrange(metar.plot.rand, exprob.plot.rand, extime.plot.rand, ncol = 1, nrow = 3, align="hv",
          labels="AUTO", label.x=0.95, label.y=0.9, font.label = list(size=7, face ="plain"),
          common.legend = TRUE, legend="bottom")
dev.off()

###TREE NETWORK WITH TOP CONNECTIVITY SELECTED###
tree_topcon <- read.csv(file="Tree_TopCon_LatinSQ.csv", header=T)
names(tree_topcon)

tree_topcon$nDisturbed <- (tree_topcon$nDisturbed/60) * 100

hist(tree_topcon$metaR) #left skew
hist(tree_topcon$ExtinctTime) #right skew with a lot of values at 500
hist(tree_topcon$Attract)

tree_topcon_center <- tree_topcon # new dataset for scaling

#METAPOPULATION GROWTH RATE#

#Scale all the variables to SD and variance
tree_topcon_center$metaR <- scale(tree_topcon_center$metaR) #mean centre and scale to SD
tree_topcon_center$nDisturbed <- scale(tree_topcon_center$nDisturbed)
tree_topcon_center$Sev <- scale(tree_topcon_center$Sev)
tree_topcon_center$Attract <- scale(tree_topcon_center$Attract)
tree_topcon_center$Disperse <- scale(tree_topcon_center$Disperse)
#hist(tree_topcon_center$nDisturbed)
#hist(tree_topcon_center$Sev)

hist(tree_topcon_center$metaR) #left skew with negative values
min(tree_topcon_center$metaR) #value to shift the data to the right to get rid of negative values

#mod1.meta <- glm((metaR + 6.513615) ~ nDisturbed * Sev * Attract * Disperse, family = Gamma(link="inverse"), data = full_rand_center)
tree.topcon.meta.mod1 <- glm(metaR ~ nDisturbed * Sev * Attract * Disperse, data = tree_topcon_center)
summary(tree.topcon.meta.mod1)
write.csv(summary(tree.topcon.meta.mod1)$coefficients, file="tree_topcon_meta_coef.csv")

tree.topcon.metacoef <- as.data.frame(cbind(summary(tree.topcon.meta.mod1)$coefficients[,1], confint.default(tree.topcon.meta.mod1)))
tree.topcon.metacoef <- tree.topcon.metacoef[-1,]
colnames(tree.topcon.metacoef) <- c("mean_MetaR", "lowerCI", "upperCI")
tree.topcon.metacoef$var <- rownames(tree.topcon.metacoef)
tree.topcon.metacoef <- tree.topcon.metacoef[order(-tree.topcon.metacoef$mean),]
tree.topcon.metacoef$var <- factor(tree.topcon.metacoef$var, levels=tree.topcon.metacoef$var) 
tree.topcon.metacoef$Network <- "Tree"
tree.topcon.metacoef$mSelect <- "Top"

min(tree.topcon.metacoef$lowerCI)
max(tree.topcon.metacoef$upperCI)

p.tree.meta.topcon <- ggplot(tree.topcon.metacoef, aes(x=var, y=mean_MetaR)) +
  geom_point() +
  geom_errorbar(tree.topcon.metacoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree Top Con - Metapopulation Growth Rate") +
  #scale_y_continuous(breaks=c(-0.0002, 0, 0.0002, 0.0004), labels = c("-0.002", "0", "0.002", "0.004"), name = "Effect Size") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('tree_topcon_meta.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.tree.meta.topcon
dev.off()


#EXTINCTION#
hist(tree_topcon_center$ExtinctTime) #right skew with lots of values at one location
max(tree_topcon$ExtinctTime)

#probability of extinction
tree_topcon_center$Extinct[tree_topcon_center$ExtinctTime == 501] <- 0 #survive = 0
tree_topcon_center$Extinct[tree_topcon_center$ExtinctTime < 501] <- 1 #extinct = 1

tree.topcon.ex.mod1 <- glm(Extinct ~ nDisturbed * Sev * Attract * Disperse, data = tree_topcon_center, family = binomial(link = "logit"))
summary(tree.topcon.ex.mod1)
write.csv(summary(tree.topcon.ex.mod1)$coefficients, file="tree_topcon_exprob_coef.csv")


tree.topcon.excoef <- as.data.frame(cbind(summary(tree.topcon.ex.mod1)$coefficients[,1], confint.default(tree.topcon.ex.mod1)))
tree.topcon.excoef <- tree.topcon.excoef[-1,]
colnames(tree.topcon.excoef) <- c("mean_Exprob", "lowerCI", "upperCI")
tree.topcon.excoef$var <- rownames(tree.topcon.excoef)
tree.topcon.excoef <- tree.topcon.excoef[order(-tree.topcon.excoef$mean),]
tree.topcon.excoef$var <- factor(tree.topcon.excoef$var, levels=tree.topcon.excoef$var) 
tree.topcon.excoef$Network <- "Tree"
tree.topcon.excoef$mSelect <- "Top"

min(tree.topcon.excoef$lowerCI)
max(tree.topcon.excoef$upperCI)

p.tree.topcon.exprob <- ggplot(tree.topcon.excoef, aes(x=var, y=mean_Exprob)) +
  geom_point() +
  geom_errorbar(tree.topcon.excoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree TopCon- Probability of Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('tree_topcon_exprob.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.tree.topcon.exprob
dev.off()


#Time to extinction
tree.topcon.ex.mod2 <- glm(ExtinctTime ~ nDisturbed * Sev * Attract * Disperse, data = tree_topcon_center)
summary(tree.topcon.ex.mod2)
write.csv(summary(tree.topcon.ex.mod2)$coefficients, file="tree_topcon_extime_coef.csv")

tree.topcon.excoef.mod2 <- as.data.frame(cbind(summary(tree.topcon.ex.mod2)$coefficients[,1], confint.default(tree.topcon.ex.mod2)))
tree.topcon.excoef.mod2 <- tree.topcon.excoef.mod2[-1,]
colnames(tree.topcon.excoef.mod2) <- c("mean_ExTime", "lowerCI", "upperCI")
tree.topcon.excoef.mod2$var <- rownames(tree.topcon.excoef.mod2)
tree.topcon.excoef.mod2 <- tree.topcon.excoef.mod2[order(-tree.topcon.excoef.mod2$mean),]
tree.topcon.excoef.mod2$var <- factor(tree.topcon.excoef.mod2$var, levels=tree.topcon.excoef.mod2$var) 
tree.topcon.excoef.mod2$Network <- "Tree"
tree.topcon.excoef.mod2$mSelect <- "Top"

p.extime.tree.topcon <- ggplot(tree.topcon.excoef.mod2, aes(x=var, y=mean_ExTime)) +
  geom_point() +
  geom_errorbar(tree.topcon.excoef.mod2, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree TopCon - Time to Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )

jpeg('tree_topcon_extime.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.extime.tree.topcon
dev.off()



###TREE NETWORK WITH Bottom CONNECTIVITY SELECTED###
tree_botcon <- read.csv(file="Tree_BotCon_LatinSQ.csv", header=T)
names(tree_botcon)

tree_botcon$nDisturbed <- (tree_botcon$nDisturbed/60) * 100

hist(tree_botcon$metaR) #left skew
hist(tree_botcon$ExtinctTime) #right skew with a lot of values at 500
hist(tree_botcon$Attract)

tree_botcon_center <- tree_botcon # new dataset for scaling

#METAPOPULATION GROWTH RATE#

#Scale all the variables to SD and variance
tree_botcon_center$metaR <- scale(tree_botcon_center$metaR) #mean centre and scale to SD
tree_botcon_center$nDisturbed <- scale(tree_botcon_center$nDisturbed)
tree_botcon_center$Sev <- scale(tree_botcon_center$Sev)
tree_botcon_center$Attract <- scale(tree_botcon_center$Attract)
tree_botcon_center$Disperse <- scale(tree_botcon_center$Disperse)
#hist(tree_botcon_center$nDisturbed)
#hist(tree_botcon_center$Sev)

hist(tree_botcon_center$metaR) #left skew with negative values
min(tree_botcon_center$metaR) #value to shift the data to the right to get rid of negative values

#mod1.meta <- glm((metaR + 6.513615) ~ nDisturbed * Sev * Attract * Disperse, family = Gamma(link="inverse"), data = full_rand_center)
tree.botcon.meta.mod1 <- glm(metaR ~ nDisturbed * Sev * Attract * Disperse, data = tree_botcon_center)
summary(tree.botcon.meta.mod1)
write.csv(summary(tree.botcon.meta.mod1)$coefficients, file="tree_botcon_meta_coef.csv")

tree.botcon.metacoef <- as.data.frame(cbind(summary(tree.botcon.meta.mod1)$coefficients[,1], confint.default(tree.botcon.meta.mod1)))
tree.botcon.metacoef <- tree.botcon.metacoef[-1,]
colnames(tree.botcon.metacoef) <- c("mean_MetaR", "lowerCI", "upperCI")
tree.botcon.metacoef$var <- rownames(tree.botcon.metacoef)
tree.botcon.metacoef <- tree.botcon.metacoef[order(-tree.botcon.metacoef$mean),]
tree.botcon.metacoef$var <- factor(tree.botcon.metacoef$var, levels=tree.botcon.metacoef$var) 
tree.botcon.metacoef$Network <- "Tree"
tree.botcon.metacoef$mSelect <- "Bottom"

min(tree.botcon.metacoef$lowerCI)
max(tree.botcon.metacoef$upperCI)

p.tree.meta.botcon <- ggplot(tree.botcon.metacoef, aes(x=var, y=mean_MetaR)) +
  geom_point() +
  geom_errorbar(tree.botcon.metacoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree Bot Con - Metapopulation Growth Rate") +
  #scale_y_continuous(breaks=c(-0.0002, 0, 0.0002, 0.0004), labels = c("-0.002", "0", "0.002", "0.004"), name = "Effect Size") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('tree_botconn_meta.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.tree.meta.botcon
dev.off()


#EXTINCTION#
hist(tree_botcon_center$ExtinctTime) #right skew with lots of values at one location
max(tree_botcon$ExtinctTime)

#probability of extinction
tree_botcon_center$Extinct[tree_botcon_center$ExtinctTime == 501] <- 0 #survive = 0
tree_botcon_center$Extinct[tree_botcon_center$ExtinctTime < 501] <- 1 #extinct = 1

tree.botcon.ex.mod1 <- glm(Extinct ~ nDisturbed * Sev * Attract * Disperse, data = tree_botcon_center, family = binomial(link = "logit"))
summary(tree.botcon.ex.mod1)
write.csv(summary(tree.botcon.ex.mod1)$coefficients, file="tree_botcon_exprob_coef.csv")


tree.botcon.excoef <- as.data.frame(cbind(summary(tree.botcon.ex.mod1)$coefficients[,1], confint.default(tree.botcon.ex.mod1)))
tree.botcon.excoef <- tree.botcon.excoef[-1,]
colnames(tree.botcon.excoef) <- c("mean_Exprob", "lowerCI", "upperCI")
tree.botcon.excoef$var <- rownames(tree.botcon.excoef)
tree.botcon.excoef <- tree.botcon.excoef[order(-tree.botcon.excoef$mean),]
tree.botcon.excoef$var <- factor(tree.botcon.excoef$var, levels=tree.botcon.excoef$var) 
tree.botcon.excoef$Network <- "Tree"
tree.botcon.excoef$mSelect <- "Bottom"

min(tree.botcon.excoef$lowerCI)
max(tree.botcon.excoef$upperCI)

p.tree.botcon.exprob <- ggplot(tree.botcon.excoef, aes(x=var, y=mean_Exprob)) +
  geom_point() +
  geom_errorbar(tree.botcon.excoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree BotCon- Probability of Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('tree_botcon_exprob.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.tree.botcon.exprob
dev.off()


#Time to extinction
tree.botcon.ex.mod2 <- glm(ExtinctTime ~ nDisturbed * Sev * Attract * Disperse, data = tree_botcon_center)
summary(tree.botcon.ex.mod2)
write.csv(summary(tree.botcon.ex.mod2)$coefficients, file="tree_botcon_extime_coef.csv")

tree.botcon.excoef.mod2 <- as.data.frame(cbind(summary(tree.botcon.ex.mod2)$coefficients[,1], confint.default(tree.botcon.ex.mod2)))
tree.botcon.excoef.mod2 <- tree.botcon.excoef.mod2[-1,]
colnames(tree.botcon.excoef.mod2) <- c("mean_ExTime", "lowerCI", "upperCI")
tree.botcon.excoef.mod2$var <- rownames(tree.botcon.excoef.mod2)
tree.botcon.excoef.mod2 <- tree.botcon.excoef.mod2[order(-tree.botcon.excoef.mod2$mean),]
tree.botcon.excoef.mod2$var <- factor(tree.botcon.excoef.mod2$var, levels=tree.botcon.excoef.mod2$var) 
tree.botcon.excoef.mod2$Network <- "Tree"
tree.botcon.excoef.mod2$mSelect <- "Bottom"

p.extime.tree.botcon <- ggplot(tree.botcon.excoef.mod2, aes(x=var, y=mean_ExTime)) +
  geom_point() +
  geom_errorbar(tree.botcon.excoef.mod2, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Tree BotCon - Time to Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )

jpeg('tree_botcon_extime.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.extime.tree.botcon
dev.off()

###FIGURE WITH ALL TREE NETWORKS

#Meta R
tree.metar <- rbind(tree.metacoef, tree.topcon.metacoef, tree.botcon.metacoef)
tree.metar$mSelect<- factor(tree.metar$mSelect, levels = c("Random", "Top", "Bottom"))
tree.metar$var <- gsub("nDisturbed", "Disturbed", tree.metar$var)
tree.metar$var <- gsub("Sev", "Severity", tree.metar$var)


metar.plot.tree <- ggplot(tree.metar, aes(x=reorder(var, -mean_MetaR), y=mean_MetaR, colour=mSelect, order=mSelect)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(tree.metar, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Mean Metapopulation Growth Rate") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Selection") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('MetaR_Tree.jpeg', width=10, height=8, units="cm", res=300)
metar.plot.tree
dev.off()

#Ex Prob
tree.exprob <- rbind(tree.excoef, tree.topcon.excoef, tree.botcon.excoef)
tree.exprob$mSelect <- factor(tree.exprob$mSelect, levels=c("Random", "Top", "Bottom"))
tree.exprob$var <- gsub("nDisturbed", "Disturbed", tree.exprob$var)
tree.exprob$var <- gsub("Sev", "Severity", tree.exprob$var)

exprob.plot.tree <- ggplot(tree.exprob, aes(x=reorder(var, mean_Exprob), y=mean_Exprob, colour=mSelect, order=mSelect)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(tree.exprob, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Extinction Probability") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Tree") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('ExProb_Tree.jpeg', width=10, height=8, units="cm", res=300)
exprob.plot.tree
dev.off()

#Time to Extinction
tree.extime <- rbind(tree.excoef.mod2, tree.topcon.excoef.mod2, tree.botcon.excoef.mod2)
tree.extime$mSelect <- factor(tree.extime$mSelect, levels=c("Random", "Top", "Bottom"))
tree.extime$var <- gsub("nDisturbed", "Disturbed", tree.extime$var)
tree.extime$var <- gsub("Sev", "Severity", tree.extime$var)

extime.plot.tree <- ggplot(tree.extime, aes(x=reorder(var, -mean_ExTime), y=mean_ExTime, colour=mSelect, order=mSelect)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(tree.extime, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Time to Extinction") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Selection") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('ExTime_Tree.jpeg', width=10, height=8, units="cm", res=300)
extime.plot.tree
dev.off()


#Writing all three figures together
tiff('Tree_ConSelection.tif', width=10, height=20, units="cm", res=300)
ggarrange(metar.plot.tree, exprob.plot.tree, extime.plot.tree, ncol = 1, nrow = 3, align="hv",
          labels="AUTO", label.x=0.95, label.y=0.9, font.label = list(size=7, face ="plain"),
          common.legend = TRUE, legend="bottom")
dev.off()


###ERDOS RENYI SELECTIONG PATCHES BASED ON CONNECTIVITY

###ERDOS-RENYI NETWORK WITH TOP CONNECTIVITY SELECTED###
erdos_topcon<- read.csv(file="Erdos_TopCon_LatinSQ.csv", header=T)
names(erdos_topcon)

erdos_topcon$nDisturbed <- (erdos_topcon$nDisturbed/60) * 100

hist(erdos_topcon$metaR) #left skew
hist(erdos_topcon$ExtinctTime) #right skew with a lot of values at 500
hist(erdos_topcon$Attract)

erdos_topcon_center <- erdos_topcon # new dataset for scaling

#METAPOPULATION GROWTH RATE#

#Scale all the variables to SD and variance
erdos_topcon_center$metaR <- scale(erdos_topcon_center$metaR) #mean centre and scale to SD
erdos_topcon_center$nDisturbed <- scale(erdos_topcon_center$nDisturbed)
erdos_topcon_center$Sev <- scale(erdos_topcon_center$Sev)
erdos_topcon_center$Attract <- scale(erdos_topcon_center$Attract)
erdos_topcon_center$Disperse <- scale(erdos_topcon_center$Disperse)
#hist(erdos_topcon_center$nDisturbed)
#hist(erdos_topcon_center$Sev)

hist(erdos_topcon_center$metaR) #left skew with negative values
min(erdos_topcon_center$metaR) #value to shift the data to the right to get rid of negative values

#mod1.meta <- glm((metaR + 6.513615) ~ nDisturbed * Sev * Attract * Disperse, family = Gamma(link="inverse"), data = full_rand_center)
erdos.topcon.meta.mod1 <- glm(metaR ~ nDisturbed * Sev * Attract * Disperse, data = erdos_topcon_center)
summary(erdos.topcon.meta.mod1)
write.csv(summary(erdos.topcon.meta.mod1)$coefficients, file="erdos_topcon_meta_coef.csv")

erdos.topcon.metacoef <- as.data.frame(cbind(summary(erdos.topcon.meta.mod1)$coefficients[,1], confint.default(erdos.topcon.meta.mod1)))
erdos.topcon.metacoef <- erdos.topcon.metacoef[-1,]
colnames(erdos.topcon.metacoef) <- c("mean_MetaR", "lowerCI", "upperCI")
erdos.topcon.metacoef$var <- rownames(erdos.topcon.metacoef)
erdos.topcon.metacoef <- erdos.topcon.metacoef[order(-erdos.topcon.metacoef$mean),]
erdos.topcon.metacoef$var <- factor(erdos.topcon.metacoef$var, levels=erdos.topcon.metacoef$var) 
erdos.topcon.metacoef$Network <- "Erdos"
erdos.topcon.metacoef$mSelect <- "Top"

min(erdos.topcon.metacoef$lowerCI)
max(erdos.topcon.metacoef$upperCI)

p.erdos.meta.topcon <- ggplot(erdos.topcon.metacoef, aes(x=var, y=mean_MetaR)) +
  geom_point() +
  geom_errorbar(erdos.topcon.metacoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos TopCon - Metapopulation Growth Rate") +
  #scale_y_continuous(breaks=c(-0.0002, 0, 0.0002, 0.0004), labels = c("-0.002", "0", "0.002", "0.004"), name = "Effect Size") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('erdos_topcon_meta.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.erdos.meta.topcon
dev.off()


#EXTINCTION#
hist(erdos_topcon_center$ExtinctTime) #right skew with lots of values at one location
max(erdos_topcon$ExtinctTime)

#probability of extinction
erdos_topcon_center$Extinct[erdos_topcon_center$ExtinctTime == 501] <- 0 #survive = 0
erdos_topcon_center$Extinct[erdos_topcon_center$ExtinctTime < 501] <- 1 #extinct = 1

erdos.topcon.ex.mod1 <- glm(Extinct ~ nDisturbed * Sev * Attract * Disperse, data = erdos_topcon_center, family = binomial(link = "logit"))
summary(erdos.topcon.ex.mod1)
write.csv(summary(erdos.topcon.ex.mod1)$coefficients, file="erdos_topcon_exprob_coef.csv")

erdos.topcon.excoef <- as.data.frame(cbind(summary(erdos.topcon.ex.mod1)$coefficients[,1], confint.default(erdos.topcon.ex.mod1)))
erdos.topcon.excoef <- erdos.topcon.excoef[-1,]
colnames(erdos.topcon.excoef) <- c("mean_Exprob", "lowerCI", "upperCI")
erdos.topcon.excoef$var <- rownames(erdos.topcon.excoef)
erdos.topcon.excoef <- erdos.topcon.excoef[order(-erdos.topcon.excoef$mean),]
erdos.topcon.excoef$var <- factor(erdos.topcon.excoef$var, levels=erdos.topcon.excoef$var) 
erdos.topcon.excoef$Network <- "Erdos"
erdos.topcon.excoef$mSelect <- "Top"

min(erdos.topcon.excoef$lowerCI)
max(erdos.topcon.excoef$upperCI)

p.erdos.topcon.exprob <- ggplot(erdos.topcon.excoef, aes(x=var, y=mean_Exprob)) +
  geom_point() +
  geom_errorbar(erdos.topcon.excoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos TopCon- Probability of Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('erdos_topcon_exprob.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.erdos.topcon.exprob
dev.off()


#Time to extinction
erdos.topcon.ex.mod2 <- glm(ExtinctTime ~ nDisturbed * Sev * Attract * Disperse, data = erdos_topcon_center)
summary(erdos.topcon.ex.mod2)
write.csv(summary(erdos.topcon.ex.mod2)$coefficients, file="erdos_topcon_extime_coef.csv")

erdos.topcon.excoef.mod2 <- as.data.frame(cbind(summary(erdos.topcon.ex.mod2)$coefficients[,1], confint.default(erdos.topcon.ex.mod2)))
erdos.topcon.excoef.mod2 <- erdos.topcon.excoef.mod2[-1,]
colnames(erdos.topcon.excoef.mod2) <- c("mean_ExTime", "lowerCI", "upperCI")
erdos.topcon.excoef.mod2$var <- rownames(erdos.topcon.excoef.mod2)
erdos.topcon.excoef.mod2 <- erdos.topcon.excoef.mod2[order(-erdos.topcon.excoef.mod2$mean),]
erdos.topcon.excoef.mod2$var <- factor(erdos.topcon.excoef.mod2$var, levels=erdos.topcon.excoef.mod2$var) 
erdos.topcon.excoef.mod2$Network <- "Erdos"
erdos.topcon.excoef.mod2$mSelect <- "Top"

p.extime.erdos.topcon <- ggplot(erdos.topcon.excoef.mod2, aes(x=var, y=mean_ExTime)) +
  geom_point() +
  geom_errorbar(erdos.topcon.excoef.mod2, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos TopCon - Time to Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )

jpeg('erdos_topcon_extime.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.extime.erdos.topcon
dev.off()



###ERDOS NETWORK WITH Bottom CONNECTIVITY SELECTED###
erdos_botcon <- read.csv(file="Erdos_BotCon_LatinSQ.csv", header=T)
names(erdos_botcon)

erdos_botcon$nDisturbed <- (erdos_botcon$nDisturbed/60) * 100

hist(erdos_botcon$metaR) #left skew
hist(erdos_botcon$ExtinctTime) #right skew with a lot of values at 500
hist(erdos_botcon$Attract)

erdos_botcon_center <- erdos_botcon # new dataset for scaling

#METAPOPULATION GROWTH RATE#

#Scale all the variables to SD and variance
erdos_botcon_center$metaR <- scale(erdos_botcon_center$metaR) #mean centre and scale to SD
erdos_botcon_center$nDisturbed <- scale(erdos_botcon_center$nDisturbed)
erdos_botcon_center$Sev <- scale(erdos_botcon_center$Sev)
erdos_botcon_center$Attract <- scale(erdos_botcon_center$Attract)
erdos_botcon_center$Disperse <- scale(erdos_botcon_center$Disperse)
#hist(erdos_botcon_center$nDisturbed)
#hist(erdos_botcon_center$Sev)

hist(erdos_botcon_center$metaR) #left skew with negative values
min(erdos_botcon_center$metaR) #value to shift the data to the right to get rid of negative values

#mod1.meta <- glm((metaR + 6.513615) ~ nDisturbed * Sev * Attract * Disperse, family = Gamma(link="inverse"), data = full_rand_center)
erdos.botcon.meta.mod1 <- glm(metaR ~ nDisturbed * Sev * Attract * Disperse, data = erdos_botcon_center)
summary(erdos.botcon.meta.mod1)
write.csv(summary(erdos.botcon.meta.mod1)$coefficients, file="erdos_botcon_meta_coef.csv")

erdos.botcon.metacoef <- as.data.frame(cbind(summary(erdos.botcon.meta.mod1)$coefficients[,1], confint.default(erdos.botcon.meta.mod1)))
erdos.botcon.metacoef <- erdos.botcon.metacoef[-1,]
colnames(erdos.botcon.metacoef) <- c("mean_MetaR", "lowerCI", "upperCI")
erdos.botcon.metacoef$var <- rownames(erdos.botcon.metacoef)
erdos.botcon.metacoef <- erdos.botcon.metacoef[order(-erdos.botcon.metacoef$mean),]
erdos.botcon.metacoef$var <- factor(erdos.botcon.metacoef$var, levels=erdos.botcon.metacoef$var) 
erdos.botcon.metacoef$Network <- "Erdos"
erdos.botcon.metacoef$mSelect <- "Bottom"

min(erdos.botcon.metacoef$lowerCI)
max(erdos.botcon.metacoef$upperCI)

p.erdos.meta.botcon <- ggplot(erdos.botcon.metacoef, aes(x=var, y=mean_MetaR)) +
  geom_point() +
  geom_errorbar(erdos.botcon.metacoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos BotCon - Metapopulation Growth Rate") +
  #scale_y_continuous(breaks=c(-0.0002, 0, 0.0002, 0.0004), labels = c("-0.002", "0", "0.002", "0.004"), name = "Effect Size") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('erdos_botconn_meta.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.erdos.meta.botcon
dev.off()


#EXTINCTION#
hist(erdos_botcon_center$ExtinctTime) #right skew with lots of values at one location
max(erdos_botcon$ExtinctTime)

#probability of extinction
erdos_botcon_center$Extinct[erdos_botcon_center$ExtinctTime == 501] <- 0 #survive = 0
erdos_botcon_center$Extinct[erdos_botcon_center$ExtinctTime < 501] <- 1 #extinct = 1

erdos.botcon.ex.mod1 <- glm(Extinct ~ nDisturbed * Sev * Attract * Disperse, data = erdos_botcon_center, family = binomial(link = "logit"))
summary(erdos.botcon.ex.mod1)
write.csv(summary(erdos.botcon.ex.mod1)$coefficients, file="erdos_botcon_exprob_coef.csv")


erdos.botcon.excoef <- as.data.frame(cbind(summary(erdos.botcon.ex.mod1)$coefficients[,1], confint.default(erdos.botcon.ex.mod1)))
erdos.botcon.excoef <- erdos.botcon.excoef[-1,]
colnames(erdos.botcon.excoef) <- c("mean_Exprob", "lowerCI", "upperCI")
erdos.botcon.excoef$var <- rownames(erdos.botcon.excoef)
erdos.botcon.excoef <- erdos.botcon.excoef[order(-erdos.botcon.excoef$mean),]
erdos.botcon.excoef$var <- factor(erdos.botcon.excoef$var, levels=erdos.botcon.excoef$var) 
erdos.botcon.excoef$Network <- "Erdos"
erdos.botcon.excoef$mSelect <- "Bottom"

min(erdos.botcon.excoef$lowerCI)
max(erdos.botcon.excoef$upperCI)

p.erdos.botcon.exprob <- ggplot(erdos.botcon.excoef, aes(x=var, y=mean_Exprob)) +
  geom_point() +
  geom_errorbar(erdos.botcon.excoef, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos BotCon- Probability of Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )


jpeg('erdos_botcon_exprob.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.erdos.botcon.exprob
dev.off()


#Time to extinction
erdos.botcon.ex.mod2 <- glm(ExtinctTime ~ nDisturbed * Sev * Attract * Disperse, data = erdos_botcon_center)
summary(erdos.botcon.ex.mod2)
write.csv(summary(erdos.botcon.ex.mod2)$coefficients, file="erdos_botcon_extime_coef.csv")

erdos.botcon.excoef.mod2 <- as.data.frame(cbind(summary(erdos.botcon.ex.mod2)$coefficients[,1], confint.default(erdos.botcon.ex.mod2)))
erdos.botcon.excoef.mod2 <- erdos.botcon.excoef.mod2[-1,]
colnames(erdos.botcon.excoef.mod2) <- c("mean_ExTime", "lowerCI", "upperCI")
erdos.botcon.excoef.mod2$var <- rownames(erdos.botcon.excoef.mod2)
erdos.botcon.excoef.mod2 <- erdos.botcon.excoef.mod2[order(-erdos.botcon.excoef.mod2$mean),]
erdos.botcon.excoef.mod2$var <- factor(erdos.botcon.excoef.mod2$var, levels=erdos.botcon.excoef.mod2$var) 
erdos.botcon.excoef.mod2$Network <- "Erdos"
erdos.botcon.excoef.mod2$mSelect <- "Bottom"

p.extime.erdos.botcon <- ggplot(erdos.botcon.excoef.mod2, aes(x=var, y=mean_ExTime)) +
  geom_point() +
  geom_errorbar(erdos.botcon.excoef.mod2, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5) +
  geom_hline(yintercept=0) +
  labs(title = "Erdos BotCon - Time to Extinction") +
  #scale_y_continuous(breaks=c(-0.3, -0.15, 0, 0.15), labels = c("-3", "-0.15", "0", "0.15"), name = "Effect Size") +
  scale_y_continuous(name = "Effect Size") +
  scale_x_discrete(name = "Factor") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6)
  )

jpeg('erdos_botcon_extime.jpeg', width=7.5, height=7.5, units="cm", res=300)
p.extime.erdos.botcon
dev.off()

###FIGURE WITH ALL Erdos NETWORKS

#Meta R
erdos.metar <- rbind(erdos.metacoef, erdos.topcon.metacoef, erdos.botcon.metacoef)
erdos.metar$mSelect <- factor(erdos.metar$mSelect, levels = c("Random", "Top", "Bottom"))
erdos.metar$var <- gsub("nDisturbed", "Disturbed", erdos.metar$var)
erdos.metar$var <- gsub("Sev", "Severity", erdos.metar$var)

metar.plot.erdos <- ggplot(erdos.metar, aes(x=reorder(var, -mean_MetaR), y=mean_MetaR, colour=mSelect, order=mSelect)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(erdos.metar, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Mean Metapopulation Growth Rate") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Selection") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('MetaR_Erdos.jpeg', width=10, height=8, units="cm", res=300)
metar.plot.erdos
dev.off()

#Ex Prob
erdos.exprob <- rbind(erdos.excoef, erdos.topcon.excoef, erdos.botcon.excoef)
erdos.exprob$mSelect <- factor(erdos.exprob$mSelect, levels = c("Random", "Top", "Bottom"))
erdos.exprob$var <- gsub("nDisturbed", "Disturbed", erdos.exprob$var)
erdos.exprob$var <- gsub("Sev", "Severity", erdos.exprob$var)

exprob.plot.erdos <- ggplot(erdos.exprob, aes(x=reorder(var, mean_Exprob), y=mean_Exprob, colour=mSelect, order=mSelect)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(erdos.exprob, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Extinction Probability") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Selection") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('ExProb_Erdos.jpeg', width=10, height=8, units="cm", res=300)
exprob.plot.erdos
dev.off()

#Time to Extinction
erdos.extime <- rbind(erdos.excoef.mod2, erdos.topcon.excoef.mod2, erdos.botcon.excoef.mod2)
erdos.extime$mSelect <- factor(erdos.extime$mSelect, levels = c("Random", "Top", "Bottom"))
erdos.extime$var <- gsub("nDisturbed", "Disturbed", erdos.extime$var)
erdos.extime$var <- gsub("Sev", "Severity", erdos.extime$var)

extime.plot.erdos <- ggplot(erdos.extime, aes(x=reorder(var, -mean_ExTime), y=mean_ExTime, colour=mSelect, order=mSelect)) +
  geom_point(size = 0.3, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_errorbar(erdos.extime, mapping=aes(ymin=lowerCI, ymax=upperCI), width=0.5, position=position_dodge2(width=0.5, reverse = TRUE)) +
  geom_hline(yintercept=0) +
  labs(title = "Time to Extinction") +
  scale_y_continuous(name= "Effect Size") +
  scale_x_discrete(name = "Factor") +
  labs(col="Selection") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size=7, hjust=-0),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6),
    legend.position = c(0.2, 0.25),
    legend.background = element_rect(size=0.5, linetype="solid", colour ="black")
  )

jpeg('ExTime_Erdos.jpeg', width=10, height=8, units="cm", res=300)
extime.plot.erdos
dev.off()

tiff('Erdos_ConSelection.tif', width=10, height=20, units="cm", res=300)
ggarrange(metar.plot.erdos, exprob.plot.erdos, extime.plot.erdos, ncol = 1, nrow = 3, align="hv",
          labels="AUTO", label.x=0.95, label.y=0.9, font.label = list(size=7, face ="plain"),
          common.legend = TRUE, legend="bottom")
dev.off()


