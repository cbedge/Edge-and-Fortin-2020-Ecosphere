###Analyses of individual effects
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(plotly)
library(plot3D)
#library(RColorBrewer)
library(ggpubr)

setwd("N:/Population Model Results") #Set directory to store outputs
setwd("Z:/cedge/Population Model Results") #Laptop

###INVESTIGATE THE EFFECTS OF DISPERSAL ON ALL THREE RESPONSE VARIABLES
###ALL THE NETWORKS AT THE SAME TIME

#Effect of dispersal
dispersal <- read.csv("Dispersal.csv", header=T)
names(dispersal)
dispersal$Extinct[dispersal$ExtinctTime == 501] <- 0 #survive = 0
dispersal$Extinct[dispersal$ExtinctTime < 501] <- 1 #extinct = 1
dispersal$Network <- factor(dispersal$Network, levels = c("Full", "Erdos", "Tree"))

#plot of dispersal V metaR
DispVR_P <- ggplot(dispersal, aes(x = Disperse, y = metaR, factor = Network)) +
  geom_point(size = 1, aes(color = factor(Network)), alpha = 1/10) +
  geom_smooth(se=T, method = "glm", aes(color = factor(Network))) +
  scale_y_continuous(name= "Mean metapopulation growth rate") +
  scale_x_continuous(name = "Dispersal") +
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
    legend.title=element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size=7),
    legend.key = element_rect(fill = "white")
  )

#plot of dispersal V Exprob
DispVExProb_P <- ggplot(dispersal, aes(x = Disperse, y = Extinct, factor = Network)) +
  geom_point(size = 1, aes(color = factor(Network)), alpha = 1/10) +
  geom_smooth(se=T, method = "glm", aes(color = factor(Network))) +
  scale_y_continuous(breaks=c(0,1), name= "Extinction Probabilty") +
  scale_x_continuous(name = "Dispersal") +
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
    legend.title=element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size=7),
    legend.key = element_rect(fill = "white")
  )

#plot of dispersal V Extime
DispVTime_P <- ggplot(dispersal, aes(x = Disperse, y = ExtinctTime, factor = Network)) +
  geom_point(size = 1, aes(color = factor(Network)), alpha = 1/10) +
  geom_smooth(se=T, method = "glm", aes(color = factor(Network))) +
  scale_y_continuous(name= "Time to extinction") +
  scale_x_continuous(name = "Dispersal") +
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
    legend.title=element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size=7),
    legend.key = element_rect(fill = "white")
  )

tiff('Dispersal.tif', width=8, height=18, units="cm", res=300)
ggarrange(DispVR_P, DispVExProb_P, DispVTime_P, ncol = 1, nrow = 3, align="hv",
          labels="AUTO", label.x=0.93, label.y=0.95, font.label = list(size=7, face ="plain"),
          common.legend = TRUE, legend="bottom")
dev.off()


###USING 25% DISPERSAL INVESTIGATE THE INTERACITON BETWEEN NUMBER OF DISTURBACNES AND SEVERITY

##ERDOS
erdos_distxsev <- read.csv("Erdos_ndistxsev.csv", header=T)

erdos_distxsev$Extinct[erdos_distxsev$ExtinctTime == 501] <- 1 #survive = 1
erdos_distxsev$Extinct[erdos_distxsev$ExtinctTime < 501] <- 0 #extinct = 0

erdos_distxsev_center <- erdos_distxsev

erdos_distxsev_center$metaR <- as.numeric(scale(erdos_distxsev_center$metaR)) #mean centre and scale to SD
erdos_distxsev_center$nDisturbed <- as.numeric(scale(erdos_distxsev_center$nDisturbed))
erdos_distxsev_center$Sev <- as.numeric(scale(erdos_distxsev_center$Sev))
erdos_distxsev_center$Attract <- as.numeric(scale(erdos_distxsev_center$Attract))
erdos_distxsev_center$Disperse <- as.numeric(scale(erdos_distxsev_center$Disperse))

E_mod1 <- glm(metaR ~ nDisturbed * Sev, data = erdos_distxsev_center)
summary(E_mod1)

E_mod2 <- glm(Extinct ~ nDisturbed * Sev, data = erdos_distxsev_center, family = binomial(link = "logit"))
summary(E_mod2)

E_mod3 <- glm(ExtinctTime ~ nDisturbed * Sev, data = erdos_distxsev_center)
summary(E_mod3)


grid.lines = 25
E.x.pred <- seq(min(erdos_distxsev_center$nDisturbed), max(erdos_distxsev_center$nDisturbed), length.out = grid.lines)
E.y.pred <- seq(min(erdos_distxsev_center$Sev), max(erdos_distxsev_center$Sev), length.out = grid.lines)
E.xy <- expand.grid(nDisturbed = E.x.pred, Sev = E.y.pred)
E.z.pred <- matrix(predict(E_mod1, newdata = E.xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
E.fitpoints <- predict(E_mod1)

scatter3D(x = erdos_distxsev_center$nDisturbed, y = erdos_distxsev_center$Sev, z = erdos_distxsev_center$metaR, 
          main = "Erdos", xlab = "Disturbed",
          ylab ="Severity", zlab = "Mean R",
          col = ramp.col(c("red", "yellow", "blue")), pch = 19, cex = 0.5, 
          phi = 10, theta = 60,
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))



##FULL
full_distxsev <- read.csv("Full_ndistxsev.csv", header=T)

full_distxsev$Extinct[full_distxsev$ExtinctTime == 501] <- 1 #survive = 1
full_distxsev$Extinct[full_distxsev$ExtinctTime < 501] <- 0 #extinct = 0

full_distxsev_center <- full_distxsev

full_distxsev_center$metaR <- as.numeric(scale(full_distxsev_center$metaR)) #mean centre and scale to SD
full_distxsev_center$nDisturbed <- as.numeric(scale(full_distxsev_center$nDisturbed))
full_distxsev_center$Sev <- as.numeric(scale(full_distxsev_center$Sev))
full_distxsev_center$Attract <- as.numeric(scale(full_distxsev_center$Attract))
full_distxsev_center$Disperse <- as.numeric(scale(full_distxsev_center$Disperse))

F_mod1 <- glm(metaR ~ nDisturbed * Sev, data = full_distxsev_center)
summary(F_mod1)

F_mod2 <- glm(Extinct ~ nDisturbed * Sev, data = full_distxsev_center, family = binomial(link = "logit"))
summary(F_mod2)

F_mod3 <- glm(ExtinctTime ~ nDisturbed * Sev, data = full_distxsev_center)
summary(F_mod3)


grid.lines = 25
F.x.pred <- seq(min(full_distxsev_center$nDisturbed), max(full_distxsev_center$nDisturbed), length.out = grid.lines)
F.y.pred <- seq(min(full_distxsev_center$Sev), max(full_distxsev_center$Sev), length.out = grid.lines)
F.xy <- expand.grid(nDisturbed = F.x.pred, Sev = F.y.pred)
F.z.pred <- matrix(predict(F_mod1, newdata = F.xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
F.fitpoints <- predict(F_mod1)

scatter3D(x = full_distxsev_center$nDisturbed, y = full_distxsev_center$Sev, z = full_distxsev_center$metaR, 
          main = "Full", xlab = "Disturbed",
          ylab ="Severity", zlab = "Mean R",
          col = ramp.col(c("red", "yellow", "blue")), pch = 19, cex = 0.5, 
          phi = 10, theta = 60,
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))



##Tree
tree_distxsev <- read.csv("Tree_ndistxsev.csv", header=T)

tree_distxsev$Extinct[tree_distxsev$ExtinctTime == 501] <- 1 #survive = 1
tree_distxsev$Extinct[tree_distxsev$ExtinctTime < 501] <- 0 #extinct = 0

tree_distxsev_center <- tree_distxsev

tree_distxsev_center$metaR <- as.numeric(scale(tree_distxsev_center$metaR)) #mean centre and scale to SD
tree_distxsev_center$nDisturbed <- as.numeric(scale(tree_distxsev_center$nDisturbed))
tree_distxsev_center$Sev <- as.numeric(scale(tree_distxsev_center$Sev))
tree_distxsev_center$Attract <- as.numeric(scale(tree_distxsev_center$Attract))
tree_distxsev_center$Disperse <- as.numeric(scale(tree_distxsev_center$Disperse))

T_mod1 <- glm(metaR ~ nDisturbed * Sev, data = tree_distxsev_center)
summary(T_mod1)

T_mod2 <- glm(Extinct ~ nDisturbed * Sev, data = tree_distxsev_center, family = binomial(link = "logit"))
summary(T_mod2)

T_mod3 <- glm(ExtinctTime ~ nDisturbed * Sev, data = tree_distxsev_center)
summary(T_mod3)


grid.lines = 25
T.x.pred <- seq(min(tree_distxsev_center$nDisturbed), max(tree_distxsev_center$nDisturbed), length.out = grid.lines)
T.y.pred <- seq(min(tree_distxsev_center$Sev), max(tree_distxsev_center$Sev), length.out = grid.lines)
T.xy <- expand.grid(nDisturbed = T.x.pred, Sev = T.y.pred)
T.z.pred <- matrix(predict(T_mod1, newdata = T.xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
T.fitpoints <- predict(T_mod1)

scatter3D(x = tree_distxsev_center$nDisturbed, y = tree_distxsev_center$Sev, z = tree_distxsev_center$metaR, 
          main = "Tree", xlab = "Disturbed",
          ylab ="Severity", zlab = "Mean R",
          col = ramp.col(c("red", "yellow", "blue")), pch = 19, cex = 0.5, 
          phi = 10, theta = 60,
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA))


plot_ly(showscale = FALSE) %>%
  add_markers(tree_distxsev_center, x = ~nDisturbed,
              y = ~Sev,
              z = ~metaR)

plot_ly(showscale = FALSE) %>%
  add_surface(z = ~F.z.pred, opacity = 0.5) %>%
  add_surface(z = ~E.z.pred, opacity = 0.5) %>%
  add_surface(z = ~T.z.pred, opacity = 0.5)




################################################################


colvar = NULL, col = "#0072B2", 
, colvar=NULL, col="#D55E00"

persp(x.pred, y.pred, z.pred, col="lightblue", phi = 10, theta = 60)


distsev_p <- plot_ly(erdos_distxsev, x = ~nDisturbed, y = ~Sev, z = ~metaR) %>%
  add_markers()
