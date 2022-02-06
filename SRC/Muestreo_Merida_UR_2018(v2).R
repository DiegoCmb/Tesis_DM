
#R Version: 4.0.3

####### PACKAGES AND DATABASES ########

########### Loading all packages needed
library(tidyverse)
library(car)
library(agricolae)
library(grid)
library(gridExtra)
library(lmtest)
library(lme4) #lmer
library(magrittr) #estimadores marginales de prob.
library(ggeffects)
library(sjmisc)
library(nlme)
library(gstat)  #coord
library(sp)    #coord
library(ade4)  #coord
library(vegan) #shannon
library(olsrr)    #
library(lmerTest) #ranova
library(MuMIn) #r.squaredGLMM
library(GGally) #for CCA
library(CCA)  #for CCA
library(effsize) #for effect size
library(factoextra) # para get_pca_var y fviz_pca_biplot

############# Loaded the databases
Muestreo_UR <- read.csv("Database/Muestreo_UR_Merida_2018.csv")[,-1]
#Sum.Muestreo_UR <- read.csv("Database/Muestreo_Summary.csv")[-1]
Sum.Muestreo_UR1 <- read.csv("Database/Muestreo_Summary1.csv")[-1]


# Changing columns to factors
Muestreo_UR$site <-as.factor(Muestreo_UR$site)
Muestreo_UR$area <-as.factor(Muestreo_UR$area)
Muestreo_UR$ambient <-as.factor(Muestreo_UR$ambient)

Muestreo_UR<-Muestreo_UR[Muestreo_UR$area!="open",]


############ FUNCTIONS ###########

##Function for Standard Errors
se<-function(x, na.rm=F){
  if(na.rm==F){
    se=sd(x)/sqrt(length(x))
    return(se)
  }
  else{
    x=x[!is.na(x)]
    se=sd(x)/sqrt(length(x))
    return(se)
  }
}


##Corralations between soil nutrient and characteristics

#extract R-squared
cor.r <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X)
  above <- row(R) < col(R)
  R[upper.tri(R)] <- NA
  R[diag(R)] <- NA
  R
}

#as matrix
cor.r.m <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X)
  above <- row(R) < col(R)
  R
}

#extract P-values
cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X)
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  
  cor.mat <- t(R)
  cor.mat[upper.tri(cor.mat)] <- NA
  cor.mat
}

#as matrix
cor.prob.m <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X)
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  
  cor.mat <- t(R)
  cor.mat[above] <- 1 - pf(Fstat, 1, dfr)
  cor.mat
}


################################################################################

############################################ 
##############                ############## 
##############   GRAFICS AND  ############## 
##############    ANALYSES    ##############
##############                ############## 
############################################ 



######################################
#######    Nutrient Analyses   #######
######################################

##ESTOS SON PARTE DE LOS RESULTADOS
#ESTOS ANALISIS SE MUESTRAN EN TABLA 1


####### PH * AREA #######

#Linear regresion
PHSOILa <- lm(pH~area, Sum.Muestreo_UR1)
summary(PHSOILa)
anova(PHSOILa)

hist(Sum.Muestreo_UR1$pH)
plot(PHSOILa,2)#the data is normaly distributed
ols_test_normality(PHSOILa)#TESTING RESIDUALS


SP<-summary(lm(pH~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))
#Graphic
pHarea<-ggplot(subset(Sum.Muestreo_UR1, !is.na(pH)), 
       aes(x=area, y=pH, fill=area))+
  geom_boxplot() + theme_classic()+
  labs(title="A) pH del suelo", x="Area") +
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural", 
                            "sidew" = "Profundo"))+
  geom_text(data=data.frame(Sum.Muestreo_UR1 %>% 
                              group_by(area) %>% 
                              summarize(Max.pH=max(pH, na.rm = T))),
            aes(x=area,y=0.05+Max.pH,
                label=c("a","b")),vjust=0)


pHarea
##ggpredict
fit_PHSOILa<-ggpredict(PHSOILa, terms=c("area"))
plot(fit_PHSOILa)




####### K * AREA #######

#Linear regresion
KSOILa <- lm(log(K)~area, Sum.Muestreo_UR1)
summary(KSOILa)
anova(KSOILa)

hist(log(Sum.Muestreo_UR1$K))
plot(KSOILa,2)#the data is NOT normaly distributed
ols_test_normality(KSOILa)#TESTING RESIDUALS


KSOILa0 <- glm(K~1, Sum.Muestreo_UR1,
              family = Gamma)
KSOILa1 <- glm(K~area, Sum.Muestreo_UR1,
              family = Gamma)

anova(KSOILa0, KSOILa1)
summary(KSOILa1)
Anova(KSOILa1)


SP<-summary(lm(K~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))

#Graphic
Karea<-ggplot(subset(Sum.Muestreo_UR1, !is.na(K)), 
       aes(x=area, y=K, fill=area))+
  geom_boxplot() + theme_classic()+
  labs(title="B) Concentración de Potasio",x="Area", y="K (mg/kg)")+
  geom_text(data=data.frame(Sum.Muestreo_UR1 %>% 
                              group_by(area) %>% 
                              summarize(Max.K=max(K, na.rm = T))),
            aes(x=area,y=0.02+Max.K,
                label=c("a","b")),vjust=0)+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural", 
                            "sidew" = "Profundo"))



####### P * AREA #######

#Linear regresion
PSOILa <- lm(log(P)~area, Sum.Muestreo_UR1)
summary(PSOILa)
anova(PSOILa)

hist(log(Sum.Muestreo_UR1$P))
plot(PSOILa,2)#the data is NOT normaly distributed
ols_test_normality(PSOILa)#TESTING RESIDUALS


PSOILa0 <- glm(P~1, Sum.Muestreo_UR1,
               family = Gamma)
PSOILa1 <- glm(P~area, Sum.Muestreo_UR1,
               family = Gamma)



anova(PSOILa0, PSOILa1)
Anova(PSOILa1)

SP<-summary(lm(P~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))




#Graphic
Parea<-ggplot(subset(Sum.Muestreo_UR1, !is.na(P)), 
       aes(x=area, y=P, fill=area))+
  geom_boxplot() + theme_classic()+
  labs(title="C) Concentración de Fósforo", 
       x="Area", y="P(mg/kg)")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural", 
                            "sidew" = "Profundo"))+
  geom_text(data=data.frame(Sum.Muestreo_UR1 %>% 
                              group_by(area) %>% 
                              summarize(Max.P=max(P, na.rm = T))),
            aes(x=area,y=5+Max.P,
                label=c("b","a")),vjust=0)


####### N * AREA #######

#Linear regresion
NSOILa <- lm(log(N)~area, Sum.Muestreo_UR1)
summary(NSOILa)
anova(NSOILa)

plot(NSOILa,2)#the data is NOT normaly distributed
ols_test_normality(NSOILa)#TESTING RESIDUALS


NSOILa0 <- glm(N~1, Sum.Muestreo_UR1,
               family = Gamma)
NSOILa1 <- glm(N~area, Sum.Muestreo_UR1,
               family = Gamma)

anova(NSOILa0, NSOILa1)
Anova(NSOILa1)
summary(NSOILa1)


SP<-summary(lm(N~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))

#Graphic
Narea<-ggplot(subset(Sum.Muestreo_UR1, !is.na(N)), 
       aes(x=area, y=N, fill=area))+
  geom_boxplot() + theme_classic()+
  labs(title="D) Concentración de Nitrógeno", x="Area", y="N(%)")+
  guides(fill=F)+
  geom_text(data=data.frame(Sum.Muestreo_UR1 %>% 
                              group_by(area) %>% 
                              summarize(Max.N=max(N, na.rm = T))),
            aes(x=area,y=0.05+Max.N,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural", 
                            "sidew" = "Profundo"))


grid.arrange(pHarea, Karea, Parea, Narea,
                   widths= c(1,1,1,1), 
                   layout_matrix= rbind(c(4,3,2,1)))



####### Ca * AREA #######

#Linear Model
CaSOILa <- lm(log(Ca)~area, Sum.Muestreo_UR1)
summary(CaSOILa)
anova(CaSOILa)

hist(log(Sum.Muestreo_UR1$Ca))
plot(CaSOILa,2)#the data is NOT normaly distributed
ols_test_normality(CaSOILa)#TESTING RESIDUALS


SP<-summary(lm(Ca~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))


#Generalized Linear Model
CaSOILa0 <- glm(Ca~1, Sum.Muestreo_UR1,
               family = Gamma)
CaSOILa1 <- glm(Ca~area, Sum.Muestreo_UR1,
               family = Gamma)

anova(CaSOILa0, CaSOILa1)
Anova(CaSOILa1)
summary(CaSOILa1)



####### Na * AREA #######


#Linear Model
NaSOILa <- lm(log(Na)~area, Sum.Muestreo_UR1)
summary(NaSOILa)
anova(NaSOILa)

hist(log(Sum.Muestreo_UR1$Na))
plot(NaSOILa,2)#the data is NOT normaly distributed
ols_test_normality(NaSOILa)#TESTING RESIDUALS


#Generalized Linear Model
NaSOILa0 <- glm(Na~1, Sum.Muestreo_UR1,
               family = Gamma)
NaSOILa1 <- glm(Na~area, Sum.Muestreo_UR1,
               family = Gamma)

anova(NaSOILa0, NaSOILa1)
Anova(NaSOILa1)
summary(NaSOILa1)



SP<-summary(lm(Na~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))




#Extract means, sd & se
Sum.Muestreo_UR1 %>% 
  group_by(area) %>% 
  summarise(
    pH_mean=mean(pH, na.rm = T),
    pH_sd=sd(pH, na.rm = T),
    pH_se=se(pH, na.rm = T),
    K_mean=mean(K, na.rm = T),
    K_sd=sd(K, na.rm = T),
    K_se=se(K, na.rm = T),
    P_mean=mean(P, na.rm = T),
    P_sd=sd(P, na.rm = T),
    P_se=se(P, na.rm = T),
    N_mean=mean(N, na.rm = T),
    N_sd=sd(N, na.rm = T),
    N_se=se(N, na.rm = T),
    
    Ca_mean=mean(Ca, na.rm = T),
    Ca_sd=sd(Ca, na.rm = T),
    Ca_se=se(Ca, na.rm = T),
    Na_mean=mean(Na, na.rm = T),
    Na_sd=sd(Na, na.rm = T),
    Na_se=se(Na, na.rm = T)
  ) -> mean_suelos

(mean_suelos<-as.data.frame(mean_suelos))



#effect size
#cohen.d(?~area, data=Sum.Muestreo_UR1, na.action = na.omit, hedges.correction=T)
#No se uso en los resultados

Sum.Muestreo_UR1 %>%  
  summarise(
    pH_effsize=cohen.d(pH~area, data=Sum.Muestreo_UR1, 
                       na.action = na.omit, 
                       hedges.correction=T)[[4]],
    
    K_effsize=cohen.d(K~area, data=Sum.Muestreo_UR1, 
                      na.action = na.omit, 
                      hedges.correction=T)[[4]],
    
    P_effsize=cohen.d(P~area, data=Sum.Muestreo_UR1, 
                      na.action = na.omit, 
                      hedges.correction=T)[[4]],
    
    N_effsize=cohen.d(N~area, data=Sum.Muestreo_UR1, 
                      na.action = na.omit, 
                      hedges.correction=T)[[4]],
    
    Ca_effsize=cohen.d(Ca~area, data=Sum.Muestreo_UR1, 
                       na.action = na.omit, 
                       hedges.correction=T)[[4]],
    
    Na_effsize=cohen.d(Na~area, data=Sum.Muestreo_UR1, 
                       na.action = na.omit, 
                       hedges.correction=T)[[4]],
  ) -> effsize_suelos

(effsize_suelos<-as.data.frame(effsize_suelos))



###################################################
#### CORRELATIONS BETWEEN SOIL CHARACTERISTICS ####
###################################################

##SE MUESTRAN EN RESULTADOS, TABLA SUPLEMENTARIA

cor.r(Sum.Muestreo_UR1[,c(25,26,27,30)])#R-squared
cor.prob(Sum.Muestreo_UR1[,c(25,26,27,30)])#P-values

cor.r(Sum.Muestreo_UR1[,c(25,26,27,30,28,29,24)])#R-squared
cor.prob(Sum.Muestreo_UR1[,c(25,26,27,30,28,29,24)])#P-values


hmap_cor<-cor(Sum.Muestreo_UR1[,c(25,26,27,30,28,29,24)])
rownames(hmap_cor)[7]<-"Col"
colnames(hmap_cor)[7]<-"Col"

heatmap(hmap_cor, scale="column", 
        col=RColorBrewer::colorRampPalette(brewer.pal(5, "RdBu"))(256))

RColorBrewer::colorRampPalette(brewer.pal(5, "RdBu"))(2)

#HEATMAP FOR CORRELATIONS
hmap_cor[upper.tri(hmap_cor)]<-NA
hmap_cor<-round(hmap_cor, 2)
melt_corr<-reshape2::melt(hmap_cor, na.rm=TRUE)


corr_heatmap<-ggplot(melt_corr, aes(x=Var1, y=Var2, fill=value))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="#0571B0",high = "#CA0020", mid="white",
                       midpoint=0, limit=c(-1,1), space="Lab",
                       name="Pearson\nCorrelation")+
  theme_minimal()+
  coord_fixed()

corr_heatmap+
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1.1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))



###################################################
####         PCA SOIL CHARACTERISTICS          ####
###################################################

## SE MUESTRAN EN RESULTADOS

subset(Sum.Muestreo_UR1, area=="sidew")

soilpca<-Sum.Muestreo_UR1[,c(25:30)]
soil_spore_m<-Sum.Muestreo_UR1[,c(25:30,42,44,24)]

row.names(soilpca)<-Sum.Muestreo_UR1$site
# Estas transformaciones son las mismas que se usaron en el Structural eq modelling
# en particular se trabaja una re escalarización de P y Ca para que las varianzas 
# sean en ordenes de magnitud similares a las de las otras variables. 
# como PCA tambien requiere normalidad se hicieron las transformaciones necesarias
# y aplicadas en Str. eq. Modl. 

soilpca<-soilpca%>%transmute(
                 t.pH = scale(log(pH)),
                 t.N = scale(log(N)),
                 t.K = scale(log(K)),
                 t.Na = scale(log(Na)),
                 t.Ca = scale(log(Ca/100)),
                 t.P = scale(log(P/100)))

PCAsoil <- prcomp(soilpca, scale = TRUE)
summary(PCAsoil)
PCAsoil$rotation # loadings. <---- MUY IMPORTANTE REPORTAR. 
PCAsoil$center # media 
PCAsoil$sd # desv estandar 

PCAsoil1<-PCAsoil[["x"]][,1]
PCAsoil2<-PCAsoil[["x"]][,2]


Sum.Muestreo_UR2<-cbind(Sum.Muestreo_UR1, PCAsoil1, PCAsoil2)

Sum.Muestreo_UR2<-Sum.Muestreo_UR2%>%mutate( # incorporando variables transformadas
       t.pH = scale(log(pH)),
       t.N = scale(log(N)),
       t.K = scale(log(K)),
       t.Na = scale(log(Na)),
       t.Ca = scale(log(Ca/100)),
       t.P = scale(log(P/100)))


fviz_pca_biplot(PCAsoil,
             geom.ind = "point", # show points only (nbut not "text")
             pointshape = 19, pointsize = 1.5,
             col.ind = Sum.Muestreo_UR1$area, # color by groups
             palette = c("#00ba38","#f8766d"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Ambiente") #Grafica Figura 2

# Analisis de del impacto de CP1 sobre la tasa de colonización ############
# Creo que no necesitamos interpretar nada sobre el CP2

PCA1soila <- lm(PCAsoil1 ~ area, data = Sum.Muestreo_UR2)
anova(PCA1soila)
plot (PCA1soila)
summary (PCA1soila)

asin(sqrt(prop01))
col_pca1_ar_rural <- lm(prop01 ~ PCAsoil1 , data = Sum.Muestreo_UR2, , subset = area == "rural")
anova(col_pca1_ar_rural)
Anova(col_pca1_ar_rural)
s1<-summary(col_pca1_ar_rural)
fit1<-col_pca1_ar_rural

col_pca1_ar_urban <- lm(prop01 ~ PCAsoil1 , data = Sum.Muestreo_UR2, subset = area == "sidew")
anova(col_pca1_ar_urban)
Anova(col_pca1_ar_urban, type="2")
s2<-summary(col_pca1_ar_urban)
fit2<-col_pca1_ar_urban

df1 <- data.frame(x=1:3, y=1:3+rnorm(3))
df2 <- data.frame(x=1:3, y=1:3+rnorm(3))
fit1 <- lm(y~x, df1)
s1 <- summary(fit1)$coefficients
fit2 <- lm(y~x, df2)
s2 <- summary(fit2)$coefficients
db <- (s2[2,1]-s1[2,1])
sd <- sqrt(s2[2,2]^2+s1[2,2]^2)
df <- (fit1$df.residual+fit2$df.residual)
td <- db/sd
2*pt(-abs(td), df)








# Grafico de anova regresion
ggplot(Sum.Muestreo_UR2, aes(x=PCAsoil1, y=prop01*100, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+ stat_smooth(method="lm")+
  labs(title="Colonizacion x Soil PCA", 
       x="Soil PCA1", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))
  # guides(color=F, fill=F)

# urban
Sum.Muestreo_UR2%>%filter(area =="sidew")%>%
ggplot(., aes(x=PCAsoil1, y=prop01*100))+
  geom_point(pch=20, size=8, color = "#f8766d")+
  theme_classic()+ stat_smooth(method="lm", color = "black", size = 0.8)+
  xlab("Soil PC1") + ylab("AMF-colonization rate") +
  scale_x_continuous(breaks = seq(0, 4, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 80, by = 10))+
  theme(text=element_text(size=16,  family="Arial"))+
  theme(axis.text.x = element_text(color="black",  # face = "bold
                                   size=14, angle=0),
        axis.text.y = element_text(color="black", 
                                   size=14, angle=0),
        axis.line = element_line(colour = "black", 
                                        size = 1, linetype = "solid"))

#~/Downloads/Tesis_final


# Rural
Sum.Muestreo_UR2%>%filter(area =="rural")%>%
  ggplot(., aes(x = PCAsoil1, y=prop01*100))+
  geom_point(pch=20, size=8, color = "#00ba38")+
  theme_classic()+ stat_smooth(method="lm", color = "black", size = 0.8)+
  xlab("Soil PC1") + ylab("AMF-colonization rate") + 
  scale_x_continuous(breaks = seq(-4, 0, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 80, by = 10))+
  theme(text=element_text(size=16,  family="Arial"))+
  theme(axis.text.x = element_text(color="black",  # face = "bold
                                   size=14, angle=0),
        axis.text.y = element_text(color="black", 
                                   size=14, angle=0),
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))




# Evaluar las regresiones entr PCA y las vaeriables que son partes del PCA
# sirve para entender su relacion e interprtear mejor las escalas del PCA
pca_p<-lm(PCAsoil1~t.P, Sum.Muestreo_UR2)
summary (pca_p)

#### subset urban and rural

# soilpca_u<-subset(Sum.Muestreo_UR1, area=="sidew")[,c(25:30)]
# soilpca_r<-subset(Sum.Muestreo_UR1, area=="rural")[,c(25:30)]
# 
# PCAsoil_u<-prcomp(soilpca_u, scale = TRUE)
# summary(PCAsoil_u)
# 
# fviz_pca_biplot(PCAsoil_u,
#                 geom.ind = "point", # show points only (nbut not "text")
#                 pointshape = 19, pointsize = 1.5,
#               ) #Grafica Figura 2
# 
# PCAsoil_r<-prcomp(soilpca_r, scale = TRUE)
# summary(PCAsoil_r)
# 
# fviz_pca_biplot(PCAsoil_r,
#                 geom.ind = "point", # show points only (nbut not "text")
#                 pointshape = 19, pointsize = 1.5,
# ) #Grafica Figura 2
# 
# 
# 
# PCAsoil1_u<-PCAsoil_u[["x"]][,1]
# PCAsoil2_u<-PCAsoil_u[["x"]][,2]
# PCAsoil1_r<-PCAsoil_r[["x"]][,1]
# PCAsoil2_r<-PCAsoil_r[["x"]][,2]
# 
# 
# c(PCAsoil1_r,PCAsoil1_u)
# c(PCAsoil2_r,PCAsoil2_u)
# 
# Sum.Muestreo_UR2<-cbind(Sum.Muestreo_UR2, 
#                         sub_PCAsoil1=c(PCAsoil1_r,PCAsoil1_u), 
#                         sub_PCAsoil2=c(PCAsoil2_r,PCAsoil2_u))
# Sum.Muestreo_UR2[,c(1,45,47,46,48)]
# 
# 
# 
# col_pca1_ar<-lm(asin(sqrt(prop01))~sub_PCAsoil1, data = Sum.Muestreo_UR2)
# col_pca1_ar<-lm(asin(sqrt(prop01))~sub_PCAsoil1*area, data = Sum.Muestreo_UR2)
# anova(col_pca1_ar)
# Anova(col_pca1_ar, type="3")
# summary(col_pca1_ar)
# 
# ggplot(Sum.Muestreo_UR2, aes(x=sub_PCAsoil1, y=prop01*100, color=area))+
#   geom_point(aes(fill=area),pch=21,size=2)+
#   theme_classic()+
#   geom_smooth(method="lm", se=F)+
#   labs(title="Colonizacion x Soil PCA", 
#        x="Soil PCA1", y="% Colonizacion",
#        color="Area", fill="Area")+
#   scale_color_manual(labels=c("Rural","Profundo"), 
#                      values=c("#00ba38","#f8766d"))+
#   scale_fill_manual(labels=c("Rural","Profundo"),
#                     values=c("#00ba38","#f8766d"))+
#   guides(color=F, fill=F)
# 
# 
# 
# ### Standardize soilPCA1
# 
# std_PCAsoil1_u<-scale(subset(Sum.Muestreo_UR2, area=="sidew")[,c(45)])
# std_PCAsoil1_r<-scale(subset(Sum.Muestreo_UR2, area=="rural")[,c(45)])
# std_PCAsoil1<-c(std_PCAsoil1_r,std_PCAsoil1_u)
# 
# Sum.Muestreo_UR2<-cbind(Sum.Muestreo_UR2, 
#                         std_PCAsoil1)
# 
# col_pca1_ar<-lm(asin(sqrt(prop01))~std_PCAsoil1*area, 
#                 data = Sum.Muestreo_UR2)
# anova(col_pca1_ar)
# 
# 
# ggplot(Sum.Muestreo_UR2, aes(x=std_PCAsoil1, y=prop01*100, color=area))+
#   geom_point(aes(fill=area),pch=21,size=2)+
#   theme_classic()+
#   geom_smooth(method="lm", se=F)+
#   labs(title="Colonizacion x Soil PCA", 
#        x="Soil PCA1", y="% Colonizacion",
#        color="Area", fill="Area")+
#   scale_color_manual(labels=c("Rural","Profundo"), 
#                      values=c("#00ba38","#f8766d"))+
#   scale_fill_manual(labels=c("Rural","Profundo"),
#                     values=c("#00ba38","#f8766d"))+
#   guides(color=F, fill=F)
# 
# 
# 
# corr_PCA<-lm(PCAsoil1~sub_PCAsoil1, data = Sum.Muestreo_UR2)
# summary(corr_PCA)
# 
# cor(Sum.Muestreo_UR2$PCAsoil1, Sum.Muestreo_UR2$sub_PCAsoil1)
# 
# 
# ggplot(Sum.Muestreo_UR2, aes(x=sub_PCAsoil1, y=PCAsoil1, color=area))+
#   geom_point(aes(fill=area),pch=21,size=2)+
#   theme_classic()+
#   geom_smooth(method="lm", se=F)+
#   labs(title="Correlation PCA ~ subset PCA", 
#        x="subset PCA", y="PCA",
#        color="Area", fill="Area")+
#   scale_color_manual(labels=c("Rural","Profundo"), 
#                      values=c("#00ba38","#f8766d"))+
#   scale_fill_manual(labels=c("Rural","Profundo"),
#                     values=c("#00ba38","#f8766d"))+
#   guides(color=F, fill=F)

###################################################
####                  NMDS                     ####
###################################################

#NO SE MUESTRAN EN RESULTADOS

head(soilpca)
mtx_soil<-as.matrix(soilpca)
set.seed(123)
s_nmds <- metaMDS(mtx_soil, distance = "bray")
s_nmds
View(s_nmds)
plot(s_nmds)


envfit(nmds, env, permutations = 999, na.rm = TRUE)


soil_data.scores <- as.data.frame(scores(s_nmds))
soil_data.scores$Site <- Sum.Muestreo_UR1$Site
soil_data.scores$area <- Sum.Muestreo_UR1$area

ggplot(soil_data.scores, aes(x = NMDS1, y = NMDS2, color=area)) + 
  geom_point(size = 4)+ 
  theme_classic() + 
  labs(x = "NMDS1", colour = "Ambiente", y = "NMDS2")  + 
  scale_color_manual(labels=c("rural" = "Rural", 
                              "sidew" = "Profundo"),
                     values=c("#00ba38","#f8766d"))

s_anosim<-anosim(mtx_soil, Sum.Muestreo_UR1$area, 
       distance = "bray", permutations = 9999)



#SOIL AND SPORE
mtx_soilsp<-as.matrix(soil_spore_m)
set.seed(123)
ssp_nmds <- metaMDS(mtx_soilsp, distance = "bray")
ssp_nmds
plot(ssp_nmds)

soilsp_data.scores <- as.data.frame(scores(ssp_nmds))
soilsp_data.scores$Site <- Sum.Muestreo_UR1$Site
soilsp_data.scores$area <- Sum.Muestreo_UR1$area

ggplot(soilsp_data.scores, aes(x = NMDS1, y = NMDS2, color=area)) + 
  geom_point(size = 4)+ 
  theme_classic() + 
  labs(x = "NMDS1", colour = "Ambient", y = "NMDS2")  + 
  scale_color_manual(labels=c("rural" = "Rural", 
                              "sidew" = "Profundo"),
                     values=c("#00ba38","#f8766d"))

sp_anosim<-anosim(mtx_soilsp, Sum.Muestreo_UR1$area, 
       distance = "bray", permutations = 9999)




#########################################################
#########################################################
###                                                   ###
###                   HMA INTERACION                  ###
###                                                   ###
#########################################################
#########################################################

#ESTOS ANALISIS SE MUESTRAN EN TABLA 1



################################
#######  SPORE DENSITY   #######
################################



#######  Spore density vs Area ####### 

# Linear regression
sporesa<- lm(log(dense)~area, Sum.Muestreo_UR1)
summary(sporesa)
anova(sporesa) #En Resutlados

hist(log(Sum.Muestreo_UR1$dense))
plot(sporesa, 2)# data is NOT distributed normally
ols_test_normality(sporesa)#TESTING RESIDUALS


# Generalized Linear Model
sporesa0<- glm(dense~1, Sum.Muestreo_UR1, family=quasipoisson)
sporesa<- glm(dense~area, Sum.Muestreo_UR1, family=quasipoisson)
anova(sporesa0, sporesa)
Anova(sporesa)
summary(sporesa)
ls_means(sporesa)


SP<-summary(lm(dense~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))



#Graphic
Sum.Muestreo_UR1 %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(.,
         aes(x=area, 
             y=dense, fill=area))+
  geom_boxplot() + theme_classic()+
  labs(x="Subambiente", y="Número de esporas / 100 g")+
  guides(fill=F)+
  #geom_text(data=data.frame(Sum.Muestreo_UR1 %>% 
   #                           group_by(area) %>% 
    #                          summarize(Max.spore=max(dense, na.rm = T))),
            #aes(x=area,y=5+Max.spore,
            #    label=c("a","a","b")),vjust=0)+
  scale_fill_manual(values=c("#10cedbff","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural", 
                            "sidew" = "Urbano"))








Sum.Muestreo_UR1 %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(dense, na.rm = T),
    sd=sd(dense, na.rm = T),
    se=se(dense, na.rm = T)
  ) -> mean_spore
(mean_sporea<-as.data.frame(mean_spore)) #En Resutlados




##################################
#######  SPORE DIVERSITY   #######
##################################

div_spAr <- lm(SHANNON~area, data=Sum.Muestreo_UR1)
summary(div_spAr)
anova(div_spAr) #En Resutlados

hist(Sum.Muestreo_UR1$SHANNON)
plot(div_spAr,2)# data is distributed normally
ols_test_normality(div_spAr)#TESTING RESIDUALS



SP<-summary(lm(SHANNON~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))



Sum.Muestreo_UR1 %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew", "open"))) %>%
  ggplot(., 
         aes(x=area, y=SHANNON, fill=area))+
  theme_classic()+ geom_boxplot()+
  labs(title="Spore Diversity per area", 
       x="Area", y="Shannon Index")+
  scale_fill_manual(values=c("#10cedbff","#f8766d","#619cff"))+
  scale_x_discrete(labels=c("rural" = "Rural", 
                            "sidew" = "Urbano"))+
  guides(fill=F)



Sum.Muestreo_UR1 %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(SHANNON, na.rm = T),
    sd=sd(SHANNON, na.rm = T),
    se=se(SHANNON, na.rm = T)
  ) -> mean_diversity
(mean_diversity<-as.data.frame(mean_diversity)) #En Resutlados



############################
####### COLONIZATION #######
############################


########  Colonization vs area ########
bin_prop <- lm(asin(sqrt(prop01))~area, data= Muestreo_UR)
summary(bin_prop)
anova(bin_prop) 

hist(asin(sqrt(Muestreo_UR$prop01)))
plot(bin_prop,2)#as expected, the distribution is NOT normal
ols_test_normality(bin_prop)#TESTING RESIDUALS

# Generalized regresion with quasibinomial distribution
bin_prop <- glm(prop01 ~ area, 
                data = Muestreo_UR, 
                family = quasibinomial, 
                weights = total_field)
summary(bin_prop)
Anova(bin_prop, test="LR")#results dont change
                          



SP<-summary(lm(prop01*100~area, Sum.Muestreo_UR1)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4]))))


#Using sites as random effects and binomial distribution
bin_prop0 <- lmer(asin(sqrt(prop01)) ~ (1|site), 
                data = Muestreo_UR)
bin_prop1 <- lmer(asin(sqrt(prop01)) ~ area + (1|site), 
                 data = Muestreo_UR) #En Resutlados

bin_prop1 <- lmer(prop01 ~ area + (1|site), 
                  data = Muestreo_UR)

anova(bin_prop0, bin_prop1)#results dont change, but P-value is larger
Anova(bin_prop1)  #ddf="Kenward-Roger"
                  #En Resutlados
ranova(bin_prop1) #significant
r.squaredGLMM(bin_prop1) #
summary(bin_prop1)
ls_means(bin_prop1)



#Extracting means, sd, & se
colon_bar<-data.frame(
  area=c("rural","sidew"),
  mean=c(mean(Muestreo_UR$prop_colonization[Muestreo_UR$area=="rural"], na.rm = T),
         mean(Muestreo_UR$prop_colonization[Muestreo_UR$area=="sidew"], na.rm = T)),
  sd=c(sd(Muestreo_UR$prop_colonization[Muestreo_UR$area=="rural"], na.rm = T),
       sd(Muestreo_UR$prop_colonization[Muestreo_UR$area=="sidew"], na.rm = T)),
  se=c(se(Muestreo_UR$prop_colonization[Muestreo_UR$area=="rural"], na.rm = T),
       se(Muestreo_UR$prop_colonization[Muestreo_UR$area=="sidew"], na.rm = T)))
colon_bar #En Resutlados


#Graphic
ggplot(Muestreo_UR, 
         aes(x=area, y=prop01*100, fill=area))+
  theme_classic()+ 
  geom_boxplot()+
  labs(x="Area", y="% Colonizacion")+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)
  geom_text(aes(x=area,y=2+mean+se,
                label=c("a","b")),vjust=0)


colon_bar %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>% 
  ggplot(., 
         aes(x=area, y=mean, fill=area))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=area, ymin=mean-se, ymax=mean+se, 
                    width=0.1))+
  labs(x="Area", y="% Colonizacion")+
  geom_text(aes(x=area,y=2+mean+se,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)



bin_prop1gg<-ggpredict(bin_prop1, terms="area")


ggplot(bin_prop1gg, aes(x, predicted, color=x))+
  geom_point(position = position_dodge(0))+
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high)
    #position = position_dodge(0.1)
    )+ theme_classic()+
  labs(x="Area", y="% Colonizacion")+
  geom_text(aes(x=x, y=conf.high+0.02,
                label=c("a","b")),
            color="black",vjust=0)+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  scale_color_manual(values=c("#00ba38","#f8766d"))+
  guides(color=F)


  
#########################################################
#########################################################
###                                                   ###
###                   COORDENADAS                     ###
###                                                   ###
#########################################################
#########################################################

## NO SE MUESTRAN EN RESULTADOS


###################
### COL & COORD ###
###################


site.dists <- dist(cbind(Sum.Muestreo_UR1$long, Sum.Muestreo_UR1$lat))
col.dists <- dist(Sum.Muestreo_UR1$prop01)
vol.dists <- dist(Sum.Muestreo_UR1$volumen_ml)
mantel.rtest(site.dists,col.dists, nrepet = 9999) #NOT Significant
mantel.rtest(site.dists,vol.dists, nrepet = 9999) #NOT Significant
  

#####################
### SPORE & COORD ###
#####################


dense.dists <- dist(Sum.Muestreo_UR1$dense)
diverse.dists <- dist(Sum.Muestreo_UR1$SHANNON)
mantel.rtest(site.dists,dense.dists,nrepet = 9999)#NOT Significant
mantel.rtest(site.dists,diverse.dists,nrepet = 9999)#NOT Significant


###################
##  SOIL & COORD ##
###################


N.dists <- dist(Sum.Muestreo_UR1$N)
K.dists <- dist(Sum.Muestreo_UR1$K)
P.dists <- dist(Sum.Muestreo_UR1$P)
pH.dists <- dist(Sum.Muestreo_UR1$pH)
mantel.rtest(site.dists, N.dists, nrepet = 9999)#Significant
mantel.rtest(site.dists, P.dists, nrepet = 9999)#NOT Significant
mantel.rtest(site.dists, K.dists, nrepet = 9999)#Significant
mantel.rtest(site.dists, pH.dists, nrepet = 9999)#NOT Significant



##### COORDENADAS

# NITROGEN

f.N1 <- formula(N~area)
N.gls <- nlme::gls(f.N1, data=Sum.Muestreo_UR1,
                   na.action = na.omit)

N.gls.A<-gls(f.N1, correlation = corSpher(form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1)
N.gls.C<-gls(f.N1, correlation = corRatio(form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1)
N.gls.D<-gls(f.N1, correlation = corGaus (form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1)
N.gls.E<-gls(f.N1, correlation = corExp  (form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1)

AIC(N.gls, N.gls.A, N.gls.C, N.gls.D, N.gls.E)
anova(N.gls, N.gls.C)

summary(N.gls)
summary(N.gls.C) #Still Significant



#further revisions
U_subset<-subset(Sum.Muestreo_UR1, area=="sidew")#subset urban sites
R_subset<-subset(Sum.Muestreo_UR1, area=="rural")#subset rural sites
U_site.dists <- dist(cbind(U_subset$long, U_subset$lat))
R_site.dists <- dist(cbind(R_subset$long, U_subset$lat))
U_N.dists <- dist(U_subset$N)
R_N.dists <- dist(R_subset$N)
mantel.rtest(U_site.dists, U_N.dists, nrepet = 9999)#  NOT significant (P=0.2)
mantel.rtest(R_site.dists, R_N.dists, nrepet = 9999)#  NOT significant  (P=0.8)
# Urban points are closer to each other
# The general analyses confuses the effect of Urban Area with distances?



# POTASSIUM

f.K1 <- formula(K~area)
K.gls <- nlme::gls(f.K1, data=Sum.Muestreo_UR1,
                        na.action = na.omit)

K.gls.A<-gls(f.K1, correlation = corSpher(form=~long+lat, nugget=T),
              data = Sum.Muestreo_UR1)
K.gls.B<-gls(f.K1, correlation = corLin  (form=~long+lat, nugget=T),
              data = Sum.Muestreo_UR1)
K.gls.C<-gls(f.K1, correlation = corRatio(form=~long+lat, nugget=T),
              data = Sum.Muestreo_UR1)
K.gls.D<-gls(f.K1, correlation = corGaus (form=~long+lat, nugget=T),
              data = Sum.Muestreo_UR1)
K.gls.E<-gls(f.K1, correlation = corExp  (form=~long+lat, nugget=T),
              data = Sum.Muestreo_UR1)

AIC(K.gls, K.gls.A, K.gls.B, 
    K.gls.C, K.gls.D, K.gls.E)
anova(K.gls, K.gls.E)

summary(K.gls)
summary(K.gls.E) #NOT longer Significant
anova(K.gls.E)


#further revisions
U_K.dists <- dist(U_subset$K)
R_K.dists <- dist(R_subset$K)
mantel.rtest(U_site.dists, U_K.dists, nrepet = 9999)#  NOT significant (P=0.8)
mantel.rtest(R_site.dists, R_K.dists, nrepet = 9999)#  NOT significant  (P=0.8)
# Are these distance analyses removing the URBAN effect?
# Results from these analyses might confuse interpretation.     

   
##### NEW

# ALL

f.all1 <- formula(prop01~N+P+K+pH) 
all.gls <- nlme::gls(f.all1, data=Sum.Muestreo_UR1,
                   na.action = na.omit) 

all.gls.A<-gls(f.all1, correlation = corSpher(form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1) 
all.gls.C<-gls(f.all1, correlation = corRatio(form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1) 
all.gls.D<-gls(f.all1, correlation = corGaus (form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1) 
all.gls.E<-gls(f.all1, correlation = corExp  (form=~long+lat, nugget=T),
             data = Sum.Muestreo_UR1) 

AIC(all.gls, all.gls.A, all.gls.C, all.gls.D, all.gls.E)
anova(all.gls, all.gls.A) #Not Significant

summary(all.gls)
summary(all.gls.C) 
anova(all.gls)
anova(all.gls.A)


##########################################
# Multiple interaction
##########################################


all_soil01<-glm(prop01~N*P*K*pH, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
summary(all_soil01)
anova(all_soil01)
Anova(all_soil01, type="2")


all_soil01ar<-glm(prop01~N*P*K*pH*area, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
summary(all_soil01ar)
anova(all_soil01ar)
Anova(all_soil01ar, type="2")


all_soil01<-glm(prop01~N*P, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
anova(all_soil01)
Anova(all_soil01, type="3")

all_soil01<-glm(prop01~P*K, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
anova(all_soil01)
Anova(all_soil01, type="3")


all_soil01<-glm(prop01~P*K, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
anova(all_soil01)
Anova(all_soil01, type="3")


#########################################################
#########################################################
###                                                   ###
###              PLANT CHARACTERISTICS                ###
###                                                   ###
#########################################################
#########################################################




############################
####### ROOT VOLUME  #######
############################

#Linear Model
vola <- lm(log(volumen_ml) ~ area, data=Muestreo_UR)
summary(vola)
anova(vola)

hist(log(Muestreo_UR$volumen_ml))
plot(vola, 2)# the distribution is NOT normal, log --> doesn't change
ols_test_normality(vola)#TESTING RESIDUALS


vola1 <- lmer(log(volumen_ml) ~ (1|site), data=Muestreo_UR)
vola2 <- lmer(log(volumen_ml) ~ area + (1|site), data=Muestreo_UR)
anova(vola1, vola2)
anova(vola2)
Anova(vola2) #En Resutlados tesis
ranova(vola2)


vola <- lmer(volumen_ml ~ area + (1|site), data=Muestreo_UR)
ls_means(vola)



ggpredict(vola, "area") %>% 
  plot()



#Gr?fico 
Muestreo_UR %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(., aes(x=area, y=r_wet, fill=area))+
  geom_boxplot() + theme_classic()+
  labs(x="Area", y="Peso de Raíz (g)")+
  guides(fill=F)+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.vol=max(r_wet, na.rm = T))),
            aes(x=area,y=0.5+Max.vol,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)




Muestreo_UR %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(volumen_ml, na.rm = T),
    sd=sd(volumen_ml, na.rm = T),
    se=se(volumen_ml, na.rm = T)
  ) -> mean_vol

(mean_vol<-as.data.frame(mean_vol))



#############################
#######  ROOT LENGTH  #######
#############################  

#Linear Model
lra <- lm(log(long_root) ~ area, data=Muestreo_UR)
summary(lra)
anova(lra)

hist(log(Muestreo_UR$long_root))
plot(lra, 2)# the distribution is NOT normal, lod --> does
ols_test_normality(lra)#TESTING RESIDUALS


# Linear Mixed Models
lra1 <- lmer(log(long_root) ~ (1|site), data=Muestreo_UR)
lra2 <- lmer(log(long_root) ~ area + (1|site), data=Muestreo_UR) #En Resutlados
anova(lra1, lra2)
Anova(lra2) #En Resutlados
ranova(lra2)


lra <- lmer(long_root ~ area + (1|site), data=Muestreo_UR)
ls_means(lra)
summary(lra)


#Gr?fico 
ggplot(Muestreo_UR, aes(x=area, y=long_root, fill=area))+
  geom_boxplot()+theme_classic()+
  labs(title="Root Length by Area", 
       x="Area", y="Length (cm)")+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.vol=max(long_root, na.rm = T))),
            aes(x=area,y=2+Max.vol,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)



Muestreo_UR %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(long_root, na.rm = T),
    sd=sd(long_root, na.rm = T),
    se=se(long_root, na.rm = T)
  ) -> mean_rootl

(mean_rootl<-as.data.frame(mean_rootl)) #En Resutlados






#################################
#######   AEREAL LENGTH   #######
#################################

#Linear Model
laa <- lm(log(long_plant) ~ area, data=Muestreo_UR)
summary(laa)
anova(laa)

hist(log(Muestreo_UR$long_plant))
plot(laa, 2)# the distribution is NOT normal, log --> doesn't change
ols_test_normality(laa)#TESTING RESIDUALS


# Generalized linear model
laa <- glm(long_plant~area, Muestreo_UR, family=Gamma)
summary(laa)
Anova(laa)


# Linear Mixed Models
laa1 <- lmer(log(long_plant) ~ (1|site), data=Muestreo_UR)
laa2 <- lmer(log(long_plant) ~ area + (1|site), data=Muestreo_UR) #En Resutlados
anova(laa1, laa2)
Anova(laa2) #En Resutlados
ranova(laa2)


laa <- lmer(long_plant ~ area + (1|site), data=Muestreo_UR)
ls_means(lra)


#Gr?fico 
ggplot(Muestreo_UR, aes(x=area, y=long_plant, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs(title="Top Length by Area", 
       x="Area", y="Length (cm)")+
 
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)



Muestreo_UR %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(long_plant, na.rm = T),
    sd=sd(long_plant, na.rm = T),
    se=se(long_plant, na.rm = T)
  ) -> mean_plantl

(mean_plantl<-as.data.frame(mean_plantl)) 




##############################
#######  TOTAL WEIGHT  #######
##############################

#Linear Model
wt2<-c(Muestreo_UR$w_above+Muestreo_UR$r_wet)[-c(1:42)]
wt1<-Muestreo_UR$w_total[1:42]
na1<-rep(NA,42)
w_total1<-c(wt1,wt2)
w_total2<-c(na1,wt2)
Muestreo_UR<-cbind(Muestreo_UR,w_total1)
Muestreo_UR<-cbind(Muestreo_UR,w_total2)

wt1s<-Muestreo_UR %>% 
  group_by(site) %>% 
  summarize(w_total1 = mean(w_total1,na.rm=T))

Sum.Muestreo_UR1<-left_join(Sum.Muestreo_UR1, wt1s)


waa <- lm(log(w_total2) ~ area, data=Muestreo_UR)
summary(waa)
anova(waa)

hist(log(Muestreo_UR$w_total2))
plot(waa, 2)# the distribution is NOT normal --> log does
ols_test_normality(waa)#TESTING RESIDUALS


# Regresi?n lineal 
waa<- glm(w_total1~area, Muestreo_UR, family=Gamma)
summary(waa)
Anova(waa, test="LR")


# Linear Mixed Models
waa1 <- lmer(log(w_total2) ~ (1|site), data=Muestreo_UR)
waa2 <- lmer(log(w_total2) ~ area + (1|site), data=Muestreo_UR) #En Resutlados
anova(waa1, waa2)
anova(waa2)
summary(waa2)
Anova(waa2)#En Resutlados
ranova(waa2)


waa <- lmer(w_total1 ~ area + (1|site), data=Muestreo_UR)
ls_means(waa)
waa <- lmer(w_total2 ~ area + (1|site), data=Muestreo_UR)
ls_means(waa)


# Generalized linear mixed model
waa1 <- glmer(w_total1 ~ (1|site), data=Muestreo_UR, family=Gamma)
waa2 <- glmer(w_total1 ~ area + (1|site), data=Muestreo_UR, family=Gamma)
anova(waa1, waa2)
Anova(waa2)
ranova(waa2)





SP<-summary(lm(w_total~area, Muestreo_UR)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4])))) #En Resutlados


#Gr?fico 
Muestreo_UR %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(., aes(x=area, y=w_total, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs( x="Area", y="Peso (g)")+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.weight=max(w_total, na.rm = T))),
            aes(x=area,y=2+Max.weight,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)




Muestreo_UR %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(w_total, na.rm = T),
    sd=sd(w_total, na.rm = T),
    se=se(w_total, na.rm = T)
  ) -> mean_plantw

(mean_plantw<-as.data.frame(mean_plantw)) 



################################
####### Num. de BRANCHES #######
################################

#Linear Model
nbrancha <- lm(log(no_branch) ~ area, data=Muestreo_UR)
summary(nbrancha)
anova(nbrancha)

hist(log(Muestreo_UR$no_branch)) #transformation makes them worse
plot(nbrancha, 2)# the distribution is NOT normal
ols_test_normality(nbrancha)#TESTING RESIDUALS


# Regresi?n lineal 
nbrancha<- glm(no_branch~area, Muestreo_UR, family=poisson)
summary(nbrancha)#no overdisperssion, nothing changes
Anova(nbrancha, test="LR")


#Mixed model
nbrancha1 <- lmer(no_branch ~ (1|site), data=Muestreo_UR)
nbrancha2 <- lmer(log(no_branch) ~ area + (1|site), data=Muestreo_UR) #En Resutlados
anova(nbrancha1, nbrancha2) #no difference
Anova(nbrancha2) #En Resutlados
ranova(nbrancha2) #random effect not important




SP<-summary(lm(no_branch~area, Muestreo_UR)) 
cbind(area=c("rural", "urban"),
      as.data.frame(cbind(mean=c(SP$coefficients[1],
                                 SP$coefficients[1]+SP$coefficients[2]),
                          se=c(SP$coefficients[3],SP$coefficients[4])))) #En Resutlados


#Gr?fico 
ggplot(Muestreo_UR, aes(x=area, y=no_branch, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs(#title="Number of Branches by Area", 
       x="Area", y="Number of Branches")+

  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)



Muestreo_UR %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(no_branch, na.rm = T),
    sd=sd(no_branch, na.rm = T),
    se=se(no_branch, na.rm = T)
  ) -> mean_nobranch

(mean_nobranch<-as.data.frame(mean_nobranch))


############################
####### Num. de ROOTS ######
############################

#Linear Model
nroota <- lm(log(no_roots) ~ area, data=Muestreo_UR)
summary(nroota)
anova(nroota)


hist(log(Muestreo_UR$no_roots))
plot(nroota, 2)# the distribution is NOT normal, log --> does
ols_test_normality(nroota)#TESTING RESIDUALS


# Generalized linear model
nroota<- glm(no_roots~area, Muestreo_UR, family=poisson)
summary(nroota)
Anova(nroota, test="LR")


# Linear Mixed Models
nroota1 <- lmer(log(no_roots) ~ (1|site), data=Muestreo_UR)
nroota2 <- lmer(log(no_roots) ~ area + (1|site), data=Muestreo_UR) #En Resutlados
anova(nroota1, nroota2)
Anova(nroota2) #En Resutlados
ranova(nroota2)#random effect is important


nroota <- lmer(no_roots ~ area + (1|site), data=Muestreo_UR)
ls_means(nroota)
summary(nroota)


# Generalized linear mixed model
nroota1 <- glmer(no_roots ~ (1|site), data=Muestreo_UR, family=poisson)
nroota2 <- glmer(no_roots ~ area + (1|site), data=Muestreo_UR, family=poisson)
anova(nroota1, nroota2)
Anova(nroota2)
ranova(nroota2)

#Gr?fico
ggplot(Muestreo_UR, aes(x=area, y=no_roots, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs(#title="Number of Roots by Area", 
       x="Area", y="Number of Roots")+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.root=max(no_roots, na.rm = T))),
            aes(x=area,y=5+Max.root,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)


Muestreo_UR %>% 
  group_by(area) %>% 
  summarise(
    mean=mean(no_roots, na.rm = T),
    sd=sd(no_roots, na.rm = T),
    se=se(no_roots, na.rm = T)
  ) -> mean_noroots

(mean_noroots<-as.data.frame(mean_noroots)) #En Resutlados




####################################
#######  SPECIFIC LEAF AREA  #######
####################################
## NO ESTA EN RESULTADOS

#Linear Model
SLAa <- lm(SLFA ~ prop01, data = Muestreo_UR, subset = area == "sidew")
summary(SLAa)
anova(SLAa)

hist(Muestreo_UR$w_s_area)
plot(SLAa, 2)# the distribution is NOT normal
ols_test_normality(SLAa)#TESTING RESIDUALS

# Linear Mixed Models
SLAa1 <- lmer(SLFA ~ (1|site), data = Muestreo_UR)
SLAa2 <- lmer(F ~ area + (1|site), data = Muestreo_UR)
anova(SLAa1, SLAa2)
Anova(SLAa2)
ranova(SLAa2)


SLAa <- lmer(w_s_area*1000 ~ area + (1|site), data = Muestreo_UR)
ls_means(SLAa)
summary(SLAa)

 
#Gr?fico 
Muestreo_UR %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(., aes(x=area, y=w_total, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs( x="Area", y="Peso (g)")+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.weight=max(w_total, na.rm = T))),
            aes(x=area,y=2+Max.weight,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)


########################################################################
#######################                   ##############################
#######################      BIOMASS      ##############################
#######################                   ##############################
########################################################################

####################################
#######     BIOMASS ROOT     #######
####################################

#Linear Model
mass_roota <- lm(log(mass_r) ~ area, data=Muestreo_UR)
summary(mass_roota)
anova(mass_roota) #significant

hist(log(Muestreo_UR$mass_r))
plot(mass_roota, 2)# the distribution is NOT normal, use log
ols_test_normality(mass_roota)#TESTING RESIDUALS, log is normal

# Linear Mixed Models
mass_roota1 <- lmer(log(mass_r) ~ (1|site), data=Muestreo_UR)
mass_roota2 <- lmer(log(mass_r) ~ area + (1|site), data=Muestreo_UR)
anova(mass_roota1, mass_roota2)
Anova(mass_roota2) #significant
ranova(mass_roota2) #no point


#Gr?fico 
Muestreo_UR %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(., aes(x=area, y=mass_r, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs( x="Environment", y="Root mass (g)")+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.mass_r=max(mass_r, na.rm = T))),
            aes(x=area,y=2+Max.mass_r,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Urban"))+
  guides(fill=F)


####################################
#######   BIOMASS SHOOT     #######
####################################

#Linear Model
mass_shoota <- lm(log(mass_a) ~ area, data=Muestreo_UR)
summary(mass_shoota)
anova(mass_shoota) #significant

hist(log(Muestreo_UR$mass_a))
plot(mass_shoota, 2)# the distribution is NOT normal, use log
ols_test_normality(mass_shoota)#TESTING RESIDUALS, log is normal

# Linear Mixed Models
mass_shoota1 <- lmer(log(mass_a) ~ (1|site), data=Muestreo_UR)
mass_shoota2 <- lmer(log(mass_a) ~ area + (1|site), data=Muestreo_UR)
anova(mass_shoota1, mass_shoota2)
Anova(mass_shoota2) #not significant
ranova(mass_shoota2) #use fix model !!


#Gr?fico 
Muestreo_UR %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(., aes(x=area, y=mass_a, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs( x="Environment", y="Shoot mass (g)")+
  #geom_text(data=data.frame(Muestreo_UR %>% 
   #                           group_by(area) %>% 
    #                          summarize(Max.mass_a=max(mass_a, na.rm = T))),
     #       aes(x=area,y=2+Max.mass_a,
      #          label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Urban"))+
  guides(fill=F)


####################################
#######    BIOMASS TOTAL     #######
####################################

#Linear Model
mass_totala <- lm(log(mass_t) ~ area, data=Muestreo_UR)
summary(mass_totala)
anova(mass_totala) #significant

hist(log(Muestreo_UR$mass_t))
plot(mass_totala, 2)# the distribution is NOT normal, use log
ols_test_normality(mass_totala)#TESTING RESIDUALS, log is normal

# Linear Mixed Models
mass_totala1 <- lmer(log(mass_t) ~ (1|site), data=Muestreo_UR)
mass_totala2 <- lmer(log(mass_t) ~ area + (1|site), data=Muestreo_UR)
anova(mass_totala1, mass_totala2)
Anova(mass_totala2) #significant
ranova(mass_totala2) #no point


#Gr?fico 
Muestreo_UR %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(., aes(x=area, y=mass_t, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs( x="Environment", y="Totalbiomass (g)")+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.mass_t=max(mass_t, na.rm = T))),
            aes(x=area,y=2+Max.mass_t,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Urban"))+
  guides(fill=F)



####################################
#######     BIOMASS SLFA     #######
####################################

#Linear Model
SLAa <- lm(SLFA ~ area, data=Muestreo_UR)
summary(SLAa)
anova(SLAa) #significant

hist(Muestreo_UR$SLFA)
plot(SLAa, 2)# the distribution is normal, YES
ols_test_normality(SLAa)#TESTING RESIDUALS

# Linear Mixed Models
SLAa1 <- lmer(SLFA ~ (1|site), data=Muestreo_UR)
SLAa2 <- lmer(SLFA ~ area + (1|site), data=Muestreo_UR)
anova(SLAa1, SLAa2)
Anova(SLAa2)  # significant
ranova(SLAa2) # use fix model !!


#Gr?fico 
Muestreo_UR %>% 
  dplyr::mutate(area = factor(area, 
                              levels = c("rural", "sidew"))) %>%
  ggplot(., aes(x=area, y=SLFA, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs( x="Area", y="Specific Leaf Area (mg/unit squared)")+
  geom_text(data=data.frame(Muestreo_UR %>% 
                              group_by(area) %>% 
                              summarize(Max.SLFA=max(SLFA, na.rm = T))),
            aes(x=area,y=2+Max.SLFA,
                label=c("a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Profundo"))+
  guides(fill=F)




###ROOT BIOMASS vs COL

biomass_col <- lm(log(mass_r) ~ prop01, data=Muestreo_UR)
plot(biomass_col, 2)# the distribution is normal, YES
ols_test_normality(biomass_col)#TESTING RESIDUALS
summary(biomass_col)
anova(biomass_col)


ggplot(Muestreo_UR, aes(x=prop01, y=mass_r))+
  geom_point(pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonizacion x Root biomass", 
       x="Colonizacion rate", y="Root biomass")

#by area
biomass_col <- lm(log(mass_r) ~ prop01*area, data=Muestreo_UR)
plot(biomass_col, 2)# the distribution is normal, YES
ols_test_normality(biomass_col)#TESTING RESIDUALS
summary(biomass_col)
Anova(biomass_col, type=3)

ggplot(Muestreo_UR, aes(x=prop01, y=mass_r, color=area))+
  geom_point(pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonizacion x Root biomass", 
       x="Colonizacion rate", y="Root biomass",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)




################# PCA PLANT CHARACTERISTICS         ####
## ESTA EN RESULTADOS


plant_root<-Muestreo_UR[,c(26,3,4,  12,14,17)] #by plant
sum.plant_root<-Sum.Muestreo_UR1[,c(1,2,  12,14,17)] #by site


planta_01_m<-Muestreo_UR[-c(1:42),c(26, 3, 12,14:15,17,19:21,27)] #w/ colnizacion by plant
#row.names(planta_01_m)<-Muestreo_UR$site


## WITH COLONIZATION
#no uasdo en resultados

planta_01_m %>% 
  group_by(site) %>% 
  summarise(long_root = mean(long_root, na.rm=T),
            no_roots = mean(no_roots , na.rm=T),
            no_branch = mean(no_branch , na.rm=T),
            volumen_ml = mean(volumen_ml , na.rm=T),
            w_above = mean(w_above , na.rm=T),
            r_wet = mean(r_wet , na.rm=T),
            w_s_area = mean(w_s_area, na.rm=T),
            prop01 = mean(prop01, na.rm=T)
  ) -> plantpca_sum

plantpca_sum<-as.data.frame(plantpca_sum)
plantpca_sum


planta_01_m[is.na(planta_01_m$prop01),]


planta_01_m[112,3]<-plantpca_sum[2,2]
planta_01_m[121,3]<-plantpca_sum[3,2]
planta_01_m[112,4]<-plantpca_sum[2,3]
planta_01_m[121,4]<-plantpca_sum[3,3]
planta_01_m[112,5]<-plantpca_sum[2,4]
planta_01_m[114,5]<-plantpca_sum[2,4]
planta_01_m[118,6]<-plantpca_sum[3,5]
planta_01_m[123,6]<-plantpca_sum[4,5]
planta_01_m[112,7]<-plantpca_sum[2,6]
planta_01_m[113,7]<-plantpca_sum[2,6]
planta_01_m[114,7]<-plantpca_sum[2,6]
planta_01_m[112,9]<-plantpca_sum[2,8]
planta_01_m[114,9]<-plantpca_sum[2,8]
planta_01_m[56,10]<-plantpca_sum[22,9]

planta_01_m

plantpca<-planta_01_m[,-10]

rootpca<-planta_01_m[,-c(5,7:10)]
root01pca<-planta_01_m[,-c(5,7:9)]



PCAplant<-prcomp(plantpca[,-c(1,2)], scale = TRUE)
summary(PCAplant)

PCAplant1<-PCAplant[["x"]][,1]


get_pca_var(PCAplant)
get_pca_ind(PCAplant)

NAv<-rep(NA,42)
NAv[43:203]<-as.vector(PCAplant1)

Muestreo_UR1<-cbind(Muestreo_UR, NAv)
names(Muestreo_UR1)[36]<-"PCAplant1"



fviz_pca_biplot(PCAplant,
             geom.ind = "point", # show points only (nbut not "text")
             pointshape = 19, pointsize = 1.5,
             col.ind = plantpca$area, # color by groups
             palette = c("#00ba38","#f8766d"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)





###### ROOT
#Usado en resutlados

plant_root<-na.omit(plant_root)


PCA_root<-prcomp(plant_root[,-c(1:3)], scale = TRUE)
summary(PCA_root)

PCA_sumroot<-prcomp(sum.plant_root[,-c(1:2)], scale = TRUE)
summary(PCA_sumroot)


par(mfrow=c(1,1))

fviz_pca_biplot(PCA_root,
                geom.ind = "point", # show points only (nbut not "text")
                pointshape = 19, pointsize = 1.5,
                col.ind = plant_root$area, # color by groups
                palette = c("#00ba38","#f8766d"),
                addEllipses = TRUE, # Concentration ellipses
                legend.title = "Ambient",
                title = NULL
)+theme(legend.position = "none") -> pcaplot #EN Resultados, figura 2


fviz_pca_biplot(PCA_sumroot,
                geom.ind = "point", # show points only (nbut not "text")
                pointshape = 19, pointsize = 1.5,
                col.ind = sum.plant_root$area, # color by groups
                palette = c("#00ba38","#f8766d"),
                addEllipses = TRUE, # Concentration ellipses
                legend.title = "Ambient",
                title = NULL
)+theme(legend.position = "none")-> pcaplotsum


pcaplotsum

grid.arrange(pcaplot,pcaplotsum,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))

#PCA1

PCA1_root<-PCA_root[["x"]][,1]
PCA1_roots<-PCA_sumroot[["x"]][,1]

plant_root1<-cbind(plant_root,PCA1_root)
sum.plant_root<-cbind(sum.plant_root,PCA1_roots)

lm_PCAroot1<-lm(PCA1_root~area,data=plant_root1)
lm_PCAroot1<-lm(PCA1_roots~area,data=sum.plant_root)
summary(lm_PCAroot1)
anova(lm_PCAroot1)


ggplot(plant_root1, aes(x=area, y=PCA1_root, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs(#title="Number of Roots by Area", 
    x="Area", y="Tamaño de Raíz (PCA1)")+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Urbano"))+
  guides(fill=F)

ggplot(sum.plant_root, aes(x=area, y=PCA1_roots, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs(#title="Number of Roots by Area", 
    x="Area", y="Tamaño de Raíz (PCA1)")+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Urbano"))+
  guides(fill=F)



#PCA2

PCA2_root<-PCA_root[["x"]][,2]
PCA2_roots<-PCA_sumroot[["x"]][,2]

plant_root1<-cbind(plant_root1,PCA2_root)
sum.plant_root<-cbind(sum.plant_root,PCA2_roots)

lm_PCAroot2<-lm(PCA2_root~area,data=plant_root1)
lm_PCAroot2<-lm(PCA2_roots~area,data=sum.plant_root)
summary(lm_PCAroot2)
anova(lm_PCAroot2)


ggplot(plant_root1, aes(x=area, y=PCA2_root, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs(#title="Number of Roots by Area", 
    x="Area", y="Forma de Raíz (PCA2)")+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Urbano"))+
  guides(fill=F)


ggplot(sum.plant_root, aes(x=area, y=PCA2_roots, fill=area))+
  geom_boxplot()+ theme_classic()+
  labs(#title="Number of Roots by Area", 
    x="Area", y="Forma de Raíz (PCA2)")+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("rural" = "Rural",
                            "sidew" = "Urbano"))+
  guides(fill=F)



######################
### PLANTS VS AREA ### 
######################
#PCA, no se uso

#Linear Model
PCAplantAr<-lm(PCAplant1~area, 
                data=Muestreo_UR1)

shapiro.test(Muestreo_UR1$PCAplant1)
plot(PCAplantAr, 2)#data is NOT normally distributed
ols_test_normality(PCAplantAr)#TESTING RESIDUALS
summary(PCAplantAr)

#to work with Gamma
Muestreo_UR1$PCAplant1_10<-Muestreo_UR1$PCAplant1+10

#Generalized Linear Model
PCAplantAr<-glm(PCAplant1_10~area, 
                 data=Muestreo_UR1,
                 family = Gamma)
Anova(PCAplantAr)

#Generalized Linear Mixed Model
PCAplantAr1<-glmer(PCAplant1_10~ (1|site), 
                data=Muestreo_UR1,
                family = Gamma)
PCAplantAr2<-glmer(PCAplant1_10~ area + (1|site), 
                  data=Muestreo_UR1,
                  family = Gamma)
anova(PCAplantAr1,PCAplantAr2)
Anova(PCAplantAr2)


ggplot(Muestreo_UR1, aes(x=area, y=PCAplant1, fill=area))+
  geom_boxplot()+
  theme_classic()+
  labs(title="Soil PCA by area", 
       x="Area", y="Plant PCA Dim1",
       color="Area", fill="Area")+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(fill=F)



#####################
### PLANTS VS COL ### 
#####################
#PCA, no se uso

#Linear Model
colPCAplant<-lm(PCAplant1~prop01, 
                data=Muestreo_UR1)


plot(colPCAplant, 2)#data is NOT normally distributed
ols_test_normality(colPCAplant)#TESTING RESIDUALS
summary(colPCAplant)

#Generalized Linear Model
colPCAplant<-glm(PCAplant1~prop01, 
                data=Muestreo_UR1,
                family = Gamma)
summary(colPCAplant)


#Generalized Linear Mixed Model
colPCAplant1<-lmer(PCAplant1~ (1|site), 
                   data=Muestreo_UR1,
                   subset = !is.na(prop01))
colPCAplant2<-lmer(PCAplant1~ prop01 + (1|site), 
                   data=Muestreo_UR1,
                   subset = !is.na(prop01))
anova(colPCAplant1,colPCAplant2)
Anova(colPCAplant2)


ggplot(Muestreo_UR1, aes(x=prop01*100, y=PCAplant1))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonizacion x Plant PCA", 
       x="Plant PCA Dim1", y="% Colonizacion")



#BY AREA


#Linear Model
colPCAplantAr<-lm(PCAplant1~prop01*area, 
                   data=Muestreo_UR1)
plot(colPCAplantAr, 2)#data is NOT normally distributed
ols_test_normality(colPCAplantAr)#TESTING RESIDUALS

summary(colPCAplantAr)
anova(colPCAplantAr)

#Generalized Linear Model
colPCAplantAr<-glm(PCAplant1~prop01*area, 
                   data=Muestreo_UR1,
                   family = Gamma)
summary(colPCAplantAr)
Anova(colPCAplantAr, type="3", test="LR")


#Generalized Linear Mixed Model
colPCAplantAr0<-lmer(PCAplant1~ 1+ (1|site), 
                     data=Muestreo_UR1,
                     subset = !is.na(prop01))
colPCAplantAr1<-lmer(PCAplant1~ prop01 +(1|site), 
                    data=Muestreo_UR1,
                    subset = !is.na(prop01))
colPCAplantAr2<-lmer(PCAplant1~ prop01 * area + (1|site), 
                    data=Muestreo_UR1,
                    subset = !is.na(prop01))
anova(colPCAplantAr0,colPCAplantAr1)
anova(colPCAplantAr1,colPCAplantAr2)

Anova(colPCAplantAr2, type="3")

ggplot(Muestreo_UR1, aes(x=prop01*100, y=PCAplant1, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonizacion x Soil PCA", 
       x="Plant PCA Dim1", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)





##########################################
# COLONIZATION VS ROOT PCA
##########################################
#PCA, no se uso

Sum.Muestreo_UR1_pca<-left_join(Sum.Muestreo_UR1 ,sum.plant_root, by="site")


# PCA1 - SIZE
col_root<-lm(PCA1_roots~prop01, data=Sum.Muestreo_UR1)
plot(col_root,2)
summary(col_root)

col_roota<-lm(PCA1_roots~prop01*area, data=Sum.Muestreo_UR1)
Anova(col_roota, type="2", test="F")


ggplot(Sum.Muestreo_UR1, aes(x=PCA1_roots, y=prop01*100))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonizacion x Soil PCA", 
       x="Root Size (PCA1)", y="% Colonizacion",
       color="Area", fill="Area")


ggplot(Sum.Muestreo_UR1, aes(x=PCA1_roots, y=prop01*100, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonizacion x Soil PCA", 
       x="Root Size (PCA1)", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

# PCA2 - FORM
col_root2<-lm(PCA2_roots~prop01, data=Sum.Muestreo_UR1)
plot(col_root2,2)
summary(col_root2)

col_root2a<-lm(PCA2_roots~prop01*area, data=Sum.Muestreo_UR1)
Anova(col_root2a, type="2", test="F")


ggplot(Sum.Muestreo_UR1, aes(x=PCA2_roots, y=prop01*100))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonizacion x Soil PCA", 
       x="Root Form (PCA2)", y="% Colonizacion",
       color="Area", fill="Area")


ggplot(Sum.Muestreo_UR1, aes(x=PCA2_roots, y=prop01*100, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonizacion x Soil PCA", 
       x="Root Form (PCA2)", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)




##########################################
# Multiple interaction
##########################################
# no se uso

all_soil01<-glm(prop01~N*P*K*pH, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
summary(all_soil01)
Anova(all_soil01, type="3")

all_soil01<-glm(prop01~PCA1_roots+N*P*K*pH, data=Sum.Muestreo_UR1_pca,
                weights = total_field, family=quasibinomial)
summary(all_soil01)
Anova(all_soil01, type="3")

all_soil01ar<-glm(prop01~N*P*K*pH*area, data=Sum.Muestreo_UR1,
                  weights = total_field, family=quasibinomial)
summary(all_soil01ar)
Anova(all_soil01ar, type="3")


all_soil01ar<-glm(prop01~PCA1_roots+N*P*K*pH*area, data=Sum.Muestreo_UR1,
                  weights = total_field, family=binomial)
summary(all_soil01ar)
Anova(all_soil01ar, type="3")



### DOUBLE

NP_soil01<-glm(prop01~N*P, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
NP_soil01<-glm(prop01~N*P*area, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
summary(NP_soil01)
Anova(NP_soil01, type="3")


PK_soil01<-glm(prop01~P*K, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
PK_soil01<-glm(prop01~P*K*area, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
summary(PK_soil01)
Anova(PK_soil01, type="3")


NK_soil01<-glm(prop01~N*K, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
NK_soil01<-glm(prop01~N*K*area, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
summary(NK_soil01)
Anova(NK_soil01, type="3")



NPK_soil01<-glm(prop01~N*P*K, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
NPK_soil01<-glm(prop01~N*P*K*area, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
NPK_soil01<-glm(prop01~PCA1_roots+N*P*K*area, data=Sum.Muestreo_UR1,
                weights = total_field, family=quasibinomial)
summary(NPK_soil01)
Anova(NPK_soil01, type="3")



pH_soil01<-glm(prop01~N*pH, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
pH_soil01<-glm(prop01~P*pH, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
pH_soil01<-glm(prop01~K*pH, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
pH_soil01<-glm(prop01~P*K*pH*area, data=Sum.Muestreo_UR1,
               weights = total_field, family=quasibinomial)
pH_soil01<-lm(prop01~PCA1_roots+N+P+K+pH+area, data=Sum.Muestreo_UR1)
pH_soil01<-lm(prop01~PCA1_roots+N*area+P*area+K*area+pH*area,
              data=Sum.Muestreo_UR1)
summary(pH_soil01)
plot(pH_soil01,2)
Anova(pH_soil01, type="3")


AIC(pH_soil01)

step(pH_soil01,directions="both")
step_all<-stepAIC(pH_soil01, directions="both")

step_all$anova


Sum.Muestreo_UR1$P_N<-Sum.Muestreo_UR1$P/Sum.Muestreo_UR1$N
Sum.Muestreo_UR1$K_N<-Sum.Muestreo_UR1$K/Sum.Muestreo_UR1$N

pH_soil01<-lm(prop01~PCA1_roots+P_N*K_N*pH*area, data=Sum.Muestreo_UR1)
Anova(pH_soil01, type="3")




#########################################################
#########################################################
###                                                   ###
###               INTERACTION ANALYSES                ###
###                                                   ###
#########################################################
#########################################################

####################################
####### COLONIZATION VS ROOT #######
####################################
# NO ESTA EN RESULTADOS

#Linear Model
col_vol <- lm(asin(sqrt(prop01))~volumen_ml, data=Muestreo_UR)
summary(col_vol)
anova(col_vol)

plot(col_vol, 2)# the distribution is NOT normal
ols_test_normality(col_vol)#TESTING RESIDUALS


#Generalized Linear regression
col_vol<-glm(prop01~volumen_ml, 
             data=Muestreo_UR, 
             family=quasibinomial, #dispersion = 26
             weights = total_field,
             subset = !is.na(volumen_ml))
summary(col_vol)
Anova(col_vol)#binomial is significant
              #quasibinomial is NOT significant


#USING BINOMIAL MAKES THESE SIGNIFICANT TOO
col_vol0<-lmer(asin(sqrt(prop01))~(1|site), 
                data=Muestreo_UR, 
                subset = !is.na(volumen_ml))

col_vol1<-lmer(asin(sqrt(prop01))~volumen_ml+(1|site), 
                data=Muestreo_UR, 
                subset = !is.na(volumen_ml))

col_vol1<-glmer(prop01~volumen_ml+(1|site), 
             data=Muestreo_UR, 
             family=binomial, #dispersion = 26
             weights = total_field,
             subset = !is.na(volumen_ml))

anova(col_vol0,col_vol1)
ranova(col_vol1)
anova(col_vol1)

summary(col_vol1)



ggplot(Muestreo_UR, aes(x=volumen_ml, y=prop_colonization))+
  geom_point()+theme_classic()+
  labs(       x="Volume (ml)", y="% Colonization")+
  geom_smooth(method="lm", se=F, color="black")




########  BY AREA
#volumen_ml
#long_root
#no_roots

#w_total
#long_root



Muestreo_UR2<-Muestreo_UR
for (i in seq(203)){
  if (is.na(Muestreo_UR2$w_total[i])){
    Muestreo_UR2$w_total[i]<-
      Muestreo_UR2$w_above[i]+Muestreo_UR2$r_wet[i]
  }
  else{
    Muestreo_UR2$w_total[i]<-NA
  }
}

#regresion linea
col_vola0<-lm(log(w_total)~prop01*area, 
              data=Muestreo_UR2)

plot(col_vola0, 2)# the distribution is NOT normal
ols_test_normality(col_vola0)#TESTING RESIDUALS

Anova(col_vola0, type = "2")#Interaction effect NOT significat


col_vola0<-lmer(log(long_root)~prop01*area+(1|site), 
              data=Muestreo_UR2)
ranova(col_vola0)
Anova(col_vola0, type = "2")#Interaction effect NOT significat: Vol Root
                            #Interaction effect NOT significat: Long Root
                            #Interaction effect NOT significat: No Root

#Interaction effect NOT significat: Total weight
#Interaction effect NOT significat: Plant height


#no_branch
col_vola0<-glm(no_branch~prop01*area, 
              data=Muestreo_UR2, 
              family=poisson)
summary(col_vola0) #dispersion = 0.89
Anova(col_vola0, type = "2")#Interaction effect NOT significat


col_vola0<-glmer(no_branch~prop01*area+(1|site), 
                data=Muestreo_UR2, 
                family=poisson)
ranova(col_vola0)
Anova(col_vola0, type = "2")

#Interaction effect NOT significat: Num of Branches



#Graficas
ggplot(Muestreo_UR, aes(x=volumen_ml, y=prop01*100, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(x="Volumen (ml)", y="% Colonizacion",
       color="Ambient", fill="Ambient")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)




########################################
####### COLONIZATION VS DENSITY  #######
########################################
#ESTA EN RESULTADOS, Tabla S5


#Linear Model
col_dense<-lm(log(dense)~prop01, 
             data=Sum.Muestreo_UR1)
col_dense<-lm(prop01~dense, 
              data=Sum.Muestreo_UR1)
summary(col_dense)
anova(col_dense)


plot(col_dense,2)
ols_test_normality(col_dense)#TESTING RESIDUALS


#Generalized lineal regresion
col_dense<-glm(dense~prop01, 
              data=Sum.Muestreo_UR1, 
              family=quasipoisson) 
summary(col_dense)
Anova(col_dense) #results not so different from lm


#BY AREA
col_densea<-lm(log(dense)~prop01*area, 
              data=Sum.Muestreo_UR1)
col_densea<-lm(prop01~dense*area, 
              data=Sum.Muestreo_UR1)
plot(col_densea,2)
ols_test_normality(col_densea)

summary(col_densea)
anova(col_densea)
Anova(col_densea, type=2, test="F") #ESTA EN RESULTADOS, Tabla S5

#GLM
col_densea0<-glm(dense~prop01, 
                data=Sum.Muestreo_UR1, 
                family=quasipoisson)
col_densea1<-glm(dense~prop01*area, 
          data=Sum.Muestreo_UR1, 
          family=quasipoisson)


anova(col_densea0,col_densea1)
Anova(col_densea1, type=3, test="LR")


#Graficas
ggplot(Sum.Muestreo_UR1, aes(x=dense, y=prop01))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonization vs density", 
       x="No. spores / 100g", y="% Colonization")

ggplot(Sum.Muestreo_UR1, aes(x=dense, y=prop01*100, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonization vs density per area",
       x="No. spores / 100g", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

  

#vol vs soil

volN<-lm(volumen_ml~N, 
         data=Sum.Muestreo_UR1) #normal
volP<-lm(volumen_ml~P, 
         data=Sum.Muestreo_UR1) #normal
volK<-lm(volumen_ml~K, 
         data=Sum.Muestreo_UR1) #normal
volpH<-lm(volumen_ml~pH, 
          data=Sum.Muestreo_UR1)#normal

summary(volN) #significant
summary(volP) #NOT significant
summary(volK) #significant
summary(volpH)#NOT significant

###by area
volNa<-lm(volumen_ml~N*area, 
         data=Sum.Muestreo_UR1) #normal
volPa<-lm(volumen_ml~P*area, 
         data=Sum.Muestreo_UR1) #normal
volKa<-lm(volumen_ml~K*area, 
         data=Sum.Muestreo_UR1) #normal
volpHa<-lm(volumen_ml~pH*area, 
          data=Sum.Muestreo_UR1)#normal

Anova(volNa, type="2") #NOT significant
Anova(volPa, type="2") #NOT significant / area*
Anova(volKa, type="2") #NOT significant
Anova(volpHa, type="2")#YES !! significant

#root length vs soil

lrootN<-lm(log(long_root)~N, 
         data=Sum.Muestreo_UR1) #log normal
lrootP<-lm(long_root~P, 
         data=Sum.Muestreo_UR1) #normal
lrootK<-lm(log(long_root)~K, 
         data=Sum.Muestreo_UR1) #log normal
lrootpH<-lm(log(long_root)~pH, 
          data=Sum.Muestreo_UR1)#log normal

summary(lrootN) #significant
summary(lrootP) #significant
summary(lrootK) #significant
summary(lrootpH)#significant


###by area
lrootNa<-lm(long_root~N*area, 
           data=Sum.Muestreo_UR1) #normal
lrootPa<-lm(long_root~P*area, 
           data=Sum.Muestreo_UR1) #normal
lrootKa<-lm(long_root~K*area, 
           data=Sum.Muestreo_UR1) #normal
lrootpHa<-lm(long_root~pH*area, 
            data=Sum.Muestreo_UR1)#normal???

plot(lrootpHa, 2)#data is normally distributed
ols_test_normality(lrootpHa)#TESTING RESIDUALS

Anova(lrootNa, type="2") #NOT significant / area***
Anova(lrootPa, type="2") #NOT significant / area***
Anova(lrootKa, type="2") #NOT significant / P*, area***
Anova(lrootpHa, type="2")#NOT significant / area***

#TOTAL WEIGHT => NONE SIGNIFICANT
#PLANT HEIGHT => NONE SIGNIFICANT

  
##########################################
####### COLONIZACION VS SOIL CHAR  #######
##########################################
#ESTA EN RESULTADOS, Tabla S4


#Linear regression
colN<-lm(prop01~N, 
         data=Sum.Muestreo_UR1)
colP<-lm(prop01~P, 
         data=Sum.Muestreo_UR1)
colK<-lm(log(prop01)~K, 
         data=Sum.Muestreo_UR1)
colpH<-lm(log(prop01)~pH, 
          data=Sum.Muestreo_UR1)
#N
plot(colN, 2)#data is normally distributed
ols_test_normality(colN)#TESTING RESIDUALS
  #statistical test correspond with visual test

#P
plot(colP, 2)#data is normally distributed
ols_test_normality(colP)#TESTING RESIDUALS
  #statistical test correspond with visual test

#K
plot(colK, 2)#data is NOT normally distributed
ols_test_normality(colK)#TESTING RESIDUALS
  #statistical test DOES NOT correspond with visual test

#pH
plot(colpH, 2)#data is NOT normally distributed
ols_test_normality(colpH)#TESTING RESIDUALS
  #statistical test DOES NOT correspond with visual test



#Generalized linear regression

  #colN<-glm(prop01~N, 
  #          data=Sum.Muestreo_UR1, 
  #          family=quasibinomial, weights = total_field)
  #colP<-glm(prop01~P, 
  #          data=Sum.Muestreo_UR1, 
  #          family=quasibinomial, weights = total_field)

colK<-glm(prop01~K, 
          data=Sum.Muestreo_UR1, 
          family=quasibinomial, weights = total_field)
colpH<-glm(prop01~pH, 
          data=Sum.Muestreo_UR1, 
          family=quasibinomial, weights = total_field)

summary(colN)
summary(colP)
summary(colK)
summary(colpH)


#Graficas
col_N<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=long_root))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonizacion por Nitrogeno", 
       x="N (%)", y="% Colonizacion")

col_P<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=long_root))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonizacion por Fosforo", 
       x="P (mg/kg)", y="% Colonizacion")

col_K<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=long_root))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonizacion por Potasio", 
       x="K (mg/kg)", y="% Colonizacion")

col_pH<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=long_root))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonization vs pH", 
       x="pH", y="% Colonizacion")


grid.arrange(col_N, col_P, col_K, col_pH,
             widths= c(1,1,1,1), 
             layout_matrix= rbind(c(1,2,3,4)))

grid.arrange(col_N, col_P, 
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))

########  BY AREA ######## 

colNar<-lm(prop01~N*area, 
           data=Sum.Muestreo_UR1)
colPar<-lm(prop01~P*area, 
           data=Sum.Muestreo_UR1) #log?
colKar<-lm(log(prop01)~K*area, 
           data=Sum.Muestreo_UR1)
colpHar<-lm(prop01~pH*area, 
            data=Sum.Muestreo_UR1)

#N
plot(colNar, 2)#data is normally distributed
ols_test_normality(colNar)#TESTING RESIDUALS

#P
plot(colPar, 2)#data is normally distributed ???
ols_test_normality(colPar)#TESTING RESIDUALS ?

#K
plot(colKar, 2)#data is normally distributed 
ols_test_normality(colKar)#TESTING RESIDUALS

#pH
plot(colpHar, 2)#data is NOT normally distributed ???
ols_test_normality(colpHar)#TESTING RESIDUALS
  #statistical test show: normal distribution


anova(colNar)
anova(colPar)
anova(colKar)
anova(colpHar)

Anova(colNar, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(colPar, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(colKar, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(colpHar, type=3, test="F")#ESTA EN RESULTADOS, Tabla S4

summary(colKar)
beta.coef(colKar)


#Generalized linear regression
colNar<-glm(prop01~N*area, 
            data=Sum.Muestreo_UR1, 
            family=quasibinomial, weights = total_field)
colPar<-glm(prop01~P*area, 
            data=Sum.Muestreo_UR1, 
            family=quasibinomial, weights = total_field)
colKar<-glm(prop01~K*area, 
            data=Sum.Muestreo_UR1, 
            family=quasibinomial, weights = total_field)
colpHar<-glm(prop01~pH*area, 
            data=Sum.Muestreo_UR1, 
            family=quasibinomial, weights = total_field)


summary(colKar)
summary(aov(colKar))
Anova(colNar, type="3", test="F")
Anova(colPar, type="3", test="F")
Anova(colKar, type="3", test="F")
Anova(colpHar, type="3", test="F")

#Graficas
colNA1<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=prop01, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonization vs Nitrogen by Area",
       x="N (%)", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

colPA1<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=prop01, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              se = FALSE, color="black")+
  
  geom_smooth(method="lm", se=F)+
  labs(title="Colonization vs Phosphorus by Area",
       x="P (mg/kg)", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

colKA1<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=prop01, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonization vs Potasium by Area",
       x="K (mg/kg)", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

colpHA1<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=prop01, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonization vs pH by Area",
       x="pH", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)


grid.arrange(colNA1, colPA1, colKA1, colpHA1,
             widths= c(1,1,1,1), 
             layout_matrix= rbind(c(1,2,3,4)))




################
### SOIL PCA ###
################
# NO SE USO EN RESULTADOS

colPCAsoil<-lm(prop01~PCAsoil1, 
               data=Sum.Muestreo_UR2)

plot(colPCAsoil, 2)#data is normally distributed
ols_test_normality(colPCAsoil)#TESTING RESIDUALS
summary(colPCAsoil)


ggplot(Sum.Muestreo_UR1, aes(x=PCAsoil1, y=prop01*100))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonizacion x Soil PCA", 
       x="Soil PCA Dim1", y="% Colonizacion")






### BY AREA
colPCAsoilar<-lm(prop01~PCAsoil1*area, 
                 data=Sum.Muestreo_UR2)

plot(colPCAsoilar, 2)#data is normally distributed
ols_test_normality(colPCAsoilar)#TESTING RESIDUALS

summary(colPCAsoilar)
anova(colPCAsoilar)


ggplot(Sum.Muestreo_UR1, aes(x=PCAsoil1, y=prop01*100, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonizacion x Soil PCA", 
       x="Soil PCA Dim1", y="% Colonizacion",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Profundo"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Profundo"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)










#########################################
####### Nitrogeno vs Volumen Raiz #######
#########################################
#NO SE USO EN RESULTADOS

#Linear regression
Vol_N<-lm(volumen_ml~N, 
         data=Sum.Muestreo_UR1)
Vol_P<-lm(volumen_ml~P, 
         data=Sum.Muestreo_UR1)
Vol_K<-lm(volumen_ml~K, 
         data=Sum.Muestreo_UR1)
Vol_pH<-lm(volumen_ml~pH, 
          data=Sum.Muestreo_UR1)
#N
plot(Vol_N, 2)#data is normally distributed
ols_test_normality(Vol_N)#TESTING RESIDUALS

#P
plot(Vol_P, 2)#data is normally distributed
ols_test_normality(Vol_P)#TESTING RESIDUALS

#K
plot(Vol_K, 2)#data is normally distributed
ols_test_normality(Vol_K)#TESTING RESIDUALS

#pH
plot(Vol_pH, 2)#data is normally distributed
ols_test_normality(Vol_pH)#TESTING RESIDUALS

summary(Vol_N)#YES**
summary(Vol_P)#NOT
summary(Vol_K)#YES*
summary(Vol_pH)#NOT


#Grafico
NVol<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=volumen_ml))+
    geom_point()+theme_classic()+
    geom_smooth(method = "lm", se=F, color="black")+
    labs(title=" ", 
         x="N (%)", y="Volume (ml)")

PVol<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=volumen_ml))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="P (mg/kg)", y="Volume (ml)")

KVol<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=volumen_ml))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="K (mg/kg)", y="Volume (ml)")

pHVol<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=volumen_ml))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="pH", y="Volume (ml)")


grid.arrange(NVol, PVol, KVol, pHVol,
             widths= c(1,1,1,1), 
             layout_matrix= rbind(c(1,2,3,4)))





########  BY AREA ######## 
#volumen_ml
#long_root
#no_roots

#Regresion lineal
Vol_N_lma<-lm(no_roots ~ N*area, data = Sum.Muestreo_UR1)
Vol_P_lma<-lm(log(no_roots) ~ P*area, data = Sum.Muestreo_UR1)
Vol_K_lma<-lm(no_roots ~ K*area, data = Sum.Muestreo_UR1)
Vol_pH_lma<-lm(log(no_roots) ~ pH*area, data = Sum.Muestreo_UR1)

#w_total
#long_plant
#no_branch
Sum.Muestreo_UR3<-Sum.Muestreo_UR1
for (i in seq(40)){
  if (is.na(Sum.Muestreo_UR3$w_total[i])){
    Sum.Muestreo_UR3$w_total[i]<-
      Sum.Muestreo_UR3$w_above[i]+Sum.Muestreo_UR3$r_wet[i]
  }
  else{
    Sum.Muestreo_UR3$w_total[i]<-NA
  }
}


Vol_N_lma<-lm(log(no_branch) ~ N*area, data = Sum.Muestreo_UR3)
Vol_P_lma<-lm(log(no_branch) ~ P*area, data = Sum.Muestreo_UR3)
Vol_K_lma<-lm(log(no_branch) ~ K*area, data = Sum.Muestreo_UR3)
Vol_pH_lma<-lm(no_branch ~ pH*area, data = Sum.Muestreo_UR3)

#N
plot(Vol_N_lma, 2)#data is normally distributed
ols_test_normality(Vol_N_lma)#TESTING RESIDUALS

#P
plot(Vol_P_lma, 2)#data is normally distributed
ols_test_normality(Vol_P_lma)#TESTING RESIDUALS

#K
plot(Vol_K_lma, 2)#data is normally distributed
ols_test_normality(Vol_K_lma)#TESTING RESIDUALS

#pH
plot(Vol_pH_lma, 2)#data is normally distributed
ols_test_normality(Vol_pH_lma)#TESTING RESIDUALS

summary(Vol_N_lma)
summary(Vol_P_lma)
summary(Vol_K_lma)
summary(Vol_pH_lma)

Anova(Vol_N_lma, type="2")#nothing
Anova(Vol_P_lma, type="2")#nothing
Anova(Vol_K_lma, type="2")#nothing
Anova(Vol_pH_lma, type="2")#YES INTERACTION





#Graphics
NVola<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=volumen_ml, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="N (%)", y="Volume (ml)")+
  scale_color_manual(values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

PVola<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=volumen_ml, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="P (mg/kg)", y="Volume (ml)")+
  scale_color_manual(values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

KVola<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=volumen_ml, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="K (mg/kg)", y="Volume (ml)")+
  scale_color_manual(values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

pHVola<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=volumen_ml, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="pH", y="Volume (ml)")+
  scale_color_manual(values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)




grid.arrange(NVola, PVola, KVola, pHVola,
             widths= c(1,1,1,1), 
             layout_matrix= rbind(c(1,2,3,4)))





#########################################
####### SUELO vs Densidad Esporas #######
#########################################
#ESTA EN RESULTADOS, Tabla S4

#Linear regression
Spore_N<-lm(dense ~ N, data = Sum.Muestreo_UR1)
Spore_P<-lm(dense ~ P, data = Sum.Muestreo_UR1)
Spore_K<-lm(dense ~ K, data = Sum.Muestreo_UR1)
Spore_pH<-lm(dense ~ pH, data = Sum.Muestreo_UR1)


#N
plot(Spore_N, 2)#data is NOT normally distributed
ols_test_normality(Spore_N)#TESTING RESIDUALS

#P
plot(Spore_P, 2)#data is NOT normally distributed
ols_test_normality(Spore_P)#TESTING RESIDUALS

#K
plot(Spore_K, 2)#data is NOT normally distributed
ols_test_normality(Spore_K)#TESTING RESIDUALS

#pH
plot(Spore_pH, 2)#data is NOT normally distributed
ols_test_normality(Spore_pH)#TESTING RESIDUALS

Spore_N<-lm(log(dense) ~ N, data = Sum.Muestreo_UR1)
Spore_P<-lm(log(dense) ~ P, data = Sum.Muestreo_UR1)
Spore_K<-lm(log(dense) ~ K, data = Sum.Muestreo_UR1)
Spore_pH<-lm(log(dense) ~ pH, data = Sum.Muestreo_UR1)#norm


#Regresion lineal
Spore_N<-glm(dense ~ N, data = Sum.Muestreo_UR1, family=poisson)
Spore_P<-glm(dense ~ P, data = Sum.Muestreo_UR1, family=poisson)
Spore_K<-glm(dense ~ K, data = Sum.Muestreo_UR1, family=poisson)
#Spore_pH<-glm(dense ~ pH, data = Sum.Muestreo_UR1, family=poisson)

summary (Spore_N, corr=T)#nothing
summary (Spore_P, corr=T)#nothing
summary (Spore_K, corr=T)#nothing
summary (Spore_pH, corr=T)#nothing





#Grafico
pHD<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=dense))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="pH", y="No. de esporas / 100g")+
  guides(fill=F)

KD<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=dense))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="Concentracion K(mg/kg)", y="No. de esporas / 100g")+
  guides(fill=F)

PD<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=dense))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="Concentracion P(mg/kg)", y="No. de esporas / 100g")+
  guides(fill=F)

ND<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=dense))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="Concentracion N(%)", y="No. de esporas / 100g")+
  guides(fill=F)




grid.arrange(pHD, KD, PD, ND,
                   widths= c(1,1,1,1), 
                   layout_matrix= rbind(c(4,3,2,1)))




####### BY AREA
#ESTA EN RESULTADOS, Tabla S4

#Linear regression
Spore_NAr<-lm(log(dense) ~ N*area, data = Sum.Muestreo_UR1)
Spore_PAr<-lm(log(dense) ~ P*area, data = Sum.Muestreo_UR1)
Spore_KAr<-lm(dense ~ K*area, data = Sum.Muestreo_UR1)
Spore_pHAr<-lm(log(dense) ~ pH*area, data = Sum.Muestreo_UR1)


#N
plot(Spore_NAr, 2)#data is NOT normally distributed
ols_test_normality(Spore_NAr)#TESTING RESIDUALS

#P
plot(Spore_PAr, 2)#data is NOT normally distributed
ols_test_normality(Spore_PAr)#TESTING RESIDUALS

#K
plot(Spore_KAr, 2)#data is normally distributed
ols_test_normality(Spore_KAr)#TESTING RESIDUALS

#pH
plot(Spore_pHAr, 2)#data is NOT normally distributed
ols_test_normality(Spore_pHAr)#TESTING RESIDUALS


anova(Spore_NAr)#nothing
anova(Spore_PAr)#nothing
anova(Spore_KAr)#nothing
anova(Spore_pHAr)#nothing

Anova(Spore_NAr, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(Spore_PAr, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(Spore_KAr, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(Spore_pHAr, type=3, test="F")#ESTA EN RESULTADOS, Tabla S4


#Generalized linear regression (part1)
Spore_<-glm(dense ~ 1, data = Sum.Muestreo_UR1, family=poisson)
Spore_N<-glm(dense ~ N, data = Sum.Muestreo_UR1, family=poisson)
anova(Spore_, Spore_N)

Spore_NA<-glm(dense ~ N+area, data = Sum.Muestreo_UR1, family=poisson)
anova(Spore_N, Spore_NA)

Spore_NAr<-glm(dense ~ N*area, data = Sum.Muestreo_UR1, family = quasipoisson)
summary(Spore_NAr)
anova(Spore_NA, Spore_NAr)


#Generalized linear regression (part 2 - analyses)
Spore_NAr1<-glm(dense ~ N*area, data = Sum.Muestreo_UR1, family = quasipoisson)
Spore_PAr1<-glm(dense ~ P*area, data = Sum.Muestreo_UR1, family = quasipoisson)
Spore_KAr1<-glm(dense ~ K*area, data = Sum.Muestreo_UR1, family = quasipoisson)
Spore_pHAr1<-glm(dense ~ pH*area, data = Sum.Muestreo_UR1, family = quasipoisson)
summary(Spore_NAr1)
summary(Spore_PAr1)
summary(Spore_KAr1)
summary(Spore_pHAr1)


Anova(Spore_NAr1 , type="3", test="LR")#  N + Area + Interaction
Anova(Spore_PAr1 , type="3", test="LR")#  ! nothing
Anova(Spore_KAr1 , type="3", test="LR")#  K + Aarea + Interaction
Anova(Spore_pHAr1, type="3", test="LR")#  Area + Interaction




#Grafico
pHa<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=dense, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="pH", y="No. of spores / 100g",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

Ka<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=dense, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="K concentration (mg/kg)", y="No. of spores / 100g",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
guides(color=F, fill=F)

Pa<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=dense, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="P concentration (mg/kg)", y="No. of spores / 100g",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
guides(color=F, fill=F)

Na<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=dense, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="N concentration N(%)", y="No. of spores / 100g",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
guides(color=F, fill=F)



print(grid.arrange(pHa, Ka, Pa, Na,
                   widths= c(1,1,1,1), 
                   layout_matrix= rbind(c(4,3,2,1))))







#########################################
#######       VS DIVERSITY       #######
#########################################
#ESTA EN RESULTADOS, Tabla S4 y S5

##### DIVERSIDAD vs COLINIZATION

#Linear regression
col_div<-lm(log(prop01)~SHANNON, 
             data=Sum.Muestreo_UR1)

plot(col_div, 2)#data is normally distributed
ols_test_normality(col_div)#TESTING RESIDUALS

summary(col_div)


### BY AREA
#ESTA EN RESULTADOS, Tabla S5

#Linear regression
col_diva<-lm(prop01~SHANNON*area, 
             data=Sum.Muestreo_UR1)

plot(col_diva, 2)#data is normally distributed
ols_test_normality(col_diva)#TESTING RESIDUALS

summary(col_diva)
anova(col_diva)
Anova(col_diva, type=2, test="F") #ESTA EN RESULTADOS, Tabla S5


#Graficas
ggplot(Sum.Muestreo_UR1, aes(x=SHANNON, y=prop01*100))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Colonization vs spore diversity", 
       x="Shannon Index", y="% Colonization")

ggplot(Sum.Muestreo_UR1, aes(x=SHANNON, y=prop01*100, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Colonization vs spore diversity per area",
       x="Shannon Index", y="% Colonizacion",
       color="Subambiente", fill="Area")+
  scale_color_manual(labels=c("Rural","Deed"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)




##### DIVERSIDAD vs DENSITY


#Linear regression
dense_div<-lm(dense~SHANNON, 
            data=Sum.Muestreo_UR1)

plot(dense_div, 2)#data is normally distributed
ols_test_normality(dense_div)#TESTING RESIDUALS

summary(dense_div)


### BY AREA

#Linear regression
dense_diva<-lm(log(dense)~SHANNON*area, 
             data=Sum.Muestreo_UR1)

plot(dense_diva, 2)#data is normally distributed
ols_test_normality(dense_diva)#TESTING RESIDUALS

summary(dense_diva)
anova(dense_diva)
Anova(dense_diva, type=2, test="F")


#Graficas
ggplot(Sum.Muestreo_UR1, aes(x=SHANNON, y=dense))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", se=F, color="black")+
  labs(title="Spore density vs spore diversity", 
       x="Shannon Index", y="Num. Spores / 100g")

ggplot(Sum.Muestreo_UR1, aes(x=SHANNON, y=dense, color=area))+
  geom_point(aes(fill=area),pch=21,size=2)+
  theme_classic()+
  geom_smooth(method="lm", se=F)+
  labs(title="Spore density vs spore diversity per area",
       x="Shannon Index", y="Num. Spores / 100g",
       color="Subambiente", fill="Area")+
  scale_color_manual(labels=c("Rural","Deed"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)



##### DIVERSIDAD vs NUTIENRES
#ESTA EN RESULTADOS, Tabla S4


#NITROGEN
Diverse_N<-lm(SHANNON~N, data = Sum.Muestreo_UR1)
plot(Diverse_N, 2)#data is normally distributed
ols_test_normality(Diverse_N)#TESTING RESIDUALS

summary(Diverse_N)#NOT significant


#PHOSPHORUS
Diverse_P<-lm(SHANNON~P, data = Sum.Muestreo_UR1)
plot(Diverse_P, 2)#data is normally distributed
ols_test_normality(Diverse_P)#TESTING RESIDUALS

summary(Diverse_P)#YES!


#POTASIUM
Diverse_K<-lm(SHANNON~K, data = Sum.Muestreo_UR1)
plot(Diverse_K, 2)#data is normally distributed
ols_test_normality(Diverse_K)#TESTING RESIDUALS

summary(Diverse_K)#NOT significant


#pH
Diverse_pH<-lm(SHANNON~pH, data = Sum.Muestreo_UR1)
plot(Diverse_pH, 2)#data is normally distributed
ols_test_normality(Diverse_pH)#TESTING RESIDUALS

summary(Diverse_pH)#NOT significant


#Grafico
pHDiv<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=SHANNON))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="pH", y="Shanon Index")+
  guides(fill=F)

KDiv<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=SHANNON))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="Concentracion K(mg/kg)", y="Shanon Index")+
  guides(fill=F)

PDiv<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=SHANNON))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="Concentracion P(mg/kg)", y="Shanon Index")+
  guides(fill=F)

NDiv<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=SHANNON))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F, color="black")+
  labs(title=" ", 
       x="Concentracion N(%)", y="Shanon Index")+
  guides(fill=F)




grid.arrange(pHDiv, KDiv, PDiv, NDiv,
             widths= c(1,1,1,1), 
             layout_matrix= rbind(c(4,3,2,1)))


#### POR AREA

#NITROGEN
Diverse_NAr<-lm(SHANNON~N*area, data = Sum.Muestreo_UR1)
plot(Diverse_NAr, 2)#data is normally distributed
ols_test_normality(Diverse_NAr)#TESTING RESIDUALS

summary(Diverse_NAr)
anova(Diverse_NAr)#nothing

#PHOSPHORUS
Diverse_PAr<-lm(SHANNON~P*area, data = Sum.Muestreo_UR1)
plot(Diverse_PAr, 2)#data is normally distributed
ols_test_normality(Diverse_PAr)#TESTING RESIDUALS

summary(Diverse_PAr)
anova(Diverse_PAr)#nothing

#POTASIUM
Diverse_KAr<-lm(SHANNON~K*area, data = Sum.Muestreo_UR1)
plot(Diverse_KAr, 2)#data is normally distributed
ols_test_normality(Diverse_KAr)#TESTING RESIDUALS

summary(Diverse_KAr)
anova(Diverse_KAr)#nothing

#pH
Diverse_pHAr<-lm(SHANNON~pH*area, data = Sum.Muestreo_UR1)
plot(Diverse_pHAr, 2)#data is normally distributed
ols_test_normality(Diverse_pHAr)#TESTING RESIDUALS

summary(Diverse_pHAr)
anova(Diverse_pHAr)#nothing


anova(Diverse_NAr)
anova(Diverse_PAr)
anova(Diverse_KAr)
anova(Diverse_pHAr)


Anova(Diverse_NAr, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(Diverse_PAr, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(Diverse_KAr, type=3, test="F") #ESTA EN RESULTADOS, Tabla S4
Anova(Diverse_pHAr, type=3, test="F")#ESTA EN RESULTADOS, Tabla S4

#Grafico
pHDiva<-ggplot(Sum.Muestreo_UR1, aes(x=pH, y=SHANNON, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="pH", y="Shannon Index",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

KDiva<-ggplot(Sum.Muestreo_UR1, aes(x=K, y=SHANNON, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="K concentration (mg/kg)", y="Shannon Index",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

PDiva<-ggplot(Sum.Muestreo_UR1, aes(x=P, y=SHANNON, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="P concentration (mg/kg)", y="Shannon Index",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)

NDiva<-ggplot(Sum.Muestreo_UR1, aes(x=N, y=SHANNON, color=area))+
  geom_point()+theme_classic()+
  geom_smooth(method = "lm", se=F)+
  labs(title=" ", 
       x="N concentration N(%)", y="Shannon Index",
       color="Area", fill="Area")+
  scale_color_manual(labels=c("Rural","Deep"), 
                     values=c("#00ba38","#f8766d"))+
  scale_fill_manual(labels=c("Rural","Deep"),
                    values=c("#00ba38","#f8766d"))+
  guides(color=F, fill=F)



print(grid.arrange(pHDiva, KDiva, PDiva, NDiva,
                   widths= c(1,1,1,1), 
                   layout_matrix= rbind(c(4,3,2,1))))





################################################
####                                        ####
####                                        ####
####          EFFECT SIZE DATASET           ####
####                                        ####
####                                        ####
################################################
# NO SE USO EN RESULTADOS


planta_01_m
plantpca
rootpca
root01pca
soilpca
soil_spore_m


cohen.d(pH~area, data=Sum.Muestreo_UR1, 
                  na.action = na.omit, 
                  hedges.correction=T)[["conf.int"]][[1]]

Sum.Muestreo_UR1 %>%  
  summarise(
    effsize=c(pH_effsize=cohen.d(pH~area, data=Sum.Muestreo_UR1, 
                       na.action = na.omit, 
                       hedges.correction=T)[[4]],
              K_effsize=cohen.d(K~area, data=Sum.Muestreo_UR1, 
                      na.action = na.omit, 
                      hedges.correction=T)[[4]],
              P_effsize=cohen.d(P~area, data=Sum.Muestreo_UR1, 
                      na.action = na.omit, 
                      hedges.correction=T)[[4]],
              N_effsize=cohen.d(N~area, data=Sum.Muestreo_UR1, 
                      na.action = na.omit, 
                      hedges.correction=T)[[4]],
              Ca_effsize=cohen.d(Ca~area, data=Sum.Muestreo_UR1, 
                       na.action = na.omit, 
                       hedges.correction=T)[[4]],
              Na_effsize=cohen.d(Na~area, data=Sum.Muestreo_UR1, 
                       na.action = na.omit, 
                       hedges.correction=T)[[4]]),
    
    magnitude=c(pH_effsize=as.character(cohen.d(pH~area, data=Sum.Muestreo_UR1, 
                                 na.action = na.omit, 
                                 hedges.correction=T)[["magnitude"]]),
              K_effsize=as.character(cohen.d(K~area, data=Sum.Muestreo_UR1, 
                                na.action = na.omit, 
                                hedges.correction=T)[["magnitude"]]),
              P_effsize=as.character(cohen.d(P~area, data=Sum.Muestreo_UR1, 
                                na.action = na.omit, 
                                hedges.correction=T)[["magnitude"]]),
              N_effsize=as.character(cohen.d(N~area, data=Sum.Muestreo_UR1, 
                                na.action = na.omit, 
                                hedges.correction=T)[["magnitude"]]),
              Ca_effsize=as.character(cohen.d(Ca~area, data=Sum.Muestreo_UR1, 
                                 na.action = na.omit, 
                                 hedges.correction=T)[["magnitude"]]),
              Na_effsize=as.character(cohen.d(Na~area, data=Sum.Muestreo_UR1, 
                                 na.action = na.omit, 
                                 hedges.correction=T)[["magnitude"]]))
  ) -> effsize_suelos

(effsize_suelos<-as.data.frame(effsize_suelos))

rownames(effsize_suelos) <- c("pH","K","P",
                             "N","Ca","Na")

Muestreo_UR %>%  
  summarise(
    effsize=c(TotalW=cohen.d(w_total~ambient, data=Muestreo_UR, 
                       na.action = na.omit, 
                       hedges.correction=T)[[4]],
              Hieght=cohen.d(long_plant~ambient, data=Muestreo_UR, 
                             na.action = na.omit, 
                             hedges.correction=T)[[4]],
              No.Branch=cohen.d(no_branch~ambient, data=Muestreo_UR, 
                                na.action = na.omit, 
                                hedges.correction=T)[[4]],
              Root_Vol=cohen.d(volumen_ml~ambient, data=Muestreo_UR, 
                               na.action = na.omit, 
                               hedges.correction=T)[[4]],
              Root_length=cohen.d(long_root~ambient, data=Muestreo_UR, 
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]],
              No.Root=cohen.d(no_roots~ambient, data=Muestreo_UR, 
                              na.action = na.omit, 
                              hedges.correction=T)[[4]]
              ),
    
    magnitude=c(TotalW=as.character(cohen.d(w_total~ambient, data=Muestreo_UR, 
                               na.action = na.omit, 
                               hedges.correction=T)[["magnitude"]]),
                Hieght=as.character(cohen.d(long_plant~ambient, data=Muestreo_UR, 
                               na.action = na.omit, 
                               hedges.correction=T)[["magnitude"]]),
                No.Branch=as.character(cohen.d(no_branch~ambient, data=Muestreo_UR, 
                                  na.action = na.omit, 
                                  hedges.correction=T)[["magnitude"]]),
                Root_Vol=as.character(cohen.d(volumen_ml~ambient, data=Muestreo_UR, 
                                 na.action = na.omit, 
                                 hedges.correction=T)[["magnitude"]]),
                Root_length=as.character(cohen.d(long_root~ambient, data=Muestreo_UR, 
                                    na.action = na.omit, 
                                    hedges.correction=T)[["magnitude"]]),
                No.Root=as.character(cohen.d(no_roots~ambient, data=Muestreo_UR, 
                                na.action = na.omit, 
                                hedges.correction=T)[["magnitude"]]))
  ) -> effsize_plant


(effsize_plant<-as.data.frame(effsize_plant))
rownames(effsize_plant) <- c("TotalW","Hieght","No.Branch",
                             "Root_Vol","Root_length","No.Root")


effsize_col<-data.frame(effsize=cohen.d(prop01~ambient, data=Muestreo_UR, 
                                na.action = na.omit, 
                                hedges.correction=F)[[3]],
          magnitude=as.character(cohen.d(prop01~ambient, data=Muestreo_UR, 
                                  na.action = na.omit, 
                                  hedges.correction=F)[["magnitude"]]))

Sum.Muestreo_UR1 %>%  
  summarise(
    effsize=c(Dense=cohen.d(dense~ambient, data=Sum.Muestreo_UR1, 
                      na.action = na.omit, 
                      hedges.correction=F)[[3]],
              Diversity=cohen.d(SHANNON~ambient, data=Sum.Muestreo_UR1, 
                                na.action = na.omit, 
                                hedges.correction=F)[[3]]),
    
    magnitude=c(Dense=as.character(cohen.d(dense~ambient, data=Sum.Muestreo_UR1, 
                              na.action = na.omit, 
                              hedges.correction=F)[["magnitude"]]),
                Diversity=as.character(cohen.d(SHANNON~ambient, data=Sum.Muestreo_UR1, 
                        na.action = na.omit, 
                        hedges.correction=F)[["magnitude"]]))
  ) -> effsize_spore


effsize_amf<-(rbind(effsize_col,effsize_spore))
rownames(effsize_amf)<-c("Colonization","Density","Diversity")

effsize_suelos
effsize_plant
effsize_amf

#### EMMEANS

df<-emmeans(#anova
  )
stuff<-eff_size(df, 
                  sigma = sqrt(mean(sigma(#anova
                    )^2)),   # RMS of sigma()
                  edf = df.residual(#anova
                    ))






treatment = rnorm(100,mean=10)
control = rnorm(100,mean=12)
d = (c(treatment,control))
f = rep(c("Treatment","Control"),each=100)
## compute Cohen's d
## treatment and control
cohen.d(treatment,control)
## data and factor
cohen.d(d,f)
## formula interface
cohen.d(d ~ f)
## compute Hedges' g
cohen.d(d,f,hedges.correction=TRUE)
