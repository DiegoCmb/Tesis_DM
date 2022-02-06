
#R Version: 4.0.3

####### PAQUETES Y BASE DE DATOS ########

########### Cargando todos los paquetes
library(tidyverse)
library(car)
library(agricolae)
library(grid)
library(gridExtra)
library(lmtest)
library(lme4)
library(magrittr)
library(ggeffects)
library(sjmisc)

library(MASS)
library(mvoutlier)
library(mvnormtest)
library(pastecs)
library(reshape)
library(reshape2)

############# Cargando las bases del projecto
Experim_UR<-read.csv("Database/Exp_Estable_UR_Merida_2020.csv")[,-1]
Experim_UR_FQ<-read.csv("Database/Exp_Estable_UR_Merida_2020(FQ).csv")[,-1]



# Convirtiendo los datos importantes a factores y numericos
Experim_UR$ORIGEN<-as.factor(Experim_UR$ORIGEN)
Experim_UR$VASO<-as.factor(Experim_UR$VASO)
Experim_UR$TIERRA<-as.factor(Experim_UR$TIERRA)
Experim_UR$ESTERIL<-as.factor(Experim_UR$ESTERIL)
Experim_UR$GERM<-as.numeric(Experim_UR$GERM)

str(Experim_UR)

Experim_UR_FQ$ORIGEN<-as.factor(Experim_UR_FQ$ORIGEN)
Experim_UR_FQ$VASO<-as.factor(Experim_UR_FQ$VASO)
Experim_UR_FQ$TIERRA<-as.factor(Experim_UR_FQ$TIERRA)
Experim_UR_FQ$ESTERIL<-as.factor(Experim_UR_FQ$ESTERIL)
Experim_UR_FQ$SEMILLA<-as.factor(Experim_UR_FQ$SEMILLA)

Experim_UR_FQ$GERM<-as.numeric(Experim_UR_FQ$GERM)

Experim_UR_FQ$PAR<-as.numeric(Experim_UR_FQ$PAR)
Experim_UR_FQ$PSII<-as.numeric(Experim_UR_FQ$PSII)
Experim_UR_FQ$ETR<-as.numeric(Experim_UR_FQ$ETR)
Experim_UR_FQ$qP<-as.numeric(Experim_UR_FQ$qP)
Experim_UR_FQ$qN<-as.numeric(Experim_UR_FQ$qN)
Experim_UR_FQ$Fv_Fm<-as.numeric(Experim_UR_FQ$Fv_Fm)
Experim_UR_FQ$NPQ<-as.numeric(Experim_UR_FQ$NPQ)

############ FUNCTIONS ###########

#Funcion para obtener Error Estandar
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





###############################################################################


Experim_UR.NE<-Experim_UR[Experim_UR$ESTERIL=="NO",]
Experim_UR.E<-Experim_UR[Experim_UR$ESTERIL=="YES",]

#############################################
########## EXPERIMENTO GERMINACION ##########
#############################################

View(Experim_UR)
head(Experim_UR)
names(Experim_UR)
str(Experim_UR)




Experim_UR.NE1<-(gather(Experim_UR.NE, "g1","g2","g3","g4","g5",
                        "g6","g7","g8","g9","g10", 
                  key="seed", value="dias"))
Experim_UR.E1<-(gather(Experim_UR.E, "g1","g2","g3","g4","g5",
                        "g6","g7","g8","g9","g10", 
                        key="seed", value="dias"))


rta.NE_<-cbind(Experim_UR.NE$GERM_, Experim_UR.NE$weights- Experim_UR.NE$GERM_)
rta.E_<-cbind(Experim_UR.E$GERM_, Experim_UR.E$weights- Experim_UR.E$GERM_)


### Peso de semillas por Origen

peso_sem<-lm(PESO_m_s~ORIGEN, data=Experim_UR)
summary(peso_sem)
anova(peso_sem)

mean(Experim_UR$PESO_m_s[Experim_UR$ORIGEN == "Rural"])/
  mean(Experim_UR$PESO_m_s[Experim_UR$ORIGEN == "Urban"])

ggplot(Experim_UR, aes(x=ORIGEN, y=PESO_m_s, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Peso de semillas por ambiente", 
       x="Ambiente de origen", y="Peso (mg)")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))




### Numero de germinación


#Numero de Germinación por Origen de semilla
hist(Experim_UR$tiemp_g)

#binomial agregado
num_sem_orsem<-lm(GERM01~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR)



num_sem_orsem<-glm(GERM01~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR,
                   family=quasibinomial, weights = weights)


num_sem_orsem<-glm(rta~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR,
                   family=binomial(logit))

num_sem_orsem<-glm(rta_~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR,
                   family=binomial(logit))


num_sem_orsem<-glmer(rta~PESO_m_s+ORIGEN*TIERRA*ESTERIL+(1|SEMILLA), data=Experim_UR,
                     family=binomial(logit))
anova(num_sem_orsem, test="Chisq")

summary(num_sem_orsem)
Anova(num_sem_orsem)



prob_fit_sk1 <- ggpredict(num_sem_orsem, terms=c("TIERRA","ORIGEN","ESTERIL"))
sk1<-plot(prob_fit_sk1)




prob_fit_sk2 <- ggpredict(tiemp_sem_orsem, terms=c("TIERRA","ORIGEN","ESTERIL"))
sk2<-plot(prob_fit_sk2)




print(grid.arrange(sk1,sk2,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))






tiemp_sem_orsem0<-lm(tiemp_g~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR)
anova(tiemp_sem_orsem0)

tiemp_sem_orsem<-lm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1)


tiemp_sem_orsem<-glm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1,
                     family=Gamma)

anova(tiemp_sem_orsem, test="F")

Anova(tiemp_sem_orsem)
summary(aov(tiemp_sem_orsem))
plot(tiemp_sem_orsem)

summary.aov(tiemp_sem_orsem, split = list(TIERRA=list(Rural=1, Urban=2)))


TukeyHSD(aov(tiemp_sem_orsem))




ggplot(Experim_UR1, aes(x=VASO, y=dias, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Tiempo de germinación de semillas por tratamiento", 
       x="Tratamiento de tierra", y="Tiempo de germ (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))










por_sem<-data.frame(
  TIERRA=c("Rural","Urban"),
  meanT=c(mean(Experim_UR$GERM01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
          mean(Experim_UR$GERM01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  sdT=c(sd(Experim_UR$GERM01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        sd(Experim_UR$GERM01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  seT=c(se(Experim_UR$GERM01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        se(Experim_UR$GERM01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  ESTERIL=c("NO","YES"),
  meanE=c(mean(Experim_UR$GERM01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
          mean(Experim_UR$GERM01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  sdE=c(sd(Experim_UR$GERM01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        sd(Experim_UR$GERM01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  seE=c(se(Experim_UR$GERM01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        se(Experim_UR$GERM01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100))


por_sem_<-data.frame(
  TIERRA=c("Rural","Urban"),
  meanT=c(mean(Experim_UR$GERM01_[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
          mean(Experim_UR$GERM01_[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  sdT=c(sd(Experim_UR$GERM01_[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        sd(Experim_UR$GERM01_[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  seT=c(se(Experim_UR$GERM01_[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        se(Experim_UR$GERM01_[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  ESTERIL=c("NO","YES"),
  meanE=c(mean(Experim_UR$GERM01_[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
          mean(Experim_UR$GERM01_[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  sdE=c(sd(Experim_UR$GERM01_[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        sd(Experim_UR$GERM01_[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  seE=c(se(Experim_UR$GERM01_[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        se(Experim_UR$GERM01_[Experim_UR$ESTERIL=="YES"], na.rm = T)*100))




plot1a<-ggplot(por_sem_, 
               aes(x=TIERRA, y=meanT, fill=TIERRA))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=TIERRA, ymin=meanT-seT, ymax=meanT+seT, 
                    width=0.1))+
  labs(x="Origen de Tierra", y="% Germinación")+
  geom_text(aes(x=TIERRA,y=10+meanT+seT,
                label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))+
  guides(fill=F)



plot2a<-ggplot(por_sem, 
               aes(x=ESTERIL, y=meanE, fill=ESTERIL))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=ESTERIL, ymin=meanE-seE, ymax=meanE+seE, 
                    width=0.1))+
  labs(x="Estado de Esterilidad", y="% Germinación")+
  geom_text(aes(x=ESTERIL,y=+meanE+seE-10,
                label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("NO" = "No Esteril", "YES" = "Esteril"))+
  guides(fill=F)


Experim_UR %>% 
  group_by(VASO,ORIGEN) %>% 
  summarize(
    mean=mean(GERM01, na.rm = T)*100,
    sd=sd(GERM01,na.rm = T)*100,
    se=se(GERM01, na.rm = T)*100
  ) -> por_sem_completo


plot3a<-ggplot(por_sem_completo, 
               aes(x=VASO, y=mean, fill=ORIGEN))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black", position="dodge") +
  geom_errorbar(aes(x=VASO, ymin=mean-se, ymax=mean+se), 
                width=0.1, position=position_dodge(width = 0.9))+
  labs(x="Estado de Esterilidad", y="% Germinación")+
  scale_fill_manual(values=c("#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)


print(grid.arrange(plot1a,plot2a,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


print(grid.arrange(plot3a, plot7b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))  



#########################  #########################
######################### #########################
#########################  #########################
######################### #########################

tiemp_germ<-data.frame(
  TIERRA=c("Rural","Urban"),
  meanT=c(mean(Experim_UR$tiemp_g[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
          mean(Experim_UR$tiemp_g[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  sdT=c(sd(Experim_UR$tiemp_g[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        sd(Experim_UR$tiemp_g[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  seT=c(se(Experim_UR$tiemp_g[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        se(Experim_UR$tiemp_g[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  ESTERIL=c("NO","YES"),
  meanE=c(mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
          mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  sdE=c(sd(Experim_UR$tiemp_g[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        sd(Experim_UR$tiemp_g[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  seE=c(se(Experim_UR$tiemp_g[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        se(Experim_UR$tiemp_g[Experim_UR$ESTERIL=="YES"], na.rm = T)*100))


Experim_UR %>% 
  group_by(VASO,ORIGEN) %>% 
  summarise(
    mean=mean(tiemp_g, na.rm = T)*100,
    sd=sd(tiemp_g,na.rm = T)*100,
    se=se(tiemp_g, na.rm = T)*100
  ) -> tiemp_germ_completo




print(grid.arrange(plot2b,plot3b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


#########################  #########################
######################### #########################
#########################  #########################
######################### #########################

num_sem_orsem<-glm(GERM01~ORIGEN, data=Experim_UR,
                   family=binomial, weights=weights)
summary(num_sem_orsem)
anova(num_sem_orsem)


plot1a<-ggplot(Experim_UR, aes(x=ORIGEN, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de semilla", 
       x="Ambiente de origen de semilla", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Numero de Germinación por Origen de Tierra
num_sem_ortier<-lm(GERM~TIERRA, data=Experim_UR)
num_sem_orsem<-glm(GERM01~TIERRA, data=Experim_UR,
                   family=binomial, weights=weights)
summary(num_sem_ortier)
anova(num_sem_ortier)
Anova(num_sem_ortier)

mean(Experim_UR$GERM[Experim_UR$TIERRA == "Urban"])/
  mean(Experim_UR$GERM[Experim_UR$TIERRA == "Rural"])

plot2a<-ggplot(Experim_UR, aes(x=TIERRA, y=GERM, fill=TIERRA))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de tierra", 
       x="Ambiente de origen de tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Numero de Germinación por Esterilidad de Tierra
num_sem_estril<-lm(GERM~ESTERIL, data=Experim_UR)
summary(num_sem_estril)
anova(num_sem_estril)

mean(Experim_UR$GERM[Experim_UR$ESTERIL == "NO"])/
  mean(Experim_UR$GERM[Experim_UR$ESTERIL == "YES"])

mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL == "YES"], na.rm = T)/
  mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL == "NO"], na.rm = T)

plot3a<-ggplot(Experim_UR, aes(x=ESTERIL, y=GERM, fill=ESTERIL))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por esterilidad de tierra", 
       x="Esterilidad de Tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))





#Numero de Germinación por Combinacion de Tierra
num_sem_combtierr<-lm(GERM~TIERRA*ESTERIL, data=Experim_UR)
num_sem_combtierr<-lm(GERM~VASO, data=Experim_UR)
summary(num_sem_combtierr)
anova(num_sem_combtierr)

TukeyHSD(aov(num_sem_combtierr))

plot4a<-ggplot(Experim_UR, aes(x=VASO, y=GERM, fill=VASO))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por tratamiento de tierra", 
       x="Tratamiento de tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))


mean(Experim_UR$GERM[Experim_UR$VASO == "UX"], na.rm = T)/
  mean(Experim_UR$GERM[Experim_UR$VASO == "RX"], na.rm = T)

mean(Experim_UR$tiemp_g[Experim_UR$VASO == "RX"], na.rm = T)/
  mean(Experim_UR$tiemp_g[Experim_UR$VASO == "UX"], na.rm = T)






###Numero de Germinación por Combinacion de Tratamientos

#ORIGEN TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~TIERRA*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot5a<-ggplot(Experim_UR, aes(x=TIERRA, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de ambiente", 
       x="Ambiente de origen de tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))


#ESTERILIDAD DE  TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~ESTERIL*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot6a<-ggplot(Experim_UR, aes(x=ESTERIL, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por esterilidad de tierra", 
       x="Esterilidad de Tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))



#TRATAMIENTOS DE TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~VASO*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot7a<-ggplot(Experim_UR, aes(x=VASO, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por tratamiento", 
       x="Tratamiento de tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))










###########################
###Tiempo de germinación###
###########################


#Tiempo de Germinación por Origen de semilla
tiemp_sem_orsem<-lm(tiemp_g~ORIGEN, data=Experim_UR)
summary(tiemp_sem_orsem)
anova(tiemp_sem_orsem)

plot1b<-ggplot(Experim_UR, aes(x=ORIGEN, y=tiemp_g, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de semilla", y="Tiempo de germ (d????as)")+
  guides(fill=F)+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Tiempo de Germinación por Origen de Tierra
tiemp_sem_ortier<-lm(tiemp_g~TIERRA, data=Experim_UR)
summary(tiemp_sem_ortier)
anova(tiemp_sem_ortier)

plot2b<-ggplot(Experim_UR, aes(x=TIERRA, y=tiemp_g, fill=TIERRA))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de tierra", y="Tiempo de germ (d????as)")+
  guides(fill=F)+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Tiempo de Germinación por Esterilidad de Tierra
tiemp_sem_estril<-lm(tiemp_g~ESTERIL, data=Experim_UR)
summary(tiemp_sem_estril)
anova(tiemp_sem_estril)

plot3b<-ggplot(Experim_UR, aes(x=ESTERIL, y=tiemp_g, fill=ESTERIL))+
  geom_boxplot()+theme_classic()+
  labs(x="Esterilidad de Tierra", y="Tiempo de germ (d????as)")+
  guides(fill=F)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))



#Tiempo de Germinación por Combinacion de Tierra
tiemp_sem_combtierr<-lm(tiemp_g~TIERRA*ESTERIL, data=Experim_UR)
tiemp_sem_combtierr<-lm(tiemp_g~VASO, data=Experim_UR)
summary(tiemp_sem_combtierr)
anova(tiemp_sem_combtierr)

TukeyHSD(aov(tiemp_sem_combtierr))

plot4b<-ggplot(Experim_UR, aes(x=VASO, y=tiemp_g, fill=VASO))+
  geom_boxplot()+theme_classic()+
  labs(x="Tratamiento de tierra", y="Tiempo de germ (d????as)")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))






print(grid.arrange(plot1a,plot1b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


grid.arrange(plot2a,plot2b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))

grid.arrange(plot3a,plot3b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))


grid.arrange(plot4a,plot4b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))







###Tiempo de Germinación por Combinacion de Tratamientos

#ORIGEN TIERRA * ORIGEN SEMILLA
tiemp_sem_combtot<-lm(tiemp_g~TIERRA*ORIGEN, data=Experim_UR)
summary(tiemp_sem_combtot)
anova(tiemp_sem_combtot)

TukeyHSD(aov(tiemp_sem_combtot))

plot5b<-ggplot(Experim_UR, aes(x=TIERRA, y=tiemp_g, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de tierra", y="Tiempo de germ (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))


#ESTERILIDAD DE  TIERRA * ORIGEN SEMILLA
tiemp_sem_combtot<-lm(tiemp_g~ESTERIL*ORIGEN, data=Experim_UR)
summary(tiemp_sem_combtot)
anova(tiemp_sem_combtot)

TukeyHSD(aov(tiemp_sem_combtot))

plot6b<-ggplot(Experim_UR, aes(x=ESTERIL, y=tiemp_g, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Esterilidad de Tierra", y="Tiempo de germ (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))


#TRATAMIENTOS DE TIERRA * ORIGEN SEMILLA
tiemp_sem_combtot<-lm(tiemp_g~VASO*ORIGEN, data=Experim_UR)
summary(tiemp_sem_combtot)
anova(tiemp_sem_combtot)

TukeyHSD(aov(tiemp_sem_combtot))

plot7b<-ggplot(Experim_UR, aes(x=VASO, y=tiemp_g, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Tratamiento de tierra", y="Tiempo de germ (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)





grid.arrange(plot5a,plot5b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))



grid.arrange(plot6a,plot6b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))



grid.arrange(plot7a,plot7b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))













#####################
#     EMERGENCIA    #
#####################



###############################################################################
###############################################################################
###############################################################################





rtaE<-cbind(Experim_UR$EMER,Experim_UR$weights-Experim_UR$EMER)

Experim_UR$EMER01_<-Experim_UR$EMER/Experim_UR$GERM
rtaE_<-cbind(Experim_UR$EMER,Experim_UR$GERM-Experim_UR$EMER)
rtaE_c<-Experim_UR$GERM-Experim_UR$EMER


rtaE_c[rtaE_c>0]<-0
rtaE_c<-rtaE_c*(-1)
Experim_UR$GERM_<-Experim_UR$GERM+rtaE_c

Experim_UR$EMER01_>1
Experim_UR$EMER[58]
Experim_UR[58,13]<-10

TF<-rtaE_[,2]<0
rtaE_[rtaE_[,2]<0,2]<-0

### Numero de germinación


#Numero de Germinación por Origen de semilla
hist(Experim_UR$tiemp_e)

#binomial agregado

num_sem_emer <- lm(EMER01~ORIGEN, data=Experim_UR)
summary(num_sem_emer)


num_srsemem_o<-glm(EMER01~ORIGEN, data=Experim_UR,
                   family=binomial, weights = weights)


tiemp_sem_oremer<-glm(EMER01~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR,
                      family=binomial, weights = weights)

num_sem_emer<-glm(rtaE~ESTERIL, data=Experim_UR,
                  family=binomial(logit))


num_sem_oremer<-glmer(rtaE_~PESO_m_s+ORIGEN*TIERRA*ESTERIL+(1|SEMILLA), data=Experim_UR,
                      family=binomial(logit))
anova(num_sem_orsem, test="Chisq")

summary(num_sem_emer)
plot(EMER~ESTERIL, data=Experim_UR)
Anova(num_sem_orsem)



prob_fit_sk1e <- ggpredict(num_sem_oremer, terms=c("TIERRA","ORIGEN","ESTERIL"))
sk1e<-plot(prob_fit_sk1e)
summary(prob_fit_sk1e)

Experim_UR1e<-(gather(Experim_UR, "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10", 
                   key="seed", value="dias"))
tiemp_sem_oremer<-glm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1e,
                      family=Gamma)


prob_fit_sk2e <- ggpredict(tiemp_sem_oremer, terms=c("TIERRA","ORIGEN","ESTERIL"))
sk2e<-plot(prob_fit_sk2e)




print(grid.arrange(sk1e,sk2e,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))






tiemp_sem_orsem0<-lm(tiemp_g~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR)
anova(tiemp_sem_orsem0)

tiemp_sem_orsem<-lm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1)


tiemp_sem_orsem<-glm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1,
                     family=Gamma)

anova(tiemp_sem_orsem, test="F")

Anova(tiemp_sem_orsem)
summary(aov(tiemp_sem_orsem))
plot(tiemp_sem_orsem)

summary.aov(tiemp_sem_orsem, split = list(TIERRA=list(Rural=1, Urban=2)))


TukeyHSD(aov(tiemp_sem_orsem))




ggplot(Experim_UR1, aes(x=VASO, y=dias, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Tiempo de germinación de semillas por tratamiento", 
       x="Tratamiento de tierra", y="Tiempo de germ (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))










por_emer<-data.frame(
  TIERRA=c("Rural","Urban"),
  meanT=c(mean(Experim_UR$EMER01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
          mean(Experim_UR$EMER01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  sdT=c(sd(Experim_UR$EMER01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        sd(Experim_UR$EMER01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  seT=c(se(Experim_UR$EMER01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        se(Experim_UR$EMER01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  ESTERIL=c("NO","YES"),
  meanE=c(mean(Experim_UR$EMER01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
          mean(Experim_UR$EMER01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  sdE=c(sd(Experim_UR$EMER01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        sd(Experim_UR$EMER01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  seE=c(se(Experim_UR$EMER01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        se(Experim_UR$EMER01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100))




e_plot1a<-ggplot(por_emer, 
                 aes(x=TIERRA, y=meanT, fill=TIERRA))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=TIERRA, ymin=meanT-seT, ymax=meanT+seT, 
                    width=0.1))+
  labs(x="Origen de Tierra", y="% Emergencia")+
  geom_text(aes(x=TIERRA,y=10+meanT+seT,
                label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))+
  guides(fill=F)



e_plot2a<-ggplot(por_emer, 
                 aes(x=ESTERIL, y=meanE, fill=ESTERIL))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=ESTERIL, ymin=meanE-seE, ymax=meanE+seE, 
                    width=0.1))+
  labs(x="Estado de Esterilidad", y="% Emergencia")+
  geom_text(aes(x=ESTERIL,y=+meanE+seE-10,
                label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("NO" = "No Esteril", "YES" = "Esteril"))+
  guides(fill=F)




Experim_UR %>% 
  group_by(VASO,ORIGEN) %>% 
  summarise(
    mean=mean(EMER01, na.rm = T)*100,
    sd=sd(EMER01,na.rm = T)*100,
    se=se(EMER01, na.rm = T)*100
  ) -> por_emer_completo


e_plot3a<-ggplot(por_emer_completo, 
                 aes(x=VASO, y=mean, fill=ORIGEN))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black", position="dodge") +
  geom_errorbar(aes(x=VASO, ymin=mean-se, ymax=mean+se), 
                width=0.1, position=position_dodge(width = 0.9))+
  labs(x="Estado de Esterilidad", y="% Emergencia")+
  scale_fill_manual(values=c("#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)


print(grid.arrange(e_plot1a, e_plot2a,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


print(grid.arrange(e_plot3a, e_plot7b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))  



#########################  #########################
######################### #########################
#########################  #########################
######################### #########################

tiemp_emerg<-data.frame(
  TIERRA=c("Rural","Urban"),
  meanT=c(mean(Experim_UR$tiemp_e[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
          mean(Experim_UR$tiemp_e[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  sdT=c(sd(Experim_UR$tiemp_e[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        sd(Experim_UR$tiemp_e[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  seT=c(se(Experim_UR$tiemp_e[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        se(Experim_UR$tiemp_e[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  ESTERIL=c("NO","YES"),
  meanE=c(mean(Experim_UR$tiemp_e[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
          mean(Experim_UR$tiemp_e[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  sdE=c(sd(Experim_UR$tiemp_e[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        sd(Experim_UR$tiemp_e[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  seE=c(se(Experim_UR$tiemp_e[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        se(Experim_UR$tiemp_e[Experim_UR$ESTERIL=="YES"], na.rm = T)*100))


Experim_UR %>% 
  group_by(VASO,ORIGEN) %>% 
  summarise(
    mean=mean(tiemp_e, na.rm = T)*100,
    sd=sd(tiemp_e,na.rm = T)*100,
    se=se(tiemp_e, na.rm = T)*100
  ) -> tiemp_emer_completo




print(grid.arrange(e_plot2b,e_plot3b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


#########################  #########################
######################### #########################
#########################  #########################
######################### #########################

num_sem_orsem<-glm(GERM01~ORIGEN, data=Experim_UR,
                   family=binomial, weights=weights)
summary(num_sem_orsem)
anova(num_sem_orsem)


plot1a<-ggplot(Experim_UR, aes(x=ORIGEN, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de semilla", 
       x="Ambiente de origen de semilla", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Numero de Germinación por Origen de Tierra
num_sem_ortier<-lm(GERM~TIERRA, data=Experim_UR)
num_sem_orsem<-glm(GERM01~TIERRA, data=Experim_UR,
                   family=binomial, weights=weights)
summary(num_sem_ortier)
anova(num_sem_ortier)
Anova(num_sem_ortier)

mean(Experim_UR$GERM[Experim_UR$TIERRA == "Urban"])/
  mean(Experim_UR$GERM[Experim_UR$TIERRA == "Rural"])

plot2a<-ggplot(Experim_UR, aes(x=TIERRA, y=GERM, fill=TIERRA))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de tierra", 
       x="Ambiente de origen de tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Numero de Germinación por Esterilidad de Tierra
num_sem_estril<-lm(GERM~ESTERIL, data=Experim_UR)
summary(num_sem_estril)
anova(num_sem_estril)

mean(Experim_UR$GERM[Experim_UR$ESTERIL == "NO"])/
  mean(Experim_UR$GERM[Experim_UR$ESTERIL == "YES"])

mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL == "YES"], na.rm = T)/
  mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL == "NO"], na.rm = T)

plot3a<-ggplot(Experim_UR, aes(x=ESTERIL, y=GERM, fill=ESTERIL))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por esterilidad de tierra", 
       x="Esterilidad de Tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))





#Numero de Germinación por Combinacion de Tierra
num_sem_combtierr<-lm(GERM~TIERRA*ESTERIL, data=Experim_UR)
num_sem_combtierr<-lm(GERM~VASO, data=Experim_UR)
summary(num_sem_combtierr)
anova(num_sem_combtierr)

TukeyHSD(aov(num_sem_combtierr))

plot4a<-ggplot(Experim_UR, aes(x=VASO, y=GERM, fill=VASO))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por tratamiento de tierra", 
       x="Tratamiento de tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))


mean(Experim_UR$GERM[Experim_UR$VASO == "UX"], na.rm = T)/
  mean(Experim_UR$GERM[Experim_UR$VASO == "RX"], na.rm = T)

mean(Experim_UR$tiemp_g[Experim_UR$VASO == "RX"], na.rm = T)/
  mean(Experim_UR$tiemp_g[Experim_UR$VASO == "UX"], na.rm = T)






###Numero de Germinación por Combinacion de Tratamientos

#ORIGEN TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~TIERRA*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot5a<-ggplot(Experim_UR, aes(x=TIERRA, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de ambiente", 
       x="Ambiente de origen de tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))


#ESTERILIDAD DE  TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~ESTERIL*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot6a<-ggplot(Experim_UR, aes(x=ESTERIL, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por esterilidad de tierra", 
       x="Esterilidad de Tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))



#TRATAMIENTOS DE TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~VASO*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot7a<-ggplot(Experim_UR, aes(x=VASO, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por tratamiento", 
       x="Tratamiento de tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))










###########################
###Tiempo de germinación###
###########################


#Tiempo de Germinación por Origen de semilla
tiemp_sem_orsem<-lm(tiemp_g~ORIGEN, data=Experim_UR)
summary(tiemp_sem_orsem)
anova(tiemp_sem_orsem)

e_plot1b<-ggplot(Experim_UR, aes(x=ORIGEN, y=tiemp_e, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de semilla", y="Tiempo de emer (dias)")+
  guides(fill=F)+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Tiempo de Germinación por Origen de Tierra
tiemp_sem_ortier<-lm(tiemp_g~TIERRA, data=Experim_UR)
summary(tiemp_sem_ortier)
anova(tiemp_sem_ortier)

e_plot2b<-ggplot(Experim_UR, aes(x=TIERRA, y=tiemp_e, fill=TIERRA))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de tierra", y="Tiempo de emer (diaas)")+
  guides(fill=F)+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Tiempo de Germinación por Esterilidad de Tierra
tiemp_sem_estril<-lm(tiemp_g~ESTERIL, data=Experim_UR)
summary(tiemp_sem_estril)
anova(tiemp_sem_estril)

e_plot3b<-ggplot(Experim_UR, aes(x=ESTERIL, y=tiemp_e, fill=ESTERIL))+
  geom_boxplot()+theme_classic()+
  labs(x="Esterilidad de Tierra", y="Tiempo de emer (diaas)")+
  guides(fill=F)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))



#Tiempo de Germinación por Combinacion de Tierra
tiemp_sem_combtierr<-lm(tiemp_g~TIERRA*ESTERIL, data=Experim_UR)
tiemp_sem_combtierr<-lm(tiemp_g~VASO, data=Experim_UR)
summary(tiemp_sem_combtierr)
anova(tiemp_sem_combtierr)

TukeyHSD(aov(tiemp_sem_combtierr))

e_plot4b<-ggplot(Experim_UR, aes(x=VASO, y=tiemp_e, fill=VASO))+
  geom_boxplot()+theme_classic()+
  labs(x="Tratamiento de tierra", y="Tiempo de emer (diaas)")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))






print(grid.arrange(e_plot1a,e_plot1b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


grid.arrange(e_plot2a,e_plot2b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))

grid.arrange(e_plot3a,e_plot3b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))


grid.arrange(e_plot4a,e_plot4b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))







###Tiempo de Germinación por Combinacion de Tratamientos

#ORIGEN TIERRA * ORIGEN SEMILLA
tiemp_sem_combtot<-lm(tiemp_g~TIERRA*ORIGEN, data=Experim_UR)
summary(tiemp_sem_combtot)
anova(tiemp_sem_combtot)

TukeyHSD(aov(tiemp_sem_combtot))

e_plot5b<-ggplot(Experim_UR, aes(x=TIERRA, y=tiemp_e, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de tierra", y="Tiempo de emer (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))


#ESTERILIDAD DE  TIERRA * ORIGEN SEMILLA
tiemp_sem_combtot<-lm(tiemp_g~ESTERIL*ORIGEN, data=Experim_UR)
summary(tiemp_sem_combtot)
anova(tiemp_sem_combtot)

TukeyHSD(aov(tiemp_sem_combtot))

e_plot6b<-ggplot(Experim_UR, aes(x=ESTERIL, y=tiemp_e, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Esterilidad de Tierra", y="Tiempo de emer (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))


#TRATAMIENTOS DE TIERRA * ORIGEN SEMILLA
tiemp_sem_combtot<-lm(tiemp_g~VASO*ORIGEN, data=Experim_UR)
summary(tiemp_sem_combtot)
anova(tiemp_sem_combtot)

TukeyHSD(aov(tiemp_sem_combtot))

e_plot7b<-ggplot(Experim_UR, aes(x=VASO, y=tiemp_e, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Tratamiento de tierra", y="Tiempo de emer (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)





grid.arrange(plot5a,plot5b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))



grid.arrange(plot6a,plot6b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))



grid.arrange(plot7a,plot7b,
             widths= c(1,1), 
             layout_matrix= rbind(c(1,2)))






















###############################################################################
###############################################################################
###############################################################################



###################
#    AUTOTROFIA   #
###################

###############################################################################
###############################################################################
###############################################################################








rtaA<-cbind(Experim_UR$AUTO,Experim_UR$weights-Experim_UR$AUTO)



Experim_UR$AUTO_<-Experim_UR$AUTO/Experim_UR$EMER
rtaA_<-cbind(Experim_UR$AUTO,Experim_UR$EMER-Experim_UR$AUTO)

Experim_UR$AUTO_>1
Experim_UR$AUTO_

### Numero de germinación


#Numero de Germinación por Origen de semilla
hist(Experim_UR$tiemp_e)

#binomial agregado

#num_sem_emer <- lm(EMER01~ORIGEN, data=Experim_UR)
#summary(num_sem_emer)


#num_srsemem_o<-glm(EMER01~ORIGEN, data=Experim_UR,
#                   family=binomial, weights = weights)


tiemp_sem_oraut<-glm(tiemp_h~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR,
                     family=binomial, weights = weights)

#num_sem_emer<-glm(rtaE~ESTERIL, data=Experim_UR,
#                  family=binomial(logit))


num_sem_oraut<-glmer(rtaA_~PESO_m_s+ORIGEN*TIERRA*ESTERIL+(1|SEMILLA), 
                     data=Experim_UR,
                     family=binomial(logit))
anova(num_sem_orsem, test="Chisq")

summary(num_sem_emer)
plot(EMER~ESTERIL, data=Experim_UR)
Anova(num_sem_orsem)



prob_fit_sk1a <- ggpredict(num_sem_oraut, terms=c("TIERRA","ORIGEN","ESTERIL"))
sk1a<-plot(prob_fit_sk1a)



Experim_UR1a<-(gather(Experim_UR, "H1","H2","H3","H4","H5","H6","H7","H8","H9","H10", 
                   key="seed", value="dias"))
tiemp_sem_oraut<-glm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1a,
                     family=Gamma)


prob_fit_sk2a <- ggpredict(tiemp_sem_oraut, terms=c("TIERRA","ORIGEN","ESTERIL"))
sk2a<-plot(prob_fit_sk2a)




print(grid.arrange(sk1a,sk2a,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))






tiemp_sem_orsem0<-lm(tiemp_g~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR)
anova(tiemp_sem_orsem0)

tiemp_sem_orsem<-lm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1)


tiemp_sem_orsem<-glm(dias~PESO_m_s+ORIGEN*TIERRA*ESTERIL, data=Experim_UR1,
                     family=Gamma)

anova(tiemp_sem_orsem, test="F")

Anova(tiemp_sem_orsem)
summary(aov(tiemp_sem_orsem))
plot(tiemp_sem_orsem)

summary.aov(tiemp_sem_orsem, split = list(TIERRA=list(Rural=1, Urban=2)))


TukeyHSD(aov(tiemp_sem_orsem))




ggplot(Experim_UR1, aes(x=VASO, y=dias, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Tiempo de germinación de semillas por tratamiento", 
       x="Tratamiento de tierra", y="Tiempo de germ (d????as)",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))










por_auto<-data.frame(
  TIERRA=c("Rural","Urban"),
  meanT=c(mean(Experim_UR$AUTO01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
          mean(Experim_UR$AUTO01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  sdT=c(sd(Experim_UR$AUTO01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        sd(Experim_UR$AUTO01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  seT=c(se(Experim_UR$AUTO01[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        se(Experim_UR$AUTO01[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  ESTERIL=c("NO","YES"),
  meanE=c(mean(Experim_UR$AUTO01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
          mean(Experim_UR$AUTO01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  sdE=c(sd(Experim_UR$AUTO01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        sd(Experim_UR$AUTO01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  seE=c(se(Experim_UR$AUTO01[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        se(Experim_UR$AUTO01[Experim_UR$ESTERIL=="YES"], na.rm = T)*100))




a_plot1a<-ggplot(por_auto, 
                 aes(x=TIERRA, y=meanT, fill=TIERRA))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=TIERRA, ymin=meanT-seT, ymax=meanT+seT, 
                    width=0.1))+
  labs(x="Origen de Tierra", y="% Hoja Verdadera")+
  geom_text(aes(x=TIERRA,y=10+meanT+seT,
                label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))+
  guides(fill=F)



a_plot2a<-ggplot(por_auto, 
                 aes(x=ESTERIL, y=meanE, fill=ESTERIL))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=ESTERIL, ymin=meanE-seE, ymax=meanE+seE, 
                    width=0.1))+
  labs(x="Estado de Esterilidad", y="% Hoja Verdadera")+
  geom_text(aes(x=ESTERIL,y=+meanE+seE-10,
                label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("NO" = "No Esteril", "YES" = "Esteril"))+
  guides(fill=F)






Experim_UR %>% 
  group_by(TIERRA) %>% 
  summarize(
    mean=mean(AUTO01, na.rm = T)*100,
    sd=sd(AUTO01, na.rm = T)*100,
    se=se(AUTO01, na.rm = T)*100
  ) -> por_auto_completo


a_plot3a<-ggplot(por_auto_completo, 
                 aes(x=VASO, y=mean, fill=ORIGEN))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black", position="dodge") +
  geom_errorbar(aes(x=VASO, ymin=mean-se, ymax=mean+se), 
                width=0.1, position=position_dodge(width = 0.9))+
  labs(x="Estado de Esterilidad", y="% Hoja Verdadera")+
  scale_fill_manual(values=c("#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff",
                             "#00b835", "#f8756b","85de9e","ebc2bfff"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)


print(grid.arrange(a_plot1a, a_plot2a,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


print(grid.arrange(a_plot3a, a_plot7b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))  



#########################  #########################
######################### #########################
#########################  #########################
######################### #########################

tiemp_auto<-data.frame(
  TIERRA=c("Rural","Urban"),
  meanT=c(mean(Experim_UR$tiemp_h[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
          mean(Experim_UR$tiemp_h[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  sdT=c(sd(Experim_UR$tiemp_h[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        sd(Experim_UR$tiemp_h[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  seT=c(se(Experim_UR$tiemp_h[Experim_UR$TIERRA=="Rural"], na.rm = T)*100,
        se(Experim_UR$tiemp_h[Experim_UR$TIERRA=="Urban"], na.rm = T)*100),
  ESTERIL=c("NO","YES"),
  meanE=c(mean(Experim_UR$tiemp_h[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
          mean(Experim_UR$tiemp_h[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  sdE=c(sd(Experim_UR$tiemp_h[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        sd(Experim_UR$tiemp_h[Experim_UR$ESTERIL=="YES"], na.rm = T)*100),
  seE=c(se(Experim_UR$tiemp_h[Experim_UR$ESTERIL=="NO"], na.rm = T)*100,
        se(Experim_UR$tiemp_h[Experim_UR$ESTERIL=="YES"], na.rm = T)*100))


Experim_UR %>% 
  group_by(VASO,ORIGEN) %>% 
  summarise(
    mean=mean(tiemp_h, na.rm = T)*100,
    sd=sd(tiemp_h,na.rm = T)*100,
    se=se(tiemp_h, na.rm = T)*100
  ) -> tiemp_auto_completo




print(grid.arrange(a_plot2b,a_plot3b,
                   widths= c(1,1), 
                   layout_matrix= rbind(c(1,2))))


#########################  #########################
######################### #########################
#########################  #########################
######################### #########################

num_sem_orsem<-glm(GERM01~ORIGEN, data=Experim_UR,
                   family=binomial, weights=weights)
summary(num_sem_orsem)
anova(num_sem_orsem)


a_plot1a<-ggplot(Experim_UR, aes(x=ORIGEN, y=AUTO, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de semilla", 
       x="Ambiente de origen de semilla", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Numero de Germinación por Origen de Tierra
num_sem_ortier<-lm(GERM~TIERRA, data=Experim_UR)
num_sem_orsem<-glm(GERM01~TIERRA, data=Experim_UR,
                   family=binomial, weights=weights)
summary(num_sem_ortier)
anova(num_sem_ortier)
Anova(num_sem_ortier)

mean(Experim_UR$GERM[Experim_UR$TIERRA == "Urban"])/
  mean(Experim_UR$GERM[Experim_UR$TIERRA == "Rural"])

a_plot2a<-ggplot(Experim_UR, aes(x=TIERRA, y=AUTO, fill=TIERRA))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de tierra", 
       x="Ambiente de origen de tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Numero de Germinación por Esterilidad de Tierra
num_sem_estril<-lm(GERM~ESTERIL, data=Experim_UR)
summary(num_sem_estril)
anova(num_sem_estril)

mean(Experim_UR$GERM[Experim_UR$ESTERIL == "NO"])/
  mean(Experim_UR$GERM[Experim_UR$ESTERIL == "YES"])

mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL == "YES"], na.rm = T)/
  mean(Experim_UR$tiemp_g[Experim_UR$ESTERIL == "NO"], na.rm = T)

a_plot3a<-ggplot(Experim_UR, aes(x=ESTERIL, y=GERM, fill=ESTERIL))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por esterilidad de tierra", 
       x="Esterilidad de Tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))





#Numero de Germinación por Combinacion de Tierra
num_sem_combtierr<-lm(GERM~TIERRA*ESTERIL, data=Experim_UR)
num_sem_combtierr<-lm(GERM~VASO, data=Experim_UR)
summary(num_sem_combtierr)
anova(num_sem_combtierr)

TukeyHSD(aov(num_sem_combtierr))

plot4a<-ggplot(Experim_UR, aes(x=VASO, y=GERM, fill=VASO))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por tratamiento de tierra", 
       x="Tratamiento de tierra", y="Numero de semillas germinadas")+
  guides(fill=F)+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))


mean(Experim_UR$GERM[Experim_UR$VASO == "UX"], na.rm = T)/
  mean(Experim_UR$GERM[Experim_UR$VASO == "RX"], na.rm = T)

mean(Experim_UR$tiemp_g[Experim_UR$VASO == "RX"], na.rm = T)/
  mean(Experim_UR$tiemp_g[Experim_UR$VASO == "UX"], na.rm = T)






###Numero de Germinación por Combinacion de Tratamientos

#ORIGEN TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~TIERRA*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot5a<-ggplot(Experim_UR, aes(x=TIERRA, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por origen de ambiente", 
       x="Ambiente de origen de tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))


#ESTERILIDAD DE  TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~ESTERIL*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot6a<-ggplot(Experim_UR, aes(x=ESTERIL, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por esterilidad de tierra", 
       x="Esterilidad de Tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))



#TRATAMIENTOS DE TIERRA * ORIGEN SEMILLA
num_sem_combtot<-lm(GERM~VASO*ORIGEN, data=Experim_UR)
summary(num_sem_combtot)
anova(num_sem_combtot)

TukeyHSD(aov(num_sem_combtot))

plot7a<-ggplot(Experim_UR, aes(x=VASO, y=GERM, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(title="Numero de semillas germinadas por tratamiento", 
       x="Tratamiento de tierra", y="Numero de semillas germinadas",
       fill="Origen de\nSemillas")+
  #guides(fill=F)
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))










###########################
###Tiempo de germinación###
###########################


#Tiempo de Germinación por Origen de semilla
tiemp_sem_orsem<-lm(tiemp_g~ORIGEN, data=Experim_UR)
summary(tiemp_sem_orsem)
anova(tiemp_sem_orsem)

a_plot1b<-ggplot(Experim_UR, aes(x=ORIGEN, y=tiemp_h, fill=ORIGEN))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de semilla", y="Tiempo de hojas (dias)")+
  guides(fill=F)+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Tiempo de Germinación por Origen de Tierra
tiemp_sem_ortier<-lm(tiemp_g~TIERRA, data=Experim_UR)
summary(tiemp_sem_ortier)
anova(tiemp_sem_ortier)

a_plot2b<-ggplot(Experim_UR, aes(x=TIERRA, y=tiemp_h, fill=TIERRA))+
  geom_boxplot()+theme_classic()+
  labs(x="Ambiente de origen de tierra", y="Tiempo de hojas (dias)")+
  guides(fill=F)+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))



#Tiempo de Germinación por Esterilidad de Tierra
tiemp_sem_estril<-lm(tiemp_g~ESTERIL, data=Experim_UR)
summary(tiemp_sem_estril)
anova(tiemp_sem_estril)

a_plot3b<-ggplot(Experim_UR, aes(x=ESTERIL, y=tiemp_h, fill=ESTERIL))+
  geom_boxplot()+theme_classic()+
  labs(x="Esterilidad de Tierra", y="Tiempo de hojas (dias)")+
  guides(fill=F)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))









###############################################################################
###############################################################################
###############################################################################



###########################
###     Colonization    ###
###########################


###############################################################################
###############################################################################
###############################################################################


Experim_UR %>% 
  group_by(VASO) %>% 
  summarize(
    mean=mean(prop01, na.rm = T)*100,
    sd=sd(prop01,na.rm = T)*100,
    se=se(prop01, na.rm = T)*100
  ) -> por_colV

Experim_UR %>% 
  group_by(TIERRA, ORIGEN) %>% 
  summarize(
    mean=mean(prop01, na.rm = T)*100,
    sd=sd(prop01,na.rm = T)*100,
    se=se(prop01, na.rm = T)*100
  ) -> por_colTO

Experim_UR %>% 
  group_by(ESTERIL, ORIGEN) %>% 
  summarize(
    mean=mean(prop01, na.rm = T)*100,
    sd=sd(prop01,na.rm = T)*100,
    se=se(prop01, na.rm = T)*100
  ) -> por_colEO

Experim_UR %>% 
  group_by(VASO,ORIGEN) %>% 
  summarize(
    mean=mean(prop01, na.rm = T)*100,
    sd=sd(prop01,na.rm = T)*100,
    se=se(prop01, na.rm = T)*100
  ) -> por_colVO








#Colonización por Origen de semilla
colexp_orsem<-lm(prop01~ORIGEN, data = Experim_UR)
plot(colexp_orsem,2)

colexp_orsem<-glm(prop01~ORIGEN, data = Experim_UR,
               family=quasibinomial, weights = Sector)
summary(colexp_orsem)
Anova(colexp_orsem)



rtaCOL <- c(Experim_UR$Col,Experim_UR$Sector - Experim_UR$Col)


col_plot1<-ggplot(por_col, 
  aes(x=ORIGEN, y=meanS, fill=ORIGEN))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=ORIGEN, ymin=meanS-seS, ymax=meanS+seS, 
                    width=0.1))+
  labs(x="Origen de Semilla", y="% Colonization")+
  #geom_text(aes(x=TIERRA,y=10+meanS+seS,
  #              label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))+
  guides(fill=F)
  





#Colonización por Origen de tierra
colexp_ortier<-lm(prop01~TIERRA, data = Experim_UR)
plot(colexp_ortier,2)

colexp_ortier<-glm(prop01~TIERRA, data = Experim_UR, 
               family=binomial, weights = Sector)
summary(colexp_ortier)
Anova(colexp_ortier)

col_plot2<-ggplot(por_col, 
  aes(x=TIERRA, y=meanT, fill=TIERRA))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=TIERRA, ymin=meanT-seT, ymax=meanT+seT, 
                    width=0.1))+
  labs(x="Origen de Tierra", y="% Colonization")+
  #geom_text(aes(x=TIERRA,y=10+meanS+seS,
  #              label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  scale_x_discrete(labels=c("Rural" = "Rural", "Urban" = "Urbano"))+
  guides(fill=F)





#Colonización por esterilidad de tierra
colexp_esteril<-lm(prop01~ESTERIL, data = Experim_UR)
plot(colexp_esteril,2)

colexp_esteril<-glm(prop01~ESTERIL, data = Experim_UR, 
              family=binomial, weights = Sector)
summary(colexp_esteril)
Anova(colexp_esteril)



col_plot3<-ggplot(por_col, 
  aes(x=ESTERIL, y=meanE, fill=TIERRA))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(x=ESTERIL, ymin=meanE-seE, ymax=meanE+seE, 
                    width=0.1))+
  labs(x="Esterilidad de Tierra", y="% Colonization")+
  #geom_text(aes(x=ESTERIL,y=10+meanE+seE,
  #              label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("NO" = "No Esteril", "YES" = "Esteril"))+
  guides(fill=F)




print(grid.arrange(col_plot1,col_plot2, col_plot3,
                   widths= c(1,1,1), 
                   layout_matrix= rbind(c(1,2,3))))

















#Colonización por tratamiento de tierra
colexp_combtier<-lm(prop01~VASO, data = Experim_UR,
                    subset = ESTERIL=="NO")
plot(colexp_combtier,2)

colexp_combtier<-glm(prop01~VASO, data = Experim_UR, 
                   family=quasibinomial, 
                   weights = Sector #,subset = Experim_UR$VASO!="RR"
                   ) #RR > RX 
summary(colexp_combtier) #RR -- RX***, UU   , UX***
                         #RX --        UU***, UX·
                         #UU --               UX***
Anova(colexp_combtier)


TukeyHSD(aov(colexp_combtier))

col_plot4<-ggplot(por_colV, aes(x=VASO, y=mean, fill=VASO))+
  theme_classic()+ 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(x=VASO, ymin=mean-se, ymax=mean+se, 
                    width=0.1))+
  labs(x="Tratamiento de tierra", y="% Colonization")+
  geom_text(aes(x=VASO,y=3+mean+se,
                label=c("a","b","a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)





mean(Experim_UR$prop01[Experim_UR$VASO == "UX"], na.rm = T)/
  mean(Experim_UR$prop01[Experim_UR$VASO == "RX"], na.rm = T)





#Colonización por origenes de tierra y semilla
colexp_origenes<-lm(prop01~TIERRA*ORIGEN, data = Experim_UR,
                    subset=ESTERIL=="NO")
plot(colexp_origenes,2)

colexp_origenes<-glm(prop01~TIERRA*ORIGEN, data = Experim_UR, 
                    family=quasibinomial,  
                    weights = Sector #,subset = Experim_UR$VASO=="RR"
                    ) #RR > RX 
summary(colexp_origenes)
anova(colexp_origenes)
Anova(colexp_origenes, type="3", test="F")

TukeyHSD(aov(colexp_origenes))

col_plot5<-ggplot(por_colTO, aes(x=TIERRA, y=mean, fill=interaction(TIERRA,ORIGEN)))+
  theme_classic()+ 
  geom_bar(position="dodge", stat="identity", color="black")+
  geom_errorbar(aes(x=TIERRA, ymin=mean-se, ymax=mean+se), 
                width=0.1, position=position_dodge(.9))+
  labs(x="Origen de tierra", y="% Colonization")+
  #geom_text(aes(x=TIERRA,y=3+mean+se,
  #              label=c("a","a","a","a")),
  #          vjust=0, position = position_dodge(.9))+
  
  scale_fill_manual(values=c("#00ba38","#00ba38", "#f8766d","#f8766d"))+
  scale_x_discrete(labels=c("Rural", "Urban"))+
  guides(fill=F)
  

  




#ESTERILIDAD DE  ESTERILIDAD * ORIGEN SEMILLA
colexp_semesteril<-lm(prop01~ESTERIL*ORIGEN, data = Experim_UR)
colexp_semesteril<-glm(prop01~ESTERIL*ORIGEN, data = Experim_UR, 
                     family=quasibinomial,  
                     weights = Sector #,subset = Experim_UR$VASO=="RR"
                     ) #RR > RX 
summary(colexp_semesteril)
Anova(colexp_semesteril, type="3", test="F")

TukeyHSD(aov(colexp_semesteril))

col_plot6<-ggplot(por_colEO, aes(x=ESTERIL, y=mean, 
                                 fill=interaction(ESTERIL,ORIGEN)))+
  theme_classic()+ 
  geom_bar(position="dodge", stat="identity", color="black")+
  geom_errorbar(aes(x=ESTERIL, ymin=mean-se, ymax=mean+se), 
                width=0.1, position=position_dodge(.9))+
  labs(x="Origen de tierra", y="% Colonization")+
  geom_text(aes(x=ESTERIL,y=3+mean+se,
                label=c("a","a","b","b")),
            vjust=0, position = position_dodge(.9))+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))+
  guides(fill=F)

  

  
  
print(grid.arrange(col_plot4,col_plot5, col_plot6,
                    widths= c(1,1,1), 
                    layout_matrix= rbind(c(1,2,3))))
  















#TRATAMIENTOS DE TIERRA * ORIGEN SEMILLA
colexp_tot<-glm(prop01~TIERRA*ESTERIL*ORIGEN, data = Experim_UR, 
                     family=quasibinomial,  
                weights = Sector #,subset = Experim_UR$VASO=="RR"
                ) #RR > RX 
summary(colexp_tot)
Anova(colexp_tot, type="3", test="F")

TukeyHSD(aov(colexp_tot))

col_plot7<-ggplot(por_colVO, aes(x=VASO, y=mean, 
                                 fill=interaction(VASO,ORIGEN)))+
  theme_classic()+ 
  geom_bar(position="dodge", stat="identity", color="black")+
  geom_errorbar(aes(x=VASO, ymin=mean-se, ymax=mean+se), 
                width=0.1, position=position_dodge(.9))+
  labs(x="Tratamiento de tierra", y="% Colonization")+
  geom_text(aes(x=VASO,y=3+mean+se,
                label=c("a","a","c","bc",
                        "a","ab","c","c")),
            vjust=0, position = position_dodge(.9))+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#00ba38", "#85df9f",
                             "#f8766d","#ecc3c0","#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)










###############################################################################
###############################################################################
###############################################################################




AUTOCOL<-lm(tiemp_h~prop01*ORIGEN, data=Experim_UR)
anova(AUTOCOL)


###############################################################################
###############################################################################
###############################################################################

lmod<-lm(Fv_Fm~ORIGEN, data = Exp.Fv_Fm)
plot(lmod,2)

summary(lmod)
Anova(lmod)

Exp.Fv_Fm$Fv_Fm[119]

Experim_UR_FQ %>% 
  group_by(SEMILLA, VASO, ORIGEN, TIERRA, ESTERIL, MED) %>% 
  summarize(Fv_Fm=mean(Fv_Fm, na.rm=T)) -> Exp.Fv_Fm

Exp.Fv_Fm_out<-Exp.Fv_Fm[!Exp.Fv_Fm$Fv_Fm<0.6,]
Exp.Fv_Fm_out<-Exp.Fv_Fm[-c(103,111,118,200),]
Exp.Fv_Fm[118,]


Exp.Fv_Fm$SEMILLA<-as.factor(Exp.Fv_Fm$SEMILLA)
Exp.Fv_Fm$VASO<- as.factor(Exp.Fv_Fm$VASO)
Exp.Fv_Fm$ORIGEN<- as.factor(Exp.Fv_Fm$ORIGEN)
Exp.Fv_Fm$TIERRA <-as.factor(Exp.Fv_Fm$TIERRA)
Exp.Fv_Fm$ESTERIL <-as.factor(Exp.Fv_Fm$ESTERIL)
Exp.Fv_Fm$MED <-as.factor(Exp.Fv_Fm$MED)

levels(Exp.Fv_Fm$ORIGEN)

hist(Exp.Fv_Fm_out$Fv_Fm)
hist(asin(sqrt(Exp.Fv_Fm_out$Fv_Fm)))






#Colonización por Origen de semilla
FVFM_orsem<-lm(Fv_Fm~ORIGEN, data = Exp.Fv_Fm_out)
plot(FVFM_orsem,2)

summary(FVFM_orsem)
Anova(FVFM_orsem)

FVFM_plot1<-ggplot(Exp.Fv_Fm_out, 
                  aes(x=ORIGEN, y=Fv_Fm, fill=ORIGEN))+
  theme_classic()+ geom_boxplot()+
  labs(x="Origen de Semilla", y="Fm/Fm")+
  #geom_text(aes(x=TIERRA,y=10+meanS+seS,
  #              label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  guides(fill=F)


FVFM_orsem<-lm(asin(sqrt(Exp.Fv_Fm_out$Fv_Fm))~ORIGEN, 
               data = Exp.Fv_Fm_out)
plot(FVFM_orsem,2)

summary(FVFM_orsem)
Anova(FVFM_orsem)


###ggpredict

FVFMfit1 <- ggpredict(FVFM_orsem, terms=c("ORIGEN"), type="re")
FVFMgg1<-plot(FVFMfit1)
View(FVFMfit1)


#Colonización por Origen de tierra
FVFM_ortier<-lm(Fv_Fm~TIERRA, data = Exp.Fv_Fm_out)
plot(FVFM_ortier,2)

summary(FVFM_ortier)
Anova(FVFM_ortier)

FVFM_plot2<-ggplot(Exp.Fv_Fm_out, 
                  aes(x=TIERRA, y=Fv_Fm, fill=TIERRA))+
  theme_classic()+ geom_boxplot()+
  labs(x="Origen de Tierra", y="Fm/Fm")+
  #geom_text(aes(x=TIERRA,y=10+meanS+seS,
  #              label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#00ba38","#f8766d"))+
  guides(fill=F)


###ggpredict

FVFMfit2 <- ggpredict(FVFM_ortier, terms=c("TIERRA"))
FVFMgg2<-plot(FVFMfit2)


#Colonización por esterilidad de tierra
FVFM_esteril<-lm(Fv_Fm~ESTERIL, data = Exp.Fv_Fm_out)
plot(FVFM_esteril,2)

summary(FVFM_esteril)
Anova(FVFM_esteril)

FVFM_plot3<-ggplot(Exp.Fv_Fm_out, 
                  aes(x=ESTERIL, y=Fv_Fm, fill=ESTERIL))+
  theme_classic()+ geom_boxplot()+
  labs(x="Esterilidad de Tierra", y="Fm/Fm")+
  #geom_text(aes(x=ESTERIL,y=10+meanE+seE,
  #              label=c("b","a")),vjust=0)+
  scale_fill_manual(values=c("#999999", "#ffffff"))+
  scale_x_discrete(labels=c("NO" = "No Esteril", "YES" = "Esteril"))+
  guides(fill=F)



###ggpredict

FVFMfit3 <- ggpredict(FVFM_esteril, terms=c("ESTERIL"))
FVFMgg3<-plot(FVFMfit3)




print(grid.arrange(FVFM_plot1, FVFM_plot2, FVFM_plot3,
                   widths= c(1,1,1), 
                   layout_matrix= rbind(c(1,2,3))))



grid.arrange(FVFMfit1, FVFMfit2, FVFMfit3,
                   widths= c(1,1,1), 
                   layout_matrix= rbind(c(1,2,3)))













#Colonización por tratamiento de tierra
FVFM_combtier<-lm(Fv_Fm~VASO, data = Exp.Fv_Fm_out)
plot(FVFM_combtier,2)

#,subset = Exp.Fv_Fm_out$VASO!="RR"
#RR > RX 
summary(FVFM_combtier) #RR -- RX***, UU   , UX***
#RX --        UU***, UX·
#UU --               UX***
Anova(FVFM_combtier)


TukeyHSD(aov(FVFM_combtier))

FVFM_combtier_<-lm(Fv_Fm~TIERRA*ESTERIL, data = Exp.Fv_Fm_out)



FVFM_plot4<-ggplot(Exp.Fv_Fm_out, aes(x=VASO, y=Fv_Fm, 
                                      fill=VASO))+
  theme_classic()+ geom_boxplot()+
  labs(x="Tratamiento de tierra", y="Fv/Fm")+
  #geom_text(aes(x=VASO,y=3+mean+se,
  #              label=c("a","b","a","b")),vjust=0)+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)


###ggpredict

FVFMfit4 <- ggpredict(FVFM_combtier_, terms=c("TIERRA","ESTERIL"))
FVFMgg4<-plot(FVFMfit4, connect.lines=T)


mean(Exp.Fv_Fm_out$Fv_Fm[Exp.Fv_Fm_out$VASO == "UX"], na.rm = T)/
  mean(Exp.Fv_Fm_out$Fv_Fm[Exp.Fv_Fm_out$VASO == "RX"], na.rm = T)







#Colonización por origenes de tierra y semilla
FVFM_origenes<-lm(Fv_Fm~TIERRA*ORIGEN, data = Exp.Fv_Fm_out)
plot(FVFM_origenes,2)

#,subset = Exp.Fv_Fm_out$VASO=="RR"
#RR > RX 
summary(FVFM_origenes)
anova(FVFM_origenes)
Anova(FVFM_origenes, type="3", test="F")

TukeyHSD(aov(FVFM_origenes))

FVFM_plot5<-ggplot(Exp.Fv_Fm_out, aes(x=TIERRA, y=Fv_Fm, 
                      fill=interaction(TIERRA,ORIGEN)))+
  theme_classic()+ geom_boxplot()+
  labs(x="Origen de tierra", y="Fv/Fm")+
  #geom_text(aes(x=TIERRA,y=3+mean+se,
  #              label=c("a","a","a","a")),
  #          vjust=0, position = position_dodge(.9))+
  scale_fill_manual(values=c("#00ba38","#00ba38", "#f8766d","#f8766d"))+
  scale_x_discrete(labels=c("Rural", "Urban"))+
  guides(fill=F)


###ggpredict

FVFMfit5 <- ggpredict(FVFM_origenes, terms=c("TIERRA","ORIGEN"))
FVFMgg5<-plot(FVFMfit5, connect.lines=T)







#ESTERILIDAD DE  ESTERILIDAD * ORIGEN SEMILLA
FVFM_semesteril<-lm(Fv_Fm~ESTERIL*ORIGEN, data = Exp.Fv_Fm_out)
#,subset = Exp.Fv_Fm_out$VASO=="RR"
#RR > RX 
summary(FVFM_semesteril)
Anova(FVFM_semesteril, type="3", test="F")

TukeyHSD(aov(FVFM_semesteril))

FVFM_plot6<-ggplot(Exp.Fv_Fm_out, aes(x=ESTERIL, y=Fv_Fm, 
                                 fill=interaction(ESTERIL,ORIGEN)))+
  theme_classic()+ geom_boxplot()+
  labs(x="Origen de tierra", y="Fv/Fm")+
  #geom_text(aes(x=ESTERIL,y=3+mean+se,
  #              label=c("a","a","b","b")),
  #          vjust=0, position = position_dodge(.9))+
  scale_fill_manual(labels=c("Rural" = "Rural", "Urban" = "Urbano"),
                    values=c("#00ba38", "#85df9f", "#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("No Esteril", "Esteril"))+
  guides(fill=F)

###ggpredict

FVFMfit6 <- ggpredict(FVFM_semesteril, terms=c("ESTERIL","ORIGEN"))
FVFMgg6<-plot(FVFMfit6, connect.lines=T)



print(grid.arrange(FVFM_plot4, FVFM_plot5, FVFM_plot6,
                   widths= c(1,1,1), 
                   layout_matrix= rbind(c(1,2,3))))




print(grid.arrange(FVFMfit4, FVFMfit5, FVFMfit6,
                   widths= c(1,1,1), 
                   layout_matrix= rbind(c(1,2,3))))











#TRATAMIENTOS DE TIERRA * ORIGEN SEMILLA
FVFM_tot<-lm(Fv_Fm~TIERRA*ESTERIL*ORIGEN, data = Exp.Fv_Fm_out) 
#,subset = Exp.Fv_Fm_out$VASO=="RR"
#RR > RX 
summary(FVFM_tot)
Anova(FVFM_tot, type="3", test="F")

TukeyHSD(aov(FVFM_tot))

FVFM_plot7<-ggplot(Exp.Fv_Fm_out, aes(x=VASO, y=Fv_Fm, 
                                 fill=interaction(VASO,ORIGEN)))+
  theme_classic()+ geom_boxplot()+
  labs(x="Tratamiento de tierra", y="Fv/Fm")+
  #geom_text(aes(x=VASO,y=3+mean+se,
  #              label=c("a","a","c","bc",
  #                      "a","ab","c","c")),
  #          vjust=0, position = position_dodge(.9))+
  scale_fill_manual(values=c("#00ba38", "#85df9f", "#00ba38", "#85df9f",
                             "#f8766d","#ecc3c0","#f8766d","#ecc3c0"))+
  scale_x_discrete(labels=c("Rural\nNo Esteril", "Rural\nEsteril", 
                            "Urbano\nNo Esteril", "Urbano\nEsteril" ))+
  guides(fill=F)




FVFM_plot7



###ggpredict

FVFMfit7 <- ggpredict(FVFM_tot, terms=c("TIERRA","ORIGEN","ESTERIL"),type="re")
structure(FVFMfit1)
FVFMgg7<-plot(FVFMfit7, connect.lines=T)


FVFMgg1
FVFMgg2
FVFMgg3
FVFMgg4
FVFMgg5
FVFMgg6
FVFMgg7

as.data.frame(FVFMgg1)

ggPredict(FVFM_tot)+
  geom_ribbon(aes(ymin = predicted - std.error, 
                  ymax = predicted + std.error), alpha = 0.5)

  
  
  
  
  
  

#########################################
#######                           #######
#######        FOTOQUIMICA        #######
#######                           #######
#########################################

source("https://dornsife.usc.edu/assets/sites/239/docs/Rallfun-v38.txt")


#############
#### ETR ####
#############


Experim_UR_ETR<-Experim_UR_FQ[,c(1:6,65,70,71)]
Experim_UR_ETR<-as.data.frame(pivot_wider(Experim_UR_ETR, 
                          names_from = MED, values_from = ETR))
Experim_UR_ETR<-Experim_UR_ETR[,-8]
names(Experim_UR_ETR)[8:15]<-c("ETR2","ETR3","ETR4","ETR5",
                               "ETR6","ETR7","ETR8","ETR9")


m_ETR<-as.matrix(Experim_UR_ETR[8:15])
pairs(m_ETR)
mshapiro.test(m_ETR[,1:4])


ETR_mlm<-lm(m_ETR~TIERRA, Experim_UR_ETR )
summary(ETR_mlm)
anova(ETR_mlm)

ETR_manova<-manova(m_ETR~ESTERIL, Experim_UR_ETR)
summary(ETR_manova, intercept=T)
summary.aov(ETR_manova)


#Plot

Experim_UR_ETR.T<-Experim_UR_ETR[,-c(1,4,5,7)]
ETR_melt<-melt(Experim_UR_ETR.T, id=c("TIERRA","ESTERIL","ORIGEN"), 
     measured=c("ETR2","ETR3","ETR4","ETR5",
                "ETR6","ETR7","ETR8","ETR9"))
names(ETR_melt)[c(4,5)]<-c("Measure","Value")


ETR_melt %>% 
  group_by(ESTERIL, Measure) %>% 
  summarize(mean=mean(Value, na.rm=T),
            sd=sd(Value, na.rm=T),
            se=se(Value, na.rm=T)) -> Sum.ETR_melt


ggplot(Sum.ETR_melt, aes(Measure, mean, group=ESTERIL, color = ESTERIL)) +
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR")+
  scale_x_discrete(labels= c("125","190","285","420",
                 "625","845","1150","1500"))




##############
#### PSII ####
##############


#names(Experim_UR_FQ1)
#Experim_UR_FQ1[,c(1:6,63:71)]
Experim_UR_PSII<-Experim_UR_FQ1[,c(1:6,64,70,71)]
Experim_UR_PSII<-as.data.frame(pivot_wider(Experim_UR_PSII, 
                                          names_from = MED, values_from = PSII))
names(Experim_UR_PSII)[8:16]<-c("PSII1","PSII2","PSII3","PSII4","PSII5",
                               "PSII6","PSII7","PSII8","PSII9")


m_PSII<-as.matrix(Experim_UR_PSII[8:16])
pairs(m_PSII)

PSII_mlm<-lm(PSII1~TIERRA, Experim_UR_PSII )
anova(PSII_mlm)

PSII_manova<-manova(m_PSII~TIERRA, Experim_UR_PSII)
summary(PSII_manova, intercept=T)
summary.aov(PSII_manova)

PSII_manova<-manova(m_PSII~ORIGEN, Experim_UR_PSII)
summary(PSII_manova, intercept=T)
summary.aov(PSII_manova)

PSII_manova<-manova(m_PSII~ESTERIL, Experim_UR_PSII)
summary(PSII_manova, intercept=T)
summary.aov(PSII_manova)

#Plot
Experim_UR_PSII.T<-Experim_UR_PSII[,-c(1,4,5,7)]
PSII_melt<-melt(Experim_UR_PSII.T, id=c("TIERRA","ESTERIL","ORIGEN"), 
              measured=c("PSII2","PSII3","PSII4","PSII5",
                         "PSII6","PSII7","PSII8","PSII9"))
names(PSII_melt)[c(4,5)]<-c("Measure","Value")


PSII_melt %>% 
  group_by(ESTERIL, Measure) %>% 
  summarize(mean=mean(Value, na.rm=T),
            sd=sd(Value, na.rm=T),
            se=se(Value, na.rm=T)) -> Sum.PSII_melt


ggplot(Sum.PSII_melt, aes(Measure, mean, group=ESTERIL, color = ESTERIL)) +
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="PSII")+
  scale_x_discrete(labels= c("125","190","285","420",
                             "625","845","1150","1500"))


##############
####  qP  ####
##############


#names(Experim_UR_FQ1)
#Experim_UR_FQ1[,c(1:6,63:71)]
Experim_UR_qP<-Experim_UR_FQ1[,c(1:6,66,70,71)]
Experim_UR_qP<-as.data.frame(pivot_wider(Experim_UR_qP, 
                                           names_from = MED, values_from = qP))
Experim_UR_qP<-Experim_UR_qP[,-8]
names(Experim_UR_qP)[8:15]<-c("qP2","qP3","qP4","qP5",
                                "qP6","qP7","qP8","qP9")


m_qP<-as.matrix(Experim_UR_qP[8:15])
pairs(m_qP)

qP_mlm<-lm(qP2~TIERRA, Experim_UR_qP )
anova(qP_mlm)

qP_manova<-manova(m_qP~ESTERIL, Experim_UR_qP)
summary(qP_manova, intercept=T)
summary.aov(qP_manova)


#Plot
Experim_UR_qP.T<-Experim_UR_qP[,-c(1,4,5,7)]
qP_melt<-melt(Experim_UR_qP.T, id=c("TIERRA","ESTERIL","ORIGEN"), 
              measured=c("qP2","qP3","qP4","qP5",
                         "qP6","qP7","qP8","qP9"))
names(qP_melt)[c(4,5)]<-c("Measure","Value")


qP_melt %>% 
  group_by(ESTERIL, Measure) %>% 
  summarize(mean=mean(Value, na.rm=T),
            sd=sd(Value, na.rm=T),
            se=se(Value, na.rm=T)) -> Sum.qP_melt


ggplot(Sum.qP_melt, aes(Measure, mean, group=ESTERIL, color = ESTERIL)) +
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="qP")+
  scale_x_discrete(labels= c("125","190","285","420",
                             "625","845","1150","1500"))


##############
####  qN  ####
##############


#names(Experim_UR_FQ1)
#Experim_UR_FQ1[,c(1:6,63:71)]
Experim_UR_qN<-Experim_UR_FQ1[,c(1:6,67,70,71)]
Experim_UR_qN<-as.data.frame(pivot_wider(Experim_UR_qN, 
                                         names_from = MED, values_from = qN))
Experim_UR_qN<-Experim_UR_qN[,-8]
names(Experim_UR_qN)[8:15]<-c("qN1","qN2","qN3","qN4","qN5",
                              "qN6","qN7","qN8","qN9")


m_qN<-as.matrix(Experim_UR_qN[8:15])
pairs(m_qN)

qN_mlm<-lm(qN1~TIERRA, Experim_UR_qN )
anova(qN_mlm)

qN_manova<-manova(m_qN~ESTERIL, Experim_UR_qN)
summary(qN_manova, intercept=T)
summary.aov(qN_manova)


#Plot
Experim_UR_qN.T<-Experim_UR_qN[,-c(1,4,5,7)]
qN_melt<-melt(Experim_UR_qN.T, id=c("TIERRA","ESTERIL","ORIGEN"), 
               measured=c("qN2","qN3","qN4","qN5",
                          "qN6","qN7","qN8","qN9"))
names(qN_melt)[c(4,5)]<-c("Measure","Value")


qN_melt %>% 
  group_by(ESTERIL, Measure) %>% 
  summarize(mean=mean(Value, na.rm=T),
            sd=sd(Value, na.rm=T),
            se=se(Value, na.rm=T)) -> Sum.qN_melt


ggplot(Sum.qN_melt, aes(Measure, mean, group=ESTERIL, color = ESTERIL)) +
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="qN")+
  scale_x_discrete(labels= c("125","190","285","420",
                             "625","845","1150","1500"))


###############
####  NPQ  ####
###############


#names(Experim_UR_FQ1)
#Experim_UR_FQ1[,c(1:6,63:71)]
Experim_UR_NPQ<-Experim_UR_FQ1[,c(1:6,69,70,71)]
Experim_UR_NPQ<-as.data.frame(pivot_wider(Experim_UR_NPQ, 
                                         names_from = MED, values_from = NPQ))
Experim_UR_NPQ<-Experim_UR_NPQ[,-8]
names(Experim_UR_NPQ)[8:15]<-c("NPQ2","NPQ3","NPQ4","NPQ5",
                              "NPQ6","NPQ7","NPQ8","NPQ9")


m_NPQ<-as.matrix(Experim_UR_NPQ[8:15])
pairs(m_NPQ)

NPQ_mlm<-lm(m_NPQ~ORIGEN, Experim_UR_NPQ )
anova(NPQ_mlm)
summary(NPQ_mlm)

NPQ_manova<-manova(m_NPQ~TIERRA*ESTERIL*ORIGEN, Experim_UR_NPQ)
summary(NPQ_manova)
summary.aov(NPQ_manova)


#Plot
Experim_UR_NPQ.T<-Experim_UR_NPQ[,-c(1,4,5,7)]
NPQ_melt<-melt(Experim_UR_NPQ.T, id=c("TIERRA","ESTERIL","ORIGEN"), 
               measured=c("NPQ2","NPQ3","NPQ4","NPQ5",
                          "NPQ6","NPQ7","NPQ8","NPQ9"))
names(NPQ_melt)[c(4,5)]<-c("Measure","Value")


NPQ_melt %>% 
  group_by(TIERRA,ESTERIL,ORIGEN, Measure) %>% 
  summarize(mean=mean(Value, na.rm=T),
            sd=sd(Value, na.rm=T),
            se=se(Value, na.rm=T)) -> Sum.NPQ_melt



ggplot(Sum.NPQ_melt, aes(Measure, mean, 
                         group = ESTERIL, 
                         color = ESTERIL)) +
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR")+
  scale_x_discrete(labels= c("125","190","285","420",
                             "625","845","1150","1500"))













ggplot(Sum.NPQ_melt, aes(Measure, mean, 
                         group = interaction(TIERRA,ESTERIL, ORIGEN), 
                         color = interaction(TIERRA,ESTERIL, ORIGEN),
                         shape = interaction(TIERRA,ORIGEN))) +
  geom_line() + geom_point(size=2) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0))+
  theme_classic()+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="Electron Transfer Rate (ETR)")+
  scale_x_discrete(labels= c("125","190","285","420",
                             "625","845","1150","1500"))+
  scale_color_manual(leyend="Tratamiento de Tierra",
    labels=c("Rural NO Esteril", "Urbano NO Esteril", "Rural Esteril",
             "Urban Esteril", "Rural NO Esteril", "Urbano NO Esteril", 
             "Rural Esteril", "Urban Esteril"),
    values=c("#00ba38", "#f8766d", "#85df9f", "#ecc3c0",
                           "#00ba38","#f8766d","#85df9f","#ecc3c0"))+
  scale_shape_manual(labels=c("Rural Sympatric", "Rural Alopatric",
                               "Urban Alopatric", "Urban Sympatric"),
                     values = c(19,1,2,17)) 

#00ba38 deep green
#85df9f pale green
#f8766d deep red
#ecc3c0 pale red
