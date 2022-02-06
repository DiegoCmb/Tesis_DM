library (tidyverse)
library (lme4)
library (lmerTest)
library (modelr)
library (car)
library (broom.mixed)
library (emmeans) # para calcular medias a partir de estimadores de modelos 
library (ggeffects)
library(effsize) #for effect size



db<-read_csv("Database/Exp_Estable_UR_Merida_2020(FQ).csv")[,-1]
fq<-db%>%mutate_if(is.character, as.factor)%>%
  select (VASO, SEMILLA, COUNT, CODIGO, TIERRA, ESTERIL, ORIGEN, Col , prop01 ,PAR:MED)%>%
  mutate(key = as.factor(paste(as.character(SEMILLA), as.character(COUNT), sep = ".")))%>%
  mutate (COUNT = as.factor(COUNT), MED = as.factor(MED))
#%>%filter(COUNT != "10") # no se porque me sale una obs 10 
# %>% filter (COUNT == c("1","2", "3")) # por si se quiere reducir la complejidad de la base para 
                                        # entender con dos ejemplos como funcionan las cosas y detectar errores

names(fq)[c(3,17)]<-c("MED","COUNT")

############# Promediando a nivel de SEMILLA


fq_mean <-  # promedia el valor de la lectura i entre los inds de cada vaso
  fq %>% group_by(SEMILLA, COUNT) %>% summarise(
    m_Col = mean (Col, na.rm = TRUE),   # no hay bronca, el prom de algo no varible es igual al dato en si. 
    m_prop01 = mean (prop01, na.rm = TRUE),
    m_par = mean(PAR, na.rm = TRUE),
    m_psii = mean(PSII, na.rm = TRUE),
    m_etr = mean(ETR, na.rm = TRUE),
    m_qP = mean(qP, na.rm = TRUE),
    m_qN = mean(qN, na.rm = TRUE),
    m_fv_fm = mean (Fv_Fm, na.rm = TRUE), # no hay bronca, el prom de algo no varible es igual al dato en si.
    m_npq = mean(NPQ, na.rm = TRUE)
  )

fq_fact <-
  fq %>% distinct_at(vars(SEMILLA), .keep_all = TRUE) %>% 
     select (VASO:ORIGEN, key) # para obtener la info de cada obs. no se si la variable key vaya a servir para algo 
                                                                  
fq_purrr <- left_join(fq_mean, fq_fact, by = "SEMILLA") # juntando bases de datos con medias y con factores
fq_purrr


fq_purrr.NE<-fq_purrr[fq_purrr$ESTERIL=="NO",]
fq_purrr.E<-fq_purrr[fq_purrr$ESTERIL=="YES",]

fq_purrr.NE9<-fq_purrr.NE[fq_purrr.NE$COUNT==9,]
fq_purrr.E9<-fq_purrr.E[fq_purrr.E$COUNT==9,]
range(na.omit(fq_purrr.NE9$m_etr))
range(na.omit(fq_purrr.E9$m_etr))

ssss<-fq_purrr.NE9[-24,]
ssss<-fq_purrr.E9[-32,]
########################### COrriendo Modelos para cada Registro (Count)  ######


##### Ejemplo con m_psii

fit <-  lm(m_psii ~ ORIGEN * TIERRA ,
    data = fq_purrr,
    na.action = na.omit,
    subset = COUNT == "1"
  )

summary (fit)
anova (fit)
emmeans(fit, specs = "ORIGEN") 
emmeans(fit, pairwise ~ ORIGEN)


###################### AHORA SI.... CORRIENDO TODO PARA PSII COMO EJEMPLO 


# CORRIENDO EJEMPLO PARA PSII
#################  Anidando bases por medicion. Crea una lista de dataframes en 
##########                     #funcion para cada Medida (COUNT)


fq_by_medida_lm <- fq_purrr %>%
  group_by(COUNT) %>%
  nest()

fq_by_medida_lm.NE <- fq_purrr.NE %>%
  group_by(COUNT) %>%
  nest()

fq_by_medida_lm.E <- fq_purrr.E %>%
  group_by(COUNT) %>%
  nest()


#################   Creando modelo 
#

anova_model <- function(data){
  lm(m_psii~ ORIGEN*TIERRA , data = data, na.action = na.omit)
}


lsmeans_model <- function(data){ # no se como conseguir estas lsmeans del modelo con interaccion 
  lm(m_psii~ ORIGEN*TIERRA , data = data, na.action = na.omit)
}

################## Corriendo modelo para cada data y anidando los resultados
### 

ds_nest_anova <- fq_by_medida_lm  %>% mutate(model = map(.x = data, .f = anova_model))
ds_nest_lsmeans <- fq_by_medida_lm  %>% mutate(model = map(.x = data, .f = lsmeans_model))

################# obteniendo outputs de los modelos 

################# Extrayendo informacion de cadas modelo 
#         almacenado en la columna model 

# anovas
gral_fit_modelos_lm <- ds_nest_anova %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)

anova_modelos_lm <- ds_nest_anova %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)

# summaries 
summary_modelos_lm <- ds_nest_anova %>%
  mutate(glance_lm = map(model, broom::tidy)) %>%
  unnest(glance_lm, .drop = TRUE)


# Obteniendo means (ls means usando emmeans)
lsmeans_modelos <- ds_nest_lsmeans %>% 
  mutate(emmeans = pmap(.l = list(object = model, specs = "TIERRA", type="response"), 
                        .f = emmeans)) 

lsmeans_modelos_lm<-lsmeans_modelos %>% mutate(emm2 = map(emmeans, data.frame)) %>% 
                unnest(emm2)






# PAIRWISE
# Obteniendo means (ls means usando emmeans)

em<-emmeans(
  ds_nest_anova.NE$model[[1]], pairwise~TIERRA|ORIGEN, type="response")
em$contrasts

pair<-data.frame()
pair<-list()
for (i in c(1:9)){
  em<-emmeans(
    ds_nest_anova.NE$model[[i]], pairwise~TIERRA|ORIGEN, type="response")
  pair[[i]]<-em$contrasts
}
pair



##### Resultado finales para graficar o checar significancias 
gral_fit_modelos_lm
anova_modelos_lm$terms <- as.factor(rep(c("Origen","Tierra", "OrigenXTierra","residuales"),9))
anova_modelos_lm <- anova_modelos_lm%>%select(COUNT:model,terms, Df:terms)
na.omit(anova_modelos_lm[anova_modelos_lm$`Pr(>F)`<0.05,])
summary_modelos_lm
lsmeans_modelos_lm

## NO ESTERIL 

ds_nest_anova.NE <- fq_by_medida_lm.NE  %>% mutate(model = map(.x = data, .f = anova_model))

# anovas
gral_fit_modelos_lm.NE <- ds_nest_anova.NE %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
gral_fit_modelos_lm.NE

anova_modelos_lm.NE <- ds_nest_anova.NE %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
anova_modelos_lm.NE
anova_modelos_lm.NE$terms <- as.factor(rep(c("Origen","Tierra", "OrigenXTierra","residuales"),9))
anova_modelos_lm.NE <- anova_modelos_lm.NE%>%select(COUNT:model,terms, Df:terms)
na.omit(anova_modelos_lm.NE[anova_modelos_lm.NE$`Pr(>F)`<0.05,])



## ESTERIL 

ds_nest_anova.E <- fq_by_medida_lm.E  %>% mutate(model = map(.x = data, .f = anova_model))

# anovas
gral_fit_modelos_lm.E <- ds_nest_anova.E %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
gral_fit_modelos_lm.E

anova_modelos_lm.E <- ds_nest_anova.E %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
anova_modelos_lm.E
anova_modelos_lm.E$terms <- as.factor(rep(c("Origen","Tierra", "OrigenXTierra","residuales"),9))
anova_modelos_lm.E <- anova_modelos_lm.E%>%select(COUNT:model,terms, Df:terms)
na.omit(anova_modelos_lm.E[anova_modelos_lm.E$`Pr(>F)`<0.05,])





plm.NE<-ds_nest_anova.NE$model[[1]]
plm.E<-ds_nest_anova.E$model[[1]]


plot(ggpredict(plm.NE, terms=c("TIERRA","ORIGEN")),
     connect.lines=T)
plot(ggpredict(plm.E, terms=c("TIERRA","ORIGEN")),
     connect.lines=T)




df<-data.frame(TIERRA=1,ORIGEN=1,emmean=1,SE=1,df=1,
               lower.CL=1,upper.CL=1,COUNT=1)
for (i in c(1:9)){
  em<-as.data.frame(emmeans(
    ds_nest_anova.NE$model[[i]], specs = c("TIERRA","ORIGEN"), type="response"))
  em$COUNT=rep(i,4)
  df<-rbind(df,em)
}
df$COUNT<-as.factor(df$COUNT)
df1<-df[-1,]



ggplot(df1[1:4,], aes(COUNT, emmean, 
                 group = interaction(TIERRA, ORIGEN), 
                 color = interaction(TIERRA, ORIGEN),
                 shape = interaction(TIERRA, ORIGEN))) +
  geom_point(size=3, position=position_dodge(0.05)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("0"))+
  scale_color_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values=c("#777777", "#000000","#777777", "#000000"))+
  scale_shape_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values = c(19,1,2,17))



ggplot(df1[-c(1:4),], aes(COUNT, emmean, 
                          group = interaction(TIERRA, ORIGEN), 
                          color = interaction(TIERRA, ORIGEN),
                          shape = interaction(TIERRA, ORIGEN))) +
  geom_line(position=position_dodge(0)) + 
  geom_point(size=3, position=position_dodge(0)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0))+
  theme_classic()+
  theme(legend.position=c(0.9,0.5))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("125","190","285","420",
                             "625","845","1150","1500"))+
  scale_color_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values=c("#777777", "#000000","#777777", "#000000"))+
  scale_shape_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values = c(19,1,2,17))







### FIN FIN FIN FIN ........... HACERLO PARA CADA UNA DE LAS VARIABLES. 

g.any_anova_model <- function(data, X){
  lm(m_etr ~ m_Col+ESTERIL , data = data, na.action = na.omit)
}

g.any_anova_model <- function(data, X){
  lm(m_etr ~ m_Col+ORIGEN*TIERRA*ESTERIL , data = data, na.action = na.omit)
}


## CHOOSE ##

### PSII

any_anova_model <- function(data, X){
  lm(m_psii ~ m_Col+ORIGEN*TIERRA , data = data, na.action = na.omit)
}

### ETR

any_anova_model <- function(data, X){
  lm(m_etr ~ m_Col+ORIGEN*TIERRA , data = data, na.action = na.omit)
}

### qN
any_anova_model <- function(data, X){
  lm(m_qN ~ m_Col+ORIGEN*TIERRA , data = data, na.action = na.omit)
}

### NPQ
any_anova_model <- function(data, X){
  lm(m_npq ~ m_Col+ORIGEN*TIERRA , data = data, na.action = na.omit)
}

### qP
any_anova_model <- function(data, X){
  lm(m_qP ~ m_Col+ORIGEN*TIERRA , data = data, na.action = na.omit)
}




## RUN ##

### ALL

any_ds_nest_anova <- fq_by_medida_lm[-1,]  %>% mutate(model = map(.x = data, .f = g.any_anova_model))



# anovas
any_gral_fit_modelos_lm <- any_ds_nest_anova %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
any_gral_fit_modelos_lm

any_anova_modelos_lm <- any_ds_nest_anova %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
any_anova_modelos_lm<-na.omit(any_anova_modelos_lm)
any_anova_modelos_lm$terms <- 
  as.factor(rep(c("Tasa_Col","Soil","Seed","SeedxSoil"),16))
any_anova_modelos_lm <- any_anova_modelos_lm.NE %>%
  select(COUNT:model,terms, Df:terms)
na.omit(any_anova_modelos_lm[any_anova_modelos_lm$`Pr(>F)`<0.05,])


#
# anovas
any_gral_fit_modelos_lm <- any_ds_nest_anova %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
any_gral_fit_modelos_lm

any_anova_modelos_lm <- any_ds_nest_anova %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
any_anova_modelos_lm<-na.omit(any_anova_modelos_lm)
any_anova_modelos_lm$terms <- 
  as.factor(rep(c("Tasa_Col","Soil","Seed","Esteril",
                  "SoilxSeed","SoilxEsteril","SeedxEsteril","SxSxE"),8))
any_anova_modelos_lm <- any_anova_modelos_lm %>%
  select(COUNT:model,terms, Df:terms)
na.omit(any_anova_modelos_lm[any_anova_modelos_lm$`Pr(>F)`<0.05,])



g.any<-data.frame()
for (i in as.numeric(any_ds_nest_anova$COUNT)){
  g.em<-as.data.frame(emmeans(
    any_ds_nest_anova$model[any_ds_nest_anova$COUNT==i][[1]], 
    specs = "ESTERIL", type="response"))
  g.em$COUNT=rep(i,2)
  g.any<-rbind(g.any,g.em)
}
g.any$COUNT<-as.factor(g.any$COUNT)




g.any<-data.frame()
for (i in as.numeric(any_ds_nest_anova$COUNT)){
  g.em<-as.data.frame(emmeans(
    any_ds_nest_anova$model[any_ds_nest_anova$COUNT==i][[1]], 
    specs = c("TIERRA","ORIGEN","ESTERIL"), type="response"))
  g.em$COUNT=rep(i,2)
  g.any<-rbind(g.any,g.em)
}
g.any$COUNT<-as.factor(g.any$COUNT)



g.any<-data.frame()
for (i in as.numeric(any_ds_nest_anova$COUNT)){
  g.em<-as.data.frame(emmeans(
    any_ds_nest_anova$model[any_ds_nest_anova$COUNT==i][[1]], 
    specs = c("TIERRA","ORIGEN"), type="response"))
  g.em$COUNT=rep(i,2)
  g.any<-rbind(g.any,g.em)
}
g.any$COUNT<-as.factor(g.any$COUNT)


g.any11<-g.any[-c(1:2),]





pair<-list()
for (i in c(as.numeric(any_ds_nest_anova$COUNT))){
  em<-emmeans(
    any_ds_nest_anova$model[any_ds_nest_anova$COUNT==i][[1]], 
    pairwise~ESTERIL, type="response")
  pair[[i]]<-em$contrasts
}
pair


pair<-list()
for (i in c(as.numeric(any_ds_nest_anova$COUNT))){
  em<-emmeans(
    any_ds_nest_anova$model[any_ds_nest_anova$COUNT==i][[1]], 
    pairwise~ORIGEN|TIERRA|ESTERIL, type="response")
  pair[[i]]<-em$contrasts
}
pair


pair<-list()
for (i in c(as.numeric(any_ds_nest_anova$COUNT))){
  em<-emmeans(
    any_ds_nest_anova$model[any_ds_nest_anova$COUNT==i][[1]], 
    pairwise~TIERRA|ORIGEN|ESTERIL, type="response")
  pair[[i]]<-em$contrasts
}
pair


pair<-list()
for (i in c(as.numeric(any_ds_nest_anova$COUNT))){
  em<-emmeans(
    any_ds_nest_anova$model[any_ds_nest_anova$COUNT==i][[1]], 
    pairwise~TIERRA|ORIGEN, type="response")
  pair[[i]]<-em$contrasts
}
pair


ggplot(g.any, aes(COUNT, emmean, 
                 group = ESTERIL, 
                 color = ESTERIL)) +
  geom_line(position=position_dodge(0)) + 
  geom_point(size=3, position=position_dodge(0)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("125","190","285","420", 
                             "625","845","1150","1500"))+
  scale_color_manual(name=c("Biotic factor"),
                     labels=c(YES="Absent",  NO="Present"),
                     values=c(YES="#777777", NO="#000000"))





ggplot(g.any, aes(COUNT, emmean, 
                  group = interaction(TIERRA, ORIGEN, ESTERIL), 
                  color = interaction(TIERRA, ORIGEN, ESTERIL),
                  shape = interaction(TIERRA, ORIGEN, ESTERIL))) +
  geom_line(position=position_dodge(0.5)) + 
  geom_point(size=3, position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=1,
                position=position_dodge(0.5))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("125","190","285","420", 
                             "625","845","1150","1500"))+
  scale_color_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural.NO="Rural-Rural-NE", 
                              Urban.Urban.NO="Urban-Urban-NE",
                              Urban.Rural.NO="Rural-Urban-NE",
                              Rural.Urban.NO="Urban-Rural-NE",
                              
                              Rural.Rural.YES="Rural-Rural-E", 
                              Urban.Urban.YES="Urban-Urban-E",
                              Urban.Rural.YES="Rural-Urban-E",
                              Rural.Urban.YES="Urban-Rural-E"),
                     
                     values=c(Rural.Rural.NO="#10cedb", 
                              Urban.Urban.NO="#f8766d",
                              Urban.Rural.NO="#10cedb",
                              Rural.Urban.NO="#f8766d",
                              
                              Rural.Rural.YES="#10cedb", 
                              Urban.Urban.YES="#f8766d",
                              Urban.Rural.YES="#10cedb",
                              Rural.Urban.YES="#f8766d"))+
    
  scale_shape_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural.NO="Rural-Rural-NE", 
                              Urban.Urban.NO="Urban-Urban-NE",
                              Urban.Rural.NO="Rural-Urban-NE",
                              Rural.Urban.NO="Urban-Rural-NE",
                              
                              Rural.Rural.YES="Rural-Rural-E", 
                              Urban.Urban.YES="Urban-Urban-E",
                              Urban.Rural.YES="Rural-Urban-E",
                              Rural.Urban.YES="Urban-Rural-E"),
                     
                     values = c(Rural.Rural.NO=1, 
                                Urban.Urban.NO=1,
                                Urban.Rural.NO=2,
                                Rural.Urban.NO=2,
                                
                                Rural.Rural.YES=1, 
                                Urban.Urban.YES=1,
                                Urban.Rural.YES=2,
                                Rural.Urban.YES=2))+
  guides(color=F, shape=F)





ggplot(g.any, aes(COUNT, emmean, 
                  group = interaction(TIERRA,  ESTERIL), 
                  color = interaction(TIERRA,  ESTERIL),
                  shape = interaction(TIERRA,  ESTERIL))) +
  geom_line(position=position_dodge(0.5)) + 
  geom_point(size=3, position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=1,
                position=position_dodge(0.5))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("125","190","285","420", 
                             "625","845","1150","1500"))+
  scale_color_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.NO="Rural-Rural-NE", 
                              Urban.NO="Urban-Urban-NE",
                              
                              Rural.YES="Rural-Rural-E", 
                              Urban.YES="Urban-Urban-E"),
                     
                     values=c(Rural.NO="#10cedb", 
                              Urban.NO="#f8766d",
                              
                              Rural.YES="#10cedb", 
                              Urban.YES="#f8766d"))+
  
  scale_shape_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.NO="Rural-Rural-NE", 
                              Urban.NO="Urban-Urban-NE",
                              
                              Rural.YES="Rural-Rural-E", 
                              Urban.YES="Urban-Urban-E"),
                     
                     values = c(Rural.NO=1, 
                                Urban.NO=1,
                                
                                Rural.YES=1, 
                                Urban.YES=1))+
guides(color=F, shape=F)

  
  
  
  


###

ggplot(g.any1[1:2,], aes(COUNT, emmean, 
                   group = ESTERIL, 
                   color = ESTERIL)) +
  geom_line(position=position_dodge(0.5)) + 
  geom_point(size=3, position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.5))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="PSII", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("0"))+
  scale_color_manual(name=c("Biotic factor"),
                     labels=c(YES="Absent",  NO="Present"),
                     values=c(YES="#777777", NO="#000000"))+
  guides(color=F)








### NO ESTERIL


any_ds_nest_anova.NE <- fq_by_medida_lm.NE[-1,]  %>% mutate(model = map(.x = data, .f = any_anova_model))



# anovas
any_gral_fit_modelos_lm.NE <- any_ds_nest_anova.NE %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
any_gral_fit_modelos_lm.NE

any_anova_modelos_lm.NE <- any_ds_nest_anova.NE %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
any_anova_modelos_lm.NE$terms <- 
  as.factor(rep(c("Tasa_Col","Origen","Tierra", "OrigenXTierra","residuales"),8))
any_anova_modelos_lm.NE <- any_anova_modelos_lm.NE %>%
  select(COUNT:model,terms, Df:terms)
na.omit(any_anova_modelos_lm.NE[any_anova_modelos_lm.NE$`Pr(>F)`<0.05,])




pair<-data.frame()
pair.NE<-list()
for (i in c(as.numeric(any_ds_nest_anova.E$COUNT))){
  emNE<-emmeans(
    any_ds_nest_anova.NE$model[[i]], pairwise~TIERRA|ORIGEN, type="response")
  pair.NE[[i]]<-emNE$contrasts
}
pair.NE[[1]]
pair.NE[[2]]
pair.NE[[3]]
pair.NE[[4]]
pair.NE[[5]]
pair.NE[[6]]
pair.NE[[7]]
pair.NE[[8]]
pair.NE[[9]]


pair<-data.frame()
pair.NE<-list()
for (i in c(as.numeric(any_ds_nest_anova.E$COUNT))){
  emNE<-emmeans(
    any_ds_nest_anova.NE$model[[i]], pairwise~ORIGEN|TIERRA, type="response")
  pair.NE[[i]]<-emNE$contrasts
}
pair.NE[[1]]
pair.NE[[2]]
pair.NE[[3]]
pair.NE[[4]]
pair.NE[[5]]
pair.NE[[6]]
pair.NE[[7]]
pair.NE[[8]]
pair.NE[[9]]




### ESTERIL


any_ds_nest_anova.E <- fq_by_medida_lm.E[-1,]  %>% mutate(model = map(.x = data, .f = any_anova_model))


# anovas
any_gral_fit_modelos_lm.E <- any_ds_nest_anova.E %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
any_gral_fit_modelos_lm.E

any_anova_modelos_lm.E <- any_ds_nest_anova.E %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
any_anova_modelos_lm.E$terms <- 
  as.factor(rep(c("Tasa_Col","Origen","Tierra", "OrigenXTierra","residuales"),8))
any_anova_modelos_lm.E <- any_anova_modelos_lm.E %>%
  select(COUNT:model,terms, Df:terms)
na.omit(any_anova_modelos_lm.E[any_anova_modelos_lm.E$`Pr(>F)`<0.05,])

etrNE9<-lm(m_etr~m_Col+TIERRA,data=fq_purrr.NE9, na.action=na.omit)
etrE9<-lm(m_etr~m_Col+TIERRA,data=fq_purrr.E9,na.action=na.omit)
Anova(etrNE9, type="3")
Anova(etrE9, type="3")


pair<-data.frame()
pair.E<-list()
for (i in c(as.numeric(any_ds_nest_anova.E$COUNT))){
  emE<-emmeans(
    ds_nest_anova.NE$model[[i]], pairwise~ORIGEN|TIERRA, type="response")
  pair.E[[i]]<-emE$contrasts
}
pair.E[[1]]
pair.E[[2]]
pair.E[[3]]
pair.E[[4]]
pair.E[[5]]
pair.E[[6]]
pair.E[[7]]
pair.E[[8]]
pair.E[[9]]


pair<-data.frame()
pair.E<-list()
for (i in c(as.numeric(any_ds_nest_anova.E$COUNT))){
  emE<-emmeans(
    ds_nest_anova.NE$model[[i]], pairwise~TIERRA|ORIGEN, type="response")
  pair.E[[i]]<-emE$contrasts
}
pair.E[[1]]
pair.E[[2]]
pair.E[[3]]
pair.E[[4]]
pair.E[[5]]
pair.E[[6]]
pair.E[[7]]
pair.E[[8]]
pair.E[[9]]




lm.NE<-any_ds_nest_anova.NE$model[[8]]
lm.E<-any_ds_nest_anova.E$model[[8]]


plot(ggpredict(lm.NE, terms=c("TIERRA","ORIGEN")),
     connect.lines=T)
plot(ggpredict(lm.E, terms=c("TIERRA","ORIGEN")),
     connect.lines=T)


any_ds_nest_anova.NE$COUNT
any<-data.frame()
for (i in seq(8)){
  em<-as.data.frame(emmeans(
    any_ds_nest_anova.E$model[[i]], specs = c("TIERRA","ORIGEN"), type="response"))
  em$COUNT=rep(i,4)
  any<-rbind(any,em)
}
any$COUNT<-as.factor(any$COUNT)
any
any11<-any[-c(1:4),]



ggplot(any11, aes(COUNT, emmean, 
                          group = interaction(TIERRA, ORIGEN), 
                          color = interaction(TIERRA, ORIGEN),
                          shape = interaction(TIERRA, ORIGEN))) +
  geom_line(position=position_dodge(0)) + 
  geom_point(size=3, position=position_dodge(0)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="PSII", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("125","190","285","420", 
                             "625","845","1150","1500"))+
  scale_color_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values=c(Rural.Rural="#00ba38", 
                              Urban.Urban="#f8766d",
                              Urban.Rural="#00ba38",
                              Rural.Urban="#f8766d"))+
  scale_shape_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values = c(Rural.Rural=19, 
                                Urban.Urban=17,
                                Urban.Rural=2,
                                Rural.Urban=1))



ggplot(any[29:32,], aes(COUNT, emmean, 
                  group = interaction(TIERRA, ORIGEN), 
                  color = interaction(TIERRA, ORIGEN),
                  shape = interaction(TIERRA, ORIGEN))) +
  geom_line(position=position_dodge(0.1)) + 
  geom_point(size=3, position=position_dodge(0.1)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.1))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="PSII", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("0"))+
  scale_color_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values=c(Rural.Rural="#00ba38", 
                              Urban.Urban="#f8766d",
                              Urban.Rural="#00ba38",
                              Rural.Urban="#f8766d"))+
  scale_shape_manual(name=c("Seed-Soil origin"),
                     labels=c(Rural.Rural="Rural-Rural", 
                              Urban.Urban="Urban-Urban",
                              Urban.Rural="Rural-Urban",
                              Rural.Urban="Urban-Rural"),
                     values = c(Rural.Rural=19, 
                                Urban.Urban=17,
                                Urban.Rural=2,
                                Rural.Urban=1))+
  guides(color=F, shape=F)










################## TIERRA



any<-data.frame()
for (i in c(as.numeric(any_ds_nest_anova.E$COUNT))){
  em<-as.data.frame(emmeans(
    any_ds_nest_anova.E$model[any_ds_nest_anova.E$COUNT==i][[1]], 
    specs = c("TIERRA"), type="response"))
  em$COUNT=rep(i,2)
  any<-rbind(any,em)
}
any$COUNT<-as.factor(any$COUNT)
any
any11<-any[-c(1:2),]



ggplot(any, aes(COUNT, emmean, 
                  group = TIERRA, 
                  color = TIERRA)) +
  geom_line(position=position_dodge(0)) + 
  geom_point(size=3, position=position_dodge(0)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="ETR (Electron Transer Rate)", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("125","190","285","420", 
                             "625","845","1150","1500"))+
  scale_color_manual(name=c("Soil origin"),
                     labels=c(Rural="Rural-Rural", 
                              Urban="Urban-Urban"),
                     values=c(Rural="#00ba38", 
                              Urban="#f8766d"))


ggplot(any[1:2,], aes(COUNT, emmean, 
                       group = TIERRA, 
                       color = TIERRA)) +
  geom_line(position=position_dodge(0.1)) + 
  geom_point(size=3, position=position_dodge(0.1)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.1))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="PSII", color="Soil Origin", shape="Seed Origin")+
  scale_x_discrete(labels= c("0"))+
  scale_color_manual(name=c("Soil origin"),
                     labels=c(Rural="Rural-Rural", 
                              Urban="Urban-Urban"),
                     values=c(Rural="#00ba38", 
                              Urban="#f8766d"))+
  guides(color=F)






any<-data.frame()
for (i in c(as.numeric(any_ds_nest_anova.E$COUNT))){
  em<-as.data.frame(emmeans(
    any_ds_nest_anova.E$model[[i]], specs = c("ORIGEN"), type="response"))
  em$COUNT=rep(i,2)
  any<-rbind(any,em)
}
any$COUNT<-as.factor(any$COUNT)
any
any11<-any[-c(1:2),]



ggplot(any11, aes(COUNT, emmean, 
                  group = ORIGEN, 
                  color = ORIGEN)) +
  geom_line(position=position_dodge(0)) + 
  geom_point(size=3, position=position_dodge(0)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="PSII", color="Seed Origin")+
  scale_x_discrete(labels= c("125","190","285","420", 
                             "625","845","1150","1500"))+
  scale_color_manual(name=c("Seed origin"),
                     labels=c(Rural="Rural-Rural", 
                              Urban="Urban-Urban"),
                     values=c(Rural="#00ba38", 
                              Urban="#f8766d"))


ggplot(any[1:2,], aes(COUNT, emmean, 
                       group = ORIGEN, 
                       color = ORIGEN)) +
  geom_line(position=position_dodge(0.1)) + 
  geom_point(size=3, position=position_dodge(0.1)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.1))+
  theme_classic()+
  #theme(legend.position=c(0.9,0.2))+
  xlab(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1," photons)")) )+
  labs(y="PSII", color="Seed Origin")+
  scale_x_discrete(labels= c("0"))+
  scale_color_manual(name=c("Seed origin"),
                     labels=c(Rural="Rural-Rural", 
                              Urban="Urban-Urban"),
                     values=c(Rural="#00ba38", 
                              Urban="#f8766d"))+
  guides(color=F)

















########

### Fv_Fm
#gfv_anova_model <- function(data, X){
#  lm(m_fv_fm ~ m_Col+ESTERIL*ORIGEN*TIERRA , data = data, na.action = na.omit)
#}
gfv_anova_model <- lm(m_fv_fm ~ m_Col, 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)

gfv_anova_model <- lm(m_fv_fm ~ m_Col+ESTERIL, 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)

gfv_anova_model <- lm(m_fv_fm ~ m_Col+ESTERIL*ORIGEN*TIERRA , 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
                      )


gfv_anova_model <- lm(m_fv_fm ~ m_Col+ORIGEN*TIERRA , 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)

gfv_anova_model <- lm(m_fv_fm ~ m_Col+ORIGEN , 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)

gfv_anova_model <- lm(m_fv_fm ~ m_Col+ORIGEN*TIERRA , 
                    data = fq_by_medida_lm$data[[1]], na.action = na.omit ,
                    subset=ESTERIL=="YES"
                    )

gfv_anova_model <- lm(m_fv_fm ~ m_Col+ORIGEN*TIERRA, 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit ,
                      subset=ESTERIL=="NO"
                      )

summary(gfv_anova_model)
Anova(gfv_anova_model, type="2")
Anova(gfv_anova_model, type="3")

emmeans(gfv_anova_model, pairwise~ORIGEN,  type="response")


fvemm<-as.data.frame(emmeans(gfv_anova_model, spec=c("TIERRA","ORIGEN","ESTERIL")))
fvemm<-as.data.frame(emmeans(gfv_anova_model, spec=c("ORIGEN")))
emmeans(gfv_anova_model, pairwise~TIERRA|ORIGEN|ESTERIL,  type="response")
emmeans(gfv_anova_model, pairwise~ORIGEN|TIERRA|ESTERIL,  type="response")

emmeans(gfv_anova_model, pairwise~ESTERIL|ORIGEN|TIERRA,  type="response")




gfv_anova_model <- lm(m_fv_fm ~ ORIGEN , 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)
fv_SP<-summary(gfv_anova_model)
anova(gfv_anova_model)
cbind(ORIGEN=c("rural", "urban"),
      as.data.frame(cbind(mean=c(fv_SP$coefficients[1],
                                 fv_SP$coefficients[1]+
                                   fv_SP$coefficients[2]),
                          se=c(fv_SP$coefficients[3],
                               fv_SP$coefficients[4]))))

FVFM_solo<-fq_by_medida_lm$data[[1]]

FVFM_solo.E<-FVFM_solo[FVFM_solo$ESTERIL=="YES",]
FVFM_solo.NE<-FVFM_solo[FVFM_solo$ESTERIL=="NO",]

#ESTERIL
with(FVFM_solo.E, tapply(m_fv_fm, TIERRA, length))
with(FVFM_solo.E, tapply(m_fv_fm, ORIGEN, length))

fvfm_E.SU<-FVFM_solo.E[FVFM_solo.E$ORIGEN=="Urban",]
fvfm_E.SR<-FVFM_solo.E[FVFM_solo.E$ORIGEN=="Rural",]

0.007217/sqrt(as.vector(with(fvfm_E.SU, tapply(m_fv_fm, TIERRA, length))))
0.007217/sqrt(as.vector(with(fvfm_E.SR, tapply(m_fv_fm, TIERRA, length))))




#NO ESTERIL
with(FVFM_solo.NE, tapply(m_fv_fm, TIERRA, length))
with(FVFM_solo.NE, tapply(m_fv_fm, ORIGEN, length))

fvfm_NE.SU<-FVFM_solo.NE[FVFM_solo.NE$ORIGEN=="Urban",]
fvfm_NE.SR<-FVFM_solo.NE[FVFM_solo.NE$ORIGEN=="Rural",]

0.007217/sqrt(as.vector(with(fvfm_NE.SU, tapply(m_fv_fm, TIERRA, length))))
0.007217/sqrt(as.vector(with(fvfm_NE.SR, tapply(m_fv_fm, TIERRA, length))))




ES<-cohen.d(m_fv_fm ~ ESTERIL, data = fq_by_medida_lm$data[[1]], 
        na.action = na.omit, 
        hedges.correction=T)
View(ES)
ES[[4]]

ES<-cohen.d(m_fv_fm ~ ORIGEN, data = fq_by_medida_lm$data[[1]], 
        na.action = na.omit, 
        hedges.correction=T)

ggplot(fvemm, aes(TIERRA, emmean, 
                      group = ORIGEN, 
                      color = ORIGEN)) +
  geom_point(size=3, position=position_dodge(0.05)) + geom_line()+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(legend.position=c(0.9,0.2))+
  labs(y="Fv / Fm", color="Origen de Semilla")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))

fvemm<-as.data.frame(emmeans(gfv_anova_model, spec=c("TIERRA","ESTERIL")))
emmeans(gfv_anova_model, pairwise~TIERRA|ESTERIL)
emmeans(gfv_anova_model, pairwise~ESTERIL|TIERRA)


ggplot(fvemm, aes(TIERRA, emmean, 
                  group = ESTERIL, 
                  color = ESTERIL)) +
  geom_point(size=3, position=position_dodge(0.05)) + geom_line()+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(legend.position=c(0.9,0.2))+
  labs(y="Fv / Fm", color="Biotic factor")+
  scale_color_manual(labels=c(NO="Present", YES="Absent"),
                     values=c(YES="#aaaaaa", NO="#000000"))







fvemm<-as.data.frame(emmeans(gfv_anova_model, spec=c("ORIGEN","ESTERIL")))
emmeans(gfv_anova_model, pairwise~ESTERIL|ORIGEN)
emmeans(gfv_anova_model, pairwise~ORIGEN|ESTERIL)



ggplot(fvemm, aes(ESTERIL, emmean, 
                  group = ORIGEN, 
                  color = ORIGEN)) +
  geom_point(size=3, position=position_dodge(0.05)) + 
  geom_line(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(legend.position=c(0.9,0.2))+
  labs(y="Fv / Fm", x="Biotic factor",
       color="Seed origin")+
  scale_x_discrete(labels= c("Present", "Absent"))+
  scale_color_manual(values=c("#00ba38", "#f8766d"))





fvemm<-as.data.frame(emmeans(gfv_anova_model, spec=c("TIERRA","ESTERIL")))
emmeans(gfv_anova_model, pairwise~TIERRA|ESTERIL)
emmeans(gfv_anova_model, pairwise~ESTERIL|TIERRA)


ggplot(fvemm, aes(TIERRA, emmean, 
                  group = ESTERIL, 
                  color = ESTERIL)) +
  geom_point(size=3, position=position_dodge(0.05)) + geom_line()+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(legend.position=c(0.9,0.2))+
  labs(y="Fv / Fm", color="Biotic factor")+
  scale_color_manual(labels=c(NO="Present", YES="Absent"),
                     values=c(YES="#aaaaaa", NO="#000000"))








fvemm<-as.data.frame(emmeans(gfv_anova_model, 
                             spec=c("ESTERIL", "TIERRA","ORIGEN")))

emmeans(gfv_anova_model, pairwise~TIERRA|ORIGEN|ESTERIL)
emmeans(gfv_anova_model, pairwise~ORIGEN|TIERRA|ESTERIL)

eff_fvemm<-emmeans(gfv_anova_model, 
                   spec=c("ESTERIL", "TIERRA","ORIGEN"))

EFF_FVFM<-eff_size(eff_fvemm, 
         sigma = sqrt(mean(sigma(gfv_anova_model)^2)),   # RMS of sigma()
         edf = df.residual(gfv_anova_model))
df.EFF_FVFM<-as.data.frame(EFF_FVFM)

df.EFF_FVFM.O<-df.EFF_FVFM[c(2,24,9,27),]
df.EFF_FVFM.T<-df.EFF_FVFM[c(4,17,11,22),]

df.EFF_FVFM.O$ESTERIL<-c("NO","NO","YES","YES")
df.EFF_FVFM.O$ORIGEN<-c("Rural", "Urban","Rural", "Urban")

df.EFF_FVFM.T$ESTERIL<-c("NO","NO","YES","YES")
df.EFF_FVFM.T$TIERRA<-c("Rural", "Urban","Rural", "Urban")



ggplot(fvemm, aes(TIERRA, emmean, 
                    group = ORIGEN, 
                    color = ORIGEN)) +
  geom_point(size=3, position=position_dodge(0.1)) + 
  geom_line(position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1,
                position=position_dodge(0.1))+
  theme_classic()+
  facet_grid(. ~ ESTERIL)+
  theme(legend.position=c(0.9,0.2),)+
  labs(y="Fv/Fm", x="Soil origin", color="Seed origin")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))



################################################################

fv_anova_model <- function(data, X){
  lm(m_fv_fm ~ ESTERIL , data = data, na.action = na.omit)
}

fv_anova_model <- function(data, X){
  lm(m_fv_fm ~ ORIGEN*TIERRA , data = data, na.action = na.omit)
}

fv_ds_nest_anova <- fq_by_medida_lm[1,]  %>% mutate(model = map(.x = data, .f = fv_anova_model))

# anovas
fv_gral_fit_modelos_lm <- fv_ds_nest_anova %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
fv_gral_fit_modelos_lm

fv_anova_modelos_lm <- fv_ds_nest_anova %>%
  mutate(glance_lm = map(model, Anova)) %>%
  unnest(glance_lm, .drop = TRUE)
fv_anova_modelos_lm$terms <- 
  as.factor(rep(c("Origen","Tierra", "OrigenXTierra","residuales"),1))
fv_anova_modelos_lm$terms <- 
  as.factor(rep(c("Esteril","residuales"),1))
fv_anova_modelos_lm <- fv_anova_modelos_lm %>%
  select(COUNT:model,terms, Df:terms)
fv_anova_modelos_lm


fv_ds_nest_anova.NE <- fq_by_medida_lm.NE[1,]  %>% mutate(model = map(.x = data, .f = fv_anova_model))


# anovas
fv_gral_fit_modelos_lm.NE <- fv_ds_nest_anova.NE %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
fv_gral_fit_modelos_lm.NE

fv_anova_modelos_lm.NE <- fv_ds_nest_anova.NE %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
fv_anova_modelos_lm.NE$terms <- 
  as.factor(rep(c("Origen","Tierra", "OrigenXTierra","residuales"),1))
fv_anova_modelos_lm.NE <- fv_anova_modelos_lm.NE %>%
  select(COUNT:model,terms, Df:terms)
fv_anova_modelos_lm.NE



fv_ds_nest_anova.E <- fq_by_medida_lm.E[1,]  %>% mutate(model = map(.x = data, .f = fv_anova_model))


# anovas
fv_gral_fit_modelos_lm.E <- fv_ds_nest_anova.E %>%
  mutate(glance_lm = map(model, broom::glance)) %>%
  unnest(glance_lm, .drop = TRUE)
fv_gral_fit_modelos_lm.E

fv_anova_modelos_lm.E <- fv_ds_nest_anova.E %>%
  mutate(glance_lm = map(model, anova)) %>%
  unnest(glance_lm, .drop = TRUE)
fv_anova_modelos_lm.E$terms <- 
  as.factor(rep(c("Origen","Tierra", "OrigenXTierra","residuales"),1))
fv_anova_modelos_lm.E <- fv_anova_modelos_lm.E %>%
  select(COUNT:model,terms, Df:terms)
fv_anova_modelos_lm.E



plot(ggpredict(fv_ds_nest_anova$model[[1]], terms=c("TIERRA","ORIGEN","ESTERIL")),
     connect.lines=T)

plot(ggpredict(fv_ds_nest_anova.NE$model[[1]], terms=c("TIERRA","ORIGEN")),
     connect.lines=T)
plot(ggpredict(fv_ds_nest_anova.E$model[[1]], terms=c("TIERRA","ORIGEN")),
     connect.lines=T)




emmeans(fv_ds_nest_anova.NE$model[[1]], specs = c("TIERRA","ORIGEN"))
emmeans(fv_ds_nest_anova.E$model[[1]], specs = c("TIERRA","ORIGEN"))





# Un ejemplo que nos ayudaria a obtener los estimados de cada variable de un modelo 
# que considera interacciones... 
# tomado de https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#contrasts
noise.lm <- lm(noise/10 ~ size * type * side, data = auto.noise)
emmeans(noise.lm, specs = c("size", "type"))
emmeans(noise.lm, pairwise ~ size)




gcol_anova_model <- lm(m_Col ~ ESTERIL*ORIGEN*TIERRA , 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)
anova(gcol_anova_model)



gcol<-as.data.frame(emmeans(gcol_anova_model, 
                            spec=c("TIERRA","ORIGEN","ESTERIL"),type="response"))
gcol
emmeans(gfv_anova_model, pairwise~TIERRA|ESTERIL)
emmeans(gfv_anova_model, pairwise~ESTERIL|TIERRA)



ggplot(gcol, aes(TIERRA, emmean, 
                  color = ORIGEN,
                  group = ORIGEN)) +
  geom_point(size=3, position=position_dodge(0.1)) + 
  geom_line(position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1,
                position=position_dodge(0.1))+
  theme_classic()+
  facet_grid(. ~ ESTERIL)+
  theme(legend.position=c(0.95,0.2),)+
  labs(y="% Colonization", x="Soil origin", fill="Seed origin")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))




gcol_anova_modelNE <- lm(m_Col ~ ORIGEN*TIERRA, 
                         data = fq_by_medida_lm$data[[1]], na.action = na.omit ,subset=ESTERIL=="NO"
)

gcol_anova_modelE <- lm(m_Col ~ ORIGEN*TIERRA , 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit ,subset=ESTERIL=="YES"
)


anova(gcol_anova_modelNE)
anova(gcol_anova_modelE)

gcolNE<-as.data.frame(emmeans(gcol_anova_modelNE, 
                            spec=c("TIERRA","ORIGEN"),type="response"))

emmeans(gcol_anova_modelNE, pairwise~ORIGEN|TIERRA)
emmeans(gcol_anova_modelNE, pairwise~TIERRA|ORIGEN)

ggplot(gcolNE, aes(TIERRA, emmean, 
                 color = ORIGEN,
                 group = ORIGEN)) +
  geom_point(size=3, position=position_dodge(0.1)) + 
  geom_line(position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1,
                position=position_dodge(0.1))+
  theme_classic()+
  theme(legend.position=c(0.95,0.2),)+
  labs(y="% Colonization", x="Soil origin", fill="Seed origin")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))




####################################
#                                  #
#            SIZE EFFECT           #
#                                  #
####################################

gfv_anova_model <- lm(m_fv_fm ~ m_prop01, 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)
summary(gfv_anova_model)
anova(gfv_anova_model)

##### Prepare the effect sizes -- Rates of events #####

fq_eff_NE<-fq_by_medida_lm.NE[1,2][[1]][[1]]
fq_eff_E<-fq_by_medida_lm.E[1,2][[1]][[1]]


fq_eff_E[fq_eff_E$TIERRA=="Urban",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[[4]]*-1,
              col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[[4]]*-1),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                               na.action = na.omit, 
                                               hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]]*-1,
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]]*-1),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]]*-1,
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]]*-1)
  ) -> fq.effsize_Local.E.TU


fq_eff_E[fq_eff_E$TIERRA=="Rural",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                   na.action = na.omit, 
                                   hedges.correction=T)[[4]],
              col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]]),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                                 na.action = na.omit, 
                                                 hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]],
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[1]]),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]],
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[2]])
  ) -> fq.effsize_Local.E.TR



fq_eff_E[fq_eff_E$ORIGEN=="Urban",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                   na.action = na.omit, 
                                   hedges.correction=T)[[4]]*-1,
              col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]]*-1),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~TIERRA,  data = . ,
                                                 na.action = na.omit, 
                                                 hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]]*-1,
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[1]]*-1),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]]*-1,
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[2]]*-1)
  ) -> fq.effsize_Local.E.SU



fq_eff_E[fq_eff_E$ORIGEN=="Rural",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                   na.action = na.omit, 
                                   hedges.correction=T)[[4]],
              col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]]),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~TIERRA,  data = . ,
                                                 na.action = na.omit, 
                                                 hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]],
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[1]]),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]],
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[2]])
  ) -> fq.effsize_Local.E.SR


fq.effsize_Local.E.TU$ZONE<-"Urban"
fq.effsize_Local.E.TR$ZONE<-"Rural"
fq.effsize_Local.E.SU$ZONE<-"Urban"
fq.effsize_Local.E.SR$ZONE<-"Rural"


names(fq.effsize_Local.E.TU)[c(3,4)]<-c("upper","lower")
names(fq.effsize_Local.E.SU)[c(3,4)]<-c("upper","lower")


fq.effsize_Local.E.T<-rbind(fq.effsize_Local.E.TU,
                         fq.effsize_Local.E.TR)

fq.effsize_Local.E.S<-rbind(fq.effsize_Local.E.SU,
                         fq.effsize_Local.E.SR)


fq.effsize_Local.E.Tfvfm<-fq.effsize_Local.E.T[c(1,3),-c(2)]
fq.effsize_Local.E.Tcol<-fq.effsize_Local.E.T[c(2,4),-c(2)]

fq.effsize_Local.E.Sfvfm<-fq.effsize_Local.E.S[c(1,3),-c(2)]
fq.effsize_Local.E.Scol<-fq.effsize_Local.E.S[c(2,4),-c(2)]


names(fq.effsize_Local.E.Tfvfm)[1:3]<-c("TIERRA","T.lower","T.upper")
names(fq.effsize_Local.E.Tcol)[1:3]<-c("TIERRA","T.lower","T.upper")

names(fq.effsize_Local.E.Sfvfm)[1:3]<-c("ORIGEN","O.lower","O.upper")
names(fq.effsize_Local.E.Scol)[1:3]<-c("ORIGEN","O.lower","O.upper")


fq.effsize_Local.Efvfm<-left_join(fq.effsize_Local.E.Tfvfm,fq.effsize_Local.E.Sfvfm, by="ZONE")
fq.effsize_Local.Ecol<-left_join(fq.effsize_Local.E.Tcol,fq.effsize_Local.E.Scol, by="ZONE")




# NO ESTERIL

fq_eff_NE[fq_eff_NE$TIERRA=="Urban",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                   na.action = na.omit, 
                                   hedges.correction=T)[[4]]*-1,
              col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]]*-1),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                                 na.action = na.omit, 
                                                 hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]]*-1,
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[1]]*-1),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]]*-1,
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[2]]*-1)
  ) -> fq.effsize_Local.NE.TU


fq_eff_NE[fq_eff_NE$TIERRA=="Rural",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                   na.action = na.omit, 
                                   hedges.correction=T)[[4]],
              col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]]),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                                 na.action = na.omit, 
                                                 hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]],
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[1]]),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~ORIGEN,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]],
            col_effsize=cohen.d(m_prop01 ~ORIGEN,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[2]])
  ) -> fq.effsize_Local.NE.TR



fq_eff_NE[fq_eff_NE$ORIGEN=="Urban",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                   na.action = na.omit, 
                                   hedges.correction=T)[[4]]*-1,
              col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]]*-1),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~TIERRA,  data = . ,
                                                 na.action = na.omit, 
                                                 hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]]*-1,
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[1]]*-1),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]]*-1,
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[2]]*-1)
  ) -> fq.effsize_Local.NE.SU



fq_eff_NE[fq_eff_NE$ORIGEN=="Rural",] %>%  
  summarise(
    effsize=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                   na.action = na.omit, 
                                   hedges.correction=T)[[4]],
              col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                  na.action = na.omit, 
                                  hedges.correction=T)[[4]]),
    
    
    magnitude=c(fvfm_effsize=as.character(cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                                  na.action = na.omit, 
                                                  hedges.correction=T)[["magnitude"]]),
                col_effsize=as.character(cohen.d(m_prop01 ~TIERRA,  data = . ,
                                                 na.action = na.omit, 
                                                 hedges.correction=T)[["magnitude"]])),
    
    
    lower=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[1]],
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[1]]),
    
    upper=c(fvfm_effsize=cohen.d(m_fv_fm ~TIERRA,  data = . ,
                                 na.action = na.omit, 
                                 hedges.correction=T)[["conf.int"]][[2]],
            col_effsize=cohen.d(m_prop01 ~TIERRA,  data = . ,
                                na.action = na.omit, 
                                hedges.correction=T)[["conf.int"]][[2]])
  ) -> fq.effsize_Local.NE.SR


fq.effsize_Local.NE.TU$ZONE<-"Urban"
fq.effsize_Local.NE.TR$ZONE<-"Rural"
fq.effsize_Local.NE.SU$ZONE<-"Urban"
fq.effsize_Local.NE.SR$ZONE<-"Rural"


names(fq.effsize_Local.NE.TU)[c(3,4)]<-c("upper","lower")
names(fq.effsize_Local.NE.SU)[c(3,4)]<-c("upper","lower")


fq.effsize_Local.NE.T<-rbind(fq.effsize_Local.NE.TU,
                            fq.effsize_Local.NE.TR)

fq.effsize_Local.NE.S<-rbind(fq.effsize_Local.NE.SU,
                            fq.effsize_Local.NE.SR)


fq.effsize_Local.NE.Tfvfm<-fq.effsize_Local.NE.T[c(1,3),-c(2)]
fq.effsize_Local.NE.Tcol<-fq.effsize_Local.NE.T[c(2,4),-c(2)]

fq.effsize_Local.NE.Sfvfm<-fq.effsize_Local.NE.S[c(1,3),-c(2)]
fq.effsize_Local.NE.Scol<-fq.effsize_Local.NE.S[c(2,4),-c(2)]


names(fq.effsize_Local.NE.Tfvfm)[1:3]<-c("TIERRA","T.lower","T.upper")
names(fq.effsize_Local.NE.Tcol)[1:3]<-c("TIERRA","T.lower","T.upper")

names(fq.effsize_Local.NE.Sfvfm)[1:3]<-c("ORIGEN","O.lower","O.upper")
names(fq.effsize_Local.NE.Scol)[1:3]<-c("ORIGEN","O.lower","O.upper")


fq.effsize_Local.NEfvfm<-left_join(fq.effsize_Local.NE.Tfvfm,fq.effsize_Local.NE.Sfvfm, by="ZONE")
fq.effsize_Local.NEcol<-left_join(fq.effsize_Local.NE.Tcol,fq.effsize_Local.NE.Scol, by="ZONE")

#

fq.effsize_Local.Efvfm$ESTERIL <- "YES"
fq.effsize_Local.Ecol$ESTERIL <- "YES"


fq.effsize_Local.NEfvfm$ESTERIL <- "NO"
fq.effsize_Local.NEcol$ESTERIL <- "NO"


fq.effsize_Local.fvfm<-rbind(fq.effsize_Local.Efvfm, fq.effsize_Local.NEfvfm)
fq.effsize_Local.col<-rbind(fq.effsize_Local.Ecol, fq.effsize_Local.NEcol)

#####

gfv_anova_model <- lm(m_fv_fm ~ m_prop01+ESTERIL*ORIGEN*TIERRA , 
                      data = fq_by_medida_lm$data[[1]], na.action = na.omit
)
anova(gfv_anova_model)

fvemm<-as.data.frame(emmeans(gfv_anova_model, 
                             spec=c("ESTERIL", "TIERRA","ORIGEN")))

emmeans(gfv_anova_model, pairwise~TIERRA|ORIGEN|ESTERIL)
emmeans(gfv_anova_model, pairwise~ORIGEN|TIERRA|ESTERIL)

eff_fvemm<-emmeans(gfv_anova_model, 
                   spec=c("ESTERIL", "TIERRA","ORIGEN"))

EFF_FVFM<-eff_size(eff_fvemm, 
                   sigma = sqrt(mean(sigma(gfv_anova_model)^2)),   # RMS of sigma()
                   edf = df.residual(gfv_anova_model))
df.EFF_FVFM<-as.data.frame(EFF_FVFM)

df.EFF_FVFM.O<-df.EFF_FVFM[c(2,9),-c(1,3,4)]
df.EFF_FVFM.T<-df.EFF_FVFM[c(4,11),-c(1,3,4)]

df.EFF_FVFM.O1<-df.EFF_FVFM[c(24,27),-c(1,3,4)]*-1
df.EFF_FVFM.T1<-df.EFF_FVFM[c(17,22),-c(1,3,4)]*-1

names(df.EFF_FVFM.O1)[c(3,2)]<-c("lower.CL","upper.CL")
names(df.EFF_FVFM.T1)[c(3,2)]<-c("lower.CL","upper.CL")

df.EFF_FVFM.O<-rbind(df.EFF_FVFM.O,df.EFF_FVFM.O1)
df.EFF_FVFM.T<-rbind(df.EFF_FVFM.T,df.EFF_FVFM.T1)


df.EFF_FVFM.O$ESTERIL<-c("NO","YES","NO","YES")
df.EFF_FVFM.O$ZONE<-c("Rural","Rural", "Urban", "Urban")

df.EFF_FVFM.T$ESTERIL<-c("NO","YES","NO","YES")
df.EFF_FVFM.T$ZONE<-c("Rural","Rural", "Urban", "Urban")

names(df.EFF_FVFM.O)<-c("ORIGEN","O.lower","O.upper","ESTERIL","ZONE")
names(df.EFF_FVFM.T)<-c("TIERRA","T.lower","T.upper","ESTERIL","ZONE")

fq.effsize_Local.fvfm<-left_join(df.EFF_FVFM.T,df.EFF_FVFM.O,by=c("ZONE","ESTERIL"))


###COLONIZATION

#tot_01<-glm(cbind(rtaCOL[,1],rtaCOL[,2])~tiemp_h+ESTERIL*TIERRA*ORIGEN, data = Experim_UR, family=quasibinomial)
tot_01.NE<-glm(rtaCOL.NE~tiemp_h+TIERRA*ORIGEN, data = Experim_UR.NE,
               family=quasibinomial)
Anova(tot_01.NE, type="2")

eff_colemm<-emmeans(tot_01.NE, 
                   spec=c("TIERRA","ORIGEN"))

EFF_COL<-eff_size(eff_colemm, 
                   sigma = sqrt(mean(sigma(tot_01.NE)^2)),   # RMS of sigma()
                   edf = df.residual(tot_01.NE))
df.EFF_COL<-as.data.frame(EFF_COL)

df.EFF_COL.O<-df.EFF_COL[c(1),-c(1,3,4)]
df.EFF_COL.T<-df.EFF_COL[c(2),-c(1,3,4)]

df.EFF_COL.O1<-df.EFF_COL[c(6),-c(1,3,4)]*-1
df.EFF_COL.T1<-df.EFF_COL[c(5),-c(1,3,4)]*-1

names(df.EFF_COL.O1)[c(2,3)]<-c("asymp.UCL", "asymp.LCL")
names(df.EFF_COL.T1)[c(2,3)]<-c("asymp.UCL", "asymp.LCL")

df.EFF_COL.O<-rbind(df.EFF_COL.O,df.EFF_COL.O1)
df.EFF_COL.T<-rbind(df.EFF_COL.T,df.EFF_COL.T1)

df.EFF_COL.O$ZONE<-c("Rural", "Urban")
df.EFF_COL.T$ZONE<-c("Rural", "Urban")

names(df.EFF_COL.O)<-c("ORIGEN","O.lower","O.upper","ZONE")
names(df.EFF_COL.T)<-c("TIERRA","T.lower","T.upper","ZONE")

fq.effsize_Local.col<-left_join(df.EFF_COL.T,df.EFF_COL.O,by=c("ZONE"))


#Data for size effects
fq.effsize_Local.fvfm
fq.effsize_Local.col



#FULL Fv/Fm
ggplot(fq.effsize_Local.fvfm, aes(ZONE, TIERRA, color=ZONE)) +
  geom_point(size=5,shape=18, position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=T.lower, ymax=T.upper), width=0.1,
                position=position_dodge(0.1))+
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "grey", size=1)+
  theme_classic()+
  facet_grid(.~ESTERIL)+
  theme(legend.position=c(0.95,0.2),)+
  labs(y="Effectsize", x="Soil Origin")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))+
  guides(color=F)

ggplot(fq.effsize_Local.fvfm, aes(ZONE,ORIGEN, color=ZONE)) +
  geom_point(size=5,shape=18, position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=O.lower, ymax=O.upper), width=0.1,
                position=position_dodge(0.1))+
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "grey", size=1)+
  theme_classic()+
  facet_grid(.~ESTERIL)+
  theme(legend.position=c(0.95,0.2),)+
  labs(y="Effectsize", x="Seed Origin")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))+
  guides(color=F)


#Colonization rate

ggplot(fq.effsize_Local.col, aes(ZONE, TIERRA, color=ZONE)) +
  geom_point(size=5,shape=18, position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=T.lower, ymax=T.upper), width=0.1,
                position=position_dodge(0.1))+
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "grey", size=1)+
  theme_classic()+
  theme(legend.position=c(0.95,0.2),)+
  labs(y="Effectsize", x="Soil Origin")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))+
  guides(color=F)

ggplot(fq.effsize_Local.col, aes(ZONE, ORIGEN, color=ZONE)) +
  geom_point(size=5,shape=18, position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=O.lower, ymax=O.upper), width=0.1,
                position=position_dodge(0.1))+
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "grey", size=1)+
  theme_classic()+
  theme(legend.position=c(0.95,0.2),)+
  labs(y="Effectsize", x="Seed Origin")+
  scale_color_manual(values=c("#00ba38", "#f8766d"))+
  guides(color=F)

#####






################################################################################
#
#                   Usando la aproximacion de modelos mixtos 
#
################################################################################



############################ Modelo mixto con SEMILLAS as random 
# CORRIENDO EJEMPLO PARA PSII
#################  Anidando bases por medicion. Crea una lista de dataframes en 
##########                     #funcion para cada Medida (COUNT)

fq_by_medida <- fq %>%
  group_by(COUNT) %>%
  nest()

#################   Creando modelo 
# Creando modelo, decidi usar a SEMILLA como random para controlar por pseudo replica.. .y si es sig el termino!
fq_model <- function(df) {
  lmer(PSII~ ORIGEN*TIERRA+ (1|SEMILLA), data = df, na.action = na.omit)
}

# paso intermedio .....
# para entender como funciona... 
#models<- map(fq_by_medida$data, fq_model)
#anova (models[[1]])

################## Corriendo modelo para cada data
### 

by_medida <- fq_by_medida %>%
  mutate( model = map(data, fq_model))

################# Extrayendo informacion de cadas modelo 
#         almacenado en la columna model 

#### Obteniendo estimados generales del modelo 

glance_estimados_modelos <- by_medida %>%
  mutate(glance = map(model, broom::glance)) %>%
  unnest(glance, .drop = TRUE)

glance_estimados_modelos


#####         obteniendo ANOVAS 
#para cada uno de los terminos del modelo para cada toma 
# de medida

glance_anovas_temrinos<-by_medida %>%
  mutate(glance = map(model, anova)) %>%   # hay que ver como logramos meter suma de cuadrados III
  unnest(glance, .drop = TRUE)  # habria qe agregar una col con rep(c("Origen","Tierra","OxT"), 9))
glance_anovas_terminos # ya esta la base con las anovas para cada termino!


# Buscando un termino significativo 
glance_anovas %>%filter(`Pr(>F)`< 0.05) # ahora queremos las filas que sean significativas ...  
                                        # y solo es una :(














