

#########--------------Structural equation modeling (SEM)-----------############

########  Loading packages  #####
library (tidyverse)
library(Hmisc)
library (lavaan) # para hacer SEM (path analysis)
library (semTools) # para correr measurementInvariance
library (tidySEM) # para graficar path diagrams de SEM
library (corrg)  # para analisiar patrones de correl 
library (PerformanceAnalytics) # grafico de correlaciones ver... https://rb.gy/akjidd
library(mvnormtest) # para probar multinormalidad
library(ggpubr) # permite incorporar varios graficos en un solo espacio de graf
library (emmeans)

####### reading database   ##############
Muestreo_summary_UR <- read_csv("Database/Muestreo_Summary1.csv")[,-1]

mur <- Muestreo_summary_UR %>% mutate_if(is.character, as.factor) # 

View (mur)


################## Creando PCA                  ################################

pca_plant <- prcomp(mur[,c("long_root", "no_roots", "volumen_ml")], scale = TRUE)
summary(pca_plant)
pca_plant

pca1_roots <- pca_plant[["x"]][, 1] # obteniendo eigen values PCA1
pca2_roots <- pca_plant[["x"]][, 2] # obteniendo eigen values PCA2

mur$pca1_roots <- pca1_roots # solo trabajaremos con PCA1 como estimador de 
# tamano de raiz


######## Carpinteria de datos ################

mur_cor_all <-
  mur %>% select (
    ambient, pca1_roots,
    volumen_ml,
    prop01,
    w_above:w_s_area,
    pH:P,
    dense,
    SHANNON) %>% mutate(P.10 = P/100, Ca.10 = Ca / 100, 
                        t.prop01 = scale(asin(sqrt(prop01))),
                        t.pH = scale(log(pH)),
                        t.N = scale(log(N)),
                        t.K = scale(log(K)),
                        t.Na = scale(log(Na)),
                        t.Ca = scale(log(Ca.10)),
                        t.P = scale(log(P.10)))  # scale(variables) a priori ver 
                                                 #. pagina 191 Manual EQS ultimo parrafo

# " however, for most nodels it is permissible in practice for the variables in 
# one sample to be rescaled so that the resulting covariance matrix is close to 
# (or, cheating somewhat, identical to) a correlation matrix, providing that the
# variables in all other samples are rescaled by the scaling constants form the 
# reference sample" Bentler EQS manual 

# Sobre estandarizaciones y no estandarizaciones ... 
# https://jslefche.github.io/sem_book/coefficients.html#unstandardized-and-standardized-coefficients

# P.10 y Ca.10 resescala variables dividiendolas x100, por lo que sus unidades cambiaran 
# y hay que pensar como debran ser leidas las nuevas unidades. Esto es porque en
# analisis subsecuentes con SEM indica que las varianzas son extremas entre variables
# Ej. P tiene una varianza de 630 mientras que Na de 0.0000798. 
# PARA PROBLEMAS CON magnitudes de varianza entre variables... 
# ver este link https://groups.google.com/g/lavaan/c/FzShzsXYNco?pli=1


# aqui muestro el efecto las diferencias en varianzas 

varianzas <- mur_cor_all %>% group_by(ambient) %>% 
  summarise_if(is.numeric, list(var))
varianzas # fijate el contraste en los ordenes de magnitud de P y P.10

                       

########### Pairwise corr para variables importantes para el SEM:   ############

# Antes del SEM, debemos de entender los patrones de covariacion entre las var
# que usaremos. Creo que por lo que hemos platicado, variables como diversidad
# o conteos de esporas quedan fuera de los analisis del SEM. Y nos concentraremos
# en variables abioticas y usaremos el tama;o de la raiz como covariable
# para eliminar efectos numericos (i.e. grandes raices tienen mas colonizacion)


### Agrupando todas las observaciones: agrupa todos los datos  

### Nota
# esta es la matriz de obs que usaremos en el SEM general 
# para separar los efectos indirectos producto de la asociacion entre tantas 
# variables. 


####### Corr. on No transformed variables  ##########
chart.Correlation(mur_cor_all[, c("prop01", "pH", "N", "K"
                                  , "Na", "Ca","P" ,"pca1_roots")],
                  histogram = TRUE, pch = 19) # sin transformar
####### Corr. on transformed variables    #############
chart.Correlation(mur_cor_all[, c("t.prop01", "t.pH","t.N","t.K"
                                             ,"t.Na","t.Ca", "P","pca1_roots")],
                  histogram = TRUE, pch = 19) # transformados... 

#### Separados por ambientes usando transformed variables 

mur_cor_urb_soil <- mur_cor_all %>% filter(ambient == "urban")
mur_cor_rur_soil <- mur_cor_all %>% filter(ambient == "rural")

# Urban
chart.Correlation(mur_cor_urb_soil[, c("t.prop01", "t.pH","t.N","t.K"
                                       ,"t.Na","t.Ca","t.P", "pca1_roots")], histogram = TRUE, pch = 19) # transformados... 
# Rural
chart.Correlation(mur_cor_rur_soil[, c("t.prop01", "t.pH","t.N","t.K"
                                       ,"t.Na","t.Ca","t.P","pca1_roots")], histogram = TRUE, pch = 19) # transformados... 


######## Testing multivariate normality -----------------------------------------
# NO SE PORQUE MARCA ERROR... 
# mshapiro.test(t(mur_cor_all[,c("prop01", "pH")))  # hay que checar como correrlo, yo diria con correrlo para la matriz general
                               # que incluye todos los ambientes. 
# no la reporte en el paper

##########                 SEM        ###############
################################################################################
#
############                  Modelado SEM    #################
#
################################################################################

# De manera general, el SEM trabaja con las matrices observadas de covariacion 
# que pueden ser observadas en los charts de pairwise correlation y trata de 
# explicar la variacion de estas covarianzas pero incorporando una direccionalidad 
# causal la cual permite a su vez separar efectos directos de indirectos. 

# en este primer modelo la matriz a explicar es de 1x1 (prop01 vs K) 
# y como el modelo solo incluye esa regresion es un modelo saturado ya que 
# el modelo que se trata de probar tiene como unico parametro el mismo parametro que 
# quiere ser evaluado


modelK<-'				
	t.prop01.res ~ t.K
'
fitk <- sem(modelK, data = mur_cor_all)
summary(fitk, fit.measures = T) # un buen modelo estadistico debe tener una CHisq no significativa! ojo!
fitMeasures(fitk, c("chisq", "df", "pvalue", "cfi", "rmsea"))
standardizedsolution(fitk)
modificationIndices(fitk, standardized = F) # mi es el aporte a CHisq de ser incluido el path sugerido
graph_sem(model = fitk)



#################           ANALISIS DEL PAPER                  ################ 

# Generales de la relacion de macro nutrientes y tasa de colonizacion 
# PLANTA obtiene fosforo y nitrogeno / hongo recibe carbono. Carbohidratos o lipidos. 
# Mayor nitrogeno mayor tasa de colinizacion pero mucho nitrogeno no necesita a AMF 
# Mayor fosforo menor tasa de colinizacion 

# Generales del efecto entre macro nutrientes afectando disponibilidad en suelo 
# N afecta negativamente disponibilidad de K
# P afecta negativamente disponibilidad de K
# No hay relacion entre N y P al menos directamente. 


# Primeros modelos introductorios 
# la logica causal que me gustaria es: suelo -> micorrizas -> planta 
# idealmente ubiera sido fitness en planta, pero creo que podemos hacer algo 
# con las variables que tenemos. 

######  PASO 1)  CONTRASTES ENTRE MODELOS ALTERNATIVOS ########## 

##### MODELO 1 #####
model1<-'				
	t.prop01 ~ t.N + t.P + t.K    
	t.K ~ t.N  # se sabe que N afecta neg. a K
  t.K ~ t.P  # se sabe que P afecta neg, a K
	t.N ~ t.pH # se sabe que pH afecta disp de cada macronutriente                  
	t.P ~ t.pH                           
  t.K ~ t.pH
'

fit1 <- sem(model1, data = mur_cor_all)
fitMeasures(fit1, c("chisq", "df", "pvalue", "cfi", "rmsea"))
modificationIndices(fit1, standardized = F, sort = T) # "mi" es el aporte a CHisq de ser incluido el path sugerido
                                                      # The EPC is the value by which a constrained relationship would change from  
                                                      # zero if it was freed to be estimated by the model.
summary(fit1, standardized = T)  #, fit.measures = T
standardizedsolution(fit1) #<---------------- REPORTADO EN TESIS Y PAPER
graph_sem(model = fit1) # para entender el grafico 
                        # ve https://cran.r-project.org/web/packages/tidySEM/vignettes/sem_graph.html


# Este modelo muestra como el efecto de K sobre la prop01 se debe a los 
# efectos indirectos de P y N sobre K. sIN EMBARGO NO TIENE UN BUEN AJUSTE


###### MODELO 2  ##################
model2<-'				
	t.prop01 ~  t.K    
	t.K ~ t.N  # se sabe que N afecta neg. a K
  t.K ~ t.P  # se sabe que P afecta neg, a K
	t.N ~ t.pH # se sabe que pH afecta disp de cada macronutriente                  
	t.P ~ t.pH                           
  t.K ~ t.pH
'

fit2 <- sem(model2, data = mur_cor_all)
fitMeasures(fit2, c("chisq", "df", "pvalue", "cfi", "rmsea"))
modificationIndices(fit2, standardized = F, sort = T) # mi es el aporte a CHisq de ser incluido el path sugerido
# The EPC is the value by which a constrained relationship would change from  
# zero if it was freed to be estimated by the model.
summary(fit2)  #, fit.measures = T
standardizedsolution(fit2)
graph_sem(model = fit2)


# Este modelo asume que todos los efectos importantes van via N y P 
# como se suele pensar que ocurre en AMF-plant 


# ESTE MODELO ES SIMILAR EN AJUSTE AL MODELO1 E INCLUSO SI SE CONSIDERA AIC ES MEJOR. 
# SIN EMBARGO, SE VUELVE MENOS INTERESANTE CUANDO SE HACE LA COMPARACION ENTRE AMBIENTES 
# YA QUE LA RELACION ENTRE T.PROP01 Y K ES INTERESANTE. 

####### MODELO 3   ############
model3<-'				
	t.prop01 ~ t.N + t.P    
	t.N ~ t.pH # se sabe que pH afecta disp de cada macronutriente                  
	t.P ~ t.pH                           
'
fit3 <- sem(model3, data = mur_cor_all)
fitMeasures(fit3, c("chisq", "df", "pvalue", "cfi", "rmsea"))
modificationIndices(fit3, standardized = F, sort = T) # mi es el aporte a CHisq de ser incluido el path sugerido
summary(fit3, standardized = T)  #, fit.measures = T
standardizedsolution(fit3)
graph_sem(model = fit3)



# en este caso no es bueno usar una anova para comparar ya que el modelo 
# fit2 no es estrictamente un modelo anidado en el modelo fit1. la estructura
# cambia mucho. No es una variante de fit1. Por esto lo mejor seria 
# utilizar los valores de RMSA y CFI para decidir cual modelo es mejor 
# valores de CFI cercanos a 1 indican mejor ajuste y rmsea bajos tambien. Por lo que 
# nos debemos de quedar con el primer modelo. 

fitMeasures(fit1, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr","AIC")) # primer modelo 
fitMeasures(fit2, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr","AIC")) # segundo modelo
fitMeasures(fit3, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr","AIC"))

anova(fit1, fit3)
lavTestLRT(fit1, fit3)


#############             MULTIGROUP TEST                        ############


###### MODELO 1 ganador y que REPORTAMOS EN EL PAPER   #######
model1<-'				
	t.prop01 ~ t.N + t.P + t.K    
	t.K ~ t.N  # se sabe que N afecta neg. a K
  t.K ~ t.P  # se sabe que P afecta neg, a K
	t.N ~ t.pH # se sabe que pH afecta disp de cada macronutriente                  
	t.P ~ t.pH                           
  t.K ~ t.pH
'

fit1_overall <- sem(model1,      
                    data = mur_cor_all, estimator ="MLM")
summary (fit1_overall)
fitMeasures(fit1_overall, c("chisq", "df", "pvalue", "cfi", "rmsea","srmr")) # <- reportado

######### Modelo general aplicable a ambos ambientes?   ######################

####### Evaluando que en ambos ambientes es aplicable el modelo 
# rural<-subset(mur_cor_all, mur_cor_all$ambient == "rural")
# urban<-subset(mur_cor_all, mur_cor_all$ambient == "urban")
# 
# fit1_rural <- sem(model1,
#                   data = rural)
# fit1_urban <- sem(model1,
#                  data = urban)
# summary (fit1_rural)
# summary (fit1_urban)
# fitMeasures(fit1_rural, c("chisq", "df", "pvalue", "cfi", "rmsea","srmr"))
# fitMeasures(fit1_urban, c("chisq", "df", "pvalue", "cfi", "rmsea","srmr"))

###################     iniciando multigroup test   ################################# 

########  FREE PARAMETERS IN BOTH ENVIRONMENTS ########## 

####### Model 1: configural invariance.  ########### 
# The same factor structure is imposed on all groups.
fit1_config_invariance <- sem(model1, group = "ambient",
                  data = mur_cor_all, estimator ="MLM")

summary(fit1_config_invariance, standardize = TRUE, rsq = T) 
graph_sem(model = fit1_config_invariance) 
modindices(fit1_config_invariance, sort = TRUE, minimum.value = 3) 
fitMeasures(fit1_config_invariance, c("chisq", "df", "pvalue", "cfi", "rmsea","srmr"))



#######  CONSTRAINING REGRESSIONS TO BE EQUAL BETWEEN ENVIRONMENTS ##########
####### Model 2: weak invariance. ###########
# The factor loadings are constrained to be equal across groups.
# Constraining regressions: como estamos interesados en las asociaciones...
# creamos un modelo donde todas las regresiones son fijadas a ser iguales entre
# ambientes. Es decir, construimos nuestro modelo nulo. 

fit1_weak_constrain <- sem(model1, group = "ambient",
                           data = mur_cor_all , estimator ="MLM",
                  group.equal = "regressions") # regresiones iguales osea, path coef 


### CONSTRAINING INTERCEPTS AND REGRESSIONS ################
####### Model 3: strong invariance. ########## 
#The factor loadings and intercepts are constrained to be equal across groups.
# # Constraining regressions and intercepts: este modelo NO lo usaremos, pero muestra
# # la posibilidad de probar que los interceptos (las medias) y las regresiones
# # son iguales en ambos ambientes. En el paper ya mostramos que hay diferencia 
# # de medias. 

fit1_strong_constrain <- sem(model1, group = "ambient",
                           data = mur_cor_all, estimator ="MLM",
                    group.equal = c("intercepts", "regressions"))

lavTestScore(fit1_strong_constrain, epc = TRUE, cumulative = TRUE)

########### Anova: free vs. constrained model   ##############

lavTestLRT(fit1_config_invariance, fit1_weak_constrain, fit1_strong_constrain)
anova (fit1_config_invariance, fit1_strong_constrain)  # <---- esta indica que difieren los modelos *
compareFit(fit1_config_invariance, fit1_weak_constrain)

# * Segun https://jslefche.github.io/sem_book/multigroup-analysis.html#introduction-to-multigroup-analysis
# Y este resultado justifica buscar que puede mejorar el modelo ... sin embargo, 
# https://lavaan.ugent.be/tutorial/groups.html (ver hasta el final) concluye directamente que 
# config invariance vs weak constrain probara si hay diferencias en las regresiones
# weak constrain vs strong constrain pribara si hay diferencias en las medias (interceptos)


################# LAGRANGE MULTIPLIER TEST  ############# 
###    Usando lagrange multiplier test para evaluar que parametro se deberia 
# liberar del modelo 

# EVALUA QUE RESTRICCIONES DEBEMOS QUITAR DEL MODELO RESTRINGIDO 
# originalmente pensado
lavTestScore(fit1_weak_constrain, epc = TRUE, cumulative = TRUE) # Yo andaba buscando esta y ya la encontre 
                                                                 # Lagrange multiplier test for releasing one 
                                                                 # or more fixed or constrained parameters and 
                                                                 # if we include epc = TRUE will return the modindices

#


# Sobre el output... 
# El univariate score test sugiere liberar la restriccion 
# p3 == p21: t.prop01  ~      t.K

fit1_weak_constrain.release.uni <- sem(model1, group = "ambient",
                           data = mur_cor_all, 
                           group.equal = "regressions",  estimator ="MLM",
                           group.partial = c("t.prop01~ t.K")) # liberando la restriccion 

summary(fit1_weak_constrain.release.uni, standardize = TRUE, rsq = T)
lavTestLRT(fit1_weak_constrain.release.uni, fit1_weak_constrain)   #  la liberacion de este path mejora el modelo sig 


# El cumulative score test sugiere liberar mas restricciones 
# p3 == p21: t.prop01  ~      t.K
# p4 == p22: t.K  ~      t.N 
# p6 == p24: t.N  ~     t.pH 


fit1_weak_constrain.release.cum <- sem(model1, group = "ambient",
                                    data = mur_cor_all, 
                                    group.equal = "regressions", estimator ="MLM",
                                    group.partial = c("t.prop01~ t.K",
                                                      "t.K ~ t.N", 
                                                      "t.N ~ t.pH")) 

lavTestLRT(fit1_weak_constrain.release.cum, fit1_weak_constrain)  

summary (fit1_weak_constrain.release.cum, rsq = T)
standardizedSolution(fit1_weak_constrain.release.cum) # no creo que sea buena idea usarlo.  Dan diferencias entre paths que se fijaron 
parameterEstimates(fit1_weak_constrain.release.cum, standardized = TRUE)
graph_sem(model = fit1_weak_constrain.release.cum, ellipses_width = 3 )
