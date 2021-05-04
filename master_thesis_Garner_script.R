library(dplyr)
library(tidyr)
library(lubridate)
library(DHARMa)
library(MuMIn)
library(glmmTMB) 


data <- read.csv("master_thesis_Garner_raw_data.csv",
                 header=TRUE,sep=";") %>%
  separate(colony_ID,c("region","management",NA),sep="_",remove = F) # add variables

#### FLIGHT DISTANCES ####

# new variable: time between start of observation and dances for nest sites
day1 <- make_difftime(86400) # 1 day in s
data1 <- data  %>% 
  mutate(observation_start_date=dmy(observation_start_date),
         observation_start_daytime=hm(observation_start_daytime),
         dance_date=dmy(dance_date),
         dance_daytime=hm(dance_daytime)) %>% 
  unite(observation_start,observation_start_date,observation_start_daytime,sep = " ",remove = F) %>% 
  unite(dance_time,dance_date,dance_daytime,sep = " ",remove = F) %>% 
  mutate(observation_start=ymd_hms(observation_start),
         dance_time=ymd_hms(dance_time),
         dance_time_diff=difftime(dance_time,observation_start,units = "secs"),
         dance_time_diff=if_else(day(dance_time)==day(observation_start)+1,difftime(dance_time-day1/2,observation_start,units = "secs"),dance_time_diff),
         dance_time_diff=if_else(day(dance_time)==day(observation_start)+2,difftime(dance_time-day1,observation_start,units = "secs"),dance_time_diff),
         dance_time_diff=if_else(day(dance_time)==day(observation_start)+3,difftime(dance_time-day1/2*3,observation_start,units = "secs"),dance_time_diff),
         dance_time_diff=as.numeric(dance_time_diff)/60)
         

# compare models via AIC
data1$flight_distance0 <- round(data1$flight_distance) 
data1$flight_distance1 <- round(if_else(data1$flight_distance<1,1,data1$flight_distance)) # flight_distance >= 1

m1 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|region),
              data=data1,family="nbinom2")
m2 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|region/observation_start_date),
              data=data1,family="nbinom2")
m3 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|region/colony_ID),
              data=data1,family="nbinom2")
m4 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|region/observation_start_date/colony_ID),
              data=data1,family="nbinom2")
m5 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|observation_start_date/colony_ID),
              data=data1,family="nbinom2")
m6 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|observation_start_date),
              data=data1,family="nbinom2")
m7 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|colony_ID),
              data=data1,family="nbinom2")

AICc.tab <- MuMIn::AICc(m1,m2,m3,m4,m5,m6,m7) # m5


m8 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|observation_start_date/colony_ID), # Poisson
              data=data1)
m9 <- glmmTMB(flight_distance0~management*dance_time_diff+(1|observation_start_date/colony_ID), # zero inflated
              data=data1,family="nbinom2",ziformula=~1)
m10 <- glmmTMB(flight_distance1~management*dance_time_diff+(1|observation_start_date)+(1|colony_ID), # dist >= 1
               data=data1,family="nbinom2")

AICc.tab <- MuMIn::AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10) # m9

simulateResiduals(fittedModel=m9,plot=T)
summary(m9)


## data frame with predicted flight distances (m9)
data_pred <- expand.grid(dance_time_diff=data1$dance_time_diff,
                         management=c("NP","WW"))
## design matrix (fixed effects)
mm <- model.matrix(delete.response(terms(m9)),data_pred)
data_pred$predicted_distance <- drop(mm %*% fixef(m9)[["cond"]])

predvar <- diag(mm %*% vcov(m9)[["cond"]] %*% t(mm))
data_pred$SE <- sqrt(predvar) 
data_pred$SE2 <- sqrt(predvar+sigma(m9)^2)

data_pred <- data_pred %>%
  mutate(y=exp(predicted_distance),
         ymin=exp(predicted_distance-1.96*SE), #..CI of the regression line
         ymax=exp(predicted_distance+1.96*SE))



#### NEST SITE LOCATIONS ####
library(rgdal)
library(raster)

# example Kellerwald
# add 75% uncertainty buffer to nest sites
data_K <- data %>% 
  filter(region=="K") %>%
  mutate(nest_site_.coordinates_E=as.numeric(nest_site_.coordinates_E),
         nest_site_.coordinates_N=as.numeric(nest_site_.coordinates_N),
         nest_site_distance=as.numeric(flight_distance)) %>%
  mutate(uncert75=33.7335867+0.3301049*flight_distance) 

# (-> QGIS -> save as shape)

nestsites <- shapefile("data_K.shp")
swarmstands <- shapefile("swarmstands_K.shp") 
forest_types <- shapefile("forest_types_K.shp") # forestry data from forestry administrations and CORINE

# buffer + intersect nest sites
nestsites_buffer <- buffer(nestsites, width=nestsites$uncert75, dissolve = FALSE)
nestsites_inters <- raster::intersect(nestsites_buffer,forest_types)
nestsites_inters$area <- area(nestsites_inters)

# define area for analysis (area with 90% of the dances / r <= 4 km)
buffer <- data_K %>% 
  group_by(colony_ID) %>% 
  summarise(Q90=quantile(flight_distance,c(.90))) %>% 
  mutate(R4=4000,
         Q90.75=Q90+(33.7335867+0.3301049*Q90),
         R4.75=4000+(33.7335867+0.3301049*4000),
         buffer = if_else(R4<Q90,R4,Q90),
         buffer.75 = if_else(R4<Q90,R4.75,Q90.75)
  )

# buffer + intersect areas around swarmstands
swarmstands@data <- left_join(swarmstands@data,buffer,by="colony_ID")
swarmstands_buffer <- buffer(swarmstands, width=swarmstands$buffer.75, dissolve = FALSE)
swarmstands_inters <- raster::intersect(swarmstands_buffer,forest_types)
swarmstands_inters$area <- area(swarmstands_inters)

# shapes to data frames
nestsites_DF <- as.data.frame(nestsites_inters)
swarmstands_DF <- as.data.frame(swarmstands_inters)


#### proportions of forest types in dance buffers 
#### -> observed proportions
nestsites.p <- nestsites_DF %>% 
  mutate(prop=area/(pi*(uncert75^2))) %>% 
  group_by(dance_ID,forest_type) %>% 
  summarise(prop=sum(prop))

# add nest sites in unforested areas
nestsites_spread <- full_join(nestsites.p,data_K %>% dplyr::select(dance_ID)) %>% 
  mutate(forest_type=if_else(is.na(forest_type),"unforested",forest_type)) %>% 
  spread(key=forest_type,value=prop,fill=0)

# add missing unforested areas in buffers
nestsites_spread$unforested <- 1-rowSums(nestsites_spread[,-1])

# transform back
nestsites_gather <- nestsites_spread %>% 
  gather(forest_type,"prop",-1)

# filter nestsites in area of analysis
nestsites_new <- full_join(data_K,nestsites_gather,by="dance_ID") 
nestsites_new <- left_join(nestsites_new,buffer,by="colony_ID") %>% 
  filter(flight_distance<=buffer)


#### proportions of forest types around the swarmstands 
#### -> expected proportions

swarmstands.p <- swarmstands_DF %>% 
  mutate(prop=area/(pi*(buffer.75^2))) %>% 
  group_by(colony_ID,forest_types) %>% 
  summarise(prop=sum(prop))

# add unforested areas in buffers
swarmstands_spread <- swarmstands.p %>%
  spread(key=forest_types,value=prop,fill=0) 
swarmstands_spread$unforested <- 1-rowSums(swarmstands_spread[,-1])

# transform back
swarmstands_new <- swarmstands_spread%>% 
  gather(forest_types,"prop",-1) 


#### Chi-square Goodness of fit test
# differ expected forest types from observed?

n.nestsites<- nestsites_new %>% # number of dances per colony
  group_by(colony_ID,dance_ID) %>% 
  summarise(n=n()) %>% 
  group_by(colony_ID) %>% 
  summarise(n=n())
nestsites.mean <- nestsites_new %>% # forest type proportion per colony = prop
  group_by(colony_ID,forest_types) %>% 
  summarise(prop=mean(prop)) 

obs.colony <- left_join(nestsites.mean,n.nestsites,by="colony_ID") %>% 
  mutate(number.of.dances = prop*n) %>% # observed = number.of.dances
  mutate(area="observed")
exp.colony <- left_join(swarmstands_new,n.nestsites,by="colony_ID") %>% 
  mutate(number.of.dances = prop*n)%>% # expected = prop
  mutate(area="expected")


# test per colony
chi <- left_join(obs.colony %>% dplyr::select(colony_ID,forest_types,number.of.dances),
                 exp.colony %>% dplyr::select(colony_ID,forest_types,prop),
                 by=c("colony_ID","forest_types")) %>% 
  filter(prop!="0")

chisq.test <- data.frame()

val <- chi %>% 
  dplyr::select(colony_ID) %>% 
  distinct(colony_ID)
col <- val$colony_ID

for (colony in col) {
  test <- chisq.test(chi$number.of.dances[chi$colony_ID==colony],p=chi$prop[chi$colony_ID==colony])
  n <- sum(chi$number.of.dances[chi$colony_ID==colony])
  K <- colony
  X <- test$statistic
  df <- test$parameter
  p <- test$p.value
  chix <- data.frame(K,df,n,X,p)
  chisq.test <- rbind.data.frame(chisq.test,chix)
}


# test for all dances
chi2 <- left_join(obs.colony %>% dplyr::select(colony_ID,forest_types,number.of.dances),
                  exp.colony %>% dplyr::select(colony_ID,forest_types,number.of.dances),
                  by=c("colony_ID","forest_types"),suffix=c(".obs",".exp")) %>% 
  group_by(forest_types) %>% 
  summarise(number.of.dances.obs=sum(number.of.dances.obs), # observed n
            number.of.dances.exp=sum(number.of.dances.exp)) %>% 
  mutate(prop.obs=number.of.dances.obs/sum(number.of.dances.obs),
         prop.exp=number.of.dances.exp/sum(number.of.dances.exp)) %>% # expected prop
  filter(prop.exp!="0")
n <- sum(chi2$number.of.dances.exp)


test <- chisq.test(chi2$number.of.dances.obs,p=chi2$prop.exp)
X <- test$statistic
df <- test$parameter
p <- test$p.value
chisq.test2 <- data.frame(df,n,X,p)