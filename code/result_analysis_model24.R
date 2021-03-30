load("model24.RData")
load("data_for_nimble.RData")
library("nimble")
source("utility.R")
summary_statistics = samples$summary$all.chains
#summary_statistics = (samples$summary$chain1 + samples$summary$chain2) / 2
#summary_statistics = samples$summary$chain1
#summary_statistics = samples$summary$chain2
row_names = rownames(summary_statistics)
summary_statistics = summary_statistics[,1]
grand_mean_ind = grep("mu", row_names)
grand_mean_effect = summary_statistics[grand_mean_ind]
fbeta0_ind = grep("fbeta0", row_names)
mbeta0_ind = grep("mbeta0",row_names)
sex_effect = c(summary_statistics[mbeta0_ind], summary_statistics[fbeta0_ind])
fhetero_ind = grep("fhetero", row_names)
fhetero_effect = summary_statistics[fhetero_ind]
mhetero_ind = grep("mhetero", row_names)
mhetero_effect = summary_statistics[mhetero_ind]
fcar_age_ind = grep("fcar_age", row_names)
fcar_age_effect = summary_statistics[fcar_age_ind]
mcar_age_ind = grep("mcar_age", row_names)
mcar_age_effect = summary_statistics[mcar_age_ind]
fcar_time_ind = grep("fcar_time", row_names)
fcar_time_effect = summary_statistics[fcar_time_ind]
mcar_time_ind = grep("mcar_time", row_names)
mcar_time_effect = summary_statistics[mcar_time_ind]




# recover space effect and save all the results----------------------------------------------------

raw_data = (samples$samples$chain1 + samples$samples$chain2 + samples$samples$chain3) / 3
#raw_data = (samples$samples$chain1 + samples$samples$chain2) / 2
u1_ind = grep("u\\[1", row_names)
u2_ind = grep("u\\[2", row_names)
omega_ind = grep("omega", row_names)
mvspace = array(NA,dim = c(2,length(u1_ind),nrow(raw_data)))
for (i in 1:nrow(raw_data)){
  temp_u = rbind(raw_data[i,u1_ind],raw_data[i,u2_ind])
  temp_omega = raw_data[i,omega_ind]
  temp_omega = matrix(temp_omega,nrow = 2)
  temp_Cov <- inverse(temp_omega)
  temp_achol <- t(chol(temp_Cov))
  # Based on matrix algebra, the mvspace will also have sum to zero constraint
  mvspace[,,i] = temp_achol%*%temp_u
}

cur_sum = matrix(0,nrow = 2, ncol = length(u1_ind))
for (i in 1:nrow(raw_data)){
  cur_sum = cur_sum + mvspace[,,i]
}
mvspace_effect = cur_sum / nrow(raw_data)

mvspace_sd = matrix(0,nrow = 2, ncol = length(u1_ind))
for(i in 1:2){
  for (j in 1:length(u1_ind)){
    mvspace_sd[i,j] = sd(mvspace[i,j,])
  }
}



# result analysis: check convergence ---------------------------------------------------------

chain1 = samples$samples$chain1
chain2 = samples$samples$chain2
chain3 = samples$samples$chain3

mu_1 = chain1[,'mu']
mu_2 = chain2[,'mu']
mu_3 = chain3[,'mu']

fbeta0_1 = chain1[,"fbeta0"]
fbeta0_2 = chain2[,"fbeta0"]
fbeta0_3 = chain3[,"fbeta0"]

mbeta0_1 = chain1[,"mbeta0"]
mbeta0_2 = chain2[,"mbeta0"]
mbeta0_3 = chain3[,"mbeta0"]

fcar_age_1 = chain1[,"fcar_age[2]"]
fcar_age_2 = chain2[,"fcar_age[2]"]
fcar_age_3 = chain3[,"fcar_age[2]"]

plot_chains(mu_1,mu_2,mu_3)
plot_chains(fbeta0_1,fbeta0_2,fbeta0_3)

plot_chains(fcar_age_1,fcar_age_2,fcar_age_3)




# result analysis: heat map -----------------------------------------------
library(readr)
library("ggmap")
library("ggplot2")
library("proj4")

#register_google(key="AIzaSyBx4Key-P9lg0FUE5S0SyZwFSeQYQz0Szg")
#map <- get_map("michigan")
#ggmap(map)-
tbdata_adj <- read_csv("~/Nutstore Files/study/research assitant project/tb project/data/tbdata_adj.csv")
location_data <- cbind(tbdata_adj$X, tbdata_adj$Y)

mspace_effect = mvspace_effect[1,]
fspace_effect = mvspace_effect[2,]
#save(mspace_effect,fspace_effect,mvspace_sd,location_data,file = "model24_for_map.RData")

plot_sptial_heatmap(mspace_effect,location_data[,1],location_data[,2],10,10)
plot_sptial_heatmap(fspace_effect,location_data[,1],location_data[,2],10,10)
plot_sptial_heatmap((fspace_effect+mspace_effect)/2,location_data[,1],location_data[,2],10,10)


# result analysis: age effect contribution to the hazard function ---------------------------------------------
# result analysis: age effect contribution to the hazard function ---------------------------------------------
summary_statistics = samples$summary$all.chains
confidence_interval = summary_statistics[,c(4,5)]
fcar_age_lower = confidence_interval[fcar_age_ind,1]
fcar_age_upper = confidence_interval[fcar_age_ind,2]
mcar_age_lower = confidence_interval[mcar_age_ind,1]
mcar_age_upper = confidence_interval[mcar_age_ind,2]
std = summary_statistics[,3]
f_std = std[fcar_age_ind]
m_std = std[mcar_age_ind]
limit = 12
fcar_age_lower = fcar_age_effect
fcar_age_lower[1:limit] = fcar_age_lower[1:limit] - f_std[1:limit]
fcar_age_upper = fcar_age_effect
fcar_age_upper[1:limit] = fcar_age_upper[1:limit] + f_std[1:limit]

mcar_age_lower = mcar_age_effect 
mcar_age_lower[1:limit] = mcar_age_lower[1:limit] - m_std[1:limit]
mcar_age_upper = mcar_age_effect
mcar_age_upper[1:limit] = mcar_age_upper[1:limit] + m_std[1:limit]



f_age_hazard = rep(NA,Tage)
m_age_hazard = rep(NA,Tage)
f_age_lower_hazard =rep(NA,Tage)
f_age_upper_hazard = rep(NA,Tage)
m_age_lower_hazard = rep(NA,Tage)
m_age_upper_hazard = rep(NA,Tage)
for (k in 1:Tage) {
  f_age_hazard[k] <- exp(grand_mean_effect + (sex_effect[2]+fcar_age_effect[k])+mean(mvspace_effect[2,])+mean(fhetero_effect)+mean(fcar_time_effect))
  m_age_hazard[k] <- exp(grand_mean_effect + (sex_effect[1]+mcar_age_effect[k])+mean(mvspace_effect[1,])+mean(mhetero_effect)+mean(mcar_time_effect))
  f_age_lower_hazard[k] <- exp(grand_mean_effect + (sex_effect[2]+fcar_age_lower[k])+mean(mvspace_effect[2,])+mean(fhetero_effect)+mean(fcar_time_effect))
  m_age_lower_hazard[k] <- exp(grand_mean_effect + (sex_effect[1]+mcar_age_lower[k])+mean(mvspace_effect[1,])+mean(mhetero_effect)+mean(mcar_time_effect))
  f_age_upper_hazard[k] <-  exp(grand_mean_effect + (sex_effect[2]+fcar_age_upper[k])+mean(mvspace_effect[2,])+mean(fhetero_effect)+mean(fcar_time_effect))
  m_age_upper_hazard[k] <-  exp(grand_mean_effect + (sex_effect[1]+mcar_age_upper[k])+mean(mvspace_effect[1,])+mean(mhetero_effect)+mean(mcar_time_effect))
  
} 

#step_hazard_plot(seq(0.5,Tage/2,0.5),c("female","male"),f_age_hazard,m_age_hazard)
truncated = 12

#age_hazard_plot = step_hazard_plot(seq(1,Tage),c("female","male"),f_age_hazard,m_age_hazard)
#age_hazard_plot+geom_errorbar(aes(x= x-0.5,ymin = c(f_age_lower_hazard,m_age_lower_hazard),ymax=c(f_age_upper_hazard,m_age_upper_hazard)),width = 0)
age_hazard_plot = step_hazard_plot(seq(0.5,truncated/2,0.5),c("female","male"),f_age_hazard[1:truncated],m_age_hazard[1:truncated])
age_hazard_plot+geom_errorbar(aes(x= x-0.25,ymin = c(f_age_lower_hazard[1:truncated],m_age_lower_hazard[1:truncated]),ymax=c(f_age_upper_hazard[1:truncated],m_age_upper_hazard[1:truncated])),width = 0)



# result analysis: time effect contribution to the hazard function --------
f_std = std[fcar_time_ind]
m_std = std[mcar_time_ind]
limit = Ttime
ftime_lower = fcar_time_effect
ftime_lower[1:limit] = ftime_lower[1:limit] - f_std[1:limit]
ftime_upper = fcar_time_effect
ftime_upper[1:limit] = ftime_upper[1:limit] + f_std[1:limit]

mtime_lower = mcar_time_effect
mtime_lower[1:limit] = mtime_lower[1:limit] - m_std[1:limit]
mtime_upper = fcar_time_effect
mtime_upper[1:limit] = mtime_upper[1:limit] + m_std[1:limit]
f_time_hazard = rep(NA,Ttime)
m_time_hazard = rep(NA,Ttime)
f_time_lower_hazard= rep(NA,Ttime)
f_time_upper_hazard= rep(NA,Ttime)
m_time_lower_hazard= rep(NA,Ttime)
m_time_upper_hazard= rep(NA,Ttime)
for (k in 1:Ttime){
  f_time_hazard[k] = exp(grand_mean_effect +mean(fcar_age_effect) + sex_effect[2] + fcar_time_effect[k]+mean(mvspace_effect[2,])+mean(fhetero_effect))
  m_time_hazard[k] = exp(grand_mean_effect +mean(mcar_age_effect) + sex_effect[1] + fcar_time_effect[k]+mean(mvspace_effect[2,])+mean(mhetero_effect))
  f_time_lower_hazard[k] = exp(grand_mean_effect +mean(fcar_age_effect) + sex_effect[2] + ftime_lower[k]+mean(mvspace_effect[2,])+mean(fhetero_effect))
  f_time_upper_hazard[k] = exp(grand_mean_effect +mean(fcar_age_effect) + sex_effect[2] + ftime_upper[k]+mean(mvspace_effect[2,])+mean(fhetero_effect))
  m_time_lower_hazard[k] = exp(grand_mean_effect +mean(mcar_age_effect) + sex_effect[1] + mtime_lower[k]+mean(mvspace_effect[1,])+mean(mhetero_effect))
  m_time_upper_hazard[k] = exp(grand_mean_effect +mean(mcar_age_effect) + sex_effect[1] + mtime_upper[k]+mean(mvspace_effect[1,])+mean(mhetero_effect))
}

start = 
  
  step_hazard_plot(seq(1,Ttime),c("female","male"),f_time_hazard,m_time_hazard)

library("lubridate")
start_time = as.Date("1996-05-27")
temp =seq(1,Ttime)
new_temp = rep(start_time,length(temp))
for (j in 1:length(temp)){
  for (i in 1:6){
   new_temp[j] = new_temp[j]  %m+% months(temp[j])
  }
}


a = step_hazard_plot(new_temp,c("female","male"),f_time_hazard,m_time_hazard)+scale_x_date(date_labels = "%Y",date_breaks = "5 year",limit=c(as.Date("1996-05-27"),as.Date("2018-05-27")))
a
a + geom_errorbar(aes(x= x-0.5,ymin = c(f_time_lower_hazard,m_time_lower_hazard),ymax=c(f_time_upper_hazard,m_time_upper_hazard)),width = 0)

error_temp = new_temp
for (i in 1:6){
  error_temp[i] = error_temp[i] %m+% months(3)
}
b = step_hazard_plot(new_temp,c("male"),m_time_hazard)+scale_x_date(date_labels = "%Y",date_breaks = "5 year",limit=c(as.Date("1996-05-27"),as.Date("2018-05-27")))

b + geom_errorbar(aes(x=error_temp,ymin = c(m_time_lower_hazard),ymax=c(m_time_upper_hazard)),width = 0)


c = step_hazard_plot(new_temp,c("female"),f_time_hazard)+scale_x_date(date_labels = "%Y",date_breaks = "5 year",limit=c(as.Date("1996-05-27"),as.Date("2018-05-27")))
c + geom_errorbar(aes(x=error_temp,ymin = c(f_time_lower_hazard),ymax=c(f_time_upper_hazard)),width = 0)

# result analysis: 1.5year old time hazard ---------------------
f_std = std[fcar_time_ind]
m_std = std[mcar_time_ind]
limit = Ttime
ftime_lower = fcar_time_effect
ftime_lower[1:limit] = ftime_lower[1:limit] - f_std[1:limit]
ftime_upper = fcar_time_effect
ftime_upper[1:limit] = ftime_upper[1:limit] + f_std[1:limit]

mtime_lower = mcar_time_effect
mtime_lower[1:limit] = mtime_lower[1:limit] - m_std[1:limit]
mtime_upper = fcar_time_effect
mtime_upper[1:limit] = mtime_upper[1:limit] + m_std[1:limit]
f_time_hazard = rep(NA,Ttime)
m_time_hazard = rep(NA,Ttime)
f_time_lower_hazard= rep(NA,Ttime)
f_time_upper_hazard= rep(NA,Ttime)
m_time_lower_hazard= rep(NA,Ttime)
m_time_upper_hazard= rep(NA,Ttime)
for (k in 1:Ttime){
  f_time_hazard[k] = exp(grand_mean_effect +fcar_age_effect[3] + sex_effect[2] + fcar_time_effect[k]+mean(mvspace_effect[2,])+mean(fhetero_effect))
  m_time_hazard[k] = exp(grand_mean_effect +mcar_age_effect[3] + sex_effect[1] + fcar_time_effect[k]+mean(mvspace_effect[2,])+mean(mhetero_effect))
  f_time_lower_hazard[k] = exp(grand_mean_effect +fcar_age_effect[3] + sex_effect[2] + ftime_lower[k]+mean(mvspace_effect[2,])+mean(fhetero_effect))
  f_time_upper_hazard[k] = exp(grand_mean_effect +fcar_age_effect[3] + sex_effect[2] + ftime_upper[k]+mean(mvspace_effect[2,])+mean(fhetero_effect))
  m_time_lower_hazard[k] = exp(grand_mean_effect +mcar_age_effect[3] + sex_effect[1] + mtime_lower[k]+mean(mvspace_effect[1,])+mean(mhetero_effect))
  m_time_upper_hazard[k] = exp(grand_mean_effect +mcar_age_effect[3] + sex_effect[1] + mtime_upper[k]+mean(mvspace_effect[1,])+mean(mhetero_effect))
}

start = 
  
  step_hazard_plot(seq(1,Ttime),c("female","male"),f_time_hazard,m_time_hazard)

library("lubridate")
start_time = as.Date("1996-05-27")
temp =seq(1,Ttime)
new_temp = rep(start_time,length(temp))
for (j in 1:length(temp)){
  for (i in 1:6){
    new_temp[j] = new_temp[j]  %m+% months(temp[j]-1)
  }
}


a = step_hazard_plot(new_temp,c("female","male"),f_time_hazard,m_time_hazard)+scale_x_date(date_labels = "%Y",date_breaks = "5 year",limit=c(as.Date("1996-05-27"),as.Date("2018-11-27")))
a
a + geom_errorbar(aes(x= x-0.5,ymin = c(f_time_lower_hazard,m_time_lower_hazard),ymax=c(f_time_upper_hazard,m_time_upper_hazard)),width = 0)

error_temp = new_temp
for (i in 1:6){
  error_temp[i] = error_temp[i] %m+% months(3)
}
b = step_hazard_plot(new_temp,c("male"),m_time_hazard)+scale_x_date(date_labels = "%Y",date_breaks = "5 year",limit=c(as.Date("1996-05-27"),as.Date("2018-05-27")))

b + geom_errorbar(aes(x=error_temp,ymin = c(m_time_lower_hazard),ymax=c(m_time_upper_hazard)),width = 0)


c = step_hazard_plot(new_temp,c("female"),f_time_hazard)+scale_x_date(date_labels = "%Y",date_breaks = "5 year",limit=c(as.Date("1996-05-27"),as.Date("2018-05-27")))
c + geom_errorbar(aes(x=error_temp,ymin = c(f_time_lower_hazard),ymax=c(f_time_upper_hazard)),width = 0)

combined_hazard = (m_time_hazard+f_time_hazard)/2
combined_lower = (m_time_lower_hazard + f_time_lower_hazard)/2
combined_upper = (m_time_upper_hazard + f_time_upper_hazard)/2

d = step_hazard_plot(new_temp,c("total population"),combined_hazard)+scale_x_date(date_labels = "%Y",date_breaks = "5 year",limit=c(as.Date("1996-05-27"),as.Date("2018-05-27")))
d + geom_errorbar(aes(x=error_temp,ymin = c(combined_lower),ymax=c(combined_upper)),width = 0)

# 1.5 year incidence rate, 1000 deer
convert_hazard_to_incidence <- function(haz,num){
  incidence = (1-exp(-haz*0.5))*num
  annual_incidence = c()
  for (i in seq(1,length(incidence),2)){
    annual_incidence = c(annual_incidence,incidence[i]+incidence[i+1])
  }
  return (annual_incidence)
}


f_time_incidence = convert_hazard_to_incidence(f_time_hazard,1000)

m_time_incidence = convert_hazard_to_incidence(m_time_hazard,1000)

f_time_lower_incidence = convert_hazard_to_incidence(f_time_lower_hazard,1000)

f_time_upper_incidence = convert_hazard_to_incidence(f_time_upper_hazard,1000)

m_time_lower_incidence = convert_hazard_to_incidence(m_time_lower_hazard,1000)

m_time_upper_incidence = convert_hazard_to_incidence(m_time_upper_hazard,1000)



start_time = as.Date("1996-05-27")

end_time = as.Date("2018-05-27")
new_temp2 = seq(start_time,end_time,"years")

# e = ggplot(data.frame(time=rep(new_temp2,2),incidence_per_1000=c(f_annual_time_incidence,m_annual_time_incidence),sex = c(rep("female",length(new_temp2)),rep("male",length(new_temp2)))),
#            aes(x = time,y=incidence_per_1000,fill = sex)) +
#   geom_bar(stat = "identity",color = "black",position = position_dodge())+
#   geom_errorbar(aes(ymin = c(f_time_lower_incidence,m_time_lower_incidence),ymax = c(f_time_upper_incidence,m_time_upper_incidence) ),width=.2,
#                 position=position_dodge(.9))

e = ggplot( data.frame(time=rep(new_temp2,2),incidence_per_1000=c(f_annual_time_incidence,m_annual_time_incidence),sex = c(rep("female",length(new_temp2)),rep("male",length(new_temp2)))),aes(x=time,y=incidence_per_1000,color=sex))


e + geom_point() +  geom_errorbar(aes(ymin = c(f_time_lower_incidence,m_time_lower_incidence),ymax = c(f_time_upper_incidence,m_time_upper_incidence) ),
                                              position=position_dodge(.9))


# result analysis: space sex correlation ----------------------------------
omega_ind = grep("omega", row_names)
cor12 = rep(NA, nrow(raw_data))
for (i in 1:nrow(raw_data)){
  temp_omega = raw_data[i,omega_ind]
  temp_omega = matrix(temp_omega,nrow = 2)
  temp_Cov <- inverse(temp_omega)
  cor12[i] <- temp_Cov[1,2]/(sqrt(temp_Cov[1,1])*sqrt(temp_Cov[2,2]))
}

# the result suggests high correlation......
median(cor12)

