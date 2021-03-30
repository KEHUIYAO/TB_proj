library(readr)
library(lubridate)
library(dplyr)
library(ggplot2)
data <- read_csv("../data/MI_TBData_FOIProj092519.csv")
data$Tag_Date = as.Date(data$Tag_Date)
outlier_date = sort(data$Tag_Date)[1:3]
data = data[!(data$Tag_Date %in% outlier_date),]

# prevalence rate
nrow(data)
sum(data$`TB+` == "Pos")
sum(data$`TB+` == "Pos") / nrow(data)
sum(data$`TB+` == "Pos" & data$Sex == "female") / sum(data$Sex == "female")
sum(data$`TB+` == "Pos" & data$Sex == "male") / sum(data$Sex == "male")
data$Tag_Date = year(data$Tag_Date)
colnames(data)[8] = "TB"
# Number of cases in each year and general disease trend across 30 --------


temp_data = data


ggplot(temp_data,aes(Tag_Date,fill = Sex))+geom_bar()




# incidence rate
# temp_data = 
#   temp_data %>% 
#   mutate(Tag_Date = as.factor(as.character(Tag_Date))) %>%
#   filter(!is.na(Tag_Date)) %>%
#   filter(Sex %in% c("male","female")) %>%
#   group_by(Tag_Date,Sex) %>% 
#   summarise(positive_case = sum(TB == "Pos"), positive_rate = sum(TB == "Pos") / n()) 
# 
# temp_data$Tag_Date = as.numeric(as.character(temp_data$Tag_Date))
# ggplot(temp_data,aes(x=Tag_Date,y=positive_rate,color = Sex)) + geom_line()



# prevalence rate across year

temp_data = data


temp_data = 
  temp_data %>% 
  mutate(Tag_Date = as.factor(as.character(Tag_Date))) %>%
  filter(!is.na(Tag_Date)) %>%
  filter(Sex %in% c("male","female")) %>%
  group_by(Tag_Date,Sex) %>% 
  summarise(positive_case = sum(TB == "Pos"), num_of_cases = n()) 

temp_data
temp_data$Tag_Date = as.numeric(as.character(temp_data$Tag_Date))
female = temp_data%>%filter(Sex == "female")
male = temp_data%>%filter(Sex == "male")
x1 = cumsum(female$positive_case)/cumsum(female$num_of_cases)
x2 = cumsum(male$positive_case)/cumsum(male$num_of_cases)
x3 = c(rep("female",length(x1)),rep("male",length(x1)))
new_data = data.frame(Tag_Date =c(unique(temp_data$Tag_Date),unique(temp_data$Tag_Date)),prevalence = c(x1,x2), sex = x3)
ggplot(new_data,aes(x=Tag_Date,y=prevalence,color = sex)) + geom_line()+xlab("Tag Date")+ylab("Prevalence rate")


# Age structure and disease pattern across ages-----------------------------------------------------------
temp_data = data%>%filter(Sex %in% c("male","female"))
newdata = temp_data%>%group_by(Age,Sex)%>%summarise(count = n())
ggplot(newdata,aes(x=reorder(Age,-count),y=count,fill=Sex))+geom_bar(stat="identity",position="dodge")+xlab("Age")+ylab("Count")


# 
# temp_data = temp_data %>%
#   group_by(Age,Sex) %>%
#   summarise(positive_case = sum(TB == "Pos"), num_of_case =n())
# temp_data$Age = as.numeric(temp_data$Age)
# female = temp_data%>%filter(Sex == "female")
# male = temp_data%>%filter(Sex == "male")
# female$prevalence = female$positive_case /sum(data$Sex == "female")
# male$prevalence = male$positive_case / sum(data$Sex == "male")
# newdata = rbind(female,male)
# ggplot(newdata,aes(x=Age,y=prevalence,fill=Sex)) + geom_bar(stat="identity",position = "dodge")+xlab("Age")+ylab("Age prevalence")

# positive rate among each age class
temp_data = data%>%filter(Sex %in% c("male","female"))
temp_data = temp_data %>%
  group_by(Age,Sex) %>%
  summarise(positive_case = sum(TB == "Pos"), num_of_case =n())
temp_data = temp_data[temp_data$Age != "Q",]
temp_data$Age = as.numeric(temp_data$Age)
temp_data = temp_data %>% arrange(Age)
temp_data$age = as.character(temp_data$Age)
female = temp_data%>%filter(Sex == "female")
male = temp_data%>%filter(Sex == "male")
male = as.data.frame(male)
female = as.data.frame(female)
female$positive_rate = female$positive_case /female$num_of_case
male$positive_rate = male$positive_case / male$num_of_case

temp = data.frame(Age = 5, Sex = 'male',positive_case = 0,num_of_case = 0, age = 5, positive_rate = 0)
male = rbind(male,temp)

female$positive_rate[female$positive_rate == 0] = 10e-5

male$positive_rate[male$positive_rate == 0] = 10e-5



newdata = rbind(female,male)
#ind = seq(1,nrow(female))
ggplot(newdata,aes(x=reorder(age,Age),y=positive_rate,fill=Sex)) + geom_bar(stat="identity",position =position_dodge(width =1))+geom_text(aes(label = round(positive_rate,2)),position = position_dodge(width = 1),vjust = -0.5,size = 3)+xlab("Age")+ylab("Age positive rate")



# spatial
load("data_for_nimble_time.RData")
head(data_cleaned)
library(readr)
library("ggmap")
library("ggplot2")
library("proj4")

#register_google(key="AIzaSyBx4Key-P9lg0FUE5S0SyZwFSeQYQz0Szg")
#map <- get_map("michigan")
#ggmap(map)
tbdata_adj <- read_csv("../data/tbdata_adj.csv")
location_data <- cbind(tbdata_adj$X, tbdata_adj$Y)
str(data_cleaned)
space_pos = data_cleaned %>% arrange(Adj_ID) %>% 
  group_by(Adj_ID)%>% 
  summarise(rate = sum(status) / n())

space_positive_rate = rep(0,N)
for (i in 1:nrow(space_pos)){
  space_positive_rate[as.numeric(space_pos[i,1])] = space_positive_rate[as.numeric(space_pos[i,1])] + as.numeric(space_pos[i,2])
}

plot_sptial_heatmap(space_positive_rate,location_data[,1],location_data[,2],10,10)
#save(space_positive_rate,location_data,file = "eda_for_map.RData")
