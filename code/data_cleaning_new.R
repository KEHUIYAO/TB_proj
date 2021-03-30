library("dplyr")
library("lubridate")

# load utility functions
source("utility.R")

# load the data
col_names = colnames(read.csv("../data/MI_TBData_FOIProj092519.csv"))
data = read.csv("../data/RFOIData_ AllAge_ 96-20.txt")
colnames(data) = col_names

# convert Age to numeric variables
data$Age = as.numeric(as.character(data$Age))

# how many NA in Age?
sum(is.na(data$Age))

# remove the row where Age column has missing value
data_filtered = data[!is.na(data$Age),]



# remove the row where TB. contains "IS", which is neither positive nor negative
data_filtered$TB. = as.character(data_filtered$TB.)
data_filtered = data_filtered[data_filtered$TB.!="IS",]

# create a new column called status, which maps Pos to 1, Neg to 0
data_filtered$status = 0
data_filtered$status[data_filtered$TB.=="Pos"] = 1


data_filtered$Sex = as.character(data_filtered$Sex)

# remove the rows where sex is neither male nor female
data_filtered = data_filtered[data_filtered$Sex=="male" | data_filtered$Sex=="female",]

# create a new column called covar, which maps Male to 1 and Female to 0
data_filtered$covar = 0
data_filtered$covar[data_filtered$Sex=="male"] = 1

# fit the model with a subset
#data_subset = data_filtered[1:100,]
#data_subset = data_filtered

# fit the loglogistic model with estimation of u
# res1 = fit_force_of_infection_model(data_subset$Age, data_subset$status, 400, 'loglogistic', 
#                              c(0.25,-0.8,0.7,-0.25), data_subset$covar)

# fit the loglogistic model with no infection-asscociated mortality but with sex effect
# res2 = fit_force_of_infection_model(data_subset$Age, data_subset$status, 400, 'loglogistic', 
#                                     c(0.0,-0.8,0.7,-0.25), fixmort=T, covars=data_subset$covar)

# fit the Bernoulli/Cohen model with a sex effect
# res3 = fit_force_of_infection_model(data_subset$Age, data_subset$status, 1, 'exponential', 
#                                     c(0.4,-0.8,-0.25), covars=data_subset$covar)

# load the dbf file and add spatial convariate
library(foreign)
shape_file = read.dbf("../data/sectionsubset.dbf")

unique(data_filtered$Town)
unique(data_filtered$Range)
unique(data_filtered$Section)

unique(shape_file$TOWN)
unique(shape_file$RANGE)
unique(shape_file$SECTION)

data_filtered$Town = as.character(data_filtered$Town)
is_in(unique(data_filtered$Town), unique(shape_file$TOWN))

# make spellings consistent
data_filtered$Town[data_filtered$Town=="28n"] = "28N"
data_filtered$Town[data_filtered$Town=="31n"] = "31N"
is_in(unique(data_filtered$Town), unique(shape_file$TOWN))

data_filtered$Range = as.character(data_filtered$Range)
is_in(unique(data_filtered$Range), unique(shape_file$RANGE))
data_filtered$Range[data_filtered$Range=="06e"] = "06E"
data_filtered$Range[data_filtered$Range=="04e"] = "04E"
is_in(unique(data_filtered$Range), unique(shape_file$RANGE))

data_filtered$Section = as.character(data_filtered$Section)
is_in(unique(data_filtered$Section), unique(shape_file$SECTION))
data_filtered$Section[data_filtered$Section=="1"] = "01"
data_filtered$Section[data_filtered$Section=="2"] = "02"
data_filtered$Section[data_filtered$Section=="3"] = "03"
data_filtered$Section[data_filtered$Section=="4"] = "04"
data_filtered$Section[data_filtered$Section=="5"] = "05"
data_filtered$Section[data_filtered$Section=="6"] = "06"
data_filtered$Section[data_filtered$Section=="7"] = "07"
data_filtered$Section[data_filtered$Section=="8"] = "08"
data_filtered$Section[data_filtered$Section=="9"] = "09"
is_in(unique(data_filtered$Section), unique(shape_file$SECTION))

# see which row contains strange factor level
row_index_to_be_removed = c(where_s_in_dataframe("Q", data_filtered, c(4,5,6)),
                            where_s_in_dataframe("31_32_or_33", data_filtered, 6),
                            where_s_in_dataframe("4 or 5 or 6", data_filtered, 6))
length(row_index_to_be_removed)

# remove the rows which contains the strange factor level
data_filtered = data_filtered[-row_index_to_be_removed,]

# merge the dataset with spatial locations
data_cleaned <- merge(data_filtered,shape_file,by.x = c("Town","Range","Section"), by.y = c("TOWN","RANGE","SECTION"),
                      all.x = T)


# nimble convert adjacency inforamtion into the required format
library(nimble)
adj_info = read.csv("../data/Adj.txt",header = F)
adj_num = read.csv("../data/Num.txt",header = F)[,1]

adj = rep(0, sum(adj_num))
count = 1
for (i in 1:nrow(adj_info)){
  for (j in 1:ncol(adj_info)){
    if (!is.na(adj_info[i,j])){
      adj[count] = adj_info[i,j]
      count = count+1
    }
  }
}

weights = rep(1, sum(adj_num))
num = adj_num

# add spatial id based on neighborhood idenfication to the original data
library(readxl)
tbdata_adj <- read_excel("../data/tbdata_adj.xlsx")
tbdata_adj$Section = as.character(tbdata_adj$SECTION)
tbdata_adj$Town = tbdata_adj$TOWN
tbdata_adj$Range = tbdata_adj$RANGE
tbdata_adj$Section[tbdata_adj$Section=="1"] = "01"
tbdata_adj$Section[tbdata_adj$Section=="2"] = "02"
tbdata_adj$Section[tbdata_adj$Section=="3"] = "03"
tbdata_adj$Section[tbdata_adj$Section=="4"] = "04"
tbdata_adj$Section[tbdata_adj$Section=="5"] = "05"
tbdata_adj$Section[tbdata_adj$Section=="6"] = "06"
tbdata_adj$Section[tbdata_adj$Section=="7"] = "07"
tbdata_adj$Section[tbdata_adj$Section=="8"] = "08"
tbdata_adj$Section[tbdata_adj$Section=="9"] = "09"
tbdata_adj = tbdata_adj %>% select(Section, Range, Town, Adj_ID)

# merge data_cleaned and tbdata_adj with key = ("Section","Range","Town")
data_cleaned <- merge(data_cleaned, tbdata_adj, all.x = T)
data_cleaned <- data_cleaned[order(data_cleaned$Adj_ID),]


# convert Tag_Date to how many 0.5 years passed from the starting time

data_cleaned$Tag_Date = as.Date(data_cleaned$Tag_Date)

# find outliers, delete it
hist(data_cleaned$Tag_Date,"years")
outlier_date = sort(data_cleaned$Tag_Date)[1:3]
data_cleaned = data_cleaned[!(data_cleaned$Tag_Date %in% outlier_date),]
head(sort(data_cleaned$Tag_Date),n = 100)
tail(sort(data_cleaned$Tag_Date),n = 100)
hist(data_cleaned$Tag_Date,"years")
#time_end = max(data_cleaned$Tag_Date,na.rm =TRUE)

# obtain the birth date for each animal, here the new Tag_Date is actually the birth date
for (i in 1:6){
  data_cleaned$Tag_Date = data_cleaned$Tag_Date  %m-% months(data_cleaned$Age*2)
}

hist(data_cleaned$Tag_Date,"years")

head(sort(data_cleaned$Tag_Date),n = 1000)
tail(sort(data_cleaned$Tag_Date),n = 1000)

# We only care about the animals born later than 1996-05-27
time_start = as.Date("1996-05-27")
data_cleaned = data_cleaned %>% filter(Tag_Date > time_start)
data_cleaned$Tag_Date = difftime(data_cleaned$Tag_Date,time_start,units="days")
data_cleaned = data_cleaned[!is.na(data_cleaned$Tag_Date),]

# Now the Tag Date become the discrete time class where the birth date of an individual falls in
time_interval = 182
data_cleaned$Tag_Date = as.numeric(as.character(round(data_cleaned$Tag_Date/time_interval)))+1









# prepare data for nimble model

# number of spatial locations
N = length(tbdata_adj$Adj_ID)



# sex effect
data_cleaned$male = ifelse(data_cleaned$Sex == "male",1,0)

# Tage, maximum age
max(data_cleaned$Age) # 15.5
unique(data_cleaned$Age)

# use 0.5 year as a unit
data_cleaned$agem = data_cleaned$Age*2

agem = data_cleaned$agem
# pos or neg
pos = data_cleaned$status

# data_cleaned
data_cleaned = data_cleaned %>% select(Adj_ID,male,agem,status,Tag_Date)

Ttime = max(data_cleaned$Tag_Date+data_cleaned$agem)

Tage = max(data_cleaned$agem)


# adj
adj =adj
# weights
weights = weights
# num
num = num


# double check
nrow(data_cleaned)
data_cleaned$male
data_cleaned$status

# save the result in an R data file
save(data_cleaned,Tage,N,adj,weights,num,Ttime,file = "data_for_nimble.RData")


