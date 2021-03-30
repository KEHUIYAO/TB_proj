library(ggplot2)
library(pammtools)



# find if all elements in list A are in list B
is_in <- function(x,y){
  for (i in x){
    if (i %in% y){
      next
    }else{
      cat(i,"not in","the set\n")
    }
  }
}

# is_in(c(1,2,3,4),c(1,2,3))

# 在dataframe的某几列里找到所有含有字符串s的行的位置
where_s_in_dataframe <- function(s,data,col_ind){
	res = rep(0,nrow(data))
	for (i in 1:nrow(data)){
		if (s %in% data[i,col_ind]){
			res = c(res, i)
		}
	}

	return(res[res != 0])
}



# hazard plot
step_hazard_plot <- function(x,name,...){
  
 
  
  dots = list(...)
  ndots = length(dots)

  data = data.frame(x=rep(x,ndots),y=unlist(dots),group = rep(name,each = length(x)))
  return(ggplot(data,aes(x=x,y=y,colour = group))+geom_stephazard())
  
}


# create a function to plot the difference of traces of one parameter in different chains
plot_chains <- function(...){
  dots = list(...)
  ndots = length(dots)
  chain_length = length(dots[[1]])
  data = data.frame(x = rep(seq(1,chain_length),ndots), y = unlist(dots), chain = rep(seq(1,ndots),each = chain_length))
  data$chain = as.factor(data$chain)
  ggplot(data,aes(x,y))+geom_line(aes(color = chain))
}


# spatial heat map of log hazard ratio
plot_sptial_heatmap <- function(effect,lat,long,lat_band=10,long_band=10){
  boundary = c(min(lat),max(lat),min(long),max(long))
  x <- seq(1,lat_band)
  y <- seq(1,long_band)
  a = (boundary[2]-boundary[1]) / lat_band
  b = (boundary[4]-boundary[3]) / long_band
  data <- expand.grid(X=as.character(x), Y=as.character(y))
  z = c()
  for(i in x){
    for (j in y){
      cur_val = c()
      cur_bound = c(boundary[1]+(i-1)*a, 
                    boundary[1]+ i*a,
                    boundary[3]+(j-1)*a,
                    boundary[3]+ j*a)
      for (k in 1:length(lat)){
        if (lat[k]<cur_bound[2] & lat[k]>cur_bound[1] & long[k]<cur_bound[4] & long[k]>cur_bound[3]){
          cur_val = c(cur_val,effect[k])
        }
      
      }
      if (length(cur_val) == 0){cur_val = 0}else{cur_val = mean(cur_val)}
      
      z = c(z,cur_val)
    }
  }
  
  data$Z = z
  return(ggplot(data,aes(X,Y,fill = Z)) + geom_tile() +
    scale_fill_gradient(low="white",high = "black"))
  
  
}

