################################################################################################
# Threshold age of H
################################################################################################

library(tidyverse)
library(readr)
library(data.table)
setwd("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Threshold of H")

# Function to calculate threshold age
age.H<- function(Age, ax, dx, lx, ex){ 
  
  ax<- ax
  dx<- dx /100000
  lx<- lx /100000
  ex<- ex
  Age <-  Age
  
  Age[Age == "110+"] <-110
  Age <- as.integer(as.character(Age))
  #e-dagger at age x
  ex.bar    <- ex +ax*(ex - c(ex[-1], ex[length(ex)]))
  ex.bar[length(ex.bar)] <- ex
  ex.bar.dx <- ex.bar * dx
  e.dag.x   <-  rev(cumsum(rev(ex.bar.dx)))/ lx
  
  #Keyfitz's entropy at age x
   Hx <- e.dag.x / ex
  
  #Cumulative hazard
   cumhaz <- -log(lx)
   
  #g(x)
   gx <- cumhaz + Hx - 1- Hx[1]
   
   #wx
   
   wx <- (dx / Hx) * gx
  
   #a.dagger
   a.dag <- ex * (1- cumhaz)
   
   #get threshold age of H
    
   return(data.frame(Age,ex, lx, e.dag.x, Hx, cumhaz, gx, wx, a.dag))}

get.a <- function(Age,gx,ex, e.dag.x, a.dag){
  e0 <-  ex[1]
  f1 <- approxfun(Age,gx,method = "linear",rule=2 )
  a.H <- uniroot(function(x) f1(x),c(0,110))$root
  
  
  f2 <- approxfun(Age,e.dag.x,method = "linear",rule=2 )
  f3 <- approxfun(Age, a.dag ,method = "linear",rule=2 )
  
  a.dagger <- uniroot(function(x) f2(x) - f3(x),c(0,110))$root
  
  result <-  data.frame(e0,a.H, a.dagger)
  return(result)
}



Dat <- data.table(read_table2("DNK.fltper_1x1.txt", skip = 3))

pr <- Dat[ ,age.H(Age = Age, ax = ax, dx = dx, lx = lx, ex = ex), by = list(Year)]


pr.a <-  subset(pr, Year >=1800)


# get thresholds for all years
as <- pr.a[ ,get.a(Age = Age, gx = gx, ex = ex,  e.dag.x = e.dag.x, a.dag =a.dag), by = list(Year)]

pr.t <- merge(pr.a, as, by = "Year")

ggplot(as)+
  geom_line(aes(Year, e0), size =1)+
  geom_line(aes(Year, a.H), colour= "red")+
  geom_line(aes(Year, a.dagger), colour= "blue")+
  theme_bw()+
  scale_y_continuous(expand=c(0,0), limits = c(20,90))+
  scale_x_continuous(expand=c(0,0))+
  theme(strip.background = element_rect(fill="white"))+
  theme(legend.title=element_blank(), 
        panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "right",
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(2, "lines"))+
  xlab("Year") +
  ylab("Threshold ages and life expectancy at birth (in years)")

pr.b <- subset(pr.t, Year==2005)

ggplot(pr.b)+
  geom_line(aes(Age, gx), colour = "blue")+
  geom_line(aes(Age, a.dag / 10), colour = "red")+
  geom_line(aes(Age, e.dag.x / 10), colour = "gray60")+
  geom_vline(xintercept = pr.b$a.H[1], linetype = "dotted", colour = "blue", size = 1)+ #threshold aH
  geom_vline(xintercept = pr.b$a.dagger[1], linetype = "dotted", colour = "red", size = 1)+ #threshold aH
  geom_vline(xintercept = pr.b$ex[1], linetype = "dashed", colour = "black")+ #threshold aH
  theme_bw()+
  scale_y_continuous(expand=c(0,0), limits = c(-1,8))+
  scale_x_continuous(expand=c(0,0), limits = c(0,110))+
  theme(strip.background = element_rect(fill="white"))+
  theme(legend.title=element_blank(), 
        panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "right",
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(2, "lines"))+
  xlab("Age") +
  ylab("Threshold ages (in years)")


