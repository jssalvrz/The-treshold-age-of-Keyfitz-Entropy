## cleaning the workspace
rm(list=ls(all=TRUE))

setwd(  "C:/Users/jmaburto.SAM/Documents/GitHub/The-treshold-age-of-Keyfitz-Entropy/Submission DR/Revision/Gompertz exercise/")

## load packages
library(data.table)
library(HMDHFDplus)
library(expint)
library(ggplot2)

# #get data from HMD

getHMDcountries()
XYZ <- c('FRATNP','SWE')
us <- "jmaburto@colmex.mx"
pw <- "kolmogorov"

HMDL <- do.call(rbind,lapply(XYZ, function(x, us, pw){
  cat(x,"\n")
  Females      <- readHMDweb(x,"fltper_1x1",username=us,password=pw)
  Females$Sex  <- "f"
  CTRY         <- Females
  CTRY$PopName <- x
  CTRY
}, us = us, pw = pw))

HMDL <- data.table(HMDL)

HMD_Counts <- do.call(rbind,lapply(XYZ, function(x, us, pw){
  cat(x,"\n")
  Deaths          <- readHMDweb(x,"Deaths_1x1",username=us,password=pw)
  Exposures       <- readHMDweb(x,"Exposures_1x1",username=us,password=pw)
  CTRY            <- cbind(Deaths, Exposures)
  CTRY$PopName <- x
  CTRY <- CTRY[,c(1,2,3,9,13)]
  CTRY
}, us = us, pw = pw))

HMD_Counts <- data.table(HMD_Counts)
names(HMD_Counts) <- c("Year","Age","Deaths","Exposures","PopName")

save(HMD_Counts,HMDL, file="Gompertz_exercise.RData")

## dimensions & data
load('Gompertz_exercise.RData')
ages <- 30:90
HMD_Counts[,lmx:= Deaths/Exposures]
HMD_Counts <- HMD_Counts[HMD_Counts$Age %in% ages,]

##--- METHOD 1: MLE + DELTA METHOD --------

## Gompertz mu and log-likelihood
GompMu <- function(pars,ages){
  ## parameters
  a <- pars[1]
  b <- pars[2]
  ## Gompertz law
  mu <- a*exp(b*ages)
  return(mu)
}
GompLL <- function(pars,ages,deaths,exposures){
  ## parameters
  a <- pars[1]
  b <- pars[2]
  ## Gompertz mu
  mu <- GompMu(pars = c(a,b),ages = ages)
  ## log-like (minus to maximize)
  llk <- -sum(deaths*log(mu)-exposures*mu)
  return(llk)
}


Results <- NULL
## optimize parameters
for (i in unique(HMD_Counts$PopName)){
    new.data <- HMD_Counts[HMD_Counts$PopName == i,]
    for (j in unique(new.data$Year)) {
      new.data2    <- new.data[new.data$Year == j,]
      fitGomp      <- optim(par=c(1e-5,0.1),fn=GompLL,ages=ages,deaths=new.data2$Deaths,exposures=new.data2$Exposures)
      pars         <- data.table(cbind(a.hat = fitGomp$par[1], b.hat = fitGomp$par[2]))
      pars$PopName <- i
      pars$Year    <- j
      Results      <- rbind(Results, pars)
      pars         <- NULL
      fitGomp      <- NULL
      new.data2    <- NULL
      }
}


Results$a.hat
Results$b.hat

#life expectancy from gompertz
Results[,e0 := (1/b.hat)*exp(a.hat/b.hat)*expint_E1(a.hat/b.hat)]

#delta from Gompertz
gama <- digamma(1)

#Results[,delta := (1/(gama + log(a.hat/b.hat) - exp(a.hat/b.hat))) + 1/exp(a.hat/b.hat)]

Results[,delta := (1/exp(a.hat/b.hat))*(1 + 1/(exp(a.hat/b.hat)*(gama + log(a.hat/b.hat)) - 1))]

Results[,delta.e0G := delta*e0]

Results[,Country:= ifelse(PopName == 'FRATNP', 'France', 'Sweden')]

#Plot


f1 <- ggplot(Results)+
  geom_line(aes(Year, delta), size =1)+
  theme_bw()+
  facet_grid(. ~ Country)+
  #scale_y_continuous(expand=c(0,0), limits = c(.6,.92))+
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
  ylab(quote(delta))
f1


pdf(file="Figure_delta.pdf",width=11,height=5,pointsize=6,useDingbats = F)
print(f1)
dev.off()


f2 <- ggplot(Results)+
  geom_line(aes(Year, delta.e0G), size =1)+
  theme_bw()+
  facet_grid(. ~ Country)+
  #scale_y_continuous(expand=c(0,0), limits = c(.6,.92))+
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
  ylab(quote(delta))
f2


data2 <- HMDL[ HMDL$Age %in% 30:100,]

data3 <- data2[, list(muertes = sum(dx)/100000), by = list(PopName,Sex,Year)]


f3 <- ggplot(data3)+
  geom_line(aes(Year, muertes), size =1)+
  theme_bw()+
  facet_grid(. ~ PopName)+
  #scale_y_continuous(expand=c(0,0), limits = c(.6,.92))+
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
  ylab(quote(delta))
f3




