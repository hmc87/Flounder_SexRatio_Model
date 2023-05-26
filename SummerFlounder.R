library(ggplot2)
library(reshape2)
rm(list=ls(all=TRUE))
 
#Population simulation parameters
Npop <- 100  #initial pop size
Tmax <- 300 #max years
Tfishing <- 100 #when we start fishing
Trecovery <- 200 #when fishing ends
Ngroup <- 2 #Number of separate population groups - IE females,males
Amax <- c(14, 12) #max age across all groups (females, males)

#Life history parameters for females and males:
initialsize = c(12, 12) #cm size at "recruitment" #females, males from Terceiro 2016
Linf = c(80.6, 63.9)  #cm, average max size at maturity for use in von Bert growth function (VBGF) from NOAA 2019
k = c(0.6, 0.5)   # VBGF growth coefficients
c = c(0.0133,0.0133) #mass-at-size coefficient (in kilograms) from Joseph Charles Desfosse dissertation 1995
b = c(3, 3)   #mass-at-size exp #Should I be optimizing this?

##Proportion female at recruitment

Pfemale <- c(0.3, .25, .2, .15, .1, .05)
#We will run the population projection model for varying levels of fishing pressure:
mu_f = c(0,0.1,0.2, 0.3, 0.4, 0.5, 0.6,0.7,0.8,0.9)

lRdiff <- length(Pfemale)

#assuming Beverton Holt recruitment, these are general parameters I adjust  by eye
alpha = 500
beta =1e-200

matsize <- c(35, 35) #Safe to assume all individuals mature at 35cm (Henderson and Fabrizio 2013)

natmort <- c(0.29, 0.54) #From Maunder and Wong 2011


pmat <- array(dim=c(Amax[1],Ngroup, lRdiff), data=0) #this will hold probability of maturation each year
#pmat[Amax[2]:(Amax[1]+1),1] <- 0 #age-dependent maturity
pmat[(Amax[2]+1):Amax[1], 2,lRdiff] <- 0

###Set up vectors to store age-specific life history characters and fishing probability
W <- array(dim=c(Amax[1],Ngroup, lRdiff), data=0) #this vector will hold individual weight-at-age
L <- array(dim=c(Amax[1],Ngroup, lRdiff), data=0) #length-at-age
mu <- array(dim=c(Amax[1],Ngroup, lRdiff), data=0) #mortality-at-age (this is constant in initial case)
select <- array(dim=c(Amax[1],Ngroup), data=0) #probability of being caught at each age, given size-at-age
Fishing <- array(dim=c(Amax[1],Ngroup), data=0) #Fishing mortality at age

N <- array(dim=c(Amax[1],Tmax,Ngroup, lRdiff), data=0)  #this matrix holds population numbers in each year class over time

E <- array(dim=c(1,Tmax, lRdiff), data=0)  #holds number of eggs at each time
P <- array(dim=c(1,Tmax, lRdiff), data=0)  #holds number of larvae at each time
Catch <- array(dim=c(Amax[1],Tmax,Ngroup, lRdiff), data=0)  #store Catch number
Yield <- array(dim=c(Amax[1],Tmax,Ngroup, lRdiff), data=0) #store yield biomass
Biomass <- array(dim=c(Amax[1],Tmax,Ngroup,lRdiff),data=0) #store biomass

TEP=matrix(0, nrow=length(mu_f), ncol=lRdiff)
SPR=matrix(0, nrow=length(mu_f), ncol=lRdiff)
TS_SPR=matrix(0, nrow=length(mu_f), ncol=lRdiff)
Total_Yield = matrix(0, nrow=length(mu_f), ncol=lRdiff)
Total_Catch = matrix(0, nrow=length(mu_f), ncol=lRdiff)
Sexratio = matrix(0, nrow=length(mu_f), ncol=lRdiff)

Nhat = matrix(0, nrow=length(mu_f), ncol=lRdiff)
Bhat = matrix(0, nrow=length(mu_f), ncol=lRdiff) 
#Age-dependent size, maturation, and mortality 
#we first need to set initial conditions and define the relationships between age, size, maturation, mortality, and fishery selectivity. 
#In this simple case, we are assuming asymptotic selectivity based on length, not mass.

for (R in 1:lRdiff ) {
  
  prop <- c(Pfemale[R], (1-Pfemale[R]))#we now define Prop inside the R loop 

for(g in 1:Ngroup) { #for both sexes, where g =1 is female and g = 2 is male
  
  N[1:5,1,g, R] = Npop*prop[g] #initial population size
  
  L[1,g,R] = initialsize[g] #initial size when entering the population model 
 
  #Next define age-specific growth, mortality, and selectivity functions 
  
  for (a in 1:(Amax[1]-1)) {
    if (a <= Amax[g]) { #this if statement allows us to vary maximum age of each sex 
      
      #GROWTH 		 
      L[a+1,g, R]= Linf[g]*(1-exp(-k[g])) + (L[a,g, R])*exp(-k[g]) #length at age
      W[a,g, R]=  c[g]*L[a,g, R]^b[g] #weight at age
      
      
    
      #Maturation
      if (L[a,g, R] > matsize[g]) pmat[ a, g, R] <- 1 else pmat[a,g, R] = 0
      
      
      #FISHERY SELECTIVITY
      # Asymptotic 
      if (L[a, g, R] > 35 )  select[a,g] = 1/(1+exp(-0.12*(L[a,g, R]-60))) else select[a,g] = 0 #asymptotic for fish over 35 cm
      
    } #endif
     
    
  }#end a loop
}#next group (sex)

#plot selectivity curve
matplot(select[,],type="l", lwd=5,las=1,ylab="Selectivity",xlab="Age", col=c('grey40','grey80'), lty=1)

matplot(L[,,R], type="l", lwd=5, las=1, ylab ="Length (cm)", xlab = "Age", col=c('grey40','grey80'), lty=1) 
matplot(pmat[,,R], type="l", lwd=5, las=1, ylab ="Probability Mature", xlab = "Age", col=c('grey40','grey80'), lty=1) 


#define female egg production as a function of her mass
eggs <- 7.9*L[, 1, R]^3.4  #length-fecundity equation from Morse 1981

#Eggs will be a vector, one egg value for each age

#################################################################################################################

###Simulate population dynamics, start from an arbitrary population size and let the population reach a stable age distribution, then start fishing. The population will reach a new, fished, steady state (stable age dist).
## *Note: For the base case, recruitment is based only on Female abundance.

##Set up matrices to store population-level metrics of interest (for every level of fishing)
for (fish in 1:length(mu_f)) { #loop over all levels of fishing mortality
  for(t in 1:(Tmax-1)) {
    E[t]=sum(N[, t,1, R]*eggs[-15]*pmat[, 1, R]) #assuming spawning occurs between 1 t and the next and depends ONLY on mature females
    P[t]= E[t]  #assumes all eggs are fertilized (for now)
    
    for(g in 1:2) { #here g is sex, where g = 1 is female, and g = 2 is male  
      #need to add something about for each group there's going to be a different prop
      N[1,t+1,g, R]= (alpha*P[t]/(1 + beta*P[t]))*prop[g] 
      #This prop function, otherwise prop is this, another if else related to prop
      #this is the N_0 class that is born and recruits to the population model in the next time step... 
      # calculate  probability of fishing mortality, given selectivity function
      age <- 1
      for (age in 1:(Amax[g]-1)) {  
        
        if( t > Tfishing  & t < Trecovery) {
          #this is where we add prop changes
          Fishing[age,g] = select[age,g]*mu_f[fish]
          Fishing[Amax[g], g] = select[Amax[g], g]*mu_f[fish] 
          
        } else {
          #else prop g
          Fishing[age,g] = 0
            
        } #end if
        
        
        Catch[,t+1,g, R] <- N[,t,g, R]*(1-exp(-Fishing[,g]))  #Note this is catch numbers, not biomass, assumes natural mortality occurs later in the year than fishing
        #Catch of all the ages, the number of individuals in every age class * fishing mortality of every age class (Fishing=Fishing mortality)
        #Calculate survival
        Yield[,t+1,g,R] <- N[,t,g, R]*(1-exp(-Fishing[,g]))*W[age,g, R] #biomass yield- think TAC
        N[age+1,t+1,g, R] <- N[age,t,g, R]*exp(-natmort[g]-Fishing[age,g]) #surviving fish in each group enter the next age class in the following year, all fish get to spawn before mortality 		 
        Biomass[age,t,g, R]<-N[age+1,t+1,g, R]*W[age,g, R]
      } #end second age loop
    } #end sex g loop
  } #end t loop

  
  TEP[fish, R] = E[Tfishing-1]/1
  SPR[fish,R] = E[Tfishing+50]/E[Tfishing-1]
  TS_SPR[fish,R] = E[Tfishing+50]/TEP[1,1]
  #All fish, all sexes
  Total_Yield[fish, R] = sum(Yield[,Tfishing+20, , R])
  Total_Catch[fish, R] = sum(Catch[,Tfishing+20, , R])
  #Matrix of age at time in steady state of females/males
  Sexratio[fish, R] = sum(N[, Tfishing+20, 1, R])/sum(N[-1, Tfishing+20, 2, R])
  #Sexratio[fish] = sum(N[, Tfishing+10, 1])/(sum(N[, Tfishing+10, 2]+sum(N[, Tfishing+10, 1])))
  
  
  # right now this is a simple sum, could also be standardized relative to unfished abunandance
  
  Nhat[fish, R] =  sum(N[, Tfishing+20,  , R]) #/ sum(N[, Tfishing-20,  , R])
  Bhat[fish, R] = sum(Biomass[, Tfishing+20,  , R]) # / sum(Biomass[, Tfishing-20,  , R])
} #end fishing mortality loop 

 
#just females
# plot(colSums(N[,,1, R]), type="l")
# 
# plot(colSums(Biomass[,,1, R]),type="l")

#total Abundance
plot(colSums(N[,,1, R])+colSums(N[,,2, R]),type="l")

#total Biomass 
plot(colSums(Biomass[,,1, R])+colSums(Biomass[,,2,R]),type="l")

#Plot age structures
age_structure <- N[, Tfishing-1, ,R ]

barplot(t(age_structure), beside=T)
x<- as.data.frame(age_structure)
names(x) <- c("females", "males")
#add age column
x$age<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
age<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
dat2 <- melt(x,id.vars="age")
 
ggplot(dat2) + aes(age) + coord_flip() +  theme_bw()   +
  xlab("Age")+
  ylab("Abundance")+
  geom_bar(data =  dat2[dat2[["variable"]]=="females",],
           aes(y = -value, fill="Females"), stat="identity", fill= "grey40")+ geom_bar(data =  dat2[dat2[["variable"]]=="males",], aes(y = value, fill="Males"), stat="identity", fill="grey80") +ylim(-10e+16,10e+16)

} #end R loop

#install.packages("fields")
library(fields)
par(mfrow = c(1, 1))
image(1:6, 1:10, t(Nhat), main = "Absolute abundance", xlab = "Rdiff", ylab = "Fishing mortality (inst)")

image(1:6, 1:10, t(Bhat), main = "Absolute biomass", xlab = "Rdiff", ylab = "Fishing mortality (inst)")

par(mfrow = c(1, 1))
 image(1:6, 1:10, t(TEP), main = "TEP", xlab = "Rdiff", ylab = "Fishing mortality (inst)")


image(1:6, 1:10, t(Total_Yield), main = "Yield", xlab = "Rdiff", ylab = "Fishing mortality (inst)")

image(1:6, 1:10, t(Total_Catch), main = "Catch", xlab = "Rdiff", ylab = "Fishing mortality (inst)")

#Calculating and visualizing catch per unit effort

#Annual harvest rate
harvest<-1-exp(-mu_f)
harvestmat<-matrix(nrow = 10, ncol = 6, data = rep(harvest,6))
#Removing the first row because cannot divide by 0
harvestmat<-harvestmat[-1,]

#CPUE with actual catch numbers
Total_Catch2 = Total_Catch[-1,]
cpue_catch=Total_Catch2/harvestmat
Total_Catch2 = Total_Catch[-1,]
cpue_catch=Total_Catch2/harvestmat
image(1:6, 1:10, t(cpue_catch), main = "CPUE", xlab = "Rdiff", ylab = "Fishing mortality (inst)", col=hcl.colors(12,"Mako", rev = TRUE))

#CPUE with actual catch mass
Total_Yield2<-Total_Yield[-1,]
cpue_yield=Total_Yield2/harvestmat
image(1:6, 1:10, t(cpue_yield), main = "CPUE", xlab = "Rdiff", ylab = "Fishing mortality (inst)", col=hcl.colors(12,"Mako", rev = TRUE))

