
##--------------------------------------------------------

# MOLECULAR SURVEY FOR PATHOGENS AND MARKERS OF PERMETHRIN RESISTANCE IN HUMAN HEAD LICE (PHTHIRAPTERA: PEDICULIDAE) FROM MADAGASCAR
# Matthew Anderson
# Eremeeva Lab
# Febuary 25, 2019

##--------------------------------------------------------

# our table needs to be stored as a matrix 
 kdr <- matrix(c(22,0,17,1,4,3,63,7,18,0,6,5,8,0,3,1,0,1), nrow = 6)
 
 # labels for our table 
 dimnames(kdr) <- list("Location" = c("Ranomafana","Tsararano","Ambatolahy","Ranovao","Ambatovory","Mangevo"), "Genotype" = c("S/S","S/R","R/R"))
 
 # exact test instead of chi-squared because we have "0" counts in our data 
 fisher.test(kdr,alternative = "two.side")

# package for hardy weinberg test 
 #install.packages("HardyWeinberg")
 library("HardyWeinberg")
 
 # vector of SS, SR, RR genotypes 
 hardy_weinberg_vector<-c(sum(kdr[,1]), sum(kdr[,2]), sum(kdr[,3]))
 
 # permutation test, computationally slow but does not depend on asymptotics 
 HWPerm(hardy_weinberg_vector)
 

 
 
# HO, HE table 
 
# for the resistant allele f(R) observed 
 
# double the third column  
R_count<-cbind(kdr[,1:2],kdr[,3]*2)

# p the sum of all alleles with R divided by the total times 2 since there are 2 alleles this is the observed heterozygosity 
p<-rowSums(R_count[,2:3]) / (2*rowSums(kdr[,1:3]))

# q the complement of p 
q<- 1-p

# hardy weinburg proporitons under expectation 
RR<-p^2
SR<-2*p*q
SS<-q^2


# observed R
observed<-kdr[,2]/rowSums(kdr)

#expected SR
expected<-SR


# table of Ho, He
HO_HE_table<-cbind(observed,expected)

###Inbredding Coeficient Fis per group
Fis<-(HO_HE_table[,2]-HO_HE_table[,1])/HO_HE_table[,2]

#print table 
cbind(HO_HE_table, Fis)


### overall Ho He 

kdr_all<-colSums(kdr)
R_count_all<-c(kdr_all[1:2],kdr_all[3]*2)
p_all<-sum(R_count_all[2:3]) / (2*sum(R_count_all))
q_all<-1-p_all
# Calcuating observed expected heterogeneity
RR_all<-p_all^2
 SR_all<-2*p_all*q_all
 SS_all<- q_all^2
 observed_all<-kdr_all[2]/sum(kdr_all)
expected_all<-SR_all

# Inbredding Coeficient Fis overall 
Fis_all<-(HO_HE_all_talbe[2]-HO_HE_all_talbe[1])/HO_HE_all_talbe[2]

#overall table 
HO_HE_all_talbe<-cbind(observed_all, expected_all)

# Inbredding Coeficient Fis overall 
Fis_all<-(HO_HE_all_talbe[2]-HO_HE_all_talbe[1])/HO_HE_all_talbe[2]

# overall table 
cbind(HO_HE_all_talbe, Fis_all)

# grand table
rbind(cbind(HO_HE_table, Fis),cbind(HO_HE_all_talbe, Fis_all))





 ##----------------  References  ----------------------------------##
 
 # Exploring Diallelic Genetic Markers: The HardyWeinberg Package, Jan Graffelman, May 29, 2018
 
 # https://cran.r-project.org/web/packages/HardyWeinberg/vignettes/HardyWeinberg.pdf
 