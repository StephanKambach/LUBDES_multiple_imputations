### Testing multiple imputation methods for estimating missing
### standard deviations in a real dataset
###
### Author: Stephan Kambach
### Date: 28.04.2015

library(mice)
#library(mi)
#library(Amelia)
#library(mitools)
#library(pan)

library(ggplot2)
library(metafor)

setwd("C:\\Users\\kambach\\Desktop\\aktuelle Arbeiten\\SESYNC_multiple_imputation")
dat.raw = read.csv("MA_Gerstner2014_complete_observations.csv",row.names=1)

str(dat.raw)

#remove very clear outlier of SD
dat.raw = dat.raw[-419,]

#calculate SE
dat.raw$SE..control. = dat.raw$SD..control. / sqrt(dat.raw$plot.number..control.)
dat.raw$SE..managed. = dat.raw$SD..managed. / sqrt(dat.raw$plot.number..managed.)


#visual inspection
#plot(SE..control. ~ Mean..control., data= dat.raw)
#plot(SD..control. ~ Mean..control., data= dat.raw)
#plot(SE..control. ~ plot.number..control., data= dat.raw)


#test lm with SD or SE
#test.model = lm(SD..control. ~ 1, data=dat.raw)
#test.model = lm(SE..control. ~ Mean..control. * Biodiversity.Aspect * plot.number..control. , data=dat.raw)
#summary(test.model)

#test.model2 = lm(SE..control. ~ Mean..control. + Biodiversity.Aspect + plot.number..control. , data=dat.raw )
#summary(test.model2)
#anova(test.model2)


#check for normal distribution
#hist(log(0.001 + dat.raw$SE..control.))
#hist(log(0.001 + dat.raw$Mean..control.))
#hist(log(dat.raw$plot.number..control.))
#-> apply log transformation everywhere
#dat.log = cbind(dat.raw[,c(1:2)],log(0.001 + dat.raw[,c(3:ncol(dat.raw))]))
#hist((dat.log$SE..control.))
#hist((dat.log$Mean..control.))
#hist((dat.log$plot.number..control.))






#package mice

#imputations.possible = c("pmm","norm","norm.nob","norm.boot","norm.predict","mean","2l.norm","2l.pan",
                         "2lonly.mean","2lonly.norm","2lonly.pmm","quadratic","logreg","logreg.boot",
                         "polyreg","polr","lda","cart","rf","ri")


#mean_control = dat.raw$Mean..control.
#n_control = dat.raw$plot.number..control. 
#sd_control = dat.raw$SE..control. 
#mean_managed = dat.raw$Mean..managed.
#n_managed  = dat.raw$plot.number..managed.
#sd_managed = dat.raw$SE..managed.
#sd_or_se = "sd"
#log_trans = "no"
#del_random = "yes"
#del.min = 0.02
#del.max = 0.5
#del.step = 0.005
#del.rate=0.5
#methods = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart")




test.mi = function(mean_control, n_control, sd_control, 
                   mean_managed, n_managed, sd_managed,
                   sd_or_se, methods, log_trans,del_random,
                   del.min, del.max, del.step,Fisher_step){
  
  df.results = data.frame("method"=NA, "del_rate"=NA,"sd_or_se"=NA,"log_trans"=NA,"del_random"=NA,
                          "grand_mean_full_RR" =NA, "grand_mean_full_RR_se"=NA,
                          "grand_mean_full_SMD"=NA, "grand_mean_full_SMD_se"=NA,
                          "R2_imp_control"=NA,"R2_imp_managed"=NA,
                          "R2_es_RR" =NA, "R2_es_SMD"=NA,
                          "grand_mean_imp_RR" =NA, "grand_mean_imp_RR_se"=NA,
                          "grand_mean_imp_SMD"=NA,"grand_mean_imp_SMD_se"=NA)[0,]
  
  #calculate full model, either with sd or se
  if(sd_or_se %in% "sd"){
  effect_sizes.rom.full = escalc(measure = "ROM", m1i = mean_managed, n1i = n_managed, sd1i = sd_managed,
                                 m2i = mean_control, n2i = n_control,sd2i = sd_control)
  effect_sizes.smd.full = escalc(measure = "SMD", m1i = mean_managed, n1i = n_managed, sd1i = sd_managed,
                                 m2i = mean_control, n2i = n_control,sd2i = sd_control )
  rma.RR.full = rma(yi=effect_sizes.rom.full$yi, vi = effect_sizes.rom.full$vi, method="REML")
                       #mods=~ dat.raw$Biodiversity.Aspect, method="REML")
  rma.smd.full = rma(yi=effect_sizes.smd.full$yi, vi = effect_sizes.smd.full$vi,method="REML")
                      #mods=~ dat.raw$Biodiversity.Aspect, method="REML")
  }else{
    effect_sizes.rom.full = escalc(measure = "ROM", m1i = mean_managed, n1i = n_managed, sd1i = sd_managed * sqrt(n_managed),
                                   m2i = mean_control, n2i = n_control,sd2i = sd_control*sqrt(n_control))
    effect_sizes.smd.full = escalc(measure = "SMD", m1i = mean_managed, n1i = n_managed, sd1i = sd_managed * sqrt(n_managed),
                                   m2i = mean_control, n2i = n_control,sd2i = sd_control * sqrt(n_control))
    rma.RR.full = rma(yi=effect_sizes.rom.full$yi, vi = effect_sizes.rom.full$vi, method="REML")
    #mods=~ dat.raw$Biodiversity.Aspect, method="REML")
    rma.smd.full = rma(yi=effect_sizes.smd.full$yi, vi = effect_sizes.smd.full$vi,method="REML")
    #mods=~ dat.raw$Biodiversity.Aspect, method="REML")
  }
  
  
  
  for(del.rate in seq(from=del.min, to = del.max, by = del.step)){
  print("start loop...")    
    
    #transform
    if(log_trans %in% "yes"){
      mean_control = log(1 + mean_control) 
      n_control= log(1 + n_control)
      sd_control = log(1 + sd_control) 
      mean_managed = log(1 + mean_managed)
      n_managed = log(1 + n_managed)
      sd_managed = log(1 + sd_managed)
    }
print("deleting...")    
    #delete according to del.rate and del_random
    sd.del.control = sd_control
    sd.del.managed = sd_managed 
      
    if(del_random == "yes"){
      sd.del.control[sample(seq(1:length(sd_control)), size = length(sd_control) * del.rate)] = NA
      sd.del.managed[which(is.na(sd.del.control))] = NA
    }else{
      sd.del.control[sample(seq(1:length(sd_control)), size = length(sd_control) * del.rate, prob=(mean_control/mean_managed)^2)] = NA
      sd.del.managed[which(is.na(sd.del.control))] = NA
    }
    
print("imputing...")
    #impute
    imputed.sd.control = list()
    imputed.sd.managed = list()
  
    #impute
    for(method in methods){
      if(method %in% "na.omit"){
        imputed.sd.control[["na.omit"]] = sd.del.control
        imputed.sd.managed[["na.omit"]] = sd.del.managed
      }else{
        imputed.sd.control[[method]] =  complete(mice(cbind(dat.raw$Biodiversity.Aspect,mean_control,n_control,sd.del.control),
                                                                            method = method, m=5, maxit =20, printFlag = FALSE))$sd.del.control
        imputed.sd.managed[[method]] =  complete(mice(cbind(dat.raw$Biodiversity.Aspect,mean_managed,n_managed,sd.del.managed),
                                                                            method = method, m=5, maxit =20, printFlag = FALSE))$sd.del.managed
      }
    }
          
    #back - transform
    if(log_trans %in% "yes"){
      mean_control = exp(mean_control) - 1
      n_control= exp(n_control) - 1
      sd_control = exp(sd_control) - 1 
      mean_managed = exp(mean_managed) - 1
      n_managed = exp(n_managed) - 1
      sd_managed = exp(sd_managed) -1
      for(i in names(imputed.sd.managed)){
      imputed.sd.control[[i]] = exp(imputed.sd.control[[i]]) - 1
      imputed.sd.managed[[i]] = exp(imputed.sd.managed[[i]]) - 1
      }
    }
   
    #if se is given, compute now sd
    if(sd_or_se %in% "se"){
      sd_control = sd_control * sqrt(n_control)
      sd_managed = sd_managed * sqrt(n_managed)
      for(method in methods){
      imputed.sd.control[[method]] = imputed.sd.control[[method]] * sqrt(n_control)
      imputed.sd.managed[[method]] = imputed.sd.managed[[method]] * sqrt(n_managed)
      }
    }

print("filter negative sds...")
    #filter out negative sds
    negative.sd = list()
    for(method in methods){
      negative.sd[[method]] = which(imputed.sd.control[[method]] <= 0 | imputed.sd.managed[[method]] <= 0 | is.na(imputed.sd.control[[method]]) | is.na(imputed.sd.managed[[method]]))
    }
    
print("calculate effect sizes...")    
    #calculate effect sizes
    effect_sizes.rom.imp = list()
    effect_sizes.smd.imp = list()
    
    for(method in methods){
    if(length(negative.sd) >= 1){
    effect_sizes.rom.imp[[method]] = escalc(measure = "ROM", m1i = mean_managed[-negative.sd[[method]]], n1i = n_managed[-negative.sd[[method]]], sd1i = imputed.sd.managed[[method]][-negative.sd[[method]]],
                                 m2i = mean_control[-negative.sd[[method]]], n2i = n_control[-negative.sd[[method]]] ,sd2i = imputed.sd.control[[method]][-negative.sd[[method]]] )
    effect_sizes.smd.imp[[method]] = escalc(measure = "SMD",vtype="UB", m1i = mean_managed[-negative.sd[[method]]], n1i = n_managed[-negative.sd[[method]]], sd1i = imputed.sd.managed[[method]][-negative.sd[[method]]],
                                      m2i = mean_control[-negative.sd[[method]]], n2i = n_control[-negative.sd[[method]]] ,sd2i = imputed.sd.control[[method]][-negative.sd[[method]]])
    }else{
      effect_sizes.rom.imp[[method]] = escalc(measure = "ROM", m1i = mean_managed, n1i = n_managed, sd1i = imputed.sd.managed[[method]],
                                        m2i = mean_control, n2i = n_control, sd2i = imputed.sd.control[[method]])
      effect_sizes.smd.imp[[method]] = escalc(measure = "SMD", vtype="UB", m1i = mean_managed, n1i = n_managed, sd1i = imputed.sd.managed[[method]],
                                        m2i = mean_control, n2i = n_control, sd2i = imputed.sd.control[[method]])
    }
    }

print("compute R squares...")
    #check R square between missing values and the correct values
    R2.lm.control = list()
    R2.lm.managed = list()
    R2.lm.rom.es = list()
    R2.lm.smd.es = list()
    
    for(method in methods){
      if(method %in% "na.omit"){
        R2.lm.control[[method]] = 0
        R2.lm.managed[[method]] = 0
        R2.lm.rom.es[[method]] = 0
        R2.lm.smd.es[[method]] = 0
      }else{    
        if(length(negative.sd[[method]]) >= 1){
          R2.lm.control[[method]] = summary(lm(imputed.sd.control[[method]][which(is.na(sd.del.control))] ~ sd_control[which(is.na(sd.del.control))]))$adj.r.squared
          R2.lm.managed[[method]] = summary(lm(imputed.sd.managed[[method]][which(is.na(sd.del.managed))] ~ sd_managed[which(is.na(sd.del.managed))]))$adj.r.squared
          R2.lm.rom.es[[method]] = summary(lm(effect_sizes.rom.imp[[method]]$yi ~ effect_sizes.rom.full$yi[-negative.sd[[method]]]))$adj.r.squared
          R2.lm.smd.es[[method]] = summary(lm(effect_sizes.smd.imp[[method]]$yi ~ effect_sizes.smd.full$yi[-negative.sd[[method]]]))$adj.r.squared
        }else{
          R2.lm.control[[method]] = summary(lm(imputed.sd.control[[method]][which(is.na(sd.del.control))] ~ sd_control[which(is.na(sd.del.managed))]))$adj.r.squared
          R2.lm.managed[[method]] = summary(lm(imputed.sd.managed[[method]][which(is.na(sd.del.managed))] ~ sd_managed[which(is.na(sd.del.managed))]))$adj.r.squared
          R2.lm.rom.es[[method]] = summary(lm(effect_sizes.rom.imputed[[method]]$yi ~ effect_sizes.rom.full$yi))$adj.r.squared
          R2.lm.smd.es[[method]] = summary(lm(effect_sizes.smd.imputed[[method]]$yi ~ effect_sizes.smd.full$yi))$adj.r.squared
        }}}
    
print("calculate grand means...")
  #calculate grand.mean
   rma.rom.del = list()
   rma.smd.del = list()
      
   for(method in methods){
    rma.rom.del[[method]] = rma(yi=as.vector(effect_sizes.rom.imp[[method]]$yi), vi = as.vector(effect_sizes.rom.imp[[method]]$vi),method ="REML",control=list(maxiter=1000,stepadj=Fisher_step))
    #mods=~ dat.raw$Biodiversity.Aspect, method="REML")
    rma.smd.del[[method]] = rma(yi=as.vector(effect_sizes.smd.imp[[method]]$yi), vi = as.vector(effect_sizes.smd.imp[[method]]$vi), method = "REML", control=list(maxiter=1000,stepadj =Fisher_step))
    #mods=~ dat.raw$Biodiversity.Aspect, method="REML")
   }
   
print("add to results table")

    #add to results.table
    for(method in methods){
       df.results = rbind(df.results,data.frame(
        "method"= method,
        "del_rate" = del.rate,
        "sd_or_se" = sd_or_se,
        "log_trans"= log_trans,
        "del_random" = del_random,
        "grand_mean_full_RR" = rma.RR.full$b[1], 
        "grand_mean_full_RR_se"= rma.RR.full$se[1],
        "grand_mean_full_SMD"= rma.smd.full$b[1], 
        "grand_mean_full_SMD_se"= rma.smd.full$se[1],
        "R2_imp_control"= R2.lm.control[[method]],
        "R2_imp_managed"= R2.lm.managed[[method]],
        "R2_es_RR"= R2.lm.rom.es[[method]],
        "R2_es_SMD"= R2.lm.smd.es[[method]],
        "grand_mean_imp_RR"= rma.rom.del[[method]]$b[1],
        "grand_mean_imp_RR_se"= rma.rom.del[[method]]$se[1],
        "grand_mean_imp_SMD"= rma.smd.del[[method]]$b[1],
        "grand_mean_imp_SMD_se"= rma.smd.del[[method]]$se[1]))
      }
    print(del.rate)
    }
  print(df.results)
}



#test and compile all methods

df.results= data.frame("method"=NA, "del_rate"=NA,"sd_or_se"=NA,"log_trans"=NA,"del_random"=NA,
                             "grand_mean_full_RR" =NA, "grand_mean_full_RR_se"=NA,
                             "grand_mean_full_SMD"=NA, "grand_mean_full_SMD_se"=NA,
                             "R2_imp_control"=NA,"R2_imp_managed"=NA,
                             "R2_es_RR" =NA, "R2_es_SMD"=NA,
                             "grand_mean_imp_RR" =NA, "grand_mean_imp_RR_se"=NA,
                             "grand_mean_imp_SMD"=NA,"grand_mean_imp_SMD_se"=NA)[0,]

#run all tests
df.results = test.mi(mean_control = dat.raw$Mean..control., 
                    n_control = dat.raw$plot.number..control., 
                    sd_control = dat.raw$SD..control., 
                    mean_managed = dat.raw$Mean..managed.,
                    n_managed  = dat.raw$plot.number..managed.,
                    sd_managed = dat.raw$SD..managed.,
                    sd_or_se = "sd",
                    method = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart"),
                    log_trans = "no",
                    del_random = "yes",
                    del.min = 0.02,
                    del.max = 0.5,
                    del.step = 0.005,
                    Fisher_step=0.5)

output.file = "Mice_sd_no_trans_MAR"
write.csv(df.results,paste("Output\\",output.file,".csv",sep=""))

df.results = test.mi(mean_control = dat.raw$Mean..control., 
                     n_control = dat.raw$plot.number..control., 
                     sd_control = dat.raw$SD..control., 
                     mean_managed = dat.raw$Mean..managed.,
                     n_managed  = dat.raw$plot.number..managed.,
                     sd_managed = dat.raw$SD..managed.,
                     sd_or_se = "sd",
                     method = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart"),
                     log_trans = "no",
                     del_random = "no",
                     del.min = 0.02,
                     del.max = 0.5,
                     del.step = 0.005,
                     Fisher_step=0.5)

output.file = "Mice_sd_no_trans_MNAR"
write.csv(df.results,paste("Output\\",output.file,".csv",sep=""))


df.results = test.mi(mean_control = dat.raw$Mean..control., 
                     n_control = dat.raw$plot.number..control., 
                     sd_control = dat.raw$SD..control., 
                     mean_managed = dat.raw$Mean..managed.,
                     n_managed  = dat.raw$plot.number..managed.,
                     sd_managed = dat.raw$SD..managed.,
                     sd_or_se = "sd",
                     method = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart"),
                     log_trans = "yes",
                     del_random = "yes",
                     del.min = 0.02,
                     del.max = 0.5,
                     del.step = 0.005,
                     Fisher_step=0.5)

output.file = "Mice_sd_log_trans_MAR"
write.csv(df.results,paste("Output\\",output.file,".csv",sep=""))

df.results = test.mi(mean_control = dat.raw$Mean..control., 
                     n_control = dat.raw$plot.number..control., 
                     sd_control = dat.raw$SD..control., 
                     mean_managed = dat.raw$Mean..managed.,
                     n_managed  = dat.raw$plot.number..managed.,
                     sd_managed = dat.raw$SD..managed.,
                     sd_or_se = "sd",
                     method = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart"),
                     log_trans = "yes",
                     del_random = "no",
                     del.min = 0.02,
                     del.max = 0.5,
                     del.step = 0.005,
                     Fisher_step=0.5)

output.file = "Mice_sd_log_trans_MNAR"
write.csv(df.results,paste("Output\\",output.file,".csv",sep=""))

df.results = test.mi(mean_control = dat.raw$Mean..control., 
                     n_control = dat.raw$plot.number..control., 
                     sd_control = dat.raw$SE..control., 
                     mean_managed = dat.raw$Mean..managed.,
                     n_managed  = dat.raw$plot.number..managed.,
                     sd_managed = dat.raw$SE..managed.,
                     sd_or_se = "se",
                     method = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart"),
                     log_trans = "no",
                     del_random = "yes",
                     del.min = 0.02,
                     del.max = 0.5,
                     del.step = 0.005,
                     Fisher_step=0.2)

output.file = "Mice_se_no_trans_MAR"
write.csv(df.results,paste("Output\\",output.file,".csv",sep=""))


df.results = test.mi(mean_control = dat.raw$Mean..control., 
                     n_control = dat.raw$plot.number..control., 
                     sd_control = dat.raw$SE..control., 
                     mean_managed = dat.raw$Mean..managed.,
                     n_managed  = dat.raw$plot.number..managed.,
                     sd_managed = dat.raw$SE..managed.,
                     sd_or_se = "se",
                     method = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart"),
                     log_trans = "no",
                     del_random = "no",
                     del.min = 0.02,
                     del.max = 0.5,
                     del.step = 0.005,
                     Fisher_step=0.2)

output.file = "Mice_se_no_trans_MNAR"
write.csv(df.results,paste("Output\\",output.file,".csv",sep=""))


df.results = test.mi(mean_control = dat.raw$Mean..control., 
                     n_control = dat.raw$plot.number..control., 
                     sd_control = dat.raw$SE..control., 
                     mean_managed = dat.raw$Mean..managed.,
                     n_managed  = dat.raw$plot.number..managed.,
                     sd_managed = dat.raw$SE..managed.,
                     sd_or_se = "se",
                     method = c("na.omit","sample","mean","pmm","norm.nob","norm.boot", "norm.predict", "norm","cart"),
                     log_trans = "yes",
                     del_random = "yes",
                     del.min = 0.02,
                     del.max = 0.5,
                     del.step = 0.005,
                     Fisher_step=0.2)

output.file = "Mice_se_log_trans_MAR"
write.csv(df.results,paste("Output\\",output.file,".csv",sep=""))


df.results = read.csv("Output\\Mice_sd_no_trans_MAR.csv",row.names=1)
#df.results = read.csv("Output\\Mice_sd_no_trans_MNAR.csv",row.names=1)
#df.results = read.csv("Output\\Mice_sd_log_trans_MAR.csv",row.names=1)
#df.results = read.csv("Output\\Mice_sd_log_trans_MNAR.csv",row.names=1)
#df.results = read.csv("Output\\Mice_se_no_trans_MAR.csv",row.names=1)


#plot all outputs
all.outputs = gsub(".csv","",list.files("Output\\",pattern=".csv"))


for(output in all.outputs){
df.results = read.csv(paste("Output\\",output,".csv",sep=""), row.names=1)
  
#plot the stuff

plot.R2_control = ggplot(data=df.results) +
  geom_point(aes(x=del_rate,y=R2_imp_control )) +
  geom_smooth(aes(x=del_rate,y=R2_imp_control )) +
  facet_grid(. ~ method) + ggtitle("R� between imputed and real values in sd of control treatment")

plot.R2_managed = ggplot(data=df.results) +
  geom_point(aes(x=del_rate,y=R2_imp_managed )) +
  geom_smooth(aes(x=del_rate,y=R2_imp_managed )) +
  facet_grid(. ~ method) + ggtitle("R� between imputed and real values in sd of managed treatment")

plot.R2_RR = ggplot(data=df.results) +
  geom_point(aes(x=del_rate,y=R2_es_RR )) +
  geom_smooth(aes(x=del_rate,y=R2_es_RR )) +
  facet_grid(. ~ method) + ggtitle("R� between Response ratios from imputed sd's and real values")

plot.R2_SMD = ggplot(data=df.results) +
  geom_point(aes(x=del_rate,y=R2_es_SMD )) +
  geom_smooth(aes(x=del_rate,y=R2_es_SMD )) +
  facet_grid(. ~ method) + ggtitle("R� between Hedge's d from imputed sd's and real values")

plot.grand.mean_RR = ggplot(data=df.results) +
  geom_rect(aes(xmin = min(df.results$del_rate),xmax = max(df.results$del_rate), ymin= df.results$grand_mean_full_RR[1] - (1.96 * df.results$grand_mean_full_RR_se[1]), ymax= df.results$grand_mean_full_RR[1] + (1.96 * df.results$grand_mean_full_RR_se[1])),
            alpha = 0.2,fill="#FFFFFC",color="#FFFFFC") +
  geom_point(aes(x=del_rate,y=grand_mean_imp_RR )) +
  geom_smooth(aes(x=del_rate,y=grand_mean_imp_RR )) +
  geom_line(aes(x=del_rate,y=grand_mean_full_RR),linetype="twodash",color="orange") +
  facet_grid(. ~ method) + ggtitle("grand mean of response ratios calculated from imputet sd's \n - line denotes grand mean of full dataset")

plot.grand.mean_SMD = ggplot(data=df.results) +
  geom_rect(aes(xmin = min(df.results$del_rate),xmax = max(df.results$del_rate), ymin= df.results$grand_mean_full_SMD[1] - (1.96 * df.results$grand_mean_full_SMD_se[1]), ymax= df.results$grand_mean_full_SMD[1] + (1.96 * df.results$grand_mean_full_SMD_se[1])),
            alpha = 0.2,fill="#FFFFFC",color="#FFFFFC") +
  geom_point(aes(x=del_rate,y=grand_mean_imp_SMD )) +
  geom_smooth(aes(x=del_rate,y=grand_mean_imp_SMD )) +
  geom_line(aes(x=del_rate,y=grand_mean_full_SMD),linetype="twodash",color="orange") +
facet_grid(. ~ method) + ggtitle("grand mean of Hedge's d calculated from imputet sd's \n - line denotes grand mean of full dataset")




#plot.R2.imp.con = ggplot(data=df.results) +
#  geom_point(aes(x=del_rate,y= )) +
#  geom_smooth(aes(x=del_rate,y= sqrt((grand_mean_full_RR - grand_mean_imp_RR)^2) / sqrt(grand_mean_full_RR^2) )) +

#  facet_grid(method ~ .)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


pdf(file=paste("Output\\",output,".pdf",sep=""),width=15,height=20)
multiplot(plot.R2_control,plot.R2_managed,plot.grand.mean_RR,plot.grand.mean_SMD)
dev.off()

png(file=paste("Output\\",output,".png",sep=""),width=1000,height=700,units="px")
multiplot(plot.R2_control,plot.R2_managed,plot.grand.mean_RR,plot.grand.mean_SMD)
dev.off()

svg(file=paste("Output\\",output,".svg",sep=""),width=15,height=20)
multiplot(plot.R2_control,plot.R2_managed,plot.grand.mean_RR,plot.grand.mean_SMD)
dev.off()
}





#combine all curves in one

all.outputs = gsub(".csv","",list.files("Output\\",pattern=".csv"))

for(output in all.outputs){
  df.results = read.csv(paste("Output\\",output,".csv",sep=""), row.names=1)
  plot.list.RR = list()
  plot.list.SMD = list()
  levels(df.results$method) = c("na.omit","sample","mean","pmm","norm.nob","norm.boot","norm.predict","norm","cart")
  
  grand_means_RR = 
    ggplot(data=df.results) +
    geom_line(aes(x=del_rate, y= df.results$grand_mean_full_RR[1] - (1.96 * df.results$grand_mean_full_RR_se[1])),linetype="twodash",color="black") +
    geom_line(aes(x=del_rate, y= df.results$grand_mean_full_RR[1] + (1.96 * df.results$grand_mean_full_RR_se[1])),linetype="twodash",color="black") +
    geom_line(aes(x=del_rate,y=grand_mean_full_RR),linetype="twodash",color="orange") +
    theme_bw()+
    geom_ribbon(aes(x=del_rate,ymin= grand_mean_imp_RR-(1.96*grand_mean_imp_RR_se),ymax=grand_mean_imp_RR+(1.96*grand_mean_imp_RR_se)),alpha=0.2) + 
    geom_point(aes(x=del_rate,y=grand_mean_imp_RR),color="#3399FF") +
    ylim(min(df.results$grand_mean_imp_RR -(1.96 * df.results$grand_mean_imp_RR_se)),max(df.results$grand_mean_imp_RR+(1.96 * df.results$grand_mean_imp_RR_se))) +
    ylab("Grand mean effect size") + xlab("% of deleted SD values") + 
    ggtitle(paste(method.to.print," - Response ratio",sep="")) +
    facet_grid(. ~ method)
  
    
  
  grand_means_SMD = 
    ggplot(data=df.results) +
    geom_line(aes(x=del_rate, y= df.results$grand_mean_full_SMD[1] - (1.96 * df.results$grand_mean_full_SMD_se[1])),linetype="twodash",color="black") +
    geom_line(aes(x=del_rate, y= df.results$grand_mean_full_SMD[1] + (1.96 * df.results$grand_mean_full_SMD_se[1])),linetype="twodash",color="black") +
    geom_line(aes(x=del_rate,y=grand_mean_full_SMD),linetype="twodash",color="orange") +
    theme_bw()+
    geom_ribbon(aes(x=del_rate,ymin= grand_mean_imp_SMD-(1.96*grand_mean_imp_SMD_se),ymax=grand_mean_imp_SMD+(1.96*grand_mean_imp_SMD_se)),alpha=0.2) + 
    geom_point(aes(x=del_rate,y=grand_mean_imp_SMD),color="#3399FF") +
    ylim(min(df.results$grand_mean_imp_SMD -(1.96 * df.results$grand_mean_imp_SMD_se)),max(df.results$grand_mean_imp_SMD+(1.96 * df.results$grand_mean_imp_SMD_se))) +
    ylab("Grand mean effect size") + xlab("% of deleted SD values") + 
    ggtitle(paste(method.to.print," - Hedge's d",sep="")) +
    facet_grid(. ~ method)
  
  

  pdf(file=paste("Output\\",output,"_grand_means.pdf",sep=""),width=15,height=5)
  multiplot(grand_means_RR,grand_means_SMD)  
  dev.off()
  
  png(file=paste("Output\\",output,"_grand_means.png",sep=""),width=1500,height=500,units="px")
  multiplot(grand_means_RR,grand_means_SMD)  
  dev.off()
  
  svg(file=paste("Output\\",output,"_grand_means.svg",sep=""),width=15,height=5)
  multiplot(grand_means_RR,grand_means_SMD)  
  dev.off()
}