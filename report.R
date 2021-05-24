## ----init1, echo=T, warning=F, message=FALSE, cache=F-------------------------
#load the R packages we need. 
library("pander")
library("coin")
library("knitr")
library("rmarkdown")
library("summarytools")
library("parallel")
library("DescTools")
library("ggplot2")
library("rcompanion")

options(stringsAsFactors=F)


#x = read.csv("data_temp1.csv")
#x = read.csv("data_temp2.csv")
#x = read.csv("data_20210502.csv")
#x = read.csv("data_20210502_v2.csv")
#x = read.csv("data_20210503.csv")
#x = read.csv("data_20210506.csv")
x = read.csv("data_20210522.csv")
#write.csv(file = "varsnames3.csv",cbind(v1=colnames(x), v2=tolower(colnames(x))), row.names=F)


#load in the variable name tranformation file
#newnames = read.csv("varsnames.csv")
newnames = read.csv("varsnames3.csv")
newnames$v2 = trimws(newnames$v2)
#change variable names
for (i in 1:nrow(newnames)){
    ix = which(colnames(x) == newnames$v1[i])
    colnames(x)[ix] = newnames$v2[i]
}




#fix the titer count by taking out commas and letters
x$spike.ab.au.ml = as.numeric(gsub(x$spike.ab.au.ml, pat="[<>,A-Za-z/]", rep=""))

#every categorical variable, change value to lower case and take out white space
for (i in 1:ncol(x)){
    if (is.character(x[,i])){
        x[,i] = trimws(tolower(x[,i]), which="both")

        #race has a wierd entry
        #Error: unexpected input in "�"
        #this should fix it
        #keep only alphanumeric cahracters
        x[,i] = gsub("[^[:alnum:][:space:]/]","",x[,i])
    }
}




#ignore anything that isn't positive or negative
x$spike.antibody.result[x$spike.antibody.result %nin% c("positive", "negative")] = NA

#ignore anything that isn't yes or no
x$on.cancer.therapy.at.the.time.of.vaccine[x$on.cancer.therapy.at.the.time.of.vaccine %nin% c("yes", "no")] = NA


unique(x$race)
table(x$race)
#combine black into african american.
#make white the reference
x$race[x$race == "black"] = "african american"
x$race = factor(x$race, levels=c("white", "hispanic", "african american", "asian", "other"))


#fix the time to vax column
time1 = as.numeric(x$Time.since.vax)
#a few wild values, most likely from bad date ranges in excel
time1[time1 > 300] = NA
x$Time.since.vax = time1



#create variables for time to last cd20 and sct
d1 = as.Date(x$Date.of.last.CD20, format="%M/%d/%Y")
d2 = as.Date(x$spike.antibody.date, format="%M/%d/%Y")
x$time.since.cd20 = d2-d1
#cap it as 0
x$time.since.cd20[x$time.since.cd20 < 0] = 0
#dichotomize to 1 year
x$time.since.cd20.365 = NA
x$time.since.cd20.365[x$time.since.cd20 <= 365] = 0
x$time.since.cd20.365[x$time.since.cd20 > 365] = 1



d1 = as.Date(x$Date.of.SCT, format="%M/%d/%Y")
d2 = as.Date(x$spike.antibody.date, format="%M/%d/%Y")
x$time.since.sct = d2-d1
#cap it as 0
x$time.since.sct[x$time.since.sct < 0] = 0
x$time.since.sct.365 = NA
x$time.since.sct.365[x$time.since.sct <= 365] = 0
x$time.since.sct.365[x$time.since.sct > 365] = 1



#fix mrn 1106341  in twice
#Asiedu charity- keep the latter "relapse" with AT as assigned fellow
x = x[-(which(x$mrn == 1106341)[1]),]
x = x[!duplicated(x$mrn),]



#patients without titer
mrns.notiters = x$mrn[is.na(x$spike.ab.au.ml)]
mrns.notiters


#I've narrowed them down: MRNs: . Two patients 3578225 and 8022063 had a negative Ab however no documentation of dates/type of vaccine so probably best to exclude them as well. 
mrns.incomplete = c(6692963, 7999793, 2453443, 3310828, 3010818, 1339883, 3578225, 8022063)
mrns.incomplete


#Sorry 2 more MRNs need to be excluded completely (1 shot only) 1055424 and 3872631
mrns.astha2 = c(1055424,3872631)

#remove the problem patients
x = x[x$mrn %nin% c(mrns.incomplete, mrns.notiters, mrns.astha2),]

#set anything other than no/yes to no
x$on.chemo.during.48.hours.before.after.vaccine.[x$on.chemo.during.48.hours.before.after.vaccine. %nin% c("yes", "no")] = NA

table(x$type.of.vaccine)
main.vacs = c("johnson and johnson", "moderna", "pfizer")
#set to NA all types except the 3 mains
x$type.of.vaccine[x$type.of.vaccine %nin% main.vacs] = NA
#remove patients that do have have one of the 3 main vaccines
#x = x[x$type.of.vaccine %in% main.vacs,]

#12 patients had positive titer after 1 dose of mRNA vaccine. We could potentially keep them in seroconversion analysis but exclude them from the titer analysis? We can mention that 1 dose still gave spike ab positivity so it might be somewhat protective?
mrn.pos.tit.dose1 = c(7282368, 2854987, 2594461, 6205577, 2627465, 1254864, 6668322, 1507031, 2002420, 2063523, 3212422, 1180340)
#remove the problem patients
x = x[x$mrn %nin% mrn.pos.tit.dose1,]


#15 patients have had an ab test <7 days, (MRNs below) which we wanted to exclude from the titer analysis to compare w controls.
mrns.ab.l7 = c(2672647, 3064487, 2403599, 2599061, 1582242, 7273792, 3102925, 2054377, 1032379, 6589544, 1225666, 3011089, 5579345, 6016098, 5849706)

#we want to use this dataset when testing titers
x.186 = x[x$mrn %nin% mrns.ab.l7,]


## ----results="asis"-----------------------------------------------------------
#run a quick summary for each variable
for (v in colnames(x)){
    cat("##", v, "\n")
    cat("\n")
    if (length(unique(x[,v])) < 40){
        cat(pander(table(x[,v], useNA="ifany")))
        cat("\n")
    }else{
        cat(pander(summary(x[,v])))
        cat("\n")
    }
}



## -----------------------------------------------------------------------------





## ---- echo=F------------------------------------------------------------------

runStats <- function(x, vars1, vars2){


    for (v1 in vars1){
        for (v2 in vars2){
            cat(paste0("\n##", v1, "  vs  ", v2), "\n")

            vals1 = x[,v1]
            vals2 = x[,v2]

            df1 = data.frame(vals1, vals2)



            #if one variable is numerical and one is categorical, draw a boxplot
            if (is.numeric(vals1) & length(unique(vals2)) < 15){
      
                cat("\n###boxplot\n")
                df1$vals2 = as.factor(df1$vals2)
                g2 = ggplot(df1, aes(x=vals2, y=vals1, fill=vals2)) + 
                    geom_boxplot(outlier.colour = NA) + geom_jitter(width=0.25, height=0.25) +
                    xlab(v2) + ylab(v1) + theme(legend.position="none")
                cat("\n\n")
                print(g2)
                cat("\n\n")
                
                #cat(pander(by(vals1, vals2, summary)))
                cat(pander(by(vals1, vals2, descr, stats = c("mean", "sd", "min", "med", "max", "n.valid"))))
                #cat(pander(t.test(vals1~vals2)))
            } else{

                cat("\n###scatterplot\n")
                g1 = ggplot(df1, aes(x=vals1, y=vals2)) + 
                    geom_jitter(width=0.05, height=0.05) +
                    xlab(v1) + ylab(v2)
                cat("\n\n")
                print(g1)
                cat("\n\n")
            }

            #Sys.sleep(1)
            #plot.new()
            #dev.off()

            tab1 = table(v1=vals1, v2=vals2, dnn=c(v1, v2))
            cat("\n###table counts\n")
            if (prod(dim(tab1)) < 40){
                cat(pander(ftable(tab1)))
                #cat(pander(tab1))
            }else{
                cat ("NA\n")
            }
            cat("\n")
            cat("\n###")
            if (prod(dim(tab1)) <= 25){
                res.f = (fisher.test(tab1))
                if (res.f$p.val < 0.05){
                    cat("*Fisher exact test\n")
                }else{
                    cat("Fisher exact test\n")
                }
                cat(pander(res.f))
                cat("\n")
            }else{
                cat("Fisher exact test \n \nNA\n")
            }

            #add a kruskal wallis test
            cat("\n###")
            resa = kruskal.test(vals1~vals2)
            resb = kruskal.test(vals2~vals1)
            if (resa$p.value < 0.05 | resb$p.value < 0.05){
                cat("*Kruskal Wallis\n")
            }else{
                cat("Kruskal Wallis\n")
            }
            cat(pander(resa))
            cat("\n")
            cat(pander(resb))
            cat("\n")


            #todo:  maybe just report the p-value
            #cor.test(exer, smoke, method="kendall")
            cat("\n###")
            res1 = cor.test(as.numeric(as.factor(vals1)), as.numeric(as.factor(vals2)), method="kendall")
            if (res1$p.value < 0.05){
                cat("*KendallTauB\n")
            }else{
                cat("KendallTauB\n")
            }
            cat(pander(res1))
            cat("\n")

            #cat("\n###KendallTauB_CI")
            res1 = KendallTauB(tab1, conf.level = 0.95)
            cat(pander(KendallTauB(tab1, conf.level = 0.95)))
            cat("\n")

            #cat("\n###MannWhitney\n")
            #if (length(unique(vals1)) == 2 & is.numeric(vals2)){
            #    cat("\n")
            #    res.w = wilcox.test(vals2~vals1)
            #    cat(pander(res.w))
            #    #print(res.w)
            #    cat("\n")
            #} else if (length(unique(vals2)) == 2 & is.numeric(vals1)){
            #    cat("\n")
            #    res.w = wilcox.test(vals1~vals2)
            #    cat(pander(res.w))
            #    #print(res.w)
            #    cat("\n")
            #} else{
            #    cat ("\nNA\n")
            #    cat("\n")
            #}



            #cat("\n###Anova\n")
            #cat(pander(summary(aov(vals1~vals2))))
            #cat(pander(summary(aov(vals2~vals1))))


            #cat("\n###T-test\n")
            #if (length(unique(vals1)) == 2){
            #    cat(pander(t.test(vals2~vals1)))
            #    cat("\n")
            #} else if (length(unique(vals2)) == 2){
            #    cat(pander(t.test(vals1~vals2)))
            #    cat("\n")
            #} else{
            #    cat ("\nNA\n")
            #    cat("\n")
            #}


            #cat("\n###Linear modelling\n")
            #cat(pander(summary(lm(vals1 ~ vals2 ))))
            #cat(pander(summary(lm(vals2 ~ vals1 ))))
            #cat("\n")





            #print(colnames(x)[i])
            cat("<br/><br/>\n")
            cat("<br/><br/>\n")
            cat("<br/><br/>\n")
        }
    }


}



## ---- results="asis"----------------------------------------------------------
#load in the comparisons

#comps = read.csv("comp1.csv", header=F)
comps = read.csv("comp_Astha_20210503.csv", header=F)
#comps = read.csv("comp_Astha_5221.csv", header=F)


for (i in 1:nrow(comps)){
    runStats(x, comps[i,1], comps[i,2])
}


#run the pairwise vaccine types
runStats(x[x$type.of.vaccine %in% c("moderna", "pfizer"),], "spike.ab.au.ml", "type.of.vaccine")
runStats(x[x$type.of.vaccine %in% c("moderna", "johnson and johnson"),], "spike.ab.au.ml", "type.of.vaccine")
runStats(x[x$type.of.vaccine %in% c("pfizer", "johnson and johnson"),], "spike.ab.au.ml", "type.of.vaccine")





## ---- results="asis"----------------------------------------------------------


for (i in 1:nrow(comps)){
    runStats(x.186, comps[i,1], comps[i,2])
}

#run the pairwise vaccine types
runStats(x.186[x.186$type.of.vaccine %in% c("moderna", "pfizer"),], "spike.ab.au.ml", "type.of.vaccine")
runStats(x.186[x.186$type.of.vaccine %in% c("moderna", "johnson and johnson"),], "spike.ab.au.ml", "type.of.vaccine")
runStats(x.186[x.186$type.of.vaccine %in% c("pfizer", "johnson and johnson"),], "spike.ab.au.ml", "type.of.vaccine")



## ---- results="asis"----------------------------------------------------------

runStats(x, "spike.ab.au.ml", "malignancy.category")


## ---- warning=F---------------------------------------------------------------
vals1 = x[,"spike.ab.au.ml"]
vals2 = x[,"malignancy.category"]
pairwise.wilcox.test(vals1, vals2, p.adjust="fdr")



## ---- results="asis"----------------------------------------------------------

runStats(x.186, "spike.ab.au.ml", "malignancy.category")


## ---- warning=F---------------------------------------------------------------
vals1 = x.186[,"spike.ab.au.ml"]
vals2 = x.186[,"malignancy.category"]
pairwise.wilcox.test(vals1, vals2, p.adjust="fdr")



## -----------------------------------------------------------------------------
#read in the control patients
controls= read.csv("data_controls_20210504.csv")
controls$Titer = as.numeric(gsub(controls$Titer, pat="[<>,A-Za-z/]", rep=""))

#only look at patients that have had their 2nd dose
controls = controls[controls$Dose..==2,]



df1a = data.frame(type=x$solid.liquid, type2 = "cancer", Titer=x$spike.ab.au.ml)
df1b = data.frame(type=x.186$solid.liquid, type2 = "cancer", Titer=x.186$spike.ab.au.ml)
df2 = data.frame(type="control", type2 = "control", Titer=controls$Titer)
df3a = rbind(df1a, df2)#201
df3b = rbind(df1b, df2)#186


kruskal.test(Titer~type, df3a)
pairwise.wilcox.test(df3a$Titer, df3a$type, p.adjust="fdr")
kruskal.test(Titer~type2, df3a)

kruskal.test(Titer~type, df3b)
pairwise.wilcox.test(df3b$Titer, df3b$type, p.adjust="fdr")
kruskal.test(Titer~type2, df3b)





## ---- results="asis"----------------------------------------------------------
cat("\n##201\n")
runStats(df3a, "Titer", "type")
runStats(df3a, "Titer", "type2")

cat("\n##186\n")
runStats(df3b, "Titer", "type")
runStats(df3b, "Titer", "type2")




## ---- results="asis"----------------------------------------------------------

runStats(x, "race", "spike.ab.au.ml")
runStats(x, "race", "spike.antibody.result")



## ---- results="asis"----------------------------------------------------------

runStats(x, "Time.since.vax", "spike.ab.au.ml")
runStats(x, "Time.since.vax", "spike.antibody.result")



## ---- results="asis"----------------------------------------------------------

#time.since.sct
#spike.antibody.result
runStats(x, "time.since.sct.365", "spike.antibody.result")
runStats(x, "time.since.sct.365", "spike.antibody.result")
runStats(x, "time.since.cd20.365", "spike.ab.au.ml")
runStats(x, "time.since.cd20.365", "spike.antibody.result")


## -----------------------------------------------------------------------------

knit_exit()

