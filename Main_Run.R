#Loading libraries
library(ape)
library(TransPhylo)
library(stringr)
library(reticulate)
library(lubridate)
library(readr)

#Load the ptree from the Neewick tree
ptree<-ptreeFromPhylo(read.tree('MUSCLE_beast_adjusted_tree.nwk'),dateLastSample=2018.222)
brLenChecker(ptree)

#set parameters for AutoSummary
sscalep<-1.5
sshape<-2
gscalep<-1.5
gshape<-2
pop<- 139

c<-autoSummary(ptree,pop,0.99,F,1,gshape,gscalep,sshape,scalep)

#Run TransReport
transReport(c,pop)

#Load TransReport output
py_run_file('predictVSKnown.py')

transmission_report <- read_delim("transmission_report.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
knowninfec<-data.frame(ID=transmission_report$`Actual Infector`)
knowninfec<-na.omit(knowninfec)
preddata<-data.frame(SAC_ID=transmission_report$`SAC ID`,ID=transmission_report$`Actual Infector`,PID=transmission_report$`Pred. Infector`,Dates=transmission_report$`Pred. Infect. Date`)
preddata<-na.omit(preddata)

#calculate correct transmissions
correct<-preddata[(preddata$ID==preddata$PID),4]
correct<-decimal_date(ymd(correct))
lenc<-length(correct)

#calculate incorrect transmissions
incorrect<-preddata[(preddata$ID!=preddata$PID),4]
incorrect<-decimal_date(ymd(incorrect))

#generate plotCTree 
plotCTree(c,TRUE,NA,NA,correct,incorrect)
