##BarBIQ_final_fitting_OD.r

#################################################################################################
#####Description of this code#####
# This code is used to estimate operational droplet (OD, labeled as "EDrop") of each sample and the standard error (labeled as "SE") by fitting
#################################################################################################
#####how to run this code#####
##command##
#Rscript --vanilla BarBIQ_final_fitting_OD.r inputfile_name Index_Overlap Outputfile log_file
##explaination##
#inputfile_name: output file of BarBIQ_final_merge_all_overlaps.pl
#Index_Overlap: a file contains the indexes need to be analyzed
#Outputfile: Outputfile name
#log_file: a file to save the log information 
####################################################################################################

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.19

#####code#######
if(!require(plotrix)){
    install.packages("plotrix")
    library(plotrix)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4)
   {
      stop("At least four argument must be supplied in the sequence as (inputfile_name Index_Overlap Outputfile log_file).n", call.=FALSE)
   }

cat("Now you are running BarBIQ_final_fitting_OD.r",file=args[4],append=TRUE, sep = "\n")
cat(paste("Your input file is:",args[1]),file=args[4],append=TRUE, sep = "\n")
cat(paste("Your index file is:",args[2]),file=args[4],append=TRUE, sep = "\n")
cat(paste("Your output file is:",args[3]),file=args[4],append=TRUE, sep = "\n")
cat(paste("Your log file is:",args[4]),file=args[4],append=TRUE, sep = "\n")
group1 <- read.table(args[1],
    header=T, sep="")

write.table(t(c("ID", "EDrop", "SE")),file=args[3],sep="\t",col.names = F, row.names = F, quote = F)

## ID: the index ID;
## EDrop: operational droplet (OD);
## SE: standard error.

### For the sample SX
index <- read.table(args[2], header=F, sep="\t")
indexes <- c(as.matrix(index))
for (i in indexes) ## "CHECK4" Indexes you want to analysis
    {
 sample<- i
 sampleA<-paste(c(sample,"A"),collapse = "_")
 sampleB<-paste(c(sample,"B"),collapse = "_")
 sampleO<-paste(c(sample,"O"),collapse = "_")
    group1 <- group1[(group1[,sampleA]>0 & group1[,sampleB]>0), ]
    x=log10(group1[,sampleA]*group1[,sampleB])
    y=group1[,sampleO]

data <- cbind(x,y)
data <- data.frame(data)
### Calculate the median
stock<-c()
for(i in seq(from=-0.4, to=10, by=0.2))
   {
     low=i;
     up=i+0.4;
     mid=low+0.2;
     data1<-data[data$x >low, ];
     data1<-data1[data1$x <= up, ];
    #  print(nrow(data1));
            if(nrow(data1)>2)
            {
             mean<-mean(data1[,2]);
             media<-median(data1[,2]);
             stock=rbind(stock,c(mid,mean,media));
            }
   }
stock <- data.frame(stock)
colnames(stock)=c("mid","mean","media")
stock$media<-log10(stock$media)
stock$mean<-log10(stock$mean)
stock[stock == "-Inf"]<- -2

select=stock[stock$media > -2, ];
x=select$mid
y=select$media

## Fitting
model <- lm(y ~ 1 + offset(x))
p <- summary(model)
predicted.intervals <- predict(model,data.frame(x=x),interval='confidence',
                                level=0.99)
      
# EDrop <-  read.table(outputname, header=F, sep="", colClasses=c("character"))
newdata<-t(c(sample, -model$coeff,p$coefficients[,2]))
# EDrop<-rbind(EDrop,newdata)
write.table(newdata,file=args[3],sep="\t",col.names = F, row.names = F, quote = F, append = T) 
    }

##end##






