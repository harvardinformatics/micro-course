## SCRIPT TO PREP DATA FOR INTRODUCTION TO R WORKSHOP

if ("tidyverse" %in% rownames(installed.packages())){
  library(tidyverse)
} else {
  install.packages("tidyverse", repos = "http://ftp.ussg.iu.edu/CRAN/")
  library(tidyverse)  
}

capwords <- function(s, strict = TRUE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


housing<-read.csv("https://raw.githubusercontent.com/datasets/house-prices-us/master/data/cities-month.csv", stringsAsFactors=F, strip.white = T)
housing=housing[c(1:(length(housing)-3),length(housing))]
airpass<-as.data.frame(t(matrix(AirPassengers,12,dimnames=list(month.abb,unique(floor(time(AirPassengers)))))))
airpass$Year = rownames(airpass)
rownames(airpass)<-NULL
airpass<-airpass[c(13,1:12)]
mms<-read.table("http://www.randomservices.org/random/data/MM.txt", header=TRUE,stringsAsFactors=F)
mms$BagID = seq(1:length(mms$Red))
mms <- mms %>% tbl_df %>% gather(color, count, -Weight, -BagID) %>% arrange(BagID,color)