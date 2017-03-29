des <- read.table("../Placenta-Ethnicity/Ethnicity data/des.txt", header = TRUE)
desMat <- read.table("../Placenta-Ethnicity/Ethnicity data/supplementary clinical info/desMat_all B4 placenta.txt", fill = TRUE, header = TRUE)
head(des)
head(desMat)

library(dplyr)
library(tidyr)
colnames(des) <- c("Sample_Name", "Sample_Group", "ga1", "sex1", "Ethnicity")
head(des)

newdes <- inner_join(des, desMat)
newdes
newdes2 <- select(newdes, Sample_Name, Sample_Group, Diagnosis, ga1, ga, sex1, sex, Ethnicity, Sample_Plate, Sentrix_ID, Sentrix_Position)
newdes2
newdes3 <- select(newdes2, -Diagnosis, -ga1, -sex1)
newdes3

write.csv(newdes3, file ="C:/Users/Toshiba/Documents/Github Repos/Placenta-Ethnicity/Ethnicity data/samplesheet.csv")

#forget everything below here:
#mdata <- filter(desMat, Sample_Name == row.names(des))
#rownum <- which(desMat$Sample_Name %in% row.names(des))

#df1 <- data.frame(x=rep(letters[1:26], 16))
#df2 <- data.frame(y=letters[1:4])

#df1
#df2
#which( df1$x %in% df2$y )
#table(df1[ which( df1$x %in% df2$y ), "x"])

#subset(desMat, select = Sample_Name == row.names(des))

#des <- rownames_to_column(des, "Sample_Name")
#sname <- des$Sample_Name
#head(des)
#head(sname)


#mdata <- filter(desMat, Sample_Name %in% sname)
#head(mdata)
#head(des)

#tail(mdata)
#mdata
#des

#mdata <- arrange(mdata, Sample_Name)
#des <- arrange(des, Sample_Name)
#mdata$Sample_Name == des$Sample_Name


#head(mdata$Sample_Name)
#head(des$Sample_Name)
#head(mdata)
#mdata$Ethnicity <- des$Ethnicity

#write.table(mdata, "C:/Users/vyuan/Downloads/Ethnicity data/metadata.txt", sep="\t")
