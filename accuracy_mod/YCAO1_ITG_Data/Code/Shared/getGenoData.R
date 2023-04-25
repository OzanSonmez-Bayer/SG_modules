# Check imputed processing tomato data

library(aws.s3)
library(jsonlite)
source('Code/Shared/misFunctions.R') 


AWS.Authorization('ycao1')


bucket <- '/genome-analytics-perm-space/'
path <- '/shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/'


s3load(object = paste0(path, 'ProcessingTomato20200615-1_MarkerMetadata.txt', collapse = ''), bucket = bucket)

s3read_using(read.delim2, object = paste0(path, 'ProcessingTomato20200615-1_MarkerMetadata.txt', collapse = ''), bucket = bucket)
