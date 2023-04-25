id_is_numeric =F #if input id is numeric True, otherwise False
chunkSize =  900 #CSW query batch size should be less than 1000
project <- "bcs-breeding-datasets" #CSW table project
h <- curl::new_handle()
curl::handle_setopt(h, http_version = 0)

options(scipen = 999)
library(bit64)
library(aws.s3)
library(httr)
library(bigrquery)
options(httr_oob_default=FALSE) 
options(gargle_quiet = FALSE)
# file.remove(".rAzureAPI.RData")

get_kv_vault <- function(vault_config){
  #########
  # Set Proper Variables from Config File
  vault_url <- vault_config$vault_url
  vault_login_url <- paste0(vault_url,"/auth/approle/login")
  vault_credentials_payload <- paste0('{"role_id":"',vault_config$role_id, '","secret_id":"', vault_config$secret_id,'"}')
  app_kv_url <- paste0(vault_url, vault_config$app_kv_path)
  
  #########
  # Authenticate in Vault
  # Get Token
  response <- POST(url=vault_login_url,body=vault_credentials_payload)
  client_token <- content(response)$auth$client_token
  
  ########
  # Capture Key-Values from Vault
  response <- GET(app_kv_url, add_headers("X-Vault-Token" = client_token))
  print(response)
  #######
  # Return key-values list
  return(content(response)$data)
}
vault_config <- list(vault_url =  "https://vault.agro.services/v1", 
                     role_id = Sys.getenv("VAULT_APP_ROLE_ID"),
                     secret_id = Sys.getenv("VAULT_APP_ROLE_SECRET"),
                     app_kv_path = "/secret/apd-in-house-imputation/DATABASE"
)

#------------------------------------------------------------------------------#
#------------------------------- Vault Setup ----------------------------------#
#------------------------------------------------------------------------------#
vault_kv <- get_kv_vault(vault_config = vault_config)
bq_auth(path = vault_kv$gcp)



###########################
#### Pedigree Name sql ####
###########################

ped_name_sql<-function(invs1){
  inputList=invs1
  if (id_is_numeric == T){
    makeQuote = function(x) paste0( "(", paste(x,collapse=","), ")")
  } else {
    makeQuote = function(x) paste0( "('", paste(x,collapse="','"), "')")
  }
  
  inbredAll = NULL
  
  for (i in 1:length(inputList)){
    queryString <- paste0("
   select distinct
     germ.germplasm_id              as input_germplasm_id, 
     germ.pedigree_name             as input_pedigree_name,
     gm.generation                  as input_generation, 
     gen                            as generation_difference, 
     gm2.generation                 as output_generation,  
     germ2.pedigree_name            as output_pedigree_name, 
     germ2.germplasm_id             as output_germplasm_id,
     inv2.barcode                   as output_barcode, 
     inv2.inventory_id              as output_inventory_id,
     inv2.ASORT1                    as output_ASORT1  
   from 
   `bcs-breeding-datasets.obs.breeding_midas_germplasm` germ
    join `bcs-breeding-datasets.obs.breeding_midas_genetic_material` gm
   on germ.germplasm_id = gm.germplasm_id
   left join `bcs-breeding-datasets.breeding_genomics.inbred_generation_mapping` #view for inbred relationship thru gmID
   on V1 = gm.genetic_material_id
   join  `bcs-breeding-datasets.obs.breeding_midas_genetic_material` gm2
   on gm2.genetic_material_id = V2
   join `bcs-breeding-datasets.obs.breeding_midas_germplasm` germ2
   on germ2.germplasm_id = gm2.germplasm_id 
   join `bcs-breeding-datasets.obs.breeding_midas_inventory` inv2
   on gm2.genetic_material_id = inv2.genetic_material_id
   where 
   1=1
   and gm2.generation <> 'F1' 
   -- and gm2.generation <> 'DH0'
   and germ.pedigree_name in ", makeQuote(inputList[[i]]))#inpupt condition
    
    tmp <- bq_table_download(bq_project_query(project, queryString, quiet = T), bigint = "character")
    names(tmp) <- toupper(names(tmp))
    inbredAll = plyr::rbind.fill(inbredAll, tmp)
    cat("query batch ", i, "out of ", length(inputList), "is done...")
  }
  
  return(inbredAll)
}




###########################
#### InVBID/BARCODE sql ####
###########################

bar_inv_sql<-function(invs){
  inputList=invs
  if (id_is_numeric == T){
    makeQuote = function(x) paste0( "(", paste(x,collapse=","), ")")
  } else {
    makeQuote = function(x) paste0( "('", paste(x,collapse="','"), "')")
  }
  
  inbredAll = NULL
  for (i in 1:length(inputList)){
    queryString <- paste0(
      paste0("
  with
  inputID as (
  select distinct germ1.germplasm_id, germ1.pedigree_name, gm1.generation, gm1.genetic_material_id,inv1.barcode from        
  `bcs-breeding-datasets.obs.breeding_midas_inventory` inv1
  join `bcs-breeding-datasets.obs.breeding_midas_genetic_material` gm1
  on inv1.genetic_material_id = gm1.genetic_material_id
  join `bcs-breeding-datasets.obs.breeding_midas_germplasm` germ1
  on germ1.germplasm_id = gm1.germplasm_id
  where inv1.barcode in", makeQuote(inputList[[i]]),")"),
      
      #where inv1.barcode in ('inputList[i]')
      "
select distinct
i.germplasm_id as input_germplasm_id,
i.pedigree_name as input_pedigree_name,
i.generation as input_generation,
i.barcode as input_barcode,
gen as generation_difference,
gm2.generation as output_generation,
germ2.pedigree_name as output_pedigree_name,
germ2.germplasm_id as output_germplasm_id,
inv2.barcode as output_barcode,
inv2.inventory_id as output_inventory_id
--,
-- inv2.ASORT1 as output_ASORT1
from
inputID i
left join `bcs-breeding-datasets.breeding_genomics.inbred_generation_mapping` #view for inbred relationship thru gmID
on V1 = i.genetic_material_id
join `bcs-breeding-datasets.obs.breeding_midas_genetic_material` gm2
on gm2.genetic_material_id = V2
join `bcs-breeding-datasets.obs.breeding_midas_germplasm` germ2
on germ2.germplasm_id = gm2.germplasm_id
join `bcs-breeding-datasets.obs.breeding_midas_inventory` inv2
on gm2.genetic_material_id = inv2.genetic_material_id
where
1=1
and gm2.generation <> 'F1'
-- and gm2.generation <> 'DH0'
    ")
    
  
  tmp <- bq_table_download(bq_project_query(project, queryString, quiet = T), bigint = "character")
  names(tmp) <- toupper(names(tmp))
  inbredAll = plyr::rbind.fill(inbredAll, tmp)
  cat("query batch ", i, "out of ", length(inputList), "is done...")
}
return(inbredAll)
}


#sql
# @ id is either pedigree name or barcode as input to sql 

sql_origin=function(miss_inv3,season,year,method){
  #method="peds"
  #miss_inv3=peds
  
  x1=as.character(miss_inv3$id)
  max=10
  y1=seq_along(miss_inv3$id)
  invs1 <- split(x1, ceiling(y1/max))
  
  if(method=="peds"){
    Allinbred1=ped_name_sql(invs1)
    nam=paste0("/ycao1/ITG/Cucumber/sql_pedigree_name_",season,year,".Rdata")
    s3save(Allinbred1, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = nam)}
  
  if(method=="barcode"){
    Allinbred3=bar_inv_sql(invs1)  
    nam=paste0("/ycao1/ITG/Cucumber/sql_invbid_name_",season,year,".Rdata")
    s3save(Allinbred3, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = nam)}
}




