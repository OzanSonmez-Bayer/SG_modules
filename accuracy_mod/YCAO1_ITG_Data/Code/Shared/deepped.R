library(httr)
library(jsonlite)
#install.packages("aws.s3")
library(aws.s3)
library(dplyr)

install.packages("azurequest", repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))
library(azurequest)


# if (length(grep('DOMINO',names(Sys.getenv())))>1){
#   source('Credential/Credential_sdi.R')
# }else{
#   source("Credential/Credential_sdi.R") # local file with credentials 
# }

### FROM AUTHENTICATOR 0.2
#########################################
createPingToken <- function (passwordFlow = FALSE, maxRetries = 3) 
{
  credentials = Sys.getenv("PATH_TO_CREDENTIALS", unset = NA)
  if (Sys.getenv("PATH_TO_CREDENTIALS") != "") {
    credentials <- read.table(Sys.getenv("PATH_TO_CREDENTIALS"), 
                              header = FALSE)
    client_secret <- as.character(credentials[1, ])
    client_id <- as.character(credentials[2, ])
    url <- as.character(credentials[3, ])
    if (passwordFlow) {
      username <- as.character(credentials[4, ])
      password <- as.character(credentials[5, ])
      if (is.na(username) || username == "" || is.na(password) || 
          password == "") {
        warning("password-grant flow requires username and password")
      }
      if (is.na(username) || username == "") {
        warning("username is not provided in your credentials file")
        username = ""
      }
      if (is.na(password) || password == "") {
        warning("password is not provided in your credentials file")
        password = ""
      }
    }
  }
  else {
    client_id <- Sys.getenv("CLIENT_ID")
    client_secret <- Sys.getenv("CLIENT_SECRET")
    url <- Sys.getenv("PING_URL")
    if (passwordFlow) {
      username <- Sys.getenv("USERNAME")
      password <- Sys.getenv("PASSWORD")
      if (is.na(username) || username == "" || is.na(password) || 
          password == "") {
        warning("password-grant flow requires username and password")
      }
      if (is.na(username) || username == "") {
        warning("Environmental variable USERNAME is not set")
        username = ""
      }
      if (is.na(password) || password == "") {
        warning("Environmental variable PASSWORD is not set")
        password = ""
      }
    }
  }
  auth_header <- paste0("Basic ", base64enc::base64encode(charToRaw(paste0(client_id, 
                                                                           ":", client_secret))))
  attempts <- 0
  repeat {
    if (passwordFlow) {
      response <- httr::POST(url, httr::add_headers(Authorization = auth_header), 
                             httr::add_headers(`Accept-Encoding` = "gzip, deflate"), 
                             httr::user_agent("authenticator"), httr::content_type("application/x-www-form-urlencoded;charset=UTF-8"), 
                             body = list(grant_type = "password", username = username, 
                                         password = password), encode = "form")
    }
    else {
      response <- httr::POST(url, httr::add_headers(Authorization = auth_header), 
                             httr::add_headers(`Accept-Encoding` = "gzip, deflate"), 
                             httr::user_agent("authenticator"), httr::content_type("application/x-www-form-urlencoded;charset=UTF-8"), 
                             body = list(grant_type = "client_credentials"), 
                             encode = "form")
    }
    attempts <- attempts + 1
    if (httr::status_code(response) == 200) {
      break
    }
    if (attempts >= maxRetries) {
      stop("Maximum number of retries to receive OAuth token exceed")
    }
  }
  BEARER_TOKEN <<- paste0("Bearer ", httr::content(response)$access_token)
  return(BEARER_TOKEN)
}

#createPingToken

pingSecuredGet <- function(url, userAgent = "authenticator", query = NULL, timeOut = 60) 
{
  url <- URLencode(url)
  response <- httr::GET(url, httr::add_headers(Authorization = BEARER_TOKEN), 
                        httr::user_agent(userAgent), httr::timeout(timeOut), 
                        query = query)
  if (response$status_code %in% c(302, 401, 403)) {
    httr::handle_reset(url)
    createPingToken()
    response <- pingSecuredGet(url, userAgent, timeOut)
  }
  return(response)
}


#########################################

#########################################

#### GET IDs table ###


#Crop='Watermelon' # now set in external script
#Crop = "Pepper"
# BaseInputBucketAddress <- "/shared/Veg_Phenotypes/Ref_Data/"
IDTableName <- paste0(Crop, "/",Crop,"_IDs.RData")
s3load(object = IDTableName, bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")



#########################################

### GET germplasm ID
get_germID <- function(pedigree){ # one pedigreemay have more than one germplasm id (!)
  id = unique(CropIDs$M.GERMPLASM.X_ID[CropIDs$M.GERMPLASM.PEDIGREE == pedigree])
  
  if (length(id) == 0){
    return(data.frame(pedigree = pedigree, 
                      germplasm_id = NA, 
                      n = 0, 
                      stringsAsFactors = F))
    print(paste("pedigree not found:", pedigree))
  } else{
    return(data.frame(pedigree = pedigree, 
                      germplasm_id = id, 
                      n = 1:length(id), 
                      stringsAsFactors = F))
  }
}
#get_germID <- function(pedigree){ # one pedigreemay have more than one germplasm id (!)
#  id_and_gen = unique(CropIDs[CropIDs$M.GERMPLASM.PEDIGREE == pedigree, 
#                              c("M.GERMPLASM.X_ID", "M.GENETICMATERIAL.GENERATION")])
#  names(id_and_gen) <- c("germplasm_id", "generation")
#  
#  if (nrow(id_and_gen) == 0){
#    return(data.frame(pedigree = pedigree, 
#                      germplasm_id = NA, 
#                      generation = NA, 
#                      n = 0, 
#                      stringsAsFactors = F))
#    print(paste("pedigree not found:", pedigree))
#  } else{
#    return(data.frame(pedigree = pedigree, 
#                      id_and_gen, 
#                      n = 1:nrow(id_and_gen), 
#                      stringsAsFactors = F))
#  }
#}

# test:
#get_germID("FD_HP_219_BS2RX2/PSR_206590//DON_PEDRO:@.1.1.")

# helper functions to parse number of inbreeding generations from strings with generation codes, e.g. 'INBRED' etc.
# specifically backcrosses:
backcross2inbreeding <- function(gen){
  f_gen = stringr::str_replace(gen, "BC[0-9]+", "")
  if(f_gen == ""){
    inbreeding = 0 # no inbreeding after backcross
  }else if(stringr::str_detect(f_gen, "F")){
    inbreeding = as.numeric(stringr::str_replace(f_gen, "F", "")) - 1
  }
  else {
    print(paste("invalid backcross gen:", gen))
    inbreeding = NA
  }
}
# general parser of generations character code -> inbreeding generations
generation2inbreeding <- function(gen, gen_assume_inbred = 7){
  inbreeding = ifelse(gen == "INBRED", gen_assume_inbred, # inbreds are assumed 7 gens of inbreeding
                      ifelse(stringr::str_detect(gen, "^F"), # e.g. F4 becomes 3 gens of inbreeding 
                             as.numeric(stringr::str_replace(gen, "^F", "")) - 1,
                             ifelse(stringr::str_detect(gen, "DH"),
                                    12, # DH lines are assumed to equivalent to 12 gen inbreeding
                                    ifelse(stringr::str_detect(gen, "BC"), # remove BC (b/c encoded in pedigree), only care about inbreeding afterwards
                                           backcross2inbreeding(gen), 
                                           NA)))) # if doesn't match INBRED or F? or DH, can't determine generations of inbreeding
  if (is.na(inbreeding)) print(paste("could not parse gen inbreeding from ", gen))
  return(inbreeding)
}

# test
#sapply(c("BC1", "BC12", "BC5F2", "F8", "INBRED", "", "DH0", "DH1"), generation2inbreeding) == c(0, 0, 1, 7, 7, NA, 12, 12)

# using helper functions above, get generations of inbreeding using germline id
get_generation <- function(germID){
  gen = unique(CropIDs$M.GENETICMATERIAL.GENERATION[CropIDs$M.GERMPLASM.X_ID == germID])
  
  if (length(gen) == 0){
    d = data.frame(germplasm_id = germID, 
                   generation = NA, 
                   inbreeding = NA,
                   stringsAsFactors = F)
    print(paste("germplasm id not found:", germID))
  } else{
    d = data.frame(germplasm_id = germID, 
                   generation = gen, 
                   stringsAsFactors = F) %>%
      mutate(inbreeding = sapply(gen, generation2inbreeding))
  }
  return(d)
}
# test
#get_generation("351876807262208")
#get_generation("notanID")
#########################################





# 2nd version searches for generation for all germplasm IDs associated with a pedigree
get_generation2 <- function(germID){
  pedigree = unique(CropIDs$M.GERMPLASM.PEDIGREE[CropIDs$M.GERMPLASM.X_ID == germID])
  gen = unique(CropIDs$M.GENETICMATERIAL.GENERATION[CropIDs$M.GERMPLASM.PEDIGREE == pedigree | CropIDs$M.GERMPLASM.X_ID == germID])
  #print(length(gen))  
  if (length(gen) == 0){
    #print(gen)
    #print(germID)
    #print(pedigree)
    d = data.frame(germplasm_id = germID, 
                   pedigree = ifelse(length(pedigree) == 0, NA, pedigree),
                   generation = NA, 
                   inbreeding = NA,
                   stringsAsFactors = F)
    print(paste("germplasm id not found:", germID))
  } else{
    d = data.frame(germplasm_id = germID,
                   pedigree = pedigree,
                   generation = gen, 
                   stringsAsFactors = F) %>%
      mutate(inbreeding = sapply(gen, generation2inbreeding)) 
  }
  return(d)
}


### GET binary parent through the API call

# get new token for ancestry pedigrees or fingerprint data
new_token <- function(fingerprint = F){
  if (!fingerprint){ #ancestry
    Sys.setenv(CLIENT_SECRET = ANCESTRY_CLIENT_SECRET)
    Sys.setenv(CLIENT_ID = ANCESTRY_CLIENT_ID)
    Sys.setenv(PING_URL = ANCESTRY_PING_URL)
  } else{ #fingerprint
    Sys.setenv(CLIENT_SECRET = FINGERPRINT_CLIENT_SECRET)
    Sys.setenv(CLIENT_ID = FINGERPRINT_CLIENT_ID)
    Sys.setenv(PING_URL = FINGERPRINT_PING_URL)
  }
  createPingToken()
}

# simple get binary parents from pedigree (will return more if inbreeding loop)
get_pedigree <- function(GermID, print_fails = F){
  
  CallString = paste0("https://product360.agro.services/ancestry/v1/germplasm/", GermID,"/binary-parents?depth=1&props=lineCode")
  #CallTable <- rbind(CallTable,c(GermID,CallString))
  #FamilyJson <- pingSecuredGet(CallString, timeOut = 40)
  # FamilyJson <- GET(CallString, add_headers(Authorization = BEARER_TOKEN))
  
  token_gen <- azure_api_client(Sys.getenv("CLIENT_ID"), Sys.getenv("CLIENT_SECRET"))
  FamilyJson <- azure_get(CallString,list(token = token_gen$token, token_type = token_gen$token_type))
  
  
  if (FamilyJson$status_code == '200') {
    family <- fromJSON(content(FamilyJson, as = "text", encoding = "UTF-8"))
    if (length(family$nodes) > 0){
      # convert to standard pedigree relationship format
      ped <- data.frame(family$relationships,
                        stringsAsFactors = F) %>%
        mutate(from = as.character(from)) %>% # IDs are characters, not numbers
        mutate(to = as.character(to)) %>%
        mutate(relative = paste(relation, parentalRole, sep = "_")) %>%
        dplyr::rename(ID = from) %>%
        dplyr::select(., c("ID", "relative", "to")) %>%
        tidyr::spread(., "relative", "to")
    } else{
      #print(paste(GermID, "- no parents found"))
      ped <- data.frame(ID = GermID, PARENT_FEMALE = NA, PARENT_MALE = NA, stringsAsFactors = F)}
  } else {
    if(print_fails) print(paste(GermID, "- Call failed"))
    ped <- data.frame(ID = GermID, PARENT_FEMALE = NA, PARENT_MALE = NA, stringsAsFactors = F)}
  return(ped)
}   


get_all_ancestors <- function(GermID, n_gen = 3, print_fails = F){ # look back 3 generations
  #https://product360.agro.services/ancestry/v1/germplasm/2392207065088/parents?until-first=binary-cross
  CallString = paste0("https://product360.agro.services/ancestry/v1/germplasm/", GermID,"/parents?depth=", n_gen, "&props=lineCode&props=parentalRole")
  #CallTable <- rbind(CallTable,c(GermID,CallString))
  FamilyJson <- GET(CallString, add_headers(Authorization = BEARER_TOKEN))
  
  if (FamilyJson$status_code == '200') {
    family <- fromJSON(content(FamilyJson, as = "text", encoding = "UTF-8"))
    if (length(family$nodes) > 0){
      # convert to standard pedigree relationship format
      ped <- data.frame(family$relationships,
                        stringsAsFactors = F) %>%
        mutate(from = as.character(from)) %>% # IDs are characters, not numbers
        mutate(to = as.character(to)) %>%
        mutate(relative = paste(relation, parentalRole, sep = "_")) %>%
        dplyr::rename(ID = from) %>%
        dplyr::select(., c("ID", "relative", "to")) %>%
        tidyr::spread(., "relative", "to")
    } else{
      #print(paste(GermID, "- no parents found"))
      ped <- data.frame(ID = GermID, PARENT_FEMALE = NA, PARENT_MALE = NA, stringsAsFactors = F)}
  } else {
    if(print_fails) print(paste(GermID, "- Call failed"))
    ped <- data.frame(ID = GermID, PARENT_FEMALE = NA, PARENT_MALE = NA, stringsAsFactors = F)}
  return(ped)
}   



# helper function to sort_pedigree()
sort_middle <- function(ped, N, count_down = 100){ # count prevents recursion forever
  L2 <- lapply(ped$ID, function(id) N[ped$PARENT_FEMALE == id | ped$PARENT_MALE == id | ped$ID == id])
  N2 <- unlist(lapply(L2, function(x) ifelse(sum(x == max(x, na.rm = T), na.rm = T) > 1, max(x, na.rm = T) + 1, max(x, na.rm = T)))) # if tied, bump up one
  if (sum(N2 > N) == 0) {
    print(paste("found sorted solution with", count_down, "iterations left"))
    return(ped %>%
             mutate(n = N) %>%
             arrange(desc(n)))
  }
  else if (count_down <= 0){
    print("could not come to a sorted solution -- increase count and/or check for inbreeding loops!")
    return(ped %>%
             mutate(n = N2) %>%
             arrange(desc(n))) # give best answer possible
  } 
  else return(sort_middle(ped = ped, N = N2, count_down = count_down - 1))
}

# sorting function to remove duplicates and always list parents before kids
sort_pedigree <- function(ped, nsorts = 100){
  if (is.null(ped)){
    ped_sorted <- ped
  } else{
    ped1 <- ped %>%
      arrange(desc(gen2go)) %>%
      filter(!duplicated(ID)) %>% # keep only lowest duplicate in the tree
      arrange(gen2go) %>%
      arrange(!(is.na(PARENT_FEMALE) & is.na(PARENT_MALE)))
    ped_no_parents <- ped1[is.na(ped1$PARENT_FEMALE) & is.na(ped1$PARENT_MALE), ] %>%
      mutate(n = NA)
    ped_kids_only <- ped1[!(ped1$ID %in% c(ped1$PARENT_FEMALE, ped1$PARENT_MALE)) & 
                            !(is.na(ped1$PARENT_FEMALE) & is.na(ped1$PARENT_MALE)), ] %>%
      mutate(n = NA)
    # now for 'middle' entries, make sure parents come before kids
    ped_middle <- ped1[!(is.na(ped1$PARENT_FEMALE) & is.na(ped1$PARENT_MALE)) &
                         (ped1$ID %in% c(ped1$PARENT_FEMALE, ped1$PARENT_MALE)), ]
    if (nrow(ped_middle) <= 1) ped_middle_sorted = ped_middle
    else ped_middle_sorted = sort_middle(ped = ped_middle, N = nrow(ped_middle):1, count_down = nsorts)
    ped_sorted <- bind_rows(ped_no_parents, ped_middle_sorted, ped_kids_only)
  }
  return(ped_sorted)
}

get_parents <- function(GermID_list, n_gen, print_fails = F){ # list of IDs and number of parent generations back to go
  if (length(GermID_list) == 0){ # no individuals to track
    peds <- NULL
  } else if (n_gen == 0){ # don't look for any more parents, just list these individuals as parentless
    peds <- data.frame(ID = GermID_list, PARENT_FEMALE = NA, PARENT_MALE = NA, gen2go = n_gen, stringsAsFactors = F)
  }
  else{
    peds0 <- do.call(rbind,
                     lapply(GermID_list, function(id) get_pedigree(id, print_fails = print_fails))) %>%
      distinct() %>% # because binary parents can return whole inbreeding loops, thus repeating entries
      mutate(gen2go = n_gen) 
    next_ids = unique(c(peds0$PARENT_FEMALE, peds0$PARENT_MALE)[!is.na(c(peds0$PARENT_FEMALE, peds0$PARENT_MALE))])
    next_gen = n_gen - 1
    peds <- bind_rows(peds0, get_parents(GermID_list = next_ids[!(next_ids %in% peds0$ID)], # don't repeat if already in pedigree
                                         n_gen = next_gen))
  }
  return(peds)
}


# test cases:
if (F){
  ids_small <- c("366279412416512", "597819221737472", "605351958151168")
  get_parents(GermID_list = ids_small[ids_small == 0], n_gen = 0) # returns NULL because given empty string
  get_parents(GermID_list = ids_small, n_gen = 0)
  get_parents(GermID_list = ids_small, n_gen = 1, sort_me = F)
  get_parents(GermID_list = ids_small, n_gen = 5, sort_me = T) # no parents found after 1st generation, so it stops
}



# make_ped_file() takes in a dataframe of phenotypic data, which MUST have the following column names: PEDIGREE_NAME, P1, P2, and LTYPE
# it finds all pedigree relationships n_gen deep in the pedigree (just parents would be n_gen = 1)
# and returns a sorted pedigree file appropriate for making an A matrix in asreml.
# To handle missing database entries:
# It finds all individuals listed in fts and their parents
# then, if a pedigree or it's parents aren't in fts, it pulls that information from the dataframe (P1 and P2 for parents)
# it looks for germplasm id's for P1 and P2
# All pedigree relationships are listed by germplasm ID
# any listed individuals or parents without germplasm id's in the system will have ID = their listed pedigree name

make_ped_file <- function(df, crop_ids = CropIDs, n_gen = 1){ # list of IDs and number of parent generations back to go
  # new_token(fingerprint = F)
  # unique individuals
  distinct_lines <- distinct(df[ , c("pedigree", "P1", "P2", "LTYPE")]) %>%
    #rename(pedigree = PEDIGREE_NAME) %>%
    left_join(., distinct(crop_ids[ , c("M.GERMPLASM.PEDIGREE", "M.GERMPLASM.X_ID")]), by = c("pedigree"="M.GERMPLASM.PEDIGREE")) %>%
    dplyr::rename(germplasm_id = M.GERMPLASM.X_ID) %>%
    arrange(is.na(germplasm_id)) %>% # put NA's at the bottom of the list
    filter(!duplicated(pedigree)) %>% # remove duplicated pedigrees
    mutate(ID = ifelse(is.na(germplasm_id), pedigree, germplasm_id))

  # get first generation back of parents
  parents <- do.call(rbind,
                     lapply(distinct_lines$germplasm_id[!is.na(distinct_lines$germplasm_id)], function(id) get_pedigree(id, print_fails = F))) %>%
    distinct() %>% # because binary parents can return whole inbreeding loops, thus repeating entries
    mutate(gen2go = n_gen) %>%
    filter(!(is.na(PARENT_MALE) & is.na(PARENT_FEMALE)))

  # then for any parents not found in fts, try to look them up from the data frame names listed for P1 and P2
  parents_not_found0 <-  distinct_lines %>%
    mutate(ID = ifelse(!is.na(germplasm_id), germplasm_id, pedigree)) %>% # if no germplasm ID is found for an individual, use their pedigree as their ID
    left_join(., parents, by = "ID") %>%
    filter(is.na(PARENT_FEMALE)) %>%
    dplyr::select(c("ID", "P1", "P2", "pedigree"))
  if (dim(parents_not_found0)[1] == 0){# deals with empty case (all parents found)
    parent_ids <- unique(c(parents$PARENT_FEMALE, parents$PARENT_MALE))
    peds_found <- bind_rows(parents, get_parents(GermID_list = parent_ids[!(parent_ids %in% parents$ID) & !is.na(parent_ids)], # don't repeat if already in pedigree
                                                 n_gen = n_gen - 1))
  } else{
    parents_not_found <- parents_not_found0 %>%
      tidyr::pivot_longer(data = ., cols = c("P1", "P2"), names_to = "which_parent", values_to = "parent_pedigree") %>%
      #filter(., !is.na(parent_pedigree)) %>% # for all individuals with missing parents, but parents listed in df, find the germplasm ID of their listed parents
      left_join(., distinct(crop_ids[ , c("M.GERMPLASM.PEDIGREE", "M.GERMPLASM.X_ID")]),
                by = c("parent_pedigree"="M.GERMPLASM.PEDIGREE")) %>%
      mutate(parent_ID = ifelse(!is.na(M.GERMPLASM.X_ID), M.GERMPLASM.X_ID, parent_pedigree)) %>% # parent ID is the germplasm id (if found) or else the pedigree (if germplasm id not found)
      dplyr::select(-parent_pedigree) %>%
      tidyr::pivot_wider(data = ., names_from = which_parent, values_from = parent_ID) %>%
      dplyr::rename(PARENT_FEMALE = P1) %>%
      mutate(PARENT_MALE = ifelse(is.na(P2), PARENT_FEMALE, P2)) %>% # if only 1 parent is listed, it's an inbreeding event and P2 = P1
      dplyr::select(ID, PARENT_FEMALE, PARENT_MALE) %>%
      mutate(gen2go = n_gen)
    parent_ids <- unique(c(parents$PARENT_FEMALE, parents$PARENT_MALE, parents_not_found$PARENT_FEMALE, parents_not_found$PARENT_MALE))
    peds_found <- bind_rows(parents, parents_not_found, get_parents(GermID_list = parent_ids[!(parent_ids %in% c(parents$ID, parents_not_found$ID)) & !is.na(parent_ids)], # don't repeat if already in pedigree
                                                                    n_gen = n_gen - 1))
  }

  peds_sorted <- sort_pedigree(peds_found) %>%
    mutate(order_id = 1:nrow(.)) %>%
    left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GENETICMATERIAL.GENERATION")]),
              by = c("ID"="M.GERMPLASM.X_ID")) %>%
    left_join(., distinct_lines, by = c("ID"="germplasm_id")) %>%
    # defer to information given about this genetic material. If it's a hybrid, the generation is F1
    mutate(gen = ifelse(is.na(LTYPE), M.GENETICMATERIAL.GENERATION, ifelse(LTYPE == "Hybrid",
                                                                           "F1",
                                                                           ifelse(LTYPE == "Inbred", "INBRED", M.GENETICMATERIAL.GENERATION)))) %>%
    mutate(inbreeding = sapply(gen, function(g) generation2inbreeding(g, gen_assume_inbred = 5))) %>% # assume 5 gen. inbreeding if unknown for an inbred
    arrange(order_id, desc(inbreeding)) %>% # choose most inbred
    filter(!duplicated(ID)) %>%
    mutate(inbreeding = ifelse(gen == "" | is.na(gen) & is.na(inbreeding), 5, inbreeding)) %>%# (!) Assume unknowns are inbred if no generations listed
    dplyr::select(ID, PARENT_FEMALE, PARENT_MALE, inbreeding, gen2go)

  return(list(ped = peds_sorted, n_distinct_lines = dim(distinct_lines)[1], n_parents_found = dim(parents)[1], distinct_lines = distinct_lines)) # also returns df of pedigree-germplasm link for samples w/ data
}


# make_ped_file <- function(df, crop_ids = CropIDs, n_gen){ # list of IDs and number of parent generations back to go
#   # new_token(ancestry = T)
#   new_token(fingerprint = F)
#   #crop_ids <- mutate_all(funs(as,character))
#   
#   if (!("female_germplasm_id" %in% colnames(df) & "male_germplasm_id" %in% colnames(df))){ # if either column doesn't exist, create them
#     df$female_germplasm_id <- NA
#     df$male_germplasm_id <- NA
#     print("No parent germplasm IDs listed in database. Will attempt to look up ancestry.")
#   }
#   if (!("parent.female_pedigree_name" %in% colnames(df) & "parent.male_pedigree_name" %in% colnames(df))){ # if either column doesn't exist, create them
#     df$parent.female_pedigree_name <- NA
#     df$parent.male_pedigree_name <- NA
#     print("No parent pedigrees listed in database. Will attempt to look up ancestry.")
#   }
#   # find unique individuals. Assumes colnames below all exist and focal individual has pedigree, germplasm ID and genetic material ID fields complete
#   distinct_lines0 <- df[ , c("bmrd.germplasm_id", "bmrd.genetic_material_id", "bmrd.pedigree_name", 
#                              "female_germplasm_id", "parent.female_pedigree_name", 
#                              "male_germplasm_id", "parent.male_pedigree_name")] %>% #(pedigree, "P1", "P2", "LTYPE")]) %>%
#     mutate_all(funs(as.character)) %>% # convert all columns to character columns
#     filter(!is.na(bmrd.germplasm_id)) %>%
#     distinct() %>% # remove duplicates
#     # get generations for focal individuals:
#     left_join(., distinct(crop_ids[ , c("M.GENETICMATERIAL.X_ID", "M.GENETICMATERIAL.GENERATION")]), by = c("bmrd.genetic_material_id"="M.GENETICMATERIAL.X_ID")) %>%
#     rename("gen" = "M.GENETICMATERIAL.GENERATION") %>%
#     mutate(inbreeding = sapply(gen, function(g) generation2inbreeding(g, gen_assume_inbred = 5))) %>% # assume 5 gen. inbreeding if unknown for an inbred
#     dplyr::select(-bmrd.genetic_material_id) %>%
#     rename(ID = bmrd.germplasm_id,
#            PARENT_FEMALE = female_germplasm_id,
#            PARENT_MALE = male_germplasm_id,
#            PEDIGREE = bmrd.pedigree_name,
#            PARENT_FEMALE_PEDIGREE = parent.female_pedigree_name,
#            PARENT_MALE_PEDIGREE = parent.male_pedigree_name) %>%
#     group_by(ID, PEDIGREE, 
#              PARENT_FEMALE, PARENT_FEMALE_PEDIGREE, 
#              PARENT_MALE, PARENT_MALE_PEDIGREE) %>%
#     summarise(# inbreeding is assumed 0 if no information. Or if multiple options, take minimum value.
#       inbreeding = ifelse(is.na(min(inbreeding, na.rm = T)), 0, min(inbreeding, na.rm = T))) %>%  # ACK! there are some cases where the same germplasm ID has different parents listed, e.g. distinct_lines[distinct_lines$bmrd.germplasm_id == "261092344061",]
#     ungroup()
#   
#   distinct_lines0[distinct_lines0 == ""] <- NA # turn any empty strings into NAs
#   
#   # if only female parent listed, it means this was selfed (so male parent = female parent)
#   distinct_lines0 <- distinct_lines0 %>%
#     mutate( # also convert any listed pedigrees and parent germplasm IDs to character strings
#       PARENT_FEMALE = as.character(PARENT_FEMALE),
#       PARENT_FEMALE_PEDIGREE = as.character(PARENT_FEMALE_PEDIGREE),
#       PARENT_MALE = ifelse(is.na(PARENT_MALE), PARENT_FEMALE, as.character(PARENT_MALE)),
#       PARENT_MALE_PEDIGREE = ifelse(is.na(PARENT_MALE_PEDIGREE), PARENT_FEMALE_PEDIGREE, as.character(PARENT_MALE_PEDIGREE))) %>%
#     filter(!is.na(ID)) # don't include any NA lines
#   
#   # find duplicated lines
#   duplicated_lines <- distinct_lines0$ID[duplicated(distinct_lines0$ID)]
#   
#   if (length(duplicated_lines > 0)){ # print error if there is nonredundant duplication (indicating underlying data errors)
#     print("Some germplasm IDs are duplicated with inconsistent metadata. Only the first entry below was used in analysis:")
#     duplicates <- distinct_lines0 %>%
#       filter(ID %in% duplicated_lines) %>%
#       # count observations to find the most common entry of several (e.g. repeated parents)
#       mutate(n_observed = sapply(1:nrow(.), function(i) sum(c(distinct_lines0$PARENT_FEMALE_PEDIGREE, 
#                                                               distinct_lines0$PARENT_MALE_PEDIGREE,
#                                                               distinct_lines0$ID) %in% 
#                                                               c(.$PARENT_FEMALE_PEDIGREE[i], 
#                                                                 .$PARENT_MALE_PEDIGREE[i], # don't count NA matches
#                                                                 .$ID[i])[!is.na(c(.$PARENT_FEMALE_PEDIGREE[i], 
#                                                                                   .$PARENT_MALE_PEDIGREE[i], 
#                                                                                   .$ID[i]))]))) %>%
#       arrange(desc(n_observed), ID) %>%
#       print()
#   }
#   
#   distinct_lines <- distinct_lines0 %>% # ignore any duplicate entries. omit least common.
#     mutate(n_observed = sapply(1:nrow(.), function(i) sum(c(distinct_lines0$PARENT_FEMALE_PEDIGREE, 
#                                                             distinct_lines0$PARENT_MALE_PEDIGREE,
#                                                             distinct_lines0$ID) %in% 
#                                                             c(.$PARENT_FEMALE_PEDIGREE[i], 
#                                                               .$PARENT_MALE_PEDIGREE[i], # don't count NA matches
#                                                               .$ID[i])[!is.na(c(.$PARENT_FEMALE_PEDIGREE[i], 
#                                                                                 .$PARENT_MALE_PEDIGREE[i], 
#                                                                                 .$ID[i]))]))) %>%
#     arrange(desc(n_observed), ID) %>%
#     filter(!duplicated(ID)) %>% 
#     mutate(gen2go = n_gen) %>%
#     dplyr::select(-n_observed)
#   
#   
#   
#   
#   #-------------------------------------------------------------------------------------------------------------------#
#   # THE FOLLOWING 100+ lines of code deal with edge cases where the parent germplasm and/or pedigree data is incomplete
#   
#   # some test pedigrees:
#   # WCSEJ11-2538+WCS-EJ17-2790
#   # EJ_(WDL-42-69/WCSEJ14-2649:@.0003.0002.0001.0001.)+WDL-EJ17-2803
#   # CRISBY
#   # CANT_PARSE_ME++
#   
#   
#   # split into cases with parents listed and cases without parents listed:
#   parents_listed_peds_complete <- distinct_lines %>%
#     filter(!is.na(PARENT_FEMALE) & !is.na(PARENT_FEMALE_PEDIGREE) & !is.na(PARENT_MALE_PEDIGREE))
#   
#   # are the pedigrees listed too for parents?
#   parents_listed_peds_missing <- distinct_lines %>%
#     filter(!is.na(PARENT_FEMALE) & (is.na(PARENT_FEMALE_PEDIGREE) | is.na(PARENT_MALE_PEDIGREE)))
#   if (nrow(parents_listed_peds_missing) == 0){ # if yes, all done
#     parents_listed = parents_listed_peds_complete
#   } else{ # if not, look up parent pedigrees in CropIDs database
#     parents_listed = parents_listed_peds_missing %>%
#       left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GERMPLASM.PEDIGREE")]), by = c("PARENT_FEMALE"="M.GERMPLASM.X_ID")) %>%
#       rename(PARENT_FEMALE_PEDIGREE = M.GERMPLASM.PEDIGREE) %>%
#       left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GERMPLASM.PEDIGREE")]), by = c("PARENT_MALE"="M.GERMPLASM.X_ID")) %>%
#       rename(PARENT_MALE_PEDIGREE = M.GERMPLASM.PEDIGREE) %>%
#       bind_rows(., parents_listed_peds_complete) # combine with parents already having pedigrees
#   }
#   
#   
#   # for any germplasm with missing parent pedigrees, look up parent pedigrees by germplasm ID
#   no_parents_listed <- distinct_lines %>%
#     filter(is.na(PARENT_FEMALE)) %>%
#     filter(!is.na(ID))
#   
#   # get first generation back of parents by looking up parents in ancestry database
#   if (nrow(no_parents_listed) == 0){ # all parents have germplasm listed
#     # parents is a complete list of all individuals and all possibly known parents, with pedigrees
#     parents <- dplyr::select(parents_listed, ID, PARENT_FEMALE, PARENT_MALE, PEDIGREE, PARENT_FEMALE_PEDIGREE, PARENT_MALE_PEDIGREE, inbreeding, gen2go)
#     
#   } else { # otherwise, find the pedigrees
#     print(paste(nrow(no_parents_listed), "have no parent germplasm listed..looking up germplasm in ancestry database"))
#     parents_from_ancestry <- do.call(rbind,
#                                      lapply(no_parents_listed$ID, function(id) get_pedigree(id, print_fails = F))) %>%
#       distinct() %>% # because binary parents can return whole inbreeding loops, thus repeating entries
#       filter(!(is.na(PARENT_FEMALE))) %>%
#       mutate(PARENT_MALE = ifelse(is.na(PARENT_MALE), PARENT_FEMALE, PARENT_MALE)) # if only female parent, means it's selfing progeny
#     if (nrow(parents_from_ancestry) == 0){
#       parents_from_ancestry = no_parents_listed # no parents were found
#     } else{
#       parents_from_ancestry <- parents_from_ancestry %>%
#         left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GERMPLASM.PEDIGREE")]), by = c("PARENT_FEMALE"="M.GERMPLASM.X_ID")) %>%
#         rename(PARENT_FEMALE_PEDIGREE = M.GERMPLASM.PEDIGREE) %>%
#         left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GERMPLASM.PEDIGREE")]), by = c("PARENT_MALE"="M.GERMPLASM.X_ID")) %>%
#         rename(PARENT_MALE_PEDIGREE = M.GERMPLASM.PEDIGREE) %>%
#         left_join(no_parents_listed[ , c("ID", "PEDIGREE", "inbreeding", "gen2go")], ., by = "ID") # add back in 
#     }
#     
#     
#     # parents found vs. no parents found in ancestry database
#     parents_found <- parents_from_ancestry %>%
#       filter(!is.na(PARENT_FEMALE))
#     parents_not_found <- parents_from_ancestry %>%
#       filter(is.na(PARENT_FEMALE))
#     
#     # parse pedigrees for any individuals whose parents are still not found
#     if (nrow(parents_not_found) == 0){ # all individuals had parents found in ancestry database
#       parents <-  bind_rows(parents_listed, parents_found) %>%
#         dplyr::select(., ID, PARENT_FEMALE, PARENT_MALE, PEDIGREE, PARENT_FEMALE_PEDIGREE, PARENT_MALE_PEDIGREE, inbreeding, gen2go)
#     } else{ # still individuals without parents
#       print(paste("no parent information found for", nrow(parents_not_found), "individuals in ancestry database. Parsing their pedigrees and looking up germplasm ids from parsed pedigree:"))
#       print(parents_not_found[ , c("ID", "PEDIGREE")])
#       if (sum(is.na(parents_not_found$PARENT_FEMALE_PEDIGREE)) > 0){ # continue if at least one parent has no pedigree listed..find the pedigree by parsing the offspring's pedigree
#         
#         parents_not_found <- parents_not_found %>%
#           mutate(ped_can_be_parsed = grepl('[+]', PEDIGREE)) %>% # is there a + in the pedigree
#           # try to parse the pedigree into 2 parts by splitting at the +
#           separate(PEDIGREE, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left') %>%
#           mutate(PARENT_FEMALE_PEDIGREE = ifelse(is.na(PARENT_FEMALE_PEDIGREE) & ped_can_be_parsed, 
#                                                  P1, 
#                                                  PARENT_FEMALE_PEDIGREE)) %>%
#           mutate(PARENT_MALE_PEDIGREE = ifelse(is.na(PARENT_MALE_PEDIGREE) & ped_can_be_parsed, 
#                                                P2,
#                                                PARENT_MALE_PEDIGREE)) %>%
#           dplyr::select(-c(P1, P2, ped_can_be_parsed)) %>%
#           # only keep newly parsed parent pedigrees if they don't match the offspring's pedigree
#           mutate(PARENT_FEMALE_PEDIGREE = ifelse(PEDIGREE %in% c(PARENT_FEMALE_PEDIGREE, PARENT_FEMALE_PEDIGREE), NA, PARENT_FEMALE_PEDIGREE),
#                  PARENT_MALE_PEDIGREE = ifelse(PEDIGREE %in% c(PARENT_FEMALE_PEDIGREE, PARENT_FEMALE_PEDIGREE), NA, PARENT_MALE_PEDIGREE))
#       }
#       
#       # Every parent should now have a pedigree if it's possible to look one up. Now look up germplasm IDs from these pedigrees:
#       parents_not_found <- parents_not_found %>%
#         dplyr::select(c("PARENT_FEMALE_PEDIGREE", "PARENT_MALE_PEDIGREE", "ID")) %>%
#         left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GERMPLASM.PEDIGREE")]) %>%
#                     filter(!is.na(M.GERMPLASM.X_ID)) %>%
#                     filter(!duplicated(M.GERMPLASM.PEDIGREE)), # keep only one non-NA possible germplasm_id match
#                   by = c("PARENT_FEMALE_PEDIGREE"="M.GERMPLASM.PEDIGREE")) %>%
#         rename(PARENT_FEMALE = M.GERMPLASM.X_ID) %>%
#         left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GERMPLASM.PEDIGREE")]) %>%
#                     filter(!is.na(M.GERMPLASM.X_ID)) %>%
#                     filter(!duplicated(M.GERMPLASM.PEDIGREE)),
#                   by = c("PARENT_MALE_PEDIGREE"="M.GERMPLASM.PEDIGREE")) %>%
#         rename(PARENT_MALE = M.GERMPLASM.X_ID) %>%
#         left_join(parents_not_found[ , c("ID", "PEDIGREE", "inbreeding", "gen2go")], ., by = "ID") %>%
#         # if absolutely no germplasm can be found in the system, use the pedigree as the germplasm ID
#         mutate(PARENT_FEMALE = ifelse(is.na(PARENT_FEMALE), PARENT_FEMALE_PEDIGREE, PARENT_FEMALE),
#                PARENT_MALE = ifelse(is.na(PARENT_MALE), PARENT_MALE_PEDIGREE, PARENT_MALE))
#       
#       # pedigree was not successfully parsed for the following individuals:
#       print("Attempted to parse pedigree but still did not find both parents for the following individuals:")
#       parents_not_found %>%
#         filter(is.na(PARENT_FEMALE) | is.na(PARENT_MALE)) %>%
#         dplyr::select(-c(inbreeding, gen2go)) %>%
#         print(.)
#       
#       # fit the model with dummy germplasm ID = pedigree. Note: there may still be some individuals with NA listed for parents
#       parents <- bind_rows(parents_listed, parents_found, parents_not_found) %>%
#         dplyr::select(., ID, PARENT_FEMALE, PARENT_MALE, PEDIGREE, PARENT_FEMALE_PEDIGREE, PARENT_MALE_PEDIGREE, inbreeding, gen2go) 
#     }
#   }
#   
#   
#   
#   #-------------------------------------------------------------------------------------------------------------------#  
#   # Now we've used the trait database and ancestry database to create as complete as possible 'parents' dataframe with 
#   # germplasm ID's and pedigrees for all focal individuals and all their parents too.
#   
#   child_ids <- parents[ , c("ID", "PEDIGREE")] # focal individuals
#   #parent_ids <- unique(c(parents$PARENT_FEMALE, parents$PARENT_MALE, parents_not_found$PARENT_FEMALE, parents_not_found$PARENT_MALE))
#   parent_ids <- data.frame(ID = c(parents$PARENT_FEMALE, parents$PARENT_MALE),
#                            PEDIGREE = c(parents$PARENT_FEMALE_PEDIGREE, parents$PARENT_MALE_PEDIGREE),
#                            stringsAsFactors = F) %>%
#     distinct() %>%
#     arrange(PEDIGREE) %>%
#     filter(!(ID %in% child_ids$ID)) %>% # if individual is a parent and child in the experiment, default to child group only
#     filter(!is.na(ID) & !duplicated(ID)) # only keep non-NA, non-duplicated parents
#   
#   # find the rest of the generations (beyong parents, back to n_gen):
#   next_gens = get_parents(GermID_list = parent_ids$ID[!(parent_ids %in% parents$ID)], # don't repeat if already in pedigree
#                           n_gen = n_gen - 1) %>%
#     filter(!duplicated(ID)) %>% # keep earliest found in the pedigree when removing duplicates
#     left_join(., CropIDs[ , c("M.GERMPLASM.X_ID", "M.GENETICMATERIAL.GENERATION")], by = c("ID"="M.GERMPLASM.X_ID")) %>%
#     rename("gen" = "M.GENETICMATERIAL.GENERATION") %>%
#     mutate(inbreeding = sapply(gen, function(g) generation2inbreeding(g, gen_assume_inbred = 5))) %>% # assume 5 gen. inbreeding if unknown for an inbred
#     dplyr::select(-gen) %>%
#     group_by(ID, PARENT_FEMALE, PARENT_MALE, gen2go) %>%
#     summarise(# inbreeding is assumed 0 if no information. Or if multiple options, take minimum value.
#       inbreeding = ifelse(is.na(max(inbreeding, na.rm = T)), 5, max(inbreeding, na.rm = T))) %>% 
#     ungroup(.)
#   
#   peds_found <- bind_rows(parents[ , c("ID", "PARENT_FEMALE", "PARENT_MALE", "inbreeding", "gen2go")], next_gens)
#   #peds_found <- bind_rows(parents[ , c("ID", "PARENT_FEMALE", "PARENT_MALE", "inbreeding", "gen2go")], get_parents(GermID_list = parent_ids$ID[!(parent_ids %in% parents$ID)][1:10], # don't repeat if already in pedigree
#   #                                             n_gen = n_gen - 1))
#   #peds_found <- bind_rows(parents, parents_not_found, get_parents(GermID_list = parent_ids[!(parent_ids %in% c(parents$ID, parents_not_found$ID)) & !is.na(parent_ids)], # don't repeat if already in pedigree
#   #                                                                  n_gen = n_gen - 1))
#   
#   peds_sorted <- sort_pedigree(as.data.frame(peds_found))
#   
#   return(list(ped = peds_sorted, child_ids = child_ids, parent_ids = parent_ids)) # also returns df of pedigree-germplasm link for samples w/ data
# }

