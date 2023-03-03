## From the attribute column, extract the classification as is
extract_raw_classification <- function(attribute_string) {
  classification <- substring(attribute_string, regexpr("*Classification=", attribute_string),regexpr(";Se", attribute_string)-1)
  classification <- unlist(str_split(classification,pattern='='))[2]
  return(classification)
}

## From the attribute column, extract the classification as is
##     depending on value, reclassify into 10 larger groups to
##     reduce the number of classifications + obtain broader classifications
extract_classification <- function(attribute_string) {
  classification <- substring(attribute_string, regexpr("*Classification=", attribute_string),regexpr(";Se", attribute_string)-1)
  classification <- unlist(str_split(classification,pattern='='))[2]
  if (grepl('LINE',classification)) {
    upd_classification <- 'LINE'
  } else if (grepl('LTR',classification)) {
    upd_classification <- 'LTR' 
  } else if (grepl('Helitron',classification)) {
    upd_classification <- 'Helitron' 
  } else if (grepl('MITE',classification)) {
    upd_classification <- 'MITE' 
  } else if (grepl('DNA',classification)) {
    upd_classification <- 'TIR' 
  } else if (grepl('knob',classification)) {
    upd_classification <- 'knob' 
  } else if (grepl('Simple',classification)) {
    upd_classification <- 'Simple Repeat' 
  } else if (grepl('Cent',classification)) {
    upd_classification <- 'Cent/CentC' 
  } else if (grepl('rDNA',classification)) {
    upd_classification <- 'rDNA/spacer' 
  } else if (grepl('Low',classification)) {
    upd_classification <- 'Low Complexity' 
  } else (
    upd_classification <- "subtelomere/4-12-1"
  )
  return(upd_classification)
}

## From the attribute column, extract the family as is
extract_family <- function(attribute_string) {
  family <- substring(attribute_string, regexpr("*Name=", attribute_string),regexpr(";C", attribute_string)-1)
  family <- unlist(str_split(family,pattern='='))[2]
  return(family)
}

## From the attribute column, extract the method as is
##    to see if the call was made via structure vs homology
extract_method <- function(attribute_string) {
  method <- substring(attribute_string, regexpr("*Method=", attribute_string))
  method <- unlist(str_split(method,pattern='='))[2]
  
  if (nchar(method) > 10) {
    new_method <- unlist(str_split(method,pattern=';'))[1]
  } else {
    new_method <- method
  }
  return(new_method)
}

## From the attribute column, extract the identity as is
extract_identity <- function(attribute_string) {
  identity <- substring(attribute_string, regexpr("*dentity=", attribute_string),regexpr(";M", attribute_string)-1)
  identity <- unlist(str_split(identity,pattern='='))[2]
  return(identity)
}