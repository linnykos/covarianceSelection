.is.listOfMatrix <- function(obj, obj.name){
  if(!is.list(obj)) stop(paste(obj.name, "is not a list"))
  
  vec <- sapply(obj, is.matrix)
  if(!all(vec)) stop(paste(obj.name, "has elements in list not matrices"))
  
  TRUE
}

.is.listofNumeric <- function(obj, obj.name){
  if(!is.list(obj)) stop(paste(obj.name, "is not a list"))
  
  vec <- sapply(obj, is.numeric)
  if(!all(vec)) stop(paste(obj.name, "has elements in list not numeric"))
  
  vec <- sapply(obj, is.vector)
  if(!all(vec)) stop(paste(obj.name, "has elements in list not vector"))
  
  TRUE
}

.is.nonNegInteger <- function(vec, name){
  if(!all(!is.na(vec))) stop(paste(name,
                                   "cannot contain NA and non-NA elements"))
  if(!is.numeric(vec)) stop(paste(name, "must be a numeric vector"))
  if(any(vec %% 1 != 0)) stop(paste(name, "must be a vector of integers"))
  if(any(vec < 0)) stop(paste(name, "must be a vector of positive integers"))
  if(any(duplicated(vec))) stop(paste(name, "must contain all unique",
                                      "elements"))
  
  TRUE
}