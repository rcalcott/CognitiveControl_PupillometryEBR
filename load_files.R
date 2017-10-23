#load a bunch of files based on a pattern,save all into a big dataframe
load_files <- function(path,pattern,sep='\t'){
  
  file.glob <- glob2rx(pattern)
  
  file_list <- list.files(path=path,pattern=file.glob,full.names=TRUE)
  print('Found files:')
  print(file_list)
  
  readf <- function(x) readr::read_delim(x,sep)
  
  alldfs <- lapply(file_list,readf)
  dataset <- dplyr::bind_rows(alldfs)
  
  
  return(dataset)
}


