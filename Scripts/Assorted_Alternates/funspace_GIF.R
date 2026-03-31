
library("purrr")
library("magick")

setwd("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Figures/funspace_GIF")

gifname <- "funspace_1myr_Timeslice"

rev(list.files(path = getwd(), pattern = "*.jpg", full.names = T)) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=5) %>% # animates, can set the Frames Per Second or number of loops
  image_write(paste0(getwd(),"/", gifname,".gif")) # write to current dir
