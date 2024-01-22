copy_from_subfolders <- function(parent_folder, new_folder, recurse=FALSE) {
     subfolders <- dir_info(parent_folder, recurse = recurse, type = "directory")$path
     
       for (subfolder in subfolders) {
           files <- dir_info(subfolder, type = "file")$path
           
             for (file in files) {
                 file_name <- basename(file)
                 destination <- file.path(new_folder, file_name)
                 file_copy(file, destination)
               }
         }
  }



parent_folder<-"E:/(A) Cell Free-NIR Measurements/NIR Measurments - Avantes/Gellan Gum - SPD"
new_folder<-"E:/(A) Cell Free-NIR Measurements/NIR Measurments - Avantes/SpdNew"
copy_from_subfolders(parent_folder,new_folder)