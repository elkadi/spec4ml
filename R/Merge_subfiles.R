#' Copying files from sub-folders in one folder
#'
#' @param parent_folder Folder with subfolders containing the subfiles to be merged.
#' @param new_folder Folder for the subfiles to be copied to.
#' @return All files copied from subfolders in the parent_folder to the new_folder.
#' @import fs
#' @examples
#' Merge_subfiles(parent_folder, new_folder, recurse=FALSE)


Merge_subfiles <- function(parent_folder, new_folder, recurse=FALSE) {
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

