#' Add_FolderName_to_Files
#'
#' Adding the folder names to their files
#'
#' @param parent_folder Folder with subfolders containing the files to be copied.
#' @param new_folder  the Folder to contain the copied files.
#' @param recurse to specify if it recurse or not (default recurse=FALSE)
#' @return All the files to be copied in the new_folder
#' @import fs
#' @export


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
