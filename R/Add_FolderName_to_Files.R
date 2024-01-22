#' Add_FolderName_to_Files
#'
#' Adding the folder names to their files
#'
#' @param parent_folder Folder with subfolders containing the files to be renamed by adding their folder name.
#' @return All files in subfolders in the parent_folder renamed by adding their folder name.
#' @examples
#' Add_FolderName_to_Files(parent_folder)


Add_FolderName_to_Files <- function(parent_folder) {
  folders <- list.dirs(parent_folder, recursive = FALSE, full.names = TRUE)
  for (folder in folders) {
    files <- list.files(folder, full.names = TRUE)
    folder_name <- basename(folder)
    for (file in files) {
      if (!file.info(file)$isdir) {
        new_filename <- paste0(folder_name, "_", basename(file))
        file.rename(file, file.path(dirname(file), new_filename))
      }
    }
  }
}
