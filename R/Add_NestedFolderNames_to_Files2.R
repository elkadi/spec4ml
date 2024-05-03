#'Add_NestedFolderNames_to_Files2
#'
#' Adding the folder names in the path to their files
#'
#' @param parent_folder Folder with nested subfolders containing the files to be renamed by adding their folder names in the path.
#' @param new_folder Folder with files with new names to be merged together.
#' @importFrom utils tail
#' @return non--All files in subfolders in the parent_folder will be renamed by adding their folder names in the path and moved to the new folder.
#' @export

Add_NestedFolderNames_to_Files2 <- function(parent_folder,new_folder) {
  folders <- list.dirs(parent_folder, recursive = TRUE, full.names = TRUE)
  parent_folder_length <- length(strsplit(parent_folder, "/")[[1]])
  for (folder in folders) {
    files <- list.files(folder, full.names = TRUE)
    folder_names <- unlist(strsplit(folder, "/"))
    relevant_names <- tail(folder_names, length(folder_names) - parent_folder_length + 1)
    for (file in files) {
      if (!file.info(file)$isdir) {
        new_filename <- paste0(paste(relevant_names, collapse = "_"), "_", basename(file))
        file.copy(file, file.path(new_folder, new_filename))
      }
    }
  }
}
