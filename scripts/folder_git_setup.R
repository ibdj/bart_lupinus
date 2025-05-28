#### set up of standard folders ####

# Define the folder names
folders <- c("scripts", "data", "outputs")

# Loop through each folder name
for (folder in folders) {
  # Check if the folder exists, and create it if it doesn't
  if (!dir.exists(folder)) {
    dir.create(folder)
    message(paste("Created folder:", folder))
  } else {
    message(paste("Folder already exists:", folder))
  }
}

#### install packages ####
install.packages("pacman")
pacman::p_load(usethis, tidyverse) 

#### setting up git and github ####
use_git_config(
  user.name = "ibdj", 
  user.email = "ibdjacobsen@gmail.com"
)

# see tokens https://github.com/settings/tokens
usethis::create_github_token()
gitcreds::gitcreds_set()

git_vaccinate() 

usethis::use_git()

use_github()

git_default_branch_rename()



