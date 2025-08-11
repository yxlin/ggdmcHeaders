q(save = "no")
rm(list = ls())


tools::showNonASCIIfile("DESCRIPTION")
devtools::check()
usethis::use_gpl_license(version = 2) # or use_gpl_license(version = 2)
devtools::check(manual = TRUE, cran = TRUE) # Generates a PDF manual if needed
devtools::build()

# install.packages("gitcreds")
library(gitcreds)

# Set the GitHub token
gitcreds_set()

gitcreds_get()
gitcreds_clear()

usethis::gh_token_help()

library(gitcreds)
gitcreds_set(url = "https://github.com/yxlin/ggdmcHeaders")

usethis::create_github_token(
    scopes = c("repo", "workflow"),
    description = "GitHub token for ggdmcHeaders package"
)


rhub::rhub_setup()
rhub::rhub_platforms()

rhub::rhub_doctor()

clang20_R_devel <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmcHeaders",
    platforms = "clang20"
)

windows_R_devel <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmcHeaders",
    platforms = "windows"
)

macos_arm64 <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmcHeaders",
    platforms = "macos-arm64"
)
