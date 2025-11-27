options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("rmarkdown")
install.packages("knitr")
rmarkdown::render("Projet.Rmd")
