
rm(list = ls())
Rcpp::compileAttributes(pkgdir="/Users/zhiz/Downloads/psbcSpeedUp/psbcSpeedUp/")
devtools::document("/Users/zhiz/Downloads/psbcSpeedUp/psbcSpeedUp")
devtools::build("/Users/zhiz/Downloads/psbcSpeedUp/psbcSpeedUp", vignettes=T)
install.packages("/Users/zhiz/Downloads/psbcSpeedUp/psbcSpeedUp_2.0.5.tar.gz",repos = NULL,type = "source",build_vignettes = T)
