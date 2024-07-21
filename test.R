

# Preambles ---------------------------------------------------------------

# install.packages("devtools")
# install.packages("usethis")
# install.packages("rmarkdown")
# install.packages("fs")

library(devtools)
library(roxygen2)
library(knitr)
library(rmarkdown)
library(usethis)
# library(fs)


use_git_config(
  user.name = "lmaowisc",
  user.email = "lmao@biostat.wisc.edu"
)

# use_git()

# usethis:::use_devtools()


# Create, edit, and load code ---------------------------------------------
## create R file in "./R"
# use_r("mi2level.R")
# use_r("misc.R")
## Test your function in the new package
### devtools::load_all()
### Ctrl+Shift+L
load_all()



# Check package -----------------------------------------------------------
## 3 types of messages
## • ERRORs: Severe problems - always fix.
## • WARNINGs: Problems that you should fix, and must fix if you’re planning to
## submit to CRAN.
## • NOTEs: Mild problems or, in a few cases, just an observation.
## • When submitting to CRAN, try to eliminate all NOTEs.
check()



# Licenses ----------------------------------------------------------------

use_ccby_license()


# The DESCRIPTION file ----------------------------------------------------

# Type: Package
# Package: poset
# Title: Analysis of Partially Ordered Data
# Version: 1.0
# Author: Lu Mao
# Maintainer: Lu Mao <lmao@biostat.wisc.edu>
#   Description: Win ratio
# License: CC BY 4.0
# URL: https://sites.google.com/view/lmaowisc/
#   Depends:
#   R (>= 3.10)
# Suggests: knitr, rmarkdown
# VignetteBuilder:
#   knitr
# Config/testthat/edition: 3
# Encoding: UTF-8
# LazyData: true
# RoxygenNote: 7.3.1


# Commit changes to git ---------------------------------------------------

# Prerequisites:
# • GitHub account
# • create_github_token() - follow instructions
# • gitcreds::gitcreds_set() - paste PAT
# • git_sitrep() - verify
# use_github()
# - push content to new repository on GitHub

# usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)

edit_r_environ()
###### GitHub #######
create_github_token()
gitcreds::gitcreds_set()
git_sitrep()

use_github()
####################

# Documentation -----------------------------------------------------------

# RStudio: Code > Insert Roxygen Skeleton
## Ctrl+Alt+Shift+R
# Special comments (#') above function
#   definition in R/*.R

#### ----------- an example ---
#' Multiplicative win ratio (MWR) regression analysis
#'
#' @description Fit a multiplicative win ratio (MWR) regression model to
#' partially ordered outcome against covariates
#' @return An object of class \code{MWRreg} with the following components:
#' \item{beta}{A vector of estimated regression coefficients.}
#' @seealso \code{\link{wprod}}
#' @export
#' @importFrom utils combn
#' @importFrom stats complete.cases
#' @aliases MWRreg
#' @keywords MWRreg
#' @references Mao, L. (2024). Win ratio for partially ordered data.
#' Under revision.
#' @examples
#' set.seed(12345)
#' n <- 200
#' Z <- cbind(rnorm(n),rnorm(n))
#' \dontrun{
#'   use_git()
#' }
####

### Steps ###
# Go to function definition: Ctrl+.(type function name)
# • Cursor in function definition
# • Insert roxygen skeleton (Ctrl+Alt+Shift+R)
# • Complete the roxygen fields
# • document() (Ctrl+Shift+D) - create .rd files in ./man
# • ?myfunction

document()




# document() updates the NAMESPACE file with directives from Roxygen
# comments in your R code about export() and import()


# Package-level documentation ---------------------------------------------

use_package_doc()
#> ✔ Writing 'R/mypackage-package.R'
#> • Modify ‘R/mypackage-package.R’
document()
# ?poset

check() # again



# Install package to your library -----------------------------------------

install()




use_pkgdown_github_pages()
use_readme_rmd()
build_readme()

# Run once to configure your package to use pkgdown
usethis::use_pkgdown()
pkgdown::build_site()

usethis::use_pkgdown_github_pages()


# Work space  ------------------------------------------------------------

## data
bladder <- read.table("bladder.txt",sep="\t",header=T)

bladder$status[bladder$status != 0] <- 3 - bladder$status[bladder$status != 0]

use_data(bladder, overwrite = TRUE)

?rccf2


?bladder

use_r("data.R")

library(survival)

# build_manual(pkg = ".", path = NULL)
data("bladder")

id <- bladder$id
time <- bladder$time
status <- bladder$status
trt <- bladder$trt
trt <- ifelse(trt == 1,  "Thiotepa", "Placebo")

obj <- rccf(id,time,status,trt)

stat=obj$stat
S=obj$S
u1=obj$u1
u2=obj$u2
t=obj$t


# test on recurrent event
QLR=stat[1]/sqrt(S[1,1])
pLR=1-pchisq(QLR^2,1)

# test on death
QD=stat[2]/sqrt(S[2,2])
pD=1-pchisq(QD^2,1)

#joint test
Q=t(stat)%*%solve(S)%*%stat
p=1-pchisq(Q,2)

par(mfrow=c(1,2))

fatal=bladder[status != 1,]

obj=survfit(Surv(time,status > 0)~trt,data=fatal)
plot(obj,lty=c(1,3),frame=F,xlim=c(0,50),lwd=2,cex.lab=1.5,
     cex.axis=1.5,xlab="Time (months)",ylab="Survival probabilities")
# legend(0.82,2.8,c("Thiotepa","Placebo"),lty=c(1,3),cex=1.5,lwd=2)
text(35,0.14,paste0("Log-rank test: p=",round(pD,3)),cex=1.2)

#plot the mean functions
stepmu1=stepfun(t,c(0,u1))
stepmu2=stepfun(t,c(0,u2))
plot(stepmu2,do.points=F,xlab="Time (months)",ylab="Cumulative tumor frequency",
     xlim=c(0,60),ylim=c(0,3), main="", frame=F,lwd=2,cex.lab=1.5,
     cex.axis=1.5,lty=3)
plot(stepmu1,add=T,lty=1,do.points=F,lwd=2)
# legend(1,2.8,c("Thiotepa","Placebo"),lty=c(1,3),cex=1.5,lwd=2)
text(40,0.4,paste0("Ghosh-Lin test: p=",round(pLR,3)),cex=1.2)


### VIGNETTES #######


# use_article("Two-level-mi")
use_vignette("bladder-tumor")


