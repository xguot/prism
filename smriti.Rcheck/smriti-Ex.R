pkgname <- "smriti"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "smriti-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('smriti')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("smriti_forest")
### * smriti_forest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: smriti_forest
### Title: missForest-Smriti Refinement Wrapper
### Aliases: smriti_forest

### ** Examples

## Not run: 
##D df <- data.frame(T1 = c(1, NA, 3, 4), T2 = c(NA, 2, 3, 4))
##D smriti_forest(df, time_cols = 1:2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("smriti_forest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("smriti_impute")
### * smriti_impute

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: smriti_impute
### Title: Smriti Automated Longitudinal Imputation
### Aliases: smriti_impute

### ** Examples

# Simulated longitudinal data with scattered missingness
df <- data.frame(
  T1 = c(1.2, NA, 2.8, 3.1),
  T2 = c(2.1, 2.5, NA, 4.0),
  T3 = c(3.0, 3.3, 4.1, NA)
)
smriti_impute(df, time_cols = 1:3, robust = FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("smriti_impute", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("smriti_mi")
### * smriti_mi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: smriti_mi
### Title: Smriti Multiple Imputation Wrapper
### Aliases: smriti_mi

### ** Examples

## Not run: 
##D df <- data.frame(
##D   T1 = c(1, 2, NA, 4),
##D   T2 = c(NA, 2, 3, 4),
##D   T3 = c(1, NA, 3, 4)
##D )
##D mi_list <- smriti_mi(df, time_cols = c("T1", "T2", "T3"), m = 5)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("smriti_mi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("smriti_mice")
### * smriti_mice

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: smriti_mice
### Title: mice-Smriti Refinement Wrapper
### Aliases: smriti_mice

### ** Examples

## Not run: 
##D df <- data.frame(T1 = c(1, NA, 3, 4), T2 = c(NA, 2, 3, 4))
##D smriti_mice(df, time_cols = 1:2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("smriti_mice", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("smriti_ranger")
### * smriti_ranger

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: smriti_ranger
### Title: missRanger-Smriti Refinement Wrapper
### Aliases: smriti_ranger

### ** Examples

## Not run: 
##D df <- data.frame(T1 = c(1, NA, 3, 4), T2 = c(NA, 2, 3, 4))
##D smriti_ranger(df, time_cols = 1:2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("smriti_ranger", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
