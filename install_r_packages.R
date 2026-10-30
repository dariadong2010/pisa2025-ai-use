# Installs the R packages used by the analysis scripts.
pkgs <- c("nanoparquet", "poLCA", "mice", "multilevLCA")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
for (p in pkgs) cat(p, as.character(packageVersion(p)), "\n")
