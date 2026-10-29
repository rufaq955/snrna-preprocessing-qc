#!/usr/bin/env Rscript
#
# Install the R dependencies required by steps 04, 05 and 06b.
#
# DoubletFinder and SeuratDisk are not on CRAN; they come from GitHub.
# SeuratDisk in particular needs a specific remotes install and can fail
# against a Seurat v5 object - see docs/troubleshooting.md.
#
# Usage: Rscript scripts/install_r_deps.R
#
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c("Seurat", "optparse", "data.table", "hdf5r", "remotes")

missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                                 quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing) > 0) {
  message("Installing from CRAN: ", paste(missing, collapse = ", "))
  install.packages(missing)
} else {
  message("All CRAN dependencies already present.")
}

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  message("Installing DoubletFinder from GitHub (chris-mcginnis-ucsf/DoubletFinder)")
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
}

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  message("Installing SeuratDisk from GitHub (mojaveazure/seurat-disk)")
  remotes::install_github("mojaveazure/seurat-disk")
}

message("\nR dependencies installed. Versions:")
for (pkg in c("Seurat", "DoubletFinder", "SeuratDisk", "optparse")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  %-14s %s", pkg, as.character(packageVersion(pkg))))
  } else {
    message(sprintf("  %-14s NOT INSTALLED", pkg))
  }
}
