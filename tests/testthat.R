# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(pairedGSEA)
data("example_se", "example_ora_results", "example_diff_result", package = "pairedGSEA")
test_check("pairedGSEA")
