#!/bin/sh

PACKAGE_NAME=$(sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PACKAGE_VERSION=$(sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PACKAGE_FILE=${PACKAGE_NAME}_${PACKAGE_VERSION}.tar.gz

R -q -e 'library("roxygen2"); roxygen2::roxygenise()'

R CMD build .
R CMD check "$PACKAGE_FILE" --no-manual
R CMD INSTALL "$PACKAGE_FILE"