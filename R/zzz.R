######################################################################
#
# zzz.R
#
# Written by John Hughes <hughesj@umn.edu>.
#
# Last Modified 05/31/12
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/ngspatial package
#
# .onLoad is run when the package is loaded with library(ngspatial)
#
######################################################################

#' @import utils

.onLoad = function(libname, pkgname)
{
    # library.dynam("ngspatial", package = pkgname, lib.loc = libname)
    temp = packageDescription("ngspatial")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2012, John Hughes, University of Minnesota\n",
#"                    Xiaohui Cui, University of Minnesota\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("ngspatial").\n')
    msg = paste(msg, 'Type help("ngspatial-package") to get started.\n')
    packageStartupMessage(msg)
}
