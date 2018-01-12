
.onAttach = function(libname, pkgname)
{
    temp = packageDescription("ngspatial")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2013-2018, John Hughes, University of Colorado Denver\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("ngspatial").\n', sep = "")
    msg = paste(msg, 'Type help(package = ngspatial) to get started.\n', sep = "")
    packageStartupMessage(msg)
}
