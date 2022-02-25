cap noi ado uninstall multe
mata: mata clear
mata: mata set matastrict on
mata: mata set mataoptimize on
cap noi mkdir build
cap noi erase build/lmanyiv.mlib
do src/mata/multe_helpers.mata
* do src/mata/multe_decomposition.mata
do src/mata/multe.mata
mata: mata mlib create lmulte, dir("build") replace
mata: mata mlib add lmulte MulTE*() multe_helper*(), dir("build") complete
net install multe, from(`c(pwd)') replace
