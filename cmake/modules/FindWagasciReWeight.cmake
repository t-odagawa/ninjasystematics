# This module defines
# WAGASCI_REWEIGHT_LIBRARY, the name of the library to link against
# WAGASCI_REWEIGHT_FOUND, if false, do not try to link to WAGASCI_REWEIGHT
# WAGASCI_REWEIGHT_INCLUDE_DIR, where to find WAGASCI_REWEIGHT headers
#
# Note: If you see an empty WAGASCI_REWEIGHT_LIBRARY_TEMP in your configuration
# and no WAGASCI_REWEIGHT_LIBRARY, it means CMake did not find your WAGASCI_REWEIGHT library
# (libWAGASCI_REWEIGHT.dylib, libWAGASCI_REWEIGHT.so, etc).
# These values are used to generate the final WAGASCI_REWEIGHT_LIBRARY variable,
# but when these values are unset, WAGASCI_REWEIGHT_LIBRARY does not get created.
#
#
# $WAGASCI_REWEIGHTDIR is an environment variable that would
# correspond to the CMAKE_INSTALL_PREFIX=$WAGASCI_REWEIGHTDIR
# used in building WAGASCI_REWEIGHT.

SET(WAGASCI_REWEIGHT_SEARCH_PATHS
        /usr/local
        /opt/local
        ${WAGASCI_REWEIGHT_PATH}
        )

FIND_PATH(WAGASCI_REWEIGHT_INCLUDE_DIR WeightTree.hpp
        HINTS
        $ENV{WAGASCI_REWEIGHTDIR}
        /usr/local
        /opt/local
        PATH_SUFFIXES include/wagasci/reweight include
        PATHS ${WAGASCI_REWEIGHT_SEARCH_PATHS}
        )

FIND_LIBRARY(WAGASCI_REWEIGHT_LIBRARY_TEMP
        NAMES WagasciReWeight
        HINTS
        $ENV{WAGASCI_REWEIGHTDIR}
        PATH_SUFFIXES lib64/wagasci/reweight lib/wagasci/reweight lib lib64
        PATHS ${WAGASCI_REWEIGHT_SEARCH_PATHS}
        )

IF(WAGASCI_REWEIGHT_LIBRARY_TEMP)
    # Set the final string here so the GUI reflects the final state.
    SET(WAGASCI_REWEIGHT_LIBRARY ${WAGASCI_REWEIGHT_LIBRARY_TEMP} CACHE STRING "Where the WAGASCI_REWEIGHT Library can be found")
    # Set the temp variable to INTERNAL so it is not seen in the CMake GUI
    SET(WAGASCI_REWEIGHT_LIBRARY_TEMP "${WAGASCI_REWEIGHT_LIBRARY_TEMP}" CACHE INTERNAL "")
ENDIF(WAGASCI_REWEIGHT_LIBRARY_TEMP)

INCLUDE(FindPackageHandleStandardArgs)

#FIND_PACKAGE_HANDLE_STANDARD_ARGS(WAGASCI_REWEIGHT REQUIRED_VARS WAGASCI_REWEIGHT_LIBRARY WAGASCI_REWEIGHT_INCLUDE_DIR)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(WagasciReWeight REQUIRED_VARS WAGASCI_REWEIGHT_LIBRARY WAGASCI_REWEIGHT_INCLUDE_DIR)

