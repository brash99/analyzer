###### Hall A SDK Main SConscript File #####

### Modify 'libname' and the list of source files 'src' as needed

from podd_util import build_library
Import ('baseenv')

# Library name (lib<libname>.so will be built)
libname = 'SBS'

# Sources (each must have a corresponding .h header file)
# This must be a space-separated string, not a Python list
src = """
FbusDetector.cxx FbusApparatus.cxx
HcalDetector.cxx HcalApparatus.cxx
CdetDetector.cxx CdetApparatus.cxx
"""

# Tell SCons to build this library from these sources.
# A ROOT dictionary, defined in <libname>_LinkDef.h, will be built as well.
# For more info, see the documention in $ANALYZER/site_scons/podd_util.py
build_library(baseenv, libname, src)
