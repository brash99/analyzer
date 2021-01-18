#!/bin/csh
#

# echo "[setupCODA]: Setting up CODA 2.6.2..."

# CODA 2.6.2 setup
# for now, we use /site
setenv CODA  /site/coda/2.6.2
source $CODA/.setup

setenv MSQL_TCP_HOST sbs2
setenv MSQL_HOME /home/adaq/rcdb
# SBS fastbus (FB)
setenv EXPID cdetfb
setenv SESSION cdetfb1
setenv KILL_CODA_SCRIPT /home/adaq/bin/kill_coda_script
 
# CODA 3 Run Control
setenv AFECS_HOME ${CODA}/afecs
setenv CMSG_HOME  ${CODA}/cMsg
setenv COOL_HOME  /home/adaq/coda/cool
setenv RCREMLOG   /home/adaq/coda/

setenv CODA_DB_DRIVER com.imaginary.sql.msql.MsqlDriver
setenv CODA_DB_URL jdbc:msql://${MSQL_TCP_HOST}:8101/$EXPID
#setenv CODA_DB_USER 
#setenv CODA_DB_PASSWORD

setenv JAVA_HOME /u/apps/java/jdk1.6.0_03
setenv PATH ${AFECS_HOME}/bin:${JAVA_HOME}/bin:~/bin:$PATH

setenv ROOTSYS /u/apps/root/PRO/root
#setenv ROOTSYS /u/apps/root/5.34.02/root
setenv SCONS ${HOME}/scons
setenv PPC_CROSS_COMPILER /home/adaq/vxworks/5.5/gnu_vxworks
setenv PATH ${SCONS}/bin:${ROOTSYS}/bin:${PPC_CROSS_COMPILER}/bin:${PATH}
#setenv LD_LIBRARY_PATH ${SCONS}/lib:${ROOTSYS}/lib:${PPC_CROSS_COMPILER}/lib:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH ${SCONS}/lib:${ROOTSYS}/lib:${PPC_CROSS_COMPILER}/lib

setenv ANALYZER /home/adaq/analyzer_updates/analyzer_jan132021
setenv REPLAY /home/adaq/analyzer_updates/analyzer_jan132021/replay
setenv DB_DIR $ANALYZER/DB
setenv DATA_DIR $ANALYZER/../data
setenv OUT_DIR $REPLAY/rootfiles

setenv PATH ${ANALYZER}:${ANALYZER}/apps:${ANALYZER}/bin:${PATH}
setenv LD_LIBRARY_PATH ${ANALYZER}/Podd:${ANALYZER}/SDK:${LD_LIBRARY_PATH}

# echo "[setupCODA]: --> Done."
