#! /bin/bash

EXPECTED_N_ARGS=1
if [ $# -ne $EXPECTED_N_ARGS ]; then
    echo "Input the correct amount of arguments!"
    echo "Usage:"
    echo "------"
    echo "source ./SourceGromPyEnv.sh \$GMXINSTALLDIR"
    echo "EXITING..."
else
    GMXINSTALLDIR=$1

    source $GMXINSTALLDIR/bin/GMXRC

    LIBCPATH=`ldd $GMXINSTALLDIR/bin/mdrun | grep libc.so.6 | awk '{print $3}'`
    GROMPY_LIBGMX=$GMXINSTALLDIR/lib/libgmx_grompy.so
    GROMPY_LIBMD=$GMXINSTALLDIR/lib/libmd_grompy.so
    GROMPY_LIBMDRUN=$GMXINSTALLDIR/lib/libmdrun_grompy.so

    export LIBCPATH
    export GROMPY_LIBGMX
    export GROMPY_LIBMD
    export GROMPY_LIBMDRUN
fi

