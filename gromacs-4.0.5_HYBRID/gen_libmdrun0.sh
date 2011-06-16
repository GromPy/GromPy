#! /bin/bash

KERNELDIR=/home/rpool/src/gromacs-4.0.5_TEST/src/kernel
LIBMDRUNDIR=$KERNELDIR/libmdrun0
LIBGROMACSDIR=/home/rpool/src/gromacs-4.0.5_TEST/lib

if [ ! -e $LIBMDRUNDIR ]; then
  mkdir $LIBMDRUNDIR
fi

# in libmdrun0/ (na compilen, extracting objects from peviously compiled dlls):
cd $KERNELDIR
cp compute_io.o glaasje.o gctio.o ionize.o do_gct.o repl_ex.o xutils.o md.o mdrun.o genalg.o $LIBMDRUNDIR
cd $LIBMDRUNDIR
cp $LIBGROMACSDIR/libgmx.a ./
cp $LIBGROMACSDIR/libmd.a ./
ar -x libmd.a
ar -x libgmx.a
gcc -shared -Wl,-soname,mdrun.so.5 -o mdrun.so.5.0.0 *.o -lnsl -lfftw3f -lm 
ln -sf mdrun.so.5.0.0 mdrun.so.5
ln -sf mdrun.so.5 mdrun.so

exit
