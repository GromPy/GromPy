#! /bin/bash

KERNELDIR=/home/rpool/src/gromacs-4.0.5_TEST/src/kernel
LIBMDRUNDIR=$KERNELDIR/libmdrun1
LIBGROMACSDIR=/home/rpool/src/gromacs-4.0.5_TEST/lib
INCLUDEGROMACSDIR=/home/rpool/src/gromacs-4.0.5_TEST/include

if [ ! -e $LIBMDRUNDIR ]; then
  mkdir $LIBMDRUNDIR
fi

# in libmdrun1/ (na compilen, linking with dlls):
cd $KERNELDIR
cp compute_io.o glaasje.o gctio.o ionize.o do_gct.o repl_ex.o xutils.o md.o mdrun.o genalg.o $LIBMDRUNDIR
cd $LIBMDRUNDIR
gcc -shared -Wl,-soname,mdrun.so.5 -o mdrun.so.5.0.0 *.o $LIBGROMACSDIR/libgmx.so $LIBGROMACSDIR/libmd.so -lnsl -lfftw3f -lm 
ln -sf mdrun.so.5.0.0 mdrun.so.5
ln -sf mdrun.so.5 mdrun.so

exit
