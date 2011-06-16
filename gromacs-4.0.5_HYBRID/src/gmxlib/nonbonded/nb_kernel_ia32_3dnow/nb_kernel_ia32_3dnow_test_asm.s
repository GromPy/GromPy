## $Id: nb_kernel_ia32_3dnow_test_asm.s,v 1.3 2006/04/19 09:04:53 lindahl Exp $
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##


.globl nb_kernel_ia32_3dnow_test_asm
.globl _nb_kernel_ia32_3dnow_test_asm
nb_kernel_ia32_3dnow_test_asm: 
_nb_kernel_ia32_3dnow_test_asm: 
        femms
        pswapd %mm0,%mm0        ## Extended 3DNow instr.
        pxor   %mm1,%mm1
        pf2iw %mm1,%mm0         ## Extended 3DNow instr.
        femms
        ret




