/* -*- mode: c; tab-width: 4; indent-tabs-mode: n; c-basic-offset: 4 -*- 
 *
 * $Id: nb_kernel314_ppc_altivec.h,v 1.1 2004/12/26 19:26:00 lindahl Exp $
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _NB_KERNEL314_ALTIVEC_H_
#define _NB_KERNEL314_ALTIVEC_H_

/*! \file  nb_kernel314_ppc_altivec.h
 *  \brief Altivec-optimized versions of nonbonded kernel 314
 *
 *  \internal
 */



/*! \brief Nonbonded kernel 314 with forces, optimized for Altivec.
 *
 *  \internal
 *
 *  <b>Coulomb interaction:</b> Tabulated <br>
 *  <b>VdW interaction:</b> Lennard-Jones <br>
 *  <b>Water optimization:</b> Pairs of TIP4P waters interaction <br>
 *  <b>Forces calculated:</b> Yes <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel.h
 */
void
nb_kernel314_ppc_altivec  (int * restrict  p_nri,        int * restrict    iinr,   int * restrict    jindex,
						   int  * restrict   jjnr,     int   * restrict  shift,  float * restrict  shiftvec,
						   float * restrict  fshift,   int   * restrict  gid,    float * restrict  pos,
						   float * restrict  faction,  float  * restrict charge, float * restrict p_facel,
						   float * restrict p_krf,        float * restrict p_crf,      float * restrict  Vc,
						   int  * restrict   type,     int * restrict  p_ntype,    float * restrict  vdwparam,
						   float * restrict  Vvdw,     float * restrict p_tabscale, float * restrict  VFtab,
						   float * restrict  invsqrta, float * restrict  dvda,   float * restrict p_gbtabscale,
						   float  * restrict GBtab,    int * restrict  nthreads, int *  restrict count,
						   void *  mtx,        int * restrict  outeriter,int * restrict  inneriter,
						   float * work);



/*! \brief Nonbonded kernel 314 without forces, optimized for Altivec.
 *
 *  \internal
 *
 *  <b>Coulomb interaction:</b> Tabulated <br>
 *  <b>VdW interaction:</b> Lennard-Jones <br>
 *  <b>Water optimization:</b> Pairs of TIP4P waters interaction <br>
 *  <b>Forces calculated:</b> No <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel.h
 */
void
nb_kernel314nf_ppc_altivec(int * restrict  p_nri,        int * restrict    iinr,   int * restrict    jindex,
						   int  * restrict   jjnr,     int   * restrict  shift,  float * restrict  shiftvec,
						   float * restrict  fshift,   int   * restrict  gid,    float * restrict  pos,
						   float * restrict  faction,  float  * restrict charge, float * restrict p_facel,
						   float * restrict p_krf,        float * restrict p_crf,      float * restrict  Vc,
						   int  * restrict   type,     int * restrict  p_ntype,    float * restrict  vdwparam,
						   float * restrict  Vvdw,     float * restrict p_tabscale, float * restrict  VFtab,
						   float * restrict  invsqrta, float * restrict  dvda,   float * restrict p_gbtabscale,
						   float  * restrict GBtab,    int * restrict  nthreads, int *  restrict count,
						   void *  mtx,        int * restrict  outeriter,int * restrict  inneriter,
						   float * work);





#endif /* _NB_KERNEL314_ALTIVEC_H_ */
