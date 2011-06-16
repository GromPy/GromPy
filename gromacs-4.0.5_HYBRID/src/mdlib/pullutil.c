/*
 * $Id: pullutil.c,v 1.32.2.1 2009/01/24 17:42:13 hess Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include "sysstuff.h"
#include "princ.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "network.h"
#include "pbc.h"
#include "pull.h"

void pull_d_pbc_dx(int npbcdim,matrix box,
		   const dvec x1, const dvec x2, dvec dx)
{
  int m,d;
  
  /* Only correct for atom pairs with a distance within
   * half of the smallest diagonal element of box.
   */
  dvec_sub(x1,x2,dx);
  for(m=npbcdim-1; m>=0; m--) {
    while (dx[m] < -0.5*box[m][m]) {
      for(d=0; d<DIM; d++)
	dx[d] += box[m][d];
    }
    while (dx[m] >=  0.5*box[m][m]) {
      for(d=0; d<DIM; d++)
	dx[d] -= box[m][d];
    }
  }
}

static void pull_set_pbcatom(t_commrec *cr, t_pullgrp *pg,
			     t_mdatoms *md, rvec *x,
			     rvec x_pbc)
{
  int a,m;

  if (cr && PAR(cr)) {
    if (DOMAINDECOMP(cr)) {
      if (cr->dd->ga2la[pg->pbcatom].cell == 0)
	a = cr->dd->ga2la[pg->pbcatom].a;
      else
	a = -1;
    } else {
      a = pg->pbcatom;
    }
    
    if (a >= md->start && a < md->start+md->homenr) {
      copy_rvec(x[a],x_pbc);
    } else {
      clear_rvec(x_pbc);
    }
  } else {
    copy_rvec(x[pg->pbcatom],x_pbc);
  }
}

static void pull_set_pbcatoms(t_commrec *cr, t_pull *pull,
			      t_mdatoms *md, rvec *x,
			      rvec *x_pbc)
{
  int g,n,m;

  n = 0;
  for(g=0; g<1+pull->ngrp; g++) {
    if ((g==0 && PULL_CYL(pull)) || pull->grp[g].pbcatom == -1) {
      clear_rvec(x_pbc[g]);
    } else {
      pull_set_pbcatom(cr,&pull->grp[g],md,x,x_pbc[g]);
      for(m=0; m<DIM; m++) {
	if (pull->dim[m] == 0) {
	  x_pbc[g][m] = 0.0;
	}
      }
      n++;
    }
  }
  
  if (cr && PAR(cr) && n > 0) {
    /* Sum over the nodes to get x_pbc from the home node of pbcatom */
    gmx_sum((1+pull->ngrp)*DIM,x_pbc[0],cr);
  }
}

/* switch function, x between r and w */
static real get_weight(real x, real r1, real r0)
{
  real weight; 

  if (x >= r0)
    weight = 0;
  else if (x <= r1)
    weight = 1;
  else
    weight = (r0 - x)/(r0 - r1);

  return weight;
}

static void make_cyl_refgrps(t_commrec *cr,t_pull *pull,t_mdatoms *md,
			     t_pbc *pbc,rvec *x,rvec *xp) 
{
  static double *dbuf=NULL;
  int g,i,ii,m,start,end;
  rvec g_x,dx;
  double r0_2,sum_z,sum_zp,dr2,mass,weight,wmass,wwmass;
  t_pullgrp *pref,*pdyna;
  gmx_ga2la_t *ga2la=NULL;

  if (dbuf == NULL) {
    snew(dbuf,pull->ngrp*4);
  }

  if (cr && DOMAINDECOMP(cr))
    ga2la = cr->dd->ga2la;
  
  start = md->start;
  end   = md->homenr + start;

  r0_2 = dsqr(pull->cyl_r0);

  /* loop over all groups to make a reference group for each*/
  pref = &pull->grp[0];
  for(g=1; g<1+pull->ngrp; g++) {
    pdyna = &pull->dyna[g];
    sum_z = 0;
    sum_zp = 0;
    wmass = 0;
    wwmass = 0;
    pdyna->nat_loc = 0;

    for(m=0; m<DIM; m++)
      g_x[m] = pull->grp[g].x[m];

    /* loop over all atoms in the main ref group */
    for(i=0; i<pref->nat; i++) {
      ii = pull->grp[0].ind[i];
      if (ga2la) {
	if (ga2la[pref->ind[i]].cell == 0)
	  ii = ga2la[pref->ind[i]].a;
	else
	  ii = -1;
      }
      if (ii >= start && ii < end) {
	/* get_distance takes pbc into account */
	pbc_dx_aiuc(pbc,x[ii],g_x,dx);
	dr2 = dx[XX]*dx[XX] + dx[YY]*dx[YY];

	if (dr2 < r0_2) {
	  /* add to index, to sum of COM, to weight array */
	  if (pdyna->nat_loc >= pdyna->nalloc_loc) {
	    pdyna->nalloc_loc = over_alloc_large(pdyna->nat_loc+1);
	    srenew(pdyna->ind_loc,pdyna->nalloc_loc);
	    srenew(pdyna->weight_loc,pdyna->nalloc_loc);
	  }
	  pdyna->ind_loc[pdyna->nat_loc] = ii;
	  mass = md->massT[ii];
	  weight = get_weight(sqrt(dr2),pull->cyl_r1,pull->cyl_r0);
	  pdyna->weight_loc[pdyna->nat_loc] = weight;
	  sum_z += mass*weight*x[ii][ZZ];
	  if (xp) {
	    sum_zp += mass*weight*xp[ii][ZZ];
	  }
	  wmass += mass*weight;
	  wwmass += mass*sqr(weight);
	  pdyna->nat_loc++;
	}
      }
    }
    dbuf[(g-1)*4+0] = wmass;
    dbuf[(g-1)*4+1] = wwmass;
    dbuf[(g-1)*4+2] = sum_z;
    dbuf[(g-1)*4+3] = sum_zp;
  }

  if (cr && PAR(cr)) {
    /* Sum the contributions over the nodes */
    gmx_sumd(pull->ngrp*4,dbuf,cr);
  }

  for(g=1; g<1+pull->ngrp; g++) {
    pdyna = &pull->dyna[g];

    wmass        = dbuf[(g-1)*4+0];
    wwmass       = dbuf[(g-1)*4+1];
    pdyna->wscale = wmass/wwmass;
    pdyna->invtm = 1.0/(pdyna->wscale*wmass);

    pdyna->x[XX] = 0;
    pdyna->x[YY] = 0;
    pdyna->x[ZZ] = dbuf[(g-1)*4+2]/wmass;
    if (xp) {
      pdyna->xp[XX] = 0;
      pdyna->xp[YY] = 0;
      pdyna->xp[ZZ] = dbuf[(g-1)*4+3]/wmass;
    }

    if (debug) {
      fprintf(debug,"Pull cylinder group %d:%8.3f%8.3f%8.3f m:%8.3f\n",
              g,pdyna->x[0],pdyna->x[1],
              pdyna->x[2],1.0/pdyna->invtm);
    }
  }
}

static double atan2_0_2pi(double y,double x)
{
  double a;

  a = atan2(y,x);
  if (a < 0) {
    a += 2.0*M_PI;
  }
  return a;
}

/* calculates center of mass of selection index from all coordinates x */
void pull_calc_coms(t_commrec *cr,
		    t_pull *pull, t_mdatoms *md, t_pbc *pbc,
		    rvec x[], rvec *xp)
{
  static rvec *rbuf=NULL;
  static dvec *dbuf=NULL;
  int  g,i,ii,m;
  real mass,w,wm,twopi_box=0;
  double wmass,wwmass,invwmass;
  dvec com,comp;
  double cm,sm,cmp,smp,ccm,csm,ssm,csw,snw;
  rvec *xx[2],x_pbc={0,0,0},dx;
  t_pullgrp *pgrp;

  if (rbuf == NULL) {
    snew(rbuf,1+pull->ngrp);
  }
  if (dbuf == NULL) {
    snew(dbuf,3*(1+pull->ngrp));
  }

  if (pull->bRefAt) {
    pull_set_pbcatoms(cr,pull,md,x,rbuf);
  }

  if (pull->cosdim >= 0) {
    for(m=pull->cosdim+1; m<pull->npbcdim; m++) {
      if (pbc->box[m][pull->cosdim] != 0) {
	gmx_fatal(FARGS,"Can not do cosine weighting for trilinic dimensions");
      }
    }
    twopi_box = 2.0*M_PI/pbc->box[pull->cosdim][pull->cosdim];
  }
  
  for (g=0; g<1+pull->ngrp; g++) {
    pgrp = &pull->grp[g];
    clear_dvec(com);
    clear_dvec(comp);
    wmass  = 0;
    wwmass = 0;
    cm  = 0;
    sm  = 0;
    cmp = 0;
    smp = 0;
    ccm = 0;
    csm = 0;
    ssm = 0;
    if (!(g==0 && PULL_CYL(pull))) {
      if (pgrp->epgrppbc == epgrppbcREFAT) {
	/* Set the pbc atom */
	copy_rvec(rbuf[g],x_pbc);
      }
      w = 1;
      for(i=0; i<pgrp->nat_loc; i++) {
	ii = pgrp->ind_loc[i];
	mass = md->massT[ii];
	if (pgrp->epgrppbc != epgrppbcCOS) {
	  if (pgrp->weight_loc) {
	    w = pgrp->weight_loc[i];
	  }
	  wm = w*mass;
	  wmass  += wm;
	  wwmass += wm*w;
	  if (pgrp->epgrppbc == epgrppbcNONE) {
	    /* Plain COM: sum the coordinates */
	    for(m=0; m<DIM; m++)
	      com[m]    += wm*x[ii][m];
	    if (xp) {
	      for(m=0; m<DIM; m++)
		comp[m] += wm*xp[ii][m];
	    }
	  } else {
	    /* Sum the difference with the reference atom */
	    pbc_dx(pbc,x[ii],x_pbc,dx);
	    for(m=0; m<DIM; m++)
	      com[m]    += wm*dx[m];
	    if (xp) {
	      pbc_dx(pbc,xp[ii],x_pbc,dx);
	      for(m=0; m<DIM; m++)
		comp[m] += wm*dx[m];
	    }
	  }
	} else {
	  /* Determine cos and sin sums */
	  csw = cos(x[ii][pull->cosdim]*twopi_box);
	  snw = sin(x[ii][pull->cosdim]*twopi_box);
	  cm  += csw*mass;
	  sm  += snw*mass;
	  ccm += csw*csw*mass;
	  csm += csw*snw*mass;
	  ssm += snw*snw*mass;

	  if (xp) {
	    csw = cos(xp[ii][pull->cosdim]*twopi_box);
	    snw = sin(xp[ii][pull->cosdim]*twopi_box);
	    cmp += csw*mass;
	    smp += snw*mass;
 	  }
	}
      }
    }

    /* Copy local sums to a buffer for global summing */
    switch (pgrp->epgrppbc) {
    case epgrppbcNONE:
    case epgrppbcREFAT:
      copy_dvec(com,dbuf[g*3]);
      copy_dvec(comp,dbuf[g*3+1]);
      dbuf[g*3+2][0] = wmass;
      dbuf[g*3+2][1] = wwmass;
      dbuf[g*3+2][2] = 0;
      break;
    case epgrppbcCOS:
      dbuf[g*3  ][0] = cm;
      dbuf[g*3  ][1] = sm;
      dbuf[g*3  ][2] = 0;
      dbuf[g*3+1][0] = ccm;
      dbuf[g*3+1][1] = csm;
      dbuf[g*3+1][2] = ssm;
      dbuf[g*3+2][0] = cmp;
      dbuf[g*3+2][1] = smp;
      dbuf[g*3+2][2] = 0;
      break;
    }
  }

  if (cr && PAR(cr)) {
    /* Sum the contributions over the nodes */
    gmx_sumd((1+pull->ngrp)*3*DIM,dbuf[0],cr);
  }
  
  for (g=0; g<1+pull->ngrp; g++) {
    pgrp = &pull->grp[g];
    if (pgrp->nat > 0 && !(g==0 && PULL_CYL(pull))) {
      if (pgrp->epgrppbc != epgrppbcCOS) {
	/* Determine the inverse mass */
	wmass  = dbuf[g*3+2][0];
	wwmass = dbuf[g*3+2][1];
	invwmass = 1/wmass;
	/* invtm==0 signals a frozen group, so then we should keep it zero */
	if (pgrp->invtm > 0) {
	  pgrp->wscale = wmass/wwmass;
	  pgrp->invtm  = 1.0/(pgrp->wscale*wmass);
	}
	/* Divide by the total mass */
	for(m=0; m<DIM; m++) {
	  pgrp->x[m]    = dbuf[g*3  ][m]*invwmass;
	  if (xp) {
	    pgrp->xp[m] = dbuf[g*3+1][m]*invwmass;
	  }
	  if (pgrp->epgrppbc == epgrppbcREFAT) {
	    pgrp->x[m]    += rbuf[g][m];
	    if (xp) {
	      pgrp->xp[m] += rbuf[g][m];
	    }
	  }
	}
      } else {
	/* Determine the optimal location of the cosine weight */
	csw = dbuf[g*3][0];
	snw = dbuf[g*3][1];
	pgrp->x[pull->cosdim] = atan2_0_2pi(snw,csw)/twopi_box;
	/* Set the weights for the local atoms */
	wmass = sqrt(csw*csw + snw*snw);
	wwmass = (dbuf[g*3+1][0]*csw*csw + dbuf[g*3+1][1]*csw*snw + dbuf[g*3+1][2]*snw*snw)/(wmass*wmass);
	pgrp->wscale = wmass/wwmass;
	pgrp->invtm  = 1.0/(pgrp->wscale*wmass);
	/* Set the weights for the local atoms */
	csw *= pgrp->invtm;
	snw *= pgrp->invtm;
	for(i=0; i<pgrp->nat_loc; i++) {
	  ii = pgrp->ind_loc[i];
	  pgrp->weight_loc[i] = csw*cos(twopi_box*x[ii][pull->cosdim]) +
				snw*sin(twopi_box*x[ii][pull->cosdim]);
	}
	if (xp) {
	  csw = dbuf[g*3+2][0];
	  snw = dbuf[g*3+2][1];
	  pgrp->xp[pull->cosdim] = atan2_0_2pi(snw,csw)/twopi_box;
	}
      }
      if (debug) {
	fprintf(debug,"Pull group %d wmass %f wwmass %f invtm %f\n",
		g,wmass,wwmass,pgrp->invtm);
      }
    }
  }
  
  if (PULL_CYL(pull)) {
    /* Calculate the COMs for the cyclinder reference groups */
    make_cyl_refgrps(cr,pull,md,pbc,x,xp);
  }  
}
