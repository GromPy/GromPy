/*
 * $Id: commrec.h,v 1.51 2008/08/25 14:19:12 hess Exp $
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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_MPI
#include <mpi.h>
#endif

#include "idef.h"

#define DD_MAXCELL  8
#define DD_MAXICELL 4

typedef struct gmx_domdec_master *gmx_domdec_master_p_t;

typedef struct {
  int  j0;       /* j-cell start               */
  int  j1;       /* j-cell end                 */
  int  cg1;      /* i-charge-group end         */
  int  jcg0;     /* j-charge-group start       */
  int  jcg1;     /* j-charge-group end         */
  ivec shift0;   /* Minimum shifts to consider */
  ivec shift1;   /* Maximum shifts to consider */
} gmx_domdec_ns_ranges_t;

typedef struct {
  int     cell;
  atom_id a;
} gmx_ga2la_t;

typedef struct gmx_reverse_top *gmx_reverse_top_p_t;

typedef struct gmx_domdec_constraints *gmx_domdec_constraints_p_t;

typedef struct gmx_domdec_specat_comm *gmx_domdec_specat_comm_p_t;

typedef struct gmx_domdec_comm *gmx_domdec_comm_p_t;

typedef struct gmx_pme_comm_n_box *gmx_pme_comm_n_box_p_t;

typedef struct {
  /* The DD particle-particle nodes only */
  /* The communication setup within the communicator all
   * defined in dd->comm in domdec.c
   */
  int  nnodes;
#ifdef GMX_MPI
  MPI_Comm mpi_comm_all;
#endif
  /* Use MPI_Sendrecv communication instead of non-blocking calls */
  bool bSendRecv2;
  /* The local DD cell index and rank */
  ivec ci;
  int  rank;
  ivec master_ci;
  int  masterrank;
  /* Communication with the PME only nodes */
  int  pme_nodeid;
  bool pme_receive_vir_ener;
  gmx_pme_comm_n_box_p_t cnb;
#ifdef GMX_MPI
  int  nreq_pme;
  MPI_Request req_pme[4];
#endif
  

  /* The communication setup, identical for each cell, cartesian index */
  ivec nc;
  int  ndim;
  ivec dim;  /* indexed by 0 to ndim */
  bool bGridJump;
  /* The bonded and non-bonded communication setup, cartesian index */
  int  ncell;
  ivec shift[DD_MAXCELL];

  /* Screw PBC? */
  bool bScrewPBC;

  /* Tells if the box is skewed for each of the three cartesian directions */
  ivec tric_dir;
  rvec skew_fac;

  /* The cell boundaries */
  rvec cell_x0;
  rvec cell_x1;
  /* The extreme sizes of the local cells needed for the neighbor searching */
  rvec cell_ns_x0;
  rvec cell_ns_x1;
  /* Forward and backward neighboring cells, indexed by 0 to ndim */
  int  neighbor[DIM][2];

  /* Only available on the master node */
  gmx_domdec_master_p_t ma;

  /* Are there inter charge group constraints */
  bool bInterCGcons;

  /* Global atom number to interaction list */
  gmx_reverse_top_p_t reverse_top;
  int  nbonded_global;
  int  nbonded_local;

  /* The number of inter charge-group exclusions */
  int  n_intercg_excl;

  /* Vsite stuff */
  int  *ga2la_vsite;
  gmx_domdec_specat_comm_p_t vsite_comm;

  /* Constraint stuff */
  gmx_domdec_constraints_p_t constraints;
  gmx_domdec_specat_comm_p_t constraint_comm;

  /* The charge group boundaries for the cells */
  int ncg_cell[DD_MAXCELL+1];

  /* The local to gobal charge group index and local cg to local atom index */
  int  ncg_home;
  int  ncg_tot;
  int  *index_gl;
  int  *cgindex;
  int  cg_nalloc;
  /* Local atom to local cg index, only for special cases */
  int  *la2lc;
  int  la2lc_nalloc;

  /* The number of home atoms */
  int  nat_home;
  /* The total number of atoms: home and received zones */
  int  nat_tot;
  /* Index from the local atoms to the global atoms */
  int  *gatindex;
  int  gatindex_nalloc;

  /* Global atom number to local atom number, -1 if not local */
  gmx_ga2la_t *ga2la;

  /* For neighborsearching */
  int  nicell;
  gmx_domdec_ns_ranges_t icell[DD_MAXICELL];

  /* Communication stuff */
  gmx_domdec_comm_p_t comm;

  /* The partioning count, to keep track of the state */
  int ddp_count;
} gmx_domdec_t;

typedef struct gmx_partdec *gmx_partdec_p_t;

typedef struct {
  int nsim;
  int sim;
#ifdef GMX_MPI
  MPI_Group mpi_group_masters;
  MPI_Comm mpi_comm_masters;
#endif
} gmx_multisim_t;

#define DUTY_PP  (1<<0)
#define DUTY_PME (1<<1)

typedef struct {
  int      bUse;
#ifdef GMX_MPI
  MPI_Comm comm_intra;
  int      rank_intra;
  MPI_Comm comm_inter;
#endif
  
} gmx_nodecomm_t;

typedef struct {
  /* The nodids in one sim are numbered sequentially from 0.
   * All communication within some simulation should happen
   * in mpi_comm_mysim, or its subset mpi_comm_mygroup.
   */
  int sim_nodeid,nnodes,npmenodes;
  int threadid,nthreads;
  /* The nodeid in the PP/PME, PP or PME group */
  int nodeid;
#ifdef GMX_MPI
  MPI_Comm mpi_comm_mysim;
  MPI_Comm mpi_comm_mygroup;
#endif
  
  gmx_nodecomm_t nc;
  
  /* For domain decomposition */
  gmx_domdec_t *dd;

  /* For particle decomposition */
  gmx_partdec_p_t pd;

  /* The duties of this node, see the defines above */
  int duty;

  gmx_multisim_t *ms;
} t_commrec;

#define MASTERNODE(cr)     ((cr)->nodeid == 0)
#define MASTERTHREAD(cr)   ((cr)->threadid == 0)
#define MASTER(cr)         (MASTERNODE(cr) && MASTERTHREAD(cr))
#define SIMMASTER(cr)      (MASTER(cr) && ((cr)->duty & DUTY_PP))
#define NODEPAR(cr)        ((cr)->nnodes > 1)
#define THREADPAR(cr)      ((cr)->nthreads > 1)
#define PAR(cr)            (NODEPAR(cr) || THREADPAR(cr))
#define RANK(cr,nodeid)    (nodeid)
#define MASTERRANK(cr)     (0)

#define DOMAINDECOMP(cr)   ((cr)->dd != NULL)
#define DDMASTER(dd)       ((dd)->rank == (dd)->masterrank)

#define PARTDECOMP(cr)     ((cr)->pd != NULL)

#define MULTISIM(cr)       ((cr)->ms)
#define MSRANK(ms,nodeid)  (nodeid)
#define MASTERSIM(ms)      ((ms)->sim == 0)

/* The master of all (the node that prints the remaining run time etc.) */
#define MULTIMASTER(cr)    (MASTER(cr) && (!MULTISIM(cr) || MASTERSIM((cr)->ms)))
