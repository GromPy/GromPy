#ifndef _topsort_h
#define _topsort_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef CPLUSPLUS 
extern "C" {
#endif

#include "typedefs.h"

  /* Returns if the are bonded interactions for free energy calculations */
  extern bool gmx_mtop_bondeds_free_energy(gmx_mtop_t *mtop);

  /* Sort all the bonded ilists in idef to have the perturbed ones at the end
   * and set nr_nr_nonperturbed in ilist.
   */
  extern void gmx_sort_ilist_fe(t_idef *idef);

#ifdef CPLUSPLUS
}
#endif

#endif
