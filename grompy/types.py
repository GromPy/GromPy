from os import environ
from ctypes import c_char_p,\
                   c_int,\
                   c_char,\
                   c_ubyte,\
                   c_ushort,\
                   Structure,\
                   Union,\
                   POINTER,\
                   c_double,\
                   c_uint,\
                   c_ulonglong,\
                   c_void_p,\
                   c_ulong,\
                   c_float,\
                   c_long

from grompy import c_real,rvec,ivec,tensor,matrix,splinevec,dvec

time_t = c_long

#typedef int atom_id;
atom_id = c_int
#typedef int t_functype;
t_functype = c_int
#typedef atom_id t_iatom;
t_iatom = atom_id

NULL  = None
TRUE  = 1
FALSE = 0

DIM = 3
MAXFORCEPARAM = 12
NR_RBDIHS = 6

STRLEN = 4096

# from simple.h
XX = 0
YY = 1
ZZ = 2

#enum {
(
  F_BONDS,
  F_G96BONDS,
  F_MORSE,
  F_CUBICBONDS,
  F_CONNBONDS,
  F_HARMONIC,
  F_FENEBONDS,
  F_TABBONDS,
  F_TABBONDSNC,
  F_ANGLES,
  F_G96ANGLES,
  F_CROSS_BOND_BONDS,
  F_CROSS_BOND_ANGLES,
  F_UREY_BRADLEY,
  F_QUARTIC_ANGLES,
  F_TABANGLES,
  F_PDIHS,
  F_RBDIHS,
  F_FOURDIHS,
  F_IDIHS,
  F_PIDIHS,
  F_TABDIHS,
  F_LJ14,
  F_COUL14,
  F_LJC14_Q,
  F_LJC_PAIRS_NB,
  F_LJ,
  F_BHAM,
  F_LJ_LR,
  F_BHAM_LR,
  F_DISPCORR,
  F_COUL_SR,
  F_COUL_LR,
  F_RF_EXCL,
  F_COUL_RECIP,
  F_DPD,
  F_POLARIZATION,
  F_WATER_POL,
  F_THOLE_POL,
  F_POSRES,
  F_DISRES,
  F_DISRESVIOL,
  F_ORIRES,
  F_ORIRESDEV,
  F_ANGRES,
  F_ANGRESZ,
  F_DIHRES,
  F_DIHRESVIOL,
  F_CONSTR,
  F_CONSTRNC,
  F_SETTLE,
  F_VSITE2,
  F_VSITE3,
  F_VSITE3FD,
  F_VSITE3FAD,
  F_VSITE3OUT,
  F_VSITE4FD,
  F_VSITE4FDN,
  F_VSITEN,
  F_COM_PULL,
  F_EQM,
  F_EPOT,
  F_EKIN,
  F_ETOT,
  F_ECONSERVED,
  F_TEMP,
  F_PRES,
  F_DVDL,
  F_DKDL,
  F_DGDL_CON,
  F_NRE) = map(c_int, xrange(71))
#/* This number is for the total number of energies	*/
#};

class t_block(Structure):
    _fields_ = [("nr", c_int),
                ("index", POINTER(atom_id)),
                ("nalloc_index", c_int)]
#  int nr;			/* The number of blocks			*/
#  atom_id *index;		/* Array of indices (dim: nr+1) 	*/
#  int nalloc_index;             /* The allocation size for index        */

class t_blocka(Structure):
    _fields_ = [("nr", c_int),
                ("index", POINTER(atom_id)),
                ("nra", c_int),
                ("a", POINTER(atom_id)),
                ("nalloc_index", c_int),
                ("nalloc_a", c_int)]
#  int nr;			/* The number of blocks			*/
#  atom_id *index;		/* Array of indices in a (dim: nr+1)	*/
#  int nra;			/* The number of atoms 			*/
#  atom_id *a;			/* Array of atom numbers in each group 	*/
#				/* (dim: nra)				*/
#				/* Block i (0<=i<nr) runs from		*/
#				/* index[i] to index[i+1]-1. There will */
#				/* allways be an extra entry in index	*/
#				/* to terminate the table		*/
#  int nalloc_index;             /* The allocation size for index        */
#  int nalloc_a;                 /* The allocation size for a            */

class t_symbuf(Structure):
    pass

t_symbuf._fields_ = [("bufsize", c_int),
                     ("buf", POINTER(c_char_p)),
                     ("next", POINTER(t_symbuf))]
#  int bufsize;
#  char **buf;
#  struct symbuf *next;

class t_symtab(Structure):
    _fields_ = [("nr", c_int),
                ("symbuf", POINTER(t_symbuf))]
#  int      nr;
#  t_symbuf *symbuf;

class t_iparams(Union):
    class bham(Structure):
        _fields_ = [("a", c_real),
                    ("b", c_real),
                    ("c", c_real)]
    class harmonic(Structure):
        _fields_ = [("rA", c_real),
                    ("krA", c_real),
                    ("rb", c_real),
                    ("krB", c_real)]
    class cubic(Structure):
        _fields_ = [("b0", c_real),
                    ("kb", c_real),
                    ("kcub", c_real)]
    class fene(Structure):
        _fields_ = [("bm", c_real),
                    ("kb", c_real)]
    class cross_bb(Structure):
        _fields_ = [("r1e", c_real),
                    ("r2e", c_real),
                    ("krr", c_real)]
    class cross_ba(Structure):
        _fields_ = [("r1e", c_real),
                    ("r2e", c_real),
                    ("r3e", c_real),
                    ("krt", c_real)]
    class u_b(Structure):
        _fields_ = [("theta", c_real),
                    ("ktheta", c_real),
                    ("r13", c_real),
                    ("kUB", c_real)]
    class qangle(Structure):
        _fields_ = [("theta", c_real),
                    ("c", c_real * 5)]
    class polarize(Structure):
        _fields_ = [("alpha", c_real)]
    class wpol(Structure):
        _fields_ = [("al_x", c_real),
                    ("al_y", c_real),
                    ("al_z", c_real),
                    ("rOH", c_real),
                    ("rHH", c_real),
                    ("rOD", c_real)]
    class thole(Structure):
        _fields_ = [("a", c_real),
                    ("alpha1", c_real),
                    ("alpha2", c_real),
                    ("rfac", c_real)]
    class lj(Structure):
        _fields_ = [("c6", c_real),
                    ("c12", c_real)]
    class lj14(Structure):
        _fields_ = [("c6A", c_real),
                    ("c12A", c_real),
                    ("c6B", c_real),
                    ("c12B", c_real)]
    class ljc14(Structure):
        _fields_ = [("fqq", c_real),
                    ("qi", c_real),
                    ("qj", c_real),
                    ("c6", c_real),
                    ("c12", c_real)]
    class ljcnb(Structure):
        _fields_ = [("qi", c_real),
                    ("qj", c_real),
                    ("c6", c_real),
                    ("c12", c_real)]
    class pdihs(Structure):
        _fields_ = [("phiA", c_real),
                    ("cpA", c_real),
                    ("mult", c_int),
                    ("phiB", c_real),
                    ("cpB", c_real)]
    class constr(Structure):
        _fields_ = [("dA", c_real),
                    ("dB", c_real)]
    class settle(Structure):
        _fields_ = [("doh", c_real),
                    ("dhh", c_real)]
    class morse(Structure):
        _fields_ = [("b0", c_real),
                    ("cb", c_real),
                    ("beta", c_real)]
    class posres(Structure):
        _fields_ = [("pos0A", c_real * DIM),
                    ("fcA", c_real * DIM),
                    ("pos0B", c_real * DIM),
                    ("fcB", c_real * DIM)]
    class rbdihs(Structure):
        _fields_ = [("rbcA", c_real * NR_RBDIHS),
                    ("rbcB", c_real * NR_RBDIHS)]
    class vsite(Structure):
        _fields_ = [("a", c_real),
                    ("b", c_real),
                    ("c", c_real),
                    ("d", c_real),
                    ("e", c_real),
                    ("f", c_real)]
    class vsiten(Structure):
        _fields_ = [("n", c_int),
                    ("a", c_real)]
    class disres(Structure):
        _fields_ = [("low", c_real),
                    ("up1", c_real),
                    ("up2", c_real),
                    ("kfac", c_real),
                    ("type", c_int),
                    ("label", c_int)]
    class dihres(Structure):
        _fields_ = [("phi", c_real),
                    ("dphi", c_real),
                    ("kfac", c_real),
                    ("label", c_int),
                    ("power", c_int)]
    class orires(Structure):
        _fields_ = [("ex", c_int),
                    ("power", c_int),
                    ("label", c_int),
                    ("c", c_real),
                    ("obs", c_real),
                    ("kfac", c_real)]
    class tab(Structure):
        _fields_ = [("table", c_int),
                    ("kA", c_real),
                    ("kB", c_real)]
    class generic(Structure):
        _fields_ = [("buf", c_real * MAXFORCEPARAM)]

    _fields_ = [("bham", bham),
                ("harmonic", harmonic),
                ("cubic", cubic),
                ("fene", fene),
                ("cross_bb", cross_bb),
                ("cross_ba", cross_ba),
                ("u_b", u_b),
                ("qangle", qangle),
                ("polarize", polarize),
                ("wpol", wpol),
                ("thole", thole),
                ("lj", lj),
                ("lj14", lj14),
                ("ljc14", ljc14),
                ("ljcnb", ljcnb),
                ("pdihs", pdihs),
                ("constr", constr),
                ("settle", settle),
                ("morse", morse),
                ("posres", posres),
                ("rbdihs", rbdihs),
                ("vsite", vsite),
                ("vsiten", vsiten),
                ("disres", disres),
                ("dihres", dihres),
                ("orires", orires),
                ("tab", tab),
                ("generic", generic)]

#  /* Some parameters have A and B values for free energy calculations.
#   * The B values are not used for regular simulations of course.
#   * Free Energy for nonbondeds can be computed by changing the atom type.
#   * The harmonic type is used for all harmonic potentials:
#   * bonds, angles and improper dihedrals
#   */
#  struct {real a,b,c;	                                   } bham;
#  struct {real rA,krA,rB,krB;           	           } harmonic;
#  /* No free energy supported for cubic bonds, FENE, WPOL or cross terms */
#  struct {real b0,kb,kcub;                                 } cubic;
#  struct {real bm,kb;                                      } fene;
#  struct {real r1e,r2e,krr;                                } cross_bb;
#  struct {real r1e,r2e,r3e,krt;                            } cross_ba;
#  struct {real theta,ktheta,r13,kUB;                       } u_b;
#  struct {real theta,c[5];                                 } qangle;
#  struct {real alpha;                                      } polarize;
#  struct {real al_x,al_y,al_z,rOH,rHH,rOD;                 } wpol;
#  struct {real a,alpha1,alpha2,rfac;                       } thole;
#  struct {real c6,c12;				           } lj;
#  struct {real c6A,c12A,c6B,c12B;		           } lj14;
#  struct {real fqq,qi,qj,c6,c12;	                   } ljc14;
#  struct {real qi,qj,c6,c12;		                   } ljcnb;
#  /* Proper dihedrals can not have different multiplicity when
#   * doing free energy calculations, because the potential would not
#   * be periodic anymore.
#   */
#  struct {real phiA,cpA;int mult;real phiB,cpB;            } pdihs;
#  struct {real dA,dB;		        	           } constr;
#  /* Settle can not be used for Free energy calculations of water bond geometry.
#   * Use shake (or lincs) instead if you have to change the water bonds.
#   */
#  struct {real doh,dhh;                                   } settle;
#  /* No free energy supported for morse bonds */
#  struct {real b0,cb,beta;                        	  } morse;
#  struct {real pos0A[DIM],fcA[DIM],pos0B[DIM],fcB[DIM];   } posres;
#  struct {real rbcA[NR_RBDIHS], rbcB[NR_RBDIHS];          } rbdihs;
#  struct {real a,b,c,d,e,f;                               } vsite;
#  struct {int  n; real a;                                 } vsiten;
#  struct {real low,up1,up2,kfac;int type,label;           } disres;
#  struct {real phi,dphi,kfac;int label,power;             } dihres;
#  struct {int  ex,power,label; real c,obs,kfac;           } orires;
#  struct {int  table;real kA;real kB;                     } tab;
#  struct {real buf[MAXFORCEPARAM];	  	          } generic; /* Conversion */

class t_ilist(Structure):
    _fields_ = [("nr", c_int),
                ("nr_nonperturbed", c_int),
                ("iatoms", POINTER(t_iatom)),
                ("nalloc", c_int)]
#  int nr;
#  int nr_nonperturbed;
#  t_iatom *iatoms;
#  int nalloc;

class t_idef(Structure):
    _fields_ = [("ntypes", c_int),
                ("atnr", c_int),
                ("t_functype", POINTER(t_functype)),
                ("t_iparams", POINTER(t_iparams)),
                ("fudgeQQ", c_real),
                ("iparams_posres", POINTER(t_iparams)),
                ("iparams_posres_nalloc", c_int),
                ("il", t_ilist * F_NRE.value),
                ("ilsort", c_int)]
#  int ntypes;
#  int atnr;
#  t_functype *functype;
#  t_iparams  *iparams;
#  real fudgeQQ;
#  t_iparams  *iparams_posres;
#  int iparams_posres_nalloc;
#  t_ilist il[F_NRE];
#  int ilsort;


class t_atom(Structure):
    _fields_ = [("m", c_real),
                ("q", c_real),
                ("mB", c_real),
                ("qB", c_real),
                ("type", c_ushort),
                ("typeB", c_ushort),
                ("ptype", c_int),
                ("resnr", c_int),
                ("atomnumber", c_int),
                ("chain", c_ubyte)]
#  real 		m,q;		/* Mass and charge			*/
#  real 		mB,qB;		/* Mass and charge for Free Energy calc */
#  unsigned short type;		/* Atom type				*/
#  unsigned short typeB;		/* Atom type for Free Energy calc	*/
#  int           ptype;		/* Particle type			*/
#  int 		resnr;		/* Residue number			*/
#  int           atomnumber;     /* Atomic Number or NOTSET              */
#  unsigned char chain;          /* chain identifier                     */

class t_pdbinfo(Structure):
    _fields_ = [("type", c_int),
                ("atomnr", c_int),
                ("altloc", c_char),
                ("atomnm", c_char * 6),
                ("pdbresnr", c_char * 6),
                ("occup", c_real),
                ("bfac", c_real),
                ("bAnisotropic", c_int),
                ("uij", c_int * 6)]
#  int  type;                    /* PDB record name                      */
#  int  atomnr;                  /* PDB atom number                      */
#  char altloc;                  /* Alternate location indicator         */
#  char atomnm[6];               /* True atom name including spaces      */
#  char pdbresnr[6];             /* PDB res number                       */
#  real occup;                   /* Occupancy                            */
#  real bfac;                    /* B-factor                             */
#  bool bAnisotropic;            /* (an)isotropic switch                 */
#  int  uij[6];                  /* Anisotropic B-factor                 */

class t_grps(Structure):
    _fields_ = [("nr", c_int),
                ("nm_ind", POINTER(c_int))]
#  int  nr;			/* Number of different groups		*/
#  int  *nm_ind;                 /* Index in the group names             */

class t_atoms(Structure):
    _fields_ = [("nr", c_int),
                ("atom", POINTER(t_atom)),
                ("atomname", POINTER(POINTER(c_char_p))),
                ("atomtype", POINTER(POINTER(c_char_p))),
                ("atomtypeB", POINTER(POINTER(c_char_p))),
                ("nres", c_int),
                ("resname", POINTER(POINTER(c_char_p))),
                ("pdbinfo", c_char_p)]
#  int           nr;             /* Nr of atoms                          */
#  t_atom	*atom;		/* Array of atoms (dim: nr)		*/
#				/* The following entries will not 	*/
#				/* allways be used (nres==0)	 	*/
#  char		***atomname;	/* Array of pointers to atom name	*/
#				/* use: (*(atomname[i]))		*/
#  char		***atomtype;	/* Array of pointers to atom types	*/
#				/* use: (*(atomtype[i]))		*/
#  char		***atomtypeB;	/* Array of pointers to B atom types	*/
#				/* use: (*(atomtypeB[i]))		*/
#  int		nres;		/* Nr of residue names			*/
#  char		***resname; 	/* Array of pointers to residue names 	*/
#				/* use: (*(resname[i]))	       	*/
#  t_pdbinfo     *pdbinfo;       /* PDB Information, such as aniso. Bfac */

class t_atomtypes(Structure):
    _fields_ = [("nr", c_int),
                ("radius", POINTER(c_real)),
                ("vol", POINTER(c_real)),
                ("surftens", POINTER(c_real)),
                ("atomnumber", POINTER(c_int))]
#  int           nr;              /* number of atomtypes                     */
#  real         *radius;         /* GBSA radius for each atomtype        */
#  real         *vol;            /* GBSA efective volume for each atomtype   */
#  real         *surftens;       /* implicit solvent surftens for each atomtype */
#  int          *atomnumber;     /* Atomic number, used for QM/MM */

class t_topology(Structure):
    _fields_ = [("name", POINTER(c_char_p)),
                ("idef", t_idef),
                ("atoms", t_atoms),
                ("atomtypes", t_atomtypes),
                ("cgs", t_block),
                ("mols", t_block),
                ("excls", t_blocka),
                ("symtab", t_symtab)]
#  char  	**name;		/* Name of the topology	       	        */
#  t_idef	idef;		/* The interaction function definition	*/
#  t_atoms	atoms;		/* The atoms		       	        */
#  t_atomtypes   atomtypes;      /* Atomtype properties                  */
#  t_block       cgs;            /* The charge groups                    */
#  t_block       mols;           /* The molecules                        */
#  t_blocka      excls;          /* The exclusions                       */
#  t_symtab	symtab;		/* The symbol table			*/

class t_cosines(Structure):
    _fields_ = [("n", c_int),
                ("a", POINTER(c_real)),
                ("phi", POINTER(c_real))]
#  int  n;        /* Number of terms                */
#  real *a;        /* Coeffients (V / nm )                      */
#  real *phi;        /* Phase angles                    */

class t_grpopts(Structure):
    _fields_ = [("ngtc", c_int),
                ("ngacc", c_int),
                ("ngfrz", c_int),
                ("ngener", c_int),
                ("nrdf", POINTER(c_real)),
                ("ref_t", POINTER(c_real)),
                ("annealing", POINTER(c_int)),
                ("anneal_npoints", POINTER(c_int)),
                ("anneal_time", POINTER(POINTER(c_real))),
                ("anneal_temp", POINTER(POINTER(c_real))),
                ("tau_t", POINTER(c_real)),
                ("acc", POINTER(rvec)),
                ("nFreeze", POINTER(ivec)),
                ("egp_flags", POINTER(c_int)),
                ("ngQM", c_int),
                ("QMmethod", POINTER(c_int)),
                ("QMbasis", POINTER(c_int)),
                ("QMcharge", POINTER(c_int)),
                ("QMmult", POINTER(c_int)),
                ("bSH", POINTER(c_int)),
                ("CASorbitals", POINTER(c_int)),
                ("CASelectrons", POINTER(c_int)),
                ("SAon", POINTER(c_real)),
                ("SAoff", POINTER(c_real)),
                ("SAsteps", POINTER(c_int)),
                ("bOPT", POINTER(c_int)),
                ("bTS", POINTER(c_int))]
#  int     ngtc;                  /* # T-Coupl groups                        */
#  int     ngacc;                 /* # Accelerate groups                     */
#  int     ngfrz;                 /* # Freeze groups                         */
#  int     ngener;             /* # Ener groups                */
#  real    *nrdf;             /* Nr of degrees of freedom in a group        */
#  real    *ref_t;             /* Coupling temperature    per group   */
#  int     *annealing;            /* No/simple/periodic SA for each group    */
#  int     *anneal_npoints;       /* Number of annealing time points per grp */
#  real    **anneal_time;         /* For ea. group: Time points              */
#  real    **anneal_temp;         /* For ea. grp: Temperature at these times */
#                                 /* Final temp after all intervals is ref_t */
#  real    *tau_t;             /* Tau coupling time                 */
#  rvec    *acc;                 /* Acceleration per group            */
#  ivec    *nFreeze;             /* Freeze the group in each direction ?    */
#  int     *egp_flags;            /* Exclusions/tables of energy group pairs */
#
#  /* QMMM stuff */
#  int     ngQM;         /* nr of QM groups                              */
#  int     *QMmethod;    /* Level of theory in the QM calculation        */
#  int     *QMbasis;     /* Basisset in the QM calculation               */
#  int     *QMcharge;    /* Total charge in the QM region                */
#  int     *QMmult;      /* Spin multiplicicty in the QM region          */
#  bool    *bSH;         /* surface hopping (diabatic hop only)          */
#  int     *CASorbitals; /* number of orbiatls in the active space       */
#  int     *CASelectrons;/* number of electrons in the active space      */
#  real    *SAon;        /* at which gap (A.U.) the SA is switched on    */
#  real    *SAoff;
#  int     *SAsteps;     /* in how many steps SA goes from 1-1 to 0.5-0.5*/
#  bool    *bOPT;
#  bool    *bTS;

class t_pullgrp(Structure):
    _fields_ = [("nat", c_int),
                ("ind", POINTER(atom_id)),
                ("nat_loc", c_int),
                ("nalloc_loc", c_int),
                ("ind_loc", POINTER(atom_id)),
                ("nweight", c_int),
                ("weight", POINTER(c_real)),
                ("weight_loc", POINTER(c_real)),
                ("epgrppbc", c_int),
                ("pbcatom", atom_id),
                ("vec", rvec),
                ("init", rvec),
                ("rate", c_real),
                ("k", c_real),
                ("kB", c_real),
                ("wscale", c_real),
                ("invtm" , c_real),
                ("x", dvec),
                ("xp", dvec),
                ("dr", dvec),
                ("f_scal", c_double),
                ("f", dvec)]
#  int        nat;      /* Number of atoms in the pull group */
#  atom_id    *ind;     /* The global atoms numbers */
#  int        nat_loc;  /* Number of local pull atoms */
#  int        nalloc_loc; /* Allocation size for ind_loc and weight_loc */
#  atom_id    *ind_loc; /* Local pull indices */
#  int        nweight;  /* The number of weights (0 or nat) */
#  real       *weight;  /* Weights (use all 1 when weight==NULL) */
#  real       *weight_loc; /* Weights for the local indices */
#  int        epgrppbc; /* The type of pbc for this pull group, see enum above */
#  atom_id    pbcatom;  /* The reference atom for pbc (global number) */
#  rvec       vec;      /* The pull vector, direction or position */
#  rvec       init;     /* Initial reference displacement */
#  real       rate;     /* Rate of motion (nm/ps) */
#  real       k;        /* force constant */
#  real       kB;       /* force constant for state B */
#  real       wscale;   /* scaling factor for the weights: sum w m/sum w w m */
#  real       invtm;    /* inverse total mass of the group: 1/wscale sum w m */
#  dvec       x;        /* center of mass before update */
#  dvec       xp;       /* center of mass after update before constraining */
#  dvec       dr;       /* The distance from the reference group */
#  double     f_scal;   /* Scalar force for directional pulling */
#  dvec       f;        /* force due to the pulling/constraining */

class t_pull(Structure):
    _fields_ = [("ngrp", c_int),
                ("eGeom", c_int),
                ("dim", ivec),
                ("cyl_r1", c_real),
                ("cyl_r0", c_real),
                ("constr_tol", c_real),
                ("nstxout",c_int),
                ("nstfout",c_int),
                ("ePBC", c_int),
                ("npbcdim", c_int),
                ("bRefAt", c_int),
                ("cosdim", c_int),
                ("grp", POINTER(t_pullgrp)),
                ("dyna", POINTER(t_pullgrp)),
                ("out_x", c_char_p),
                ("out_f", c_char_p)]
#  int        ngrp;        /* number of groups */
#  int        eGeom;       /* pull geometry */
#  ivec       dim;         /* used to select components for constraint */
#  real       cyl_r1;      /* radius of cylinder for dynamic COM */
#  real       cyl_r0;      /* radius of cylinder including switch length */
#  real       constr_tol;  /* absolute tolerance for constraints in (nm) */
#  int        nstxout;     /* Output frequency for pull x */
#  int        nstfout;     /* Output frequency for pull f */
#  int        ePBC;        /* the boundary conditions */
#  int        npbcdim;     /* do pbc in dims 0 <= dim < npbcdim */
#  bool       bRefAt;      /* do we need reference atoms for a group COM ? */
#  int        cosdim;      /* dimension for cosine weighting, -1 if none */
#  t_pullgrp  *grp;        /* groups to pull/restrain/etc/ */
#  t_pullgrp  *dyna;       /* dynamic groups for use with local constraints */
#  FILE       *out_x;      /* output file for pull data */
#  FILE       *out_f;      /* output file for pull data */

class t_inputrec(Structure):
    _fields_ = [("eI", c_int),
                ("nsteps", c_int),
                ("simulation_part", c_int),
                ("init_step", c_int),
                ("ns_type", c_int),
                ("nstlist", c_int),
                ("ndelta", c_int),
                ("nstcomm", c_int),
                ("comm_mode", c_int),
                ("nstcheckpoint", c_int),
                ("nstlog", c_int),
                ("nstxout", c_int),
                ("nstvout", c_int),
                ("nstfout", c_int),
                ("nstenergy", c_int),
                ("nstxtcout", c_int),
                ("init_t", c_real),
                ("delta_t", c_real),
                ("xtcprec", c_real),
                ("nkx", c_int),
                ("nky", c_int),
                ("nkz", c_int),
                ("pme_order", c_int),
                ("ewald_rtol", c_real),
                ("ewald_geometry", c_int),
                ("epsilon_surface", c_real),
                ("bOptFFT", c_int),
                ("ePBC", c_int),
                ("bPeriodicMols", c_int),
                ("bContinuation", c_int),
                ("etc", c_int),
                ("epc", c_int),
                ("epct", c_int),
                ("tau_p", c_real),
                ("ref_p", tensor),
                ("compress", tensor),
                ("refcoord_scaling", c_int),
                ("posres_com", rvec),
                ("posres_comB", rvec),
                ("andersen_seed", c_int),
                ("rlist", c_real),
                ("rtpi", c_real),
                ("coulombtype", c_int),
                ("rcoulomb_switch", c_real),
                ("rcoulomb", c_real),
                ("epsilon_r", c_real),
                ("epsilon_rf", c_real),
                ("implicit_solvent", c_int),
                ("gb_algorithm", c_int),
                ("nstgbradii", c_int),
                ("rgbradii", c_real),
                ("gb_saltconc", c_real),
                ("gb_epsilon_solvent", c_real),
                ("gb_obc_alpha", c_real),
                ("gb_obc_beta", c_real),
                ("gb_obc_gamma", c_real),
                ("sa_surface_tension", c_real),
                ("vdwtype", c_int),
                ("rvdw_switch", c_real),
                ("rvdw", c_real),
                ("eDispCorr", c_int),
                ("tabext", c_real),
                ("shake_tol", c_real),
                ("efep", c_int),
                ("init_lambda", c_real),
                ("delta_lambda", c_real),
                ("sc_alpha", c_real),
                ("sc_power", c_int),
                ("sc_sigma", c_real),
                ("eDisre", c_int),
                ("dr_fc", c_real),
                ("eDisreWeighting", c_int),
                ("bDisreMixed", c_int),
                ("nstdisreout", c_int),
                ("dr_tau", c_real),
                ("orires_fc", c_real),
                ("orires_tau", c_real),
                ("nstorireout", c_int),
                ("dihre_fc", c_real),
                ("em_stepsize", c_real),
                ("em_tol", c_real),
                ("niter", c_int),
                ("fc_stepsize", c_real),
                ("nstcgsteep", c_int),
                ("nbfgscorr", c_int),
                ("eConstrAlg", c_int),
                ("nProjOrder", c_int),
                ("LincsWarnAngle", c_real),
                ("nLincsIter", c_int),
                ("bShakeSOR", c_int),
                ("bd_fric", c_real),
                ("ld_seed", c_int),
                ("nwall", c_int),
                ("wall_type", c_int),
                ("wall_r_linpot", c_real),
                ("wall_atomtype", c_int*2),
                ("wall_density", c_real*2),
                ("wall_ewald_zfac", c_real),
                ("ePull", c_int),
                ("pull", POINTER(t_pull)),
                ("cos_accel", c_real),
                ("deform", tensor),
                ("userint1", c_int),
                ("userint2", c_int),
                ("userint3", c_int),
                ("userint4", c_int),
                ("userreal1", c_real),
                ("userreal2", c_real),
                ("userreal3", c_real),
                ("userreal4", c_real),
                ("opts", t_grpopts),
                ("ex", t_cosines*3),
                ("et", t_cosines*3),
                ("bQMMM", c_int),
                ("QMconstraints", c_int),
                ("QMMMscheme", c_int),
                ("scalefactor", c_real)]
#  int  eI;              /* Integration method                 */
#  int  nsteps;        /* number of steps to be taken            */
#  int  simulation_part; /* Used in checkpointing to separate chunks */
#  int  init_step;    /* start at a stepcount >0 (used w. tpbconv)    */
#  int  ns_type;        /* which ns method should we use?               */
#  int  nstlist;        /* number of steps before pairlist is generated    */
#  int  ndelta;        /* number of cells per rlong            */
#  int  nstcomm;        /* number of steps after which center of mass    */
#                        /* motion is removed                */
#  int  comm_mode;       /* Center of mass motion removal algorithm      */
#  int nstcheckpoint;    /* checkpointing frequency                      */
#  int nstlog;        /* number of steps after which print to logfile    */
#  int nstxout;        /* number of steps after which X is output    */
#  int nstvout;        /* id. for V                    */
#  int nstfout;        /* id. for F                    */
#  int nstenergy;    /* number of steps after which energies printed */
#  int nstxtcout;    /* id. for compressed trj (.xtc)        */
#  real init_t;        /* initial time (ps)                 */
#  real delta_t;        /* time step (ps)                */
#  real xtcprec;         /* precision of xtc file                        */
#  int  nkx,nky,nkz;     /* number of k vectors in each spatial dimension*/
#                        /* for fourier methods for long range electrost.*/
#  int  pme_order;       /* interpolation order for PME                  */
#  real ewald_rtol;      /* Real space tolerance for Ewald, determines   */
#                        /* the real/reciprocal space relative weight    */
#  int  ewald_geometry;  /* normal/3d ewald, or pseudo-2d LR corrections */
#  real epsilon_surface; /* Epsilon for PME dipole correction            */
#  bool bOptFFT;         /* optimize the fft plan at start               */
#  int  ePBC;        /* Type of periodic boundary conditions        */
#  int  bPeriodicMols;   /* Periodic molecules                           */
#  bool bContinuation;   /* Continuation run: starting state is correct    */
#  int  etc;        /* temperature coupling                 */
#  int  epc;        /* pressure coupling                            */
#  int  epct;        /* pressure coupling type            */
#  real tau_p;        /* pressure coupling time (ps)            */
#  tensor ref_p;        /* reference pressure (kJ/(mol nm^3))        */
#  tensor compress;    /* compressability ((mol nm^3)/kJ)         */
#  int  refcoord_scaling;/* How to scale absolute reference coordinates  */
#  rvec posres_com;      /* The COM of the posres atoms                  */
#  rvec posres_comB;     /* The B-state COM of the posres atoms          */
#  int  andersen_seed;   /* Random seed for Andersen thermostat.         */
#  real rlist;        /* short range pairlist cut-off (nm)        */
#  real rtpi;            /* Radius for test particle insertion           */
#  int  coulombtype;    /* Type of electrostatics treatment             */
#  real rcoulomb_switch; /* Coulomb switch range start (nm)        */
#  real rcoulomb;        /* Coulomb cutoff (nm)                        */
#  real epsilon_r;       /* relative dielectric constant                 */
#  real epsilon_rf;      /* relative dielectric constant of the RF       */
#  int  implicit_solvent;/* No (=explicit water), or GBSA solvent models */
#  int  gb_algorithm;    /* Algorithm to use for calculation Born radii  */
#  int  nstgbradii;      /* Frequency of updating Generalized Born radii */
#  real rgbradii;        /* Cutoff for GB radii calculation              */
#  real gb_saltconc;     /* Salt concentration (M) for GBSA models       */
#  real gb_epsilon_solvent; /* dielectric coeff. of implicit solvent     */
#  real gb_obc_alpha;    /* 1st scaling factor for Bashford-Case GB      */
#  real gb_obc_beta;     /* 2nd scaling factor for Bashford-Case GB      */
#  real gb_obc_gamma;    /* 3rd scaling factor for Bashford-Case GB      */
#  real sa_surface_tension; /* Energy factor for SA part of GBSA */
#  int  vdwtype;         /* Type of Van der Waals treatment              */
#  real rvdw_switch;     /* Van der Waals switch range start (nm)        */
#  real rvdw;            /* Van der Waals cutoff (nm)                */
#  int  eDispCorr;       /* Perform Long range dispersion corrections    */
#  real tabext;          /* Extension of the table beyond the cut-off,   *
#                          * as well as the table length for 1-4 interac. */
#  real shake_tol;    /* tolerance for shake                */
#  int  efep;           /* free energy interpolation no/yes        */
#  real init_lambda;    /* initial value for perturbation variable    */
#  real delta_lambda;    /* change of lambda per time step (1/dt)    */
#  real sc_alpha;        /* free energy soft-core parameter              */
#  int  sc_power;        /* lambda power for soft-core interactions      */
#  real sc_sigma;        /* free energy soft-core sigma when c6 or c12=0 */
#  int  eDisre;          /* Type of distance restraining                 */
#  real dr_fc;            /* force constant for ta_disre            */
#  int  eDisreWeighting; /* type of weighting of pairs in one restraints    */
#  bool bDisreMixed;     /* Use comb of time averaged and instan. viol's    */
#  int  nstdisreout;     /* frequency of writing pair distances to enx   */
#  real dr_tau;            /* time constant for memory function in disres     */
#  real orires_fc;        /* force constant for orientational restraints  */
#  real orires_tau;        /* time constant for memory function in orires     */
#  int  nstorireout;     /* frequency of writing tr(SD) to enx           */
#  real dihre_fc;        /* force constant for dihedral restraints    */
#  real em_stepsize;        /* The stepsize for updating            */
#  real em_tol;            /* The tolerance                */
#  int  niter;           /* Number of iterations for convergence of      */
#                        /* steepest descent in relax_shells             */
#  real fc_stepsize;     /* Stepsize for directional minimization        */
#                        /* in relax_shells                              */
#  int  nstcgsteep;      /* number of steps after which a steepest       */
#                        /* descents step is done while doing cg         */
#  int  nbfgscorr;       /* Number of corrections to the hessian to keep */
#  int  eConstrAlg;      /* Type of constraint algorithm                 */
#  int  nProjOrder;      /* Order of the LINCS Projection Algorithm      */
#  real LincsWarnAngle;  /* If bond rotates more than %g degrees, warn   */
#  int  nLincsIter;      /* Number of iterations in the final Lincs step */
#  bool bShakeSOR;       /* Use successive overrelaxation for shake      */
#  real bd_fric;         /* Friction coefficient for BD (amu/ps)         */
#  int  ld_seed;         /* Random seed for SD and BD                    */
#  int  nwall;           /* The number of walls                          */
#  int  wall_type;       /* The type of walls                            */
#  real wall_r_linpot;   /* The potentail is linear for r<=wall_r_linpot */
#  int  wall_atomtype[2];/* The atom type for walls                      */
#  real wall_density[2]; /* Number density for walls                     */
#  real wall_ewald_zfac; /* Scaling factor for the box for Ewald         */
#  int  ePull;           /* Type of pulling: no, umbrella or constraint  */
#  t_pull *pull;         /* The data for center of mass pulling          */
#  real cos_accel;       /* Acceleration for viscosity calculation       */
#  tensor deform;        /* Triclinic deformation velocities (nm/ps)     */
#  int  userint1;        /* User determined parameters                   */
#  int  userint2;
#  int  userint3;
#  int  userint4;
#  real userreal1;
#  real userreal2;
#  real userreal3;
#  real userreal4;
#  t_grpopts opts;    /* Group options                */
#  t_cosines ex[DIM];    /* Electric field stuff    (spatial part)        */
#  t_cosines et[DIM];    /* Electric field stuff    (time part)        */
#  bool bQMMM;           /* QM/MM calculation                            */
#  int  QMconstraints;   /* constraints on QM bonds                      */
#  int  QMMMscheme;      /* Scheme: ONIOM or normal                      */
#  real scalefactor;     /* factor for scaling the MM charges in QM calc.*/

# enum {
(
  estLAMBDA,
  estBOX,
  estBOX_REL,
  estBOXV,
  estPRES_PREV,
  estNH_XI,
  estTC_INT,
  estX,
  estV,
  estSDX,
  estCGP,
  estLD_RNG,
  estLD_RNGI,
  estDISRE_INITF,
  estDISRE_RM3TAV,
  estORIRE_INITF,
  estORIRE_DTAV,
  estNR) = map(c_int,xrange(18))
# };
#/* These enums are used in flags as (1<<est...).
# * The order of these enums should not be changed,
# * since that affects the checkpoint (.cpt) file format.
# */

est_names = POINTER(c_char_p) * estNR.value
# extern const char *est_names[estNR];
#/* The names of the state entries, defined in src/gmxib/checkpoint.c */

class history_t(Structure):
    _fields_ = [("disre_initf", c_real),
                ("ndisrepairs", c_int),
                ("disre_rm3tav", POINTER(c_real)),
                ("orire_initf", c_real),
                ("norire_Dtav", c_int),
                ("orire_Dtav", POINTER(c_real))]
#  real disre_initf;    /* The scaling factor for initializing the time av. */
#  int  ndisrepairs;    /* The number of distance restraints                */
#  real *disre_rm3tav;  /* The r^-3 time averaged pair distances            */
#  real orire_initf;    /* The scaling factor for initializing the time av. */
#  int  norire_Dtav;    /* The number of matrix element in dtav (npair*5)   */
#  real *orire_Dtav;    /* The time averaged orientation tensors            */

class ekinstate_t(Structure):
    _fields_ = [("bUpToDate", c_int),
                ("ekinh_n", c_int),
                ("ekinh", POINTER(tensor)),
                ("dekindl", c_real),
                ("mvcos", c_real)]
#/* Struct used for checkpointing only.
# * This struct would not be required with unlimited precision.
# * But because of limited precision, the COM motion removal implementation
# * can cause the kinetic energy in the MD loop to differ by a few bits from
# * the kinetic energy one would determine from state.v.
# */
#  bool     bUpToDate;
#  int      ekinh_n;
#  tensor * ekinh;
#  real     dekindl;
#  real     mvcos;

class energyhistory_t(Structure):
    _fields_ = [("ener_ave", POINTER(c_double)),
                ("ener_sum", POINTER(c_double)),
                ("nener", c_int)]
#  double *   ener_ave;  /* Energy term history sum to get fluctuations      */
#  double *   ener_sum;  /* Energy term history sum to get fluctuations      */
#  int        nener;     /* Number of energy terms in two previous arrays    */

class t_state(Structure):
    _fields_ = [("natoms", c_int),
                ("ngtc", c_int),
                ("nrng", c_int),
                ("nrngi", c_int),
                ("flags", c_int),
                ("Lambda", c_real),
                ("box", matrix),
                ("box_rel", matrix),
                ("boxv", matrix),
                ("pres_prev", matrix),
                ("nosehoover_xi", POINTER(c_real)),
                ("therm_integral", POINTER(c_double)),
                ("nalloc", c_int),
                ("x", POINTER(rvec)),
                ("v", POINTER(rvec)),
                ("sd_X", POINTER(rvec)),
                ("cg_p", POINTER(rvec)),
                ("ld_rng", POINTER(c_uint)),
                ("ld_rngi", POINTER(c_int)),
                ("hist", history_t),
                ("ekinstate", ekinstate_t),
                ("enerhist", energyhistory_t),
                ("ddp_count", c_int),
                ("ddp_count_cg_gl", c_int),
                ("ncg_gl", c_int),
                ("cg_gl", POINTER(c_int)),
                ("cg_gl_nalloc", c_int)]
#  int           natoms;
#  int           ngtc;
#  int           nrng;
#  int           nrngi;
#  int           flags;  /* Flags telling which entries are present      */
#  real          lambda; /* the free energy switching parameter          */
#  matrix         box;    /* box vector coordinates                          */
#  matrix         box_rel; /* Relitaive box vectors to preserve shape        */
#  matrix         boxv;   /* box velocitites for Parrinello-Rahman pcoupl */
#  matrix        pres_prev; /* Pressure of the previous step for pcoupl  */
#  real          *nosehoover_xi;  /* for Nose-Hoover tcoupl (ngtc)       */
#  double        *therm_integral; /* for N-H/V-rescale tcoupl (ngtc)     */
#  int           nalloc; /* Allocation size for x, v and sd_x when !=NULL*/
#  rvec          *x;     /* the coordinates (natoms)                     */
#  rvec          *v;     /* the velocities (natoms)                      */
#  rvec          *sd_X;  /* random part of the x update for stoch. dyn.  */
#  rvec          *cg_p;  /* p vector for conjugate gradient minimization */
#
#  unsigned int  *ld_rng;  /* RNG random state                           */
#  int           *ld_rngi; /* RNG index                                  */
#
#  history_t     hist;   /* Time history for restraints                  */
#
#  ekinstate_t   ekinstate; /* The state of the kinetic energy data      */
#
#  energyhistory_t  enerhist; /* Energy history for statistics           */
#
#  int           ddp_count; /* The DD partitioning count for this state  */
#  int           ddp_count_cg_gl; /* The DD part. count for index_gl     */
#  int           ncg_gl; /* The number of local charge groups            */
#  int           *cg_gl; /* The global cg number of the local cgs        */
#  int           cg_gl_nalloc; /* Allocation size of cg_gl;              */


#/* The nonbonded kernels are documented in gmxlib/nonbonded_kernels,
# * but here's a lazy version of the numbering. The first position
# * is the Coulomb interaction (0 for none), second is Van der Waals
# * (again, 0 means no interaction), and the third is the water optimization
# * (0 meaning no water optimization = standard atom-atom loop)
# *
# *                                     value
# * pos                 1                   2           3              4
# * 1st Coul        Normal,1/r       Reaction-field  Table            Generalized born
# * 2nd Vdw         Lennard-Jones    Buckingham      Table             n/a
# * 3rd Water. opt  SPC-other atom   SPC-SPC         TIP4p-other at.  TIP4p-TIP4p
# */

eNR_NBKERNEL_NONE = -1
#define eNR_NBKERNEL_NONE -1

#enum{
(
  eNR_NBKERNEL010,
  eNR_NBKERNEL020,
  eNR_NBKERNEL030,
  eNR_NBKERNEL100,
  eNR_NBKERNEL101,
  eNR_NBKERNEL102,
  eNR_NBKERNEL103,
  eNR_NBKERNEL104,
  eNR_NBKERNEL110,
  eNR_NBKERNEL111,
  eNR_NBKERNEL112,
  eNR_NBKERNEL113,
  eNR_NBKERNEL114,
  eNR_NBKERNEL120,
  eNR_NBKERNEL121,
  eNR_NBKERNEL122,
  eNR_NBKERNEL123,
  eNR_NBKERNEL124,
  eNR_NBKERNEL130,
  eNR_NBKERNEL131,
  eNR_NBKERNEL132,
  eNR_NBKERNEL133,
  eNR_NBKERNEL134,
  eNR_NBKERNEL200,
  eNR_NBKERNEL201,
  eNR_NBKERNEL202,
  eNR_NBKERNEL203,
  eNR_NBKERNEL204,
  eNR_NBKERNEL210,
  eNR_NBKERNEL211,
  eNR_NBKERNEL212,
  eNR_NBKERNEL213,
  eNR_NBKERNEL214,
  eNR_NBKERNEL220,
  eNR_NBKERNEL221,
  eNR_NBKERNEL222,
  eNR_NBKERNEL223,
  eNR_NBKERNEL224,
  eNR_NBKERNEL230,
  eNR_NBKERNEL231,
  eNR_NBKERNEL232,
  eNR_NBKERNEL233,
  eNR_NBKERNEL234,
  eNR_NBKERNEL300,
  eNR_NBKERNEL301,
  eNR_NBKERNEL302,
  eNR_NBKERNEL303,
  eNR_NBKERNEL304,
  eNR_NBKERNEL310,
  eNR_NBKERNEL311,
  eNR_NBKERNEL312,
  eNR_NBKERNEL313,
  eNR_NBKERNEL314,
  eNR_NBKERNEL320,
  eNR_NBKERNEL321,
  eNR_NBKERNEL322,
  eNR_NBKERNEL323,
  eNR_NBKERNEL324,
  eNR_NBKERNEL330,
  eNR_NBKERNEL331,
  eNR_NBKERNEL332,
  eNR_NBKERNEL333,
  eNR_NBKERNEL334,
  eNR_NBKERNEL400,
  eNR_NBKERNEL410,
  eNR_NBKERNEL430,
  eNR_NBKERNEL010NF,
  eNR_NBKERNEL020NF,
  eNR_NBKERNEL030NF,
  eNR_NBKERNEL100NF,
  eNR_NBKERNEL101NF,
  eNR_NBKERNEL102NF,
  eNR_NBKERNEL103NF,
  eNR_NBKERNEL104NF,
  eNR_NBKERNEL110NF,
  eNR_NBKERNEL111NF,
  eNR_NBKERNEL112NF,
  eNR_NBKERNEL113NF,
  eNR_NBKERNEL114NF,
  eNR_NBKERNEL120NF,
  eNR_NBKERNEL121NF,
  eNR_NBKERNEL122NF,
  eNR_NBKERNEL123NF,
  eNR_NBKERNEL124NF,
  eNR_NBKERNEL130NF,
  eNR_NBKERNEL131NF,
  eNR_NBKERNEL132NF,
  eNR_NBKERNEL133NF,
  eNR_NBKERNEL134NF,
  eNR_NBKERNEL200NF,
  eNR_NBKERNEL201NF,
  eNR_NBKERNEL202NF,
  eNR_NBKERNEL203NF,
  eNR_NBKERNEL204NF,
  eNR_NBKERNEL210NF,
  eNR_NBKERNEL211NF,
  eNR_NBKERNEL212NF,
  eNR_NBKERNEL213NF,
  eNR_NBKERNEL214NF,
  eNR_NBKERNEL220NF,
  eNR_NBKERNEL221NF,
  eNR_NBKERNEL222NF,
  eNR_NBKERNEL223NF,
  eNR_NBKERNEL224NF,
  eNR_NBKERNEL230NF,
  eNR_NBKERNEL231NF,
  eNR_NBKERNEL232NF,
  eNR_NBKERNEL233NF,
  eNR_NBKERNEL234NF,
  eNR_NBKERNEL300NF,
  eNR_NBKERNEL301NF,
  eNR_NBKERNEL302NF,
  eNR_NBKERNEL303NF,
  eNR_NBKERNEL304NF,
  eNR_NBKERNEL310NF,
  eNR_NBKERNEL311NF,
  eNR_NBKERNEL312NF,
  eNR_NBKERNEL313NF,
  eNR_NBKERNEL314NF,
  eNR_NBKERNEL320NF,
  eNR_NBKERNEL321NF,
  eNR_NBKERNEL322NF,
  eNR_NBKERNEL323NF,
  eNR_NBKERNEL324NF,
  eNR_NBKERNEL330NF,
  eNR_NBKERNEL331NF,
  eNR_NBKERNEL332NF,
  eNR_NBKERNEL333NF,
  eNR_NBKERNEL334NF,
  eNR_NBKERNEL400NF,
  eNR_NBKERNEL410NF,
  eNR_NBKERNEL430NF,
  eNR_NBKERNEL_NR,
  eNR_NBKERNEL_OUTER,
  eNR_NB14,
  eNR_WEIGHTS,
  eNR_SPREADQ,
  eNR_SPREADQBSP,
  eNR_GATHERF,
  eNR_GATHERFBSP,
  eNR_FFT,
  eNR_CONV,
  eNR_SOLVEPME,eNR_NS,
  eNR_RESETX,
  eNR_SHIFTX,
  eNR_CGCM,
  eNR_FSUM,
  eNR_BONDS,
  eNR_G96BONDS,
  eNR_FENEBONDS,
  eNR_TABBONDS,
  eNR_ANGLES,
  eNR_G96ANGLES,
  eNR_QANGLES,
  eNR_TABANGLES,
  eNR_PROPER,
  eNR_IMPROPER,
  eNR_RB,
  eNR_FOURDIH,
  eNR_TABDIHS,
  eNR_DISRES,
  eNR_ORIRES,
  eNR_DIHRES,
  eNR_POSRES,
  eNR_ANGRES,
  eNR_ANGRESZ,
  eNR_MORSE,
  eNR_CUBICBONDS,
  eNR_WALLS,
  eNR_WPOL,
  eNR_THOLE,
  eNR_VIRIAL,
  eNR_UPDATE,
  eNR_EXTUPDATE,
  eNR_STOPCM,
  eNR_PCOUPL,
  eNR_EKIN,
  eNR_LINCS,
  eNR_LINCSMAT,
  eNR_SHAKE,
  eNR_CONSTR_V,
  eNR_SHAKE_RIJ,
  eNR_CONSTR_VIR,
  eNR_SETTLE,
  eNR_VSITE2,
  eNR_VSITE3,
  eNR_VSITE3FD,
  eNR_VSITE3FAD,
  eNR_VSITE3OUT,
  eNR_VSITE4FD,
  eNR_VSITE4FDN,
  eNR_VSITEN,
  eNRNB) = map(c_int,xrange(194))
eNR_NBKERNEL_FREE_ENERGY = eNR_NBKERNEL_NR
#};


class t_nrnb(Structure):
    _fields_ = [("n",c_double * eNRNB.value)]
#  double n[eNRNB];

gmx_cycles_t = c_ulonglong
#from gmx_cyclecounter.h:
#typedef unsigned long long
#gmx_cycles_t;

class gmx_wallcycle(Structure):
    _fields_ = [("n", c_int),
                ("c", gmx_cycles_t),
                ("start", gmx_cycles_t),
                ("last", gmx_cycles_t)]
#from gmx_wallcycle.c:
#typedef struct gmx_wallcycle {
#  int          n;
#  gmx_cycles_t c;
#  gmx_cycles_t start;
#  gmx_cycles_t last;
#} gmx_wallcycle_t_t;

gmx_wallcycle_t = POINTER(gmx_wallcycle)
#typedef struct gmx_wallcycle *gmx_wallcycle_t;

(
  egcTC,
  egcENER,
  egcACC,
  egcFREEZE,
  egcUser1,
  egcUser2,
  egcVCM,
  egcXTC,
  egcORFIT,
  egcQMMM,
  egcNR
) = map(c_int,xrange(11))
#enum {
#  egcTC,    egcENER,   egcACC, egcFREEZE,
#  egcUser1, egcUser2,  egcVCM, egcXTC,
#  egcORFIT, egcQMMM,
#  egcNR
#};

class gmx_groups_t(Structure):
    _fields_ =[("grps", t_grps * egcNR.value),
               ("ngrpname", c_int),
               ("grpname", POINTER(POINTER(c_char_p))),
               ("ngrpnr", c_int * egcNR.value),
               ("grpnr", POINTER(c_ubyte) * egcNR.value)]
#  t_grps         grps[egcNR];   /* Groups of things                     */
#  int            ngrpname;      /* Number of groupnames                 */
#  char           ***grpname;    /* Names of the groups                */
#  int            ngrpnr[egcNR];
#  unsigned char  *grpnr[egcNR]; /* Group numbers or NULL        */

class gmx_molblock_t(Structure):
    _fields_ = [("type", c_int),
                ("nmol", c_int),
                ("natoms_mol", c_int),
                ("nposres_xA", c_int),
                ("posres_xA", POINTER(rvec)),
                ("nposres_xB", c_int),
                ("posres_xB", POINTER(rvec))]
#  int           type;           /* The molcule type index in mtop.moltype */
#  int           nmol;           /* The number of molecules in this block  */
#  int           natoms_mol;     /* The number of atoms in one molecule    */
#  int           nposres_xA;     /* The number of posres coords for top A  */
#  rvec          *posres_xA;     /* The posres coords for top A            */
#  int           nposres_xB;     /* The number of posres coords for top B  */
#  rvec          *posres_xB;     /* The posres coords for top B            */

class gmx_moltype_t(Structure):
    _fields_ = [("name", POINTER(c_char_p)),
                ("atoms", t_atoms),
                ("ilist", t_ilist * F_NRE.value),
                ("cgs", t_block),
                ("excls", t_blocka)]
#  char          **name;         /* Name of the molecule type              */
#  t_atoms    atoms;        /* The atoms                           */
#  t_ilist       ilist[F_NRE];
#  t_block       cgs;            /* The charge groups                    */
#  t_blocka      excls;          /* The exclusions                       */

class gmx_ffparams_t(Structure):
    _fields_ =[("ntypes", c_int),
               ("atnr", c_int),
               ("functype", POINTER(t_functype)),
               ("iparams", POINTER(t_iparams)),
               ("fudgeQQ", c_real)]
#  int        ntypes;
#  int        atnr;
#  t_functype *functype;
#  t_iparams  *iparams;
#  real       fudgeQQ;

class gmx_mtop_t(Structure):
    _fields_ = [("name", POINTER(c_char_p)),
                ("ffparams", gmx_ffparams_t),
                ("nmoltype", c_int),
                ("moltype", POINTER(gmx_moltype_t)),
                ("nmolblock", c_int),
                ("molblock", POINTER(gmx_molblock_t)),
                ("natoms", c_int),
                ("atomtypes", t_atomtypes),
                ("mols", t_block),
                ("groups", gmx_groups_t),
                ("symtab", t_symtab)]
#/* The global topology struct, based on molecule types */
#  char           **name;    /* Name of the topology                       */
#  gmx_ffparams_t ffparams;
#  int            nmoltype;
#  gmx_moltype_t  *moltype;
#  int            nmolblock;
#  gmx_molblock_t *molblock;
#  int            natoms;
#  t_atomtypes    atomtypes;     /* Atomtype properties                  */
#  t_block        mols;          /* The molecules                        */
#  gmx_groups_t   groups;
#  t_symtab     symtab;        /* The symbol table            */

class t_mdatoms(Structure):
    _fields_ = [("tmassA", c_real),
                ("tmassB", c_real),
                ("tmass", c_real),
                ("nr", c_int),
                ("nalloc", c_int),
                ("nenergrp", c_int),
                ("bVCMgrps", c_int),
                ("nPerturbed", c_int),
                ("nMassPerturbed", c_int),
                ("nChargePerturbed", c_int),
                ("bOrires", c_int),
                ("massA", POINTER(c_real)),
                ("massB", POINTER(c_real)),
                ("massT", POINTER(c_real)),
                ("invmass", POINTER(c_real)),
                ("chargeA", POINTER(c_real)),
                ("chargeB", POINTER(c_real)),
                ("bPerturbed", POINTER(c_int)),
                ("typeA", POINTER(c_int)),
                ("typeB", POINTER(c_int)),
                ("ptype", POINTER(c_ushort)),
                ("cTC", POINTER(c_ushort)),
                ("cENER", POINTER(c_ushort)),
                ("cACC", POINTER(c_ushort)),
                ("cFREEZE", POINTER(c_ushort)),
                ("cVCM", POINTER(c_ushort)),
                ("cU1", POINTER(c_ushort)),
                ("cU2", POINTER(c_ushort)),
                ("cORF", POINTER(c_ushort)),
                ("bQM", POINTER(c_int)),
                ("start", c_int),
                ("homenr", c_int),
                ("Lambda", c_real)]
#  real          tmassA,tmassB,tmass;
#  int           nr;
#  int           nalloc;
#  int           nenergrp;
#  bool          bVCMgrps;
#  int           nPerturbed;
#  int           nMassPerturbed;
#  int           nChargePerturbed;
#  bool          bOrires;
#  real          *massA,*massB,*massT,*invmass;
#  real          *chargeA,*chargeB;
#  bool          *bPerturbed;
#  int           *typeA,*typeB;
#  unsigned short        *ptype;
#  unsigned short        *cTC,*cENER,*cACC,*cFREEZE,*cVCM;
#  unsigned short        *cU1,*cU2,*cORF;
#  /* for QMMM, atomnumber contains atomic number of the atoms */
#  bool          *bQM;
#  /* The range of home atoms */
#  int           start;
#  int           homenr;
#  /* The lambda value used to create the contents of the struct */
#  real          lambda;

(
  eNL_VDWQQ,
  eNL_VDW,
  eNL_QQ,
  eNL_VDWQQ_FREE,
  eNL_VDW_FREE,
  eNL_QQ_FREE,
  eNL_VDWQQ_WATER,
  eNL_QQ_WATER,
  eNL_VDWQQ_WATERWATER,
  eNL_QQ_WATERWATER,
  eNL_NR
) = map(c_int,xrange(11))
#enum { eNL_VDWQQ, eNL_VDW, eNL_QQ,
#       eNL_VDWQQ_FREE, eNL_VDW_FREE, eNL_QQ_FREE,
#       eNL_VDWQQ_WATER, eNL_QQ_WATER,
#       eNL_VDWQQ_WATERWATER, eNL_QQ_WATERWATER,
#       eNL_NR };

class t_nblist(Structure):
    _fields_ = [("il_code", c_int),
                ("icoul", c_int),
                ("ivdw", c_int),
                ("free_energy", c_int),
                ("nltype", c_int),
                ("nri", c_int),
                ("maxnri", c_int),
                ("nrj", c_int),
                ("maxnrj", c_int),
                ("maxlen", c_int),
                ("iinr", POINTER(c_int)),
                ("gid", POINTER(c_int)),
                ("shift", POINTER(c_int)),
                ("jindex", POINTER(c_int)),
                ("jjnr", POINTER(c_int)),
                ("count", c_int),
                ("mtx", c_void_p)]
#  int             il_code;      /* Innerloop index from nrnb.h, used     */
#                                /* for flop accounting.                  */
#  int             icoul;        /* Coulomb loop type index for kernels   */
#  int             ivdw;         /* VdW loop type index for kernels       */
#  int             free_energy;  /* Free energy setting for this list     */
#  int             nltype;       /* Atom, water, or water-water list      */
#
#  int             nri,maxnri;   /* Current/max number of i particles     */
#  int             nrj,maxnrj;   /* Current/max number of j particles     */
#  int             maxlen;       /* maxnr of j atoms for a single i atom  */
#  int *           iinr;            /* The i-elements                    */
#  int *           gid;          /* Index in energy arrays                */
#  int *           shift;        /* Shift vector index                    */
#  int *           jindex;       /* Index in jjnr                         */
#  int *           jjnr;            /* The j-atom list                       */
#  int             count;        /* counter to multithread the innerloops */
#  void *          mtx;          /* mutex to lock the counter             */

class t_forcetable(Structure):
    _fields_ =[("r", c_real),
               ("n", c_int),
               ("scale", c_real),
               ("scale_exp", c_real),
               ("tab", POINTER(c_real))]
#  real r;         /* range of the table */
#  int  n;         /* n+1 is the number of points */
#  real scale;     /* distance between two points */
#  real scale_exp; /* distance for exponential Buckingham table */
#  real *tab;      /* the actual tables, per point there are  4 numbers for
#           * Coulomb, dispersion and repulsion (in total 12 numbers)
#           */

class t_nblists(Structure):
    _fields_ = [("tab", t_forcetable),
                ("coultab", POINTER(c_real)),
                ("vdwtab", POINTER(c_real)),
                ("nlist_sr", t_nblist * eNL_NR.value),
                ("nlist_lr", t_nblist * eNL_NR.value)]
#  t_forcetable tab;
#  /* We duplicate tables for cache optimization purposes */
#  real *coultab;      /* Coul only */
#  real *vdwtab;       /* Vdw only   */
#  /* The actual neighbor lists, short and long range, see enum above
#   * for definition of neighborlist indices.
#   */
#  t_nblist nlist_sr[eNL_NR];
#  t_nblist nlist_lr[eNL_NR];


#From ./src/mdlib/gmx_fft_fftw3.c:
if environ.has_key("GROMPYDOUBLE"):
    def FFTWPREFIX(name):
        return "fftw_"+name
else:
    def FFTWPREFIX(name):
        return "fftwf_"+name
#ifdef GMX_DOUBLE
#define FFTWPREFIX(name) fftw_ ## name
#else
#define FFTWPREFIX(name) fftwf_ ## name
#endif

class gmx_fft(Structure):
    _fields_ = [("plan", c_void_p *2*2*2), # [2][2][2] array
                # ("plan", FFTWPREFIX("plan")), # [2][2][2] array
                ("real_transform", c_int),
                ("ndim", c_int)]
#From ./src/mdlib/gmx_fft_fftw3.c:
#struct gmx_fft
#{
#    /* Three alternatives (unaligned/aligned, out-of-place/in-place, forward/backward)
#     * results in 8 different FFTW plans. Keep track of them with 3 array indices:
#     * first index:   0=unaligned, 1=aligned
#     * second index:  0=out-of-place, 1=in-place
#     * third index:   0=backward, 1=forward
#     */
#    FFTWPREFIX(plan)         plan[2][2][2];
#    /* Catch user mistakes */
#    int                      real_transform;
#    int                      ndim;
#};

gmx_fft_t = POINTER(gmx_fft)
#typedef struct gmx_domdec_comm *gmx_domdec_comm_p_t;
#From gmx_fft.h:
#/*! \brief Datatype for FFT setup
# *
# *  The gmx_fft_t type contains all the setup information, e.g. twiddle
# *  factors, necessary to perform an FFT. Internally it is mapped to
# *  whatever FFT library we are using, or the built-in FFTPACK if no fast
# *  external library is available.
# *
# *  Since some of the libraries (e.g. MKL) store work array data in their
# *  handles this datatype should only be used for one thread at a time, i.e.
# *  they should allocate one instance each when executing in parallel.
# */
#typedef struct gmx_fft *gmx_fft_t;

class t_fftgrid(Structure):
    _fields_ = [("ptr", POINTER(c_real)),
                ("bParallel", c_int),
                ("workspace", POINTER(c_real)),
                ("nx", c_int),
                ("ny", c_int),
                ("nz", c_int),
                ("la2r", c_int),
                ("la2c", c_int),
                ("la12r", c_int),
                ("la12c", c_int),
                ("nptr", c_int),
                ("nxyz", c_int),
                ("fft_setup", gmx_fft_t)
             ##ifdef GMX_MPI
             #  ("mpi_fft_setup", gmx_parallel_3dfft_t),
             #  ("pfft", t_parfft)
             ##endif
               ]
#    real *                 ptr;
#    bool                   bParallel;
#    real *                 workspace;
#    int                    nx,ny,nz,la2r,la2c,la12r,la12c;
#    int                    nptr,nxyz;
#    gmx_fft_t              fft_setup;
##ifdef GMX_MPI
#    gmx_parallel_3dfft_t   mpi_fft_setup;
#    t_parfft               pfft;
##endif

class pme_grid_comm_t(Structure):
    _fields_ = [("snd0", c_int),
                ("snds", c_int),
                ("rcv0", c_int),
                ("rcvs", c_int)]
#From ./src/mdlib/pme.c:
#/* Internal datastructures */
#    int snd0;
#    int snds;
#    int rcv0;
#    int rcvs;

class pme_overlap_t(Structure):
    _fields_ = [
             ##ifdef GMX_MPI
             #  ("mpi_comm",MPI_Comm),
             ##endif
                ("nslab", c_int),
                ("s2g", POINTER(c_int)),
                ("nleftbnd", c_int),
                ("nrightbnd", c_int),
                ("nodeid", c_int),
                ("leftid", POINTER(c_int)),
                ("rightid", POINTER(c_int)),
                ("leftc", POINTER(pme_grid_comm_t)),
                ("rightc", POINTER(pme_grid_comm_t))]
#From ./src/mdlib/pme.c:
#typedef struct {
##ifdef GMX_MPI
#    MPI_Comm mpi_comm;
##endif
#    int  nslab;
#    int  *s2g;
#    int  nleftbnd,nrightbnd;  /* The number of nodes to communicate with */
#    int  nodeid,*leftid,*rightid;
#    pme_grid_comm_t *leftc,*rightc;
#} pme_overlap_t;

class pme_atomcomm_t(Structure):
    _fields_ = [
           ##ifdef GMX_MPI
           #    ("mpi_comm", MPI_Comm),
           ##endif
                ("nslab", c_int),
                ("nodeid", c_int),
                ("node_dest", POINTER(c_int)),
                ("node_src", POINTER(c_int)),
                ("buf_index", POINTER(c_int)),
                ("npd", c_int),
                ("pd_nalloc", c_int),
                ("pd", POINTER(c_int)),
                ("count", POINTER(c_int)),
                ("rcount", POINTER(c_int)),
                ("n", c_int),
                ("nalloc", c_int),
                ("x", POINTER(rvec)),
                ("q", POINTER(c_real)),
                ("f", POINTER(rvec)),
                ("bSpread", c_int),
                ("pme_order", c_int),
                ("theta", splinevec),
                ("dtheta", splinevec),
                ("idx", POINTER(ivec)),
                ("fractx", POINTER(rvec))]
#From ./src/mdlib/pme.c:
##ifdef GMX_MPI
#    MPI_Comm mpi_comm;
##endif
#    int  nslab;
#    int  nodeid;
#
#    int  *node_dest;        /* The nodes to send x and q to with DD */
#    int  *node_src;         /* The nodes to receive x and q from with DD */
#    int  *buf_index;        /* Index for commnode into the buffers */
#
#    int  npd;
#    int  pd_nalloc;
#    int  *pd;
#    int  *count;            /* The number of atoms to send to each node */
#    int  *rcount;           /* The number of atoms to receive */
#
#    int  n;
#    int  nalloc;
#    rvec *x;
#    real *q;
#    rvec *f;
#    bool bSpread;           /* These coordinates are used for spreading */
#    int  pme_order;
#    splinevec theta,dtheta;
#    ivec *idx;
#    rvec *fractx;            /* Fractional coordinate relative to the
#                              * lower cell boundary
#                              */

class gmx_pme(Structure):
    _fields_ = [("ndecompdim", c_int),
                ("nodeid", c_int),
                ("nnodes", c_int),
            # #ifdef GMX_MPI
            #   ("mpi_comm", MPI_Comm),
            # #endif
                ("bPPnode", c_int),
                ("bFEP", c_int),
                ("nkx", c_int),
                ("nky", c_int),
                ("nkz", c_int),
                ("pme_order", c_int),
                ("epsilon_r", c_real),
                ("gridA", POINTER(t_fftgrid)),
                ("gridB", POINTER(t_fftgrid)),
                ("nnx", POINTER(c_int)),
                ("nny", POINTER(c_int)),
                ("nnz", POINTER(c_int)),
                ("atc", pme_atomcomm_t*2),
                ("recipbox", matrix),
                ("bsp_mod", splinevec),
                ("overlap", pme_overlap_t*2),
                ("bufv", POINTER(rvec)),
                ("bufr", POINTER(c_real)),
                ("buf_nalloc", c_int),
                ("maxkz", c_int),
                ("work_mhz", POINTER(c_real)),
                ("work_m2", POINTER(c_real)),
                ("work_denom", POINTER(c_real)),
                ("work_tmp1", POINTER(c_real)),
                ("work_m2inv", POINTER(c_real))]
#from /home/rpool/src/gromacs-4.0.5/src/mdlib/pme.c:
#    int  ndecompdim",         /* The number of decomposition dimensions */
#    int  nodeid;             /* Our nodeid in mpi->mpi_comm */
#    int  nnodes;             /* The number of nodes doing PME */
##ifdef GMX_MPI
#    MPI_Comm mpi_comm;
##endif
#
#    bool bPPnode;            /* Node also does particle-particle forces */
#    bool bFEP;               /* Compute Free energy contribution */
#    int nkx,nky,nkz;         /* Grid dimensions */
#    int pme_order;
#    real epsilon_r;
#    t_fftgrid *gridA,*gridB;
#    int  *nnx,*nny,*nnz;
#
#    pme_atomcomm_t atc[2];
#    matrix    recipbox;
#    splinevec bsp_mod;
#
#    pme_overlap_t overlap[2];
#
#    rvec *bufv;             /* Communication buffer */
#    real *bufr;             /* Communication buffer */
#    int  buf_nalloc;        /* The communication buffer size */
#
#    /* work data for solve_pme */
#    int      maxkz;
#    real *   work_mhz;
#    real *   work_m2;
#    real *   work_denom;
#    real *   work_tmp1;
#    real *   work_m2inv;

gmx_pme_t = POINTER(gmx_pme)
#typedef struct gmx_domdec_comm *gmx_domdec_comm_p_t;
#/* Abstract type for PME that is defined only in the routine that use them. */
#typedef struct gmx_pme *gmx_pme_t;

t_excl = c_ulong
#typedef unsigned long t_excl;

class t_grid(Structure):
    _fields_ = [("nr", c_int),
                ("ndim", c_int),
                ("ncg_ideal", c_int),
                ("n", ivec),
                ("ncells", c_int),
                ("cells_nalloc", c_int),
                ("ncpddc", ivec),
                ("cell_size", rvec),
                ("cell_offset", rvec),
                ("cell_index", POINTER(c_int)),
                ("index", POINTER(c_int)),
                ("nra", POINTER(c_int)),
                ("a", POINTER(c_int)),
                ("nr_alloc", c_int),
                ("dcx2", POINTER(c_real)),
                ("dcy2", POINTER(c_real)),
                ("dcz2", POINTER(c_real)),
                ("dc_nalloc", POINTER(c_int))]
#  int     nr;        /* Total number of charge groups    */
#  int    ndim;          /* The dimensionality of the grid       */
#  int    ncg_ideal;     /* The ideal number of cg's per cell    */
#  ivec     n;            /* The dimension of the grid        */
#  int    ncells;    /* Total number of cells        */
#  int    cells_nalloc;    /* Allocation size of index and nra        */
#  ivec   ncpddc;        /* The number of cells per DD cell      */
#  rvec   cell_size;     /* The size of the cells                */
#  rvec   cell_offset;   /* The offset of the cell (0,0,0)       */
#  int     *cell_index;    /* The cell number of each cg        */
#  int    *index;    /* The index into a for each cell    */
#            /* The location of the cell in the index*/
#            /* array can be found by calling xyz2ci    */
#  int    *nra;        /* The number of entries in a cell    */
#  int    *a;        /* The grid of cgs            */
#  int    nr_alloc;      /* Allocation size of cell_index and a  */
#  real   *dcx2;         /* Squared distance from atom to j-cell */
#  real   *dcy2;         /* Squared distance from atom to j-cell */
#  real   *dcz2;         /* Squared distance from atom to j-cell */
#  int    dc_nalloc;     /* Allocation size of dcx2, dyc2, dcz2  */

MAX_CG=1024
class t_ns_buf(Structure):
    _fields_ = [("ncg", c_int),
                ("nj", c_int),
                ("jcg", atom_id*MAX_CG)]
##define MAX_CG 1024
#  int     ncg;
#  int     nj;
#  atom_id jcg[MAX_CG];

class gmx_ns_t(Structure):
    _fields_ = [("simple_aaj",POINTER(atom_id)),
                ("grid", POINTER(t_grid)),
                ("bexcl", POINTER(t_excl)),
                ("bHaveVdW", POINTER(c_int)),
                ("ns_buf", POINTER(POINTER(t_ns_buf))),
                ("bExcludeAlleg", POINTER(c_int)),
                ("nra_alloc", c_int),
                ("cg_alloc", c_int),
                ("nl_sr", POINTER(POINTER(atom_id))),
                ("nsr", POINTER(c_int)),
                ("nl_lr_ljc", POINTER(POINTER(atom_id))),
                ("nl_lr_one", POINTER(POINTER(atom_id))),
                ("nlr_ljc", POINTER(c_int)),
                ("nlr_one", POINTER(c_int))]
#  atom_id  *simple_aaj;
#  t_grid   *grid;
#  t_excl   *bexcl;
#  bool     *bHaveVdW;
#  t_ns_buf **ns_buf;
#  bool     *bExcludeAlleg;
#  int      nra_alloc;
#  int      cg_alloc;
#  atom_id  **nl_sr;
#  int      *nsr;
#  atom_id  **nl_lr_ljc;
#  atom_id  **nl_lr_one;
#  int      *nlr_ljc;
#  int      *nlr_one;

class t_MMrec(Structure):
    _fields_ = [("nrMMatoms", c_int),
                ("xMM", POINTER(rvec)),
                ("indexMM", POINTER(c_int)),
                ("MMcharges", POINTER(c_real)),
                ("shiftMM", POINTER(c_int)),
                ("MMatomtype", POINTER(c_int)),
                ("scalefactor", c_real),
                ("c6", POINTER(c_real)),
                ("c12", POINTER(c_real))]
#  int           nrMMatoms;      /* nr of MM atoms, updated every step*/
#  rvec          *xMM;           /* shifted to center of box          */
#  int           *indexMM;       /* atom i = atom indexMM[I] in mdrun */
#  real          *MMcharges;     /* MM point charges in std QMMM calc.*/
#  int           *shiftMM;
#  int           *MMatomtype;    /* only important for semi-emp.      */
#  real          scalefactor;
#  /* gaussian specific stuff */
#  real          *c6;
#  real          *c12;

class t_QMrec(Structure):
    _fields_ = [
          ("nrQMatoms", c_int),
          ("xQM", POINTER(rvec)),
          ("indexQM", POINTER(c_int)),
          ("atomicnumberQM", POINTER(c_int)),
          ("QMcharges", POINTER(c_real)),
          ("shiftQM", POINTER(c_int)),
          ("QMcharge", c_int),
          ("multiplicity", c_int),
          ("QMmethod", c_int),
          ("QMbasis", c_int),
          ("nelectrons", c_int),
          ("bTS", c_int),
          ("bOPT", c_int),
          ("frontatoms", POINTER(c_int)),
          ("nQMcpus", c_int),
          ("QMmem", c_int),
          ("accuracy", c_int),
          ("cpmcscf", c_int),
          ("gauss_dir", c_char_p),
          ("gauss_exe", c_char_p),
          ("devel_dir", c_char_p),
          ("c6", POINTER(c_real)),
          ("c12", POINTER(c_real)),
          ("bSH", c_int),
          ("SAon", c_real),
          ("SAoff", c_real),
          ("SAsteps", c_int),
          ("SAstep", c_int),
          ("CIdim", c_int),
          ("CIvec1", POINTER(c_real)),
          ("CIvec2", POINTER(c_real)),
          ("CIvec1old", POINTER(c_real)),
          ("CIvec2old", POINTER(c_real)),
          ("SHbasis", ivec),
          ("CASelectrons", c_int),
          ("CASorbitals", c_int)
  ]
#  int           nrQMatoms;      /* total nr of QM atoms              */
#  rvec          *xQM;           /* shifted to center of box          */
#  int           *indexQM;       /* atom i = atom indexQM[i] in mdrun */
#  int           *atomicnumberQM;/* atomic numbers of QM atoms        */
#  real          *QMcharges;     /* atomic charges of QM atoms(ONIOM) */
#  int           *shiftQM;
#  int           QMcharge;       /* charge of the QM system           */
#  int           multiplicity;   /* multipicity (no of unpaired eln)  */
#  int           QMmethod;       /* see enums.h for all methods       */
#  int           QMbasis;        /* see enums.h for all bases         */
#  int           nelectrons;     /* total number of elecs in QM region*/
#  bool          bTS;            /* Optimize a TS, only steep, no md  */
#  bool          bOPT;          /* Optimize QM subsys, only steep, no md  */
#  bool          *frontatoms;   /* qm atoms on the QM side of a QM-MM bond */
#  /* Gaussian specific stuff */
#  int           nQMcpus;        /* no. of CPUs used for the QM calc. */
#  int           QMmem;          /* memory for the gaussian calc.     */
#  int           accuracy;       /* convergence criterium (E(-x))     */
#  bool          cpmcscf;        /* using cpmcscf(l1003)*/
#  char          *gauss_dir;
#  char          *gauss_exe;
#  char          *devel_dir;
#  real          *c6;
#  real          *c12;
#  /* Surface hopping stuff */
#  bool          bSH;            /* surface hopping (diabatic only)   */
#  real          SAon;           /* at which energy gap the SA starts */
#  real          SAoff;          /* at which energy gap the SA stops  */
#  int           SAsteps;        /* stepwise switchinng on the SA     */
#  int           SAstep;         /* current state of SA               */
#  int           CIdim;
#  real          *CIvec1;
#  real          *CIvec2;
#  real          *CIvec1old;
#  real          *CIvec2old;
#  ivec          SHbasis;
#  int           CASelectrons;
#  int           CASorbitals;

class t_QMMMrec(Structure):
    _fields_ = [("QMMMscheme", c_int),
                ("nrQMlayers", c_int),
                ("**qm", POINTER(POINTER(t_QMrec))),
                ("*mm", POINTER(t_MMrec))]
#  int           QMMMscheme; /* ONIOM (multi-layer) or normal          */
#  int           nrQMlayers; /* number of QM layers (total layers +1 (MM)) */
#  t_QMrec       **qm;        /* atoms and run params for each QM group */
#  t_MMrec       *mm;        /* there can only be one MM subsystem !   */

class t_forcerec(Structure):
    _fields_ = [("bDomDec", c_int),                ("ePBC", c_int),                ("bMolPBC", c_int),                ("rc_scaling", c_int),                ("posres_com", rvec),                ("posres_comB", rvec),                ("rlist",  c_real),
                ("rlistlong",  c_real),                ("zsquare", c_real),
                ("temp", c_real),                ("epsilon_r", c_real),
                ("epsilon_rf", c_real),
                ("epsfac", c_real),                ("kappa", c_real),
                ("k_rf", c_real),
                ("c_rf", c_real),                ("qsum", c_double*2),                ("enershiftsix", c_real),                ("enershifttwelve", c_real),                ("enerdiffsix", c_real),                ("enerdifftwelve", c_real),                ("virdiffsix", c_real),                ("virdifftwelve", c_real),                ("avcsix", c_real*2),                ("avctwelve", c_real*2),                ("fudgeQQ", c_real),                ("bcoultab", c_int),                ("bvdwtab", c_int),                ("tab14", t_forcetable),                ("rcoulomb_switch", c_real),
                ("rcoulomb", c_real),                ("phi", POINTER(c_real)),                ("rvdw_switch", c_real),
                ("rvdw", c_real),                ("bham_b_max", c_real),                ("efep", c_int),                ("sc_alpha", c_real),                ("sc_power", c_int),                ("sc_sigma6", c_real),                ("bSepDVDL", c_int),                ("eeltype", c_int),                ("vdwtype", c_int),                ("cg0", c_int),
                ("hcg", c_int),                ("solvent_opt", c_int),                ("nWatMol", c_int),                ("bGrid", c_int),                ("cginfo_global", POINTER(c_int)),                ("cginfo", POINTER(c_int)),                ("cg_cm", POINTER(rvec)),                ("cg_nalloc", c_int),                ("shift_vec", POINTER(rvec)),                ("nnblists", c_int),                ("gid2nblists", POINTER(c_int)),                ("nblists", POINTER(t_nblists)),                ("nwall", c_int),                ("wall_tab", POINTER(POINTER(t_forcetable))),                ("bTwinRange", c_int),                ("nlr", c_int),                ("f_twin_n", c_int),                ("f_twin_nalloc", c_int),                ("f_twin", POINTER(rvec)),                ("fshift_twin", POINTER(rvec)),                ("bF_NoVirSum", c_int),                ("f_novirsum_n", c_int),                ("f_novirsum_nalloc", c_int),                ("f_novirsum", POINTER(rvec)),#                ("pmedata", gmx_pme_t),
                ("pmedata", c_void_p),                ("vir_el_recip", tensor),                ("bEwald", c_int),                ("ewaldcoeff", c_real),                ("fshift", POINTER(rvec)),                ("vir_diag_posres", rvec),                ("vir_wall_z", dvec),                ("ntype", c_int),                ("bBHAM", c_int),                ("nbfp", POINTER(c_real)),                ("egp_flags", POINTER(c_int)),                ("fc_stepsize", c_real),                ("atype_radius", POINTER(c_real)),                ("atype_vol", POINTER(c_real)),                ("atype_surftens", POINTER(c_real)),                ("n_tpi", c_int),                ("ns", gmx_ns_t),                ("bQMMM", c_int),                ("qr", POINTER(t_QMMMrec)),                ("QMMMlist_sr", t_nblist),                ("QMMMlist_lr", t_nblist),                ("print_force", c_real),                ("userint1", c_int),                ("userint2", c_int),                ("userint3", c_int),                ("userint4", c_int),                ("userreal1", c_real),                ("userreal2", c_real),                ("userreal3", c_real),                ("userreal4", c_real)]#  /* Domain Decomposition */
#  bool bDomDec;
#
#  /* PBC stuff */
#  int  ePBC;
#  bool bMolPBC;
#  int  rc_scaling;
#  rvec posres_com;
#  rvec posres_comB;
#
#  /* Cut-Off stuff */
#  real rlist,rlistlong;
#
#  /* Dielectric constant resp. multiplication factor for charges */
#  real zsquare,temp;
#  real epsilon_r,epsilon_rf,epsfac;
#
#  /* Constants for reaction fields */
#  real kappa,k_rf,c_rf;
#
#  /* Charge sum for topology A/B ([0]/[1]) for Ewald corrections */
#  double qsum[2];
#
#  /* The shift of the shift or user potentials */
#  real enershiftsix;
#  real enershifttwelve;
#  /* Integrated differces for energy and virial with cut-off functions */
#  real enerdiffsix;
#  real enerdifftwelve;
#  real virdiffsix;
#  real virdifftwelve;
#  /* Constant for long range dispersion correction (average dispersion)
#   * for topology A/B ([0]/[1]) */
#  real avcsix[2];
#  /* Constant for long range repulsion term. Relative difference of about
#   * 0.1 percent with 0.8 nm cutoffs. But hey, it's cheap anyway...
#   */
#  real avctwelve[2];
#
#  /* Fudge factors */
#  real fudgeQQ;
#
#  /* Table stuff */
#  bool bcoultab;
#  bool bvdwtab;
#  /* The normal tables are in the nblists struct(s) below */
#  t_forcetable tab14; /* for 1-4 interactions only */
#
#  /* PPPM & Shifting stuff */
#  real rcoulomb_switch,rcoulomb;
#  real *phi;
#
#  /* VdW stuff */
#  real rvdw_switch,rvdw;
#  real bham_b_max;
#
#  /* Free energy ? */
#  int  efep;
#  real sc_alpha;
#  int  sc_power;
#  real sc_sigma6;
#  bool bSepDVDL;
#
#  /* NS Stuff */
#  int  eeltype;
#  int  vdwtype;
#  int  cg0,hcg;
#  /* solvent_opt contains the enum for the most common solvent
#   * in the system, which will be optimized.
#   * It can be set to esolNO to disable all water optimization */
#  int  solvent_opt;
#  int  nWatMol;
#  bool bGrid;
#  int  *cginfo_global;
#  int  *cginfo;
#  rvec *cg_cm;
#  int  cg_nalloc;
#  rvec *shift_vec;
#
#  /* The neighborlists including tables */
#  int  nnblists;
#  int  *gid2nblists;
#  t_nblists *nblists;
#
#  /* The wall tables (if used) */
#  int  nwall;
#  t_forcetable **wall_tab;
#
#  /* This mask array of length nn determines whether or not this bit of the
#   * neighbourlists should be computed. Usually all these are true of course,
#   * but not when shells are used. During minimisation all the forces that
#   * include shells are done, then after minimsation is converged the remaining
#   * forces are computed.
#   */
#  /* bool *bMask; */
#
#  /* Twin Range stuff. */
#  bool bTwinRange;
#  int  nlr;
#  int  f_twin_n;
#  int  f_twin_nalloc;
#  rvec *f_twin;
#  rvec *fshift_twin;
#
#  /* Forces that should not enter into the virial summation:
#   * PPPM/PME/Ewald/posres
#   */
#  bool bF_NoVirSum;
#  int  f_novirsum_n;
#  int  f_novirsum_nalloc;
#  rvec *f_novirsum;
#
#  /* Long-range forces and virial for PPPM/PME/Ewald */
#  gmx_pme_t pmedata;
#  tensor    vir_el_recip;
#
#  /* PME/Ewald stuff */
#  bool bEwald;
#  real ewaldcoeff;
#
#  /* Virial Stuff */
#  rvec *fshift;
#  rvec vir_diag_posres;
#  dvec vir_wall_z;
#
#  /* Non bonded Parameter lists */
#  int  ntype; /* Number of atom types */
#  bool bBHAM;
#  real *nbfp;
#
#  /* Energy group pair flags */
#  int *egp_flags;
#
#  /* xmdrun flexible constraints */
#  real fc_stepsize;
#
#  /* Generalized born stuff */
#  /* VdW radius for each atomtype (dim is thus ntype) */
#  real *atype_radius;
#  /* Effective radius (derived from effective volume) for each type */
#  real *atype_vol;
#  /* Implicit solvent - surface tension for each atomtype */
#  real *atype_surftens;
#
#  /* If > 0 signals Test Particle Insertion,
#   * the value is the number of atoms of the molecule to insert
#   * Only the energy difference due to the addition of the last molecule
#   * should be calculated.
#   */
#  bool n_tpi;
#
#  /* Neighbor searching stuff */
#  gmx_ns_t ns;
#
#  /* QMMM stuff */
#  bool         bQMMM;
#  t_QMMMrec    *qr;
#
#  /* QM-MM neighborlists */
#  t_nblist QMMMlist_sr;
#  t_nblist QMMMlist_lr; /* not needed, one QMMM list suffices */
#
#  /* Limit for printing large forces, negative is don't print */
#  real print_force;
#
#  /* User determined parameters, copied from the inputrec */
#  int  userint1;
#  int  userint2;
#  int  userint3;
#  int  userint4;
#  real userreal1;
#  real userreal2;
#  real userreal3;
#  real userreal4;

rvec5 = c_real*5
#typedef real rvec5[5];

class t_disresdata(Structure):
    _fields_ = [("dr_weighting", c_int),
                ("dr_bMixed", c_int),
                ("dr_fc", c_real),
                ("dr_tau", c_real),
                ("ETerm", c_real),
                ("ETerm1", c_real),
                ("exp_min_t_tau", c_real),
                ("nres", c_int),
                ("npair", c_int),
                ("sumviol",c_real),
                ("rt", POINTER(c_real)),
                ("rm3tav", POINTER(c_real)),
                ("Rtl_6", POINTER(c_real)),
                ("Rt_6", POINTER(c_real)),
                ("Rtav_6", POINTER(c_real))]
#/* Distance restraining stuff */
#  int  dr_weighting;  /* Weighting of pairs in one restraint              */
#  bool dr_bMixed;     /* Use sqrt of the instantaneous times              *
#               * the time averaged violation                      */
#  real dr_fc;          /* Force constant for disres,                       *
#               * which is multiplied by a (possibly)              *
#               * different factor for each restraint              */
#  real dr_tau;          /* Time constant for disres                  */
#  real ETerm;         /* multiplication factor for time averaging         */
#  real ETerm1;        /* 1 - ETerm1                                       */
#  real exp_min_t_tau; /* Factor for slowly switching on the force         */
#  int  nres;          /* The number of distance restraints                */
#  int  npair;         /* The number of distance restraint pairs           */
#  real sumviol;       /* The sum of violations                            */
#  real *rt;           /* The calculated instantaneous distance (npr)      */
#  real *rm3tav;       /* The calculated time averaged distance (npr)      */
#  real *Rtl_6;        /* The calculated instantaneous r^-6 (nr)           */
#  real *Rt_6;         /* The calculated inst. ens. averaged r^-6 (nr)     */
#  real *Rtav_6;       /* The calculated time and ens. averaged r^-6 (nr)  */

class t_oriresdata(Structure):
    _fields_ = [("fc", c_real),
                ("edt", c_real),
                ("edt1", c_real),
                ("exp_min_t_tau", c_real),
                ("nr", c_int),
                ("nex", c_int),
                ("nref", c_int),
                ("mref", POINTER(c_real)),
                ("xref", POINTER(rvec)),
                ("xtmp", POINTER(rvec)),
                ("R", matrix),
                ("S", POINTER(tensor)),
                ("Dinsl", POINTER(rvec5)),
                ("Dins", POINTER(rvec5)),
                ("Dtav", POINTER(rvec5)),
                ("oinsl", POINTER(c_real)),
                ("oins", POINTER(c_real)),
                ("otav", POINTER(c_real)),
                ("rmsdev", POINTER(c_real)),
                ("*tmp", POINTER(rvec5)),
                ("TMP", POINTER(POINTER(POINTER(c_real)))),
                ("eig", POINTER(c_real))]
#/* Orientation restraining stuff */
#  real   fc;          /* Force constant for the restraints                  */
#  real   edt;         /* Multiplication factor for time averaging           */
#  real   edt1;        /* 1 - edt                                            */
#  real   exp_min_t_tau; /* Factor for slowly switching on the force         */
#  int    nr;          /* The number of orientation restraints               */
#  int    nex;         /* The number of experiments                          */
#  int  nref;          /* The number of atoms for the fit                    */
#  real *mref;         /* The masses of the reference atoms                  */
#  rvec *xref;         /* The reference coordinates for the fit (nref)       */
#  rvec *xtmp;         /* Temporary array for fitting (nref)                 */
#  matrix R;           /* Rotation matrix to rotate to the reference coor.   */
#  tensor *S;          /* Array of order tensors for each experiment (nexp)  */
#  rvec5  *Dinsl;      /* The order matrix D for all restraints (nr x 5)     */
#  rvec5  *Dins;       /* The ensemble averaged D (nr x 5)                   */
#  rvec5  *Dtav;       /* The time and ensemble averaged D (nr x 5)          */
#  real   *oinsl;      /* The calculated instantaneous orientations          */
#  real   *oins;       /* The calculated emsemble averaged orientations      */
#  real   *otav;       /* The calculated time and ensemble averaged orient.  */
#  real   rmsdev;      /* The weighted (using kfac) RMS deviation            */
#  rvec5  *tmp;        /* An array of temporary 5-vectors (nex);             */
#  real   ***TMP;      /* An array of temporary 5x5 matrices (nex);          */
#  real   *eig;        /* Eigenvalues/vectors, for output only (nex x 12)    */

class bondedtable_t(Structure):
    _fields_ = [("n", c_int),
                ("scale", c_real),
                ("tab", POINTER(c_real))]
#  int  n;         /* n+1 is the number of points */
#  real scale;     /* distance between two points */
#  real *tab;      /* the actual tables, per point there are  4 numbers */

class t_fcdata(Structure):
    _fields_ = [("bondtab", POINTER(bondedtable_t)),
                ("angletab", POINTER(bondedtable_t)),
                ("dihtab", POINTER(bondedtable_t)),
                ("disres", t_disresdata),
                ("orires",  t_oriresdata),
                ("dihre_fc", c_real)]
#/*
# * Data struct used in the force calculation routines
# * for storing the tables for bonded interactions and
# * for storing information which is needed in following steps
# * (for instance for time averaging in distance retraints)
# * or for storing output, since force routines only return the potential.
# */
#  bondedtable_t *bondtab;
#  bondedtable_t *angletab;
#  bondedtable_t *dihtab;
#
#  t_disresdata disres;
#  t_oriresdata orires;
#  real         dihre_fc;

##define DD_MAXCELL  8
##define DD_MAXICELL 4
DD_MAXCELL  = 8
DD_MAXICELL = 4

class gmx_pme_comm_n_box(Structure):
    _fields_ = [("natoms", c_int),
                ("box", matrix),
                ("maxshift", c_int),
                ("Lambda", c_real),
                ("flags", c_int)]
#from pme_pp.c:
#  int    natoms;
#  matrix box;
#  int    maxshift;
#  real   lambda;
#  int    flags;

gmx_pme_comm_n_box_p_t = POINTER(gmx_pme_comm_n_box)
#typedef struct gmx_pme_comm_n_box *gmx_pme_comm_n_box_p_t;

class gmx_cgsort_t(Structure):
    _fields_ = [("nsc", c_int),
                ("ind_gl", c_int),
                ("ind", c_int)]
#from domdec.c:
#    int  nsc;
#    int  ind_gl;
#    int  ind;

class gmx_domdec_sort_t(Structure):
    _fields_ = [("sort1", POINTER(gmx_cgsort_t)),
                ("sort2", POINTER(gmx_cgsort_t)),
                ("sort_nalloc", c_int),
                ("sort_new", POINTER(gmx_cgsort_t)),
                ("sort_new_nalloc", c_int),
                ("vbuf", POINTER(rvec)),
                ("vbuf_nalloc", c_int),
                ("ibuf", POINTER(c_int)),
                ("ibuf_nalloc", c_int)]
#from domdec.c:
#    gmx_cgsort_t *sort1,*sort2;
#    int  sort_nalloc;
#    gmx_cgsort_t *sort_new;
#    int  sort_new_nalloc;
#    rvec *vbuf;
#    int  vbuf_nalloc;
#    int  *ibuf;
#    int  ibuf_nalloc;

class gmx_ddpme_t(Structure):
    _fields_ = [("dim", c_int),
                ("nslab", c_int),
                ("pp_min", POINTER(c_int)),
                ("pp_max", POINTER(c_int)),
                ("maxshift",c_int)]
#from domdec.c:
#    int dim;      /* The dimension                                          */
#    int nslab;    /* The number of PME slabs in this dimension              */
#    int *pp_min;  /* The minimum pp node location, size nslab               */
#    int *pp_max;  /* The maximum pp node location,size nslab                */
#    int maxshift; /* The maximum shift for coordinate redistribution in PME */

(ddnatHOME, ddnatZONE, ddnatVSITE, ddnatCON, ddnatNR) = map(c_int,xrange(5))
#/* This enum determines the order of the coordinates.
# * ddnatHOME and ddnatZONE should be first and second,
# * the others can be ordered as wanted.
# */
#enum { ddnatHOME, ddnatZONE, ddnatVSITE, ddnatCON, ddnatNR };

(ddCyclStep, ddCyclPPduringPME, ddCyclF, ddCyclPME, ddCyclNr) = map(c_int,xrange(5))
#enum { ddCyclStep, ddCyclPPduringPME, ddCyclF, ddCyclPME, ddCyclNr };

class gmx_ddzone_t(Structure):
    _fields_ =[("min0", c_real),
               ("max1", c_real),
               ("mch0", c_real),
               ("mch1", c_real),
               ("p1_0", c_real),
               ("p1_1", c_real)]
#from domdec.c:
#    real min0;    /* The minimum bottom of this zone                        */
#    real max1;    /* The maximum top of this zone                           */
#    real mch0;    /* The maximum bottom communicaton height for this zone   */
#    real mch1;    /* The maximum top communicaton height for this zone      */
#    real p1_0;    /* The bottom value of the first cell in this zone        */
#    real p1_1;    /* The top value of the first cell in this zone           */

class gmx_domdec_ind_t(Structure):
    _fields_ =[("nsend", c_int*(DD_MAXICELL+2)),
               ("nrecv", c_int*(DD_MAXICELL+2)),
               ("index", POINTER(c_int)),
               ("nalloc", c_int),
               ("cell2at0", c_int*DD_MAXICELL),
               ("cell2at1", c_int*DD_MAXICELL)]
#from domdec.c:
#    /* The numbers of charge groups to send and receive for each cell
#     * that requires communication, the last entry contains the total
#     * number of atoms that needs to be communicated.
#     */
#    int nsend[DD_MAXICELL+2];
#    int nrecv[DD_MAXICELL+2];
#    /* The charge groups to send */
#    int *index;
#    int nalloc;
#    /* The atom range for non-in-place communication */
#    int cell2at0[DD_MAXICELL];
#    int cell2at1[DD_MAXICELL];

class gmx_domdec_comm_dim_t(Structure):
    _fields_ = [("np", c_int),
                ("np_dlb", c_int),
                ("ind", POINTER(gmx_domdec_ind_t)),
                ("np_nalloc", c_int),
                ("bInPlace", c_int)]
#from domdec.c:
#    int  np;                   /* Number of grid pulses in this dimension */
#    int  np_dlb;               /* For dlb, for use with edlbAUTO          */
#    gmx_domdec_ind_t *ind;     /* The indices to communicate, size np     */
#    int  np_nalloc;
#    bool bInPlace;             /* Can we communicate in place?            */

class gmx_domdec_root_t(Structure):
    _fields_ = [("cell_size", POINTER(c_real)),
                ("bCellMin", POINTER(c_int)),
                ("cell_f", POINTER(c_real)),
                ("old_cell_f", POINTER(c_real)),
                ("cell_f_max0", POINTER(c_real)),
                ("cell_f_min1", POINTER(c_real)),
                ("bound_min", POINTER(c_real)),
                ("bound_max", POINTER(c_real)),
                ("bLimited", c_int)]
#from domdec.c:
#    real *cell_size;
#    bool *bCellMin;
#    real *cell_f;
#    real *old_cell_f;
#    real *cell_f_max0;
#    real *cell_f_min1;
#    real *bound_min;
#    real *bound_max;
#    bool bLimited;

class gmx_domdec_load_t(Structure):
    _fields_ = [("nload", c_int),
                ("load", POINTER(c_float)),
                ("sum", c_float),
                ("max", c_float),
                ("sum_m", c_float),
                ("cvol_min", c_float),
                ("mdf", c_float),
                ("pme", c_float),
                ("flags",c_int)]
#from domdec.c:
#    int  nload;
#    float *load;
#    float sum;
#    float max;
#    float sum_m;
#    float cvol_min;
#    float mdf;
#    float pme;
#    int   flags;

class gmx_domdec_comm(Structure):
    _fields_ = [("npmenodes", c_int),
                ("bCartesianPP_PME", c_int),
                ("ntot", ivec),
                ("cartpmedim", c_int),
                ("pmenodes", POINTER(c_int)),
                ("ddindex2simnodeid", POINTER(c_int)),
                ("ddpme", gmx_ddpme_t*2),
                ("bCartesianPP", c_int),
                ("ddindex2ddnodeid", POINTER(c_int)),
                ("cgs_gl", t_block),
                ("nstSortCG", c_int),
                ("sort", POINTER(gmx_domdec_sort_t)),
                ("bFilled_nsgrid_home", c_int),
                ("bInterCGBondeds", c_int),
                ("bInterCGMultiBody", c_int),
                ("bBondComm", c_int),
                ("cglink", POINTER(t_blocka)),
                ("bLocalCG", c_char_p),
                ("eDLB", c_int),
                ("bDynLoadBal", c_int),
                ("slb_frac", POINTER(POINTER(c_real))),
                ("pme_dim_f", POINTER(c_real)),
                ("cutoff_mbody", c_real),
                ("cutoff", c_real),
                ("cellsize_min", rvec),
                ("cellsize_min_dlb", rvec),
                ("cellsize_limit", c_real),
                ("v", rvec*3*3),
                ("normal", rvec*3),
                ("old_cell_x0", rvec),
                ("old_cell_x1", rvec),
                ("zone_d1", gmx_ddzone_t*2),
                ("zone_d2", gmx_ddzone_t*2*2),
                ("cd", gmx_domdec_comm_dim_t*3),
                ("maxpulse", c_int),
                ("master_cg_ddp_count", c_int),
                ("cell_ncg1", c_int*DD_MAXCELL),
                ("nat", c_int*ddnatNR.value),
                ("buf_int", POINTER(c_int)),
                ("nalloc_int", c_int),
                ("buf_int2", POINTER(c_int)),
                ("nalloc_int2", c_int),
                ("buf_vr2", POINTER(rvec)),
                ("nalloc_vr2", c_int),
                ("cggl_flag",  POINTER(POINTER(c_int))),
                ("cggl_flag_nalloc",  c_int*6),
                ("cgcm_state", POINTER(POINTER(rvec))),
                ("cgcm_state_nalloc", c_int*6),
                ("buf_vr", POINTER(rvec)),
                ("nalloc_vr", c_int),
                ("root", POINTER(POINTER(gmx_domdec_root_t))),
                ("cell_f_row",  POINTER(c_real)),
                ("cell_f0", c_real*3),
                ("cell_f1", c_real*3),
                ("cell_f_max0", c_real*3),
                ("cell_f_min1", c_real*3),
                ("bRecordLoad", c_int),
                ("load", POINTER(gmx_domdec_load_t)),
            ###########################################
            # LEAVING THE FOLLOWING OUT OF THE CLASS! #
            ###########################################
            ##ifdef GMX_MPI
            #    ("mpi_comm_load", POINTER(MPI_Comm)),
            ##endif
                ("cycl", c_float*ddCyclNr.value),
                ("cycl_n", c_int*ddCyclNr.value),
                ("eFlop", c_int),
                ("flop", c_double),
                ("flop_n", c_int),
                ("n_load_have", c_int),
                ("n_load_collect", c_int),
                ("sum_nat", c_double*(ddnatNR.value-ddnatZONE.value)),
                ("ndecomp", c_int),
                ("nload", c_int),
                ("load_step", c_double),
                ("load_sum", c_double),
                ("load_max", c_double),
                ("load_lim", ivec),
                ("load_mdf", c_double),
                ("load_pme", c_double)]
#from domdec.c
#    /* All arrays are indexed with 0 to dd->ndim (not Cartesian indexing),
#     * unless stated otherwise.
#     */
#
#    /* The number of nodes doing PME (PP/PME or only PME) */
#    int  npmenodes;
#    /* The communication setup including the PME only nodes */
#    bool bCartesianPP_PME;
#    ivec ntot;
#    int  cartpmedim;
#    int  *pmenodes;          /* size npmenodes                         */
#    int  *ddindex2simnodeid; /* size npmenodes, only with bCartesianPP
#                              * but with bCartesianPP_PME              */
#    gmx_ddpme_t ddpme[2];
#
#    /* The DD particle-particle nodes only */
#    bool bCartesianPP;
#    int  *ddindex2ddnodeid; /* size npmenode, only with bCartesianPP_PME */
#
#    /* The global charge groups */
#    t_block cgs_gl;
#
#    /* Should we sort the cgs */
#    int  nstSortCG;
#    gmx_domdec_sort_t *sort;
#    bool bFilled_nsgrid_home;
#
#    /* Are there bonded and multi-body interactions between charge groups? */
#    bool bInterCGBondeds;
#    bool bInterCGMultiBody;
#
#    /* Data for the optional bonded interaction atom communication range */
#    bool bBondComm;
#    t_blocka *cglink;
#    char *bLocalCG;
#
#    /* The DLB option */
#    int  eDLB;
#    /* Are we actually using DLB? */
#    bool bDynLoadBal;
#
#    /* Cell sizes for static load balancing, first index cartesian */
#    real **slb_frac;
#    /* Cell sizes for determining the PME communication with SLB */
#    real *pme_dim_f;
#
#    /* The width of the communicated boundaries */
#    real cutoff_mbody;
#    real cutoff;
#    /* The minimum cell size (including triclinic correction) */
#    rvec cellsize_min;
#    /* For dlb, for use with edlbAUTO */
#    rvec cellsize_min_dlb;
#    /* The lower limit for the DD cell size with DLB */
#    real cellsize_limit;
#
#    /* Orthogonal vectors for triclinic cells, Cartesian index */
#    rvec v[DIM][DIM];
#    /* Normal vectors for the cells walls */
#    rvec normal[DIM];
#
#    /* The old location of the cell boundaries, to check cg displacements */
#    rvec old_cell_x0;
#    rvec old_cell_x1;
#
#    /* The zone limits for DD dimensions 1 and 2 (not 0), determined from
#     * cell boundaries of neighboring cells for dynamic load balancing.
#     */
#    gmx_ddzone_t zone_d1[2];
#    gmx_ddzone_t zone_d2[2][2];
#
#    /* The coordinate/force communication setup and indices */
#    gmx_domdec_comm_dim_t cd[DIM];
#    /* The maximum number of cells to communicate with in one dimension */
#    int  maxpulse;
#
#    /* Which cg distribution is stored on the master node */
#    int master_cg_ddp_count;
#
#    /* The number of cg's received from the direct neighbors */
#    int  cell_ncg1[DD_MAXCELL];
#
#    /* The atom counts, the range for each type t is nat[t-1] <= at < nat[t] */
#    int  nat[ddnatNR];
#
#    /* Communication buffer for general use */
#    int  *buf_int;
#    int  nalloc_int;
#
#    /* Communication buffers only used with multiple grid pulses */
#    int  *buf_int2;
#    int  nalloc_int2;
#    rvec *buf_vr2;
#    int  nalloc_vr2;
#
#    /* Communication buffers for local redistribution */
#    int  **cggl_flag;
#    int  cggl_flag_nalloc[DIM*2];
#    rvec **cgcm_state;
#    int  cgcm_state_nalloc[DIM*2];
#    rvec *buf_vr;
#    int  nalloc_vr;
#
#    /* Cell sizes for dynamic load balancing */
#    gmx_domdec_root_t **root;
#    real *cell_f_row;
#    real cell_f0[DIM];
#    real cell_f1[DIM];
#    real cell_f_max0[DIM];
#    real cell_f_min1[DIM];
#
#    /* Stuff for load communication */
#    bool bRecordLoad;
#    gmx_domdec_load_t *load;
##ifdef GMX_MPI
#    MPI_Comm *mpi_comm_load;
##endif
#    /* Cycle counters */
#    float cycl[ddCyclNr];
#    int   cycl_n[ddCyclNr];
#    /* Flop counter (0=no,1=yes,2=with (eFlop-1)*5% noise */
#    int eFlop;
#    double flop;
#    int    flop_n;
#    /* Have often have did we have load measurements */
#    int    n_load_have;
#    /* Have often have we collected the load measurements */
#    int    n_load_collect;
#
#    /* Statistics */
#    double sum_nat[ddnatNR-ddnatZONE];
#    int    ndecomp;
#    int    nload;
#    double load_step;
#    double load_sum;
#    double load_max;
#    ivec   load_lim;
#    double load_mdf;
#    double load_pme;

gmx_domdec_comm_p_t = POINTER(gmx_domdec_comm)
#typedef struct gmx_domdec_comm *gmx_domdec_comm_p_t;

class gmx_domdec_master(Structure):
    _fields_ = [("cell_x", POINTER(POINTER(c_real))),
                ("ncg", POINTER(c_int)),
                ("index", POINTER(c_int)),
                ("cg", POINTER(c_int)),
                ("nat", POINTER(c_int)),
                ("ibuf", POINTER(c_int)),
                ("vbuf", POINTER(rvec))]
#from domdec.c:
#    /* The cell boundaries */
#    real **cell_x;
#    /* The global charge group division */
#    int  *ncg;     /* Number of home charge groups for each node */
#    int  *index;   /* Index of nnodes+1 into cg */
#    int  *cg;      /* Global charge group index */
#    int  *nat;     /* Number of home atoms for each node. */
#    int  *ibuf;    /* Buffer for communication */
#    rvec *vbuf;    /* Buffer for state scattering and gathering */

gmx_domdec_master_p_t = POINTER(gmx_domdec_master)
#typedef struct gmx_domdec_master *gmx_domdec_master_p_t;

class gmx_reverse_ilist_t(Structure):
    _fields_ = [("index", POINTER(c_int)),
                ("il", POINTER(c_int))]
#from domdec_top.c:
#    int  *index;  /* Index for each atom into il                  */
#    int  *il;     /* ftype|type|a0|...|an|ftype|...               */

class gmx_molblock_ind_t(Structure):
    _fields_ = [("a_start", c_int),
                ("a_end", c_int),
                ("natoms_mol", c_int),
                ("type", c_int)]
#from domdec_top.c:
#    int  a_start;
#    int  a_end;
#    int  natoms_mol;
#    int  type;

class gmx_reverse_top(Structure):
    _fields_ = [("bExclRequired", c_int),
                ("bConstr", c_int),
                ("bBCheck", c_int),
                ("bMultiCGmols", c_int),
                ("ril_mt", POINTER(gmx_reverse_ilist_t)),
                ("ril_mt_tot_size", c_int),
                ("ilsort", c_int),
                ("mbi", POINTER(gmx_molblock_ind_t))]
#from domdec_top.c:
#    bool bExclRequired; /* Do we require all exclusions to be assigned? */
#    bool bConstr;       /* Are there constraints in this revserse top?  */
#    bool bBCheck;       /* All bonded interactions have to be assigned? */
#    bool bMultiCGmols;  /* Are the multi charge-group molecules?        */
#    gmx_reverse_ilist_t *ril_mt; /* Reverse ilist for all moltypes      */
#    int  ril_mt_tot_size;
#    int  ilsort;        /* The sorting state of bondeds for free energy */
#    gmx_molblock_ind_t *mbi;

gmx_reverse_top_p_t = POINTER(gmx_reverse_top)
#typedef struct gmx_reverse_top *gmx_reverse_top_p_t;

class gmx_specatsend_t(Structure):
    _fields_ = [("nsend", c_int),
                ("a", POINTER(c_int)),
                ("a_nalloc", c_int),
                ("nrecv", c_int)]
#from domdec_con:
#    int nsend;
#    int *a;
#    int a_nalloc;
#    int nrecv;

class gmx_domdec_specat_comm(Structure):
    _fields_ = [("nind_req", c_int),
                ("ind_req", POINTER(c_int)),
                ("ind_req_nalloc", c_int),
                ("nreq", c_int*3*2*2),
                ("spas", gmx_specatsend_t*3*2),
                ("bSendAtom", POINTER(c_int)),
                ("bSendAtom_nalloc", c_int),
                ("ibuf", POINTER(c_int)),
                ("ibuf_nalloc", c_int),
                ("vbuf", POINTER(rvec)),
                ("vbuf_nalloc", c_int),
                ("vbuf2", POINTER(rvec)),
                ("vbuf2_nalloc", c_int),
                ("at_start", c_int),
                ("at_end", c_int)]
#from domdec_con.c:
#    /* The atom indices we need from the surrounding cells */
#    int  nind_req;
#    int  *ind_req;
#    int  ind_req_nalloc;
#    /* The number of indices to receive during the setup */
#    int  nreq[DIM][2][2];
#    /* The atoms to send */
#    gmx_specatsend_t spas[DIM][2];
#    bool *bSendAtom;
#    int   bSendAtom_nalloc;
#    /* Send buffers */
#    int  *ibuf;
#    int  ibuf_nalloc;
#    rvec *vbuf;
#    int  vbuf_nalloc;
#    rvec *vbuf2;
#    int  vbuf2_nalloc;
#    /* The range in the local buffer(s) for received atoms */
#    int  at_start;
#    int  at_end;

gmx_domdec_specat_comm_p_t = POINTER(gmx_domdec_specat_comm)
#typedef struct gmx_domdec_specat_comm *gmx_domdec_specat_comm_p_t;

class gmx_domdec_constraints(Structure):
    _fields_ = [("molb_con_offset", POINTER(c_int)),
                ("molb_ncon_mol", POINTER(c_int)),
                ("ncon", c_int),
                ("con_gl", POINTER(c_int)),
                ("con_nlocat", POINTER(c_int)),
                ("con_nalloc", c_int),
                ("gc_req", c_char_p),
                ("ga2la",POINTER(c_int))]
#from domdec_con:
#    int  *molb_con_offset;
#    int  *molb_ncon_mol;
#    /* The fully local and connected constraints */
#    int  ncon;
#    /* The global constraint number, only required for clearing gc_req */
#    int  *con_gl;
#    int  *con_nlocat;
#    int  con_nalloc;
#    /* Boolean that tells if a global constraint index has been requested */
#    char *gc_req;
#    /* Global to local communicated constraint atom only index */
#    int  *ga2la;

gmx_domdec_constraints_p_t = POINTER(gmx_domdec_constraints)
#typedef struct gmx_domdec_constraints *gmx_domdec_constraints_p_t;

class gmx_ga2la_t(Structure):
    _fields_ = [("cell",c_int),
                ("a",   atom_id)]
#  int     cell;
#  atom_id a;

class gmx_domdec_ns_ranges_t(Structure):
    _fields_ = [("j0", c_int),
                ("j1", c_int),
                ("cg1", c_int),
                ("jcg0", c_int),
                ("jcg1", c_int),
                ("shift0", ivec),
                ("shift1", ivec)]
#  int  j0;       /* j-cell start               */
#  int  j1;       /* j-cell end                 */
#  int  cg1;      /* i-charge-group end         */
#  int  jcg0;     /* j-charge-group start       */
#  int  jcg1;     /* j-charge-group end         */
#  ivec shift0;   /* Minimum shifts to consider */
#  ivec shift1;   /* Maximum shifts to consider */

class gmx_partdec(Structure):
    _fields_ = [("neighbor", c_int*2),
                ("cgindex", POINTER(c_int)),
                ("index", POINTER(c_int)),
                ("shift", c_int),
                ("bshift", c_int)]
#from partdec.c:
#  int  neighbor[2];             /* The nodeids of left and right neighb */
#  int  *cgindex;                /* The charge group boundaries,         */
#                /* size nnodes+1,                       */
#                                /* only allocated with particle decomp. */
#  int  *index;                  /* The home particle boundaries,        */
#                /* size nnodes+1,                       */
#                                /* only allocated with particle decomp. */
#  int  shift,bshift;        /* Coordinates are shifted left for     */
#                                /* 'shift' systolic pulses, and right   */
#                /* for 'bshift' pulses. Forces are      */
#                /* shifted right for 'shift' pulses     */
#                /* and left for 'bshift' pulses         */
#                /* This way is not necessary to shift   */
#                /* the coordinates over the entire ring */

gmx_partdec_p_t = POINTER(gmx_partdec)
#typedef struct gmx_partdec *gmx_partdec_p_t;

class gmx_multisim_t(Structure):
    _fields_ = [("nsim", c_int),
                ("sim", c_int)#,
              ###########################################
              # LEAVING THE FOLLOWING OUT OF THE CLASS! #
              ###########################################
              ##ifdef GMX_MPI
              #  ("mpi_group_masters", MPI_Group),
              #  ("mpi_comm_masters",   MPI_Comm)
              ##endif
               ]
#  int nsim;
#  int sim;
##ifdef GMX_MPI
#  MPI_Group mpi_group_masters;
#  MPI_Comm mpi_comm_masters;
##endif

class gmx_domdec_t(Structure):
    _fields_ = [("nnodes", c_int),
              ###########################################
              # LEAVING THE FOLLOWING OUT OF THE CLASS! #
              ###########################################
              ##ifdef GMX_MPI
              #  ("mpi_comm_all", MPI_Comm),
              ##endif
                ("bSendRecv2", c_int),
                ("ci", ivec),
                ("rank", c_int),
                ("master_ci", ivec),
                ("masterrank", c_int),
                ("pme_nodeid", c_int),
                ("pme_receive_vir_ener", c_int),
                ("cnb", gmx_pme_comm_n_box_p_t),
              ###########################################
              # LEAVING THE FOLLOWING OUT OF THE CLASS! #
              ###########################################
              ##ifdef GMX_MPI
              #  ("nreq_pme", c_int),
              #  ("req_pme", MPI_Request*4),
              ##endif
                ("nc", ivec),
                ("ndim", c_int),
                ("dim", ivec),
                ("bGridJump", c_int),
                ("ncell", c_int),
                ("shift", ivec*DD_MAXCELL),
                ("bScrewPBC", c_int),
                ("tric_dir", ivec),
                ("skew_fac", rvec),
                ("cell_x0", rvec),
                ("cell_x1", rvec),
                ("cell_ns_x0", rvec),
                ("cell_ns_x1", rvec),
                ("neighbor", c_int*3*2),
                ("ma", gmx_domdec_master_p_t),
                ("bInterCGcons", c_int),
                ("reverse_top", gmx_reverse_top_p_t),
                ("nbonded_global", c_int),
                ("nbonded_local", c_int),
                ("n_intercg_excl", c_int),
                ("ga2la_vsite", POINTER(c_int)),
                ("vsite_comm", gmx_domdec_specat_comm_p_t),
                ("constraints", gmx_domdec_constraints_p_t),
                ("constraint_comm", gmx_domdec_specat_comm_p_t),
                ("ncg_cell", c_int*(DD_MAXCELL+1)),
                ("ncg_home", c_int),
                ("ncg_tot", c_int),
                ("index_gl", POINTER(c_int)),
                ("cgindex", POINTER(c_int)),
                ("cg_nalloc", c_int),
                ("la2lc", POINTER(c_int)),
                ("la2lc_nalloc", c_int),
                ("nat_home", c_int),
                ("nat_tot", c_int),
                ("gatindex", POINTER(c_int)),
                ("gatindex_nalloc", c_int),
                ("ga2la", POINTER(gmx_ga2la_t)),
                ("nicell", c_int),
                ("icell", gmx_domdec_ns_ranges_t*DD_MAXICELL),
                ("comm", gmx_domdec_comm_p_t),
                ("ddp_count", c_int)]
#  /* The DD particle-particle nodes only */
#  /* The communication setup within the communicator all
#   * defined in dd->comm in domdec.c
#   */
#  int  nnodes;
##ifdef GMX_MPI
#  MPI_Comm mpi_comm_all;
##endif
#  /* Use MPI_Sendrecv communication instead of non-blocking calls */
#  bool bSendRecv2;
#  /* The local DD cell index and rank */
#  ivec ci;
#  int  rank;
#  ivec master_ci;
#  int  masterrank;
#  /* Communication with the PME only nodes */
#  int  pme_nodeid;
#  bool pme_receive_vir_ener;
#  gmx_pme_comm_n_box_p_t cnb;
##ifdef GMX_MPI
#  int  nreq_pme;
#  MPI_Request req_pme[4];
##endif
#
#
#  /* The communication setup, identical for each cell, cartesian index */
#  ivec nc;
#  int  ndim;
#  ivec dim;  /* indexed by 0 to ndim */
#  bool bGridJump;
#  /* The bonded and non-bonded communication setup, cartesian index */
#  int  ncell;
#  ivec shift[DD_MAXCELL];
#
#  /* Screw PBC? */
#  bool bScrewPBC;
#
#  /* Tells if the box is skewed for each of the three cartesian directions */
#  ivec tric_dir;
#  rvec skew_fac;
#
#  /* The cell boundaries */
#  rvec cell_x0;
#  rvec cell_x1;
#  /* The extreme sizes of the local cells needed for the neighbor searching */
#  rvec cell_ns_x0;
#  rvec cell_ns_x1;
#  /* Forward and backward neighboring cells, indexed by 0 to ndim */
#  int  neighbor[DIM][2];
#
#  /* Only available on the master node */
#  gmx_domdec_master_p_t ma;
#
#  /* Are there inter charge group constraints */
#  bool bInterCGcons;
#
#  /* Global atom number to interaction list */
#  gmx_reverse_top_p_t reverse_top;
#  int  nbonded_global;
#  int  nbonded_local;
#
#  /* The number of inter charge-group exclusions */
#  int  n_intercg_excl;
#
#  /* Vsite stuff */
#  int  *ga2la_vsite;
#  gmx_domdec_specat_comm_p_t vsite_comm;
#
#  /* Constraint stuff */
#  gmx_domdec_constraints_p_t constraints;
#  gmx_domdec_specat_comm_p_t constraint_comm;
#
#  /* The charge group boundaries for the cells */
#  int ncg_cell[DD_MAXCELL+1];
#
#  /* The local to gobal charge group index and local cg to local atom index */
#  int  ncg_home;
#  int  ncg_tot;
#  int  *index_gl;
#  int  *cgindex;
#  int  cg_nalloc;
#  /* Local atom to local cg index, only for special cases */
#  int  *la2lc;
#  int  la2lc_nalloc;
#
#  /* The number of home atoms */
#  int  nat_home;
#  /* The total number of atoms: home and received zones */
#  int  nat_tot;
#  /* Index from the local atoms to the global atoms */
#  int  *gatindex;
#  int  gatindex_nalloc;
#
#  /* Global atom number to local atom number, -1 if not local */
#  gmx_ga2la_t *ga2la;
#
#  /* For neighborsearching */
#  int  nicell;
#  gmx_domdec_ns_ranges_t icell[DD_MAXICELL];
#
#  /* Communication stuff */
#  gmx_domdec_comm_p_t comm;
#
#  /* The partioning count, to keep track of the state */
#  int ddp_count;

class gmx_nodecomm_t(Structure):
    _fields_ =[("bUse", c_int)#,
           ###########################################
           # LEAVING THE FOLLOWING OUT OF THE CLASS! #
           ###########################################
           ##ifdef GMX_MPI
           #   ("comm_intra", MPI_Comm),
           #   ("rank_intra", c_int),
           #   ("comm_inter", MPI_Comm)
           ##endif
               ]
#  int      bUse;
##ifdef GMX_MPI
#  MPI_Comm comm_intra;
#  int      rank_intra;
#  MPI_Comm comm_inter;
##endif

class t_commrec(Structure):
    _fields_ = [("sim_nodeid", c_int),
                ("nnodes", c_int),
                ("npmenodes", c_int),
                ("threadid", c_int),
                ("nthreads", c_int),
                ("nodeid", c_int),
                ###########################################
                # LEAVING THE FOLLOWING OUT OF THE CLASS! #
                ###########################################
                ##ifdef GMX_MPI
                #("mpi_comm_mysim", MPI_Comm),
                #("mpi_comm_mygroup", MPI_Comm),
                ##endif
                ("nc", gmx_nodecomm_t),
                ("dd", POINTER(gmx_domdec_t)),
                ("pd", gmx_partdec_p_t),
                ("duty", c_int),
                ("ms",  POINTER(gmx_multisim_t))]
#  /* The nodids in one sim are numbered sequentially from 0.
#   * All communication within some simulation should happen
#   * in mpi_comm_mysim, or its subset mpi_comm_mygroup.
#   */
#  int sim_nodeid,nnodes,npmenodes;
#  int threadid,nthreads;
#  /* The nodeid in the PP/PME, PP or PME group */
#  int nodeid;
##ifdef GMX_MPI
#  MPI_Comm mpi_comm_mysim;
#  MPI_Comm mpi_comm_mygroup;
##endif
#
#  gmx_nodecomm_t nc;
#
#  /* For domain decomposition */
#  gmx_domdec_t *dd;
#
#  /* For particle decomposition */
#  gmx_partdec_p_t pd;
#
#  /* The duties of this node, see the defines above */
#  int duty;
#
#  gmx_multisim_t *ms;

(
  efMDP,
  efGCT,
  efTRX,
  efTRO,
  efTRN,
  efTRR,
  efTRJ,
  efXTC,
  efG87,
  efENX,
  efEDR,
  efENE,
  efSTX,
  efSTO,
  efGRO,
  efG96,
  efPDB,
  efBRK,
  efENT,
  efESP,
  efPQR,
  efCPT,
  efLOG,
  efXVG,
  efOUT,
  efNDX,
  efTOP,
  efITP,
  efTPX,
  efTPS,
  efTPR,
  efTPA,
  efTPB,
  efTEX,
  efRTP,
  efATP,
  efHDB,
  efDAT,
  efDLG,
  efMAP,
  efEPS,
  efMAT,
  efM2P,
  efMTX,
  efEDI,
  efEDO,
  efHAT,
  efXPM,
  efNR
) = map(c_int,xrange(49))
#/* this enum should correspond to the array deffile in gmxlib/filenm.c */
#enum {
#  efMDP, efGCT,
#  efTRX, efTRO, efTRN, efTRR, efTRJ, efXTC, efG87,
#  efENX, efEDR, efENE,
#  efSTX, efSTO, efGRO, efG96, efPDB, efBRK, efENT, efESP, efPQR,
#  efCPT,
#  efLOG, efXVG, efOUT,
#  efNDX,
#  efTOP, efITP,
#  efTPX, efTPS, efTPR, efTPA, efTPB,
#  efTEX, efRTP, efATP, efHDB,
#  efDAT, efDLG,
#  efMAP, efEPS, efMAT, efM2P,
#  efMTX,
#  efEDI, efEDO,
#  efHAT,
#  efXPM,
#  efNR
#};

ffSET       = (1<<0)
ffREAD      = (1<<1)
ffWRITE     = (1<<2)
ffOPT       = (1<<3)
ffLIB       = (1<<4)
ffMULT      = (1<<5)
ffRW        = (ffREAD | ffWRITE)
ffOPTRD     = (ffREAD | ffOPT)
ffOPTWR     = (ffWRITE| ffOPT)
ffOPTRW     = (ffRW   | ffOPT)
ffLIBRD     = (ffREAD | ffLIB)
ffLIBOPTRD  = (ffOPTRD | ffLIB)
ffRDMULT    = (ffREAD  | ffMULT)
ffOPTRDMULT = (ffRDMULT | ffOPT)
ffWRMULT    = (ffWRITE  | ffMULT)
ffOPTWRMULT = (ffWRMULT | ffOPT)
ffSET       = c_ulong(ffSET)
ffREAD      = c_ulong(ffREAD)
ffWRITE     = c_ulong(ffWRITE)
ffOPT       = c_ulong(ffOPT)
ffLIB       = c_ulong(ffLIB)
ffMULT      = c_ulong(ffMULT)
ffRW        = c_ulong(ffRW)
ffOPTRD     = c_ulong(ffOPTRD)
ffOPTWR     = c_ulong(ffOPTWR)
ffOPTRW     = c_ulong(ffOPTRW)
ffLIBRD     = c_ulong(ffLIBRD)
ffLIBOPTRD  = c_ulong(ffLIBOPTRD)
ffRDMULT    = c_ulong(ffRDMULT)
ffOPTRDMULT = c_ulong(ffOPTRDMULT)
ffWRMULT    = c_ulong(ffWRMULT)
ffOPTWRMULT = c_ulong(ffOPTWRMULT)
##define ffSET   1<<0
##define ffREAD  1<<1
##define ffWRITE 1<<2
##define ffOPT   1<<3
##define ffLIB   1<<4
##define ffMULT  1<<5
##define ffRW    (ffREAD | ffWRITE)
##define ffOPTRD (ffREAD | ffOPT)
##define ffOPTWR (ffWRITE| ffOPT)
##define ffOPTRW (ffRW   | ffOPT)
##define ffLIBRD (ffREAD | ffLIB)
##define ffLIBOPTRD (ffOPTRD | ffLIB)
##define ffRDMULT   (ffREAD  | ffMULT)
##define ffOPTRDMULT   (ffRDMULT | ffOPT)
##define ffWRMULT   (ffWRITE  | ffMULT)
##define ffOPTWRMULT   (ffWRMULT | ffOPT)

class t_filenm(Structure):
    _fields_ = [("ftp", c_int),
               ("opt", c_char_p),
               ("fn", c_char_p),
               ("flag", c_ulong),
               ("nfiles", c_int),
               ("fns", POINTER(c_char_p))]
#  int  ftp;             /* File type (see enum above)           */
#  char *opt;            /* Command line option                  */
#  char *fn;             /* File name (as set in source code)    */
#  unsigned long flag;   /* Flag for all kinds of info (see defs)*/
#  int  nfiles;          /* number of files                      */
#  char **fns;           /* File names                           */

(
  etINT,
  etREAL,
  etTIME,
  etSTR,
  etBOOL,
  etRVEC,
  etENUM,
  etNR
) = map(c_int,xrange(8))
#from ./include/gromacs/readinp.h:
#/* This structure is used for parsing arguments off the comand line */
#enum {
#  etINT, etREAL, etTIME, etSTR,    etBOOL, etRVEC,   etENUM, etNR
#};


class u_u(Union):
    _fields_ = [("v", c_void_p),
                ("i", POINTER(c_int)),
                ("r", POINTER(c_real)),
                ("c", POINTER(c_char_p)),
                ("b",POINTER(c_int)),
                ("rv",POINTER(rvec))]

class t_pargs(Structure):
    _fields_ = [("option", c_char_p),
                ("bSet", c_int),
                ("type", c_int),
                ("u", u_u),
                ("desc", c_char_p)]
#from ./include/gromacs/readinp.h:
#  char *option;
#  bool bSet;
#  int  type;
#  union {
#    void *v;   /* This is a nasty workaround, to be able to use initialized */
#    int  *i;   /* arrays */
#    real *r;
#    char **c;  /* Must be pointer to string (when type == etSTR)         */
#               /* or null terminated list of enums (when type == etENUM) */
#    bool *b;
#    rvec *rv;
#  } u;
#  char *desc;

#################################
# Don't know if this is correct #
#################################
FILE_p = c_void_p

class gmx_edx(Structure):
    _fields_ = [("nr", c_int),
                ("nr_loc", c_int),
                ("anrs", POINTER(c_int)),
                ("anrs_loc", POINTER(c_int)),
                ("c_ind", POINTER(c_int)),
                ("x", POINTER(rvec)),
                ("x_old", POINTER(rvec)),
                ("m", POINTER(c_real)),
                ("mtot", c_real),
                ("sqrtm", POINTER(c_real))]
#from ./src/mdlib/edsam.c:
#/* This type is for the average, reference, target, and origin structure    */
#typedef struct gmx_edx
#{
#    int           nr;             /* number of atoms this structure contains  */
#    int           nr_loc;         /* number of atoms on local node            */
#    int           *anrs;          /* atom index numbers                       */
#    int           *anrs_loc;      /* local atom index numbers                 */
#    int           *c_ind;         /* at which position of the whole anrs
#                                   * array is a local atom?, i.e.
#                                   * c_ind[0...nr_loc-1] gives the atom index
#                                   * with respect to the collective
#                                   * anrs[0...nr-1] array                     */
#    rvec          *x;             /* coordinates for this structure           */
#    rvec          *x_old;         /* used to keep track of the shift vectors
#                                     such that the ED molecule can always be
#                                     made whole in the parallel case          */
#    real          *m;             /* masses                                   */
#    real          mtot;           /* total mass (only used in sref)           */
#    real          *sqrtm;         /* sqrt of the masses used for mass-
#                                   * weighting of analysis (only used in sav) */
#} t_gmx_edx;

class t_eigvec(Structure):
    _fields_ = [("neig", c_int),
                ("ieig", POINTER(c_int)),
                ("stpsz", POINTER(c_real)),
                ("vec", POINTER(POINTER(rvec))),
                ("xproj", POINTER(c_real)),
                ("fproj", POINTER(c_real)),
                ("refproj", POINTER(c_real)),
                ("radius", c_real)]
#from ./src/mdlib/edsam.c:
#    int    neig;     /* nr of eigenvectors             */
#    int   *ieig;     /* index nrs of eigenvectors      */
#    real  *stpsz;    /* stepsizes (per eigenvector)    */
#    rvec  **vec;     /* eigenvector components         */
#    real  *xproj;    /* instantaneous x projections    */
#    real  *fproj;    /* instantaneous f projections    */
#    real  *refproj;  /* starting or target projecions  */
#    real  radius;    /* instantaneous radius           */

class t_edvecs(Structure):
    _fields_ = [("mon", t_eigvec),
                ("linfix", t_eigvec),
                ("linacc", t_eigvec),
                ("radfix", t_eigvec),
                ("radacc", t_eigvec),
                ("radcon", t_eigvec)]
#from ./src/mdlib/edsam.c:
#typedef struct
#{
#    t_eigvec      mon;            /* only monitored, no constraints       */
#    t_eigvec      linfix;         /* fixed linear constraints             */
#    t_eigvec      linacc;         /* acceptance linear constraints        */
#    t_eigvec      radfix;         /* fixed radial constraints (exp)       */
#    t_eigvec      radacc;         /* acceptance radial constraints (exp)  */
#    t_eigvec      radcon;         /* acceptance rad. contraction constr.  */
#} t_edvecs;

class t_edflood(Structure):
    _fields_ = [("deltaF0", c_real),
                ("bHarmonic", c_int),
                ("tau", c_real),
                ("deltaF", c_real),
                ("Efl", c_real),
                ("kT", c_real),
                ("Vfl", c_real),
                ("dt", c_real),
                ("constEfl", c_real),
                ("alpha2", c_real),
                ("flood_id", c_int),
                ("forces_cartesian", POINTER(rvec)),
                ("vecs",t_eigvec)]
#from ./src/mdlib/edsam.c:
#typedef struct
#{
#    real deltaF0;
#    bool bHarmonic;
#    real tau;
#    real deltaF;
#    real Efl;
#    real kT;
#    real Vfl;
#    real dt;
#    real constEfl;
#    real alpha2;
#    int flood_id;
#    rvec *forces_cartesian;
#    t_eigvec vecs;         /* use flooding for these */
#} t_edflood;

class t_do_edfit(Structure):
    _fields_ = [("omega", POINTER(POINTER(c_double))),
                ("om", POINTER(POINTER(c_double)))]
#from ./src/mdlib/edsam.c:
#struct t_do_edfit {
#    double **omega;
#    double **om;
#};

class t_do_edsam(Structure):
    _fields_ = [("old_rotmat", matrix),
                ("oldrad", c_real),
                ("old_transvec", rvec),
                ("older_transvec", rvec),
                ("transvec_compact", rvec),
                ("xcoll", POINTER(rvec)),
                ("xc_ref", POINTER(rvec)),
                ("shifts_xcoll", POINTER(ivec)),
                ("extra_shifts_xcoll", POINTER(ivec)),
                ("shifts_xc_ref", POINTER(ivec)),
                ("extra_shifts_xc_ref", POINTER(ivec)),
                ("bUpdateShifts", c_int)]
#from ./src/mdlib/edsam.c:
#struct t_do_edsam
#{
#    matrix old_rotmat;
#    real oldrad;
#    rvec old_transvec,older_transvec,transvec_compact;
#    rvec *xcoll;         /* Coordinates from all nodes, this is the collective set of coords we work on.
#                          * These are the coordinates of atoms with average structure indices */
#    rvec *xc_ref;        /* same but with reference structure indices */
#    ivec *shifts_xcoll;        /* Shifts for xcoll  */
#    ivec *extra_shifts_xcoll;  /* xcoll shift changes since last NS step */
#    ivec *shifts_xc_ref;       /* Shifts for xc_ref */
#    ivec *extra_shifts_xc_ref; /* xc_ref shift changes since last NS step */
#    bool bUpdateShifts;        /* TRUE in NS steps to indicate that the ED shifts
#                                * for this ED dataset need to be updated */
#};

class t_do_radcon(Structure):
    _fields_  = [("proj", POINTER(c_real))]
#from ./src/mdlib/edsam.c:
#struct t_do_radcon {
#    real *proj;
#};

###############################################################################################
# t_fitit and t_remove_pbc_effect not, or ill defined in GMX, making empty classes of them... #
###############################################################################################
class t_fitit(Structure):
    pass

class t_remove_pbc_effect(Structure):
    pass

class t_ed_buffer(Structure):
    _fields_ = [("fitit", POINTER(t_fitit)),
                ("do_edfit", POINTER(t_do_edfit)),
                ("remove_pbc_effect", POINTER(t_remove_pbc_effect)),
                ("do_edsam", POINTER(t_do_edsam)),
                ("do_radcon", POINTER(t_do_radcon))]
#from ./src/mdlib/edsam.c:
#/* definition of ED buffer structure */
#struct t_ed_buffer
#{
#    struct t_fitit *                fitit;
#    struct t_do_edfit *             do_edfit;
#    struct t_remove_pbc_effect *    remove_pbc_effect;
#    struct t_do_edsam *             do_edsam;
#    struct t_do_radcon *            do_radcon;
#};

########################################################
# WORKAROUND TO AVOID THE POINTER TO SELF SITUATION... #
########################################################
class edpar(Structure):
    pass

class t_edpar(Structure):
    _fields_ = [("nini", c_int),
                ("fitmas", c_int),
                ("pcamas", c_int),
                ("presteps", c_int),
                ("outfrq", c_int),
                ("maxedsteps", c_int),
                ("sref", gmx_edx),
                ("bRefEqAv", c_int),
                ("sav", gmx_edx),
                ("star", gmx_edx),
                ("sori", gmx_edx),
                ("vecs", t_edvecs),
                ("slope", c_real),
                ("bNeedDoEdsam", c_int),
                ("flood", t_edflood),
                ("buf", POINTER(t_ed_buffer)),
                ("next_edi", POINTER(edpar))]
#from ./src/mdlib/edsam.c:
#    int            nini;           /* total Nr of atoms                    */
#    bool           fitmas;         /* true if trans fit with cm            */
#    bool           pcamas;         /* true if mass-weighted PCA            */
#    int            presteps;       /* number of steps to run without any
#                                    *    perturbations ... just monitoring */
#    int            outfrq;         /* freq (in steps) of writing to edo    */
#    int            maxedsteps;     /* max nr of steps per cycle            */
#
#    /* all gmx_edx datasets are copied to all nodes in the parallel case    */
#    struct gmx_edx sref;           /* reference positions, to these fitting
#                                    * will be done                         */
#    bool           bRefEqAv;       /* If true, reference & average indices
#                                    * are the same. Used for optimization  */
#    struct gmx_edx sav;            /* average positions                    */
#    struct gmx_edx star;           /* target positions                     */
#    struct gmx_edx sori;           /* origin positions                     */
#
#    t_edvecs       vecs;           /* eigenvectors                         */
#    real           slope;          /* minimal slope in acceptance radexp   */
#
#    bool           bNeedDoEdsam;   /* if any of the options mon, linfix, ...
#                                    * is used (i.e. apart from flooding)   */
#    t_edflood      flood;          /* parameters especially for flooding   */
#    struct t_ed_buffer *buf;       /* handle to local buffers              */
#    struct edpar   *next_edi;      /* Pointer to another ed dataset        */

#class t_edpar(Structure):
#    _fields_ = [("nini", c_int),
#                ("fitmas", c_int),
#                ("pcamas", c_int),
#                ("presteps", c_int),
#                ("outfrq", c_int),
#                ("maxedsteps", c_int),
#                ("sref", gmx_edx),
#                ("bRefEqAv", c_int),
#                ("sav", gmx_edx),
#                ("star", gmx_edx),
#                ("sori", gmx_edx),
#                ("vecs", t_edvecs),
#                ("slope", c_real),
#                ("bNeedDoEdsam", c_int),
#                ("flood", t_edflood),
#                ("buf", POINTER(t_ed_buffer)),
#                ("next_edi", POINTER(edpar))]
##from ./src/mdlib/edsam.c:
##    int            nini;           /* total Nr of atoms                    */
##    bool           fitmas;         /* true if trans fit with cm            */
##    bool           pcamas;         /* true if mass-weighted PCA            */
##    int            presteps;       /* number of steps to run without any
##                                    *    perturbations ... just monitoring */
##    int            outfrq;         /* freq (in steps) of writing to edo    */
##    int            maxedsteps;     /* max nr of steps per cycle            */
##
##    /* all gmx_edx datasets are copied to all nodes in the parallel case    */
##    struct gmx_edx sref;           /* reference positions, to these fitting
##                                    * will be done                         */
##    bool           bRefEqAv;       /* If true, reference & average indices
##                                    * are the same. Used for optimization  */
##    struct gmx_edx sav;            /* average positions                    */
##    struct gmx_edx star;           /* target positions                     */
##    struct gmx_edx sori;           /* origin positions                     */
##
##    t_edvecs       vecs;           /* eigenvectors                         */
##    real           slope;          /* minimal slope in acceptance radexp   */
##
##    bool           bNeedDoEdsam;   /* if any of the options mon, linfix, ...
##                                    * is used (i.e. apart from flooding)   */
##    t_edflood      flood;          /* parameters especially for flooding   */
##    struct t_ed_buffer *buf;       /* handle to local buffers              */
##    struct edpar   *next_edi;      /* Pointer to another ed dataset        */

class gmx_edsam(Structure):
    _fields_ = [("eEDtype", c_int),
                ("edinam", c_char_p),
                ("edonam", c_char_p),
                ("edo",    FILE_p),
                ("edpar",  POINTER(t_edpar))]
#from ./src/mdlib/edsam.c:
#    int           eEDtype;        /* Type of ED: see enums above          */
#    char          *edinam;        /* name of ED sampling input file       */
#    char          *edonam;        /*                     output           */
#    FILE          *edo;           /* output file pointer                  */
#    t_edpar       *edpar;

gmx_edsam_t = POINTER(gmx_edsam)
#typedef struct gmx_edsam *gmx_edsam_t;

class t_comm_vsites(Structure):
    _fields_ = [("nprevvsite", c_int),
                ("nnextvsite", c_int),
                ("idxprevvsite", POINTER(c_int)),
                ("idxnextvsite", POINTER(c_int)),
                ("nprevconstr", c_int),
                ("nnextconstr", c_int),
                ("idxprevconstr", POINTER(c_int)),
                ("idxnextconstr", POINTER(c_int))]
#from include/vsite.h:
#typedef struct {
#  int nprevvsite; /* how many virtual sites are nonlocal */
#  int nnextvsite;
#  int *idxprevvsite; /* index of nonlocal vsite particles */
#  int *idxnextvsite;
#  int nprevconstr; /* how many constr. atoms are nonlocal */
#  int nnextconstr;
#  int *idxprevconstr; /* indices of nonlocal constructing atoms */
#  int *idxnextconstr;
#} t_comm_vsites;

class gmx_vsite_t(Structure):
    _fields_ = [("n_intercg_vsite", c_int),
                ("nvsite_pbc_molt", c_int),
                ("vsite_pbc_molt", POINTER(POINTER(POINTER(c_int)))),
                ("vsite_pbc_loc", POINTER(POINTER(c_int))),
                ("vsite_pbc_loc_nalloc", POINTER(c_int)),
                ("bPDvsitecomm", c_int),
                ("vsitecomm", POINTER(t_comm_vsites))]
#from include/vsite.h:
#typedef struct {
#  int  n_intercg_vsite;       /* The number of inter charge group vsites */
#  int  nvsite_pbc_molt;       /* The array size of vsite_pbc_molt        */
#  int  ***vsite_pbc_molt;     /* The pbc atoms for intercg vsites        */
#  int  **vsite_pbc_loc;       /* The local pbc atoms                     */
#  int  *vsite_pbc_loc_nalloc;
#  bool bPDvsitecomm;          /* Do we need vsite communication with PD? */
#  t_comm_vsites *vsitecomm;   /* The PD vsite communication struct       */
#} gmx_vsite_t;

class gmx_lincsdata(Structure):
    _fields_ = [("ncg", c_int),
                ("ncg_flex", c_int),
                ("ncg_triangle", c_int),
                ("nIter", c_int),
                ("nOrder", c_int),
                ("nc", c_int),
                ("nc_alloc", c_int),
                ("ncc", c_int),
                ("ncc_alloc", c_int),
                ("matlam", c_real),
                ("bllen0", POINTER(c_real)),
                ("ddist", POINTER(c_real)),
                ("bla", POINTER(c_int)),
                ("blc", POINTER(c_real)),
                ("blnr", POINTER(c_int)),
                ("blbnb", POINTER(c_int)),
                ("ntriangle", c_int),
                ("triangle", POINTER(c_int)),
                ("tri_bits", POINTER(c_int)),
                ("ncc_triangle", c_int),
                ("blmf", POINTER(c_real)),
                ("bllen", POINTER(c_real)),
                ("tmpv", POINTER(rvec)),
                ("tmpncc", POINTER(c_real)),
                ("tmp1", POINTER(c_real)),
                ("tmp2", POINTER(c_real)),
                ("tmp3", POINTER(c_real)),
                ("Lambda", POINTER(c_real)),
                ("rmsd_data", c_real*3)]
#from src/mdlib/clincs.c:
#typedef struct gmx_lincsdata {
#    int  ncg;         /* the global number of constraints */
#    int  ncg_flex;    /* the global number of flexible constraints */
#    int  ncg_triangle;/* the global number of constraints in triangles */
#    int  nIter;       /* the number of iterations */
#    int  nOrder;      /* the order of the matrix expansion */
#    int  nc;          /* the number of constraints */
#    int  nc_alloc;    /* the number we allocated memory for */
#    int  ncc;         /* the number of constraint connections */
#    int  ncc_alloc;   /* the number we allocated memory for */
#    real matlam;      /* the FE lambda value used for filling blc and blmf */
#    real *bllen0;     /* the reference distance in topology A */
#    real *ddist;      /* the reference distance in top B - the r.d. in top A */
#    int  *bla;        /* the atom pairs involved in the constraints */
#    real *blc;        /* 1/sqrt(invmass1 + invmass2) */
#    int  *blnr;       /* index into blbnb and blmf */
#    int  *blbnb;      /* list of constraint connections */
#    int  ntriangle;   /* the local number of constraints in triangles */
#    int  *triangle;   /* the list of triangle constraints */
#    int  *tri_bits;   /* the bits tell if the matrix element should be used */
#    int  ncc_triangle;/* the number of constraint connections in triangles */
#    real *blmf;       /* matrix of mass factors for constraint connections */
#    real *bllen;      /* the reference bond length */
#    /* arrays for temporary storage in the LINCS algorithm */
#    rvec *tmpv;
#    real *tmpncc;
#    real *tmp1;
#    real *tmp2;
#    real *tmp3;
#    real *lambda;  /* the Lagrange multipliers */
#    /* storage for the constraint RMS relative deviation output */
#    real rmsd_data[3];
#} t_gmx_lincsdata;

gmx_lincsdata_t = POINTER(gmx_lincsdata)
#from include/types/constr.h:
#/* Abstract type for LINCS that is defined only in the file that uses it */
#typedef struct gmx_lincsdata *gmx_lincsdata_t;

class gmx_constr(Structure):
    _fields_ = [("ncon_tot", c_int),
                ("nflexcon", c_int),
                ("n_at2con_mt", c_int),
                ("at2con_mt", POINTER(t_blocka)),
                ("lincsd", gmx_lincsdata_t),
                ("nblocks", c_int),
                ("sblock", POINTER(c_int)),
                ("sblock_nalloc", c_int),
                ("lagr", POINTER(c_real)),
                ("lagr_nalloc",c_int),
                ("maxwarn", c_int),
                ("warncount_lincs", c_int),
                ("warncount_settle", c_int),
                ("ed", gmx_edsam_t),
                ("warn_mtop", POINTER(gmx_mtop_t))]
#from src/mdlib/constr.c:
#typedef struct gmx_constr {
#  int             ncon_tot;     /* The total number of constraints    */
#  int             nflexcon;     /* The number of flexible constraints */
#  int             n_at2con_mt;  /* The size of at2con = #moltypes     */
#  t_blocka        *at2con_mt;   /* A list of atoms to constraints     */
#  gmx_lincsdata_t lincsd;       /* LINCS data                         */
#  int             nblocks;      /* The number of SHAKE blocks         */
#  int             *sblock;      /* The SHAKE blocks                   */
#  int             sblock_nalloc;/* The allocation size of sblock      */
#  real            *lagr;        /* Lagrange multipliers for SHAKE     */
#  int             lagr_nalloc;  /* The allocation size of lagr        */
#  int             maxwarn;      /* The maximum number of warnings     */
#  int             warncount_lincs;
#  int             warncount_settle;
#  gmx_edsam_t     ed;           /* The essential dynamics data        */
#
#  gmx_mtop_t      *warn_mtop;   /* Only used for printing warnings    */
#} t_gmx_constr;

gmx_constr_t = POINTER(gmx_constr)
#from include/types/constr.h
#/* Abstract type for constraints */
#typedef struct gmx_constr *gmx_constr_t;

class t_energy(Structure):
    _fields_= [("e", c_real),
               ("eav", c_double),
               ("esum", c_double),
               ("e2sum", c_real)]
#from energy.h:
#typedef struct {
#  real e;    /* The current energy.                    */
#  double eav;     /* The running average                       */
#  double esum;    /* The sum of energies until now.            */
#  real e2sum;    /* The sum of the square of energies until now        */
#} t_energy;

class t_ebin(Structure):
    _fields_ = [("nener", c_int),
                ("enm", POINTER(c_char_p)),
                ("e", POINTER(t_energy))]
#from ebin.h:
#typedef struct {
#  int      nener;
#  char     **enm;
#  t_energy *e;
#} t_ebin;

class t_mdebin(Structure):
    _fields_ = [("ebin", POINTER(t_ebin)),
                ("ie", c_int),
                ("iconrmsd", c_int),
                ("ib", c_int),
                ("isvir", c_int),
                ("ifvir", c_int),
                ("ipres", c_int),
                ("ivir", c_int),
                ("isurft", c_int),
                ("ipc", c_int),
                ("itemp", c_int),
                ("itc", c_int),
                ("iu", c_int),
                ("imu", c_int),
                ("ivcos", c_int),
                ("ivisc", c_int),
                ("nE", c_int),
                ("nEg", c_int),
                ("nEc", c_int),
                ("nTC", c_int),
                ("nU", c_int),
                ("igrp", POINTER(c_int))]
#from mdebin.h:
#typedef struct {
#  t_ebin *ebin;
#  int    ie,iconrmsd,ib,isvir,ifvir,ipres,ivir,isurft,ipc,itemp,itc,iu,imu;
#  int    ivcos,ivisc;
#  int    nE,nEg,nEc,nTC,nU;
#  int    *igrp;
#} t_mdebin;

# A new struct for the communication between the mdrunner parts, defined in mdrun.h
class mdrunner_comm (Structure):
    _fields_ = [("inputrec", POINTER(t_inputrec)),
                ("state", POINTER(t_state)),
                ("box", matrix),
                ("buf", POINTER(rvec)),
                ("f", POINTER(rvec)),
                ("nrnb", POINTER(t_nrnb)),
                ("mtop", POINTER(gmx_mtop_t)),
                ("mdatoms", POINTER(t_mdatoms)),
                ("fr", POINTER(t_forcerec)),
                ("fcd", POINTER(t_fcdata)),
                ("ewaldcoeff", c_real),
                ("pmedata", c_void_p), # not using its contents anyway...
            #    ("pmedata", POINTER(gmx_pme_t)),
                ("vsite", POINTER(gmx_vsite_t)),
                ("wcycle", gmx_wallcycle_t),
                ("nsteps_done",  c_int),
                ("mdebin", POINTER(t_mdebin)),
                ("xsave", POINTER(rvec)),
                ("vsave", POINTER(rvec))]
#typedef struct {
#    t_inputrec      *inputrec;
#    t_state         *state;
#    matrix          box;
#    rvec            *buf;
#    rvec            *f;
#    t_nrnb          *nrnb;
#    gmx_mtop_t      *mtop;
#    t_mdatoms       *mdatoms;
#    t_forcerec      *fr;
#    t_fcdata        *fcd;
#    real            ewaldcoeff;
#    gmx_pme_t       *pmedata;
#    gmx_vsite_t     *vsite;
#    gmx_wallcycle_t wcycle;
#    int             nsteps_done;
#    t_mdebin        *mdebin;
#} mdrunner_comm;