#from physics.h
M_PI      = 3.14159265358979323846
KILO      = 1.0e3
NANO	  = 1.0e-9
PICO	  = 1.0e-12
AMU       = 1.6605402e-27
BOLTZMANN = 1.380658e-23
AVOGADRO  = 6.0221367e23
RGAS      = BOLTZMANN*AVOGADRO
BOLTZ     = RGAS/KILO
PLANCK1   = 6.6262e-34
PLANCK    = 6.6262e-34*AVOGADRO/(PICO*KILO)

#from physics.h
#define M_PI        3.14159265358979323846
#define KILO 		 (1e3)			/* Thousand	*/
#define NANO		 (1e-9)			/* A Number	*/
#define PICO		 (1e-12)		/* A Number	*/
#define AMU              (1.6605402e-27)        /* kg           */
#define BOLTZMANN     (1.380658e-23)        /* (J/K)    */
#define AVOGADRO     (6.0221367e23)        /* ()        */
#define RGAS             (BOLTZMANN*AVOGADRO)   /* (J/(mol K))  */
#define BOLTZ            (RGAS/KILO)            /* (kJ/(mol K)) */
#define PLANCK1          (6.6262e-34)           /* J s */
#define PLANCK           (6.6262e-34*AVOGADRO/(PICO*KILO)) /* (kJ/mol) ps */