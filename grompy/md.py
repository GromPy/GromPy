import sys, signal
from ctypes import c_char_p,\
                   c_double,\
                   POINTER,\
                   c_int,\
                   addressof,\
                   byref,\
                   c_ulong,\
                   c_char,\
                   sizeof,\
                   c_size_t,\
                   c_uint,\
                   c_void_p
from grompy import libmdrun,matrix,rvec,c_real,ivec,stderr,libc,stdout
import grompy.types as gt
import grompy.commrec as commrec
import grompy.statutil as statutil
import grompy.enums as enums
import grompy.inputrec as inptrec
import grompy.vec as vec
import grompy.names as names
import grompy.sim_util as su
import grompy.gmx_wallcycle as wc
#from grompy.types import t_inputrec,\
#                         t_state,\
#                         t_nrnb,\
#                         gmx_mtop_t,\
#                         t_mdatoms,\
#                         t_forcerec,\
#                         t_fcdata,\
#                         t_commrec

# from mdrun.h
MD_GLAS         = (1<<1)
MD_IONIZE       = (1<<3)
MD_RERUN        = (1<<4)
MD_RERUN_VSITE  = (1<<5)
MD_SEPPOT       = (1<<7)
MD_PARTDEC      = (1<<9)
MD_DDBONDCHECK  = (1<<10)
MD_DDBONDCOMM   = (1<<11)
MD_CONFOUT      = (1<<12)
MD_NOGSTAT      = (1<<13)
MD_REPRODUCIBLE = (1<<14)
MD_READ_RNG     = (1<<15)
MD_APPENDFILES  = (1<<16)
MD_READ_EKIN    = (1<<17)
MD_STARTFROMCPT = (1<<18)

# from mdrun.h (octals)
LIST_SCALARS  = 0001
LIST_INPUTREC = 0002
LIST_TOP      = 0004
LIST_X        = 0010
LIST_V        = 0020
LIST_F        = 0040
LIST_LOAD     = 0100

def mdrun():
    (ddnoSEL,ddnoINTERLEAVE,ddnoPP_PME,ddnoCARTESIAN,ddnoNR) = map(c_int,xrange(5))

    descarray = c_char_p*191
    desc      = descarray(
        c_char_p("The mdrun program is the main computational chemistry engine"),
        c_char_p("within GROMACS. Obviously, it performs Molecular Dynamics simulations,"),
        c_char_p("but it can also perform Stochastic Dynamics, Energy Minimization,"),
        c_char_p("test particle insertion or (re)calculation of energies."),
        c_char_p("Normal mode analysis is another option. In this case mdrun"),
        c_char_p("builds a Hessian matrix from single conformation."),
        c_char_p("For usual Normal Modes-like calculations, make sure that"),
        c_char_p("the structure provided is properly energy-minimized."),
        c_char_p("The generated matrix can be diagonalized by g_nmeig.[PAR]"),
        c_char_p("The mdrun program reads the run input file ([TT]-s[tt])"),
        c_char_p("and distributes the topology over nodes if needed."),
        c_char_p("mdrun produces at least four output files."),
        c_char_p("A single log file ([TT]-g[tt]) is written, unless the option"),
        c_char_p("[TT]-seppot[tt] is used, in which case each node writes a log file."),
        c_char_p("The trajectory file ([TT]-o[tt]), contains coordinates, velocities and"),
        c_char_p("optionally forces."),
        c_char_p("The structure file ([TT]-c[tt]) contains the coordinates and"),
        c_char_p("velocities of the last step."),
        c_char_p("The energy file ([TT]-e[tt]) contains energies, the temperature,"),
        c_char_p("pressure, etc, a lot of these things are also printed in the log file."),
        c_char_p("Optionally coordinates can be written to a compressed trajectory file"),
        c_char_p("([TT]-x[tt]).[PAR]"),
        c_char_p("The option [TT]-dgdl[tt] is only used when free energy perturbation is"),
        c_char_p("turned on.[PAR]"),
        c_char_p("When mdrun is started using MPI with more than 1 node, parallelization"),
        c_char_p("is used. By default domain decomposition is used, unless the [TT]-pd[tt]"),
        c_char_p("option is set, which selects particle decomposition.[PAR]"),
        c_char_p("With domain decomposition, the spatial decomposition can be set"),
        c_char_p("with option [TT]-dd[tt]. By default mdrun selects a good decomposition."),
        c_char_p("The user only needs to change this when the system is very inhomogeneous."),
        c_char_p("Dynamic load balancing is set with the option [TT]-dlb[tt],"),
        c_char_p("which can give a significant performance improvement,"),
        c_char_p("especially for inhomogeneous systems. The only disadvantage of"),
        c_char_p("dynamic load balancing is that runs are no longer binary reproducible,"),
        c_char_p("but in most cases this is not important."),
        c_char_p("By default the dynamic load balancing is automatically turned on"),
        c_char_p("when the measured performance loss due to load imbalance is 5% or more."),
        c_char_p("At low parallelization these are the only important options"),
        c_char_p("for domain decomposition."),
        c_char_p("At high parallelization the options in the next two sections"),
        c_char_p("could be important for increasing the performace."),
        c_char_p("[PAR]"),
        c_char_p("When PME is used with domain decomposition, separate nodes can"),
        c_char_p("be assigned to do only the PME mesh calculation;"),
        c_char_p("this is computationally more efficient starting at about 12 nodes."),
        c_char_p("The number of PME nodes is set with option [TT]-npme[tt],"),
        c_char_p("this can not be more than half of the nodes."),
        c_char_p("By default mdrun makes a guess for the number of PME"),
        c_char_p("nodes when the number of nodes is larger than 11 or performance wise"),
        c_char_p("not compatible with the PME grid x dimension."),
        c_char_p("But the user should optimize npme. Performance statistics on this issue"),
        c_char_p("are written at the end of the log file."),
        c_char_p("For good load balancing at high parallelization, the PME grid x and y"),
        c_char_p("dimensions should be divisible by the number of PME nodes"),
        c_char_p("(the simulation will run correctly also when this is not the case)."),
        c_char_p("[PAR]"),
        c_char_p("This section lists all options that affect the domain decomposition."),
        c_char_p("[BR]"),
        c_char_p("Option [TT]-rdd[tt] can be used to set the required maximum distance"),
        c_char_p("for inter charge-group bonded interactions."),
        c_char_p("Communication for two-body bonded interactions below the non-bonded"),
        c_char_p("cut-off distance always comes for free with the non-bonded communication."),
        c_char_p("Atoms beyond the non-bonded cut-off are only communicated when they have"),
        c_char_p("missing bonded interactions; this means that the extra cost is minor"),
        c_char_p("and nearly indepedent of the value of [TT]-rdd[tt]."),
        c_char_p("With dynamic load balancing option [TT]-rdd[tt] also sets"),
        c_char_p("the lower limit for the domain decomposition cell sizes."),
        c_char_p("By default [TT]-rdd[tt] is determined by mdrun based on"),
        c_char_p("the initial coordinates. The chosen value will be a balance"),
        c_char_p("between interaction range and communication cost."),
        c_char_p("[BR]"),
        c_char_p("When inter charge-group bonded interactions are beyond"),
        c_char_p("the bonded cut-off distance, mdrun terminates with an error message."),
        c_char_p("For pair interactions and tabulated bonds"),
        c_char_p("that do not generate exclusions, this check can be turned off"),
        c_char_p("with the option [TT]-noddcheck[tt]."),
        c_char_p("[BR]"),
        c_char_p("When constraints are present, option [TT]-rcon[tt] influences"),
        c_char_p("the cell size limit as well."),
        c_char_p("Atoms connected by NC constraints, where NC is the LINCS order plus 1,"),
        c_char_p("should not be beyond the smallest cell size. A error message is"),
        c_char_p("generated when this happens and the user should change the decomposition"),
        c_char_p("or decrease the LINCS order and increase the number of LINCS iterations."),
        c_char_p("By default mdrun estimates the minimum cell size required for P-LINCS"),
        c_char_p("in a conservative fashion. For high parallelization it can be useful"),
        c_char_p("to set the distance required for P-LINCS with the option [TT]-rcon[tt]."),
        c_char_p("[BR]"),
        c_char_p("The [TT]-dds[tt] option sets the minimum allowed x, y and/or z scaling"),
        c_char_p("of the cells with dynamic load balancing. mdrun will ensure that"),
        c_char_p("the cells can scale down by at least this factor. This option is used"),
        c_char_p("for the automated spatial decomposition (when not using [TT]-dd[tt])"),
        c_char_p("as well as for determining the number of grid pulses, which in turn"),
        c_char_p("sets the minimum allowed cell size. Under certain circumstances"),
        c_char_p("the value of [TT]-dds[tt] might need to be adjusted to account for"),
        c_char_p("high or low spatial inhomogeneity of the system."),
        c_char_p("[PAR]"),
        c_char_p("The option [TT]-nosum[tt] can be used to only sum the energies"),
        c_char_p("at every neighbor search step and energy output step."),
        c_char_p("This can improve performance for highly parallel simulations"),
        c_char_p("where this global communication step becomes the bottleneck."),
        c_char_p("For a global thermostat and/or barostat the temperature"),
        c_char_p("and/or pressure will also only be updated every nstlist steps."),
        c_char_p("With this option the energy file will not contain averages and"),
        c_char_p("fluctuations over all integration steps.[PAR]"),
        c_char_p("With [TT]-rerun[tt] an input trajectory can be given for which "),
        c_char_p("forces and energies will be (re)calculated. Neighbor searching will be"),
        c_char_p("performed for every frame, unless [TT]nstlist[tt] is zero"),
        c_char_p("(see the [TT].mdp[tt] file).[PAR]"),
        c_char_p("ED (essential dynamics) sampling is switched on by using the [TT]-ei[tt]"),
        c_char_p("flag followed by an [TT].edi[tt] file."),
        c_char_p("The [TT].edi[tt] file can be produced using options in the essdyn"),
        c_char_p("menu of the WHAT IF program. mdrun produces a [TT].edo[tt] file that"),
        c_char_p("contains projections of positions, velocities and forces onto selected"),
        c_char_p("eigenvectors.[PAR]"),
        c_char_p("When user-defined potential functions have been selected in the"),
        c_char_p("[TT].mdp[tt] file the [TT]-table[tt] option is used to pass mdrun"),
        c_char_p("a formatted table with potential functions. The file is read from"),
        c_char_p("either the current directory or from the GMXLIB directory."),
        c_char_p("A number of pre-formatted tables are presented in the GMXLIB dir,"),
        c_char_p("for 6-8, 6-9, 6-10, 6-11, 6-12 Lennard Jones potentials with"),
        c_char_p("normal Coulomb."),
        c_char_p("When pair interactions are present a separate table for pair interaction"),
        c_char_p("functions is read using the [TT]-tablep[tt] option.[PAR]"),
        c_char_p("When tabulated bonded functions are present in the topology,"),
        c_char_p("interaction functions are read using the [TT]-tableb[tt] option."),
        c_char_p("For each different tabulated interaction type the table file name is"),
        c_char_p("modified in a different way: before the file extension an underscore is"),
        c_char_p("appended, then a b for bonds, an a for angles or a d for dihedrals"),
        c_char_p("and finally the table number of the interaction type.[PAR]"),
        c_char_p("The options [TT]-px[tt] and [TT]-pf[tt] are used for writing pull COM"),
        c_char_p("coordinates and forces when pulling is selected"),
        c_char_p("in the [TT].mdp[tt] file.[PAR]"),
        c_char_p("With [TT]-multi[tt] multiple systems are simulated in parallel."),
        c_char_p("As many input files are required as the number of systems."),
        c_char_p("The system number is appended to the run input and each output filename,"),
        c_char_p("for instance topol.tpr becomes topol0.tpr, topol1.tpr etc."),
        c_char_p("The number of nodes per system is the total number of nodes"),
        c_char_p("divided by the number of systems."),
        c_char_p("One use of this option is for NMR refinement: when distance"),
        c_char_p("or orientation restraints are present these can be ensemble averaged"),
        c_char_p("over all the systems.[PAR]"),
        c_char_p("With [TT]-replex[tt] replica exchange is attempted every given number"),
        c_char_p("of steps. The number of replicas is set with the [TT]-multi[tt] option,"),
        c_char_p("see above."),
        c_char_p("All run input files should use a different coupling temperature,"),
        c_char_p("the order of the files is not important. The random seed is set with"),
        c_char_p("[TT]-reseed[tt]. The velocities are scaled and neighbor searching"),
        c_char_p("is performed after every exchange.[PAR]"),
        c_char_p("Finally some experimental algorithms can be tested when the"),
        c_char_p("appropriate options have been given. Currently under"),
        c_char_p("investigation are: polarizability, glass simulations"),
        c_char_p("and X-Ray bombardments."),
        c_char_p("[PAR]"),
        c_char_p("The option [TT]-pforce[tt] is useful when you suspect a simulation"),
        c_char_p("crashes due to too large forces. With this option coordinates and"),
        c_char_p("forces of atoms with a force larger than a certain value will"),
        c_char_p("be printed to stderr."),
        c_char_p("[PAR]"),
        c_char_p("Checkpoints containing the complete state of the system are written"),
        c_char_p("at regular intervals (option [TT]-cpt[tt]) to the file [TT]-cpo[tt],"),
        c_char_p("unless option [TT]-cpt[tt] is set to -1."),
        c_char_p("A simulation can be continued by reading the full state from file"),
        c_char_p("with option [TT]-cpi[tt]. This option is intelligent in the way that"),
        c_char_p("if no checkpoint file is found, Gromacs just assumes a normal run and"),
        c_char_p("starts from the first step of the tpr file."),
        c_char_p("The simulation part number is added to all output files,"),
        c_char_p("unless [TT]-append[tt] or [TT]-noaddpart[tt] are set."),
        c_char_p("[PAR]"),
        c_char_p("With checkpointing you can also use the option [TT]-append[tt] to"),
        c_char_p("just continue writing to the previous output files. This is not"),
        c_char_p("enabled by default since it is potentially dangerous if you move files,"),
        c_char_p("but if you just leave all your files in place and restart mdrun with"),
        c_char_p("exactly the same command (with options [TT]-cpi[tt] and [TT]-append[tt])"),
        c_char_p("the result will be the same as from a single run. The contents will"),
        c_char_p("be binary identical (unless you use dynamic load balancing),"),
        c_char_p("but for technical reasons there might be some extra energy frames when"),
        c_char_p("using checkpointing (necessary for restarts without appending)."),
        c_char_p("[PAR]"),
        c_char_p("With option [TT]-maxh[tt] a simulation is terminated and a checkpoint"),
        c_char_p("file is written at the first neighbor search step where the run time"),
        c_char_p("exceeds [TT]-maxh[tt]*0.99 hours."),
        c_char_p("[PAR]"),
        c_char_p("When mdrun receives a TERM signal, it will set nsteps to the current"),
        c_char_p("step plus one. When mdrun receives a USR1 signal, it will stop after"),
        c_char_p("the next neighbor search step (with nstlist=0 at the next step)."),
        c_char_p("In both cases all the usual output will be written to file."),
        c_char_p("When running with MPI, a signal to one of the mdrun processes"),
        c_char_p("is sufficient, this signal should not be sent to mpirun or"),
        c_char_p("the mdrun process that is the parent of the others."),
        c_char_p("[PAR]"),
        c_char_p("When mdrun is started with MPI, it does not run niced by default.")
        )

    pt_commrec = POINTER(gt.t_commrec)
    cr         = pt_commrec

    fnmarray = gt.t_filenm * 27
    fnm      = fnmarray(gt.t_filenm(ftp=gt.efTPX, opt=c_char_p(),          fn=c_char_p(),           flag=gt.ffREAD , nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efTRN, opt=c_char_p("-o"),      fn=c_char_p(),           flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXTC, opt=c_char_p("-x"),      fn=c_char_p(),           flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efCPT, opt=c_char_p("-cpi"),    fn=c_char_p(),           flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efCPT, opt=c_char_p("-cpo"),    fn=c_char_p(),           flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efSTO, opt=c_char_p("-c"),      fn=c_char_p("confout"),  flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efENX, opt=c_char_p("-e"),      fn=c_char_p("ener"),     flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efLOG, opt=c_char_p("-g"),      fn=c_char_p("md"),       flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-dgdl"),   fn=c_char_p("dgdl"),     flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-field"),  fn=c_char_p("field"),    flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-table"),  fn=c_char_p("table"),    flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tablep"), fn=c_char_p("tablep"),   flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tableb"), fn=c_char_p("table"),    flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efTRX, opt=c_char_p("-rerun"),  fn=c_char_p("rerun"),    flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tpi"),    fn=c_char_p("tpi"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tpid"),   fn=c_char_p("tpidist"),  flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efEDI, opt=c_char_p("-ei"),     fn=c_char_p("sam"),      flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efEDO, opt=c_char_p("-eo"),     fn=c_char_p("sam"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efGCT, opt=c_char_p("-j"),      fn=c_char_p("wham"),     flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efGCT, opt=c_char_p("-jo"),     fn=c_char_p("bam"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-ffout"),  fn=c_char_p("gct"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-devout"), fn=c_char_p("deviatie"), flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-runav"),  fn=c_char_p("runaver"),  flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-px"),     fn=c_char_p("pullx"),    flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-pf"),     fn=c_char_p("pullf"),    flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efMTX, opt=c_char_p("-mtx"),    fn=c_char_p("nm"),       flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efNDX, opt=c_char_p("-dn"),     fn=c_char_p("dipole"),   flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()))
    NFILE =  c_int(len(fnm))

    bCart         = c_int(gt.FALSE)
    bPPPME        = c_int(gt.FALSE)
    bPartDec      = c_int(gt.FALSE)
    bDDBondCheck  = c_int(gt.TRUE)
    bDDBondComm   = c_int(gt.TRUE)
    bSumEner      = c_int(gt.TRUE)
    bVerbose      = c_int(gt.FALSE)
    bCompact      = c_int(gt.TRUE)
    bSepPot       = c_int(gt.FALSE)
    bRerunVSite   = c_int(gt.FALSE)
    bGlas         = c_int(gt.FALSE)
    bIonize       = c_int(gt.FALSE)
    bConfout      = c_int(gt.TRUE)
    bReproducible = c_int(gt.FALSE)

    npme         = c_int(-1)
    nmultisim    = c_int(0)
    repl_ex_nst  = c_int(0)
    repl_ex_seed = c_int(-1)
    nstepout     = c_int(100)
    nthreads     = c_int(1)

    realddxyz  = rvec(c_real(0.0),c_real(0.0),c_real(0.0))
    ddnoarray  = c_char_p*(ddnoNR.value+1)
    ddno_opt   = ddnoarray(c_char_p(gt.NULL),
                           c_char_p("interleave"),
                           c_char_p("pp_pme"),
                           c_char_p("cartesian"),
                           c_char_p(gt.NULL))
    ddndlbarray = c_char_p*5
    dddlb_opt   = ddndlbarray(c_char_p(gt.NULL),
                              c_char_p("auto"),
                              c_char_p("no"),
                              c_char_p("yes"),
                              c_char_p(gt.NULL))

    rdd          = c_real(0.0)
    rconstr      = c_real(0.0)
    dlb_scale    = c_real(0.8)
    pforce       = c_real(-1)
    ddcsx        = c_char_p(gt.NULL)
    ddcsy        = c_char_p(gt.NULL)
    ddcsz        = c_char_p(gt.NULL)
    cpt_period   = c_real(15.0)
    max_hours    = c_real(-1)
    bAppendFiles = c_int(gt.FALSE)
    bAddPart     = c_int(gt.TRUE)
    t_pargsarray = gt.t_pargs*32
    pa           = t_pargsarray(
                        gt.t_pargs(option=c_char_p("-pd"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bPartDec)),
                                   desc=c_char_p("Use particle decompostion")),
                        gt.t_pargs(option=c_char_p("-dd"),
                                   bSet=gt.FALSE,
                                   type=gt.etRVEC,
                                   u=gt.u_u(addressof(realddxyz)),
                                   desc=c_char_p("Domain decomposition grid, 0 is optimize")),
                        gt.t_pargs(option=c_char_p("-nt"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(nthreads)),
                                   desc=c_char_p("HIDDENNumber of threads to start on each node")),
                        gt.t_pargs(option=c_char_p("-npme"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(npme)),
                                   desc=c_char_p("Number of separate nodes to be used for PME, -1 is guess")),
                        gt.t_pargs(option=c_char_p("-ddorder"),
                                   bSet=gt.FALSE,
                                   type=gt.etENUM,
                                   u=gt.u_u(addressof(ddno_opt)),
                                   desc=c_char_p("DD node order")),
                        gt.t_pargs(option=c_char_p("-ddcheck"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bDDBondCheck)),
                                   desc=c_char_p("Check for all bonded interactions with DD")),
                        gt.t_pargs(option=c_char_p("-ddbondcomm"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bDDBondComm)),
                                   desc=c_char_p("HIDDENUse special bonded atom communication when -rdd > cut-off")),
                        gt.t_pargs(option=c_char_p("-rdd"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(rdd)),
                                   desc=c_char_p("The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates")),
                        gt.t_pargs(option=c_char_p("-rcon"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(rconstr)),
                                   desc=c_char_p("Maximum distance for P-LINCS (nm), 0 is estimate")),
                        gt.t_pargs(option=c_char_p("-dlb"),
                                   bSet=gt.FALSE,
                                   type=gt.etENUM,
                                   u=gt.u_u(addressof(dddlb_opt)),
                                   desc=c_char_p("Dynamic load balancing (with DD)")),
                        gt.t_pargs(option=c_char_p("-dds"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(dlb_scale)),
                                   desc=c_char_p("Minimum allowed dlb scaling of the DD cell size")),
                        gt.t_pargs(option=c_char_p("-ddcsx"),
                                   bSet=gt.FALSE,
                                   type=gt.etSTR,
                                   u=gt.u_u(addressof(ddcsx)),
                                   desc=c_char_p("HIDDENThe DD cell sizes in x")),
                        gt.t_pargs(option=c_char_p("-ddcsy"),
                                   bSet=gt.FALSE,
                                   type=gt.etSTR,
                                   u=gt.u_u(addressof(ddcsy)),
                                   desc=c_char_p("HIDDENThe DD cell sizes in y")),
                        gt.t_pargs(option=c_char_p("-ddcsz"),
                                   bSet=gt.FALSE,
                                   type=gt.etSTR,
                                   u=gt.u_u(addressof(ddcsz)),
                                   desc=c_char_p("HIDDENThe DD cell sizes in z")),
                        gt.t_pargs(option=c_char_p("-sum"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bSumEner)),
                                   desc=c_char_p("Sum the energies at every step")),
                        gt.t_pargs(option=c_char_p("-v"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bVerbose)),
                                   desc=c_char_p("Be loud and noisy")),
                        gt.t_pargs(option=c_char_p("-compact"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bCompact)),
                                   desc=c_char_p("Write a compact log file")),
                        gt.t_pargs(option=c_char_p("-seppot"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bSepPot)),
                                   desc=c_char_p("Write separate V and dVdl terms for each interaction type and node to the log file(s)")),
                        gt.t_pargs(option=c_char_p("-pforce"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(pforce)),
                                   desc=c_char_p("Print all forces larger than this (kJ/mol nm)")),
                        gt.t_pargs(option=c_char_p("-reprod"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bReproducible)),
                                   desc=c_char_p("Try to avoid optimizations that affect binary reproducibility")),
                        gt.t_pargs(option=c_char_p("-cpt"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(cpt_period)),
                                   desc=c_char_p("Checkpoint interval (minutes)")),
                        gt.t_pargs(option=c_char_p("-append"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bAppendFiles)),
                                   desc=c_char_p("Append to previous output files when continuing from checkpoint")),
                        gt.t_pargs(option=c_char_p("-addpart"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bAddPart)),
                                   desc=c_char_p("Add the simulation part number to all output files when continuing from checkpoint")),
                        gt.t_pargs(option=c_char_p("-maxh"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(max_hours)),
                                   desc=c_char_p("Terminate after 0.99 times this time (hours)")),
                        gt.t_pargs(option=c_char_p("-multi"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(nmultisim)),
                                   desc=c_char_p("Do multiple simulations in parallel")),
                        gt.t_pargs(option=c_char_p("-replex"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(repl_ex_nst)),
                                   desc=c_char_p("Attempt replica exchange every # steps")),
                        gt.t_pargs(option=c_char_p("-reseed"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(repl_ex_seed)),
                                   desc=c_char_p("Seed for replica exchange, -1 is generate a seed")),
                        gt.t_pargs(option=c_char_p("-rerunvsite"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bRerunVSite)),
                                   desc=c_char_p("HIDDENRecalculate virtual site coordinates with -rerun")),
                        gt.t_pargs(option=c_char_p("-glas"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bGlas)),
                                   desc=c_char_p("Do glass simulation with special long range corrections")),
                        gt.t_pargs(option=c_char_p("-ionize"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bIonize)),
                                   desc=c_char_p("Do a simulation including the effect of an X-Ray bombardment on your system")),
                        gt.t_pargs(option=c_char_p("-confout"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bConfout)),
                                   desc=c_char_p("HIDDENWrite the last configuration with -c")),
                        gt.t_pargs(option=c_char_p("-stepout"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(nstepout)),
                                   desc=c_char_p("HIDDENFrequency of writing the remaining runtime"))
                        )

    ed             = gt.gmx_edsam_t()
    Flags          = c_ulong()
    PCA_Flags      = c_ulong()
    ddxyz          = ivec()
    dd_node_order  = c_int()
    HaveCheckpoint = c_int()
    fplog          = gt.FILE_p()
    fptes          = gt.FILE_p()
    sim_part       = c_int()
    suffixarray    = c_char*gt.STRLEN
    suffix         = suffixarray()

############################
#   Quick hack: No options #
############################
    argv_array = c_char_p*1
    argv       = argv_array(c_char_p("mdrun"))
    argc       = c_int(1)

    init_par         = libmdrun.init_par
    init_par.restype = pt_commrec
    cr               = (init_par(byref(argc),byref(argv)))

    if(commrec.MASTER(cr)):
        libmdrun.CopyRight(stderr,argv[0])

    PCA_Flags = (statutil.PCA_KEEP_ARGS | statutil.PCA_NOEXIT_ON_ARGS | statutil.PCA_CAN_SET_DEFFNM | (0 if commrec.MASTER(cr) else statutil.PCA_QUIET))
    if (not c_int.in_dll(libmdrun,"gmx_parallel_env").value):
        PCA_Flags = (PCA_Flags | statutil.PCA_BE_NICE)
    PCA_Flags = c_ulong(PCA_Flags)

    libmdrun.parse_common_args(byref(argc),\
                               argv,\
                               PCA_Flags,\
                               NFILE,\
                               fnm,\
                               c_int(len(pa)),\
                               pa,\
                               c_int(len(desc)),\
                               desc,\
                               c_int(0),\
                               POINTER(c_char_p)())

    dd_node_order = c_int(libmdrun.nenum(ddno_opt))
    cr.contents.npmenodes = npme

    if(nthreads.value>1):
        print "GROMACS compiled without threads support - can only use one thread"
        sys.exit(1)

    if(repl_ex_nst.value!=0 and nmultisim.value<2):
        print "Need at least two replicas for replica exchange (option -multi)"
        sys.exit(1)

    if(nmultisim.value>1 and commrec.PAR(cr)):
        libmdrun.init_multisystem(cr,nmultisim,NFILE,fnm,c_int(gt.TRUE))

#  /* Check if there is ANY checkpoint file available */
    sim_part = c_int(1)
    if(libmdrun.opt2bSet(c_char_p("-cpi"),NFILE,fnm)):
        sim_part = libmdrun.read_checkpoint_simulation_part(libmdrun.opt2fn("-cpi",NFILE,fnm)) + 1
#       /* sim_part will now be 1 if no checkpoint file was found */
        if(sim_part==1 and commrec.MASTER(cr)):
            print "No previous checkpoint file present, assuming this is a new run.\n"

    if(sim_part.value<=1):
        bAppendFiles = c_int(gt.FALSE)

    if((not bAppendFiles.value) and bAddPart.value and sim_part.value>1):
#        /* This is a continuation run, rename trajectory output files (except checkpoint files) */
#        /* create new part name first (zero-filled) */
        if(sim_part.value<10):
            libc.sprintf(suffix,"part000%d",sim_part)
        elif(sim_part.value<100):
            libc.sprintf(suffix,"part00%d",sim_part)
        elif(sim_part.value<1000):
            libc.sprintf(suffix,"part0%d",sim_part)
        else:
            libc.sprintf(suffix,"part%d",sim_part)
        libmdrun.add_suffix_to_output_names(fnm,NFILE,suffix)
        libc.fprintf(stdout,"Checkpoint file is from part %d, new output files will be suffixed %s.\n",c_int(sim_part.value-1),suffix)

    Flags = MD_RERUN if(libmdrun.opt2bSet(c_char_p("-rerun"),NFILE,fnm)) else 0
    Flags = Flags | (MD_SEPPOT if(bSepPot.value) else 0)
    Flags = Flags | (MD_IONIZE if(bIonize.value) else 0)
    Flags = Flags | (MD_GLAS if(bGlas.value) else 0)
    Flags = Flags | (MD_PARTDEC if(bPartDec.value) else 0)
    Flags = Flags | (MD_DDBONDCHECK if(bDDBondCheck.value) else 0)
    Flags = Flags | (MD_DDBONDCOMM if(bDDBondComm.value) else 0)
    Flags = Flags | (MD_CONFOUT if(bConfout.value) else 0)
    Flags = Flags | (MD_NOGSTAT if(not bSumEner.value) else 0)
    Flags = Flags | (MD_RERUN_VSITE if(bRerunVSite.value) else 0)
    Flags = Flags | (MD_REPRODUCIBLE if(bReproducible.value) else 0)
    Flags = Flags | (MD_APPENDFILES if(bAppendFiles.value) else 0)
    Flags = Flags | (MD_STARTFROMCPT if(sim_part.value>1) else 0)
    Flags = c_ulong(Flags)

#  /* We postpone opening the log file if we are appending, so we can first truncate
#   * the old log file and append to the correct position there instead.
#   */
    if(commrec.MASTER(cr) and (not bAppendFiles.value)):
        fplog = gt.FILE_p(libmdrun.gmx_log_open(c_char_p(libmdrun.ftp2fn(gt.efLOG,NFILE,fnm)),cr,c_int(not bSepPot.value),Flags))
        libmdrun.CopyRight(fplog,argv[0])
        libmdrun.please_cite(fplog,c_char_p("Hess2008b"))
        libmdrun.please_cite(fplog,c_char_p("Spoel2005a"))
        libmdrun.please_cite(fplog,c_char_p("Lindahl2001a"))
        libmdrun.please_cite(fplog,c_char_p("Berendsen95a"))
    else:
        fplog = gt.FILE_p()

#  /* Essential dynamics */
    if (libmdrun.opt2bSet(c_char_p("-ei"),NFILE,fnm)):
#   /* Open input and output files, allocate space for ED data structure */
        ed = libmdrun.ed_open(NFILE,fnm,cr)
    else:
        ed = gt.gmx_edsam_t()

    ddxyz[gt.XX] = c_int(int(realddxyz[gt.XX] + 0.5))
    ddxyz[gt.YY] = c_int(int(realddxyz[gt.YY] + 0.5))
    ddxyz[gt.ZZ] = c_int(int(realddxyz[gt.ZZ] + 0.5))

#    mc = gt.mdrunner_comm()
#    mc.state       = gt.NULL
#    mc.buf         = gt.NULL
#    mc.f           = gt.NULL
#    mc.mtop        = gt.NULL
#    mc.mdatoms     = gt.NULL
#    mc.fr          = gt.NULL
#    mc.fcd         = gt.NULL
#    mc.ewaldcoeff  = c_real(0)
#    mc.pmedata     = gt.NULL
#    mc.vsite       = gt.NULL
#    mc.nsteps_done = c_int(0)
#    libmdrun.mdrunner_initialize(fplog,cr,NFILE,fnm,bVerbose,
#                                 ddxyz,dd_node_order,rdd,rconstr,
#                                 c_char_p(dddlb_opt[0]),dlb_scale,ddcsx,ddcsy,ddcsz,
#                                 pforce,Flags,
#                                 byref(mc))
#    mdrunner_integrate          = libmdrun.mdrunner_integrate
#    mdrunner_integrate.restypte = gt.time_t
#    start_t                     = gt.time_t(mdrunner_integrate(fplog,cr,NFILE,fnm,bVerbose,bCompact,
#                                                              dlb_scale,
#                                                              nstepout,ed,repl_ex_nst,repl_ex_seed,
#                                                              cpt_period,max_hours,Flags,
#                                                              byref(mc)))
#    libmdrun.mdrunner_finalize(fplog,cr,NFILE,fnm,Flags,
#                               start_t,byref(mc))

#    libmdrun.mdrunner(fplog,cr,NFILE,fnm,bVerbose,bCompact,
#                      ddxyz,dd_node_order,rdd,rconstr,
#                      c_char_p(dddlb_opt[0]),dlb_scale,ddcsx,ddcsy,ddcsz,
#                      nstepout,ed,repl_ex_nst,repl_ex_seed,pforce,
#                      cpt_period,max_hours,Flags)

    mdrunner(fplog,cr,NFILE,fnm,bVerbose,bCompact,
             ddxyz,dd_node_order,rdd,rconstr,
             c_char_p(dddlb_opt[0]),dlb_scale,ddcsx,ddcsy,ddcsz,
             nstepout,ed,repl_ex_nst,repl_ex_seed,pforce,
             cpt_period,max_hours,Flags)
    if(c_int.in_dll(libmdrun,"gmx_parallel_env").value):
        libmdrun.gmx_finalize(cr)

    if(commrec.MULTIMASTER(cr)):
        libmdrun.thanx(stderr)
#    /* Log file has to be closed in mdrunner if we are appending to it (fplog not set here) */
        if(commrec.MASTER(cr) and (not bAppendFiles.value)):
            libmdrun.gmx_log_close(fplog)

    return


do_md            = libmdrun.do_md
do_md.restype    = gt.time_t
do_steep         = libmdrun.do_steep
do_steep.restype = gt.time_t
do_cg            = libmdrun.do_cg
do_cg.restype    = gt.time_t
do_nm            = libmdrun.do_nm
do_nm.restype    = gt.time_t
do_lbfgs         = libmdrun.do_lbfgs
do_lbfgs.restype = gt.time_t
do_tpi           = libmdrun.do_tpi
do_tpi.restype   = gt.time_t
#/* The array should match the eI array in include/types/enums.h */
#gmx_intp_t_array = gmx_intp_t * enums.eiNR.value
integrator       = [do_md,
                    do_steep,
                    do_cg,
                    do_md,
                    do_md,
                    do_nm,
                    do_lbfgs,
                    do_tpi,
                    do_tpi,
                    do_md]
assert len(integrator)==enums.eiNR.value, "len(integrator)!=enums.eiNR.value !!"
#from md.c
#typedef struct {
#  gmx_integrator_t *func;
#} gmx_intp_t;
#/* The array should match the eI array in include/types/enums.h */
#const gmx_intp_t integrator[eiNR] = { {do_md}, {do_steep}, {do_cg}, {do_md}, {do_md}, {do_nm}, {do_lbfgs}, {do_tpi}, {do_tpi}, {do_md} };

def mdrunner(fplog,cr,nfile,fnm,bVerbose,bCompact,
             ddxyz,dd_node_order,rdd,rconstr,
             dddlb_opt,dlb_scale,ddcsx,ddcsy,ddcsz,
             nstepout,ed,repl_ex_nst,repl_ex_seed,pforce,
             cpt_period,max_hours,Flags):

    nodetime         = c_double(0.0)
    realtime         = c_double
    inputrec         = POINTER(gt.t_inputrec)
    state            = POINTER(gt.t_state)()
    box              = matrix
    buf              = POINTER(rvec)()
    f                = POINTER(rvec)()
    tmpr1            = c_real
    tmpr2            = c_real
    nrnb             = POINTER(gt.t_nrnb)
    mtop             = POINTER(gt.gmx_mtop_t)()
    mdatoms          = POINTER(gt.t_mdatoms)()
    fr               = POINTER(gt.t_forcerec)()
    fcd              = POINTER(gt.t_fcdata)()
    ewaldcoeff       = c_real(0.0)
    pmedata          = c_void_p()
#    pmedata          = POINTER(gt.gmx_pme_t)
    start_t          = gt.time_t(0)
    vsite            = POINTER(gt.gmx_vsite_t)()
    constr           = gt.gmx_constr_t
    i                = c_int
    m                = c_int
    nChargePerturbed = c_int(-1)
    status           = c_int
    nalloc           = c_int
    gro              = c_char_p
    wcycle           = gt.gmx_wallcycle_t
    bReadRNG         = c_int()
    bReadEkin        = c_int()
    list             = c_int
    nsteps_done      = c_int()

    snew         = libmdrun.save_calloc
    snew.restype = POINTER(gt.t_inputrec)
    inputrec     = snew(c_char_p("inputrec"),\
                        c_char_p("__FILE__"),\
                        c_int(0),\
                        c_uint(1),\
                        c_size_t(sizeof(gt.t_inputrec)))
    snew.restype = POINTER(gt.gmx_mtop_t)
    mtop         = snew(c_char_p("mtop"),\
                        c_char_p("__FILE__"),\
                        c_int(0),\
                        c_uint(1),\
                        c_size_t(sizeof(gt.gmx_mtop_t)))

    if(bVerbose.value and commrec.SIMMASTER(cr)):
        libc.fprintf(stderr,"Getting Loaded...\n")

    if(Flags.value & MD_APPENDFILES):
        fplog = gt.FILE_p()

    if(commrec.PAR(cr)):
#        /* The master thread on the master node reads from disk,
#        * then distributes everything to the other processors.
#        */
        list = c_int((LIST_SCALARS | LIST_INPUTREC) if (commrec.SIMMASTER(cr) and (not (Flags.value & MD_APPENDFILES))) else 0)
        snew.restype = POINTER(gt.t_state)
        state        = snew(c_char_p("state"),\
                            c_char_p("__FILE__"),\
                            c_int(0),\
                            c_uint(1),\
                            c_size_t(sizeof(gt.t_state)))
        ftp2fn         = libmdrun.ftp2fn
        ftp2fn.restype = c_char_p
        libmdrun.init_parallel(fplog,ftp2fn(gt.efTPX,nfile,fnm),\
                               cr,inputrec,mtop,state,list)
    else:
#        /* Read a file for a single processor */
        snew.restype = POINTER(gt.t_state)
        state        = snew(c_char_p("state"),\
                            c_char_p("__FILE__"),\
                            c_int(0),\
                            c_uint(1),\
                            c_size_t(sizeof(gt.t_state)))
        ftp2fn         = libmdrun.ftp2fn
        ftp2fn.restype = c_char_p
        libmdrun.init_single(fplog,inputrec,ftp2fn(gt.efTPX,nfile,fnm),\
                             mtop,state)
    if((not enums.EEL_PME(inputrec.contents.coulombtype)) or (Flags & MD_PARTDEC)):
        cr.contents.npmenodes = c_int(0)

#  /* NMR restraints must be initialized before load_checkpoint,
#   * since with time averaging the history is added to t_state.
#   * For proper consistency check we therefore need to extend
#   * t_state here.
#   * So the PME-only nodes (if present) will also initialize
#   * the distance restraints.
#   */
    snew.restype = POINTER(gt.t_fcdata)
    fcd          = snew(c_char_p("fcd"),\
                        c_char_p("__FILE__"),\
                        c_int(0),\
                        c_uint(1),\
                        c_size_t(sizeof(gt.t_fcdata)))

#    /* This needs to be called before read_checkpoint to extend the state */
    libmdrun.init_disres(fplog,mtop,inputrec,cr,c_int(Flags.value & MD_PARTDEC),fcd,state)

    if(libmdrun.gmx_mtop_ftype_count(mtop,gt.F_ORIRES)>0):
        if(commrec.PAR(cr) and  (not (Flags & MD_PARTDEC))):
            print "Orientation restraints do not work (yet) with domain decomposition, use particle decomposition (mdrun option -pd)"
            sys.exit(1)
#        /* Orientation restraints */
        if(commrec.MASTER(cr)):
            libmdrun.init_orires(fplog,mtop,state.contents.x,inputrec,cr.contents.ms,byref(fcd.contents.orires),state)

    if(inptrec.DEFORM(inputrec.contents)):
#        /* Store the deform reference box before reading the checkpoint */
        if(commrec.SIMMASTER(cr)):
            box = vec.copy_mat(state.contents.box)
        if(commrec.PAR(cr)):
            libmdrun.gmx_bcast(sizeof(box),box,cr)
        libmdrun.set_deform_reference_box(inputrec.contents.init_step,box);

    if(libmdrun.opt2bSet(c_char_p("-cpi"),nfile,fnm)):
#      /* Check if checkpoint file exists before doing continuation.
#       * This way we can use identical input options for the first and subsequent runs...
#       */
        if(libmdrun.fexist(c_int(libmdrun.opt2fn(c_char_p("-cpi"),nfile,fnm)))):
            libmdrun.load_checkpoint(c_int(libmdrun.opt2fn(c_char_p("-cpi"),nfile,fnm)),fplog,
                                     cr,c_int(Flags.value & MD_PARTDEC),ddxyz,
                                     inputrec,state,byref(bReadRNG),byref(bReadEkin),
                                     c_int(Flags.value & MD_APPENDFILES))
            if(bReadRNG.value):
                Flags = c_ulong(Flags.value | MD_READ_RNG)
            if(bReadEkin.value):
                Flags = c_ulong(Flags.value | MD_READ_EKIN)

    if(commrec.MASTER(cr) and (Flags.value & MD_APPENDFILES)):
        fplog = libmdrun.gmx_log_open(libmdrun.ftp2fn(gt.efLOG,nfile,fnm),\
                                      cr, c_int(not (Flags.value & MD_SEPPOT)),Flags);
    if(commrec.SIMMASTER(cr)):
        box = vec.copy_mat(state.contents.box)

    if(commrec.PAR(cr)):
        libmdrun.gmx_bcast(c_int(sizeof(box)),box,cr)
#
    if(bVerbose.value and commrec.SIMMASTER(cr)):
        libc.fprintf(stderr,"Loaded with Money\n\n")

    if (commrec.PAR(cr) and (not ((Flags.value & MD_PARTDEC) or enums.EI_TPI(inputrec.contents.eI)))):
        init_domain_decomposition = libmdrun.init_domain_decomposition
        init_domain_decomposition.restype = POINTER(gt.gmx_domdec_t)
        cr.contents.dd = init_domain_decomposition(fplog,cr,Flags,ddxyz,rdd,rconstr,
                                                   dddlb_opt,dlb_scale,
                                                   ddcsx,ddcsy,ddcsz,
                                                   mtop,inputrec,
                                                   box,state.contents.x)
        libmdrun.make_dd_communicators(fplog,cr,dd_node_order)
#        /* Set overallocation to avoid frequent reallocation of arrays */
        libmdrun.set_over_alloc_dd(c_int(gt.TRUE))
    else:
        cr.contents.duty = c_int(commrec.DUTY_PP | commrec.DUTY_PME)
        if(inputrec.contents.ePBC==enums.epbcSCREW.value):
            print "pbc=%s is only implemented with domain decomposition" % (names.epbc_names[inputrec.contents.ePBC])
            sys.exit(1)

    if(commrec.PAR(cr)):
#        /* After possible communicator splitting in make_dd_communicators.
#        * we can set up the intra/inter node communication.
#        */
        libmdrun.gmx_setup_nodecomm(fplog,cr)

    wallcycle_init         = libmdrun.wallcycle_init
    wallcycle_init.restype = gt.gmx_wallcycle_t
    wcycle                 = wallcycle_init(fplog,cr)

    snew.restype = POINTER(gt.t_nrnb)
    nrnb         = snew(c_char_p("nrnb"),\
                        c_char_p("__FILE__"),\
                        c_int(0),\
                        c_uint(1),\
                        c_size_t(sizeof(gt.t_nrnb)))

    if(cr.contents.duty & commrec.DUTY_PP):
#        /* For domain decomposition we allocate dynamically
#        * in dd_partition_system.
#        */
        if(commrec.DOMAINDECOMP(cr)):
            libmdrun.bcast_state_setup(cr,state)
        else:
            if(commrec.PAR(cr)):
                if(not (commrec.MASTER(cr))):
                    snew.restype = POINTER(gt.t_state)
                    state        = snew(c_char_p("state"),\
                                        c_char_p("__FILE__"),\
                                        c_int(0),\
                                        c_uint(1),\
                                        c_size_t(sizeof(gt.t_state)))
                libmdrun.bcast_state(cr,state,c_int(gt.TRUE))
            snew.restype = POINTER(rvec)
            buf          = snew(c_char_p("buf"),\
                                c_char_p("__FILE__"),\
                                c_int(0),\
                                c_uint(mtop.contents.natoms),\
                                c_size_t(sizeof(rvec)))
            snew.restype = POINTER(rvec)
            f            = snew(c_char_p("f"),\
                                c_char_p("__FILE__"),\
                                c_int(0),\
                                c_uint(mtop.contents.natoms),\
                                c_size_t(sizeof(rvec)))
#       /* Dihedral Restraints */
        if(libmdrun.gmx_mtop_ftype_count(mtop,gt.F_DIHRES) > 0):
            libmdrun.init_dihres(fplog,mtop,inputrec,fcd)
#       /* Initiate forcerecord */
        mk_forcerec         = libmdrun.mk_forcerec
        mk_forcerec.restype = POINTER(gt.t_forcerec)
        fr                  = mk_forcerec()
        libmdrun.init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,c_int(gt.FALSE),
                               libmdrun.opt2fn(c_char_p("-table"),nfile,fnm),libmdrun.opt2fn(c_char_p("-tablep"),nfile,fnm),
                               libmdrun.opt2fn(c_char_p("-tableb"),nfile,fnm),c_int(gt.FALSE),pforce)
        fr.contents.bSepDVDL = c_int((Flags.value & MD_SEPPOT)==MD_SEPPOT)

#       /* Initialize QM-MM */
        if(fr.contents.bQMMM):
            libmdrun.init_QMMMrec(cr,box,mtop,inputrec,fr)

#       /* Initialize the mdatoms structure.
#       * mdatoms is not filled with atom data,
#       * as this can not be done now with domain decomposition.
#       */
        init_mdatoms         = libmdrun.init_mdatoms
        init_mdatoms.restype = POINTER(gt.t_mdatoms)
        mdatoms              = init_mdatoms(fplog,mtop,(inputrec.contents.efep != enums.efepNO.value))
        print "mdatoms ",mdatoms

#       /* Initialize the virtual site communication */
        init_vsite         = libmdrun.init_vsite
        init_vsite.restype = POINTER(gt.gmx_vsite_t)
        vsite              = init_vsite(mtop,cr)

        libmdrun.calc_shifts(box,fr.contents.shift_vec)

#       /* With periodic molecules the charge groups should be whole at start up
#       * and the virtual sites should not be far from their proper positions.
#       */
        if((not inputrec.contents.bContinuation and commrec.MASTER(cr)) and
           (not (inputrec.contents.ePBC != enums.epbcNONE and inputrec.contents.bPeriodicMols))):
#           /* Make molecules whole at start of run */
            if(fr.contents.ePBC != enums.epbcNONE.value):
                libmdrun.do_pbc_first_mtop(fplog,inputrec.contents.ePBC,box,mtop,state.contents.x)
            if(bool(vsite)):
#              /* Correct initial vsite positions are required
#               * for the initial distribution in the domain decomposition
#               * and for the initial shell prediction.
#               */
                libmdrun.construct_vsites_mtop(fplog,vsite,mtop,state.contents.x)

#       /* Initiate PPPM if necessary */
#        LOGNAME = c_char_p()
#        getenv = libc.getenv
#        getenv.restype = c_char_p
        if(fr.contents.eeltype == enums.eelPPPM.value):
            if(mdatoms.contents.nChargePerturbed):
                print "Free energy with %s is not implemented", (names.eel_names[fr.contents.eeltype])
                sys.exit(1)
            gmx_pppm_init         = libmdrun.gmx_pppm_init
            gmx_pppm_init.restype = c_int
            status                = gmx_pppm_init(fplog,cr,c_int(gt.FALSE),c_int(gt.TRUE),
                                                  box,c_char_p(libc.getenv("GMXGHAT")),
                                                  inputrec,
                                                  (Flags.value & MD_REPRODUCIBLE))
            if(status.value != 0):
                print "Error %d initializing PPPM" % status.value
                sys.exit(1)

        if(enums.EEL_PME(fr.contents.eeltype)):
            ewaldcoeff = c_real(fr.contents.ewaldcoeff)
            pmedata    = c_void_p(fr.contents.pmedata)
        else:
            pmedata = c_void_p()
    else:
#       /* This is a PME only node */
#       /* We don't need the state */
        libmdrun.done_state(state)

        calc_ewaldcoeff         = libmdrun.calc_ewaldcoeff
        calc_ewaldcoeff.restype = c_real
        ewaldcoeff              = c_real(calc_ewaldcoeff(inputrec.contents.rcoulomb, inputrec.contents.ewald_rtol))
        snew.restype = c_void_p
        pmdeata      = snew(c_char_p("pmedata"),\
                            c_char_p("__FILE__"),\
                            c_int(0),\
                            c_uint(1),\
                            c_size_t(sizeof(c_void_p)))

#   /* Initiate PME if necessary,
#   * either on all nodes or on dedicated PME nodes only. */
    if(enums.EEL_PME(inputrec.contents.coulombtype)):
        if(bool(mdatoms)):
            nChargePerturbed = c_int(mdatoms.contents.nChargePerturbed)
        if (cr.contents.npmenodes > 0):
#           /* The PME only nodes need to know nChargePerturbed */
            libmdrun.gmx_bcast_sim(c_int(sizeof(nChargePerturbed)),byref(nChargePerturbed),cr)
        if(cr.contents.duty & commrec.DUTY_PME):
            gmx_pme_init         = libmdrun.gmx_pme_init
            gmx_pme_init.restype = c_int
            status               = gmx_pme_init(pmedata,cr,inputrec,
                                                c_int(mtop.contents.natoms if bool(mtop) else 0),nChargePerturbed,
                                                c_int(Flags.value & MD_REPRODUCIBLE))
            if(status.value != 0):
                print "Error %d initializing PME",status.value
                sys.exit(1)

    if(integrator[inputrec.contents.eI] == do_md):
#      /* Turn on signal handling on all nodes */
#      /*
#       * (A user signal from the PME nodes (if any)
#       * is communicated to the PP nodes.
#       */
        debug = gt.FILE_p.in_dll(libmdrun,"debug")
        if(c_char_p(libc.getenv("GMX_NO_TERM")).value == gt.NULL):
            if(bool(debug)):
                libc.fprintf(debug,"Installing signal handler for SIGTERM\n")
#            Don't konw how to do the following...'
#            libc.signal(c_int(signal.SIGTERM),libmdrun.signal_handler)
        if (c_char_p(libc.getenv("GMX_NO_USR1")).value == gt.NULL):
            if (bool(debug)):
                libc.fprintf(debug,"Installing signal handler for SIGUSR1\n")
#            Don't konw how to do the following...'
#            libc.signal(c_int(signal.SIGUSR1),libmdrun.signal_handler)

    if(cr.contents.duty & commrec.DUTY_PP):
        if(inputrec.contents.ePull!=enums.epullNO.value):
#           /* Initialize pull code */
            libmdrun.init_pull(fplog,inputrec,nfile,fnm,mtop,cr,
                               c_int(enums.EI_DYNAMICS(inputrec.contents.eI) and commrec.MASTER(cr)),
                               Flags)

        init_constraints         = libmdrun.init_constraints
        init_constraints.restype = gt.gmx_constr_t
        constr                   = init_constraints(fplog,mtop,inputrec,ed,state,cr)

        if(commrec.DOMAINDECOMP(cr)):
            libmdrun.dd_init_bondeds(fplog,cr.contents.dd,mtop,vsite,constr,inputrec,
                                     c_int(Flags.value & MD_DDBONDCHECK),fr.contents.cginfo_global)
            libmdrun.set_dd_parameters(fplog,cr.contents.dd,dlb_scale,inputrec,fr,box)
            libmdrun.setup_dd_grid(fplog,cr.contents.dd)

#       /* Now do whatever the user wants us to do (how flexible...) */
        start_t = integrator[inputrec.contents.eI](fplog,cr,nfile,fnm,
                                                   bVerbose,bCompact,
                                                   vsite,constr,
                                                   nstepout,inputrec,mtop,
                                                   fcd,state,f,buf,
                                                   mdatoms,nrnb,wcycle,ed,fr,
                                                   repl_ex_nst,repl_ex_seed,
                                                   cpt_period,max_hours,
                                                   Flags,
                                                   byref(nsteps_done))
        if(inputrec.contents.ePull!=enums.epullNO.value):
            libmdrun.finish_pull(fplog,inputrec.contents.pull)
    else:
#       /* do PME only */
        libmdrun.gmx_pmeonly(pmedata.contents,cr,nrnb,wcycle,ewaldcoeff,c_int(gt.FALSE))

#   /* Some timing stats */
#    print c_double.in_dll(libmdrun,"nodetime")
    if(commrec.MASTER(cr)):
        realtime = su.difftime(gt.time_t(libc.time(gt.NULL)),gt.time_t(start_t))
        nodetime = c_double(libmdrun.node_time())
        if(nodetime.value == 0.0):
#           nodetime accounting is not accessible from library calls :-(
            nodetime = realtime
    else:
        realtime=c_double(0.0)

    libmdrun.wallcycle_stop(wcycle,wc.ewcRUN)

#   /* Finish up, write some stuff
#    * if rerunMD, don't write last frame again
#    */
    libmdrun.finish_run(fplog,cr,c_char_p(libmdrun.ftp2fn(gt.efSTO,nfile,fnm)),
                        inputrec,nrnb,wcycle,nodetime,realtime,nsteps_done,
                        c_int(enums.EI_DYNAMICS(c_int(inputrec.contents.eI)) and (not bool(commrec.MULTISIM(cr)))))

#   /* Does what it says */
    libmdrun.print_date_and_time(fplog,cr.contents.nodeid,"Finished mdrun")

#   /* Close logfile already here if we were appending to it */
    if(commrec.MASTER(cr) and (Flags.value & MD_APPENDFILES)):
        libmdrun.gmx_log_close(fplog)

    return

def CheckSizesOfStructs():
    print "sizeof(gt.t_block)                    ",sizeof(gt.t_block)
    print "sizeof(gt.t_blocka)                   ",sizeof(gt.t_blocka)
    print "sizeof(gt.t_symbuf)                   ",sizeof(gt.t_symbuf)
    print "sizeof(gt.t_symtab)                   ",sizeof(gt.t_symtab)
    print "sizeof(gt.t_iparams)                  ",sizeof(gt.t_iparams)
#    print "sizeof(gt.t_iparams.bham)             ",sizeof(gt.t_iparams.bham)
#    print "sizeof(gt.t_iparams.harmonic)         ",sizeof(gt.t_iparams.harmonic)
#    print "sizeof(gt.t_iparams.cubic)            ",sizeof(gt.t_iparams.cubic)
#    print "sizeof(gt.t_iparams.fene)             ",sizeof(gt.t_iparams.fene)
#    print "sizeof(gt.t_iparams.cross_bb)         ",sizeof(gt.t_iparams.cross_bb)
#    print "sizeof(gt.t_iparams.cross_ba)         ",sizeof(gt.t_iparams.cross_ba)
#    print "sizeof(gt.t_iparams.u_b)              ",sizeof(gt.t_iparams.u_b)
#    print "sizeof(gt.t_iparams.qangle)           ",sizeof(gt.t_iparams.qangle)
#    print "sizeof(gt.t_iparams.polarize)         ",sizeof(gt.t_iparams.polarize)
#    print "sizeof(gt.t_iparams.wpol)             ",sizeof(gt.t_iparams.wpol)
#    print "sizeof(gt.t_iparams.thole)            ",sizeof(gt.t_iparams.thole)
#    print "sizeof(gt.t_iparams.lj)               ",sizeof(gt.t_iparams.lj)
#    print "sizeof(gt.t_iparams.lj14)             ",sizeof(gt.t_iparams.lj14)
#    print "sizeof(gt.t_iparams.ljc14)            ",sizeof(gt.t_iparams.ljc14)
#    print "sizeof(gt.t_iparams.ljcnb)            ",sizeof(gt.t_iparams.ljcnb)
#    print "sizeof(gt.t_iparams.pdihs)            ",sizeof(gt.t_iparams.pdihs)
#    print "sizeof(gt.t_iparams.constr)           ",sizeof(gt.t_iparams.constr)
#    print "sizeof(gt.t_iparams.settle)           ",sizeof(gt.t_iparams.settle)
#    print "sizeof(gt.t_iparams.morse)            ",sizeof(gt.t_iparams.morse)
#    print "sizeof(gt.t_iparams.posres)           ",sizeof(gt.t_iparams.posres)
#    print "sizeof(gt.t_iparams.rbdihs)           ",sizeof(gt.t_iparams.rbdihs)
#    print "sizeof(gt.t_iparams.vsite)            ",sizeof(gt.t_iparams.vsite)
#    print "sizeof(gt.t_iparams.vsiten)           ",sizeof(gt.t_iparams.vsiten)
#    print "sizeof(gt.t_iparams.disres)           ",sizeof(gt.t_iparams.disres)
#    print "sizeof(gt.t_iparams.dihres)           ",sizeof(gt.t_iparams.dihres)
#    print "sizeof(gt.t_iparams.orires)           ",sizeof(gt.t_iparams.orires)
#    print "sizeof(gt.t_iparams.tab)              ",sizeof(gt.t_iparams.tab)
#    print "sizeof(gt.t_iparams.generic)          ",sizeof(gt.t_iparams.generic)
    print "sizeof(gt.t_ilist)                    ",sizeof(gt.t_ilist)
    print "sizeof(gt.t_idef)                     ",sizeof(gt.t_idef)
    print "sizeof(gt.t_atom)                     ",sizeof(gt.t_atom)
    print "sizeof(gt.t_pdbinfo)                  ",sizeof(gt.t_pdbinfo)
    print "sizeof(gt.t_grps)                     ",sizeof(gt.t_grps)
    print "sizeof(gt.t_atoms)                    ",sizeof(gt.t_atoms)
    print "sizeof(gt.t_atomtypes)                ",sizeof(gt.t_atomtypes)
    print "sizeof(gt.t_topology)                 ",sizeof(gt.t_topology)
    print "sizeof(gt.t_cosines)                  ",sizeof(gt.t_cosines)
    print "sizeof(gt.t_grpopts)                  ",sizeof(gt.t_grpopts)
    print "sizeof(gt.t_pullgrp)                  ",sizeof(gt.t_pullgrp)
    print "sizeof(gt.t_pull)                     ",sizeof(gt.t_pull)
    print "sizeof(gt.t_inputrec)                 ",sizeof(gt.t_inputrec)
    print "sizeof(gt.est_names)                  ",sizeof(gt.est_names)
    print "sizeof(gt.history_t)                  ",sizeof(gt.history_t)
    print "sizeof(gt.ekinstate_t)                ",sizeof(gt.ekinstate_t)
    print "sizeof(gt.energyhistory_t)            ",sizeof(gt.energyhistory_t)
    print "sizeof(gt.t_state)                    ",sizeof(gt.t_state)
    print "sizeof(gt.t_nrnb)                     ",sizeof(gt.t_nrnb)
    print "sizeof(gt.gmx_cycles_t)               ",sizeof(gt.gmx_cycles_t)
    print "sizeof(gt.gmx_wallcycle)              ",sizeof(gt.gmx_wallcycle)
    print "sizeof(gt.gmx_wallcycle_t)            ",sizeof(gt.gmx_wallcycle_t)
    print "sizeof(gt.gmx_groups_t)               ",sizeof(gt.gmx_groups_t)
    print "sizeof(gt.gmx_molblock_t)             ",sizeof(gt.gmx_molblock_t)
    print "sizeof(gt.gmx_moltype_t)              ",sizeof(gt.gmx_moltype_t)
    print "sizeof(gt.gmx_ffparams_t)             ",sizeof(gt.gmx_ffparams_t)
    print "sizeof(gt.gmx_mtop_t)                 ",sizeof(gt.gmx_mtop_t)
    print "sizeof(gt.t_mdatoms)                  ",sizeof(gt.t_mdatoms)
    print "sizeof(gt.t_nblist)                   ",sizeof(gt.t_nblist)
    print "sizeof(gt.t_forcetable)               ",sizeof(gt.t_forcetable)
    print "sizeof(gt.t_nblists)                  ",sizeof(gt.t_nblists)
    print "sizeof(gt.gmx_fft)                    ",sizeof(gt.gmx_fft)
    print "sizeof(gt.gmx_fft_t)                  ",sizeof(gt.gmx_fft_t)
    print "sizeof(gt.t_fftgrid)                  ",sizeof(gt.t_fftgrid)
    print "sizeof(gt.pme_grid_comm_t)            ",sizeof(gt.pme_grid_comm_t)
    print "sizeof(gt.pme_overlap_t)              ",sizeof(gt.pme_overlap_t)
    print "sizeof(gt.pme_atomcomm_t)             ",sizeof(gt.pme_atomcomm_t)
    print "sizeof(gt.gmx_pme)                    ",sizeof(gt.gmx_pme)
    print "sizeof(gt.gmx_pme_t)                  ",sizeof(gt.gmx_pme_t)
    print "sizeof(gt.t_excl)                     ",sizeof(gt.t_excl)
    print "sizeof(gt.t_grid)                     ",sizeof(gt.t_grid)
    print "sizeof(gt.t_ns_buf)                   ",sizeof(gt.t_ns_buf)
    print "sizeof(gt.gmx_ns_t)                   ",sizeof(gt.gmx_ns_t)
    print "sizeof(gt.t_MMrec)                    ",sizeof(gt.t_MMrec)
    print "sizeof(gt.t_QMrec)                    ",sizeof(gt.t_QMrec)
    print "sizeof(gt.t_QMMMrec)                  ",sizeof(gt.t_QMMMrec)
    print "sizeof(gt.t_forcerec)                 ",sizeof(gt.t_forcerec)
    print "sizeof(gt.rvec5)                      ",sizeof(gt.rvec5)
    print "sizeof(gt.t_disresdata)               ",sizeof(gt.t_disresdata)
    print "sizeof(gt.t_oriresdata)               ",sizeof(gt.t_oriresdata)
    print "sizeof(gt.bondedtable_t)              ",sizeof(gt.bondedtable_t)
    print "sizeof(gt.t_fcdata)                   ",sizeof(gt.t_fcdata)
    print "sizeof(gt.gmx_pme_comm_n_box)         ",sizeof(gt.gmx_pme_comm_n_box)
    print "sizeof(gt.gmx_pme_comm_n_box_p_t)     ",sizeof(gt.gmx_pme_comm_n_box_p_t)
    print "sizeof(gt.gmx_cgsort_t)               ",sizeof(gt.gmx_cgsort_t)
    print "sizeof(gt.gmx_domdec_sort_t)          ",sizeof(gt.gmx_domdec_sort_t)
    print "sizeof(gt.gmx_ddpme_t)                ",sizeof(gt.gmx_ddpme_t)
    print "sizeof(gt.gmx_ddzone_t)               ",sizeof(gt.gmx_ddzone_t)
    print "sizeof(gt.gmx_domdec_ind_t)           ",sizeof(gt.gmx_domdec_ind_t)
    print "sizeof(gt.gmx_domdec_comm_dim_t)      ",sizeof(gt.gmx_domdec_comm_dim_t)
    print "sizeof(gt.gmx_domdec_root_t)          ",sizeof(gt.gmx_domdec_root_t)
    print "sizeof(gt.gmx_domdec_load_t)          ",sizeof(gt.gmx_domdec_load_t)
    print "sizeof(gt.gmx_domdec_comm)            ",sizeof(gt.gmx_domdec_comm)
    print "sizeof(gt.gmx_domdec_comm_p_t)        ",sizeof(gt.gmx_domdec_comm_p_t)
    print "sizeof(gt.gmx_domdec_master)          ",sizeof(gt.gmx_domdec_master)
    print "sizeof(gt.gmx_domdec_master_p_t)      ",sizeof(gt.gmx_domdec_master_p_t)
    print "sizeof(gt.gmx_reverse_ilist_t)        ",sizeof(gt.gmx_reverse_ilist_t)
    print "sizeof(gt.gmx_molblock_ind_t)         ",sizeof(gt.gmx_molblock_ind_t)
    print "sizeof(gt.gmx_reverse_top)            ",sizeof(gt.gmx_reverse_top)
    print "sizeof(gt.gmx_reverse_top_p_t)        ",sizeof(gt.gmx_reverse_top_p_t)
    print "sizeof(gt.gmx_specatsend_t)           ",sizeof(gt.gmx_specatsend_t)
    print "sizeof(gt.gmx_domdec_specat_comm)     ",sizeof(gt.gmx_domdec_specat_comm)
    print "sizeof(gt.gmx_domdec_specat_comm_p_t) ",sizeof(gt.gmx_domdec_specat_comm_p_t)
    print "sizeof(gt.gmx_domdec_constraints)     ",sizeof(gt.gmx_domdec_constraints)
    print "sizeof(gt.gmx_domdec_constraints_p_t) ",sizeof(gt.gmx_domdec_constraints_p_t)
    print "sizeof(gt.gmx_ga2la_t)                ",sizeof(gt.gmx_ga2la_t)
    print "sizeof(gt.gmx_domdec_ns_ranges_t)     ",sizeof(gt.gmx_domdec_ns_ranges_t)
    print "sizeof(gt.gmx_partdec)                ",sizeof(gt.gmx_partdec)
    print "sizeof(gt.gmx_partdec_p_t)            ",sizeof(gt.gmx_partdec_p_t)
    print "sizeof(gt.gmx_multisim_t)             ",sizeof(gt.gmx_multisim_t)
    print "sizeof(gt.gmx_domdec_t)               ",sizeof(gt.gmx_domdec_t)
    print "sizeof(gt.gmx_nodecomm_t)             ",sizeof(gt.gmx_nodecomm_t)
    print "sizeof(gt.t_commrec)                  ",sizeof(gt.t_commrec)
    print "sizeof(gt.t_filenm)                   ",sizeof(gt.t_filenm)
    print "sizeof(gt.t_pargs)                    ",sizeof(gt.t_pargs)
    print "sizeof(gt.gmx_edx)                    ",sizeof(gt.gmx_edx)
    print "sizeof(gt.FILE_p)                     ",sizeof(gt.FILE_p)
    print "sizeof(gt.t_eigvec)                   ",sizeof(gt.t_eigvec)
    print "sizeof(gt.t_edvecs)                   ",sizeof(gt.t_edvecs)
    print "sizeof(gt.t_edflood)                  ",sizeof(gt.t_edflood)
    print "sizeof(gt.t_do_edfit)                 ",sizeof(gt.t_do_edfit)
    print "sizeof(gt.t_do_edsam)                 ",sizeof(gt.t_do_edsam)
    print "sizeof(gt.t_do_radcon)                ",sizeof(gt.t_do_radcon)
    print "sizeof(gt.t_fitit)                    ",sizeof(gt.t_fitit)
    print "sizeof(gt.t_remove_pbc_effect)        ",sizeof(gt.t_remove_pbc_effect)
    print "sizeof(gt.t_ed_buffer)                ",sizeof(gt.t_ed_buffer)
    print "sizeof(gt.edpar)                      ",sizeof(gt.edpar)
    print "sizeof(gt.t_edpar)                    ",sizeof(gt.t_edpar)
    print "sizeof(gt.gmx_edsam)                  ",sizeof(gt.gmx_edsam)
    print "sizeof(gt.gmx_edsam_t)                ",sizeof(gt.gmx_edsam_t)
    print "sizeof(gt.t_comm_vsites)              ",sizeof(gt.t_comm_vsites)
    print "sizeof(gt.gmx_vsite_t)                ",sizeof(gt.gmx_vsite_t)
    print "sizeof(gt.gmx_lincsdata)              ",sizeof(gt.gmx_lincsdata)
    print "sizeof(gt.gmx_lincsdata_t)            ",sizeof(gt.gmx_lincsdata_t)
    print "sizeof(gt.gmx_constr)                 ",sizeof(gt.gmx_constr)

    return