from ctypes import c_int

(
  ewcRUN,
  ewcSTEP,
  ewcPPDURINGPME,
  ewcDOMDEC,
  ewcVSITECONSTR,
  ewcPP_PMESENDX,
  ewcMOVEX,
  ewcNS,
  ewcFORCE,
  ewcMOVEF,
  ewcPMEMESH,
  ewcPMEMESH_SEP,
  ewcPMEWAITCOMM,
  ewcPP_PMEWAITRECVF,
  ewcVSITESPREAD,
  ewcTRAJ,
  ewcUPDATE,
  ewcCONSTR,
  ewcMoveE,
  ewcTEST,
  ewcNR
) = map(c_int,xrange(21))
#from gmx_wallcycle.h
#enum { ewcRUN, ewcSTEP, ewcPPDURINGPME, ewcDOMDEC, ewcVSITECONSTR, ewcPP_PMESENDX, ewcMOVEX, ewcNS, ewcFORCE, ewcMOVEF, ewcPMEMESH, ewcPMEMESH_SEP, ewcPMEWAITCOMM, ewcPP_PMEWAITRECVF, ewcVSITESPREAD, ewcTRAJ, ewcUPDATE, ewcCONSTR, ewcMoveE, ewcTEST, ewcNR };