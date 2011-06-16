from ctypes import c_int
#from enums.h


(
  epbcXYZ,
  epbcNONE,
  epbcXY,
  epbcSCREW,
  epbcNR
) = map(c_int,xrange(5))

(
  eelCUT,
  eelRF,
  eelGRF,
  eelPME,
  eelEWALD,
  eelPPPM,
  eelPOISSON,
  eelSWITCH,
  eelSHIFT,
  eelUSER,
  eelGB,
  eelRF_NEC,
  eelENCADSHIFT,
  eelPMEUSER,
  eelPMESWITCH,
  eelPMEUSERSWITCH,
  eelRF_ZERO,
  eelNR
) = map(c_int,xrange(18))

def EEL_PME(e = c_int):
    result = ((e==eelPME) or\
              (e==eelPMESWITCH) or\
              (e==eelPMEUSER) or\
              (e==eelPMEUSERSWITCH))
    return result
#define EEL_PME(e)  ((e) == eelPME || (e) == eelPMESWITCH || (e) == eelPMEUSER || (e) == eelPMEUSERSWITCH)


(
  eiMD,
  eiSteep,
  eiCG,
  eiBD,
  eiSD2,
  eiNM,
  eiLBFGS,
  eiTPI,
  eiTPIC,
  eiSD1,
  eiNR
) = map(c_int,xrange(11))

def EI_SD(e = c_int):
    return (e==eiSD1.value or e==eiSD2.value)
#define EI_SD(e) ((e) == eiSD1 || (e) == eiSD2)
def EI_DYNAMICS(e = c_int):
    return (e.value==eiMD.value or EI_SD(e) or e.value==eiBD.value)
#define EI_DYNAMICS(e) ((e) == eiMD || EI_SD(e) || (e) == eiBD)
def EI_TPI(e = c_int):
    return ((e == eiTPI) or (e==eiTPIC))
#define EI_TPI(e) ((e) == eiTPI || (e) == eiTPIC)

(
  efepNO,
  efepYES,
  efepNR
) = map(c_int,xrange(3))

(
  epullNO,
  epullUMBRELLA,
  epullCONSTRAINT,
  epullCONST_F,
  epullNR
) = map(c_int,xrange(5))