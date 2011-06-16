from ctypes import c_char_p
import grompy.enums as enums
# from names.h

c_char_array = c_char_p * (enums.epbcNR.value + 1)
epbc_names   = c_char_array(c_char_p("xyz"),\
                            c_char_p("no"),\
                            c_char_p("xy"),\
                            c_char_p("screw"),\
                            c_char_p())

c_char_array = c_char_p * (enums.eelNR.value + 1)
eel_names    = c_char_array(c_char_p("Cut-off"),\
                            c_char_p("Reaction-Field"),\
                            c_char_p("Generalized-Reaction-Field"),\
                            c_char_p("PME"),\
                            c_char_p("Ewald"),\
                            c_char_p("PPPM"),\
                            c_char_p("Poisson"),\
                            c_char_p("Switch"),\
                            c_char_p("Shift"),\
                            c_char_p("User"),\
                            c_char_p("Generalized-Born"),\
                            c_char_p("Reaction-Field-nec"),\
                            c_char_p("Encad-shift"),\
                            c_char_p("PME-User"),\
                            c_char_p("PME-Switch"),\
                            c_char_p("PME-User-Switch"),\
                            c_char_p("Reaction-Field-zero"),\
                            c_char_p())
