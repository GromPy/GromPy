;#
;# $Id: nb_kernel202_ia32_sse.intel_syntax.s,v 1.5 2008/04/06 19:43:57 lindahl Exp $
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as



	
.globl nb_kernel202_ia32_sse
.globl _nb_kernel202_ia32_sse
nb_kernel202_ia32_sse:	
_nb_kernel202_ia32_sse:	
.equiv          nb202_p_nri,            8
.equiv          nb202_iinr,             12
.equiv          nb202_jindex,           16
.equiv          nb202_jjnr,             20
.equiv          nb202_shift,            24
.equiv          nb202_shiftvec,         28
.equiv          nb202_fshift,           32
.equiv          nb202_gid,              36
.equiv          nb202_pos,              40
.equiv          nb202_faction,          44
.equiv          nb202_charge,           48
.equiv          nb202_p_facel,          52
.equiv          nb202_argkrf,           56
.equiv          nb202_argcrf,           60
.equiv          nb202_Vc,               64
.equiv          nb202_type,             68
.equiv          nb202_p_ntype,          72
.equiv          nb202_vdwparam,         76
.equiv          nb202_Vvdw,             80
.equiv          nb202_p_tabscale,       84
.equiv          nb202_VFtab,            88
.equiv          nb202_invsqrta,         92
.equiv          nb202_dvda,             96
.equiv          nb202_p_gbtabscale,     100
.equiv          nb202_GBtab,            104
.equiv          nb202_p_nthreads,       108
.equiv          nb202_count,            112
.equiv          nb202_mtx,              116
.equiv          nb202_outeriter,        120
.equiv          nb202_inneriter,        124
.equiv          nb202_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb202_ixO,              0
.equiv          nb202_iyO,              16
.equiv          nb202_izO,              32
.equiv          nb202_ixH1,             48
.equiv          nb202_iyH1,             64
.equiv          nb202_izH1,             80
.equiv          nb202_ixH2,             96
.equiv          nb202_iyH2,             112
.equiv          nb202_izH2,             128
.equiv          nb202_jxO,              144
.equiv          nb202_jyO,              160
.equiv          nb202_jzO,              176
.equiv          nb202_jxH1,             192
.equiv          nb202_jyH1,             208
.equiv          nb202_jzH1,             224
.equiv          nb202_jxH2,             240
.equiv          nb202_jyH2,             256
.equiv          nb202_jzH2,             272
.equiv          nb202_dxOO,             288
.equiv          nb202_dyOO,             304
.equiv          nb202_dzOO,             320
.equiv          nb202_dxOH1,            336
.equiv          nb202_dyOH1,            352
.equiv          nb202_dzOH1,            368
.equiv          nb202_dxOH2,            384
.equiv          nb202_dyOH2,            400
.equiv          nb202_dzOH2,            416
.equiv          nb202_dxH1O,            432
.equiv          nb202_dyH1O,            448
.equiv          nb202_dzH1O,            464
.equiv          nb202_dxH1H1,           480
.equiv          nb202_dyH1H1,           496
.equiv          nb202_dzH1H1,           512
.equiv          nb202_dxH1H2,           528
.equiv          nb202_dyH1H2,           544
.equiv          nb202_dzH1H2,           560
.equiv          nb202_dxH2O,            576
.equiv          nb202_dyH2O,            592
.equiv          nb202_dzH2O,            608
.equiv          nb202_dxH2H1,           624
.equiv          nb202_dyH2H1,           640
.equiv          nb202_dzH2H1,           656
.equiv          nb202_dxH2H2,           672
.equiv          nb202_dyH2H2,           688
.equiv          nb202_dzH2H2,           704
.equiv          nb202_qqOO,             720
.equiv          nb202_qqOH,             736
.equiv          nb202_qqHH,             752
.equiv          nb202_vctot,            768
.equiv          nb202_fixO,             784
.equiv          nb202_fiyO,             800
.equiv          nb202_fizO,             816
.equiv          nb202_fixH1,            832
.equiv          nb202_fiyH1,            848
.equiv          nb202_fizH1,            864
.equiv          nb202_fixH2,            880
.equiv          nb202_fiyH2,            896
.equiv          nb202_fizH2,            912
.equiv          nb202_fjxO,             928
.equiv          nb202_fjyO,             944
.equiv          nb202_fjzO,             960
.equiv          nb202_fjxH1,            976
.equiv          nb202_fjyH1,            992
.equiv          nb202_fjzH1,            1008
.equiv          nb202_fjxH2,            1024
.equiv          nb202_fjyH2,            1040
.equiv          nb202_fjzH2,            1056
.equiv          nb202_fjzH2b,           1060
.equiv          nb202_fjzH2c,           1064
.equiv          nb202_fjzH2d,           1068
.equiv          nb202_half,             1072
.equiv          nb202_three,            1088
.equiv          nb202_rsqOO,            1104
.equiv          nb202_rsqOH1,           1120
.equiv          nb202_rsqOH2,           1136
.equiv          nb202_rsqH1O,           1152
.equiv          nb202_rsqH1H1,          1168
.equiv          nb202_rsqH1H2,          1184
.equiv          nb202_rsqH2O,           1200
.equiv          nb202_rsqH2H1,          1216
.equiv          nb202_rsqH2H2,          1232
.equiv          nb202_rinvOO,           1248
.equiv          nb202_rinvOH1,          1264
.equiv          nb202_rinvOH2,          1280
.equiv          nb202_rinvH1O,          1296
.equiv          nb202_rinvH1H1,         1312
.equiv          nb202_rinvH1H2,         1328
.equiv          nb202_rinvH2O,          1344
.equiv          nb202_rinvH2H1,         1360
.equiv          nb202_rinvH2H2,         1376
.equiv          nb202_two,              1392
.equiv          nb202_krf,              1408
.equiv          nb202_crf,              1424
.equiv          nb202_is3,              1440
.equiv          nb202_ii3,              1444
.equiv          nb202_innerjjnr,        1448
.equiv          nb202_innerk,           1452
.equiv          nb202_n,                1456
.equiv          nb202_nn1,              1460
.equiv          nb202_nri,              1464
.equiv          nb202_nouter,           1468
.equiv          nb202_ninner,           1472
.equiv          nb202_salign,           1476
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1480		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb202_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb202_p_nri]
	mov ecx, [ecx]
	mov [esp + nb202_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb202_nouter], eax
	mov [esp + nb202_ninner], eax


	mov esi, [ebp + nb202_argkrf]
	mov edi, [ebp + nb202_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb202_krf], xmm5
	movaps [esp + nb202_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb202_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb202_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb202_p_facel]
	movss xmm6, [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb202_qqOO], xmm3
	movaps [esp + nb202_qqOH], xmm4
	movaps [esp + nb202_qqHH], xmm5

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb202_half], eax
	movss xmm1, [esp + nb202_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb202_half],  xmm1
	movaps [esp + nb202_two],  xmm2
	movaps [esp + nb202_three],  xmm3

.nb202_threadloop:
        mov   esi, [ebp + nb202_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb202_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb202_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb202_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb202_n], eax
        mov [esp + nb202_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb202_outerstart
        jmp .nb202_end
	
.nb202_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb202_nouter]
	mov [esp + nb202_nouter], ebx

.nb202_outer:
	mov   eax, [ebp + nb202_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb202_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb202_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb202_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb202_pos]    ;# eax = base of pos[]  
	mov   [esp + nb202_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb202_ixO], xmm3
	movaps [esp + nb202_iyO], xmm4
	movaps [esp + nb202_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb202_ixH1], xmm0
	movaps [esp + nb202_iyH1], xmm1
	movaps [esp + nb202_izH1], xmm2
	movaps [esp + nb202_ixH2], xmm3
	movaps [esp + nb202_iyH2], xmm4
	movaps [esp + nb202_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb202_vctot], xmm4
	movaps [esp + nb202_fixO], xmm4
	movaps [esp + nb202_fiyO], xmm4
	movaps [esp + nb202_fizO], xmm4
	movaps [esp + nb202_fixH1], xmm4
	movaps [esp + nb202_fiyH1], xmm4
	movaps [esp + nb202_fizH1], xmm4
	movaps [esp + nb202_fixH2], xmm4
	movaps [esp + nb202_fiyH2], xmm4
	movaps [esp + nb202_fizH2], xmm4
	
	mov   eax, [ebp + nb202_jindex]
	mov   ecx, [eax + esi*4]     		;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]		;# jindex[n+1] 
	sub   edx, ecx               		;# number of innerloop atoms 

	mov   esi, [ebp + nb202_pos]
	mov   edi, [ebp + nb202_faction]	
	mov   eax, [ebp + nb202_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb202_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb202_ninner]
	mov   [esp + nb202_ninner], ecx
	add   edx, 0
	mov   [esp + nb202_innerk], edx    ;# number of innerloop atoms 
	jge   .nb202_unroll_loop
	jmp   .nb202_single_check
.nb202_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb202_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb202_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb202_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + nb202_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [esp + nb202_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + nb202_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + nb202_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + nb202_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + nb202_jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + nb202_jzO],  xmm0
	movaps [esp + nb202_jzH1],  xmm1
	movaps [esp + nb202_jzH2],  xmm2

	movaps xmm0, [esp + nb202_ixO]
	movaps xmm1, [esp + nb202_iyO]
	movaps xmm2, [esp + nb202_izO]
	movaps xmm3, [esp + nb202_ixO]
	movaps xmm4, [esp + nb202_iyO]
	movaps xmm5, [esp + nb202_izO]
	subps  xmm0, [esp + nb202_jxO]
	subps  xmm1, [esp + nb202_jyO]
	subps  xmm2, [esp + nb202_jzO]
	subps  xmm3, [esp + nb202_jxH1]
	subps  xmm4, [esp + nb202_jyH1]
	subps  xmm5, [esp + nb202_jzH1]
	movaps [esp + nb202_dxOO], xmm0
	movaps [esp + nb202_dyOO], xmm1
	movaps [esp + nb202_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb202_dxOH1], xmm3
	movaps [esp + nb202_dyOH1], xmm4
	movaps [esp + nb202_dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb202_rsqOO], xmm0
	movaps [esp + nb202_rsqOH1], xmm3

	movaps xmm0, [esp + nb202_ixO]
	movaps xmm1, [esp + nb202_iyO]
	movaps xmm2, [esp + nb202_izO]
	movaps xmm3, [esp + nb202_ixH1]
	movaps xmm4, [esp + nb202_iyH1]
	movaps xmm5, [esp + nb202_izH1]
	subps  xmm0, [esp + nb202_jxH2]
	subps  xmm1, [esp + nb202_jyH2]
	subps  xmm2, [esp + nb202_jzH2]
	subps  xmm3, [esp + nb202_jxO]
	subps  xmm4, [esp + nb202_jyO]
	subps  xmm5, [esp + nb202_jzO]
	movaps [esp + nb202_dxOH2], xmm0
	movaps [esp + nb202_dyOH2], xmm1
	movaps [esp + nb202_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb202_dxH1O], xmm3
	movaps [esp + nb202_dyH1O], xmm4
	movaps [esp + nb202_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb202_rsqOH2], xmm0
	movaps [esp + nb202_rsqH1O], xmm3

	movaps xmm0, [esp + nb202_ixH1]
	movaps xmm1, [esp + nb202_iyH1]
	movaps xmm2, [esp + nb202_izH1]
	movaps xmm3, [esp + nb202_ixH1]
	movaps xmm4, [esp + nb202_iyH1]
	movaps xmm5, [esp + nb202_izH1]
	subps  xmm0, [esp + nb202_jxH1]
	subps  xmm1, [esp + nb202_jyH1]
	subps  xmm2, [esp + nb202_jzH1]
	subps  xmm3, [esp + nb202_jxH2]
	subps  xmm4, [esp + nb202_jyH2]
	subps  xmm5, [esp + nb202_jzH2]
	movaps [esp + nb202_dxH1H1], xmm0
	movaps [esp + nb202_dyH1H1], xmm1
	movaps [esp + nb202_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb202_dxH1H2], xmm3
	movaps [esp + nb202_dyH1H2], xmm4
	movaps [esp + nb202_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb202_rsqH1H1], xmm0
	movaps [esp + nb202_rsqH1H2], xmm3

	movaps xmm0, [esp + nb202_ixH2]
	movaps xmm1, [esp + nb202_iyH2]
	movaps xmm2, [esp + nb202_izH2]
	movaps xmm3, [esp + nb202_ixH2]
	movaps xmm4, [esp + nb202_iyH2]
	movaps xmm5, [esp + nb202_izH2]
	subps  xmm0, [esp + nb202_jxO]
	subps  xmm1, [esp + nb202_jyO]
	subps  xmm2, [esp + nb202_jzO]
	subps  xmm3, [esp + nb202_jxH1]
	subps  xmm4, [esp + nb202_jyH1]
	subps  xmm5, [esp + nb202_jzH1]
	movaps [esp + nb202_dxH2O], xmm0
	movaps [esp + nb202_dyH2O], xmm1
	movaps [esp + nb202_dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + nb202_dxH2H1], xmm3
	movaps [esp + nb202_dyH2H1], xmm4
	movaps [esp + nb202_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb202_rsqH2O], xmm0
	movaps [esp + nb202_rsqH2H1], xmm4

	movaps xmm0, [esp + nb202_ixH2]
	movaps xmm1, [esp + nb202_iyH2]
	movaps xmm2, [esp + nb202_izH2]
	subps  xmm0, [esp + nb202_jxH2]
	subps  xmm1, [esp + nb202_jyH2]
	subps  xmm2, [esp + nb202_jzH2]
	movaps [esp + nb202_dxH2H2], xmm0
	movaps [esp + nb202_dyH2H2], xmm1
	movaps [esp + nb202_dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + nb202_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202_half] ;# rinvH2H2 
	mulps   xmm7, [esp + nb202_half] ;# rinvH2H1 
	movaps  [esp + nb202_rinvH2H2], xmm3
	movaps  [esp + nb202_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + nb202_rsqOO]
	rsqrtps xmm5, [esp + nb202_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb202_rsqOO]
	mulps   xmm5, [esp + nb202_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202_half] 
	mulps   xmm7, [esp + nb202_half]
	movaps  [esp + nb202_rinvOO], xmm3
	movaps  [esp + nb202_rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + nb202_rsqOH2]
	rsqrtps xmm5, [esp + nb202_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb202_rsqOH2]
	mulps   xmm5, [esp + nb202_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202_half] 
	mulps   xmm7, [esp + nb202_half]
	movaps  [esp + nb202_rinvOH2], xmm3
	movaps  [esp + nb202_rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + nb202_rsqH1H1]
	rsqrtps xmm5, [esp + nb202_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb202_rsqH1H1]
	mulps   xmm5, [esp + nb202_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202_half] 
	mulps   xmm7, [esp + nb202_half]
	movaps  [esp + nb202_rinvH1H1], xmm3
	movaps  [esp + nb202_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + nb202_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb202_three]
	mulps   xmm1, [esp + nb202_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb202_half] 
	movaps  [esp + nb202_rinvH2O], xmm3

	;# start with OO interaction 
	movaps xmm0, [esp + nb202_rinvOO]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]
	mulps  xmm0, xmm0	;# xmm0=rinvsq 

	mulps  xmm5, [esp + nb202_rsqOO] ;# xmm5=krsq 
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb202_crf]
	mulps  xmm6, [esp + nb202_qqOO] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	mulps xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqOO] ;# xmm7 = coul part of fscal 
	
	addps  xmm6, [esp + nb202_vctot] ;# local vctot summation variable 
	mulps  xmm0, xmm7
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb202_dxOO]
	mulps xmm1, [esp + nb202_dyOO]
	mulps xmm2, [esp + nb202_dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixO]
	addps xmm1, [esp + nb202_fiyO]
	addps xmm2, [esp + nb202_fizO]
	movaps [esp + nb202_fjxO], xmm3
	movaps [esp + nb202_fjyO], xmm4
	movaps [esp + nb202_fjzO], xmm5
	movaps [esp + nb202_fixO], xmm0
	movaps [esp + nb202_fiyO], xmm1
	movaps [esp + nb202_fizO], xmm2

	;# O-H1 interaction 
	movaps xmm0, [esp + nb202_rinvOH1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqOH1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps  xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqOH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqOH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH1  
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb202_dxOH1]
	mulps xmm1, [esp + nb202_dyOH1]
	mulps xmm2, [esp + nb202_dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixO]
	addps xmm1, [esp + nb202_fiyO]
	addps xmm2, [esp + nb202_fizO]
	movaps [esp + nb202_fjxH1], xmm3
	movaps [esp + nb202_fjyH1], xmm4
	movaps [esp + nb202_fjzH1], xmm5
	movaps [esp + nb202_fixO], xmm0
	movaps [esp + nb202_fiyO], xmm1
	movaps [esp + nb202_fizO], xmm2

	;# O-H2 interaction  
	movaps xmm0, [esp + nb202_rinvOH2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqOH2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqOH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqOH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + nb202_dxOH2]
	mulps xmm1, [esp + nb202_dyOH2]
	mulps xmm2, [esp + nb202_dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixO]
	addps xmm1, [esp + nb202_fiyO]
	addps xmm2, [esp + nb202_fizO]
	movaps [esp + nb202_fjxH2], xmm3
	movaps [esp + nb202_fjyH2], xmm4
	movaps [esp + nb202_fjzH2], xmm5
	movaps [esp + nb202_fixO], xmm0
	movaps [esp + nb202_fiyO], xmm1
	movaps [esp + nb202_fizO], xmm2

	;# H1-O interaction 
	movaps xmm0, [esp + nb202_rinvH1O]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqH1O] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqOH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqOH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb202_fjxO]
	movaps xmm4, [esp + nb202_fjyO]
	movaps xmm5, [esp + nb202_fjzO]
	mulps xmm0, [esp + nb202_dxH1O]
	mulps xmm1, [esp + nb202_dyH1O]
	mulps xmm2, [esp + nb202_dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixH1]
	addps xmm1, [esp + nb202_fiyH1]
	addps xmm2, [esp + nb202_fizH1]
	movaps [esp + nb202_fjxO], xmm3
	movaps [esp + nb202_fjyO], xmm4
	movaps [esp + nb202_fjzO], xmm5
	movaps [esp + nb202_fixH1], xmm0
	movaps [esp + nb202_fiyH1], xmm1
	movaps [esp + nb202_fizH1], xmm2

	;# H1-H1 interaction 
	movaps xmm0, [esp + nb202_rinvH1H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqH1H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb202_fjxH1]
	movaps xmm4, [esp + nb202_fjyH1]
	movaps xmm5, [esp + nb202_fjzH1]
	mulps xmm0, [esp + nb202_dxH1H1]
	mulps xmm1, [esp + nb202_dyH1H1]
	mulps xmm2, [esp + nb202_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixH1]
	addps xmm1, [esp + nb202_fiyH1]
	addps xmm2, [esp + nb202_fizH1]
	movaps [esp + nb202_fjxH1], xmm3
	movaps [esp + nb202_fjyH1], xmm4
	movaps [esp + nb202_fjzH1], xmm5
	movaps [esp + nb202_fixH1], xmm0
	movaps [esp + nb202_fiyH1], xmm1
	movaps [esp + nb202_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [esp + nb202_rinvH1H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqH1H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb202_fjxH2]
	movaps xmm4, [esp + nb202_fjyH2]
	movaps xmm5, [esp + nb202_fjzH2]
	mulps xmm0, [esp + nb202_dxH1H2]
	mulps xmm1, [esp + nb202_dyH1H2]
	mulps xmm2, [esp + nb202_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixH1]
	addps xmm1, [esp + nb202_fiyH1]
	addps xmm2, [esp + nb202_fizH1]
	movaps [esp + nb202_fjxH2], xmm3
	movaps [esp + nb202_fjyH2], xmm4
	movaps [esp + nb202_fjzH2], xmm5
	movaps [esp + nb202_fixH1], xmm0
	movaps [esp + nb202_fiyH1], xmm1
	movaps [esp + nb202_fizH1], xmm2

	;# H2-O interaction 
	movaps xmm0, [esp + nb202_rinvH2O]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqH2O] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqOH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqOH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb202_fjxO]
	movaps xmm4, [esp + nb202_fjyO]
	movaps xmm5, [esp + nb202_fjzO]
	mulps xmm0, [esp + nb202_dxH2O]
	mulps xmm1, [esp + nb202_dyH2O]
	mulps xmm2, [esp + nb202_dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixH2]
	addps xmm1, [esp + nb202_fiyH2]
	addps xmm2, [esp + nb202_fizH2]
	movaps [esp + nb202_fjxO], xmm3
	movaps [esp + nb202_fjyO], xmm4
	movaps [esp + nb202_fjzO], xmm5
	movaps [esp + nb202_fixH2], xmm0
	movaps [esp + nb202_fiyH2], xmm1
	movaps [esp + nb202_fizH2], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [esp + nb202_rinvH2H1]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqH2H1] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + nb202_fjxH1]
	movaps xmm4, [esp + nb202_fjyH1]
	movaps xmm5, [esp + nb202_fjzH1]
	mulps xmm0, [esp + nb202_dxH2H1]
	mulps xmm1, [esp + nb202_dyH2H1]
	mulps xmm2, [esp + nb202_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixH2]
	addps xmm1, [esp + nb202_fiyH2]
	addps xmm2, [esp + nb202_fizH2]
	movaps [esp + nb202_fjxH1], xmm3
	movaps [esp + nb202_fjyH1], xmm4
	movaps [esp + nb202_fjzH1], xmm5
	movaps [esp + nb202_fixH2], xmm0
	movaps [esp + nb202_fiyH2], xmm1
	movaps [esp + nb202_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [esp + nb202_rinvH2H2]
	movaps xmm7, xmm0	;# xmm7=rinv 
	movaps xmm5, [esp + nb202_krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + nb202_rsqH2H2] ;# xmm5=krsq 
	movaps xmm4, xmm5
	addps  xmm4, xmm7	;# xmm4=r inv+ krsq 
	subps  xmm4, [esp + nb202_crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + nb202_qqHH] ;# xmm4=voul=qq*(rinv+ krsq-crf) 
	mulps  xmm5, [esp + nb202_two]
	subps  xmm7, xmm5	;# xmm7=rinv-2*krsq 
	mulps  xmm7, [esp + nb202_qqHH] ;# xmm7 = coul part of fscal 
	addps  xmm6, xmm4	;# add to local vctot 
	mulps xmm0, xmm7	;# fsOH2 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm1, xmm0
	movaps [esp + nb202_vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + nb202_fjxH2]
	movaps xmm4, [esp + nb202_fjyH2]
	movaps xmm5, [esp + nb202_fjzH2]
	mulps xmm0, [esp + nb202_dxH2H2]
	mulps xmm1, [esp + nb202_dyH2H2]
	mulps xmm2, [esp + nb202_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + nb202_fixH2]
	addps xmm1, [esp + nb202_fiyH2]
	addps xmm2, [esp + nb202_fizH2]
	movaps [esp + nb202_fjxH2], xmm3
	movaps [esp + nb202_fjyH2], xmm4
	movaps [esp + nb202_fjzH2], xmm5
	movaps [esp + nb202_fixH2], xmm0
	movaps [esp + nb202_fiyH2], xmm1
	movaps [esp + nb202_fizH2], xmm2

	mov edi, [ebp + nb202_faction]
		
	;# Did all interactions - now update j forces 
	;# At this stage forces are still on the stack, in positions:
	;# fjxO, fjyO, fjzO, ... , fjzH2.
	;# Each position is a quadruplet of forces for the four 
	;# corresponding j waters, so we need to transpose them before
	;# adding to the memory positions.
	;# 
	;# This _used_ to be a simple transpose, but the resulting high number
	;# of unaligned 128-bit load/stores might trigger a possible hardware 
	;# bug on Athlon and Opteron chips, so I have worked around it
	;# to use 64-bit load/stores instead. The performance hit should be
	;# very modest, since the 128-bit unaligned memory instructions were
	;# slow anyway. 
	
		
	;# 4 j waters with three atoms each - first do Oxygen X & Y forces for 4 j particles 
	movaps xmm0, [esp + nb202_fjxO] ;# xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm2, [esp + nb202_fjyO] ;# xmm1= fjyOa  fjyOb  fjyOc  fjyOd
	movlps xmm3, [edi + eax*4]
	movlps xmm4, [edi + ecx*4]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	unpckhps xmm1, xmm2        ;# xmm1= fjxOc  fjyOc  fjxOd  fjyOd
	movhps xmm3, [edi + ebx*4]
	movhps xmm4, [edi + edx*4]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4], xmm3
	movlps [edi + ecx*4], xmm4
	movhps [edi + ebx*4], xmm3
	movhps [edi + edx*4], xmm4
	
	;# Oxygen Z & first hydrogen X forces for 4 j particles 
	movaps xmm0, [esp + nb202_fjzO]  ;# xmm0= fjzOa   fjzOb   fjzOc   fjzOd 
	movaps xmm2, [esp + nb202_fjxH1] ;# xmm1= fjxH1a  fjxH1b  fjxH1c  fjxH1d
	movlps xmm3, [edi + eax*4 + 8]
	movlps xmm4, [edi + ecx*4 + 8]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjzOa  fjxH1a  fjzOb  fjxH1b
	unpckhps xmm1, xmm2        ;# xmm1= fjzOc  fjxH1c  fjzOd  fjxH1d
	movhps xmm3, [edi + ebx*4 + 8]
	movhps xmm4, [edi + edx*4 + 8]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 8], xmm3
	movlps [edi + ecx*4 + 8], xmm4
	movhps [edi + ebx*4 + 8], xmm3
	movhps [edi + edx*4 + 8], xmm4

	
	;# First hydrogen Y & Z forces for 4 j particles 
	movaps xmm0, [esp + nb202_fjyH1]  ;# xmm0= fjyH1a  fjyH1b  fjyH1c  fjyH1d 
	movaps xmm2, [esp + nb202_fjzH1] ;# xmm1= fjzH1a  fjzH1b  fjzH1c  fjzH1d
	movlps xmm3, [edi + eax*4 + 16]
	movlps xmm4, [edi + ecx*4 + 16]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjyH1a  fjzH1a  fjyH1b  fjzH1b
	unpckhps xmm1, xmm2		;# xmm1= fjyH1c  fjzH1c  fjyH1d  fjzH1d
	movhps xmm3, [edi + ebx*4 + 16]
	movhps xmm4, [edi + edx*4 + 16]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 16], xmm3
	movlps [edi + ecx*4 + 16], xmm4
	movhps [edi + ebx*4 + 16], xmm3
	movhps [edi + edx*4 + 16], xmm4

	
	;# Second hydrogen X & Y forces for 4 j particles 
	movaps xmm0, [esp + nb202_fjxH2]  ;# xmm0= fjxH2a  fjxH2b  fjxH2c  fjxH2d 
	movaps xmm2, [esp + nb202_fjyH2] ;# xmm1= fjyH2a  fjyH2b  fjyH2c  fjyH2d
	movlps xmm3, [edi + eax*4 + 24]
	movlps xmm4, [edi + ecx*4 + 24]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjxH2a  fjyH2a  fjxH2b  fjyH2b
	unpckhps xmm1, xmm2		;# xmm1= fjxH2c  fjyH2c  fjxH2d  fjyH2d
	movhps xmm3, [edi + ebx*4 + 24]
	movhps xmm4, [edi + edx*4 + 24]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps [edi + eax*4 + 24], xmm3
	movlps [edi + ecx*4 + 24], xmm4
	movhps [edi + ebx*4 + 24], xmm3
	movhps [edi + edx*4 + 24], xmm4

	
	;# Second hydrogen Z forces for 4 j particles 
	;# Just load the four Z coords into one reg. each
	movss xmm4, [edi + eax*4 + 32]
	movss xmm5, [edi + ebx*4 + 32]
	movss xmm6, [edi + ecx*4 + 32]
	movss xmm7, [edi + edx*4 + 32]
	;# add what we have on the stack
	addss xmm4, [esp + nb202_fjzH2] 
	addss xmm5, [esp + nb202_fjzH2b] 
	addss xmm6, [esp + nb202_fjzH2c] 
	addss xmm7, [esp + nb202_fjzH2d]
	;# store back
	movss [edi + eax*4 + 32], xmm4
	movss [edi + ebx*4 + 32], xmm5
	movss [edi + ecx*4 + 32], xmm6
	movss [edi + edx*4 + 32], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb202_innerk],  4
	jl    .nb202_single_check
	jmp   .nb202_unroll_loop
.nb202_single_check:
	add dword ptr [esp + nb202_innerk],  4
	jnz   .nb202_single_loop
	jmp   .nb202_updateouterdata
.nb202_single_loop:
	mov   edx, [esp + nb202_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb202_innerjjnr],  4	

	mov esi, [ebp + nb202_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [esi + eax*4]		;# jxO  -  -  -
	movss xmm4, [esi + eax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [esi + eax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + nb202_ixO]     
	movaps  xmm1, [esp + nb202_iyO]
	movaps  xmm2, [esp + nb202_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzO   0   jzH1 jzH2
	
	;# store all j coordinates in jO  
	movaps [esp + nb202_jxO], xmm3
	movaps [esp + nb202_jyO], xmm4
	movaps [esp + nb202_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + nb202_dxOO], xmm0
	movaps [esp + nb202_dyOO], xmm1
	movaps [esp + nb202_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 

	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [esp + nb202_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6     ;# xmm7=krsq 
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb202_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + nb202_half] ;# rinv iO - j water 


	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	mulps   xmm7, [esp + nb202_two]
	subps   xmm6, [esp + nb202_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	subps   xmm3, xmm7	;# xmm3=rinv-2*krsq 
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [esp + nb202_qqOO]
	movhps  xmm4, [esp + nb202_qqOH]

	mulps xmm6, xmm4	;# vcoul  
	mulps xmm3, xmm4	;# coul part of fs  


	addps   xmm6, [esp + nb202_vctot]
	mulps   xmm0, xmm3	;# total fscal 
	movaps  [esp + nb202_vctot], xmm6	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + nb202_dxOO]
	mulps   xmm1, [esp + nb202_dyOO]
	mulps   xmm2, [esp + nb202_dzOO]
	
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb202_fjxO], xmm3
	movaps  [esp + nb202_fjyO], xmm4
	movaps  [esp + nb202_fjzO], xmm5
	addps   xmm0, [esp + nb202_fixO]
	addps   xmm1, [esp + nb202_fiyO]
	addps   xmm2, [esp + nb202_fizO]
	movaps  [esp + nb202_fixO], xmm0
	movaps  [esp + nb202_fiyO], xmm1
	movaps  [esp + nb202_fizO], xmm2

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb202_ixH1]
	movaps  xmm1, [esp + nb202_iyH1]
	movaps  xmm2, [esp + nb202_izH1]	
	movaps  xmm3, [esp + nb202_ixH2] 
	movaps  xmm4, [esp + nb202_iyH2] 
	movaps  xmm5, [esp + nb202_izH2] 
	subps   xmm0, [esp + nb202_jxO]
	subps   xmm1, [esp + nb202_jyO]
	subps   xmm2, [esp + nb202_jzO]
	subps   xmm3, [esp + nb202_jxO]
	subps   xmm4, [esp + nb202_jyO]
	subps   xmm5, [esp + nb202_jzO]
	movaps [esp + nb202_dxH1O], xmm0
	movaps [esp + nb202_dyH1O], xmm1
	movaps [esp + nb202_dzH1O], xmm2
	movaps [esp + nb202_dxH2O], xmm3
	movaps [esp + nb202_dyH2O], xmm4
	movaps [esp + nb202_dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + nb202_half] ;# rinv H2 - j water  

	mulps xmm0, [esp + nb202_krf] ;# krsq 
	mulps xmm4, [esp + nb202_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [esp + nb202_qqOH]
	movhps  xmm6, [esp + nb202_qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [esp + nb202_crf]
	subps   xmm4, [esp + nb202_crf]
	mulps   xmm1, [esp + nb202_two]
	mulps   xmm5, [esp + nb202_two]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [esp + nb202_vctot]
	movaps  [esp + nb202_vctot], xmm4
	movaps  xmm0, xmm3
	movaps  xmm4, xmm7
	mulps   xmm3, xmm3
	mulps   xmm7, xmm7
	subps   xmm0, xmm1
	subps   xmm4, xmm5
	mulps   xmm0, xmm6
	mulps   xmm4, xmm6
	mulps   xmm0, xmm3	;# fscal 
	mulps   xmm7, xmm4	;# fscal 
	
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + nb202_dxH1O]
	mulps   xmm1, [esp + nb202_dyH1O]
	mulps   xmm2, [esp + nb202_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + nb202_fjxO]
	movaps  xmm4, [esp + nb202_fjyO]
	movaps  xmm5, [esp + nb202_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + nb202_fjxO], xmm3
	movaps  [esp + nb202_fjyO], xmm4
	movaps  [esp + nb202_fjzO], xmm5
	addps   xmm0, [esp + nb202_fixH1]
	addps   xmm1, [esp + nb202_fiyH1]
	addps   xmm2, [esp + nb202_fizH1]
	movaps  [esp + nb202_fixH1], xmm0
	movaps  [esp + nb202_fiyH1], xmm1
	movaps  [esp + nb202_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + nb202_dxH2O]
	mulps   xmm1, [esp + nb202_dyH2O]
	mulps   xmm2, [esp + nb202_dzH2O]
	movaps  xmm3, [esp + nb202_fjxO]
	movaps  xmm4, [esp + nb202_fjyO]
	movaps  xmm5, [esp + nb202_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + nb202_faction]
	movaps  [esp + nb202_fjxO], xmm3
	movaps  [esp + nb202_fjyO], xmm4
	movaps  [esp + nb202_fjzO], xmm5
	addps   xmm0, [esp + nb202_fixH2]
	addps   xmm1, [esp + nb202_fiyH2]
	addps   xmm2, [esp + nb202_fizH2]
	movaps  [esp + nb202_fixH2], xmm0
	movaps  [esp + nb202_fiyH2], xmm1
	movaps  [esp + nb202_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + nb202_fjxO]
	movaps  xmm4, [esp + nb202_fjyO]
	movaps  xmm5, [esp + nb202_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# constant 00000010
	shufps  xmm7, xmm7, 3 ;# constant 00000011
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec dword ptr [esp + nb202_innerk]
	jz    .nb202_updateouterdata
	jmp   .nb202_single_loop
.nb202_updateouterdata:
	mov   ecx, [esp + nb202_ii3]
	mov   edi, [ebp + nb202_faction]
	mov   esi, [ebp + nb202_fshift]
	mov   edx, [esp + nb202_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb202_fixO]
	movaps xmm1, [esp + nb202_fiyO] 
	movaps xmm2, [esp + nb202_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	
 
	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb202_fixH1]
	movaps xmm1, [esp + nb202_fiyH1]
	movaps xmm2, [esp + nb202_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb202_fixH2]
	movaps xmm1, [esp + nb202_fiyH2]
	movaps xmm2, [esp + nb202_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb202_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb202_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb202_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb202_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb202_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb202_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb202_n], esi
        jmp .nb202_outer
.nb202_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb202_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb202_end
        ;# non-zero, do one more workunit
        jmp   .nb202_threadloop
.nb202_end:
	emms

	mov eax, [esp + nb202_nouter]
	mov ebx, [esp + nb202_ninner]
	mov ecx, [ebp + nb202_outeriter]
	mov edx, [ebp + nb202_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb202_salign]
	add esp, eax
	add esp, 1480
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret







	
.globl nb_kernel202nf_ia32_sse
.globl _nb_kernel202nf_ia32_sse
nb_kernel202nf_ia32_sse:	
_nb_kernel202nf_ia32_sse:	
.equiv          nb202nf_p_nri,          8
.equiv          nb202nf_iinr,           12
.equiv          nb202nf_jindex,         16
.equiv          nb202nf_jjnr,           20
.equiv          nb202nf_shift,          24
.equiv          nb202nf_shiftvec,       28
.equiv          nb202nf_fshift,         32
.equiv          nb202nf_gid,            36
.equiv          nb202nf_pos,            40
.equiv          nb202nf_faction,        44
.equiv          nb202nf_charge,         48
.equiv          nb202nf_p_facel,        52
.equiv          nb202nf_argkrf,         56
.equiv          nb202nf_argcrf,         60
.equiv          nb202nf_Vc,             64
.equiv          nb202nf_type,           68
.equiv          nb202nf_p_ntype,        72
.equiv          nb202nf_vdwparam,       76
.equiv          nb202nf_Vvdw,           80
.equiv          nb202nf_p_tabscale,     84
.equiv          nb202nf_VFtab,          88
.equiv          nb202nf_invsqrta,       92
.equiv          nb202nf_dvda,           96
.equiv          nb202nf_p_gbtabscale,   100
.equiv          nb202nf_GBtab,          104
.equiv          nb202nf_p_nthreads,     108
.equiv          nb202nf_count,          112
.equiv          nb202nf_mtx,            116
.equiv          nb202nf_outeriter,      120
.equiv          nb202nf_inneriter,      124
.equiv          nb202nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb202nf_ixO,            0
.equiv          nb202nf_iyO,            16
.equiv          nb202nf_izO,            32
.equiv          nb202nf_ixH1,           48
.equiv          nb202nf_iyH1,           64
.equiv          nb202nf_izH1,           80
.equiv          nb202nf_ixH2,           96
.equiv          nb202nf_iyH2,           112
.equiv          nb202nf_izH2,           128
.equiv          nb202nf_jxO,            144
.equiv          nb202nf_jyO,            160
.equiv          nb202nf_jzO,            176
.equiv          nb202nf_jxH1,           192
.equiv          nb202nf_jyH1,           208
.equiv          nb202nf_jzH1,           224
.equiv          nb202nf_jxH2,           240
.equiv          nb202nf_jyH2,           256
.equiv          nb202nf_jzH2,           272
.equiv          nb202nf_qqOO,           288
.equiv          nb202nf_qqOH,           304
.equiv          nb202nf_qqHH,           320
.equiv          nb202nf_vctot,          336
.equiv          nb202nf_half,           352
.equiv          nb202nf_three,          368
.equiv          nb202nf_rsqOO,          384
.equiv          nb202nf_rsqOH1,         400
.equiv          nb202nf_rsqOH2,         416
.equiv          nb202nf_rsqH1O,         432
.equiv          nb202nf_rsqH1H1,        448
.equiv          nb202nf_rsqH1H2,        464
.equiv          nb202nf_rsqH2O,         480
.equiv          nb202nf_rsqH2H1,        496
.equiv          nb202nf_rsqH2H2,        512
.equiv          nb202nf_rinvOO,         528
.equiv          nb202nf_rinvOH1,        544
.equiv          nb202nf_rinvOH2,        560
.equiv          nb202nf_rinvH1O,        576
.equiv          nb202nf_rinvH1H1,       592
.equiv          nb202nf_rinvH1H2,       608
.equiv          nb202nf_rinvH2O,        624
.equiv          nb202nf_rinvH2H1,       640
.equiv          nb202nf_rinvH2H2,       656
.equiv          nb202nf_krf,            672
.equiv          nb202nf_crf,            688
.equiv          nb202nf_is3,            704
.equiv          nb202nf_ii3,            708
.equiv          nb202nf_innerjjnr,      712
.equiv          nb202nf_innerk,         716
.equiv          nb202nf_n,              720
.equiv          nb202nf_nn1,            724
.equiv          nb202nf_nri,            728
.equiv          nb202nf_nouter,         732
.equiv          nb202nf_ninner,         736
.equiv          nb202nf_salign,         740
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 744		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb202nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb202nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb202nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb202nf_nouter], eax
	mov [esp + nb202nf_ninner], eax


	mov esi, [ebp + nb202nf_argkrf]
	mov edi, [ebp + nb202nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb202nf_krf], xmm5
	movaps [esp + nb202nf_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb202nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb202nf_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb202nf_p_facel]
	movss xmm6, [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb202nf_qqOO], xmm3
	movaps [esp + nb202nf_qqOH], xmm4
	movaps [esp + nb202nf_qqHH], xmm5
	

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb202nf_half], eax
	movss xmm1, [esp + nb202nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb202nf_half],  xmm1
	movaps [esp + nb202nf_three],  xmm3

.nb202nf_threadloop:
        mov   esi, [ebp + nb202nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb202nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb202nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb202nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb202nf_n], eax
        mov [esp + nb202nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb202nf_outerstart
        jmp .nb202nf_end
	
.nb202nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb202nf_nouter]
	mov [esp + nb202nf_nouter], ebx

.nb202nf_outer:
	mov   eax, [ebp + nb202nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb202nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb202nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb202nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb202nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb202nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb202nf_ixO], xmm3
	movaps [esp + nb202nf_iyO], xmm4
	movaps [esp + nb202nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb202nf_ixH1], xmm0
	movaps [esp + nb202nf_iyH1], xmm1
	movaps [esp + nb202nf_izH1], xmm2
	movaps [esp + nb202nf_ixH2], xmm3
	movaps [esp + nb202nf_iyH2], xmm4
	movaps [esp + nb202nf_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb202nf_vctot], xmm4
	
	mov   eax, [ebp + nb202nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb202nf_pos]	
	mov   eax, [ebp + nb202nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb202nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb202nf_ninner]
	mov   [esp + nb202nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb202nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb202nf_unroll_loop
	jmp   .nb202nf_single_check
.nb202nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb202nf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	
	add dword ptr [esp + nb202nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb202nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + nb202nf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [esp + nb202nf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + nb202nf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + nb202nf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + nb202nf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + nb202nf_jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + nb202nf_jzO],  xmm0
	movaps [esp + nb202nf_jzH1],  xmm1
	movaps [esp + nb202nf_jzH2],  xmm2

	movaps xmm0, [esp + nb202nf_ixO]
	movaps xmm1, [esp + nb202nf_iyO]
	movaps xmm2, [esp + nb202nf_izO]
	movaps xmm3, [esp + nb202nf_ixO]
	movaps xmm4, [esp + nb202nf_iyO]
	movaps xmm5, [esp + nb202nf_izO]
	subps  xmm0, [esp + nb202nf_jxO]
	subps  xmm1, [esp + nb202nf_jyO]
	subps  xmm2, [esp + nb202nf_jzO]
	subps  xmm3, [esp + nb202nf_jxH1]
	subps  xmm4, [esp + nb202nf_jyH1]
	subps  xmm5, [esp + nb202nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb202nf_rsqOO], xmm0
	movaps [esp + nb202nf_rsqOH1], xmm3

	movaps xmm0, [esp + nb202nf_ixO]
	movaps xmm1, [esp + nb202nf_iyO]
	movaps xmm2, [esp + nb202nf_izO]
	movaps xmm3, [esp + nb202nf_ixH1]
	movaps xmm4, [esp + nb202nf_iyH1]
	movaps xmm5, [esp + nb202nf_izH1]
	subps  xmm0, [esp + nb202nf_jxH2]
	subps  xmm1, [esp + nb202nf_jyH2]
	subps  xmm2, [esp + nb202nf_jzH2]
	subps  xmm3, [esp + nb202nf_jxO]
	subps  xmm4, [esp + nb202nf_jyO]
	subps  xmm5, [esp + nb202nf_jzO]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb202nf_rsqOH2], xmm0
	movaps [esp + nb202nf_rsqH1O], xmm3

	movaps xmm0, [esp + nb202nf_ixH1]
	movaps xmm1, [esp + nb202nf_iyH1]
	movaps xmm2, [esp + nb202nf_izH1]
	movaps xmm3, [esp + nb202nf_ixH1]
	movaps xmm4, [esp + nb202nf_iyH1]
	movaps xmm5, [esp + nb202nf_izH1]
	subps  xmm0, [esp + nb202nf_jxH1]
	subps  xmm1, [esp + nb202nf_jyH1]
	subps  xmm2, [esp + nb202nf_jzH1]
	subps  xmm3, [esp + nb202nf_jxH2]
	subps  xmm4, [esp + nb202nf_jyH2]
	subps  xmm5, [esp + nb202nf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + nb202nf_rsqH1H1], xmm0
	movaps [esp + nb202nf_rsqH1H2], xmm3

	movaps xmm0, [esp + nb202nf_ixH2]
	movaps xmm1, [esp + nb202nf_iyH2]
	movaps xmm2, [esp + nb202nf_izH2]
	movaps xmm3, [esp + nb202nf_ixH2]
	movaps xmm4, [esp + nb202nf_iyH2]
	movaps xmm5, [esp + nb202nf_izH2]
	subps  xmm0, [esp + nb202nf_jxO]
	subps  xmm1, [esp + nb202nf_jyO]
	subps  xmm2, [esp + nb202nf_jzO]
	subps  xmm3, [esp + nb202nf_jxH1]
	subps  xmm4, [esp + nb202nf_jyH1]
	subps  xmm5, [esp + nb202nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + nb202nf_rsqH2O], xmm0
	movaps [esp + nb202nf_rsqH2H1], xmm4

	movaps xmm0, [esp + nb202nf_ixH2]
	movaps xmm1, [esp + nb202nf_iyH2]
	movaps xmm2, [esp + nb202nf_izH2]
	subps  xmm0, [esp + nb202nf_jxH2]
	subps  xmm1, [esp + nb202nf_jyH2]
	subps  xmm2, [esp + nb202nf_jzH2]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + nb202nf_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202nf_half] ;# rinvH2H2 
	mulps   xmm7, [esp + nb202nf_half] ;# rinvH2H1 
	movaps  [esp + nb202nf_rinvH2H2], xmm3
	movaps  [esp + nb202nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + nb202nf_rsqOO]
	rsqrtps xmm5, [esp + nb202nf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb202nf_rsqOO]
	mulps   xmm5, [esp + nb202nf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202nf_half] 
	mulps   xmm7, [esp + nb202nf_half]
	movaps  [esp + nb202nf_rinvOO], xmm3
	movaps  [esp + nb202nf_rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + nb202nf_rsqOH2]
	rsqrtps xmm5, [esp + nb202nf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb202nf_rsqOH2]
	mulps   xmm5, [esp + nb202nf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202nf_half] 
	mulps   xmm7, [esp + nb202nf_half]
	movaps  [esp + nb202nf_rinvOH2], xmm3
	movaps  [esp + nb202nf_rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + nb202nf_rsqH1H1]
	rsqrtps xmm5, [esp + nb202nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + nb202nf_rsqH1H1]
	mulps   xmm5, [esp + nb202nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202nf_half] 
	mulps   xmm7, [esp + nb202nf_half]
	movaps  [esp + nb202nf_rinvH1H1], xmm3
	movaps  [esp + nb202nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + nb202nf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb202nf_three]
	mulps   xmm1, [esp + nb202nf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + nb202nf_half] 
	movaps  [esp + nb202nf_rinvH2O], xmm3

	;# start with OO interaction 
	movaps xmm6, [esp + nb202nf_krf]
	mulps  xmm6, [esp + nb202nf_rsqOO] ;# xmm5=krsq 
	addps  xmm6, [esp + nb202nf_rinvOO]	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb202nf_crf]
	mulps  xmm6, [esp + nb202nf_qqOO] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [esp + nb202nf_vctot] ;# local vctot summation variable 

	;# O-H interactions 
	movaps xmm0, [esp + nb202nf_krf]
	movaps xmm1, [esp + nb202nf_krf]
	movaps xmm2, [esp + nb202nf_krf]
	movaps xmm3, [esp + nb202nf_krf]
	mulps  xmm0, [esp + nb202nf_rsqOH1] ;# krsq 
	mulps  xmm1, [esp + nb202nf_rsqOH2] ;# krsq 
	mulps  xmm2, [esp + nb202nf_rsqH1O] ;# krsq 
	mulps  xmm3, [esp + nb202nf_rsqH2O] ;# krsq 
	addps  xmm0, [esp + nb202nf_rinvOH1]	;# rinv+ krsq 
	addps  xmm1, [esp + nb202nf_rinvOH2]	;# rinv+ krsq 
	addps  xmm2, [esp + nb202nf_rinvH1O]	;# rinv+ krsq 
	addps  xmm3, [esp + nb202nf_rinvH2O]	;# rinv+ krsq 
	subps  xmm0, [esp + nb202nf_crf]
	subps  xmm1, [esp + nb202nf_crf]
	subps  xmm2, [esp + nb202nf_crf]
	subps  xmm3, [esp + nb202nf_crf]
	mulps  xmm0, [esp + nb202nf_qqOH] ;# voul=qq*(rinv+ krsq-crf) 
	mulps  xmm1, [esp + nb202nf_qqOH] ;# voul=qq*(rinv+ krsq-crf) 
	mulps  xmm2, [esp + nb202nf_qqOH] ;# voul=qq*(rinv+ krsq-crf) 
	mulps  xmm3, [esp + nb202nf_qqOH] ;# voul=qq*(rinv+ krsq-crf) 
	addps xmm6, xmm0
	addps xmm1, xmm2
	addps xmm6, xmm3
	addps xmm6, xmm1
		
	;# H-H interactions 
	movaps xmm0, [esp + nb202nf_krf]
	movaps xmm1, [esp + nb202nf_krf]
	movaps xmm2, [esp + nb202nf_krf]
	movaps xmm3, [esp + nb202nf_krf]
	mulps  xmm0, [esp + nb202nf_rsqH1H1] ;# krsq 
	mulps  xmm1, [esp + nb202nf_rsqH1H2] ;# krsq 
	mulps  xmm2, [esp + nb202nf_rsqH2H1] ;# krsq 
	mulps  xmm3, [esp + nb202nf_rsqH2H2] ;# krsq 
	addps  xmm0, [esp + nb202nf_rinvH1H1]	;# rinv+ krsq 
	addps  xmm1, [esp + nb202nf_rinvH1H2]	;# rinv+ krsq 
	addps  xmm2, [esp + nb202nf_rinvH2H1]	;# rinv+ krsq 
	addps  xmm3, [esp + nb202nf_rinvH2H2]	;# rinv+ krsq 
	subps  xmm0, [esp + nb202nf_crf]
	subps  xmm1, [esp + nb202nf_crf]
	subps  xmm2, [esp + nb202nf_crf]
	subps  xmm3, [esp + nb202nf_crf]
	mulps  xmm0, [esp + nb202nf_qqHH] ;# voul=qq*(rinv+ krsq-crf) 
	mulps  xmm1, [esp + nb202nf_qqHH] ;# voul=qq*(rinv+ krsq-crf) 
	mulps  xmm2, [esp + nb202nf_qqHH] ;# voul=qq*(rinv+ krsq-crf) 
	mulps  xmm3, [esp + nb202nf_qqHH] ;# voul=qq*(rinv+ krsq-crf) 
	addps xmm6, xmm0
	addps xmm1, xmm2
	addps xmm6, xmm3
	addps xmm6, xmm1
	movaps [esp + nb202nf_vctot], xmm6
		
	;# should we do one more iteration? 
	sub dword ptr [esp + nb202nf_innerk],  4
	jl    .nb202nf_single_check
	jmp   .nb202nf_unroll_loop
.nb202nf_single_check:
	add dword ptr [esp + nb202nf_innerk],  4
	jnz   .nb202nf_single_loop
	jmp   .nb202nf_updateouterdata
.nb202nf_single_loop:
	mov   edx, [esp + nb202nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb202nf_innerjjnr],  4	

	mov esi, [ebp + nb202nf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [esi + eax*4]		;# jxO  -  -  -
	movss xmm4, [esi + eax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [esi + eax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + nb202nf_ixO]     
	movaps  xmm1, [esp + nb202nf_iyO]
	movaps  xmm2, [esp + nb202nf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzO   0   jzH1 jzH2
	
	;# store all j coordinates in jO  
	movaps [esp + nb202nf_jxO], xmm3
	movaps [esp + nb202nf_jyO], xmm4
	movaps [esp + nb202nf_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 

	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [esp + nb202nf_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + nb202nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + nb202nf_half] ;# rinv iO - j water 
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	subps   xmm6, [esp + nb202nf_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	subps   xmm3, xmm7	;# xmm3=rinv-2*krsq 
	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [esp + nb202nf_qqOO]
	movhps  xmm4, [esp + nb202nf_qqOH]

	mulps xmm6, xmm4	;# vcoul  

	addps   xmm6, [esp + nb202nf_vctot]
	movaps  [esp + nb202nf_vctot], xmm6
	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + nb202nf_ixH1]
	movaps  xmm1, [esp + nb202nf_iyH1]
	movaps  xmm2, [esp + nb202nf_izH1]	
	movaps  xmm3, [esp + nb202nf_ixH2] 
	movaps  xmm4, [esp + nb202nf_iyH2] 
	movaps  xmm5, [esp + nb202nf_izH2] 
	subps   xmm0, [esp + nb202nf_jxO]
	subps   xmm1, [esp + nb202nf_jyO]
	subps   xmm2, [esp + nb202nf_jzO]
	subps   xmm3, [esp + nb202nf_jxO]
	subps   xmm4, [esp + nb202nf_jyO]
	subps   xmm5, [esp + nb202nf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + nb202nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + nb202nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + nb202nf_half] ;# rinv H2 - j water  

	mulps xmm0, [esp + nb202nf_krf] ;# krsq 
	mulps xmm4, [esp + nb202nf_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [esp + nb202nf_qqOH]
	movhps  xmm6, [esp + nb202nf_qqHH]
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [esp + nb202nf_crf]
	subps   xmm4, [esp + nb202nf_crf]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [esp + nb202nf_vctot]
	movaps  [esp + nb202nf_vctot], xmm4
	
	dec dword ptr [esp + nb202nf_innerk]
	jz    .nb202nf_updateouterdata
	jmp   .nb202nf_single_loop
.nb202nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb202nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb202nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb202nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb202nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb202nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb202nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb202nf_n], esi
        jmp .nb202nf_outer
.nb202nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb202nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb202nf_end
        ;# non-zero, do one more workunit
        jmp   .nb202nf_threadloop
.nb202nf_end:
	emms

	mov eax, [esp + nb202nf_nouter]
	mov ebx, [esp + nb202nf_ninner]
	mov ecx, [ebp + nb202nf_outeriter]
	mov edx, [ebp + nb202nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb202nf_salign]
	add esp, eax
	add esp, 744
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


		
