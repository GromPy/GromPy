;#
;# $Id: nb_kernel230_ia32_sse.intel_syntax.s,v 1.5 2008/04/06 19:43:57 lindahl Exp $
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




.globl nb_kernel230_ia32_sse
.globl _nb_kernel230_ia32_sse
nb_kernel230_ia32_sse:	
_nb_kernel230_ia32_sse:	
.equiv          nb230_p_nri,            8
.equiv          nb230_iinr,             12
.equiv          nb230_jindex,           16
.equiv          nb230_jjnr,             20
.equiv          nb230_shift,            24
.equiv          nb230_shiftvec,         28
.equiv          nb230_fshift,           32
.equiv          nb230_gid,              36
.equiv          nb230_pos,              40
.equiv          nb230_faction,          44
.equiv          nb230_charge,           48
.equiv          nb230_p_facel,          52
.equiv          nb230_argkrf,           56
.equiv          nb230_argcrf,           60
.equiv          nb230_Vc,               64
.equiv          nb230_type,             68
.equiv          nb230_p_ntype,          72
.equiv          nb230_vdwparam,         76
.equiv          nb230_Vvdw,             80
.equiv          nb230_p_tabscale,       84
.equiv          nb230_VFtab,            88
.equiv          nb230_invsqrta,         92
.equiv          nb230_dvda,             96
.equiv          nb230_p_gbtabscale,     100
.equiv          nb230_GBtab,            104
.equiv          nb230_p_nthreads,       108
.equiv          nb230_count,            112
.equiv          nb230_mtx,              116
.equiv          nb230_outeriter,        120
.equiv          nb230_inneriter,        124
.equiv          nb230_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb230_ix,               0
.equiv          nb230_iy,               16
.equiv          nb230_iz,               32
.equiv          nb230_iq,               48
.equiv          nb230_dx,               64
.equiv          nb230_dy,               80
.equiv          nb230_dz,               96
.equiv          nb230_c6,               112
.equiv          nb230_c12,              128
.equiv          nb230_tsc,              144
.equiv          nb230_fstmp,            160
.equiv          nb230_vctot,            176
.equiv          nb230_Vvdwtot,          192
.equiv          nb230_fix,              208
.equiv          nb230_fiy,              224
.equiv          nb230_fiz,              240
.equiv          nb230_half,             256
.equiv          nb230_three,            272
.equiv          nb230_two,              288
.equiv          nb230_krf,              304
.equiv          nb230_crf,              320
.equiv          nb230_is3,              336
.equiv          nb230_ii3,              340
.equiv          nb230_ntia,             344
.equiv          nb230_innerjjnr,        348
.equiv          nb230_innerk,           352
.equiv          nb230_n,                356
.equiv          nb230_nn1,              360
.equiv          nb230_nri,              364
.equiv          nb230_facel,            368
.equiv          nb230_ntype,            372
.equiv          nb230_nouter,           376
.equiv          nb230_ninner,           380
.equiv          nb230_salign,           384
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  400		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb230_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb230_p_nri]
	mov esi, [ebp + nb230_p_facel]
	mov edi, [ebp + nb230_p_ntype]
	mov ecx, [ecx]
	mov esi, [esi]
	mov edi, [edi]
	mov [esp + nb230_nri], ecx
	mov [esp + nb230_facel], esi
	mov [esp + nb230_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb230_nouter], eax
	mov [esp + nb230_ninner], eax
	
	mov eax, [ebp + nb230_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb230_tsc], xmm3

	mov esi, [ebp + nb230_argkrf]
	mov edi, [ebp + nb230_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]	
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb230_krf], xmm5
	movaps [esp + nb230_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb230_half], eax
	movss xmm1, [esp + nb230_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb230_half],  xmm1
	movaps [esp + nb230_two],  xmm2
	movaps [esp + nb230_three],  xmm3

.nb230_threadloop:
        mov   esi, [ebp + nb230_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb230_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb230_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb230_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb230_n], eax
        mov [esp + nb230_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb230_outerstart
        jmp .nb230_end

.nb230_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb230_nouter]
	mov [esp + nb230_nouter], ebx

.nb230_outer:
	mov   eax, [ebp + nb230_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb230_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb230_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb230_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb230_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb230_facel]
	shufps xmm3, xmm3, 0

    	mov   edx, [ebp + nb230_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb230_ntype]
    	shl   edx, 1
    	mov   [esp + nb230_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb230_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb230_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb230_ix], xmm0
	movaps [esp + nb230_iy], xmm1
	movaps [esp + nb230_iz], xmm2

	mov   [esp + nb230_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb230_vctot], xmm4
	movaps [esp + nb230_Vvdwtot], xmm4
	movaps [esp + nb230_fix], xmm4
	movaps [esp + nb230_fiy], xmm4
	movaps [esp + nb230_fiz], xmm4
	
	mov   eax, [ebp + nb230_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb230_pos]
	mov   edi, [ebp + nb230_faction]	
	mov   eax, [ebp + nb230_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb230_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb230_ninner]
	mov   [esp + nb230_ninner], ecx
	add   edx, 0
	mov   [esp + nb230_innerk], edx    ;# number of innerloop atoms 
	jge   .nb230_unroll_loop
	jmp   .nb230_finish_inner
.nb230_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb230_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb230_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb230_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + nb230_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + nb230_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb230_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb230_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb230_c6], xmm4
	movaps [esp + nb230_c12], xmm6
	
	mov esi, [ebp + nb230_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [esp + nb230_ix]
	movaps xmm5, [esp + nb230_iy]
	movaps xmm6, [esp + nb230_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb230_dx], xmm4
	movaps [esp + nb230_dy], xmm5
	movaps [esp + nb230_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 
	
	movaps xmm7, [esp + nb230_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb230_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb230_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 	
	movaps xmm1, xmm0
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb230_crf]
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	mulps  xmm7, [esp + nb230_two]  ;# 2*krsq
	addps  xmm6, [esp + nb230_vctot]
	movaps [esp + nb230_vctot], xmm6

	subps  xmm0, xmm7  ;# rinv-2*krsq
	mulps  xmm0, xmm3  ;# qq*(rinv-2*krsq)
	mulps  xmm0, xmm1  ;# qq*(rinv-2*krsq)*rinv 
	movaps [esp + nb230_fstmp], xmm0
	
	;# LJ table
	mulps  xmm4, xmm1  ;# r
	mulps  xmm4, [esp + nb230_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + nb230_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [esp + nb230_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [esp + nb230_c6]
	mulps  xmm7, xmm4	 ;# fijD 
	mulps  xmm5, xmm4	 ;# Vvdw6 
	movaps xmm3, [esp + nb230_fstmp]
	mulps  xmm7, [esp + nb230_tsc]
	subps  xmm3, xmm7
	
	;# put scalar force on stack Update Vvdwtot directly 
	addps  xmm5, [esp + nb230_Vvdwtot]
	movaps [esp + nb230_fstmp], xmm3
	movaps [esp + nb230_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [esp + nb230_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [esp + nb230_c12]
	mulps  xmm7, xmm4 ;# fijR 
	mulps  xmm5, xmm4 ;# Vvdw12 
	movaps xmm3, [esp + nb230_fstmp]
	mulps  xmm7, [esp + nb230_tsc]
	subps  xmm3, xmm7
	
	addps  xmm5, [esp + nb230_Vvdwtot]
	movaps [esp + nb230_Vvdwtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm3, xmm0    

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3


	movaps xmm0, [esp + nb230_dx]
	movaps xmm1, [esp + nb230_dy]
	movaps xmm2, [esp + nb230_dz]

	mov    edi, [ebp + nb230_faction]
	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb230_fix]
	movaps xmm4, [esp + nb230_fiy]
	movaps xmm5, [esp + nb230_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb230_fix], xmm3
	movaps [esp + nb230_fiy], xmm4
	movaps [esp + nb230_fiz], xmm5
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 136  ;# constant 10001000
	shufps xmm4, xmm6, 221  ;# constant 11011101			      

	;# now xmm3-xmm5 contains fjx, fjy, fjz 
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	;# unpack them back so we can store them - first x & y in xmm3/xmm4 

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	;# xmm6(l)=x & y for j1, (h) for j2 
	;# xmm3(l)=x & y for j3, (h) for j4 
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;# and the z forces 
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 229  ;# constant 11100101
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 234  ;# constant 11101010
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 255  ;# constant 11111111
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb230_innerk],  4
	jl    .nb230_finish_inner
	jmp   .nb230_unroll_loop
.nb230_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb230_innerk],  4
	mov   edx, [esp + nb230_innerk]
	and   edx, 2
	jnz   .nb230_dopair
	jmp   .nb230_checksingle
.nb230_dopair:	
	mov esi, [ebp + nb230_charge]

    mov   ecx, [esp + nb230_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb230_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 12 ;# constant 00001100 
	shufps xmm3, xmm3, 88 ;# constant 01011000 ;# xmm3(0,1) has the charges 

	mov esi, [ebp + nb230_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + nb230_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb230_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + nb230_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# constant 00001000 	
	shufps xmm6, xmm6, 13 ;# constant 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + nb230_c6], xmm4
	movaps [esp + nb230_c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb230_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	mov    edi, [ebp + nb230_faction]
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb230_ix]
	movaps xmm5, [esp + nb230_iy]
	movaps xmm6, [esp + nb230_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb230_dx], xmm4
	movaps [esp + nb230_dy], xmm5
	movaps [esp + nb230_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movaps xmm7, [esp + nb230_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb230_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb230_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 	
	movaps xmm1, xmm0
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb230_crf]
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	mulps  xmm7, [esp + nb230_two]  ;# 2*krsq
	addps  xmm6, [esp + nb230_vctot]
	movaps [esp + nb230_vctot], xmm6

	subps  xmm0, xmm7  ;# rinv-2*krsq
	mulps  xmm0, xmm3  ;# qq*(rinv-2*krsq)
	mulps  xmm0, xmm1  ;# qq*(rinv-2*krsq)*rinv 
	movaps [esp + nb230_fstmp], xmm0

	;# LJ table
	mulps  xmm4, xmm1  ;# r
	mulps  xmm4, [esp + nb230_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb230_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movhps xmm5, [esi + ebx*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [esp + nb230_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [esp + nb230_c6]
	mulps  xmm7, xmm4	 ;# fijD 
	mulps  xmm5, xmm4	 ;# Vvdw6 
	movaps xmm3, [esp + nb230_fstmp]
	mulps  xmm7, [esp + nb230_tsc]
	subps  xmm3, xmm7
	
	;# put scalar force on stack Update Vvdwtot directly 
	addps  xmm5, [esp + nb230_Vvdwtot]
	movaps [esp + nb230_fstmp], xmm3
	movaps [esp + nb230_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [esp + nb230_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [esp + nb230_c12]
	mulps  xmm7, xmm4 ;# fijR 
	mulps  xmm5, xmm4 ;# Vvdw12 
	movaps xmm3, [esp + nb230_fstmp]
	mulps  xmm7, [esp + nb230_tsc]
	subps  xmm3, xmm7
	
	addps  xmm5, [esp + nb230_Vvdwtot]
	movaps [esp + nb230_Vvdwtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm3, xmm0    

	movd eax, mm0	
	movd ebx, mm1

	movaps xmm0, [esp + nb230_dx]
	movaps xmm1, [esp + nb230_dy]
	movaps xmm2, [esp + nb230_dz]

	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb230_fix]
	movaps xmm4, [esp + nb230_fiy]
	movaps xmm5, [esp + nb230_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb230_fix], xmm3
	movaps [esp + nb230_fiy], xmm4
	movaps [esp + nb230_fiz], xmm5
	;# update the fj's 
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 225  ;# constant 11100001
	shufps  xmm1, xmm1, 225  ;# constant 11100001
	shufps  xmm2, xmm2, 225  ;# constant 11100001

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.nb230_checksingle:				
	mov   edx, [esp + nb230_innerk]
	and   edx, 1
	jnz    .nb230_dosingle
	jmp    .nb230_updateouterdata
.nb230_dosingle:			
	mov esi, [ebp + nb230_charge]
	mov edi, [ebp + nb230_pos]
	mov   ecx, [esp + nb230_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 	

	mov esi, [ebp + nb230_type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + nb230_vdwparam]
	shl ecx, 1
	add ecx, [esp + nb230_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# constant 11111100	
	shufps xmm6, xmm6, 253  ;# constant 11111101	
			
	movaps [esp + nb230_c6], xmm4
	movaps [esp + nb230_c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb230_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb230_ix]
	movaps xmm5, [esp + nb230_iy]
	movaps xmm6, [esp + nb230_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb230_dx], xmm4
	movaps [esp + nb230_dy], xmm5
	movaps [esp + nb230_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movss xmm7, [esp + nb230_krf]
	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movss xmm2, xmm5
	mulss xmm5, xmm5
	movss xmm1, [esp + nb230_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movss xmm0, [esp + nb230_half]
	mulss  xmm7, xmm4	;# xmm7=krsq 
	subss xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 	
	movss xmm1, xmm0
	movss xmm6, xmm0
	addss  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subss  xmm6, [esp + nb230_crf]
	mulss  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	mulss  xmm7, [esp + nb230_two]  ;# 2*krsq
	addss  xmm6, [esp + nb230_vctot]
	movss [esp + nb230_vctot], xmm6

	subss  xmm0, xmm7  ;# rinv-2*krsq
	mulss  xmm0, xmm3  ;# qq*(rinv-2*krsq)
	mulss  xmm0, xmm1  ;# qq*(rinv-2*krsq)*rinv 
	movss [esp + nb230_fstmp], xmm0
	
	;# LJ table
	mulss  xmm4, xmm1  ;# r
	mulss  xmm4, [esp + nb230_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	

	mov  esi, [ebp + nb230_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [esp + nb230_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss  xmm4, [esp + nb230_c6]
	mulss  xmm7, xmm4	 ;# fijD 
	mulss  xmm5, xmm4	 ;# Vvdw6 
	movss  xmm3, [esp + nb230_fstmp]
	mulss  xmm7, [esp + nb230_tsc]
	subss  xmm3, xmm7
	
	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [esp + nb230_Vvdwtot]
	movss [esp + nb230_fstmp], xmm3
	movss [esp + nb230_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [esp + nb230_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss  xmm4, [esp + nb230_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	movss  xmm3, [esp + nb230_fstmp]
	mulss  xmm7, [esp + nb230_tsc]
	subss  xmm3, xmm7
			
	addss  xmm5, [esp + nb230_Vvdwtot]
	movss [esp + nb230_Vvdwtot], xmm5

	mulss xmm3, xmm0    

	movd eax, mm0	

	movaps xmm0, [esp + nb230_dx]
	movaps xmm1, [esp + nb230_dy]
	movaps xmm2, [esp + nb230_dz]

	mulss  xmm0, xmm3
	mulss  xmm1, xmm3
	mulss  xmm2, xmm3
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb230_fix]
	movaps xmm4, [esp + nb230_fiy]
	movaps xmm5, [esp + nb230_fiz]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movaps [esp + nb230_fix], xmm3
	movaps [esp + nb230_fiy], xmm4
	movaps [esp + nb230_fiz], xmm5
	;# update fj 
	mov edi, [ebp + nb230_faction]
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.nb230_updateouterdata:
	mov   ecx, [esp + nb230_ii3]
	mov   edi, [ebp + nb230_faction]
	mov   esi, [ebp + nb230_fshift]
	mov   edx, [esp + nb230_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb230_fix]
	movaps xmm1, [esp + nb230_fiy]
	movaps xmm2, [esp + nb230_fiz]

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

	;# increment fshift force  
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;# get n from stack
	mov esi, [esp + nb230_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb230_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb230_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb230_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb230_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb230_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb230_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb230_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb230_n], esi
        jmp .nb230_outer
.nb230_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb230_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb230_end
        ;# non-zero, do one more workunit
        jmp   .nb230_threadloop
.nb230_end:
	emms

	mov eax, [esp + nb230_nouter]
	mov ebx, [esp + nb230_ninner]
	mov ecx, [ebp + nb230_outeriter]
	mov edx, [ebp + nb230_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb230_salign]
	add esp, eax
	add esp,  400
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret







.globl nb_kernel230nf_ia32_sse
.globl _nb_kernel230nf_ia32_sse
nb_kernel230nf_ia32_sse:	
_nb_kernel230nf_ia32_sse:	
.equiv          nb230nf_p_nri,            8
.equiv          nb230nf_iinr,             12
.equiv          nb230nf_jindex,           16
.equiv          nb230nf_jjnr,             20
.equiv          nb230nf_shift,            24
.equiv          nb230nf_shiftvec,         28
.equiv          nb230nf_fshift,           32
.equiv          nb230nf_gid,              36
.equiv          nb230nf_pos,              40
.equiv          nb230nf_faction,          44
.equiv          nb230nf_charge,           48
.equiv          nb230nf_p_facel,          52
.equiv          nb230nf_argkrf,           56
.equiv          nb230nf_argcrf,           60
.equiv          nb230nf_Vc,               64
.equiv          nb230nf_type,             68
.equiv          nb230nf_p_ntype,          72
.equiv          nb230nf_vdwparam,         76
.equiv          nb230nf_Vvdw,             80
.equiv          nb230nf_p_tabscale,       84
.equiv          nb230nf_VFtab,            88
.equiv          nb230nf_invsqrta,         92
.equiv          nb230nf_dvda,             96
.equiv          nb230nf_p_gbtabscale,     100
.equiv          nb230nf_GBtab,            104
.equiv          nb230nf_p_nthreads,       108
.equiv          nb230nf_count,            112
.equiv          nb230nf_mtx,              116
.equiv          nb230nf_outeriter,        120
.equiv          nb230nf_inneriter,        124
.equiv          nb230nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb230nf_ix,               0
.equiv          nb230nf_iy,               16
.equiv          nb230nf_iz,               32
.equiv          nb230nf_iq,               48
.equiv          nb230nf_dx,               64
.equiv          nb230nf_dy,               80
.equiv          nb230nf_dz,               96
.equiv          nb230nf_c6,               112
.equiv          nb230nf_c12,              128
.equiv          nb230nf_tsc,              144
.equiv          nb230nf_vctot,            160
.equiv          nb230nf_Vvdwtot,          176
.equiv          nb230nf_half,             192
.equiv          nb230nf_three,            208
.equiv          nb230nf_krf,              224
.equiv          nb230nf_crf,              240
.equiv          nb230nf_is3,              256
.equiv          nb230nf_ii3,              260
.equiv          nb230nf_ntia,             264
.equiv          nb230nf_innerjjnr,        268
.equiv          nb230nf_innerk,           272
.equiv          nb230nf_n,                276
.equiv          nb230nf_nn1,              280
.equiv          nb230nf_nri,              284
.equiv          nb230nf_facel,            288
.equiv          nb230nf_ntype,            292
.equiv          nb230nf_nouter,           296
.equiv          nb230nf_ninner,           300
.equiv          nb230nf_salign,           304
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  304		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb230nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb230nf_p_nri]
	mov esi, [ebp + nb230nf_p_facel]
	mov edi, [ebp + nb230nf_p_ntype]
	mov ecx, [ecx]
	mov esi, [esi]
	mov edi, [edi]
	mov [esp + nb230nf_nri], ecx
	mov [esp + nb230nf_facel], esi
	mov [esp + nb230nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb230nf_nouter], eax
	mov [esp + nb230nf_ninner], eax
	
	mov eax, [ebp + nb230nf_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb230nf_tsc], xmm3

	mov esi, [ebp + nb230nf_argkrf]
	mov edi, [ebp + nb230nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]	
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb230nf_krf], xmm5
	movaps [esp + nb230nf_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb230nf_half], eax
	movss xmm1, [esp + nb230nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb230nf_half],  xmm1
	movaps [esp + nb230nf_three],  xmm3

.nb230nf_threadloop:
        mov   esi, [ebp + nb230nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb230nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb230nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb230nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb230nf_n], eax
        mov [esp + nb230nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb230nf_outerstart
        jmp .nb230nf_end

.nb230nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb230nf_nouter]
	mov [esp + nb230nf_nouter], ebx

.nb230nf_outer:
	mov   eax, [ebp + nb230nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb230nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb230nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb230nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb230nf_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb230nf_facel]
	shufps xmm3, xmm3, 0

    	mov   edx, [ebp + nb230nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb230nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb230nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb230nf_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb230nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb230nf_ix], xmm0
	movaps [esp + nb230nf_iy], xmm1
	movaps [esp + nb230nf_iz], xmm2

	mov   [esp + nb230nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb230nf_vctot], xmm4
	movaps [esp + nb230nf_Vvdwtot], xmm4
		
	mov   eax, [ebp + nb230nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb230nf_pos]
	mov   eax, [ebp + nb230nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb230nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb230nf_ninner]
	mov   [esp + nb230nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb230nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb230nf_unroll_loop
	jmp   .nb230nf_finish_inner
.nb230nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb230nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb230nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb230nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + nb230nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + nb230nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb230nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb230nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb230nf_c6], xmm4
	movaps [esp + nb230nf_c12], xmm6
	
	mov esi, [ebp + nb230nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [esp + nb230nf_ix]
	movaps xmm5, [esp + nb230nf_iy]
	movaps xmm6, [esp + nb230nf_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 
	
	movaps xmm7, [esp + nb230nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb230nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb230nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 	
	movaps xmm1, xmm0
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb230nf_crf]
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 

	addps  xmm6, [esp + nb230nf_vctot]
	movaps [esp + nb230nf_vctot], xmm6
	
	;# LJ table
	mulps  xmm4, xmm1  ;# r
	mulps  xmm4, [esp + nb230nf_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	mov  esi, [ebp + nb230nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [esp + nb230nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 
	
	;# Update Vvdwtot directly 
	addps  xmm5, [esp + nb230nf_Vvdwtot]
	movaps [esp + nb230nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [esp + nb230nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	
	addps  xmm5, [esp + nb230nf_Vvdwtot]
	movaps [esp + nb230nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb230nf_innerk],  4
	jl    .nb230nf_finish_inner
	jmp   .nb230nf_unroll_loop
.nb230nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb230nf_innerk],  4
	mov   edx, [esp + nb230nf_innerk]
	and   edx, 2
	jnz   .nb230nf_dopair
	jmp   .nb230nf_checksingle
.nb230nf_dopair:	
	mov esi, [ebp + nb230nf_charge]

    mov   ecx, [esp + nb230nf_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb230nf_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 12 ;# constant 00001100 
	shufps xmm3, xmm3, 88 ;# constant 01011000 ;# xmm3(0,1) has the charges 

	mov esi, [ebp + nb230nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + nb230nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb230nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + nb230nf_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# constant 00001000 	
	shufps xmm6, xmm6, 13 ;# constant 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + nb230nf_c6], xmm4
	movaps [esp + nb230nf_c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb230nf_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb230nf_ix]
	movaps xmm5, [esp + nb230nf_iy]
	movaps xmm6, [esp + nb230nf_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movaps xmm7, [esp + nb230nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb230nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb230nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 	
	movaps xmm1, xmm0
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb230nf_crf]
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [esp + nb230nf_vctot]
	movaps [esp + nb230nf_vctot], xmm6

	;# LJ table
	mulps  xmm4, xmm1  ;# r
	mulps  xmm4, [esp + nb230nf_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	mov  esi, [ebp + nb230nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movhps xmm5, [esi + ebx*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [esp + nb230nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 
	
	;# Update Vvdwtot directly 
	addps  xmm5, [esp + nb230nf_Vvdwtot]
	movaps [esp + nb230nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [esp + nb230nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	
	addps  xmm5, [esp + nb230nf_Vvdwtot]
	movaps [esp + nb230nf_Vvdwtot], xmm5

.nb230nf_checksingle:				
	mov   edx, [esp + nb230nf_innerk]
	and   edx, 1
	jnz    .nb230nf_dosingle
	jmp    .nb230nf_updateouterdata
.nb230nf_dosingle:			
	mov esi, [ebp + nb230nf_charge]
	mov edi, [ebp + nb230nf_pos]
	mov   ecx, [esp + nb230nf_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 	

	mov esi, [ebp + nb230nf_type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + nb230nf_vdwparam]
	shl ecx, 1
	add ecx, [esp + nb230nf_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# constant 11111100	
	shufps xmm6, xmm6, 253  ;# constant 11111101	
			
	movaps [esp + nb230nf_c6], xmm4
	movaps [esp + nb230nf_c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb230nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb230nf_ix]
	movaps xmm5, [esp + nb230nf_iy]
	movaps xmm6, [esp + nb230nf_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movss xmm7, [esp + nb230nf_krf]
	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movss xmm2, xmm5
	mulss xmm5, xmm5
	movss xmm1, [esp + nb230nf_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movss xmm0, [esp + nb230nf_half]
	mulss  xmm7, xmm4	;# xmm7=krsq 
	subss xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 	
	movss xmm1, xmm0
	movss xmm6, xmm0
	addss  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subss  xmm6, [esp + nb230nf_crf]
	mulss  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	addss  xmm6, [esp + nb230nf_vctot]
	movss [esp + nb230nf_vctot], xmm6
	
	;# LJ table
	mulss  xmm4, xmm1  ;# r
	mulss  xmm4, [esp + nb230nf_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	

	mov  esi, [ebp + nb230nf_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss  xmm4, [esp + nb230nf_c6]
	mulss  xmm5, xmm4	 ;# Vvdw6 
	
	;# Update Vvdwtot directly 
	addss  xmm5, [esp + nb230nf_Vvdwtot]
	movss [esp + nb230nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss  xmm4, [esp + nb230nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 
	
	addss  xmm5, [esp + nb230nf_Vvdwtot]
	movss [esp + nb230nf_Vvdwtot], xmm5


.nb230nf_updateouterdata:

	;# get n from stack
	mov esi, [esp + nb230nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb230nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb230nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb230nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb230nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb230nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb230nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb230nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb230nf_n], esi
        jmp .nb230nf_outer
.nb230nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb230nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb230nf_end
        ;# non-zero, do one more workunit
        jmp   .nb230nf_threadloop
.nb230nf_end:
	emms

	mov eax, [esp + nb230nf_nouter]
	mov ebx, [esp + nb230nf_ninner]
	mov ecx, [ebp + nb230nf_outeriter]
	mov edx, [ebp + nb230nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb230nf_salign]
	add esp, eax
	add esp,  304
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret





