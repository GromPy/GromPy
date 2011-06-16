##
## $Id: nb_kernel311_x86_64_sse2.s,v 1.8 2006/09/22 08:50:28 lindahl Exp $
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##





.globl nb_kernel311_x86_64_sse2
.globl _nb_kernel311_x86_64_sse2
nb_kernel311_x86_64_sse2:       
_nb_kernel311_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb311_fshift, 16
.set nb311_gid, 24
.set nb311_pos, 32
.set nb311_faction, 40
.set nb311_charge, 48
.set nb311_p_facel, 56
.set nb311_argkrf, 64
.set nb311_argcrf, 72
.set nb311_Vc, 80
.set nb311_type, 88
.set nb311_p_ntype, 96
.set nb311_vdwparam, 104
.set nb311_Vvdw, 112
.set nb311_p_tabscale, 120
.set nb311_VFtab, 128
.set nb311_invsqrta, 136
.set nb311_dvda, 144
.set nb311_p_gbtabscale, 152
.set nb311_GBtab, 160
.set nb311_p_nthreads, 168
.set nb311_count, 176
.set nb311_mtx, 184
.set nb311_outeriter, 192
.set nb311_inneriter, 200
.set nb311_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb311_ixO, 0
.set nb311_iyO, 16
.set nb311_izO, 32
.set nb311_ixH1, 48
.set nb311_iyH1, 64
.set nb311_izH1, 80
.set nb311_ixH2, 96
.set nb311_iyH2, 112
.set nb311_izH2, 128
.set nb311_iqO, 144
.set nb311_iqH, 160
.set nb311_dxO, 176
.set nb311_dyO, 192
.set nb311_dzO, 208
.set nb311_dxH1, 224
.set nb311_dyH1, 240
.set nb311_dzH1, 256
.set nb311_dxH2, 272
.set nb311_dyH2, 288
.set nb311_dzH2, 304
.set nb311_qqO, 320
.set nb311_qqH, 336
.set nb311_rinvO, 352
.set nb311_rinvH1, 368
.set nb311_rinvH2, 384
.set nb311_rO, 400
.set nb311_rH1, 416
.set nb311_rH2, 432
.set nb311_tsc, 448
.set nb311_two, 464
.set nb311_c6, 480
.set nb311_c12, 496
.set nb311_six, 512
.set nb311_twelve, 528
.set nb311_vctot, 544
.set nb311_Vvdwtot, 560
.set nb311_fixO, 576
.set nb311_fiyO, 592
.set nb311_fizO, 608
.set nb311_fixH1, 624
.set nb311_fiyH1, 640
.set nb311_fizH1, 656
.set nb311_fixH2, 672
.set nb311_fiyH2, 688
.set nb311_fizH2, 704
.set nb311_fjx, 720
.set nb311_fjy, 736
.set nb311_fjz, 752
.set nb311_half, 768
.set nb311_three, 784
.set nb311_is3, 800
.set nb311_ii3, 804
.set nb311_nri, 808
.set nb311_iinr, 816
.set nb311_jindex, 824
.set nb311_jjnr, 832
.set nb311_shift, 840
.set nb311_shiftvec, 848
.set nb311_facel, 856
.set nb311_innerjjnr, 864
.set nb311_ntia, 872
.set nb311_innerk, 876
.set nb311_n, 880
.set nb311_nn1, 884
.set nb311_nouter, 888
.set nb311_ninner, 892
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $904,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb311_nouter(%rsp)
        movl %eax,nb311_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb311_nri(%rsp)
        movq %rsi,nb311_iinr(%rsp)
        movq %rdx,nb311_jindex(%rsp)
        movq %rcx,nb311_jjnr(%rsp)
        movq %r8,nb311_shift(%rsp)
        movq %r9,nb311_shiftvec(%rsp)
        movq nb311_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb311_facel(%rsp)

        movq nb311_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb311_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb311_half(%rsp)
        movl %ebx,nb311_half+4(%rsp)
        movsd nb311_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## six
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## twelve
        movapd %xmm1,nb311_half(%rsp)
        movapd %xmm2,nb311_two(%rsp)
        movapd %xmm3,nb311_three(%rsp)
        movapd %xmm4,nb311_six(%rsp)
        movapd %xmm5,nb311_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb311_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb311_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4
        movq nb311_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb311_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb311_iqO(%rsp)
        movapd %xmm4,nb311_iqH(%rsp)

        movq  nb311_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb311_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb311_ntia(%rsp)
_nb_kernel311_x86_64_sse2.nb311_threadloop: 
        movq  nb311_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel311_x86_64_sse2.nb311_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311_x86_64_sse2.nb311_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311_n(%rsp)
        movl %ebx,nb311_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311_x86_64_sse2.nb311_outerstart
        jmp _nb_kernel311_x86_64_sse2.nb311_end

_nb_kernel311_x86_64_sse2.nb311_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311_nouter(%rsp),%ebx
        movl %ebx,nb311_nouter(%rsp)

_nb_kernel311_x86_64_sse2.nb311_outer: 
        movq  nb311_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb311_is3(%rsp)      ## store is3 

        movq  nb311_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb311_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb311_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb311_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb311_ixO(%rsp)
        movapd %xmm4,nb311_iyO(%rsp)
        movapd %xmm5,nb311_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb311_ixH1(%rsp)
        movapd %xmm1,nb311_iyH1(%rsp)
        movapd %xmm2,nb311_izH1(%rsp)
        movapd %xmm3,nb311_ixH2(%rsp)
        movapd %xmm4,nb311_iyH2(%rsp)
        movapd %xmm5,nb311_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb311_vctot(%rsp)
        movapd %xmm4,nb311_Vvdwtot(%rsp)
        movapd %xmm4,nb311_fixO(%rsp)
        movapd %xmm4,nb311_fiyO(%rsp)
        movapd %xmm4,nb311_fizO(%rsp)
        movapd %xmm4,nb311_fixH1(%rsp)
        movapd %xmm4,nb311_fiyH1(%rsp)
        movapd %xmm4,nb311_fizH1(%rsp)
        movapd %xmm4,nb311_fixH2(%rsp)
        movapd %xmm4,nb311_fiyH2(%rsp)
        movapd %xmm4,nb311_fizH2(%rsp)

        movq  nb311_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb311_pos(%rbp),%rsi
        movq  nb311_faction(%rbp),%rdi
        movq  nb311_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb311_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb311_ninner(%rsp),%ecx
        movl  %ecx,nb311_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb311_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel311_x86_64_sse2.nb311_unroll_loop
        jmp   _nb_kernel311_x86_64_sse2.nb311_checksingle
_nb_kernel311_x86_64_sse2.nb311_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb311_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb311_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb311_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311_iqO(%rsp),%xmm3
        mulpd  nb311_iqH(%rsp),%xmm4
        movapd  %xmm3,nb311_qqO(%rsp)
        movapd  %xmm4,nb311_qqH(%rsp)

        movq nb311_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb311_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb311_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb311_c6(%rsp)
        movapd %xmm6,nb311_c12(%rsp)

        movq nb311_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb311_ixO(%rsp),%xmm0
    subpd nb311_iyO(%rsp),%xmm1
    subpd nb311_izO(%rsp),%xmm2
    subpd nb311_ixH1(%rsp),%xmm3
    subpd nb311_iyH1(%rsp),%xmm4
    subpd nb311_izH1(%rsp),%xmm5
    subpd nb311_ixH2(%rsp),%xmm6
    subpd nb311_iyH2(%rsp),%xmm7
    subpd nb311_izH2(%rsp),%xmm8

        movapd %xmm0,nb311_dxO(%rsp)
        movapd %xmm1,nb311_dyO(%rsp)
        movapd %xmm2,nb311_dzO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb311_dxH1(%rsp)
        movapd %xmm4,nb311_dyH1(%rsp)
        movapd %xmm5,nb311_dzH1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb311_dxH2(%rsp)
        movapd %xmm7,nb311_dyH2(%rsp)
        movapd %xmm8,nb311_dzH2(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb311_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb311_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvO
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2  

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb311_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb311_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvO
        mulpd   %xmm15,%xmm10 ##   rinvH1
    mulpd   %xmm15,%xmm11 ##   rinvH2

        movapd  %xmm9,nb311_rinvO(%rsp)
        movapd  %xmm10,nb311_rinvH1(%rsp)
        movapd  %xmm11,nb311_rinvH2(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb311_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb311_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,%xmm12 ## epsO
    movapd    %xmm3,%xmm13 ## epsH1
    movapd    %xmm6,%xmm14 ## epsH2

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11



    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1 ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  nb311_qqO(%rsp),%xmm1     ## VV*qq = vcoul
    mulpd  nb311_qqH(%rsp),%xmm5
    mulpd  nb311_qqH(%rsp),%xmm9
    mulpd  nb311_qqO(%rsp),%xmm3      ## FF*qq = fij
    mulpd  nb311_qqH(%rsp),%xmm7
    mulpd  nb311_qqH(%rsp),%xmm11

    ## calculate LJ
    movapd nb311_rinvO(%rsp),%xmm12
    mulpd  %xmm12,%xmm12 ## rinvsq
    movapd %xmm12,%xmm13 ## rinvsq
    mulpd  %xmm12,%xmm12 ## rinv4
    mulpd  %xmm13,%xmm12 ## rinv6
    movapd %xmm12,%xmm13 ## rinv6
    mulpd  %xmm12,%xmm12 ## rinv12
        mulpd  nb311_c6(%rsp),%xmm13
        mulpd  nb311_c12(%rsp),%xmm12
    movapd %xmm12,%xmm14
    subpd  %xmm13,%xmm14

        addpd  nb311_Vvdwtot(%rsp),%xmm14
        mulpd  nb311_six(%rsp),%xmm13
        mulpd  nb311_twelve(%rsp),%xmm12
        movapd %xmm14,nb311_Vvdwtot(%rsp)
    subpd  %xmm13,%xmm12 ## LJ fscal    
    mulpd  nb311_rinvO(%rsp),%xmm12

    ## accumulate vctot
    addpd  nb311_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb311_vctot(%rsp)

    movapd nb311_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## FF*tabscale
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    subpd  %xmm3,%xmm12
    mulpd  nb311_rinvO(%rsp),%xmm12
    mulpd  nb311_rinvH1(%rsp),%xmm7
    mulpd  nb311_rinvH2(%rsp),%xmm10

    xorpd  %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd  %xmm7,%xmm8
    subpd  %xmm10,%xmm11
    ## move j forces to xmm0-xmm2
    movq nb311_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm12,%xmm3
    movapd %xmm12,%xmm4
    movapd %xmm12,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb311_dxO(%rsp),%xmm3
        mulpd nb311_dyO(%rsp),%xmm4
        mulpd nb311_dzO(%rsp),%xmm5
        mulpd nb311_dxH1(%rsp),%xmm7
        mulpd nb311_dyH1(%rsp),%xmm8
        mulpd nb311_dzH1(%rsp),%xmm9
        mulpd nb311_dxH2(%rsp),%xmm10
        mulpd nb311_dyH2(%rsp),%xmm11
        mulpd nb311_dzH2(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb311_fixO(%rsp),%xmm3
    addpd nb311_fiyO(%rsp),%xmm4
    addpd nb311_fizO(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb311_fixH1(%rsp),%xmm7
    addpd nb311_fiyH1(%rsp),%xmm8
    addpd nb311_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb311_fixH2(%rsp),%xmm10
    addpd nb311_fiyH2(%rsp),%xmm11
    addpd nb311_fizH2(%rsp),%xmm12

    movapd %xmm3,nb311_fixO(%rsp)
    movapd %xmm4,nb311_fiyO(%rsp)
    movapd %xmm5,nb311_fizO(%rsp)
    movapd %xmm7,nb311_fixH1(%rsp)
    movapd %xmm8,nb311_fiyH1(%rsp)
    movapd %xmm9,nb311_fizH1(%rsp)
    movapd %xmm10,nb311_fixH2(%rsp)
    movapd %xmm11,nb311_fiyH2(%rsp)
    movapd %xmm12,nb311_fizH2(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb311_innerk(%rsp)
        jl    _nb_kernel311_x86_64_sse2.nb311_checksingle
        jmp   _nb_kernel311_x86_64_sse2.nb311_unroll_loop
_nb_kernel311_x86_64_sse2.nb311_checksingle: 
        movl  nb311_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel311_x86_64_sse2.nb311_dosingle
        jmp   _nb_kernel311_x86_64_sse2.nb311_updateouterdata
_nb_kernel311_x86_64_sse2.nb311_dosingle: 
        movq  nb311_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb311_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311_iqO(%rsp),%xmm3
        mulpd  nb311_iqH(%rsp),%xmm4

        movapd  %xmm3,nb311_qqO(%rsp)
        movapd  %xmm4,nb311_qqH(%rsp)

        movq nb311_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb311_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb311_ntia(%rsp),%edi
        addl %edi,%r8d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb311_c6(%rsp)
        movapd %xmm6,nb311_c12(%rsp)

        movq nb311_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb311_ixO(%rsp),%xmm4
        subsd nb311_iyO(%rsp),%xmm5
        subsd nb311_izO(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb311_dxO(%rsp)
        movapd %xmm5,nb311_dyO(%rsp)
        movapd %xmm6,nb311_dzO(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move j coords to xmm4-xmm6 
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm5
        movapd %xmm2,%xmm6

        ## calc dr 
        subsd nb311_ixH1(%rsp),%xmm4
        subsd nb311_iyH1(%rsp),%xmm5
        subsd nb311_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb311_dxH1(%rsp)
        movapd %xmm5,nb311_dyH1(%rsp)
        movapd %xmm6,nb311_dzH1(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move j coords to xmm3-xmm5
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## calc dr 
        subsd nb311_ixH2(%rsp),%xmm3
        subsd nb311_iyH2(%rsp),%xmm4
        subsd nb311_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb311_dxH2(%rsp)
        movapd %xmm4,nb311_dyH2(%rsp)
        movapd %xmm5,nb311_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb311_rinvO(%rsp)         ## rinvO in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb311_rO(%rsp)    ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb311_rinvH1(%rsp)         ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb311_rH1(%rsp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb311_rinvH2(%rsp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb311_rH2(%rsp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb311_rinvO(%rsp),%xmm0
        mulsd   nb311_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb311_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7  
    mulsd  %xmm1,%xmm6      ## xmm6=Geps 
    mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
    addsd  %xmm6,%xmm5
    addsd  %xmm7,%xmm5      ## xmm5=Fp        
    mulsd  nb311_two(%rsp),%xmm7         ## two*Heps2 
    movapd nb311_qqO(%rsp),%xmm0
    addsd  %xmm6,%xmm7
    addsd  %xmm5,%xmm7 ## xmm7=FF 
    mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addsd  %xmm4,%xmm5 ## xmm5=VV 
    mulsd  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulsd  %xmm7,%xmm0 ## fijC=FF*qq 

        ## do nontable L-J 
        movapd nb311_rinvO(%rsp),%xmm2
        mulsd  %xmm2,%xmm2

    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addsd  nb311_vctot(%rsp),%xmm5
    movlpd %xmm5,nb311_vctot(%rsp)

        movapd %xmm2,%xmm1
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb311_c6(%rsp),%xmm1
        mulsd  nb311_c12(%rsp),%xmm4
        movapd %xmm4,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulsd  nb311_six(%rsp),%xmm1
        mulsd  nb311_twelve(%rsp),%xmm4
        subsd  %xmm1,%xmm4
        addsd  nb311_Vvdwtot(%rsp),%xmm3
        mulsd  nb311_rinvO(%rsp),%xmm4
        mulsd  nb311_tsc(%rsp),%xmm0
        subsd  %xmm0,%xmm4
        movlpd %xmm3,nb311_Vvdwtot(%rsp)
        mulsd  nb311_rinvO(%rsp),%xmm4

        movapd nb311_dxO(%rsp),%xmm0
        movapd nb311_dyO(%rsp),%xmm1
        movapd nb311_dzO(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2      ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb311_fixO(%rsp),%xmm3
        movapd nb311_fiyO(%rsp),%xmm4
        movapd nb311_fizO(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb311_fixO(%rsp)
        movlpd %xmm4,nb311_fiyO(%rsp)
        movlpd %xmm7,nb311_fizO(%rsp)
        ## update j forces with water O 
        movlpd %xmm0,nb311_fjx(%rsp)
        movlpd %xmm1,nb311_fjy(%rsp)
        movlpd %xmm2,nb311_fjz(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb311_rH1(%rsp),%xmm7
        mulpd nb311_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb311_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb311_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb311_qqH(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb311_vctot(%rsp),%xmm5
        mulsd  nb311_rinvH1(%rsp),%xmm3
    movlpd %xmm5,nb311_vctot(%rsp)
        mulsd  nb311_tsc(%rsp),%xmm3
        subsd %xmm3,%xmm4

        movapd nb311_dxH1(%rsp),%xmm0
        movapd nb311_dyH1(%rsp),%xmm1
        movapd nb311_dzH1(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb311_fixH1(%rsp),%xmm3
        movapd nb311_fiyH1(%rsp),%xmm4
        movapd nb311_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb311_fixH1(%rsp)
        movlpd %xmm4,nb311_fiyH1(%rsp)
        movlpd %xmm7,nb311_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb311_fjx(%rsp),%xmm0
        addsd  nb311_fjy(%rsp),%xmm1
        addsd  nb311_fjz(%rsp),%xmm2
        movlpd %xmm0,nb311_fjx(%rsp)
        movlpd %xmm1,nb311_fjy(%rsp)
        movlpd %xmm2,nb311_fjz(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb311_rH2(%rsp),%xmm7
        mulsd   nb311_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb311_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb311_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb311_qqH(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb311_vctot(%rsp),%xmm5
        mulsd  nb311_rinvH2(%rsp),%xmm3
    movlpd %xmm5,nb311_vctot(%rsp)
        mulsd  nb311_tsc(%rsp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb311_dxH2(%rsp),%xmm0
        movapd nb311_dyH2(%rsp),%xmm1
        movapd nb311_dzH2(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb311_fixH2(%rsp),%xmm3
        movapd nb311_fiyH2(%rsp),%xmm4
        movapd nb311_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb311_fixH2(%rsp)
        movlpd %xmm4,nb311_fiyH2(%rsp)
        movlpd %xmm7,nb311_fizH2(%rsp)

        movq nb311_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb311_fjx(%rsp),%xmm0
        addsd  nb311_fjy(%rsp),%xmm1
        addsd  nb311_fjz(%rsp),%xmm2

        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%rax,8),%xmm0
        addsd 8(%rdi,%rax,8),%xmm1
        addsd 16(%rdi,%rax,8),%xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)

_nb_kernel311_x86_64_sse2.nb311_updateouterdata: 
        movl  nb311_ii3(%rsp),%ecx
        movq  nb311_faction(%rbp),%rdi
        movq  nb311_fshift(%rbp),%rsi
        movl  nb311_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb311_fixO(%rsp),%xmm0
        movapd nb311_fiyO(%rsp),%xmm1
        movapd nb311_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb311_fixH1(%rsp),%xmm0
        movapd nb311_fiyH1(%rsp),%xmm1
        movapd nb311_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%rdi,%rcx,8),%xmm3
        movsd  32(%rdi,%rcx,8),%xmm4
        movsd  40(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,24(%rdi,%rcx,8)
        movsd  %xmm4,32(%rdi,%rcx,8)
        movsd  %xmm5,40(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb311_fixH2(%rsp),%xmm0
        movapd nb311_fiyH2(%rsp),%xmm1
        movapd nb311_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  48(%rdi,%rcx,8),%xmm3
        movsd  56(%rdi,%rcx,8),%xmm4
        movsd  64(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,48(%rdi,%rcx,8)
        movsd  %xmm4,56(%rdi,%rcx,8)
        movsd  %xmm5,64(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%rsi,%rdx,8),%xmm3
        movhpd 8(%rsi,%rdx,8),%xmm3
        movsd  16(%rsi,%rdx,8),%xmm4
        subpd  %xmm6,%xmm3
        subsd  %xmm7,%xmm4
        movlpd %xmm3,(%rsi,%rdx,8)
        movhpd %xmm3,8(%rsi,%rdx,8)
        movsd  %xmm4,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb311_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb311_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb311_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb311_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb311_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb311_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb311_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311_x86_64_sse2.nb311_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311_n(%rsp)
        jmp _nb_kernel311_x86_64_sse2.nb311_outer
_nb_kernel311_x86_64_sse2.nb311_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311_x86_64_sse2.nb311_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311_x86_64_sse2.nb311_threadloop
_nb_kernel311_x86_64_sse2.nb311_end: 
        movl nb311_nouter(%rsp),%eax
        movl nb311_ninner(%rsp),%ebx
        movq nb311_outeriter(%rbp),%rcx
        movq nb311_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $904,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel311nf_x86_64_sse2
.globl _nb_kernel311nf_x86_64_sse2
nb_kernel311nf_x86_64_sse2:     
_nb_kernel311nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb311nf_fshift, 16
.set nb311nf_gid, 24
.set nb311nf_pos, 32
.set nb311nf_faction, 40
.set nb311nf_charge, 48
.set nb311nf_p_facel, 56
.set nb311nf_argkrf, 64
.set nb311nf_argcrf, 72
.set nb311nf_Vc, 80
.set nb311nf_type, 88
.set nb311nf_p_ntype, 96
.set nb311nf_vdwparam, 104
.set nb311nf_Vvdw, 112
.set nb311nf_p_tabscale, 120
.set nb311nf_VFtab, 128
.set nb311nf_invsqrta, 136
.set nb311nf_dvda, 144
.set nb311nf_p_gbtabscale, 152
.set nb311nf_GBtab, 160
.set nb311nf_p_nthreads, 168
.set nb311nf_count, 176
.set nb311nf_mtx, 184
.set nb311nf_outeriter, 192
.set nb311nf_inneriter, 200
.set nb311nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb311nf_ixO, 0
.set nb311nf_iyO, 16
.set nb311nf_izO, 32
.set nb311nf_ixH1, 48
.set nb311nf_iyH1, 64
.set nb311nf_izH1, 80
.set nb311nf_ixH2, 96
.set nb311nf_iyH2, 112
.set nb311nf_izH2, 128
.set nb311nf_iqO, 144
.set nb311nf_iqH, 160
.set nb311nf_qqO, 176
.set nb311nf_qqH, 192
.set nb311nf_rinvO, 208
.set nb311nf_rinvH1, 224
.set nb311nf_rinvH2, 240
.set nb311nf_rO, 256
.set nb311nf_rH1, 272
.set nb311nf_rH2, 288
.set nb311nf_tsc, 304
.set nb311nf_c6, 320
.set nb311nf_c12, 336
.set nb311nf_vctot, 352
.set nb311nf_Vvdwtot, 368
.set nb311nf_half, 384
.set nb311nf_three, 400
.set nb311nf_is3, 416
.set nb311nf_ii3, 420
.set nb311nf_nri, 424
.set nb311nf_iinr, 432
.set nb311nf_jindex, 440
.set nb311nf_jjnr, 448
.set nb311nf_shift, 456
.set nb311nf_shiftvec, 464
.set nb311nf_facel, 472
.set nb311nf_innerjjnr, 480
.set nb311nf_ntia, 488
.set nb311nf_innerk, 492
.set nb311nf_n, 496
.set nb311nf_nn1, 500
.set nb311nf_nouter, 504
.set nb311nf_ninner, 508
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $520,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb311nf_nouter(%rsp)
        movl %eax,nb311nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb311nf_nri(%rsp)
        movq %rsi,nb311nf_iinr(%rsp)
        movq %rdx,nb311nf_jindex(%rsp)
        movq %rcx,nb311nf_jjnr(%rsp)
        movq %r8,nb311nf_shift(%rsp)
        movq %r9,nb311nf_shiftvec(%rsp)
        movq nb311nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb311nf_facel(%rsp)

        movq nb311nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb311nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb311nf_half(%rsp)
        movl %ebx,nb311nf_half+4(%rsp)
        movsd nb311nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb311nf_half(%rsp)
        movapd %xmm3,nb311nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb311nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb311nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4
        movq nb311nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb311nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb311nf_iqO(%rsp)
        movapd %xmm4,nb311nf_iqH(%rsp)

        movq  nb311nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb311nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb311nf_ntia(%rsp)
_nb_kernel311nf_x86_64_sse2.nb311nf_threadloop: 
        movq  nb311nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel311nf_x86_64_sse2.nb311nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311nf_x86_64_sse2.nb311nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311nf_n(%rsp)
        movl %ebx,nb311nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311nf_x86_64_sse2.nb311nf_outerstart
        jmp _nb_kernel311nf_x86_64_sse2.nb311nf_end

_nb_kernel311nf_x86_64_sse2.nb311nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311nf_nouter(%rsp),%ebx
        movl %ebx,nb311nf_nouter(%rsp)

_nb_kernel311nf_x86_64_sse2.nb311nf_outer: 
        movq  nb311nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb311nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb311nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb311nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb311nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb311nf_ixO(%rsp)
        movapd %xmm4,nb311nf_iyO(%rsp)
        movapd %xmm5,nb311nf_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb311nf_ixH1(%rsp)
        movapd %xmm1,nb311nf_iyH1(%rsp)
        movapd %xmm2,nb311nf_izH1(%rsp)
        movapd %xmm3,nb311nf_ixH2(%rsp)
        movapd %xmm4,nb311nf_iyH2(%rsp)
        movapd %xmm5,nb311nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb311nf_vctot(%rsp)
        movapd %xmm4,nb311nf_Vvdwtot(%rsp)

        movq  nb311nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb311nf_pos(%rbp),%rsi
        movq  nb311nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb311nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb311nf_ninner(%rsp),%ecx
        movl  %ecx,nb311nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb311nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel311nf_x86_64_sse2.nb311nf_unroll_loop
        jmp   _nb_kernel311nf_x86_64_sse2.nb311nf_checksingle
_nb_kernel311nf_x86_64_sse2.nb311nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb311nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb311nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb311nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311nf_iqO(%rsp),%xmm3
        mulpd  nb311nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb311nf_qqO(%rsp)
        movapd  %xmm4,nb311nf_qqH(%rsp)

        movq nb311nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb311nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb311nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movlpd (%rsi,%rbx,8),%xmm7      ## c6b
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        movhpd 8(%rsi,%rbx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb311nf_c6(%rsp)
        movapd %xmm6,nb311nf_c12(%rsp)

        movq nb311nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb311nf_ixO(%rsp),%xmm4
        movapd nb311nf_iyO(%rsp),%xmm5
        movapd nb311nf_izO(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb311nf_ixH1(%rsp),%xmm4
        movapd nb311nf_iyH1(%rsp),%xmm5
        movapd nb311nf_izH1(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb311nf_ixH2(%rsp),%xmm3
        movapd nb311nf_iyH2(%rsp),%xmm4
        movapd nb311nf_izH2(%rsp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb311nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb311nf_rinvO(%rsp)       ## rinvO in xmm4 
        mulpd   %xmm4,%xmm7
        movapd  %xmm7,nb311nf_rO(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb311nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH1(%rsp)       ## rinvH1 
        mulpd  %xmm4,%xmm6
        movapd %xmm6,nb311nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb311nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb311nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb311nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH2(%rsp)   ## rinv 
        mulpd %xmm4,%xmm5
        movapd %xmm5,nb311nf_rH2(%rsp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb311nf_rinvO(%rsp),%xmm0
        mulpd   nb311nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movq nb311nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7  
    mulpd  %xmm1,%xmm6      ## xmm6=Geps 
    mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
    addpd  %xmm6,%xmm5
    addpd  %xmm7,%xmm5      ## xmm5=Fp        
    movapd nb311nf_qqO(%rsp),%xmm0
    mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addpd  %xmm4,%xmm5 ## xmm5=VV 
    mulpd  %xmm0,%xmm5 ## vcoul=qq*VV  

        ## do nontable L-J 
        movapd nb311nf_rinvO(%rsp),%xmm2
        mulpd  %xmm2,%xmm2

    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addpd  nb311nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb311nf_vctot(%rsp)

        movapd %xmm2,%xmm1
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulpd  nb311nf_c6(%rsp),%xmm1
        mulpd  nb311nf_c12(%rsp),%xmm4
        movapd %xmm4,%xmm3
        subpd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addpd  nb311nf_Vvdwtot(%rsp),%xmm3
        movapd %xmm3,nb311nf_Vvdwtot(%rsp)


        ## Done with O interactions - now H1! 
        movapd nb311nf_rH1(%rsp),%xmm7
        mulpd nb311nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb311nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb311nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb311nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb311nf_vctot(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb311nf_rH2(%rsp),%xmm7
        mulpd   nb311nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb311nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb311nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb311nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb311nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb311nf_innerk(%rsp)
        jl    _nb_kernel311nf_x86_64_sse2.nb311nf_checksingle
        jmp   _nb_kernel311nf_x86_64_sse2.nb311nf_unroll_loop
_nb_kernel311nf_x86_64_sse2.nb311nf_checksingle: 
        movl  nb311nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel311nf_x86_64_sse2.nb311nf_dosingle
        jmp   _nb_kernel311nf_x86_64_sse2.nb311nf_updateouterdata
_nb_kernel311nf_x86_64_sse2.nb311nf_dosingle: 
        movq  nb311nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb311nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb311nf_iqO(%rsp),%xmm3
        mulpd  nb311nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb311nf_qqO(%rsp)
        movapd  %xmm4,nb311nf_qqH(%rsp)

        movq nb311nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb311nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb311nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb311nf_c6(%rsp)
        movapd %xmm6,nb311nf_c12(%rsp)

        movq nb311nf_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb311nf_ixO(%rsp),%xmm4
        movapd nb311nf_iyO(%rsp),%xmm5
        movapd nb311nf_izO(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb311nf_ixH1(%rsp),%xmm4
        movapd nb311nf_iyH1(%rsp),%xmm5
        movapd nb311nf_izH1(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb311nf_ixH2(%rsp),%xmm3
        movapd nb311nf_iyH2(%rsp),%xmm4
        movapd nb311nf_izH2(%rsp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb311nf_rinvO(%rsp)       ## rinvO in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb311nf_rO(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH1(%rsp)       ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb311nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb311nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb311nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb311nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb311nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb311nf_rinvH2(%rsp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb311nf_rH2(%rsp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        movapd nb311nf_rinvO(%rsp),%xmm0
        mulsd   nb311nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb311nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7  
    mulsd  %xmm1,%xmm6      ## xmm6=Geps 
    mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
    addsd  %xmm6,%xmm5
    addsd  %xmm7,%xmm5      ## xmm5=Fp        
    movapd nb311nf_qqO(%rsp),%xmm0
    mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addsd  %xmm4,%xmm5 ## xmm5=VV 
    mulsd  %xmm0,%xmm5 ## vcoul=qq*VV  

        ## do nontable L-J 
        movapd nb311nf_rinvO(%rsp),%xmm2
        mulsd  %xmm2,%xmm2

    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    addsd  nb311nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb311nf_vctot(%rsp)

        movapd %xmm2,%xmm1
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb311nf_c6(%rsp),%xmm1
        mulsd  nb311nf_c12(%rsp),%xmm4
        movapd %xmm4,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addsd  nb311nf_Vvdwtot(%rsp),%xmm3
        movlpd %xmm3,nb311nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb311nf_rH1(%rsp),%xmm7
        mulpd nb311nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb311nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb311nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addsd  nb311nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb311nf_vctot(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb311nf_rH2(%rsp),%xmm7
        mulsd   nb311nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb311nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb311nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addsd  nb311nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb311nf_vctot(%rsp)

_nb_kernel311nf_x86_64_sse2.nb311nf_updateouterdata: 
        ## get n from stack
        movl nb311nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb311nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb311nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb311nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb311nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb311nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb311nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311nf_x86_64_sse2.nb311nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311nf_n(%rsp)
        jmp _nb_kernel311nf_x86_64_sse2.nb311nf_outer
_nb_kernel311nf_x86_64_sse2.nb311nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311nf_x86_64_sse2.nb311nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311nf_x86_64_sse2.nb311nf_threadloop
_nb_kernel311nf_x86_64_sse2.nb311nf_end: 
        movl nb311nf_nouter(%rsp),%eax
        movl nb311nf_ninner(%rsp),%ebx
        movq nb311nf_outeriter(%rbp),%rcx
        movq nb311nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $520,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

