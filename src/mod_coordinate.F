!-----------------------------------------------------------------------------------------
!
!   mod_immersed.F
!
!
!   This file is part of NHWAVE.
!
!
!------------------------------------------------------------------------------------------
!
!   BSD 2-Clause License
!
!   Copyright (c) 2019, NHWAVE Development Group
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice, this
!     list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!----------------------------------------------------------------------------------------

MODULE COORDINATE_MODULE
    USE GLOBAL, ONLY: Zero,DelxSc,DelySc,D,DelxH,DelyH, &
                       DelxEta,DelyEta,sigc,Mloc,Nloc,Kloc, &
                       DelxSl,DelySl,Kloc1,sig, &
                       Ibeg,Iend,Jbeg,Jend,&
                       Iend1,Jend1,&
                       Grd_R, Ivgrd, Kbeg,Kend,Kend1,Kglob,Nghost,dsig,sigc,sig, &
                       u_bar,v_bar,GAMMA_s,U,V,W,DX,DY,Omega,MASK


                !       Eta,DelxEta,DelyEta,DelxU,DelxV,DelyU,DelyV, &
                !       DU,DelxDU,DV,DelxDV,DelyDU,DelyDV,&
                !       EtaxL,EtaxR,EtayL,EtayR,UxL,UxR,VxL,VxR,UyL,UyR,VyL,VyR, &
                !       DUxL,DUxR,DVxL,DVxR,DUyL,DUyR,DVyL,DVyR, &
                !       DxL,DxR,ExL,ExR,DyL,DyR,Ex,Ey

    IMPLICIT NONE
    SAVE

CONTAINS

    SUBROUTINE COORDINATE_INITIALIZATION

    USE INPUT_Util
    IMPLICIT NONE
    INTEGER :: i,j,k,iglob,jglob,lglob

!
!    vertical grid, this part was from initialize.F subroutine generate_grid
!
     if (Ivgrd==1) then
     do k = 1,Kloc
     dsig(k) = 1.0/float(Kglob)
     enddo
     elseif(Ivgrd==2) then
     dsig(Kbeg) = (Grd_R-1.0)/(Grd_R**float(Kglob)-1.0)
     do k = Kbeg+1,Kend
     dsig(k) = dsig(k-1)*Grd_R
     enddo

     do k = 1,Nghost
     dsig(Kbeg-k) = dsig(Kbeg+k-1)
     enddo

     do k = 1,Nghost
     dsig(Kend+k) = dsig(Kend-k+1)
     enddo
     endif

     sig(Kbeg) = Zero
     do k = Kbeg+1,Kloc1
     sig(k) = sig(k-1)+dsig(k-1)
     sigc(k-1) = sig(k-1)+0.5*dsig(k-1)
     enddo
     do k = Kbeg-1,1,-1
     sig(k) = sig(k+1)-dsig(k)
     sigc(k) = sig(k+1)-0.5*dsig(k)
     enddo

      RETURN
    END SUBROUTINE COORDINATE_INITIALIZATION

     SUBROUTINE SIGMA_TRANSFORM
!--------------------------------------------------- 
!    Calculate sigma transformation coefficient
!    Called by       
!--------------------------------------------------

     IMPLICIT NONE
     INTEGER :: i,j,k

! standard sigma-coordinate has
!      d_sigma/dx = 1/D*(1-sigma)dh/dx - 1/D*sigma d_eta/dx
! for general sigma-coordinate
!      d sigma/dx= -dz/dx / dz/d_sigma
!


     DelxSc = Zero
     DelySc = Zero

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DelxSc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelxH(i,j)-sigc(k)/D(i,j)*DelxEta(i,j)
       DelySc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelyH(i,j)-sigc(k)/D(i,j)*DelyEta(i,j)
     enddo
     enddo
     enddo


     ! added by Morteza
     DelxSl = Zero
     DelySl = Zero
     do k = 1,Kloc1
     do j = 1,Nloc
     do i = 1,Mloc
       DelxSl(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
       DelySl(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)
     enddo
     enddo
     enddo

     END SUBROUTINE SIGMA_TRANSFORM


     SUBROUTINE GET_OMEGA(R1)

!-----------------------------------------------
!    Obtain vertical velocity in sigma-coord.
!    Called by 
!       eval_duvw
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), dimension(Mloc,Nloc), intent(in) :: R1
     REAL(SP) :: D_x,D_y,U_at_s,V_at_s,W_at_s

! use the omega description directly, no TVD scheme here because it is not
! conservative form

     Omega = ZERO
     DO K=Kbeg+1,Kend1-1
     DO I=Ibeg,Iend
     DO J=Jbeg,Jend    
       D_x=(D(I+1,J)-D(I-1,J))/DX/2.0_SP
       D_y=(D(I,J+1)-D(I,J-1))/DY/2.0_SP
       U_at_s = (U(I,J,K-1)*dsig(K)+U(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
       V_at_s = (V(I,J,K-1)*dsig(K)+V(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
       W_at_s = (W(I,J,K-1)*dsig(K)+W(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
       Omega(I,J,K)= W_at_s - sig(K)*(R1(I,J)+U_at_s*D_x+V_at_s*D_y)
       if(Mask(i,j)==0) Omega(i,j,k) = zero
     ENDDO
     ENDDO
     ENDDO

     return

     END SUBROUTINE GET_OMEGA

     SUBROUTINE MEAN_VELOCITY

       IMPLICIT NONE
       INTEGER :: i,j,k

       U_bar = ZERO
       V_bar = ZERO
       DO J = 1,Nloc
       DO I = 1,Mloc
       DO K = Kbeg,Kend
         U_bar(I,J)=U_bar(I,J)+U(I,J,K)*dsig(K)
         V_bar(I,J)=V_bar(I,J)+V(I,J,K)*dsig(K)
       ENDDO
       ENDDO
       ENDDO
         
     END SUBROUTINE MEAN_VELOCITY


     SUBROUTINE GAMMAS
!--------------------------------------------------- 
!    Calculate sigma transformation coefficient
!    Called by       
!--------------------------------------------------

       IMPLICIT NONE
       INTEGER :: i,j,k

       DO K = Ibeg,Iend      
       DO J = Jbeg,Jend
       DO I = Kbeg,Kend 
         GAMMA_s(I,J,K) = (D(I+1,J)*(U(I+1,J,K)-U_bar(I+1,J))-D(I-1,J)*(U(I-1,J,K)-U_bar(I-1,J))) & 
                          /DX/2.0 &
                        + (D(I,J+1)*(V(I,J+1,K)-U_bar(I,J+1))-D(I,J-1)*(U(I,J-1,K)-U_bar(I,J-1))) & 
                          /DY/2.0 &
                        + (Omega(I,J,K+1)-Omega(I,J,K))/dsig(K)/2.0
       ENDDO
       ENDDO
       ENDDO


     END SUBROUTINE GAMMAS



END MODULE COORDINATE_MODULE


