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
                       DelxSl,DelySl,Kloc1, &
                       Ibeg,Iend,Jbeg,Jend,&
                       Iend1,Jend1,&
                       Grd_R, Ivgrd, Kbeg,Kend,Kend1,Kglob,Nghost,dsig,sigc,sig, &
                       u_bar,v_bar,GAMMA_s,U,V,W,DX,DY,Omega,MASK,SP,HC,SMALL, &
                       Eta, DelxS,DelyS,DelzS, DeltS, DeltSc,z_face,z_face0, &
                       z_center, z_center0,DT
                       


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
    REAL :: Sigma_original(Kloc1),Dsig_original(Kloc)



    ! s-function, 1) kirby double-power
    ! s = c_upper * sigma^A_upper + c_lower *sigma^A_lower
    !  c_upper = 1, c_lower = 0, A_upper = 0.33

    REAL :: C_upper, C_lower, A_upper, A_lower
      C_upper = 1.0
      C_lower = ZERO
      A_upper = 0.5
      A_lower = ZERO

     do k = 1,Kloc
     Dsig_original(k) = 1.0/float(Kglob)
     enddo

     Sigma_original(Kbeg) = Zero
     do k = Kbeg+1,Kend1
     Sigma_original(k) = Sigma_original(k-1)+Dsig_original(k-1)
     enddo

!     do k = Kbeg-1,1,-1
!     Sigma_original(k) = Sigma_original(k+1)-Dsig_original(k)
!     enddo

!   new sigma (s) 

     DO K=Kbeg,Kend1
       SIG(K) = C_upper * Sigma_original(K)**(A_upper) + C_lower * Sigma_original(K)**(A_lower)
     ENDDO

!  make constant for ghost cells
     DO K=Kend1+1,Kloc1
        SIG(K)=SIG(Kend1)+(K-Kend1)*(SIG(Kend1)-SIG(Kend1-1))
     ENDDO     
     DO K=Kbeg-1,1,-1
        SIG(K)=SIG(Kbeg)+(K-Kbeg)*(SIG(Kbeg+1)-SIG(Kbeg))
     ENDDO

     DO K=1,Kloc
        DSIG(K) = SIG(K+1)-SIG(K)
     ENDDO

     DO K=1,Kloc
       SIGC(K)= 0.5*(SIG(K)+SIG(K+1))
     ENDDO

!     print*, SIG
!     print*, SIGC
!     print*,DSIG
!    STOP

!
!    vertical grid, this part was from initialize.F subroutine generate_grid
!

!    D=10.0
!    Hc=10.0

! initialize z_face and z_center

     DO J=1,Nloc
     DO I=1,Mloc

   
        DO K=1,Kloc1
          z_face(I,J,K)=-Hc(I,J)+Hc(I,J)*SIG(K)
        ENDDO

        DO K=1,Kloc
          z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
        ENDDO

     ENDDO
     ENDDO


!     print*,Mloc,Nloc,Kloc,Nghost
!     print*,Hc(50,4),D(50,4)
!     print*, z_face(50,4,1),z_face(50,4,2),z_face(50,4,3),z_face(50,4,4),z_face(50,4,5),z_face(50,4,6), &
!             z_face(50,4,7),z_face(50,4,8),z_face(50,4,9),z_face(50,4,10),z_face(50,4,11)

!    STOP


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

!     REAL(SP),DIMENSION(Mloc,Nloc,Kloc) :: z_center
!     REAL(SP),DIMENSION(Mloc,Nloc,Kloc1) :: z_face

     z_face0=z_face
     z_center0 = z_center


     DO J=1,Nloc
     DO I=1,Mloc

   
        DO K=1,Kloc1
          z_face(I,J,K)=-Hc(I,J)+D(I,J)*SIG(K)
        ENDDO

        DO K=1,Kloc
          z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
        ENDDO

     ENDDO
     ENDDO


     DO J=1,Nloc
     DO I=1,Mloc
        DO K=Kbeg,Kend1
          z_face(I,J,K)=-Hc(I,J)+(K-Kbeg)*D(I,J)/(Kend1-Kbeg)
        ENDDO

   ! ghost cells
        DO K=1,Kbeg-1
          z_face(I,J,K) = z_face(I,K,Kbeg) - D(I,J)/(Kend1-Kbeg)*(Kbeg-K)
        ENDDO
        DO K=Kend1+1,Kloc1
          z_face(I,J,K) = z_face(I,K,Kend) + D(I,J)/(Kend1-Kbeg)*(K-Kend1)
        ENDDO

        DO K=1,Kloc
          z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
        ENDDO
     ENDDO
     ENDDO

! center points

     DO K=1,Kloc
     DO J=1,Nloc
     Do I=1,Mloc
       DelxSc(I,J,K)=-(z_center(I+1,J,K)-z_center(I-1,J,K))/dx/2.0_SP &
                     /(z_face(I,J,K+1)-z_face(I,J,K))*dsig(K)
       DelySc(I,J,K)=-(z_center(I,J+1,K)-z_center(I,J-1,K))/dy/2.0_SP &
                     /(z_face(I,J,K+1)-z_face(I,J,K))*dsig(K)
       DeltSc(I,J,K) = - (z_center(I,J,K)-z_center0(I,J,K))/DT &
                     /(z_face(I,J,K+1)-z_face(I,J,K))*dsig(K)
     ENDDO
     ENDDO
     ENDDO

! face points

     DO K=Kbeg,Kend1
     DO J=Jbeg,Jend
     Do I=Ibeg,Iend
       DelxS(I,J,K)=-(z_face(I+1,J,K)-z_face(I-1,J,K))/DX/2.0_SP &
                     /(z_center(I,J,K)-z_center(I,J,K-1))*(dsig(K-1)+dsig(K))*0.5_SP
       DelyS(I,J,K)=-(z_face(I,J+1,K)-z_face(I,J-1,K))/DY/2.0_SP &
                     /(z_center(I,J,K)-z_center(I,J,K-1))*(dsig(K-1)+dsig(K))*0.5_SP
       DelzS(I,J,K) = (dsig(K-1)+dsig(K))*0.5_SP/(z_center(I,J,K)-z_center(I,J,K-1))

       DeltS(I,J,K) = - (z_face(I,J,K)-z_face0(I,J,K))/DT &
                     *DelzS(I,J,K)
     ENDDO
     ENDDO
     ENDDO


! old code using sigma coordinate

!     DelxSc = Zero
!     DelySc = Zero

!     do k = 1,Kloc
!     do j = 1,Nloc
!     do i = 1,Mloc
!       DelxSc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelxH(i,j)-sigc(k)/D(i,j)*DelxEta(i,j)
!       DelySc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelyH(i,j)-sigc(k)/D(i,j)*DelyEta(i,j)
!     enddo
!     enddo
!     enddo


     ! added by Morteza, it is at face ponts, in fact, it is duplicated in pressure.F
     ! DelxSl and DelySl are used in diffusion terms implemented by Morteza, 
     ! subroutine kepsilon_3D and subroutine diffusion to deal with breaking waves
     ! Here, I keep the variables in case they still need them, ROXSI
     ! 

     DelxSl = DelxS
     DelySl = DelyS

!     DelxSl = Zero
!     DelySl = Zero
!     do k = 1,Kloc1
!     do j = 1,Nloc
!     do i = 1,Mloc
!       DelxSl(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
!       DelySl(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)
!     enddo
!     enddo
!     enddo

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

! for sigma-coordinate
!     DO K=Kbeg+1,Kend1-1
!     DO I=Ibeg,Iend
!     DO J=Jbeg,Jend    
!       D_x=(D(I+1,J)-D(I-1,J))/DX/2.0_SP
!       D_y=(D(I,J+1)-D(I,J-1))/DY/2.0_SP
!       U_at_s = (U(I,J,K-1)*dsig(K)+U(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
!       V_at_s = (V(I,J,K-1)*dsig(K)+V(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
!       W_at_s = (W(I,J,K-1)*dsig(K)+W(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
!       Omega(I,J,K)= W_at_s - sig(K)*(R1(I,J)+U_at_s*D_x+V_at_s*D_y)
!       if(Mask(i,j)==0) Omega(i,j,k) = zero
!     ENDDO
!     ENDDO
!     ENDDO

! for generic s-coordinate

     DO K=Kbeg+1,Kend1-1
     DO I=Ibeg,Iend
     DO J=Jbeg,Jend    
       U_at_s = (U(I,J,K-1)*dsig(K)+U(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
       V_at_s = (V(I,J,K-1)*dsig(K)+V(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
       W_at_s = (W(I,J,K-1)*dsig(K)+W(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
       Omega(I,J,K)= D(I,J) *(DeltS(I,J,K) + U_at_s*DelxS(I,J,K) &
                                             + V_at_s*DelyS(I,J,K) &
                                             + W_at_s*DelzS(I,K,K) )
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


