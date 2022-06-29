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
                       Eta, DelxS,DelyS,DelzS,DelzSc,DeltS, DeltSc,z_face,z_face0, &
                       z_center, z_center0,DT
# if defined (PARALLEL)
    USE GLOBAL, ONLY: myid,ier
# endif
                       


                !       Eta,DelxEta,DelyEta,DelxU,DelxV,DelyU,DelyV, &
                !       DU,DelxDU,DV,DelxDV,DelyDU,DelyDV,&
                !       EtaxL,EtaxR,EtayL,EtayR,UxL,UxR,VxL,VxR,UyL,UyR,VyL,VyR, &
                !       DUxL,DUxR,DVxL,DVxR,DUyL,DUyR,DVyL,DVyR, &
                !       DxL,DxR,ExL,ExR,DyL,DyR,Ex,Ey

    IMPLICIT NONE
    SAVE
    REAL,DIMENSION(:),ALLOCATABLE :: S_in_sigma,Ds_in_sigma


CONTAINS

    SUBROUTINE COORDINATE_INITIALIZATION

    USE INPUT_Util
    IMPLICIT NONE
    character(len=80) :: FILE_NAME
    integer :: line,ierr
    INTEGER :: i,j,k,iglob,jglob,lglob
    LOGICAL :: EXPONENTIAL, HYPERBOLIC 
    REAL(SP) :: C_upper, C_lower, A_upper, A_lower,HYPER_alpha,tmp1,tmp2

    ALLOCATE(S_in_sigma(Kloc1),Ds_in_sigma(Kloc))

     FILE_NAME='input.txt'

     CALL GET_LOGICAL_VAL(EXPONENTIAL,FILE_NAME,'EXPONENTIAL',line)
     IF(EXPONENTIAL)THEN
       CALL GET_Float_VAL(C_upper,FILE_NAME,'C_upper',line)
       CALL GET_Float_VAL(C_lower,FILE_NAME,'C_lower',line)
       CALL GET_Float_VAL(A_upper,FILE_NAME,'A_upper',line)
       CALL GET_Float_VAL(A_lower,FILE_NAME,'A_lower',line)
     ENDIF
     CALL GET_LOGICAL_VAL(HYPERBOLIC,FILE_NAME,'HYPERBOLIC',line)
     IF(HYPERBOLIC)THEN
       CALL GET_Float_VAL(HYPER_alpha,FILE_NAME,'HYPER_alpha',line)
     ENDIF


     do k = 1,Kloc
     Dsig(k) = 1.0/float(Kglob)
     enddo

     Sig(Kbeg) = Zero
     do k = Kbeg+1,Kend1
     Sig(k) = Sig(k-1)+Dsig(k-1)
     enddo

     do k = Kbeg-1,1,-1
     Sig(k) = Sig(k+1)-Dsig(k)
     enddo

     DO K=1,Kloc
       SIGC(K)= 0.5*(SIG(K)+SIG(K+1))
     ENDDO

   IF(EXPONENTIAL)THEN

!   s in sigma coordinate
    ! s-function, 1) kirby double-power
    ! s = c_upper * sigma^A_upper + c_lower *sigma^A_lower

     DO K=Kbeg,Kend1
       S_in_sigma(K) = C_upper * Sig(K)**(A_upper) + C_lower * Sig(K)**(A_lower)
     ENDDO

   ELSE IF (HYPERBOLIC) THEN

    ! s-function, 2) kirby tanh
    ! s = tanh(alpha*sigma) / tanh(alpha)

     DO K=Kbeg,Kend1
       tmp1=DTANH(HYPER_alpha)
       tmp2=DTANH(HYPER_alpha*Sig(K))
       S_in_sigma(K) = tmp2/tmp1
     ENDDO

   ELSE ! original sigma

     DO K=Kbeg,Kend1
       S_in_sigma(K)=Sig(K)
     ENDDO

   ENDIF

!  make constant for ghost cells
     DO K=Kend1+1,Kloc1
        S_in_sigma(K)=S_in_sigma(Kend1)+(K-Kend1)*(S_in_sigma(Kend1)-S_in_sigma(Kend1-1))
     ENDDO     
     DO K=Kbeg-1,1,-1
        S_in_sigma(K)=S_in_sigma(Kbeg)+(K-Kbeg)*(S_in_sigma(Kbeg+1)-S_in_sigma(Kbeg))
     ENDDO

     DO K=1,Kloc
        Ds_in_sigma(K) = S_in_sigma(K+1)-S_in_sigma(K)
     ENDDO

!     if(myid==0)then
!     do k=1,kloc
!     print*, k-nghost,S_in_sigma(K),DTANH(HYPER_alpha),DTANH(HYPER_alpha*Sig(K))
!     enddo
!     endif
!     stop

!
!    vertical grid, this part was from initialize.F subroutine generate_grid
!


! initialize z_face and z_center

     DO J=1,Nloc
     DO I=1,Mloc

   
        DO K=1,Kloc1
          z_face(I,J,K)=-Hc(I,J)+Hc(I,J)*S_in_sigma(K)
        ENDDO

        DO K=1,Kloc
          z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
        ENDDO

     ENDDO
     ENDDO


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
          z_face(I,J,K)=-Hc(I,J)+D(I,J)*S_in_sigma(K)
        ENDDO

        DO K=1,Kloc
          z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
        ENDDO

     ENDDO
     ENDDO


!  it was not removed last version
!     DO J=1,Nloc
!     DO I=1,Mloc
!        DO K=Kbeg,Kend1
!          z_face(I,J,K)=-Hc(I,J)+(K-Kbeg)*D(I,J)/(Kend1-Kbeg)
!        ENDDO

   ! ghost cells
!        DO K=1,Kbeg-1
!         z_face(I,J,K) = z_face(I,K,Kbeg) - D(I,J)/(Kend1-Kbeg)*(Kbeg-K)
!        ENDDO
!        DO K=Kend1+1,Kloc1
!          z_face(I,J,K) = z_face(I,K,Kend) + D(I,J)/(Kend1-Kbeg)*(K-Kend1)
!        ENDDO

!        DO K=1,Kloc
!          z_center(I,J,K)=0.5*(z_face(I,J,K)+z_face(I,J,K+1))
!        ENDDO
!     ENDDO
!     ENDDO

! center points

     DO K=1,Kloc
     DO J=1,Nloc
     Do I=1,Mloc
       DelxSc(I,J,K)=-(z_center(I+1,J,K)-z_center(I-1,J,K))/dx/2.0_SP &
                     /(z_face(I,J,K+1)-z_face(I,J,K))*dsig(K)

       DelySc(I,J,K)=-(z_center(I,J+1,K)-z_center(I,J-1,K))/dy/2.0_SP &
                     /(z_face(I,J,K+1)-z_face(I,J,K))*dsig(K)

       DelzSc(I,J,K) = dsig(K)/(z_face(I,J,K+1)-z_face(I,J,K))

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

! use Ds_in_sigma for interpolation

     DO K=Kbeg+1,Kend1-1
     DO I=Ibeg,Iend
     DO J=Jbeg,Jend    
!       U_at_s = (U(I,J,K-1)*dsig(K)+U(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
!       V_at_s = (V(I,J,K-1)*dsig(K)+V(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))
!       W_at_s = (W(I,J,K-1)*dsig(K)+W(I,J,K)*dsig(K-1))/(dsig(K)+dsig(K-1))

       U_at_s = (U(I,J,K-1)*Ds_in_sigma(K)+U(I,J,K)*Ds_in_sigma(K-1)) &
                  /(Ds_in_sigma(K)+Ds_in_sigma(K-1))
       V_at_s = (V(I,J,K-1)*Ds_in_sigma(K)+V(I,J,K)*Ds_in_sigma(K-1)) &
                  /(Ds_in_sigma(K)+Ds_in_sigma(K-1))
       W_at_s = (W(I,J,K-1)*Ds_in_sigma(K)+W(I,J,K)*Ds_in_sigma(K-1)) &
                  /(Ds_in_sigma(K)+Ds_in_sigma(K-1))

       Omega(I,J,K)= D(I,J) *(DeltS(I,J,K) + U_at_s*DelxS(I,J,K) &
                                             + V_at_s*DelyS(I,J,K) &
                                             + W_at_s*DelzS(I,K,K) )
       if(Mask(i,j)==0) Omega(i,j,k) = zero
     ENDDO
     ENDDO
     ENDDO

     return

     END SUBROUTINE GET_OMEGA


END MODULE COORDINATE_MODULE


