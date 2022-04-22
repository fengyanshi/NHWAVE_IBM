!-----------------------------------------------------------------------------------------
!
!   mod_immersed.F
!
!
!   This file is part of NHWAVE.
!
!   Subroutines in this file:
!
!        (1) obsvel
!        (2) obsflag
!        (3) imm_obs
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

! modulize IBM ----------------------------fyshi

# if defined (OBSTACLE)

MODULE IBM_MODULE
    USE GLOBAL
    IMPLICIT NONE
     real(SP) :: obs_u,obs_v,obs_w
     real(SP), dimension(:,:,:), allocatable :: ObsForceX,ObsForceY,ObsForceZ
     real(SP), dimension(:,:,:), allocatable :: set_dist_x,set_dist_y,set_dist_z
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: Mask_Obs,Rx,Ry,Rz,Tmp_obs
     integer, dimension(:,:,:), allocatable :: set_flag,set_flag_old
    SAVE

CONTAINS

    SUBROUTINE OBS_INITIALIZATION

    USE INPUT_Util
    IMPLICIT NONE
    INTEGER :: i,j,k,iglob,jglob,lglob

     ALLOCATE(set_flag(Mloc,Nloc,Kloc),ObsForceX(Mloc,Nloc,Kloc),ObsForceY(Mloc,Nloc,Kloc), &
              ObsForceZ(Mloc,Nloc,Kloc),set_dist_x(Mloc,Nloc,Kloc),set_dist_y(Mloc,Nloc,Kloc), &
              set_dist_z(Mloc,Nloc,Kloc),set_flag_old(Mloc,Nloc,Kloc))

     ALLOCATE(Tmp_obs(Mglob,Nglob,Kglob))

!,Mask_Obs(Mloc,Nloc,Kglob),Rx(Mloc,Nloc,Kglob),  &
!               Ry(Mloc,Nloc,Kglob),Rz(Mloc,Nloc,Kglob))

     ObsForceX=Zero;ObsForceY=Zero;ObsForceZ=Zero
	 set_dist_x=Zero;set_dist_y=Zero;set_dist_z=Zero
	 set_flag=0;set_flag_old=0;

! read obstacle file

! mask ---
     OPEN(7,FILE=TRIM(OBS_MASK_FILE),status='old')
         do k = 1,Kglob
         do j = 1,Nglob
           read(7,*) (Tmp_Obs(i,j,k),i=1,Mglob)
         enddo
         enddo
     CLOSE(7)
# if defined (PARALLEL)
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           lglob = k-Nghost
           set_flag(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
         enddo
# else
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = i-Nghost
           jglob = j-Nghost
           lglob = k-Nghost
           set_flag(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
# endif

     ! collect data into ghost cells 
     call phi_3D_exch(set_flag)


     IF(OBS_Require_Rxyz) THEN

! Rx ---
     OPEN(7,FILE=TRIM(OBS_Rx_FILE),status='old')
         do k = 1,Kglob
         do j = 1,Nglob
           read(7,*) (Tmp_Obs(i,j,k),i=1,Mglob)
         enddo
         enddo
     CLOSE(7)
# if defined (PARALLEL)
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           lglob = k-Nghost
           set_dist_x(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
         enddo
# else
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = i-Nghost
           jglob = j-Nghost
           lglob = k-Nghost
           set_dist_x(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
# endif

     ! collect data into ghost cells 
     call phi_3D_exch(set_dist_x)


! Ry ---------
     OPEN(7,FILE=TRIM(OBS_Ry_FILE),status='old')
         do k = 1,Kglob
         do j = 1,Nglob
           read(7,*) (Tmp_Obs(i,j,k),i=1,Mglob)
         enddo
         enddo
     CLOSE(7)
# if defined (PARALLEL)
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           lglob = k-Nghost
           set_dist_y(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
         enddo
# else
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = i-Nghost
           jglob = j-Nghost
           lglob = k-Nghost
           set_dist_y(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
# endif

     ! collect data into ghost cells 
     call phi_3D_exch(set_dist_y)

! Rz -----------------------
     OPEN(7,FILE=TRIM(OBS_Ry_FILE),status='old')
         do k = 1,Kglob
         do j = 1,Nglob
           read(7,*) (Tmp_Obs(i,j,k),i=1,Mglob)
         enddo
         enddo
     CLOSE(7)
# if defined (PARALLEL)
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           lglob = k-Nghost
           set_dist_z(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
         enddo
# else
         do k = Kbeg,Kend
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = i-Nghost
           jglob = j-Nghost
           lglob = k-Nghost
           set_dist_z(i,j,k) = Tmp_obs(iglob,jglob,lglob)
         enddo
         enddo
# endif

     ! collect data into ghost cells 
     call phi_3D_exch(set_dist_z)

    ELSE

!    if there are no Rx Ry and Rz files, set zeros

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend

      if(set_flag(i,j,k)==0) then
          ! obstacle on the right
          if(set_flag(i+1,j,k)==1) then
            set_dist_x(i,j,k) = zero
          endif

          ! obstacle on the left
          if(set_flag(i-1,j,k)==1) then
            set_dist_x(i,j,k) = zero
          endif

          ! obstacle on the front
          if(set_flag(i,j-1,k)==1) then
            set_dist_y(i,j,k) = zero
          endif
		  
          ! obstacle on the back                                                                                               
          if(set_flag(i,j+1,k)==1) then
            set_dist_y(i,j,k) = zero
          endif
		  
		  ! obstacle under the fluid
          if(set_flag(i,j,k-1)==1) then 
            set_dist_z(i,j,k) = zero
          endif
      endif

    enddo
    enddo
    enddo
      
    ENDIF

!  Gangfeng set void value as large number, to keep consistence I re-assign my 
!  negative values to large numbers

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
    
    IF(set_dist_x(i,j,k)<0.0) set_dist_x(i,j,k) = 1.e+20
    IF(set_dist_y(i,j,k)<0.0) set_dist_y(i,j,k) = 1.e+20
    IF(set_dist_z(i,j,k)<0.0) set_dist_z(i,j,k) = 1.e+20
    
    enddo
    enddo
    enddo


    END SUBROUTINE OBS_INITIALIZATION

    subroutine set_obsvel
!---------------------------------------------------------------
!   Specify or calculate obstacle velocities
!   by Gangfeng Ma, 17/08/2013
!---------------------------------------------------------------
!    use global, only: Zero,obs_u,obs_v,obs_w
    implicit none

    obs_u = Zero
    obs_v = Zero
    obs_w = Zero

    return
    end subroutine set_obsvel

    subroutine imm_obs
!---------------------------------------------------------------------
!   Immersed-Boundary Method
!   Calculate forcing at the obstacle boundary
!   By Gangfeng Ma, 17/08/2013
!--------------------------------------------------------------------
!    use global
    implicit none
    real(SP),dimension(:),allocatable :: target_vel,target_dist
    integer :: i,j,k,l,count
    real(SP) :: numer,denom,force_vel,dista,vel_grad

    ! calculate forcing at boundaries
    ObsForceX = Zero
    ObsForceY = Zero
    ObsForceZ = Zero

    allocate(target_vel(1:6))
    allocate(target_dist(1:6))

    ! x-direction
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for fluid cells
      if(set_flag(i,j,k)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j,k) is fluid, (i-1,j,k) is obstacle
        if(set_flag(i-1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i+1,j,k) to obstacle/fluid interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (U(i+1,j,k)-obs_u)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i+1,j,k) is obstacle
        if(set_flag(i+1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) 
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i-1,j,k) to obstacle/fluid interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_u-U(i-1,j,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_u-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j-1,k) is obstacle
        if(Nloc>=2.and.set_flag(i,j-1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j+1,k) to obstacle/fluid interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (U(i,j+1,k)-obs_u)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j+1,k) is obstacle 
        if((Nloc>=2).and.set_flag(i,j+1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j-1,k) to obstacle/fluid interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_u-U(i,j-1,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(yc(j)-yc(j-1))
          target_vel(count) = obs_u-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k-1) is obstacle
        if(set_flag(i,j,k-1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) -- (adjust according to obstacle)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k+1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k+1)-sigc(k))
          vel_grad = (U(i,j,k+1)-obs_u)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k+1)-sigc(k))
          target_vel(count) = obs_u+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k+1) is obstacle 
        if(set_flag(i,j,k+1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k-1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k)-sigc(k-1))
          vel_grad = (obs_u-U(i,j,k-1))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k)-sigc(k-1))
          target_vel(count) = obs_u-vel_grad*dista
          target_dist(count) = dista
        endif

        ! Search for target distances that are zero, to avoid divide-by-zero problems.
        if(count>=1) then
          numer = 0.0
          denom = 0.0
          do l = 1,count
            if (target_dist(l)<1.e-16) then
              force_vel = target_vel(l)
              goto 1
            else
              numer = numer+target_vel(l)*1.0/target_dist(l)
              denom = denom+1.0/target_dist(l)
            endif
          enddo
          force_vel = numer/denom
 1        continue
          ObsForceX(i,j,k) = (D(i,j)*force_vel-DU(i,j,k))/dt
        endif
      endif
    enddo
    enddo
    enddo

    ! y-direction
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for fluid cells
      if(set_flag(i,j,k)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j,k) is fluid, (i-1,j,k) is obstacle
        if(set_flag(i-1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i+1,j,k) to obstacle/fluid interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (V(i+1,j,k)-obs_v)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i+1,j,k) is obstacle
        if(set_flag(i+1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) 
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i-1,j,k) to obstacle/fluid interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_v-V(i-1,j,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)                                                                              
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_v-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j-1,k) is obstacle
        if(Nloc>=2.and.set_flag(i,j-1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j+1,k) to obstacle/fluid interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (V(i,j+1,k)-obs_v)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j+1,k) is obstacle 
        if((Nloc>=2).and.set_flag(i,j+1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j-1,k) to obstacle/fluid interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_v-V(i,j-1,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(yc(j)-yc(j-1))
          target_vel(count) = obs_v-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k-1) is obstacle
        if(set_flag(i,j,k-1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) -- (adjust according to obstacle)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k+1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k+1)-sigc(k))
          vel_grad = (V(i,j,k+1)-obs_v)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k+1)-sigc(k))
          target_vel(count) = obs_v+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k+1) is obstacle 
        if(set_flag(i,j,k+1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k-1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k)-sigc(k-1))
          vel_grad = (obs_v-V(i,j,k-1))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k)-sigc(k-1))
          target_vel(count) = obs_v-vel_grad*dista
          target_dist(count) = dista
        endif

        ! Search for target distances that are zero, to avoid divide-by-zero problems.
        if(count>=1) then
          numer = 0.0
          denom = 0.0
          do l = 1,count
            if (target_dist(l)<1.e-16) then
              force_vel = target_vel(l)
              goto 2
            else
              numer = numer+target_vel(l)*1.0/target_dist(l)
              denom = denom+1.0/target_dist(l)
            endif
          enddo
          force_vel = numer/denom
 2        continue
          ObsForceY(i,j,k) = (D(i,j)*force_vel-DV(i,j,k))/dt
        endif
      endif
    enddo
    enddo
    enddo
 
    ! z-direction
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! for fluid cells
      if(set_flag(i,j,k)==0) then
        count = 0
        target_vel = 0.0
        target_dist = 0.0

        ! i.e. (i,j,k) is fluid, (i-1,j,k) is obstacle
        if(set_flag(i-1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i+1,j,k) to obstacle/fluid interface
          dista = dista+(xc(i+1)-xc(i))
          vel_grad = (W(i+1,j,k)-obs_w)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(xc(i+1)-xc(i))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i+1,j,k) is obstacle
        if(set_flag(i+1,j,k)==1.and.set_dist_x(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) 
          dista = set_dist_x(i,j,k)
          ! compute velocity gradient from (i-1,j,k) to obstacle/fluid interface  
          dista = dista+(xc(i)-xc(i-1))
          vel_grad = (obs_w-W(i-1,j,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)                                                                              
          dista = dista-(xc(i)-xc(i-1))
          target_vel(count) = obs_w-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j-1,k) is obstacle
        if(Nloc>=2.and.set_flag(i,j-1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j+1,k) to obstacle/fluid interface 
          dista = dista+(yc(j+1)-yc(j))
          vel_grad = (W(i,j+1,k)-obs_w)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k) 
          dista = dista-(yc(j+1)-yc(j))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j+1,k) is obstacle 
        if((Nloc>=2).and.set_flag(i,j+1,k)==1.and.set_dist_y(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_y(i,j,k)
          ! compute velocity gradient from (i,j-1,k) to obstacle/fluid interface 
          dista = dista+(yc(j)-yc(j-1))
          vel_grad = (obs_w-W(i,j-1,k))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-(yc(j)-yc(j-1))
          target_vel(count) = obs_w-vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k-1) is obstacle
        if(set_flag(i,j,k-1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k) -- (adjust according to obstacle)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k+1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k+1)-sigc(k))
          vel_grad = (W(i,j,k+1)-obs_w)/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k+1)-sigc(k))
          target_vel(count) = obs_w+vel_grad*dista
          target_dist(count) = dista
        endif

        ! i.e. (i,j,k) is fluid, (i,j,k+1) is obstacle 
        if(set_flag(i,j,k+1)==1.and.set_dist_z(i,j,k)<1.e+10) then
          count = count+1
          ! Determine distance of interface to point (i,j,k)
          dista = set_dist_z(i,j,k)
          ! compute velocity gradient from (i,j,k-1) to obstacle/fluid interface 
          dista = dista+D(i,j)*(sigc(k)-sigc(k-1))
          vel_grad = (obs_w-W(i,j,k-1))/dista
          ! interpolate velocity from obstacle/fluid interface to (i,j,k)
          dista = dista-D(i,j)*(sigc(k)-sigc(k-1))
          target_vel(count) = obs_w-vel_grad*dista
          target_dist(count) = dista
        endif

        ! Search for target distances that are zero, to avoid divide-by-zero problems.
        if(count>=1) then
          numer = 0.0
          denom = 0.0
          do l = 1,count
            if (target_dist(l)<1.e-16) then
              force_vel = target_vel(l)
              goto 3
            else
              numer = numer+target_vel(l)*1.0/target_dist(l)
              denom = denom+1.0/target_dist(l)
            endif
          enddo
          force_vel = numer/denom
 3        continue
          ObsForceZ(i,j,k) = (D(i,j)*force_vel-DW(i,j,k))/dt
        endif
      endif
    enddo
    enddo
    enddo

    deallocate(target_vel)
    deallocate(target_dist)
 1000 continue

    return
    end subroutine imm_obs

END MODULE IBM_MODULE

# endif


