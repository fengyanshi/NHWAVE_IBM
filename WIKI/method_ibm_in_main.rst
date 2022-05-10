Modification in main program
------------------------------

1) add input file names in mod_global.F

  .. code-block:: rest

     logical :: OBS_Require_Rxyz=.true.
     character(len=80) :: OBS_MASK_FILE
     character(len=80) :: OBS_Rx_FILE
     character(len=80) :: OBS_Ry_FILE
     character(len=80) :: OBS_Rz_FILE

2) read input file names and logical parameter in initialize.F

  .. code-block:: rest

       CALL GET_LOGICAL_VAL(OBS_Require_Rxyz,FILE_NAME,'OBS_Require_Rxyz',line)
       CALL GET_STRING_VAL(OBS_MASK_FILE,FILE_NAME,'OBS_MASK_FILE',line,ierr)
       IF(OBS_Require_Rxyz)THEN
         CALL GET_STRING_VAL(OBS_Rx_FILE,FILE_NAME,'OBS_Rx_FILE',line,ierr) 
         CALL GET_STRING_VAL(OBS_Ry_FILE,FILE_NAME,'OBS_Ry_FILE',line,ierr) 
         CALL GET_STRING_VAL(OBS_Rz_FILE,FILE_NAME,'OBS_Rz_FILE',line,ierr)      
       ENDIF

3) make calls in the main program

  .. code-block:: rest

     CALL OBS_INITIALIZATION

  .. code-block:: rest

         ! obstacle velocity 
         call set_obsvel

         ! set obstacle flag                                                                                                   
         ! call set_obsflag ! change to ibm module

  .. code-block:: rest

     subroutine eval_duvw(ISTEP)
     ...

           if(set_flag(i,j,k)==1) then
             DU(i,j,k) = D(i,j)*obs_u
           else
             DU(i,j,k) = DU(i,j,k)+dt*R2(i,j,k)
           endif
     ...
           if(set_flag(i,j,k)==1) then
             DV(i,j,k) = D(i,j)*obs_v
           else
             DV(i,j,k) = DV(i,j,k)+dt*R3(i,j,k)
           endif
     ...
           if(set_flag(i,j,k)==1) then
             DW(i,j,k) = D(i,j)*obs_w
           else
             DW(i,j,k) = DW(i,j,k)+dt*R4(i,j,k)
           endif

     ...
     ! update velocities for calculating IB forces
     call get_UVW     

     ! calculate forcing at obstacle boundary                                                                                    
     call imm_obs

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

       if(set_flag(i,j,k)==1) then
         DU(i,j,k) = D(i,j)*obs_u
         DV(i,j,k) = D(i,j)*obs_v
         DW(i,j,k) = D(i,j)*obs_w
       else
         DU(i,j,k) = DU(i,j,k)+dt*ObsForceX(i,j,k)
         DV(i,j,k) = DV(i,j,k)+dt*ObsForceY(i,j,k)
         DW(i,j,k) = DW(i,j,k)+dt*ObsForceZ(i,j,k)
       endif
     enddo
     enddo
     enddo




