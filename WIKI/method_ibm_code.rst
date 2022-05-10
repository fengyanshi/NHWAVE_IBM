IBM module
------------------------------

1) The code IBM module was in mod\_immersed.F

  .. code-block:: rest

    MODULE IBM_MODULE
    USE GLOBAL
    IMPLICIT NONE
     real(SP) :: obs_u,obs_v,obs_w
     real(SP), dimension(:,:,:), allocatable :: ObsForceX,ObsForceY,ObsForceZ
     real(SP), dimension(:,:,:), allocatable :: set_dist_x,set_dist_y,set_dist_z
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: Mask_Obs,Rx,Ry,Rz,Tmp_obs
     integer, dimension(:,:,:), allocatable :: set_flag,set_flag_old
    SAVE

    END MODULE IBM_MODULE

    

2) The module contains 

  .. code-block:: rest

    SUBROUTINE OBS_INITIALIZATION
      ! used to read obstacle files

     OPEN(7,FILE=TRIM(OBS_MASK_FILE),status='old')
       ...
     CLOSE(7)

     IF(OBS_Require_Rxyz) THEN
       OPEN(7,FILE=TRIM(OBS_Rx_FILE),status='old')     
         ...
       CLOSE(7)
       OPEN(7,FILE=TRIM(OBS_Ry_FILE),status='old')     
         ...
       CLOSE(7)
       OPEN(7,FILE=TRIM(OBS_Rz_FILE),status='old')     
         ...
       CLOSE(7)
     ENDIF
    END SUBROUTINE OBS_INITIALIZATION

  .. code-block:: rest

    SUBROUTINE SET_OBSVEL
     ! set speed of obstacles (not used in this application)
      ...
    END SUBROUTINE SET_OBSVEL

  .. code-block:: rest

    SUBROUTINE IMM_OBS
      ! Calculate forcing at the obstacle boundary
      ! Gangfeng original subroutine
    END SUBROUTINE IMM_OBS


