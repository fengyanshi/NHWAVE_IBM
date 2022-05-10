Implementation of s-coordinate 
======================================================

1) The code coordinate module was in mod\_coordinate.F. The module contains the original sigma-coordinate system in Ma et al 2012 and modification in Derakhti et al. 2016 a,b. 

  .. code-block:: rest

    MODULE COORDINATE_MODULE
    USE GLOBAL, ONLY: Zero,DelxSc,DelySc,D,DelxH,DelyH, &
                       DelxEta,DelyEta,sigc,Mloc,Nloc,Kloc, &
                       DelxSl,DelySl,Kloc1,sig
    IMPLICIT NONE
    SAVE

    CONTAINS

      SUBROUTINE COORDINATE_INITIALIZATION

      END SUBROUTINE COORDINATE_INITIALIZATION

      SUBROUTINE SIGMA_TRANSFORM

      END SUBROUTINE SIGMA_TRANSFORM

    END MODULE COORDINATE_MODULE

2) 