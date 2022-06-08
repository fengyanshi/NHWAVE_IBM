Implementation of s-coordinate 
======================================================

**1) Module structure**

The code coordinate module was in mod\_coordinate.F. The module contains the original sigma-coordinate system in Ma et al 2012 and modification in Derakhti et al. 2016 a,b., the following is the model structure:

  .. code-block:: rest

    MODULE COORDINATE_MODULE
    USE GLOBAL, ONLY: Zero,DelxSc,DelySc,D,DelxH,DelyH, &
                       DelxEta,DelyEta,sigc,Mloc,Nloc,Kloc, &
                       DelxSl,DelySl,Kloc1,sig
    IMPLICIT NONE
    SAVE

    CONTAINS

      SUBROUTINE COORDINATE_INITIALIZATION
        | ...
      END SUBROUTINE COORDINATE_INITIALIZATION

      SUBROUTINE SIGMA_TRANSFORM
        | ...
      END SUBROUTINE SIGMA_TRANSFORM

      SUBROUTINE GAMMA_s
        | ...
      END SUBROUTINE GAMMA_s

      SUBROUTINE MEAN_VELOCITY
        | ...
      END SUBROUTINE MEAN_VELOCITY

      SUBROUTINE GET_OMEGA(R1)
        | ...
      END SUBROUTINE GET_OMEGA(R1)


    END MODULE COORDINATE_MODULE

**2) Procedures**

a. move the calculation of :math:`dsig` to  SUBROUTINE COORDINATE\_INITIALIZATION to generalize s-coordinate. NOTE: we still keep :math:`dsig` but for a generic :math:`\Delta s`. In this way, the code for solving the continuity equation does not need to change. 

b. move subroutine get_omega to SUBROUTINE GET\_OMEGA in MOD\_COORDINATE.F to make :math:`\omega` to the generic :math:`\tilde{\omega}`

c. calculate depth-averaged velocity :math:`\bar{u}_\alpha` in SUBROUTINE MEAN\_VELOCITY
  
d. calculate :math:`\Gamma_s` in SUBROUTINE GAMMA\_s based on

.. math::
   \Gamma_s = [D(u_\alpha - \bar{u}_\alpha)]_{,\alpha} + \tilde{\omega}_{,s} 

e. in SUBROUTINE eval\_duvw, add :math:`\Gamma_s` as an extra term in the momentum equations. Specifically, add :math:`u_{\alpha} \Gamma_s` and :math:`w \Gamma_s` in :math:`R2`, which is combination of advection, Coriolis, horizontal diffusion and extra forcing terms.   

f. add :math:`Gamma_s` term in the 3D mass conservation equation, resulting in an extra term on RHS of the Poisson equation. Following the previous approach to treat :math:`u, v, w` terms, we use explicit scheme. The :math:`Gamma_s` is implemented in SUBROUTINE generate\_coef\_rhs in pressure.F


