.. _implement_s_coor:

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


    END MODULE COORDINATE_MODULE

**2) Procedures**

* add GAMMA\_s in mod_global.F and initialize.F

* move the calculation of :math:`dsig` to  SUBROUTINE COORDINATE\_INITIALIZATION to generalize s-coordinate. NOTE: we still keep :math:`dsig` but for a generic :math:`\Delta s`. In this way, the code for solving all equations containing :math:`dsig` and associated derivative terms does not need any change. 

* in SUBTOUTINE SIGMA\_TRANSFORM in mod_coordinate.F, calculate :math:`\frac{\partial s}{\partial x^*}`,  :math:`\frac{\partial s}{\partial y^*}`, :math:`\frac{\partial s}{\partial z^*}`, and :math:`\frac{\partial s}{\partial t^*}` using :eq:`sx`, :eq:`sy`, :eq:`sz`, and :eq:`st`.

* remove get_omega(R1) in the original nhwave.F, in which :math:`\omega` is solved using the 3D mass conservation equation (1). 

* create the subroutine get_omega(R1) in mod_coordinate.F. The :math:`\omega` will be the generic :math:`\tilde{\omega}`. NOTE: :math:`\tilde{\omega}` here is calculated using :eq:`omega`, not :eq:`mass`. The residual term :math:`\Gamma_s` is calculated using :eq:`mass`.
  
* calculate :math:`\Gamma_s` in SUBROUTINE CALCULATE_GAMMAS in nhwave.F based on

.. math::

   \Gamma_s = D_{,t} + ( D u_\alpha )_{,\alpha} + \tilde{\omega}_{,s}  

where :math:`D_{,t}` has been calculated in subroutine eval\_duvw, named R1. 

* in SUBROUTINE eval\_duvw, add :math:`\Gamma_s` as extra terms in the momentum equations. Specifically, add :math:`u_{\alpha} \Gamma_s` and :math:`w \Gamma_s` in :math:`R2`, :math:`R3` and :math:`R4`, which are combination of advection, Coriolis, horizontal diffusion and extra forcing terms.   

* in pressure.F, modify the pressure Poisson solver based on :eq:`pressure`. Specifically, in SUBROUTINE\_GENERATE\_COEF\_RHS, remove DelxS, DelyS, and DelzS because they have been calculated using the generic s-coordinate in SUBTOUTINE SIGMA\_TRANSFORM. 



