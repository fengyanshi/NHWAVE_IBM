New s-Coordinate 
======================================================

According to Kirby et al. (2022, in preparation), the governing equations in a generalized :math:`s`-system :math:`(x,y,s,t)`, with a range :math:`0 \leq s \leq 1` that preserves the boundary locations from the :math:`\sigma` coordinate mapping, can be expressed by

.. math::
   :label: mass

   D_{,t} + ( D u_\alpha )_{,\alpha} + \tilde{\omega}_{,s}  =  \Gamma_s 

and

.. math::
   :label: momx
 
   (D u_\alpha)_{,t} + (D u_\alpha u_\beta + \frac{1}{2} gD^2 \delta_{\alpha\beta} )_{,\beta} + (u_\alpha \tilde{\omega})_{,s} =  S_{q\alpha} + S_{h\alpha} + u_\alpha \Gamma_s

.. math::
   :label: momz

   (Dw)_{,t} + (D u_\beta w)_{,\beta} + (w\tilde{\omega})_{,s}  =  S_{qz} + w \Gamma_s

The modified normal velocity :math:`\tilde{\omega}` is given by

.. math::
   :label: omega

    \tilde{\omega} = D ( s_{,t} + u_\alpha s_{,\alpha} + w s_{,z} )

and the inhomogeneous term :math:`\Gamma_s` is given by

.. math::
   :label: gamma_s

   \Gamma_s = D_{,t} + u_\alpha D_{,\alpha} + D [ (s_{,t})_{,s} + u_\alpha (s_{,\alpha})_{,s} + w (s_{,z})_{,s} ]

or 

.. math::
   :label: gamma_s1

   \Gamma_s = & D_{,t} + u_\alpha D_{,\alpha} + \tilde{\omega}_{,s} - D (u_{\alpha,s} s_{,\alpha} + w_{,s} s_{,z}) \\
       =& D_{,t} + (u_\alpha D)_{,\alpha} + \tilde{\omega}_{,s}

Which is the extended forms of :eq:`mass`. The current version of NHWAVE calculates the normal velocity :math:`\omega` (:math:`\tilde{\omega}` here) using :eq:`mass` (without :math:`\Gamma_s`), following the method used in typical hydrostatic ocean models such as ROMS and FVCOM. The method won't work because :math:`\Gamma_s = 0` may not hold for a generic s-coordinate. In this study, we calculate :math:`\tilde{\omega}` using :eq:`omega` and solve :math:`\Gamma_s` using :eq:`mass`.  


To get the pressure Poisson equation, we substitute the pressure-corrected velocities 

.. math::

   u^{(k)} = u^* - \frac{\Delta t}{\rho} (\frac{\partial p}{\partial x} + \frac{\partial p}{\partial s}\frac{\partial s}{\partial x^*})^{(k)}

.. math::

   v^{(k)} = v^* - \frac{\Delta t}{\rho} (\frac{\partial p}{\partial y} + \frac{\partial p}{\partial s}\frac{\partial s}{\partial y^*})^{(k)}

.. math::
   w^{(k)} = w^* - \frac{\Delta t}{\rho} \frac{\partial s}{\partial z^*} (\frac{\partial p}{\partial s})^{(k)} 

into the 3D mass conservation equation 

.. math::

   \frac{\partial u}{\partial x} + \frac{\partial u}{\partial s} \frac{\partial s}{\partial x^*} + \frac{\partial v}{\partial y} + \frac{\partial v}{\partial s} \frac{\partial s}{\partial y^*} + \frac{\partial s}{\partial z^*} \frac{\partial w}{\partial s} = 0

results in

.. math::
  :label: pressure

  \frac{\partial}{\partial x} \left [ \frac{\partial p}{\partial x}+\frac{\partial p}{\partial s} \frac{\partial s}{\partial x^*} \right ] & +  \frac{\partial}{\partial y} \left [ \frac{\partial p}{\partial y}+\frac{\partial p}{\partial s} \frac{\partial s}{\partial y^*} \right ] + \frac{\partial }{\partial s} \left ( \frac{\partial p}{\partial x} \right) \frac{\partial s}{\partial x^*} + \frac{\partial }{\partial s} \left ( \frac{\partial p}{\partial y} \right) \frac{\partial s}{\partial y^*} \\ & +  \left[  \left (\frac{\partial s}{\partial x^*} \right )^2+\left ( \frac{\partial s}{\partial y^*} \right )^2+ \underline{(\frac{\partial s}{\partial z^*})^2}  \right] \frac{\partial }{\partial s} \left( \frac{\partial p}{\partial s} \right) \\ &= \frac{\rho}{\Delta t} \left( \frac{\partial u^*}{\partial x} +\frac{\partial u^*}{\partial s} \frac{\partial s}{\partial x^*} + \frac{\partial v^*}{\partial y} +\frac{\partial v^*}{\partial s} \frac{\partial s}{\partial y^*} + \underline{\frac{\partial s}{\partial z^*}} \frac{\partial w^*}{\partial s} \right)

:math:`\underline{()}` represents modifications needed for the new s-coordinate. 


To get :math:`\frac{\partial s}{\partial x^*}`,  :math:`\frac{\partial s}{\partial y^*}`, :math:`\frac{\partial s}{\partial z^*}`, and :math:`\frac{\partial s}{\partial t^*}` needed by the momentum equations and the pressure Poisson equation, we introduce :math:`z_s` function which describes :math:`z`-location in the s-coordinate system :math:`z_s(x,y,s,t)`, 


The derivatives of :math:`s` with respect to :math:`(x^*,y^*,z^*,t^*)` can be derived as, for example, 

.. math::

   \frac{\partial z_s}{\partial x^*} = \frac{\partial z_s}{\partial x} + \frac{\partial z_s}{\partial s} \frac{\partial s}{\partial x^*} = 0

or

.. math::

    \frac{\partial s}{\partial x^*} = - \frac{\partial z_s}{\partial x} / \frac{\partial z_s}{\partial s}

Similarly,

.. math::

    \frac{\partial s}{\partial y^*} = - \frac{\partial z_s}{\partial y} / \frac{\partial z_s}{\partial s}

and 

.. math::

    \frac{\partial s}{\partial t^*} = - \frac{\partial z_s}{\partial t} / \frac{\partial z_s}{\partial s}

The discretized forms are

.. math::
   :label: sx

   \frac{\partial s}{\partial x^*} = - \frac{\frac{\Delta z_s}{\Delta x}|_{s}} {\frac{\Delta z_s}{\Delta s}}

.. math::
   :label: sy

   \frac{\partial s}{\partial y^*} = - \frac{\frac{\Delta z_s}{\Delta y}|_{s}} {\frac{\Delta z_s}{\Delta s}}

.. math::
   :label: sz

   \frac{\partial s}{\partial z^*} = \frac{1}{\frac{\Delta z_s}{\Delta s}}

.. math::
   :label: st

   \frac{\partial s}{\partial t^*} = - \frac{\frac{\Delta z_s}{\Delta t}|_{s}} {\frac{\Delta z_s}{\Delta s}}

for given :math:`z_s(x,y,s,t)`. :math:`|_s` is used to emphasize the derivatives are calculated with constant :math:`s`. 


