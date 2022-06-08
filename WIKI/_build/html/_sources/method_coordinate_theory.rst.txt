New s-Coordinate 
======================================================

According to Kirby et al. (2022, in preparation), the governing equations in a generalized :math:`s`-system :math:`(x,y,s,t)`, with a range :math:`0 \leq s \leq 1` that preserves the boundary locations from the :math:`\sigma` coordinate mapping, can be expressed by

.. math::

   D_{,t} + ( D u_\alpha )_{,\alpha} + \tilde{\omega}_{,s}  &=  \Gamma_s \\ 
  (D u_\alpha)_{,t} + (D u_\alpha u_\beta + \frac{1}{2} gD^2 \delta_{\alpha\beta} )_{,\beta} + (u_\alpha \tilde{\omega})_{,s}
 &=  S_{q\alpha} + S_{h\alpha} + u_\alpha \Gamma_s \\
  (Dw)_{,t} + (D u_\beta w)_{,\beta} + (w\tilde{\omega})_{,s}  &=  S_{qz} + w \Gamma_s

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

To evaluate :math:`\Gamma_s` at the same time level, replace :math:`D_{,t}` with the divergence of mass flux based on the continuity equation :math:`D_{,t} = - D \bar{u}_{\alpha,\alpha}`, where :math:`\bar{u}_\alpha` represents the depth averaged velocity,

.. math::
   \Gamma_s = [D(u_\alpha - \bar{u}_\alpha)]_{,\alpha} + \tilde{\omega}_{,s} 

There are two choices of :math:`s` coordinate system in the model.

**1) The original** :math:`\sigma` **coordinate**

In the original :math:`\sigma` coordinate system, i.e., :math:`s = \sigma`, 

.. math::

   \Gamma_s = 0

and

.. math::

   \tilde{\omega} = \omega = D (\sigma_{,t} + u_\alpha \sigma_{,\alpha} + w \sigma_{,z}).

**2) General cases**

For general cases, the details of the transformation and the resulting expressions for :math:`\tilde{\omega}` and :math:`\Gamma_s` must be worked out for each individual case. 

