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

   \Gamma_s = D_{,t} + u_\alpha D_{,\alpha} + D ( (s_{,t})_{,s} + u_\alpha (s_{,\alpha})_{,s} + w (s_{,z})_{,s} )


There are three choices of :math:`s` coordinate system in the model.

**1) The original** :math:`\sigma` **coordinate**

In the original :math:`\sigma` coordinate system, i.e., :math:`s = \sigma`, 

.. math::

   \Gamma_s = 0

and

.. math::

   \tilde{\omega} = \omega = D (\sigma_{,t} + u_\alpha \sigma_{,\alpha} + w \sigma_{,z}).

**2)** :math:`s = s(\sigma)`

The choice that :math:`s` is a function of :math:`\sigma` only is convenient in terms of  minimizing the impact on the structure of the governing equations. In this case, the governing equations, retaining :math:`s` derivatives and expressions in terms of :math:`s` rather than :math:`\sigma`, are given by

.. math::
   D_{,t} + ( D u_\alpha )_{,\alpha} + s_{,\sigma}(\sigma_{,s}\tilde{\omega})_{,s} & =  0 \\
   (D u_\alpha)_{,t} + (D u_\alpha u_\beta + \frac{1}{2} gD^2 \delta_{\alpha\beta} )_{,\beta} + s_{,\sigma}(\sigma_{,s} u_\alpha \tilde{\omega})_{,s} & =  S_{q\alpha} + S_{h\alpha} \\
  (Dw)_{,t} + (D u_\beta w)_{,\beta} + s_{,\sigma}(\sigma_{,s} w \tilde{\omega})_{,s} & = S_{qz}

Use of the transformation and retention of :math:`s` rather than :math:`\sigma` in all expressions leads at most to the introduction of a variable coefficient outside the :math:`s`-derivative terms in each equation. 

**3) General cases**

For cases where :math:`s` or :math:`\sigma` can't be expressed in terms of each other alone, the details of the transformation and the resulting expressions for :math:`\tilde{\omega}` and :math:`\Gamma_s` must be worked out for each individual case. 

