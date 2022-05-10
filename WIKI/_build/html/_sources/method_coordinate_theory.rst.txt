New s-Coordinate 
======================================================

Let :math:`s(x'_\alpha,z',t')` be the replacement for :math:`z'` in a new coordinate system :math:`(x,y,s,t)`, with a range :math:`0 \leq s \leq 1` that preserves the boundary locations from the :math:`\sigma` coordinate mapping.  Following the same procedure as used in the :math:`\sigma` transformation leads to the modified system of equations

.. math::

   D_{,t} + ( D u_\alpha )_{,\alpha} + \tilde{\omega}_{,s}  =  \Gamma_s  

.. math::

  (D u_\alpha)_{,t} + (D u_\alpha u_\beta + \frac{1}{2} gD^2 \delta_{\alpha\beta} )_{,\beta} + (u_\alpha \tilde{\omega})_{,s}
 =  S_{q\alpha} + S_{h\alpha} + u_\alpha \Gamma_s 

.. math::

  (Dw)_{,t} + (D u_\beta w)_{,\beta} + (w\tilde{\omega})_{,s}  =  S_{qz} + w \Gamma_s

in coordinates :math:`(x,y,s, t)`.  The modified normal velocity :math:`\tilde{\omega}` is given by

.. math::

    \tilde{\omega} = D ( s_{,t} + u_\alpha s_{,\alpha} + w s_{,z} )

and the inhomogeneous term :math:`\Gamma_s` is given by

.. math::

   \Gamma_s = D_{,t} + u_\alpha D_{,\alpha} + D ( (s_{,t})_{,s} + u_\alpha (s_{,\alpha})_{,s} + w (s_{,z})_{,s} )


The term :math:`\Gamma_s` complicates the set of equations unless we are able to find transformations which reduce :math:`\Gamma_s` to zero, as occurs for the original choice :math:`s = \sigma`.  Such a choice does not appear to be possible.  A convenient choice which minimizes the impact on the structure of the governing equations is to let :math:`s` be specified entirely as a function of :math:`\sigma`, considered further in the cases later.  For :math:`s(\sigma)` only, we obtain the derivatives

.. math::

   s_{,t} = s_{,\sigma} \sigma_{,t}; s_{,\alpha} = s_{,\sigma} \sigma_{,\alpha}; s_{,z} = s_{\sigma} \sigma_{,z}

In this case, the normal velocity is given by :math:`\tilde{\omega} = s_{,\sigma} \omega` with :math:`\omega` given by (\ref{eq10}), and thus :math:`\tilde{\omega}` satisfies the same boundary conditions as :math:`\omega` at the upper and lower boundaries.  Using the substitution for the derivatives of :math:`s` in (\ref{eq17}) allows :math:`\Gamma_s` to be written as

.. math::

   \Gamma_s = \Gamma_\sigma - (\sigma_{,s})_{,\sigma}\tilde{\omega} \label{eq19}

where :math:`\Gamma_\sigma = 0` as stated below (\ref{eq12}). Further, the combination of the :math:`s` derivative term and the residual reduces to

.. math::
   \tilde{\omega}_{,s} - \Gamma_s = s_{,\sigma} ( \sigma_{,s} \tilde{\omega} )_{,s}

which is equivalent to the term :math:`\omega_{,\sigma}` in the original equation.  The full set of equations, retaining $s$ derivatives and expressions in terms of :math:`s` rather than $\sigma$, are given by

.. math::
   D_{,t} + ( D u_\alpha )_{,\alpha} + s_{,\sigma}(\sigma_{,s}\tilde{\omega})_{,s}  =  0 

.. math::

   (D u_\alpha)_{,t} + (D u_\alpha u_\beta + \frac{1}{2} gD^2 \delta_{\alpha\beta} )_{,\beta} + s_{,\sigma}(\sigma_{,s} u_\alpha \tilde{\omega})_{,s}  =  S_{q\alpha} + S_{h\alpha}

.. math::

   (Dw)_{,t} + (D u_\beta w)_{,\beta} + s_{,\sigma}(\sigma_{,s} w \tilde{\omega})_{,s}  = S_{qz} 

Use of the transformation and retention of :math:`s` rather than :math:`\sigma` in all expressions leads at most to the introduction of a variable coefficient outside the :math:`s`-derivative terms in each equation.  For cases where :math:`s` or :math:`\sigma` can't be expressed in terms of each other alone, the details of the transformation and the resulting expressions for :math:`\tilde{\omega}` and :math:`\Gamma_s` must be worked out for each individual case, as shown in the examples in the Section next.
