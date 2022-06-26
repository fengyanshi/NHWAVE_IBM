.. _input_s_coor:

Model input for s-coordinate 
======================================================

* Double-power s-coordinate

.. math::

   s = C_1 \sigma^{A_1} + C_2 \sigma^{A_2}

where :math:`C_1` and :math:`A_1` make finer upper layers, and :math:`C_2` and :math:`A_2` make finer lower layers. In input.txt, we use 'upper' and 'lower' represent '1' and '2', respectively. 

  .. code-block:: rest

    EXPONENTIAL = T
    C_upper = 1.0
    C_lower = 0.0
    A_upper = 0.7
    A_lower = 0.0

* Hyperbolic s-coordinate

.. math::

   s = \tanh(\alpha \sigma) / \tanh (\alpha)

where :math:`\alpha` is a parameter. 

  .. code-block:: rest

    HYPERBOLIC = T
    HYPER_alpha = 1.5

