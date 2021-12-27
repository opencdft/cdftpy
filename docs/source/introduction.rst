.. role:: raw-latex(raw)
    :format: latex

Introduction
==============

Classical density-functional theory (CDFT) provides statistical mechanics description of
molecular many-body systems  in terms of the density averages of atomic sites in the molecule:


.. math::
    \rho_{\alpha}(\mathbf{r})=\left < \sum_{i=1}^{N} \delta\left(\mathbf{r}-\mathbf{r}_{i\alpha}\right)\right >

where :math:`i, \alpha` refer to molecule and atomic site indices correspondingly.

From practical point of view, the main advantage of
CDFT is that it provides direct approximation of phase space averages,
thus avoiding numerical expense of direct numerical sampling with
molecular dynamics or Monte-Carlo methods.


