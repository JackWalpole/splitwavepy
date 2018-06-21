.. _background:

*******************
Covariance
*******************

The covariance :math:`{cov} (\mathbf{x},\mathbf{y})` of the energy measured on traces :math:`(\mathbf{x},\mathbf{y})` representing the two orthogonal components in the plane of the shear wave is fundamental to the grid search methods for finding the best fitting inverse splitting operator.  These methods are the Silver and Chan method which seeks to miniminise the energy on the transverse component, and the rotation correlation method which seeks to maximise the correlation between fast and slow wavelets.


The covariance between two traces :math:`(\mathbf{x}, \mathbf{y})` each containing :math:`N` samples is

.. math:: {cov} (\mathbf{x}, \mathbf{y})=\frac {\sum(\mathbf{x}-\mu _{x})(\mathbf{y}-\mu _{y})}{N},

where :math:`\mu_{x}` is the mean of :math:`\mathbf{x}`.  Note that the covariance of a trace with itself is simply its variance

.. math:: {cov} (\mathbf{x},\mathbf{x}) = {var} (\mathbf{x}) = \sigma_{x}^{2}.

The covariance matrix :math:`\Sigma(\mathbf{x},\mathbf{y})` for this can be viewed as

.. math:: \Sigma (\mathbf{x},\mathbf{y})) = \begin{bmatrix}
			                {var} (\mathbf{x}) & {cov} (\mathbf{x},\mathbf{y}) \\
									  & {var} (\mathbf{y})
			            \end{bmatrix}.


The covariance can be normalised to produce the Pearson correlation coefficient which is the statistic which the rotation correlation method seeks to maximise.

.. math:: \rho(\mathbf{x},\mathbf{y}) = \frac{ {cov} (\mathbf{x},\mathbf{y})}{\sigma_{x}\sigma_{y}}.

where :math:`\sigma_x` is the standard deviation of :math:`\mathbf{x}`.  The value of :math:`\rho` can be extracted from the covariance matrix.

.. math:: \Sigma (X,Y) = \begin{bmatrix}
			                \sigma_x^2 & \sigma_x \sigma_y \rho_{xy} \\
									  & \sigma_y^2
			            \end{bmatrix}.

Now consider the rotation correlation method which seeks to maximise the correlation between fast and slow traces.

.. math:: (X\star Y)[n]\ {\stackrel {\mathrm {def} }{=}}\sum _{m=-\infty }^{\infty }X[m]\ Y[m+n].

The correlation coefficient at lag :math:`n` can be normalised to recover the Pearson correlation coefficient

.. math:: \rho(X,Y)[n] = \frac{(X\star Y)[n]}{\sigma_{X}\sigma_{Y}}.

Where :math:`\sigma_{X} = \sqrt{{var} (X)}` is the standard deviation on trace :math:`X`.  The rotation correlation method seeks to find the rotation and lag that maximises :\math:`\rho`.

The rotation exchanges energy between traces

.. math:: (X,Y)_{\theta} = (X,Y)_{\theta} \begin{bmatrix}
	 								cos(\theta) & -sin(\theta)
									sin(\theta) & cos(\theta)
									\end{bmatrix}.
									
The covariance matrix for test angle :math:`\theta` at lag :math:`n` is then.

.. math:: \Sigma (X,Y)_{\theta} [n] = \begin{bmatrix}}
	 								{var} (X) &  cov \\
										& {var} (Y)
										\end{bmatrix}.
									
The lag can be easily converted to a meaningful delay time if the sampling interval :math:`\Delta` of the data is known

.. math:: \delta t = n \Delta.