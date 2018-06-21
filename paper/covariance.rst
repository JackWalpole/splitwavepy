.. _covariance:

*******************
Covariance
*******************

The covariance :math:`{cov} (\mathbf{x},\mathbf{y})` of the energy measured on traces :math:`(\mathbf{x},\mathbf{y})` representing the two orthogonal components in the plane of the shear wave is fundamental to the grid search methods for finding the best fitting inverse splitting operator.  These methods are the Silver and Chan method which seeks to miniminise the energy on the transverse component, and the rotation correlation method which seeks to maximise the correlation between fast and slow wavelets.


The covariance between two traces :math:`(\mathbf{x}, \mathbf{y})` each containing :math:`N` samples is

.. math:: {cov} (\mathbf{x}, \mathbf{y})=\frac {1}{N} \sum _{m=0}^{N-1}(\mathbf{x}[m] -\mu _{x})(\mathbf{y}[m]-\mu _{y}),

where :math:`\mu_{x}` is the mean of :math:`\mathbf{x}`.  Note that the covariance of a trace with itself is simply its variance

.. math:: 

	{cov} (\mathbf{x},\mathbf{x}) = \sigma_{x}^{2},
	
	{cov} (\mathbf{x}, \mathbf{y}) = {cov} (\mathbf{y}, \mathbf{x}),
	
	{cov} (\mathbf{x}, \mathbf{y}) = \sigma_x \sigma_y \rho_{\mathbf{x},\mathbf{y}}.
	
Where :math:`\sigma_x` is the standard deviation of :math:`\mathbf{x}`, and :math:`\rho` is the normalised correlation coefficient, also known as Pearson's r. The covariance matrix :math:`\mathbf{C}(\mathbf{x},\mathbf{y})` for this can be viewed as

.. math:: \mathbf{C}(\mathbf{x},\mathbf{y}) = \begin{bmatrix}
		{cov} (\mathbf{x},\mathbf{x}) & {cov} (\mathbf{x},\mathbf{y}) \\
									  & {cov} (\mathbf{y},\mathbf{y})
			            			  		\end{bmatrix} 
										  = \begin{bmatrix}
			                \sigma_x^2 & \sigma_x \sigma_y \rho_{\mathbf{x},\mathbf{y}} \\
									  & \sigma_y^2
			            			  		\end{bmatrix}.

Pearson's r is the statistic that is maximised in the rotation correlation method.  It can be simply plucked from the covariance matrix.

.. .. math:: \rho_{x,y} = \frac{C_{x,y}}{\sqrt{C_{x,x} C_{y,y}}}

.. .. math:: \rho_{x,y} = \frac{C_{x,y}}{(C_{x,x} C_{y,y})^{1/2}}

.. math:: \rho_{\mathbf{x},\mathbf{y}}  = \frac{C_{1,2}}{(C_{1,1} C_{2,2})^{1/2}}


If the mean of each trace is zero then the covariance matrix can be cheaply calculated.

.. math:: \mathbf{C}(\mathbf{x},\mathbf{y}) = \frac{1}{N}
								\begin{bmatrix}
						\sum{x_i^2} & \sum{x_i y_i} \\
				  	  			  & \sum{y_i^2}
			            		  \end{bmatrix} 

The energy of a discrete-time signal is equal to :math:`\sum_{i=0}^{N-1}x_{i}^2`.
If :math:`x` can be guaranteed to point in the *polarisation* direction, then :math:`y` will point in the *transverse* direction.  The Silver and Chan Method seeks to minimise the energy on the transverse component, provided :math:`x` and :math:`y` are suitably rotated then this statistic is simply the lower right entry of the covariance matrix.

If the polarisation of the shear wave is not known then a guess can be made that this principal direction is the principal eigenvector of the covariance matrix.  This is the Silver and Chan eigenvalue method which seeks to minimise the small eigenvalue of the covariance matrix :math:`\lambda_2`, which is equivalent to :math:`\sigma_y^2` if the polarisation really does math the principal eigenvector.

================
Grid Search
================

Now that we know which statistics from the covariance matrix are useful for measuring shear wave splitting we can consider the grid search.

The rotation exchanges energy between traces

.. math:: \mathbf{R}(\theta) = \begin{bmatrix}
	 								cos(\theta) & -sin(\theta) \\
									sin(\theta) & cos(\theta)
									\end{bmatrix}.

Notice that the term in the upper right corner has similar form to a cross-correlation with zero lag (:math:`j=0`).

The lag can be easily converted to a meaningful delay time if the sampling interval :math:`\Delta` of the data is known

.. math:: \delta t = n \Delta.

.. math:: (\mathbf{x} \star \mathbf{y})_j\ {\stackrel {\mathrm {def} }{=}}\sum _{m=-\infty }^{\infty }x_i y_{i+j}.


Furthermore, if energy in the window is conserved for all lags, then the terms in the upper left and lower right of the covariance matrix are constant.

.. math:: \mathbf{C}(\mathbf{x},\mathbf{y})_j = \frac{1}{N}
								\begin{bmatrix}
						\sum{x_i^2}  & (\mathbf{x} \star \mathbf{y})_j \\
				  	  			  & \sum{y_i^2} 
			            		  \end{bmatrix} 
									
As a consequence of the convolution theorem a computational speedup is achieved by calculating the cross-correlation in the frequency domain. 



.. math::
    
	
	\mathbf{X} = \mathcal{F} \{ \mathbf{x} \}
	
	\mathbf{Y} = \mathcal{F} \{ \mathbf{y} \}
	
	(\mathbf{x} \star \mathbf{y})_j = (\mathcal{F}^{-1} \{ \mathbf{X} \mathbf{Y}^* \})_j

By Parseval's theorem the sum of squares in the time domain signal is equal to the sum of squares on its Fourier transformed frequency domain companion.

.. math:: \sum _{i=0}^{N-1}x_{i}^2 = {\frac {1}{N}} \sum _{i=0}^{N-1}|X_i|^{2}

This leads to the efficient calculation of multiple covariance matrices at lag :math:`j` via

.. math:: \mathbf{C}(\mathbf{x},\mathbf{y})_j = \frac{1}{N}
								\begin{bmatrix}
								\frac {1}{N} \sum _{i=0}^{N-1}|X_i|^{2}  & 
								(\mathcal{F}^{-1} \{ \mathbf{X} \mathbf{Y}^* \})_j \\
				  	  			  & \frac {1}{N} \sum _{i=0}^{N-1}|Y_i|^{2}
			            		  \end{bmatrix} 
								  
=====================
Practical Issues
=====================

									
