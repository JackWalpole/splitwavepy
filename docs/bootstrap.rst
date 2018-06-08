.. _bootstrap:

************************
Bootstrapping
************************

Principal component analysis (PCA) fits an ellipse to time series data with the long and short axes determined by the variance of the data along principal directions found by eigen analysis of the covariance matrix.  If known the orientation of the principal axes can be set if the source polarisation of the input wave is known, this stabilises the result by adding a priori information.  The aspect ratio of the ellipse is a measure of the linearity.  This method does not rely on the order of the points in the trace.  Similarly, the cross-correlation method, which finds the Pearson correlation-coefficient of the point cloud, does not rely on the order of points.

A grid search of the statistic value is performed over the splitting correcting the a range of splitting parameters.  The point with the highest 