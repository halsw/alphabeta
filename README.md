#Alpha Beta Filters for Arduino and Teensy
A library with template classes for the Alpha(G), Alpha Beta(GH) and Alpha Beta Gamma(GHK) filters

##Process Noise
Process noise is the expected deviation of the filter input between iterations which must be defined in the filters as either [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) in the constructor or [variance](https://en.wikipedia.org/wiki/Variance) by calling the method **setProcessVAR()**

##Measurement Noise
Measurement noise is the [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) of the error between the measured and actual values which must be defined in the filters as either [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) in the constructor or [variance](https://en.wikipedia.org/wiki/Variance) by calling the method **setMeasurementVAR()**

##Sampling Period
Sampling period is the time between consecutive filter updates defined in the filters through the constructor or by calling the method **setPeriod()**. The sapmpling period defined in the filter should match the period the filter gets updated and for optimal operation of the filter the sampling intervals should be as close to uniform as possible, like shown in the example alphabeta.ino