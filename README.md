# [Alpha Beta Filters](https://github.com/halsw/alphabeta) for Arduino and Teensy
A library with template classes for the α(G), α-β(GH) and α-β-γ(GHK) filters

The α-β-γ filters are quite useful because of their fast execution and low memory usage that make them valuable for real time applications needing accurate estimates from noisy readings. To implement the alpha filter in a class it may seem an overkill but it provides an easy way to experiment with different filters and types to find the best suited for an application.

## Process Noise
Process noise is the expected deviation of the filter input between iterations which must be defined in the filters as either [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) in the constructor or [variance](https://en.wikipedia.org/wiki/Variance) by calling the method **setProcessVAR()**

## Measurement Noise
Measurement noise is the [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) of the error between the measured and actual values which must be defined in the filters as either [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) in the constructor or [variance](https://en.wikipedia.org/wiki/Variance) by calling the method **setMeasurementVAR()**

## Sampling Period
Sampling period is the time between consecutive filter updates defined in the filters through the constructor or by calling the method **setPeriod()**. The sampling period defined in the filter should match the period the filter gets updated and for optimal operation of the filter the sampling intervals should be as close to uniform as possible, like shown in the example alphabeta.ino