# Wavelet Block-Thresholding Denoising 
=================================================================================

This repository contains MATLAB scripts and sample seismic data for appying denoising algorithm proposed in:

Mousavi S. M., and C. A. Langston (2016). Hybrid Seismic denoising Using Higher Order Statistics and Improved Wavelet Block Thresholding, Bulletin of the Seismological Society of America,106 (4), 1380-1393,doi:10.1785/0120150345

`demo.m` includes all info you need to know for running the code. 

you need `MATLAB statistics and signal processing toolboxes` to run this code.

## Paper
(https://www.researchgate.net/publication/303849872_Hybrid_Seismic_Denoising_Using_Higher-Order_Statistics_and_Improved_Wavelet_Block_Thresholding)

## Talk 
(https://earthquake.usgs.gov/contactus/menlo/seminars/1093)

## Abstract 
We introduce a nondiagonal seismic denoising method based on the continuous wavelet transform with hybrid block thresholding (BT). 
Parameters for the BT step are adaptively adjusted to the inferred signal property by minimizing the unbiased risk estimate of 
Stein (1980). The efficiency of the denoising for seismic data has been improved by adapting the wavelet thresholding and 
adding a preprocessing step based on a higher-order statistical analysis and a postprocessing step based on Wiener filtering. 
Application of the proposed method on synthetic and real seismic data shows the effectiveness of the method for denoising and
improving the signal-to-noise ratio of local microseismic, regional, and ocean bottom seismic data.

## A Short Description 
Seismic data recorded by surface arrays are often contaminated by unwanted noise. In many conventional seismic methods, 
the reliability of the seismic data and accuracy of parameter extraction, such as onset time, polarity, and amplitude, 
are directly affected by the background noise level. As a result, the accuracy of event location and other attributes 
derived from seismic traces are also influenced by the noise content. Therefore, there is a great need for developing 
suitable procedures that improve signal-to-noise ratios allowing for robust seismic processing. In this presentation, 
I introduce four different methods for automatic denoising of seismic data. These methods are based on the time-frequency 
thresholding approach. The efficiency and performance of the thresholding-based method for seismic data have been improved 
significantly. Proposed methods are automatic and data driven in the sense that all the filter parameters for denoising are 
dynamically adjusted to the characteristics of the signal and noise. These algorithms are applied to single channel data 
analysis and do not require large arrays of seismometers or coherency of arrivals across an array. Hence, they can be applied
to every type of seismic data and can be combined with other array based methods. Results show these methods can improve 
detection of small magnitude events and accuracy of arrival time picking.

![a-Induced microearthquake, b-local earthquake recorded by oceanic bottom seismometer, and c-regional earthquake. 
Each major panel shows the original time- series data in the upper left panel and its CWT to the right. Below are
the denoised seismogram and its CWT for comparison.](Fig.png)

a)Induced microearthquake, b)local earthquake recorded by oceanic bottom seismometer, and c)regional earthquake. 
Each major panel shows the original time- series data in the upper left panel and its CWT to the right. Below are
the denoised seismogram and its CWT for comparison.

