# EEGPreprocessing_Not_so_Basic

This repository was created for OHBM Education 2020. The code presents the different aspects discussed in the course. It combines different elements into a coherent analysis illustrating the choices to make when analysis EEG data, here from a simple ERP experiment.

The rendering for the Matlab/EEGLAB code (filtering,ICA denoising,referencing) can be seen [here](https://cpernet.github.io/EEGPreprocessing_Not_so_Basic/).
The code and rendering for the Pytjon/MNE Python can be seen [here](https://github.com/CPernet/EEGPreprocessing_Not_so_Basic/blob/master/code/demo_cov_phantom.ipynb)

## The data

The data used in this tutorial come from Wakeman and Henson (2015). In this experiment, MEG-EEG data were collected while subjects viewed famous, unfamiliar and scrambled faces. Each image was repeated twice (immediately in 50% of cases and 5–10 stimuli apart for the other 50%) and subjects pressed one of two keys with their left or right index finger indicating how symmetric they regarded each image relative to the average.
The data were prepared (i.e. EEG extracted, timing corrected, electrode position re-oriented, event recorded) by Dung Truong, Ramon Martinez & Arnaud Delorme and can be downloaded from [OpenNeuro](10.18112/openneuro.ds002718.v1.0.2).
/data folder contains sub011 used to illustrate changes related to preprocessing.  
**You must download the .fdt file and unzip it to re-run the analysis**

For the source analysis, this is done using phantom data available here: https://neuroimage.usc.edu/brainstorm/Tutorials/PhantomElekta

## Code and dependencies

## Data filtering, ICA denoising and referencing
Data filtering, ICA denoising and referencing are illustrated using Matlab and EEGLAB along with plug-ins (use the EEGLAB plug-in manager):

1. Zapline from the [NoiseToolbox](http://audition-backend.ens.fr/adc/NoiseTools/) already in the code folder
2. [ViewProps](https://sccn.ucsd.edu/wiki/Viewprops)
3. [ICLabel](https://sccn.ucsd.edu/wiki/ICLabel) to automated noise removal (labelling of ICA)
4. [REST](http://www.neuro.uestc.edu.cn/rest/Down.html) for referencing
5. [LIMO MEEG](https://limo-eeg-toolbox.github.io/limo_meeg/) for the statistical analyses

Also recommended:
[Automagic](https://github.com/methlabUZH/automagic) to automate the preprocessing steps illustrated in the script

## Source denoising
This method is demonstrated using Python and MNE python with simulated data (because we know the truth).
This can be run in the Jupyter notebook.

/code folder contains all the code necessary to run the analyses
