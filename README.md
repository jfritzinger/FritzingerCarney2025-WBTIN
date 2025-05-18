This is the code for Fritzinger and Carney (2025),
"Mechanisms of Tone-in-Noise Encoding in the Inferior Colliculus" published in the Journal of Neuroscience.

## Data 

Data is available for download at: https://osf.io/p4r82/. This zip file contains data used for all analyses, 
a PDF including rate responses and characteristics for all neurons included in the manuscript, and code. 

Data include a folder called Neural_Data which contains data for each neuron in the manuscript. 
Preprocessed data used in figure generation are included as  excel and .mat files. 

## Code Workflow 

### Generating Figures 

Run ```generate_figs.m``` with an input of the figure to be generated. For convienence, preprocessed data is included 
in the data folder due to lengthy compute times for some analyses. Functions for each figure, including 
supplemental material, are found in the _"figures"_ folder. 

Scripts to run the analyses for preprocessed data are found in the _"analysis"_ folder. Other functions used for
analysis and figure generation are found in the _"helper-functions"_ folder. 

### Running Models
#### Broad Inhibition Model 

The _"model-lat-inh"_ folder contains precompiled code necessary to run the broad inhibition model on Mac and 
Windows operating systems. This folder also contains an example script simulating a response to amplitude
modulated noise using this model. 

The _"model-fitting"_ folder includes functions that fit the broad inhibition model parameters to the example
neurons used in the manuscript. 

#### SFIE Model 

The _"UR_EAR_2022a"_ folder is the release of the AN and SFIE model used in this manuscript. The _"model-SFIE"_ 
folder contains wrapper functions for the AN and SFIE models.

#### Energy Model 

The _"model-energy"_ folder contains gammatone filters used in the manuscript. This folder includes two example scripts
that run the model and output a figure either using a population of neurons or a single neuron with the sliding 
WB-TIN stimulus. 
