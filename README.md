# ConcentratorStudies
 
## Project Status
The propagation of the track using the $\phi$ angle information and the clustering algorithm are completed, to be optimized.


 ## Table of Contents
* [General Info](#general-information)
* [Structure of the project](#structure-of-the-project)
* [Usage](#usage)
* [Output](#output)
* [Dependencies](#dependencies)


## General information

We want to study the distribution of low-quality trigger primitives (tps) in the barrel ($\eta$ < 0.8). 

The first step was to start with high-quality trigger primitives to propagate the track to other stations and match the prediction with tps in that chamber. 
To study the ghost distribution,  we grouped closeby tps into objects named clusters, separating the BestTP (highest quality at the right BX), the in-time ghost (lower qualities at the right BX) and the out-of-time ghosts (at the wrong BX). 

To have a complete description the cluster also holds information about the reconstructed segment and the digi in the same position in the chamber. If there is more than one segment reconstructed we choose the one with the most hits.

Furthermore, to study all possible types of ghosts, we also build a cluster with no tps around one or more segments. 

While to have a cluster made of digis, we ask for at least 10 close to each other in a superlayer. 

If we have a segment or a tp in the cluster, we also try to match the cluster information with the muon track extrapolated using the inner tracker information. When available the segment is preferred.


## Structure of the project
The project is divided into several files:
* TriggerPrimitives.h and TriggerPrimitives.cpp, in which the trigger primitive class is defined. Each trigger primitive is defined by an index in the event, the object stores the position and allows to predict the track expected phi in the other stations;
* Segment.h and Segment.cpp, in which the segment class is defined. It only has attributes defining the segment with index and position and the constructor;
* Digi.h and Digi.cpp, in which the digi class is defined. To compute the digis position we make use of a geometry file *DTGeom.txt* that allows converting the information from the local coordinates that are in the NTuple to the global reference used for the segments and trigger primitives.
* Analiser.h and Analiser.cpp, where the Ntuple content is analysed.

## Usage

For now, the NTuple file name is hard-coded at https://github.com/GiuliaPaggi/ConcentratorStudies/blob/main/Analiser.cpp#L73 and https://github.com/GiuliaPaggi/ConcentratorStudies/blob/main/Analiser.h#L658. 
To run the analiser, once in ROOT use:  
      
      .L TriggerPrimitive.cpp
      .L Segment.cpp
      .L Digi.cpp
      .L Cluster.cpp
      .L Analiser.cpp
      Analiser a
      a.Loop()


## Output

At the end of the analysis, the most useful information is plotted and shown on screen but also saved as a .png in a hard-coded directory and written to a .root file.

## Dependencies
The project is being developed usign :
* the standard c++ libraries *array*, *vector*, *iostream*, *fstream*, *string*, *algorithm*;
* ROOT libraries with v6.26/04. 
