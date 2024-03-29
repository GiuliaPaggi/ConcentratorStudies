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
* Cluster.h and Cluster.cpp, in which the cluster class is defined. The cluster constructor analyses all the TPs and segments to cluster within 10cm from one another,  while for digis the requirement to make a cluster is to have at least 10 in a layer in 20cm.
* Analiser.h and Analiser.cpp, where the Ntuple content is analysed.

In addition, the plotting is done with a dedicated Python program plotter.py configured using .json file. 
## Usage
To run the clustering and cluters analisis
 
      cmake -B build -S . -DCMAKE_BUILD_TYPE=Debug
      cmake --build build/
      ./build/analisis NTupleFile.root

To plot the results

       python plotter.py -b configPlotterBase.json configPlotter.json

## Output

At the end of the analysis, the results are written to a .root file.

## Dependencies
The project is being developed using:
* the standard c++ libraries *array*, *vector*, *iostream*, *fstream*, *string*, *algorithm*;
* ROOT libraries with v6.26/04. 
