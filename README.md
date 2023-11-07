# **Nearest Neighbour personal microbiome identification**

## **Authors**: 

* Markus Hittmeir
* Rudolf Mayer  (<mailto:rmayer@sba-research.org>)

## **Background**
This repository contains the code and the experiments described in the paper

_Markus Hittmeir, Rudolf Mayer, and Andreas Ekelhart. Distance-based Techniques for Personal Microbiome Identification. In International Conference on Availability, Reliability and Security (ARES), Vienna, Austria, August 2022. DOI: http://dx.doi.org/10.1145/3538969.3538985_

which aims to identify matching pairs of microbiome samples, i.e. samples from the same individual, from two different databases, by an approach on thresholded nearest neighbours, and its extended version:

_Rudolf Mayer. Markus Hittmeir, and Andreas Ekelhart. Distance-based linkage of personal microbiome records for identification and its privacy implications. In Computers & Security, Volume 136, January 2024, DOI: http://dx.doi.org/10.1016/j.cose.2023.103538_

which adds additional experiments and hypothesis testing.


## **Data**
This project uses the same data as published by

_Franzosa EA, Katherine H, Meadow JF, Gevers D, Lemon KP, Bohannan BJM, Huttenhower C. Identifying personal microbiomes using metagenomic codes. Proceedings of the National Academy of Sciences (2015): 201423854._

which are available from https://huttenhower.sph.harvard.edu/idability. They are included in this repository, in the folder [data](data).

## **Contribution**
Our solution provides a stronger method on matching the underlying individuals.

## **Main Script**

- The main script to run from this repository is [nearest_neighbour-all_experiments.py](nearest_neighbour-all_experiments.py), which runs the nearest-neighbour based matching algorithm on all provided datasets.
The actual matching code is provided in [nearest_neighbour.py](nearest_neighbour.py), which is called by the previous script
- To run the original idability scripts, run [idability-all_experiments.py](idability-all_experiments.py)
- To perform the statistical test of differences, [run hypothesis_test_idability_vs_nearest_neighbour.py](hypothesis_test_idability_vs_nearest_neighbour.py)
