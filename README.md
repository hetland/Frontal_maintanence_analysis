# Frontal maintanence analysis
Analysis for paper "Frontal Maintenance in Submesoscale Flows"

This repository contains two notebooks that will reproduce the figures in the paper "Frontal Maintenance in Submesoscale Flows" by Robert Hetland, Lixin Qu, and Kyle Hinson, submitted to the Journal of Physical oceanography.

Two notebooks are in the `Notebooks` directory. The notebook `Frontal_maintenance_idealized_channel_flow_analysis.ipynb` contains scripts that will load in the six idealized channel flow configurations found in the the [zenodo archive](https://zenodo.org/records/11069051) under the `reduced_channel_files` directory. This directory should be placed in the same diretory as the notebook.

The notebook `Frontal_maintenance_realistic_model_analysis.ipynb` reproduces the figures for the realistic TXLA model. Data used in reproduction of the figure are in the same [zenodo archive](https://zenodo.org/records/11069051) as a zarr dataset directory `surface_child.zarr`. Again, place this diretory in the same diretory as the notebook to reproduce the figures.
