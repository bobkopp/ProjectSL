ProjectSL

README file last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 26 July 2019

## Citation

This code was used to produce the results of

	R. E. Kopp, R. M. Horton, C. M. Little, J. X. Mitrovica, M. Oppenheimer,
	D. J. Rasmussen, B. H. Strauss, and C. Tebaldi (2014). Probabilistic 21st
	and 22nd century sea-level projections at a global network of tide	gauge
	sites. Earth's Future 2: 287–306, doi:10.1002/2014EF000239. 

Please cite that paper when using any results generated with this code.

## Overview

This code contains two directories, slr and lib. slr contains code for analyzing tide gauge data and generating sea-level rise projections. lib contains supporting files.

This code requires MATLAB to run. It uses the Optimization and Mapping toolkits, though some of the functionality should be available without those toolkits. It also requires one of two routines to invert the spherical harmonics of the sea-level fingerprint files. By default, the code uses ssht (https://astro-informatics.github.io/ssht/index.html), although a slower MATLAB-based routine by Frederik Simons is included and can be used as a fallback. 

Needed input files, which go in the IFILES directory (specified in configureSLRProjections). These include:

	* the CSIRO GSL reconstruction,
	* the ICE5G-VM290 GIA model (NetCDF),
	* land ice static-equilibrium fingerprints,
	* Marzeion et al. 2012 glacier and ice cap projections,
	* PSMSL tide gauge data.

## Sea level rise projections

runTrainGPSLModel.m will generate a set of parameter files with optimized hyperparameters for each of the regions described in the coastlines.txt parameter files. (It is recommended that you use the default specifications, which are stored in slr/PARAMS; the training process is slow).

runSLRProjections.m will generate the sea-level rise projections, using the configuration specified in configureSLRProjections.m and generating the output files specified in outputSLRProjections.m. You will need the slr/, slr/MFILES, and slr/MFILES/scripts directories in your path.

You will need to modify the paths in configureSLRProjections.m to match your system.

----

    Copyright (C) 2019 by Robert E. Kopp

   
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.