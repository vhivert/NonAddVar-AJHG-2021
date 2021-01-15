/*
  This file is part of PhenoCpp Simulator package developped
  by Valentin Hivert Copyright (C) 2021  University of Quuensland
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>. 
  
  Date: 08/01/2021
*/
#ifndef _SIMUL_H
#define _SIMUL_H
#include <string.h>
#include <math.h>
#include <omp.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <random>

#include "defs.h"
#include "rand.h"
using namespace std;


void ComputeGeneticValues_omp(data_struct &data, long int nbr_loci, double h, int model);
void BuildPheno(data_struct &data, double H, bool verbose);

#endif
