/*
  This file is part of PhenoCpp Simulator package developped
  by Valentin Hivert Copyright (C) 2021 The University of Queensland
  
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
#ifndef _RAND_H
#define _RAND_H
#define PI         3.14159265   // The value of pi
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
#include <algorithm>

#include "defs.h"
using namespace std;
//Declare Seed for RNG
unsigned int Seed(unsigned int seed);
// Declare engine - single instance for the whole code
extern std::mt19937 gen;
//extern std::mt19937_64 my_rng;

//Declare distributions:
extern std::uniform_real_distribution<double> unif_real_dist;
extern std::normal_distribution<double> norm;

void SampleWithoutReplacement(int populationSize, int sampleSize, bool *samples);
std::vector<int> SampleWithoutReplacement2( int populationSize, int sampleSize);
void SampleCausalsWithAscert(int nbr_loci, int Ncausals, int AscertCausals, vector<int> IndexAscert, bool *samples);
#endif
