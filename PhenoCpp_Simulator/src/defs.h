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
#include <iostream>
#include <string>
#include <string.h>
#include <random>
using namespace std;

#ifndef _DEFS_H
#define _DEFS_H

#ifndef MAIN_UNIT
#define GLOBAL_VARIABLE extern
#else
#define GLOBAL_VARIABLE
#endif

#define VERSION              "1.1.0"
#define repeat for (;;)
#define PACK_DENSITY 4
#define MASK0   3 // 3 << 2 * 0
#define MASK1  12 // 3 << 2 * 1
#define MASK2  48 // 3 << 2 * 2
#define MASK3 192 // 3 << 2 * 3

// ---------------------
// structure definitions
// ---------------------

struct data_struct {
	long int nbr_loci;    // Number of loci
	long int nbr_ind;     // Number of sampled individuals
	int Ncausals;         // Number of causal variants
	int AscertCausals;         // Number of ascertained causal variants
	double *freq;		  // Minor allele frequency for each loci	
	double *genvalue;
	double **SNPeffect;   //SNP effects for additive, dominance and additive-by-additive
	double **geno;        // genotypes per ind and loci (0 1 2 before scaling)
	string *info_loci;
	string *IID;
	string *FID;
	bool *index_causal;
	vector<int> IndexAscert; 
	double *Pheno;
	double *PhenoScaled;
};
#endif
