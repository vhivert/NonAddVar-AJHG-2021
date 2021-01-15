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
#include "mem.h"

void init_data(data_struct &data){
  long int i;
  //Allocate all memory to data
  data.FID = new string[data.nbr_ind];
  data.IID = new string[data.nbr_ind];
  data.info_loci = new string[data.nbr_loci];
  data.index_causal = new bool[data.nbr_loci];
  data.freq = new double[data.Ncausals];		  // Minor allele frequency for each loci
  data.geno = (double **) malloc(data.nbr_ind * sizeof(double *));
	for (i = 0; i < data.nbr_ind; ++i) {
		data.geno[i]=(double *) calloc(data.Ncausals,sizeof(double));
	}
  data.SNPeffect = (double **) malloc(data.Ncausals * sizeof(double *));
	for (i = 0; i < data.Ncausals; ++i) {
		data.SNPeffect[i]=(double *) calloc(3,sizeof(double));
	}
  data.genvalue = new double[data.nbr_ind]; //tmp storage of individual genetic values
  data.Pheno = new double[data.nbr_ind];
  data.PhenoScaled = new double[data.nbr_ind];
}

void free_data(data_struct &data){
  long int i;
  //Free memory to data
  delete [] data.FID;
  delete [] data.IID;
  delete [] data.info_loci;
  delete [] data.genvalue;
  delete [] data.freq;
  for(int i = 0; i < data.nbr_ind; ++i) {
        delete [] data.geno[i];   
  }
  delete [] data.geno;
  for(int i = 0; i < data.Ncausals; ++i) {
        delete [] data.SNPeffect[i];   
  }
  delete [] data.SNPeffect;
  delete [] data.Pheno;
  delete [] data.PhenoScaled;
  
}
