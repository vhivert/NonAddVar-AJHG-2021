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
#include "rand.h"

std::mt19937 gen {}; // Defines an engine

std::uniform_real_distribution<double> unif_real_dist(0., 1.); //Define uniform real distribution over [0,1)
std::normal_distribution<double> norm(0.,1.);

// Function to seed the random number generator from main file
// useful if you want the seed from a parameter file
// a negative value for seed gets you a random seed
// outputs the seed itself
unsigned int Seed(unsigned int seed)
{
  if (seed == 0) {
	random_device rd;  
    unsigned int rseed=rd();
    cout<<"Randomizing random generator, seed is "<<rseed<<"/n"<<endl;
    
    gen.seed(rseed);
    return rseed;
  } else {
    cout<<"User-provided seed is "<<seed<<"/n"<<endl;
    gen.seed(seed);
    return seed;
  }
}

void SampleWithoutReplacement
(
    int populationSize,    // size of set sampling from
    int sampleSize,        // size of each sample
    bool *samples  // output, zero-offset indicies to selected items, initialement : vector<int> & samples 
)
{
    // Use Knuth's variable names
    int n = sampleSize;
    int N = populationSize;
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;

    while (m < n)
    {
        u = unif_real_dist(gen); // call a uniform(0,1) random number generator
        if ( (N - t)*u >= n - m )
        {
            samples[t] = false;
            t++; 
        }
        else
        {
            samples[t] = true;
            t++; m++;
        }
    }
}


std::vector<int> SampleWithoutReplacement2
(
    int populationSize,    // size of set sampling from
    int sampleSize        // size of the sample
)
{
    // Use Knuth's variable names
    std::vector<int> sample;
    int n = sampleSize;
    int N = populationSize;
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;

    while (m < n)
    {
        u = unif_real_dist(gen); // call a uniform(0,1) random number generator
        if ( (N - t)*u >= n - m )
        {
            t++; 
        }
        else
        {
            sample.push_back(t);
            t++; m++;
        }
    }
    return sample;
}



void SampleCausalsWithAscert
(
	int nbr_loci,
	int Ncausals,
    int AscertCausals,
    vector<int> IndexAscert,
    bool *samples
)
{
    int i=0;
    std::vector<int> IndexSampCausals;
    
		if(IndexAscert.size()==0){cerr << "Error : Not enough SNPs for the ascertainment sampling"<<endl;
           exit(1);}
           
		std::vector<int> AscertSample = SampleWithoutReplacement2(IndexAscert.size(), AscertCausals);
		std::vector<int> RandSample = SampleWithoutReplacement2((nbr_loci-AscertCausals),(Ncausals-AscertCausals));
		std::vector<int> tmp_loci(nbr_loci);
		for(i=0;i<nbr_loci;i++){
			tmp_loci[i]=i;
		}
		
		int idx=0;
		for(i=0;i<AscertSample.size();i++){
			tmp_loci[IndexAscert[AscertSample[i]]]=-1;;
			samples[IndexAscert[AscertSample[i]]]=true;
			idx++;
		}
		cout<<idx<<" ascertained SNPs"<<endl;
		tmp_loci.erase(std::remove(tmp_loci.begin(), tmp_loci.end(), -1), tmp_loci.end());
		cout<<"Number of possible SNP for random sample = "<<tmp_loci.size()<<endl;
		
		
		idx=0;
		for(i=0;i<RandSample.size();i++){
			if(samples[tmp_loci[RandSample[i]]]==true){
			cout<<"Error : loci "<<tmp_loci[RandSample[i]]<<" already sampled from ascertainment sampling"<<endl;
			exit(1);	
			}
			samples[tmp_loci[RandSample[i]]]=true;
		}
		for(i=0;i<nbr_loci;i++){
			if(samples[i]==true){idx++;}
		}
		cout<<"Number of total causal variants = "<<idx<<endl;
}
