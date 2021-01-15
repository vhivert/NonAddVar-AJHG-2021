/*
  PhenoCpp Simulator 
  
  This program simulate phenotypic data with additive, dominance 
  and additive-by-additive variance from bed files
  
  Copyright (C) 2021 The University of Queensland
  
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
   
  Update Version 1.1.0: We include the option to input the LD score file output from GCTA and use it to force 
  the sampling of some causal variants in a specific range of MAF/LD 
  
  Author: Valentin Hivert
  Date: 08/01/2021
*/
#define MAIN_UNIT

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
#include <ctime>

#include "defs.h"
#include "rand.h"
#include "read.h"
#include "write.h"
#include "simul.h"
#include "mem.h"

using namespace std;

int main(int argc, char *argv[]){
  // Input arguments
  string bfile       = "none";
  string prefix      = "none";
  bool verbose       =  true;
  double h2          =     0;      // Narrow sense heritability parameter
  double d2          =     0;      // Dominance heritability parameter
  double e2          =     0;      // Additive-by-additive heritability parameter
  int Ncausals       =     0;      // Number of causal variants
  int nthreads		 =	   1;	   // Number of threads used for computation
  unsigned int seed			 =     0;
  // Options for ascertained causal variants
  int AscertCausals       =     0;      
  double MinMaf	=	0;
  double MaxMaf	=	0;
  double MinLD	=	0;
  double MaxLD	=	0;
  string LDfile       = "none";
  
  // Indices
  int i;
  string sw;
  
  // Declaring argument for time() 
  time_t t;
  double start,end,elapsed; 
  // Declaring variable to store return value of 
  // localtime() 
  struct tm *local; 
  
  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  
  // Read arguments
  sw = argv[1];
  if (sw == "--help"){
    cerr<<"\t--bfile      : Binary PLINK format for genotypes."<<endl;
    cerr<<"\t--Ncausals          : Number of causal variants."<<endl;
    cerr<<"\t--h2         : SNP additive heritability."<<endl;
    cerr<<"\t--d2         : SNP dominance heritability."<<endl;
    cerr<<"\t--e2         : SNP additive-by-additive heritability."<<endl;
    cerr<<"\t--out        : prefix for output file: [prefix].sim."<<endl;
    cerr<<"\t--nthreads   : Number of threads used for computaion."<<endl;
    cerr<<"\t--seed   : seed to initiate RNG (by default: random seed)."<<endl;
    cerr<<"\t[Note] Missing values are imputed to the mean genotype."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  for(i = 1; i<argc;i++){
    sw = argv[i];
    if (sw == "--bfile"){
      bfile = argv[i + 1];
    }
    if (sw == "--LDfile"){
      LDfile = argv[i + 1];
    }
    if (sw == "--Ncausals"){
      Ncausals = atoi(argv[i + 1]);
    }
    if (sw == "--AscertCausals"){
      AscertCausals = atoi(argv[i + 1]);
    }
    if (sw == "--silent"){
      verbose = false;
    }
    if (sw == "--h2"){
      h2 = atof(argv[i + 1]);
    }
    if (sw == "--d2"){
      d2 = atof(argv[i + 1]);
    }
    if (sw == "--e2"){
      e2 = atof(argv[i + 1]);
    }
    if (sw == "--out"){
      prefix = argv[i + 1];
    }
    if (sw == "--nthreads"){
      nthreads = atoi(argv[i + 1]);
    }
    if (sw == "--seed"){
	  //seed = std::stoul (argv[i + 1],nullptr,0);
      seed = atoi(argv[i + 1]);
    }
    if (sw == "--MinMaf"){
      MinMaf = atof(argv[i + 1]);
    }
    if (sw == "--MaxMaf"){
      MaxMaf = atof(argv[i + 1]);
    }
    if (sw == "--MinLD"){
      MinLD = atof(argv[i + 1]);
    }
    if (sw == "--MaxLD"){
      MaxLD = atof(argv[i + 1]);
    }
  }
  

  if(bfile=="none"){
    cerr<<"\tA prefix must be specified for files [prefix].bed, [prefix].bim and [prefix].fam"<<endl;
    exit(1);
  }
  if(h2>1. or h2<0.){
    cerr<<"\tHeritability parameter must be between 0 and 1. Use option --h2."<<endl;
    exit(1);
  }
  if(d2>1. or d2<0.){
    cerr<<"\tHeritability parameter must be between 0 and 1. Use option --d2."<<endl;
    exit(1);
  }
  if(e2>1. or e2<0.){
    cerr<<"\tHeritability parameter must be between 0 and 1. Use option --e2."<<endl;
    exit(1);
  }
  if( (h2+d2+e2)>1. ){
    cerr<<"\tBroad sense heritability parameter (h2 + d2 + e2) cannot exceed 1. Use option --h2, --d2 and --e2."<<endl;
    exit(1);
  }
  string bedfile = bfile+".bed";
  string bimfile = bfile+".bim";
  string famfile = bfile+".fam";
  
  
  t = time(NULL);
  local = localtime(&t);
  start = omp_get_wtime();
  cout<<"------------------------------------------------------------------------------------"<<endl;
  printf("%s",asctime(local));
  printf("------------------------------------------------------------------------------------\n\n");
  
  printf("This simulation was performed using PhenoCppSimulator (version %s)\n",VERSION);
  printf("Author : Valentin Hivert\n\n");
  printf("Copyright (C) 2021  The University of Queensland\nThis program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.\n\nThis program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.\n\n");
  
  printf("------------------------------------------------------------------------------------\n"); 
  if(verbose){
    cout<<">>> Simulation of phenotype <<<"<<endl;
    cout<<"The program will use "<<nthreads<<" threads.\n";
    cout<<"Input files: "<<endl;
    cout<<"\tBED: "<<bfile<<".bed\n";
    cout<<"\tBIM: "<<bfile<<".bim\n";
    cout<<"\tFAM: "<<bfile<<".fam\n";
    if(LDfile!="none"){cout<<"/tLD: "<<LDfile<<"\n";}
    cout<<endl;
    if(AscertCausals>0){
		cout<<"Simulation uses "<<Ncausals<<" causal variants ("<<Ncausals-AscertCausals<<" random variants and"<<AscertCausals<<" ascertained variants)\n";
		cout<<"Ascertained variants are selected among those with "<<MinMaf<<"<=MAF<="<<MaxMaf<<" and "<<MinLD<<"<=LDscore<="<<MaxLD<<"\n";
		}else{
			cout<<"Simulation uses "<<Ncausals<<" causal variants\n";
		}
    
    cout<<"Specified additive heritability h2 = "<<h2<<"\n";
    cout<<"Specified dominance heritability d2 = "<<d2<<"\n";
    cout<<"Specified additive-by-additive heritability e2 = "<<e2<<"\n";
    cout<<"Output file name: "<<prefix<<".out\n";
    cout<<endl;
  }
  
  seed = Seed(seed);
  
  if (nthreads == 0) {
    nthreads = omp_get_max_threads();                                          // omp_get_max_threads() returns the same value whether executing from a serial or parallel region
  }
  omp_set_num_threads(nthreads);
    
  // Set data struct
  
  data_struct data;
  data.Ncausals=Ncausals;
  data.AscertCausals=AscertCausals;
  read_data(data, bimfile, famfile, LDfile, MinMaf, MaxMaf, MinLD, MaxLD, verbose);  
  
  // Sample the causal variants
  if(AscertCausals>0){
	  SampleCausalsWithAscert(data.nbr_loci, data.Ncausals, data.AscertCausals, data.IndexAscert, data.index_causal);
  }else{
	  SampleWithoutReplacement(data.nbr_loci, data.Ncausals, data.index_causal);
  }
  
  ///////////////////////
  
  // Reading Genotypes and computing Allele Frequencies  
  int numBytes   = (int)ceil((double)data.nbr_ind / PACK_DENSITY);
  char* packed   = new char[numBytes];
  char* unpacked = new char[numBytes * PACK_DENSITY];

  int ndropped = 0;
  if((Ncausals>0 && h2>0.0) || (Ncausals>0 && d2>0.0) || (Ncausals>0 && e2>0.0)){
    ndropped = readGenotypes(data,bedfile,data.nbr_ind,data.nbr_loci,packed,unpacked,numBytes);
  }
  int nc=data.Ncausals - ndropped;
  
  if(verbose){
    cout<<"The actual number of causal variants used is m = "<<nc<<".\n\n";
  }
  
  if(h2>0){
	cout<<"Simulate additive effects";  
	ComputeGeneticValues_omp(data, nc, h2, 1);
	cout<<"...done"<<endl; 
  }
  if(d2>0){
	 cout<<"Simulate dominance effects";  
	ComputeGeneticValues_omp(data, nc, d2, 2);
	cout<<"...done"<<endl; 
  }
  if(e2>0){
	 cout<<"Simulate additive-by-additive effects";  
	ComputeGeneticValues_omp(data, nc, e2, 3);
	cout<<"...done/n"<<endl; 
  }
  
  //Build Phenotype
  BuildPheno(data,(h2+d2+e2), verbose);
  
  //Write output
  string outfile = prefix+".out";
  Write_Pheno(data, outfile);
  
  outfile = prefix+".info";
  Write_CausalInfo(data, outfile);

  if(verbose){
    cout<<"[done]: analysis finished."<<endl;
  } 
  
  end = omp_get_wtime();
  elapsed = end - start;
  t = time(NULL);
  local = localtime(&t);
  printf("------------------------------------------------------------------------------------\n");
  printf("%s------------------------\n",asctime(local));
  printf("\tTotal computing time elapsed           = %s\n",print_time(elapsed));
  printf("------------------------------------------------------------------------------------\n\n");
  printf("The program has successfully terminated.\n\n");
  
  
  delete [] packed;
  delete [] unpacked;    
  
  free_data(data);
  
  return EXIT_SUCCESS;
}


