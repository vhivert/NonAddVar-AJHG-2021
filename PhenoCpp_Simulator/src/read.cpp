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
#include "read.h"

void read_data(data_struct &data, string bimfile, string famfile,string LDfile, double MinMaf, double MaxMaf, double MinLD, double MaxLD, bool verbose){
  long int i,k;
  
  string line = "";
  string tok  = ""; 
  ifstream tmpStream;
  if(verbose){
    cout<<"1-Counting the number of SNPs..."<<endl;
  }  
  long int p = -1;
  tmpStream.open(bimfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    p++;
  }
  tmpStream.close();
  if(verbose){
    cout<<"1-[done]: found "<<p<<" SNPs."<<endl;
  }
  data.nbr_loci=p;
  if(p<data.Ncausals){
	 cerr<<"WARNING : The number of SNPs is lower than the number of SNPs in the dataset."<<endl;
	 exit(1); 
  }

  // Get sample size
  if(verbose){
    cout<<"2-Counting the number of samples..."<<endl;
  }    
  long int n = -1; 
  tmpStream.open(famfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    n++;
  }
  tmpStream.close();  
  if(verbose){
    cout<<"2-[done]: found "<<n<<" samples."<<endl;
  }
  data.nbr_ind=n;
  
  
  if(verbose){
    cout<<"3-: Allocating memory...";
  }
  init_data(data);
  if(verbose){
    cout<<"done\n"<<endl;
  }
  
  //Get individuals ID
  if(verbose){
    cout<<"4-: Getting Individuals IDs...";
  }
  int nColFamFile = 6;  
  tmpStream.open(famfile.c_str());
  for(i=0;i<n;i++){
    for(k=0;k<nColFamFile;k++){
      tmpStream >> tok;
      if(k==0){
        data.FID[i] = tok;
      }
      if(k==1){
        data.IID[i] = tok;
      }    
    }
  }
  tmpStream.close();
  if(verbose){
    cout<<"done"<<endl;
  }
  
  //Get SNPs ID
  if(verbose){
    cout<<"5-: Getting SNPs IDs...";
  }
  int nColBimFile = 6;  
  tmpStream.open(bimfile.c_str());
  for(i=0;i<p;i++){
    for(k=0;k<nColBimFile;k++){
      tmpStream >> tok;
      if(k==1){
        data.info_loci[i] = tok;
      }  
    }
  }
  tmpStream.close();
  if(verbose){
    cout<<"done"<<endl;
  }
  //Get LD score (ASSUMING THAT THE FILE CONTAINS THE SAME SNPS THAN THE BED FILE AND IN THE SAME ORDER)
  if(LDfile != "none"){
	if(verbose){
		cout<<"6-: Getting LD Scores...";
	}
	double tmpfrq = 0;
	double ldscore = 0;
	  
	tmpStream.open(LDfile.c_str());
	for(i=0;i<p;i++){
		tmpStream >> tok >> tmpfrq >> ldscore;
		if(tmpfrq>=MinMaf & tmpfrq<=MaxMaf & ldscore>=MinLD & ldscore<=MaxLD){data.IndexAscert.push_back(i); }
	}
  cout<<"Number of SNPs satisfying the conditions for ascertained sampling = "<<data.IndexAscert.size()<<endl;
  tmpStream.close();
  }
  if(verbose){
    cout<<"done"<<endl;
  }
}



void decode_plink(char *output, const char *input, const int lengthInput){
  int i, k;
  char tmp, geno;
  int a1, a2;
  
  for(i=0;i<lengthInput;++i){
    tmp = input[i];
    k   = PACK_DENSITY * i;
    geno      = (tmp & MASK0);
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK1) >> 2; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK2) >> 4; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK3) >> 6; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
  }
}


int readGenotypes(data_struct& data,
                  string bedfile, long int n, long int p,
                  char *packed, char *unpacked, int numBytes){
  ifstream influx;
  influx.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if(!influx){
    cerr << "[readGenotypes] Error reading file "<<bedfile<<endl;
    exit(1);
  }
  influx.seekg(0, ifstream::end);
  influx.seekg(3, ifstream::beg);
  cout<<"Start reading genotypes";
  long int i,j;
  double x;
  int ndropped = 0;
  int idx=0;
  vector<int> missingGeno;   
  for(j=0;j<p;j++){
	if(j==round(p/3) or j==(2*round(p/3))){cout<<".";}
    influx.read((char*)packed, sizeof(char) * numBytes);
    decode_plink(unpacked, packed, numBytes);
    if(data.index_causal[j]){
      data.freq[idx] = 0.0;
      for(i=0;i<n;i++){
        x = (double)((int) unpacked[i]);
        if(x==3.0){
          data.geno[i][idx] = (double) -1.0;
          missingGeno.push_back(i); 
        }else{
          data.geno[i][idx] = x;
          data.freq[idx]+= data.geno[i][idx];
        }
      }
      data.freq[idx]  = (double) 0.5 * data.freq[idx]/(n-missingGeno.size());
      missingGeno.clear();
      if(data.freq[idx] == 0.0){
        ndropped++;
      }
      idx++;
    }
  }
  cout<<".done"<<endl;
  influx.close();
  return ndropped;
}
