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
#include "simul.h"

void ComputeGeneticValues_omp(data_struct &data, long int nbr_loci, double h, int model)

{	
	int i,j;
	double var,mean;
	mean=0.0;
	var=0.0;
	long int nbr_ind=data.nbr_ind;
	
	if(model == 1){ //ComputeGenetic values for additive model
	
		double *coeff = new double[nbr_loci];
		for (j=0; j< nbr_loci; ++j) {
			coeff[j] = norm(gen);
		}
	
		#pragma omp parallel for schedule(guided) private(j)
		for (i = 0; i < nbr_ind; ++i) {
			for (j=0; j< nbr_loci; ++j) {
				if( data.geno[i][j]!=(-1.0) ){
					data.genvalue[i] += (data.geno[i][j]-2.0*data.freq[j]) / sqrt(2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
				}
			}
		}
		
		for (i = 0; i < nbr_ind; ++i) {
			mean+=data.genvalue[i];
		}
		mean=mean/nbr_ind;
		for (i = 0; i < nbr_ind; ++i) {
			var+=pow(data.genvalue[i]-mean,2);
		}
		var=var/nbr_ind;
		//Adjust SNP effect to reach the expected heritability
		#pragma omp parallel for schedule(guided)
		for (j=0; j< nbr_loci; ++j) {
			coeff[j] = coeff[j] * sqrt(h/var);
			data.SNPeffect[j][0]=coeff[j];
		}
	
		#pragma omp parallel for schedule(guided) private(j)
		for (i = 0; i < nbr_ind; ++i) {
			for (j=0; j< nbr_loci; ++j) {
				
					if(j==0){
						if( data.geno[i][j]!=(-1.0) ){
							data.genvalue[i] = (data.geno[i][j]-2.0*data.freq[j]) / sqrt(2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
						} else {data.genvalue[i] = 0;}
					} else {
						if( data.geno[i][j]!=(-1.0) ){
							data.genvalue[i] += (data.geno[i][j]-2.0*data.freq[j]) / sqrt(2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
						}
					}
					
			}
		}
		 
		///////////////////////////////////////////////////////////////////////////////
			mean=0;
			var=0;
			for (i = 0; i < nbr_ind; ++i) {
				mean+=data.genvalue[i];
			}
			mean=mean/nbr_ind;
			for (i = 0; i < nbr_ind; ++i) {
			var+=pow(data.genvalue[i]-mean,2);
		}
			var=var/nbr_ind;
			cout<<"variance Additive effects = "<<var<<endl;
		//////////////////////////////////////////////////////////////////////////////
		#pragma omp parallel for schedule(guided)
		for (i = 0; i < nbr_ind; ++i) {
					data.Pheno[i]+=data.genvalue[i];
					data.genvalue[i]=0.0;
			}	  
		delete [] coeff;
	} else if(model == 2){ //Dominance model
		
		
		double *coeff = new double[nbr_loci];
		for (j=0; j< nbr_loci; ++j) {
			coeff[j] = norm(gen);
		}
	
		#pragma omp parallel for schedule(guided) private(j)
		for (i = 0; i < nbr_ind; ++i) {
			for (j=0; j< nbr_loci; ++j) {
				if(data.geno[i][j] == 0.0){
					data.genvalue[i] += (0-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
				}else if(data.geno[i][j] == 1.0){
					data.genvalue[i] += ((2.0*data.freq[j])-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
				} else if(data.geno[i][j]==2.0){
					data.genvalue[i] += ((4.0*data.freq[j]-2.0)-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
				} else if(data.geno[i][j]!=(-1.0)){
					cerr<<"WARNING : One genotype is different from -1 (missing), 0, 1 or 2."<<endl;
					exit(1);
				}
			}
		}
		
		for (i = 0; i < nbr_ind; ++i) {
			mean+=data.genvalue[i];
		}
		mean=mean/nbr_ind;
		for (i = 0; i < nbr_ind; ++i) {
			var+=pow(data.genvalue[i]-mean,2);
		}
		var=var/nbr_ind;
		//Adjust SNP effect to reach the expected heritability
		#pragma omp parallel for schedule(guided)
		for (j=0; j< nbr_loci; ++j) {
			coeff[j] = coeff[j] * sqrt(h/var);
			data.SNPeffect[j][1]=coeff[j];
		}
	
		#pragma omp parallel for schedule(guided) private(j)
		for (i = 0; i < nbr_ind; ++i) {
			for (j=0; j< nbr_loci; ++j) {
				if(j==0){
					if(data.geno[i][j] == 0.0){
						data.genvalue[i] = (0-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
					}else if(data.geno[i][j] == 1.0){
						data.genvalue[i] = ((2.0*data.freq[j])-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
					} else if(data.geno[i][j]==2.0){
						data.genvalue[i] = ((4.0*data.freq[j]-2.0)-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
					} else if(data.geno[i][j]==(-1.0)){
						data.genvalue[i] = 0;
					}
				} else {
					if(data.geno[i][j] == 0.0){
						data.genvalue[i] += (0-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
					}else if(data.geno[i][j] == 1.0){
						data.genvalue[i] += ((2.0*data.freq[j])-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
					} else if(data.geno[i][j]==2.0){
						data.genvalue[i] += ((4.0*data.freq[j]-2.0)-2.0*pow(data.freq[j],2)) / (2.0 * data.freq[j]*(1.0-data.freq[j])) * coeff[j];
					}
				}	
			}
		} 
		
		///////////////////////////////////////////////////////////////////////////////
			mean=0;
			var=0;
			for (i = 0; i < nbr_ind; ++i) {
				mean+=data.genvalue[i];
			}
			mean=mean/nbr_ind;
			for (i = 0; i < nbr_ind; ++i) {
			var+=pow(data.genvalue[i]-mean,2);
		}
			var=var/nbr_ind;
			cout<<"variance dominance effects = "<<var<<endl;
		//////////////////////////////////////////////////////////////////////////////
		
		
		#pragma omp parallel for schedule(guided)
		for (i = 0; i < nbr_ind; ++i) {
					data.Pheno[i]+=data.genvalue[i];
					data.genvalue[i]=0.0;
			}	  
		
		
		delete [] coeff;
		
	} else if (model == 3){	//Additive-by-additive model
		long int k, idx;
		long int ninter=((pow(nbr_loci,2)-nbr_loci)/2);
		double tmpg1, tmpg2;
		k=0;
		idx=0;
		tmpg1=0.0;
		tmpg2=0.0;
		double *coeff = new double[ninter];
		for (j=0; j< ninter; ++j) {
			coeff[j] = norm(gen);
		}
		#pragma omp parallel for schedule(guided) private(j,idx,k,tmpg1,tmpg2)
		for (i = 0; i < nbr_ind; ++i) {
				idx=0;
				for (j=0; j< (nbr_loci-1); ++j) {
					if(data.geno[i][j]!=(-1.0)){
						tmpg1 = (data.geno[i][j]-2.0*data.freq[j]) / sqrt(2.0 * data.freq[j]*(1.0-data.freq[j]));
					} else {tmpg1 = 0;}
						for (k=(j+1); k< nbr_loci; ++k) {
							if(data.geno[i][k]!=(-1.0)){
								tmpg2 = (data.geno[i][k]-2.0*data.freq[k]) / sqrt(2.0 * data.freq[k]*(1.0-data.freq[k]));
								data.genvalue[i] += tmpg1 * tmpg2 * coeff[idx];
							}
							idx++;
						}
				}
			}
		
		for (i = 0; i < nbr_ind; ++i) {
			mean+=data.genvalue[i];
		}
		mean=mean/nbr_ind;
		
		
		for (i = 0; i < nbr_ind; ++i) {
			var+=pow(data.genvalue[i]-mean,2);
		}
		var=var/nbr_ind;
		//Adjust SNP effect to reach the expected heritability
		#pragma omp parallel for schedule(guided)
		for (j=0; j<ninter; ++j) {
				coeff[j] = coeff[j] * sqrt(h/var);
			}
		
		#pragma omp parallel for schedule(guided) private(j,k,idx,tmpg1,tmpg2)
		for (i = 0; i < nbr_ind; ++i) {
				idx=0;
				for (j=0; j< (nbr_loci-1); ++j) {
					if(data.geno[i][j]!=(-1.0)){
						tmpg1 = (data.geno[i][j]-2.0*data.freq[j]) / sqrt(2.0 * data.freq[j]*(1.0-data.freq[j]));
					} else {tmpg1 = 0.0;}
						for (k=(j+1); k< nbr_loci; ++k) {
								if(idx==0){
									if(data.geno[i][k]!=(-1.0)){
										tmpg2 = (data.geno[i][k]-2.0*data.freq[k]) / sqrt(2.0 * data.freq[k]*(1.0-data.freq[k]));
										data.genvalue[i] = tmpg1 * tmpg2 * coeff[idx];
									} else {
										data.genvalue[i] = 0;
										}
								}else{
									if(data.geno[i][j]!=(-1.0) && data.geno[i][k]!=(-1.0)){
										tmpg2 = (data.geno[i][k]-2.0*data.freq[k]) / sqrt(2.0 * data.freq[k]*(1.0-data.freq[k]));
										data.genvalue[i] += tmpg1 * tmpg2 * coeff[idx];
									}
								}
								idx++;
							}	
				}
			}
			
		///////////////////////////////////////////////////////////////////////////////
			mean=0;
			var=0;
			for (i = 0; i < nbr_ind; ++i) {
				mean+=data.genvalue[i];
			}
			mean=mean/nbr_ind;
			for (i = 0; i < nbr_ind; ++i) {
			var+=pow(data.genvalue[i]-mean,2);
		}
			var=var/nbr_ind;
			cout<<"variance additive-by-additive effects = "<<var<<endl;
		//////////////////////////////////////////////////////////////////////////////
			
		#pragma omp parallel for schedule(guided)
		for (i = 0; i < nbr_ind; ++i) {
					data.Pheno[i]+=data.genvalue[i];
					data.genvalue[i]=0.0;
			} 
		  delete [] coeff;
	}
	
}

void BuildPheno(data_struct &data, double H, bool verbose){
  long int i;
  double Ve = 1.0-(H);
  double mean=0;
  double sd =0;
  double var=0;
  vector<double> residual (data.nbr_ind); 
  for(i=0;i<data.nbr_ind;i++){
	residual[i]=norm(gen);  
    mean += residual[i];
  }
  mean/=data.nbr_ind;
  for (i = 0; i < data.nbr_ind; ++i) {
	var+=pow(residual[i]-mean,2);
  }
  var/=data.nbr_ind;
  mean=0;
  for(i=0;i<data.nbr_ind;i++){
	data.Pheno[i] += (residual[i] * sqrt(Ve/var));
	mean+=data.Pheno[i];
  }
  mean/=data.nbr_ind;
  var=0;
  for (i = 0; i < data.nbr_ind; ++i) {
		sd += pow(data.Pheno[i]-mean,2);
	}
  if(verbose){cout<<"var pheno = "<<sd/data.nbr_ind<<endl;}
  sd=sqrt(sd/data.nbr_ind);
 
	
  for (i = 0; i < data.nbr_ind; ++i) {
	data.PhenoScaled[i] = (data.Pheno[i] - mean) / sd;
  }
  
  

}
