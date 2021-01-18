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
#include "write.h"

void Write_Pheno(data_struct data, string outfile){
  long int i;
  //Write Phenotype
  ofstream fileOut(outfile.c_str());
  fileOut<<"FID\tIID\tPheno\tPhenoScaled"<<endl;
  for(i=0;i<data.nbr_ind;i++){
    fileOut<<data.FID[i]<<"\t"<<data.IID[i]<<"\t"<<data.Pheno[i]<<"\t"<<data.PhenoScaled[i]<<endl;
  }
  fileOut.close();
  
}

void Write_CausalInfo(data_struct data, string outfile){
  long int j;
  int idx=0;
  //Write Causals Id and effects
  ofstream fileOut(outfile.c_str());
  fileOut<<"SNP_ID\tAll_freq\tAdditive_effect\tDominance_effect\t"<<endl;
  for(j=0;j<data.nbr_loci;j++){
	if(data.index_causal[j]){
		fileOut<<data.info_loci[j]<<"\t"<<data.freq[idx]<<"\t"<<data.SNPeffect[idx][0]<<"\t"<<data.SNPeffect[idx][1]<<endl;
		idx++;
	}  
  }
  fileOut.close();
  
}

/*
 This function is part of the SelEstim software code: Vitalis R, Gautier M, Dawson KJ and Beaumont MA (2014) Detecting and measuring selection from gene frequency data. Genetics 196: 799-817
*/

char *print_time(double time)

{
  int total;
  int seconds = 0;
  int minutes = 0;
  int hours = 0;
  int days = 0;
  int years = 0;
  char number[8] = "\0";
  static char formatted_time[256] = "\0";
  
  total = (int) time;
  if (total < 60) {
    seconds = total;
    sprintf(number,"%d",seconds);
    strcpy(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else if (total < 3600) {
    minutes = floor(total / 60);
    seconds = total % 60;
    sprintf(number,"%d",minutes);
    strcpy(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else if (total < 86400) {
    hours = floor(total / 3600);
    minutes = floor((total - (hours * 3600)) / 60);
    seconds = total % 60;
    sprintf(number,"%d",hours);
    strcpy(formatted_time,number);
    strcat(formatted_time," hr");
    (hours > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",minutes);
    strcat(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else if (total < 31536000) {
    days = floor(total / 86400);
    hours = floor((total - (days * 86400)) / 3600);
    minutes = floor((total - (days * 86400) - (hours * 3600)) / 60);
    seconds = total % 60;
    sprintf(number,"%d",days);
    strcpy(formatted_time,number);
    strcat(formatted_time," day");
    (days > 1) ? strcat(formatted_time,"s ") : strcat(formatted_time," ");
    sprintf(number,"%d",hours);
    strcat(formatted_time,number);
    strcat(formatted_time," hr");
    (hours > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",minutes);
    strcat(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else {
    years = floor(total / 31536000);
    days = floor((total - (years * 31536000)) / 86400);
    hours = floor((total - (years * 31536000) - (days * 86400)) / 3600);
    minutes = floor((total - (years * 31536000) - (days * 86400) - (hours * 3600)) / 60);
    seconds = total % 60;
    sprintf(number,"%d",years);
    strcpy(formatted_time,number);
    strcat(formatted_time," year");
    (years > 1) ? strcat(formatted_time,"s ") : strcat(formatted_time," ");
    sprintf(number,"%d",days);
    strcat(formatted_time,number);
    strcat(formatted_time," day");
    (days > 1) ? strcat(formatted_time,"s ") : strcat(formatted_time," ");
    sprintf(number,"%d",hours);
    strcat(formatted_time,number);
    strcat(formatted_time," hr");
    (hours > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",minutes);
    strcat(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.\0") : strcat(formatted_time,".");
  }
  strcat(formatted_time,"\0");
  return(formatted_time);
}
