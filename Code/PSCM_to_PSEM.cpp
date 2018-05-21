//AUTHOR: Helge Roider
//INSTITUTION: Max Planck Institute for Molecular Genetics - Berlin
//DATE: 20/11/2006

//takes a count matrix in TRANSFAC format, adds pseudocounts and removes uninformative positions
//and then converts the modified count matrix into an energy matrix

//Modified by Marcel Schulz and Florian Schmidt
//Cluster of Excellence, Saarland University, 2016

#include<stdio.h>
#include<cmath>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<string.h>
#include<iomanip>

using namespace std;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;


int main(int argc, char *argv[]){

  //--------------------------------------------------------------
  //PARAMETERS
  //--------------------------------------------------------------

  double lambda = 0.7;
  //double slope = 0.66; // slope * effectivelength[f]
  //double intercept = -6.6;
  double slope = 0.584; 
  double intercept = -5.66;

  double lnR0; //effectifelength * slope + intercept

  //GLOBAL GENOMIC BACKGROUND
  const double gc_content = 0.43;
  const double at_content = 1 - gc_content;

  cout << "#/Lambda=" << lambda << "\t/Regression_slope=" << slope << "\t/Regression_intercept=" << intercept << "\t/GC_content=" << gc_content << "\t/Pseudocount=maximal_row_counts/counts_in_row\n";

  //---------------------------------------------------------------
  //VARIABLES
  //---------------------------------------------------------------

  const int numoffactors = 16000;
  const int numofpositions = 50;

  double *** pwm;
  double *** complement;
  pwm = new double ** [numoffactors];
  complement = new double ** [numoffactors];
  for(int i = 0; i < numoffactors; i++){
    pwm[i] = new double * [numofpositions];
    complement[i] = new double * [numofpositions];
    for(int j = 0; j < numofpositions; j++){
      pwm[i][j] = new double[4];
      complement[i][j] = new double[4];
    }
  }

  //  double pwm[numoffactors][numofpositions][4];    // energy matrix forward strand
  //double complement[numoffactors][numofpositions][4]; // energy matrix reverse strand
  int corelength[numoffactors];
  int factors = -1;
  string factornames[numoffactors];


  //----------------------------------------------------------------
  //READ TRANSFAC FILE
  //----------------------------------------------------------------

  ifstream transfac(argv[1]);

  if(!transfac){
    cout << "TRANSFAC file not opened\n";
    exit(1);
  }

  //matrix variables

  double ** entropy;
  double *** matrix;
  matrix = new double ** [numoffactors];
  entropy = new double * [numoffactors];
  for(int i = 0; i < numoffactors; i++){
    matrix[i] = new double * [numofpositions];
    entropy[i] = new double [numofpositions];
    for(int j = 0; j < numofpositions; j++){
      matrix[i][j] = new double[4];
    }
  }

  //  double matrix[numoffactors][numofpositions][4]; // read in matrix
  //double entropy[numoffactors][numofpositions];   // entropy of positions
  double Pseudocount; //count dependent
  
  int position;
  int motiflength[numoffactors];
  int effectivelength;
  
  double maxcount[numoffactors];
  double rowsum[numoffactors][numofpositions];
  
  //reading file variables
  string word[5]; //elements in row
  double max = 0; //consensus base count
  string row; //transfac file rows
  string delimiters = " \t"; //word seperators in each line
  int reading;
  int start, end;
  
  while(!transfac.eof()){
    getline(transfac,row);
    start = row.find_first_not_of(delimiters);
    int i = 0;
    
    while((start != string::npos)&&(i<5)){ //split row into tokens - word[]
      end = row.find_first_of(delimiters,start+1);
      if(end == string::npos){
	end = row.length();
      }
      word[i] = row.substr(start,end-start);
      i++;
      start = row.find_first_not_of(delimiters,end+1);
    }
    
    if((reading == 1)&&(word[0] == "XX")){ //end of matrix is reached
      reading = 0;
    }
    
    if(reading == 1){ //generate matrix
      matrix[factors][position][A] = strtod(word[1].c_str(),NULL);
      matrix[factors][position][C] = strtod(word[2].c_str(),NULL);
      matrix[factors][position][G] = strtod(word[3].c_str(),NULL);
      matrix[factors][position][T] = strtod(word[4].c_str(),NULL);
      
      rowsum[factors][position] = strtod(word[1].c_str(),NULL)+strtod(word[2].c_str(),NULL)+strtod(word[3].c_str(),NULL)+strtod(word[4].c_str(),NULL);
      
      if(rowsum[factors][position] > maxcount[factors]){
	maxcount[factors] = rowsum[factors][position];
      }
      motiflength[factors]++;
      position++;
    }
    
    if(word[0] == "P0"){ //start of matrix
      factors++;
      reading = 1;
      position = 0;
      motiflength[factors] = 0;
      maxcount[factors] = 0;
    }
    
    if(word[0] == "ID"){
      factornames[factors + 1] = word[1]+"\t"+word[2];
    }
  }
  
  transfac.close();
  
  
  //------------------------------------------------------------
  //MAKE FORWARD AND REVERSE MATRIX
  //------------------------------------------------------------
  
  for(int f = 0; f <= factors; f++){
    
    for(int p = 0; p < motiflength[f]; p++){
      
      //      Pseudocount = maxcount[f] / rowsum[f][p];
      Pseudocount = 1;

      matrix[f][p][A] = matrix[f][p][A] + Pseudocount;
      matrix[f][p][C] = matrix[f][p][C] + Pseudocount;
      matrix[f][p][G] = matrix[f][p][G] + Pseudocount;
      matrix[f][p][T] = matrix[f][p][T] + Pseudocount;
      
      entropy[f][p] = 0;
      double prob, de;
      for(int b = 0; b < 4; b++){
	prob = matrix[f][p][b] / (rowsum[f][p] + 4 * Pseudocount);
	if(prob > 0){de = prob * log(4 * prob) / log(2.00);}
	else{de = 0;}
	entropy[f][p] = entropy[f][p] + de;
      }
      
      double maxAT, maxCG;
      
      if(matrix[f][p][A] > matrix[f][p][T]){maxAT = matrix[f][p][A];}
      else{maxAT = matrix[f][p][T];}
      if(matrix[f][p][C] > matrix[f][p][G]){maxCG = matrix[f][p][C];}
      else{maxCG = matrix[f][p][G];}
      
      if(maxAT > maxCG){
	matrix[f][p][A] = log(maxAT/matrix[f][p][A]) / lambda;
	matrix[f][p][T] = log(maxAT/matrix[f][p][T]) / lambda;
	matrix[f][p][C] = log((maxAT/at_content)*(gc_content/matrix[f][p][C])) / lambda;
	matrix[f][p][G] = log((maxAT/at_content)*(gc_content/matrix[f][p][G])) / lambda;
      }
      if(maxAT < maxCG){
	matrix[f][p][C] = log(maxCG/matrix[f][p][C]) / lambda;
	matrix[f][p][G] = log(maxCG/matrix[f][p][G]) / lambda;
	matrix[f][p][A] = log((maxCG/gc_content)*(at_content/matrix[f][p][A])) / lambda;
	matrix[f][p][T] = log((maxCG/gc_content)*(at_content/matrix[f][p][T])) / lambda;
      }
      if(maxAT == maxCG){
	matrix[f][p][A] = log(maxAT/matrix[f][p][A]) / lambda;
	matrix[f][p][C] = log(maxAT/matrix[f][p][C]) / lambda;
	matrix[f][p][G] = log(maxAT/matrix[f][p][G]) / lambda;
	matrix[f][p][T] = log(maxAT/matrix[f][p][T]) / lambda;
      }
    }
    
    
    
    //DETERMINE INFORMATIVE CORE
    int n = 0;
    int core = 0;
    int pos = 0;
    corelength[f] = 0;
    
    for(int p = 0; p < motiflength[f]; p++){
	pwm[f][n][A] = matrix[f][p][A];
	pwm[f][n][C] = matrix[f][p][C];
	pwm[f][n][G] = matrix[f][p][G];
	pwm[f][n][T] = matrix[f][p][T];
	corelength[f]++;
	n++;
      }
    
    n = motiflength[f] - 1;
    
    effectivelength=corelength[f];
	lnR0 = effectivelength * slope + intercept;
    //    if(lnR0 < 1e-5){lnR0 = 0;}

    //COMPLEMENT MATRIX

    if(corelength[f] > 0){
      cout << ">" << factornames[f] << "\tlnR0: " << lnR0 << "\n"; 
      for(int p = 0; p < corelength[f]; p++){
	cout << pwm[f][p][A] << "\t" << pwm[f][p][C] << "\t" << pwm[f][p][G] << "\t" << pwm[f][p][T] << "\n";
      }
    }

  }//loop over factors

  for(int i = 0; i < numoffactors; i++){
    for(int j = 0; j < numofpositions; j++){
      delete [] matrix[i][j];
      delete [] complement[i][j];
      delete [] pwm[i][j];
    }
    delete [] matrix[i];
    delete [] entropy[i];
    delete [] complement[i];
    delete [] pwm[i];
  }
  delete [] matrix;
  delete [] entropy;
  delete [] complement;
  delete [] pwm;

}
