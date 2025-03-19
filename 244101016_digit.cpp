// 244101016_digit.cpp : Defines the entry point for the console application.
//


#include "StdAfx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include<fstream>
#include <iomanip>
#include <cmath>
#include "string"
#include "string"
#include "iomanip"
#include "cmath"
#include "conio.h"
#include <Windows.h>
using namespace std;

#define LPC_ORDER 12  // Order of LPC (Linear Predictive Coding)
#define PI 3.1428  // Value of PI
#define FRAME_SIZE 320  // Number of samples per frame

int TOTAL_FRAMES;  // Number of frames to process

// Global arrays to store signal data and computed values
int totalSamples;
double rawSignalData[65][FRAME_SIZE];  // Signal data for each frame
double autocorrelation[65][LPC_ORDER + 1];  // Autocorrelation coefficients
double lpcCoefficients[65][LPC_ORDER + 1];  // LPC coefficients
double cepstralCoefficients[65][LPC_ORDER + 1];  // Cepstral coefficients
double rawSignal[22000] = {0.0};  // Raw signal data
 double universe[15207][12]; 
//for copying to universe
//#define VALUES_PER_ROW 12
#define MAX_LINE_LENGTH 1024
#define TOTAL_FILES 400
#define FILENAME_SIZE 30


 

// Function to read raw signal data from a file
void ReadSignalDataFromFile(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Skip the first 10 lines of the file
    char tempStr[150];
    for(int skip = 0; skip < 10; skip++){
        fscanf(file, "%s", tempStr);
    }

    // Read the signal data
    while (!feof(file)) {
        if (fscanf(file, "%lf", &rawSignal[totalSamples]) && totalSamples < 22000) {
            totalSamples++;
        } else {
            break;
        }
    }

    fclose(file);
}

// Function to correct DC offset in the raw signal
void CorrectDCOffsetInSignal() {
    // Calculate the mean of the signal
    double mean = 0.0;
    for (int i = 0; i < totalSamples; i++) {
        mean += rawSignal[i];
    }
    mean /= totalSamples;

    // Subtract the mean from each sample to correct DC offset
    for (int i = 0; i < totalSamples; i++) {
        rawSignal[i] -= mean;
    }
}

// Function to apply Hamming window to a signal frame
void ApplyHammingWindowToFrame(int frameIndex) {
    // Apply the Hamming window function to each sample in the frame
    for (int sample = 0; sample < FRAME_SIZE; sample++) {
        rawSignalData[frameIndex][sample] *= 0.54 - 0.46 * cos(2 * PI * sample / (FRAME_SIZE - 1));
    }
}




// Function to normalize signal values to a specified range
void NormalizeSignalValues() {
    // Find the minimum and maximum values in the signal
    double minVal = rawSignal[0];
    double maxVal = rawSignal[0];

    for (int i = 1; i < totalSamples; i++) {
        if (rawSignal[i] < minVal) {
            minVal = rawSignal[i];
        }
        if (rawSignal[i] > maxVal) {
            maxVal = rawSignal[i];
        }
    }

    // Normalize the signal to the range [-5000, 5000]
    for (int i = 0; i < totalSamples; i++) {
        rawSignal[i] = -5000 + (((rawSignal[i] - minVal) / (maxVal - minVal)) * 10000);
    }
}


// Function to select steady state frames from the raw signal
void SelectSteadyStateFramesFromSignal() {
    for (int frame = 0; frame < TOTAL_FRAMES; frame++) {
        for (int sample = 0; sample < FRAME_SIZE; sample++) {
            rawSignalData[frame][sample] = rawSignal[frame*320 + sample];
        }
    }
}

// Function to compute autocorrelation coefficients for a frame
void ComputeAutoCorrelationForFrame(int frameIndex) {
    // Compute autocorrelation for lags from 0 to LPC_ORDER
    for (int lag = 0; lag <= LPC_ORDER; lag++) {
        autocorrelation[frameIndex][lag] = 0.0;
        // Calculate autocorrelation value for the given lag
        for (int sample = 0; sample < FRAME_SIZE - lag; sample++) {
            autocorrelation[frameIndex][lag] += rawSignalData[frameIndex][sample] * rawSignalData[frameIndex][sample + lag];
        }
    }
}

// Function to compute LPC coefficients using Levinson-Durbin recursion for a frame
void ComputeLPCoefficientsForFrame(int frameIndex) {
    double error[LPC_ORDER + 1] = {0.0};  // Prediction error values
    double reflectionCoeffs[LPC_ORDER + 1] = {0.0};  // Reflection coefficients
    double predictionErrors[LPC_ORDER + 1][LPC_ORDER + 1] = {0.0};  // Prediction errors matrix

    // Initialize the error for lag 0
    error[0] = autocorrelation[frameIndex][0];

    // Levinson-Durbin recursion to compute LPC coefficients
    for (int i = 1; i <= LPC_ORDER; i++) {
        double sum = 0.0;

        // Compute reflection coefficient
        for (int j = 1; j < i; j++) {
            sum += predictionErrors[j][i - 1] * autocorrelation[frameIndex][i - j];
        }
        reflectionCoeffs[i] = (autocorrelation[frameIndex][i] - sum) / error[i - 1];

        // Update prediction errors
        predictionErrors[i][i] = reflectionCoeffs[i];
        for (int j = 1; j < i; j++) {
            predictionErrors[j][i] = predictionErrors[j][i - 1] - reflectionCoeffs[i] * predictionErrors[i - j][i - 1];
        }

        // Update error value
        error[i] = (1 - reflectionCoeffs[i] * reflectionCoeffs[i]) * error[i - 1];
    }

    // Store the LPC coefficients
    for (int i = 1; i <= LPC_ORDER; i++) {
        lpcCoefficients[frameIndex][i] = predictionErrors[i][LPC_ORDER];
    }
}

// Function to compute Cepstral Coefficients from LPC coefficients for a frame
void ComputeCepstralCoefficientsForFrame(int frameIndex) {
    cepstralCoefficients[frameIndex][0] = 0.0;
   // weightedCepstralCoeffs[frameIndex][0] = 0.0;

    // Compute cepstral coefficients
    for (int n = 1; n <= LPC_ORDER; n++) {
        cepstralCoefficients[frameIndex][n] = lpcCoefficients[frameIndex][n];
        double sum = 0.0;
        for (int k = 1; k < n; k++) {
            if (n - k >= 0) {
                sum += k * lpcCoefficients[frameIndex][n - k] * cepstralCoefficients[frameIndex][k];
            }
        }
        cepstralCoefficients[frameIndex][n] += sum / n;
       // weightedCepstralCoeffs[frameIndex][n] += cepstralCoefficients[frameIndex][n] * (1 + (LPC_ORDER / 2) * sin(PI * n / LPC_ORDER));
    }
}
// Function to write weighted cepstral coefficients to a file
void WriteWeightedCepstralCoeffsToFile(const char* filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Write cepstral coefficients for each frame to the file
    for (int frame = 0; frame < TOTAL_FRAMES; frame++) {
        for (int coef = 1; coef <= LPC_ORDER; coef++) {
           // cepstralCoefficients[frame][coef] /= 20;  // Normalize coefficients
            fprintf(file, "%lf ", cepstralCoefficients[frame][coef]);
            cepstralCoefficients[frame][coef] = 0;  // Reset for next use
        }
        fprintf(file, "\n");
    }

    fclose(file);
}


void mergeFiles(const char *outputFile) {
    FILE *outFile = fopen(outputFile, "w");
    if (outFile == NULL) {
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }

    char fileName[100]="";
    for (int i = 1; i <= TOTAL_FILES; i++) {
        // Generate each file name dynamically
        sprintf(fileName, "./OUTPUT_CEPSTRAL/244101016_%d.txt", i);

        FILE *inFile = fopen(fileName, "r");
        if (inFile == NULL) {
            perror("Error opening input file");
            fclose(outFile);
            exit(EXIT_FAILURE);
        }

        char line[MAX_LINE_LENGTH];
        while (fgets(line, sizeof(line), inFile) != NULL) {
            fputs(line, outFile);  // Append each line to the output file
        }

        fclose(inFile);
    }

    fclose(outFile);
    printf("Files merged successfully into '%s'.\n", outputFile);
}




const int K = 32;              // Codebook size
const int DIMENSION = 12;      // Vector dimensions
const double DELTA = 0.00001;  // Convergence threshold
const int MAX_VECTORS = 15207;  // Maximum number of vectors
const double EPSILON = 0.03;   // Epsilon for centroid splitting
double codebook[K][DIMENSION];
// Function to load the universe vectors from a CSV file
void store_all_value_of_universe(const char *filename, double universe[][DIMENSION]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file.\n");
        exit(1);
    }
    int k = 0;
    while (!feof(file)) {
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
               &universe[k][0], &universe[k][1], &universe[k][2], &universe[k][3],
               &universe[k][4], &universe[k][5], &universe[k][6], &universe[k][7],
               &universe[k][8], &universe[k][9], &universe[k][10], &universe[k][11]);
        k++;
    }
	/*  while (!feof(file)) {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,",
               &universe[k][0], &universe[k][1], &universe[k][2], &universe[k][3],
               &universe[k][4], &universe[k][5], &universe[k][6], &universe[k][7],
               &universe[k][8], &universe[k][9], &universe[k][10], &universe[k][11]);
        k++;
    }*/
    fclose(file);
}

// Function to initialize the codebook with the centroid of the universe (single vector)
void initialize_codebook(double codebook[1][DIMENSION], double universe[][DIMENSION], int M) {
    for (int i = 0; i < DIMENSION; i++) {
        double sum = 0;
        for (int j = 0; j < M; j++) {
            sum += universe[j][i];
        }
        codebook[0][i] = sum / M;
    }
}

// Function to double the codebook by splitting each centroid
void split_codebook(double codebook[][DIMENSION], int *current_size) {
    int new_size = *current_size * 2;
    for (int i = 0; i < *current_size; i++) {
        // Create the new centroid by modifying the current one
        for (int j = 0; j < DIMENSION; j++) {
            codebook[i + *current_size][j] = codebook[i][j] * (1 + EPSILON);  // Split with + epsilon
            codebook[i][j] = codebook[i][j] * (1 - EPSILON);  // Split with - epsilon
        }
    }
    *current_size = new_size;
}

// Function to calculate the Euclidean distance between two vectors
double euclidean_distance(double vec1[DIMENSION], double vec2[DIMENSION]) {
    double dist = 0;
    for (int i = 0; i < DIMENSION; i++) {
        dist += pow(vec1[i] - vec2[i], 2);
    }
    return sqrt(dist);
}

// Function to calculate the Tokhura distance between two vectors
double tokhura_distance(double vec1[DIMENSION], double vec2[DIMENSION]) {
    double weights[DIMENSION] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0}; 
    double dist = 0;
    for (int i = 0; i < DIMENSION; i++) {
        dist += weights[i] * pow(vec1[i] - vec2[i], 2);
    }
    return dist;
}

// Function to assign each vector to the nearest centroid and calculate distortion
void form_cluster(double universe[][DIMENSION], int M, double codebook[][DIMENSION], int current_size, int region[], double *dist) {
    *dist = 0;
    for (int i = 0; i < M; i++) {
        double min_euc_dist = DBL_MAX;
        int nearest_centroid = -1;

        // Find the nearest centroid using Euclidean distance
        for (int j = 0; j < current_size; j++) {
            double euc_dist = euclidean_distance(universe[i], codebook[j]);
            if (euc_dist < min_euc_dist) {
                min_euc_dist = euc_dist;
                nearest_centroid = j;
            }
        }

        // Assign the vector to the nearest centroid
        region[i] = nearest_centroid;

        // Calculate the distortion using Tokhura distance (for the assigned centroid)
        *dist += tokhura_distance(universe[i], codebook[nearest_centroid]);
    }
}

// Function to update the centroids based on the assigned vectors
void update_centroids(double universe[][DIMENSION], int M, double codebook[][DIMENSION], int current_size, int region[]) {
	double store_sum[32][DIMENSION]={0};
	int current_cluster_size[32]={0};

    // Sum vectors assigned to each region
    for (int i = 0; i < M; i++) {
        int reg = region[i];
        for (int j = 0; j < DIMENSION; j++) {
            store_sum[reg][j] += universe[i][j];
        }
        current_cluster_size[reg]++;
    }

    // Compute the new centroids
    for (int i = 0; i < current_size; i++) {
        if (current_cluster_size[i] > 0) {
            for (int j = 0; j < DIMENSION; j++) {
                codebook[i][j] = store_sum[i][j] / current_cluster_size[i];
            }
        }
    }
}


///////////////////////////////////////////////////////////step 3
double framestore[100][12];
int frameonedigit;
int ansstore[100];
void readfileforstep3(const char* filename3){
	FILE *file1 = fopen(filename3, "r");
    if (!file1) {
        printf("Error opening file.here\n");
        exit(1);
    }
   frameonedigit = 0;
    while (!feof(file1)) {
        fscanf(file1, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &framestore[frameonedigit][0], &framestore[frameonedigit][1], &framestore[frameonedigit][2], &framestore[frameonedigit][3],
               &framestore[frameonedigit][4], &framestore[frameonedigit][5], &framestore[frameonedigit][6], &framestore[frameonedigit][7],
               &framestore[frameonedigit][8], &framestore[frameonedigit][9], &framestore[frameonedigit][10], &framestore[frameonedigit][11]);
        frameonedigit++;
    }
    fclose(file1);

}

// Function to calculate the Tokhura distance between two vectors
double tokhura(double vec1[DIMENSION], double vec2[DIMENSION]) {
    double weights[DIMENSION] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0}; 
    double dist = 0;
    for (int i = 0; i < DIMENSION; i++) {
        dist += weights[i] * pow(vec1[i] - vec2[i], 2);
    }
    return dist;
}

void writeobservation(const char* filename2){
	  FILE *file2 = fopen(filename2, "w");
    if (file2 == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

        for (int c = 0; c <frameonedigit ; c++) {
            fprintf(file2, "%d ", ansstore[c]);
            ansstore[c] = 0;  // Reset for next use
			if(c != frameonedigit-1)
			   fprintf(file2, "\n");
        }

    fclose(file2);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////step 4
#define N 5  // Number of states
#define O 32 // Number of observation symbols
#define T_MAX 1000 // Maximum number of observations

// Declare matrices with higher precision
long double A[N+1][N+1];  // State transition matrix
long double B[N+1][O+1];  // Observation probability matrix
long double A_comp[N+1][N+1];
long double B_comp[N+1][O+1];
long double pi[N+1]; // Initial state probabilities
long double pi_comp[N+1];
long int obs[1000]={0};
// Declare static arrays for alpha, beta, gamma, and xi
long double alpha[T_MAX][N+1];  // Forward probabilities
long double beta[T_MAX][N+1];   // Backward probabilities
long double gamma[T_MAX][N+1];// State probabilities
long double delta[T_MAX][N+1];
long double zai[T_MAX][N+1][N+1];  // Joint probabilities
int shai[T_MAX][N+1];
int qstar[T_MAX];

int T;

void readfileforstep4(const char* filename4){

	FILE *file1 = fopen(filename4, "r");
    if (!file1) {
        printf("Error opening file.here\n");
        exit(1);
    }
	T=0;
	 while (!feof(file1)) {
	     T++;
        fscanf(file1, "%ld",&obs[T]);
    }
	T--;
	//T--;
    fclose(file1);
}


long double forwardmethod()
{
	//Step 1: Initialisation
	for(int i=1;i<=N;i++)
		alpha[1][i]=pi[i]*B[i][obs[1]];
	//Step 2: Induction
	for(int t=1;t<=T-1;t++)
	{
		for(int j=1;j<=N;j++)
		{
			long double temp=0;
			for(int i=1;i<=N;i++)
				temp+=alpha[t][i]*A[i][j];
			alpha[t+1][j]=temp*B[j][obs[t+1]];
		}
	}
	//Step 3: Termination
	long double probability=0;
	for(int i=1;i<=N;i++)
		probability+=alpha[T][i];
	return probability;
}



void backwardmethod()
{
	//Initialisation
	for(int i=1;i<=N;i++)
		beta[T][i]=1;
	//Induction
	for(int t=T-1;t>=1;t--)
	{
		for(int i=1;i<=N;i++)
		{
			long double temp=0;
			for(int j=1;j<=N;j++)
			{
				temp+=A[i][j]*B[j][obs[t+1]]*beta[t+1][j];
			}
			beta[t][i]=temp;
		}
	}
}

void HMM2()
{
	for(int i=1;i<=N;i++)
	{
		for(int t=1;t<=T;t++)
		{
			long double temp=0;
			for(int j=1;j<=N;j++)
			{
				temp+=alpha[t][j]*beta[t][j];
			}
			gamma[t][i]=(alpha[t][i]*beta[t][i])/temp;
		}
	}
}
long double viterbi_algo()
{
	//Step 1:Initialization
	for(int i=1;i<=N;i++)
	{
		delta[1][i]=pi[i]*B[i][obs[1]];
		shai[1][i]=0;
	}
	//Step 2: Recursion
	long double max,arg_max,temp;
	int state=-1;
	for(int t=2;t<=T;t++)
	{
		for(int j=1;j<=N;j++)
		{
			max=0;
			arg_max=0;
			for(int i=1;i<=N;i++)
			{
				temp=delta[t-1][i]*A[i][j];
				if(temp>arg_max)
				{
					arg_max=temp;
					state=i;
				}
				temp*=B[j][obs[t]];
				if(temp>max)
					max=temp;
			}
			delta[t][j]=max;
			shai[t][j]=state;
		}
	}
	//Step 3:Termination
	long double Pstar;
	max=0,arg_max=0;
	for(int i=1;i<=N;i++)
	{
		if(delta[T][i]>max)
			max=delta[T][i];
		if(delta[T][i]>arg_max)
		{
			arg_max=delta[T][i];
			state=i;
		}
	}
	Pstar=max;
	qstar[T]=state;
	return Pstar;
}


void re_estimation()
{
	long double numerator,denominator;
	//new pi array
	for(int i=1;i<=N;i++)
	{
		pi[i]=gamma[1][i];
	}
	//new A matrix
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			numerator=0,denominator=0;
			for(int t=1;t<=T-1;t++)
			{
				numerator+=zai[t][i][j];
				denominator+=gamma[t][i];
			}
			A[i][j]=numerator/denominator;
		}
	}
	//new B matrix
	for(int j=1;j<=N;j++)
	{
		denominator=0;
		for(int t=1;t<=T;t++)
		{
			denominator+=gamma[t][j];
		}
		for(int k=1;k<=O;k++)
		{
			numerator=0;
			for(int t=1;t<=T;t++)
			{
				if(obs[t]==k)
					numerator+=gamma[t][j];
			}
			B[j][k]=numerator/denominator;
		}
	}
}

void maintain_stochastic()
{
	long double row_sum=0,max=0;
	int index;
	//pi values
	for(int i=1;i<=N;i++)
	{
		if(pi[i]>max)
		{
			max=pi[i];
			index=i;
		}
		row_sum+=pi[i];
	}
	if(row_sum<1)
		pi[index]+=1-row_sum;
	else if(row_sum>1)
		pi[index]-=row_sum-1;
	//For A matrix
	for(int i=1;i<=N;i++)
	{
		row_sum=0,max=0;
		for(int j=1;j<=N;j++)
		{
			if(A[i][j]>max)
			{
				max=A[i][j];
				index=j;
			}
			row_sum+=A[i][j];		
		}
		if(row_sum<1)
			A[i][index]+=1-row_sum;
		else if(row_sum>1)
			A[i][index]-=row_sum-1;
	}
	//For B matrix
	for(int i=1;i<=N;i++)
	{
		row_sum=0,max=0;
		for(int j=1;j<=O;j++)
		{
			if(B[i][j]<(1e-30))
				B[i][j]=1e-30;
			if(B[i][j]>max)
			{
				max=B[i][j];
				index=j;
			}
			row_sum+=B[i][j];		
		}
		if(row_sum<1)
			B[i][index]+=1-row_sum;
		else if(row_sum>1)
			B[i][index]-=row_sum-1;
	}
}

void HMM3()
{
	long double numerator,denominator;
	for(int t=1;t<=T-1;t++)
	{
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
			{
				numerator=alpha[t][i]*A[i][j]*B[j][obs[t+1]]*beta[t+1][j];
				denominator=0;
				for(int x=1;x<=N;x++)
				{
					for(int y=1;y<=N;y++)
					{
						denominator+=alpha[t][x]*A[x][y]*B[y][obs[t+1]]*beta[t+1][y];
					}
				}
				zai[t][i][j]=numerator/denominator;
			}
		}
	}
	re_estimation(); 
	maintain_stochastic();
}

void modeltraining()
{
	long double old_Pstar=-1,new_Pstar=0;
	int iter=0;
	do
	{
		iter++;
		old_Pstar=new_Pstar;
		forwardmethod();	
		backwardmethod();
		HMM2();
		new_Pstar=viterbi_algo();
		HMM3();
		cout<<"iteration"<< iter<<endl;
	}while(iter!=100 && old_Pstar<new_Pstar);
}


void store_in_file(string file)
{
	ofstream store(file);

	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=N;val++)
		{
				store << A[var][val] << " ";
		}
		store<<endl;
	}
	store<<endl<<endl<<endl;
	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=O;val++)
		{
			store << B[var][val] << " ";
		}
		store<<endl;
	}
}


void adding_values()
{
	for(int i=1;i<=N;i++)
		pi_comp[i]+=gamma[1][i];
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			A_comp[i][j]+=A[i][j];
		}
	}
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=O;j++)
		{
			B_comp[i][j]+=B[i][j];
		}
	}
}


void avg_values()
{
	for(int i=1;i<=N;i++){
		pi_comp[i]/=30.0;
	    pi[i]=pi_comp[i];
	}
	for(int var=1;var<=N;var++)
    {
		for(int val=1;val<=N;val++)
		{
				A_comp[var][val] /= 30.0;
				A[var][val]=A_comp[var][val];
		}
	}

	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=O;val++)
		{
				B_comp[var][val] /= 30.0;
				B[var][val]=B_comp[var][val];
		}
	}
	maintain_stochastic();
}



/////////////////////////////////////////////////////////////////testing

long int observation[100];
long double Anew[N+1][N+1];  // State transition matrix
long double Bnew[N+1][O+1];
long double alphanew[T_MAX][N+1];
void taketestingobservation(const char* filename4){

	FILE *file1 = fopen(filename4, "r");
    if (!file1) {
        printf("Error opening file.here\n");
        exit(1);
    }
	T=0;
	 while (!feof(file1)) {
		 T++;
        fscanf(file1, "%ld",&observation[T]);
    }
	// T--;
    fclose(file1);
}

long double probability = 0;
long double prob[10]={0};
int correct,wrong;

void testing(int val)
{
	int idx=-1;
	long double maxi=LDBL_MIN;
	for(int i=0;i<=9;i++)
	{
		if(maxi<=prob[i])
		{
			idx=i;
			maxi=prob[i];
		}
	}
	
	if(val==idx)
	{
	    correct++;
		cout<<val<<"-->"<<"correct"<<endl;
	}
	else{
		wrong++;
		cout<<val <<"-->"<<"wrong"<<endl;
	}
}

void answer(int val)
{
	for(long long int k=0;k<=9;k++)
	{
		char extrafilepath5[100]= "";
		sprintf(extrafilepath5, "./OUTPUT_a_b_final/244101016_%d.txt",k);
		ifstream fp(extrafilepath5);
		int i=1,j=1;
		while(!fp.eof())
		{
			fp >> A[i][j];
			j++;
			if(j==6)
			{
				i++;
				j=1;
			}
			if(i==6){
			  break;
			}
		}
		string get;
		i=1,j=1;
		while(!fp.eof())
		{
			fp >> B[i][j];
			j++;
			if(j==33)
			{
				i++;
				j=1;
			}
			if(i==6)
			{
				break;
			}
		}
		probability = 0;
    for(int i=1;i<=N;i++){
        alpha[1][i]=pi[i]*B[i][observation[1]];
    }

    for(int t=1;t<=T-1;t++)
    {
        for(int j=1;j<=N;j++)
        {
            long double sum=0;
            for(int i=1;i<=N;i++)
            {
                sum = sum + alpha[t][i]*A[i][j];
            }
            alpha[t+1][j] = sum * B[j][observation[t+1]];
        }
    }
    
    for(int i=1;i<=N;i++)
    {
       probability = probability + alpha[T][i];
    }
	    prob[k] = probability;
		fp.close();
	}
	//for(int i =0;i<10;i++){
	//	cout<<prob[i]<<endl;
	//}
	probability=0;
	testing(val);
	for(int i=0;i<10;i++){
		prob[i]=0;
	}
}



///////////////////////////////////////////////recording phase
#pragma comment(lib, "winmm.lib")
#define DATA 16025
#define BATCH 300

short int waveIn[32050];
int main_arr[DATA];

void recordvaluetofile(const char* filename2){
	  FILE *file2 = fopen(filename2, "w");
    if (file2 == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

        for (int c = 5000; c <27000 ; c++) {
            fprintf(file2, "%d ", main_arr[c]);
           // ansstore[c] = 0;  // Reset for next use
			if(c != 26999)
			   fprintf(file2, "\n");
        }

    fclose(file2);
}

void PlayRecord()
{
    const int NUMPTS = 16025 * 2;   // 2 seconds
    int sampleRate = 16025;
    HWAVEOUT hWaveOut;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;     // simple, uncompressed format
    pFormat.nChannels = 1;                    //  1=mono, 2=stereo
    pFormat.nSamplesPerSec = sampleRate;      // 44100
    pFormat.nAvgBytesPerSec = sampleRate * 2; // = nSamplesPerSec * n.Channels * wBitsPerSample/8
    pFormat.nBlockAlign = 2;                  // = n.Channels * wBitsPerSample/8
    pFormat.wBitsPerSample = 16;              //  16 for high quality, 8 for telephone-grade
    pFormat.cbSize = 0;

    // Specify playback parameters
    waveOutOpen(&hWaveOut, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
    WAVEHDR WaveInHdr;

    // Set up and prepare header for output
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS * 2;
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
    waveOutPrepareHeader(hWaveOut, &WaveInHdr, sizeof(WAVEHDR));

    cout << "playing the sound u have spoken..." << endl;
    waveOutWrite(hWaveOut, &WaveInHdr, sizeof(WaveInHdr)); // Playing the data
    Sleep(1000*2); //Sleep for as long as there was recorded

    waveOutClose(hWaveOut);
}
void startrecord(){
	 const int NUMPTS = 16025 *2;   // 2 seconds
    int sampleRate = 16025;  
    HWAVEIN hWaveIn;
    MMRESULT result;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;     // simple, uncompressed format
    pFormat.nChannels = 1;                    //  1=mono, 2=stereo
    pFormat.nSamplesPerSec = sampleRate;      // 8.0 kHz, 11.025 kHz, 22.05 kHz, and 44.1 kHz
    pFormat.nAvgBytesPerSec = sampleRate * 2; // =  nSamplesPerSec × nBlockAlign
    pFormat.nBlockAlign = 2;                  // = (nChannels × wBitsPerSample) / 8
    pFormat.wBitsPerSample = 16;              //  16 for high quality, 8 for telephone-grade
    pFormat.cbSize = 0;

    // Specify recording parameters
    result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
    WAVEHDR WaveInHdr;

    // Set up and prepare header for input
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS * 2;
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    // Insert a wave input buffer
    result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    // Commence sampling input
    result = waveInStart(hWaveIn);
    cout << "recording for 2 seconds..." << endl;
    Sleep(2 * 1000);

    // Wait until finished recording
    waveInClose(hWaveIn);

    // Convert the recorded data into an integer array
    for (int i = 0; i < NUMPTS; i++) {
        main_arr[i] = (int)waveIn[i];
    }

    PlayRecord();
	char outputFilePathforstep6[100]= "";
	
	int A =1;
	sprintf(outputFilePathforstep6, "./RECORDING/244101016_%d.txt",A);
     recordvaluetofile(outputFilePathforstep6);
	 char inputFilePath6[100] = "";
	 sprintf(inputFilePath6, "./RECORDING/244101016_%d.txt",A); 

            // Read the signal data from the input file
     ReadSignalDataFromFile(inputFilePath6);
	 CorrectDCOffsetInSignal();
     NormalizeSignalValues();
	 TOTAL_FRAMES = 22000/320;
	  SelectSteadyStateFramesFromSignal();
            for (int j = 0; j < TOTAL_FRAMES; j++) {
                ApplyHammingWindowToFrame(j);
            }
			 for (int j = 0; j < TOTAL_FRAMES; j++) {
                ComputeAutoCorrelationForFrame(j);
                ComputeLPCoefficientsForFrame(j);
                ComputeCepstralCoefficientsForFrame(j);
            }
			 char outputFilePath[100] = "";
			 sprintf(outputFilePath, "./OUTPUTRECORDING/244101016_%d.txt",A);
              WriteWeightedCepstralCoeffsToFile(outputFilePath);
			  char inputfileforfindingsequence[100]="";
			  sprintf(inputfileforfindingsequence, "./OUTPUTRECORDING/244101016_%d.txt",A); 
	          readfileforstep3(inputfileforfindingsequence);
			  for(int t =0;t<frameonedigit;t++){
		 double min_dist = DBL_MAX;
		 int nearest =0;
		for(int u=0;u<32;u++){
			double tok_dist = tokhura(framestore[t], codebook[u]);
			if (tok_dist < min_dist) {
                min_dist = tok_dist;
                nearest = u;
            }
		}
		ansstore[t]=nearest+1;
		nearest =0;
	}
char outputFilePathforstep8[100] = "";
    sprintf(outputFilePathforstep8, "./RECORDINGSEQUENCE/244101016_%d.txt",A);
     writeobservation(outputFilePathforstep8);
	 char inputFilePathforstep5[100] = "";
		   sprintf(inputFilePathforstep5, "./RECORDINGSEQUENCE/244101016_%d.txt",A);
	       taketestingobservation(inputFilePathforstep5);
		   answer(1);

}


//main function started
int main(int argc, char* argv[])
{
	
	char inputFilePath[100] = "";
    char outputFilePath[100] = "";
    for (int v = 0; v < 10; v++) {
        for (int i = 1; i < 41; i++) {
            totalSamples = 0;
			TOTAL_FRAMES=0;

            // Prepare file names for input and output
            sprintf(inputFilePath, "./244101016_dataset/English/txt/244101016_E_%d_%d.txt",v, i); 

            // Read the signal data from the input file
            ReadSignalDataFromFile(inputFilePath);

            // Process the signal: correct DC offset, normalize, select frames, and apply Hamming window
            CorrectDCOffsetInSignal();
            NormalizeSignalValues();


			TOTAL_FRAMES = totalSamples/FRAME_SIZE;


            SelectSteadyStateFramesFromSignal();
            for (int j = 0; j < TOTAL_FRAMES; j++) {
                ApplyHammingWindowToFrame(j);
            }

            // Compute autocorrelation, LPC coefficients, and Cepstral coefficients for each frame
            for (int j = 0; j < TOTAL_FRAMES; j++) {
                ComputeAutoCorrelationForFrame(j);
                ComputeLPCoefficientsForFrame(j);
                ComputeCepstralCoefficientsForFrame(j);
            }
			int q = 40*v+i;
			  sprintf(outputFilePath, "./OUTPUT_CEPSTRAL/244101016_%d.txt",q);
              WriteWeightedCepstralCoeffsToFile(outputFilePath);
        }
        // Write the computed weighted cepstral coefficients to the output file
       // sprintf(outputFilePath, "./OUTPUT_CEPSTRAL/244101016_%c.txt", v);
        //WriteWeightedCepstralCoeffsToFile(outputFilePath);
    }

	 mergeFiles("UNIVERSE.csv");

	printf("Step 1 done");




	//step 2 start
	int M = MAX_VECTORS;  // Universe size
    
      // Codebook with maximum size K
    int region[MAX_VECTORS];  // To store the region (nearest centroid) for each universe vector

    // Load the universe from a CSV file
    store_all_value_of_universe("./UNIVERSE.csv", universe);

    // Step 1: Initialize the codebook with one centroid (centroid of the universe)
    initialize_codebook(codebook, universe, M);
    int current_size = 1;

    // Repeat the process until the codebook size reaches K
    while (current_size < K) {
        // Step 2: Double the codebook size by splitting each centroid
        split_codebook(codebook, &current_size);

        double distortion = 0, dist, avg_dist;
        int m = 0;

        do {
            if (m >= 1)
                distortion = avg_dist;

            // Step 3: Assign vectors to regions using Euclidean distance
            form_cluster(universe, M, codebook, current_size, region, &dist);

            // Step 4: Compute average distortion using Tokhura distance
            avg_dist = dist / M;
           // printf("Iteration %d, Average Distortion: %lf\n", m, avg_dist);

            // Step 5: Update centroids
            update_centroids(universe, M, codebook, current_size, region);

            m++;
        } while (fabs(avg_dist - distortion) > DELTA);
    }
//	for (int i = 0; i < K; i++) {
  //          printf("Centroid %d: ", i);
    //        for (int j = 0; j < DIMENSION; j++) {
      //          printf("%lf ", codebook[i][j]);
        //    }
          //  printf("\n");
     // }

    FILE *file = fopen("codebook.csv", "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Write cepstral coefficients for each frame to the file
    for (int q = 0; q < 32; q++) {
        for (int w = 0; w<12; w++) {
           // cepstralCoefficients[frame][coef] /= 20;  // Normalize coefficients
            fprintf(file, "%lf ", codebook[q][w]);
          //  cepstralCoefficients[frame][coef] = 0;  // Reset for next use
        }
        fprintf(file, "\n");
    }

    fclose(file);

	printf("step 2 done");


	//////////////////////////////////////////////// step3

	int nearest=0;
	char inputFilePathforstep3[100] = "";
	char outputFilePathforstep3[100]= "";

	for(int q=1;q<401;q++){
		printf("hello");
	 sprintf(inputFilePathforstep3, "./OUTPUT_CEPSTRAL/244101016_%d.txt",q); 
	readfileforstep3(inputFilePathforstep3);
	
	//int ansstore[100];
	for(int t =0;t<frameonedigit;t++){
		 double min_dist = DBL_MAX;
		for(int u=0;u<32;u++){
			double tok_dist = tokhura(framestore[t], codebook[u]);
			if (tok_dist < min_dist) {
                min_dist = tok_dist;
                nearest = u;
            }
		}
		ansstore[t]=nearest+1;
		nearest =0;
	}
    sprintf(outputFilePathforstep3, "./OUTPUT_OBSERVATION/244101016_%d.txt",q);
     writeobservation(outputFilePathforstep3);
	 for(int p=0;p<100;p++){
	    ansstore[p]=0;
	 }
	 frameonedigit=0;
	}

   







	///////////////////////////////////////////////////////////// training part step 4
	char inputFilePathforstep4[100] = "";
	char outputFilePathforstep4[100]= "";
	char extrafilepath4[100]= "";
	for(int z=0;z<10;z++){
		for(int y=1;y<=30;y++){
			int num = z*40+y;
        sprintf(inputFilePathforstep4, "./OUTPUT_OBSERVATION/244101016_%d.txt",num); 
	    readfileforstep4(inputFilePathforstep4);

		//A matrix initialize
		for(int i=1;i<=N;i++)
	    {
		       for(int j=1;j<=N;j++)
		       {
			       if(j==i && j+1<=N)
			       {
				        A[i][j]=0.8;
				        A[i][j+1]=0.2;
			       }
			       else if(j==i && j+1>N)
			        	A[i][j]=1;
		       }
	      }

		//b matrix initialize
		for (int i = 1; i <= N; i++)
          for (int j = 1; j <= O; j++)
             B[i][j]=1.0/32;

		for(int i=1;i<=N;i++)
	    {
		  if(i==1)
			 pi[i]=1.0;
		  else
			 pi[i]=0.0;
	    }
		modeltraining();
		adding_values();
		sprintf(outputFilePathforstep4, "./OUTPUT_a_b/244101016_%d.txt",num);
		store_in_file(outputFilePathforstep4);
	  }
		avg_values();
		sprintf(extrafilepath4, "./OUTPUT_a_b_final/244101016_%d.txt",z);
		store_in_file(extrafilepath4);
	}


	/////////////////////////////////////////////////testing step 5
	
	for(long long int i=0;i<=9;i++)
	{
	   for(long long int k=31;k<=40;k++)
	   {
		   int num = i*40+k;
		   char inputFilePathforstep5[100] = "";
		   sprintf(inputFilePathforstep5, "./OUTPUT_OBSERVATION/244101016_%d.txt",num);
	       taketestingobservation(inputFilePathforstep5);
		   answer(i); 
	   }
	   cout<<endl;
	}
	cout<<"accuracy -->"<<correct<<"%";

	/////////////////////////////////////////////////recording testing phase
	Sleep(3 * 1000);

	startrecord();
	return 0;
}
