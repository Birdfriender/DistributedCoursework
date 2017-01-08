#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"

const int SEND_TOP_ROW = 1;
const int SEND_BOTTOM_ROW = 2;
const int SEND_ALL_ROWS = 3;

const double precision = 0.001;

void calcStartEndRow(int thisRank, int numActiveRows, int numProc, 
	int* startRow, int* endRow)
{
	int allocPerCore = numActiveRows/numProc;
	if (thisRank < numActiveRows % numProc)
	{	
		// +1 because the top row isnt used
		*startRow = (thisRank * (allocPerCore + 1)) + 1;
		*endRow = *startRow + allocPerCore;
	}
	else
	{
		*startRow = ((numActiveRows % numProc) * (allocPerCore + 1)) + 
		((thisRank - (numActiveRows % numProc)) * allocPerCore) + 1; 
		*endRow = *startRow + allocPerCore - 1;
	}
}

//Arguments
//Size - the width/height of the array
int main(int argc, char* argv[])
{
	int thisRank, numProc;
	int numRows;
	double **array, **resultArray;
	double *mem, *resultMem;
	

	//Do all the mpi stuff
	int sc = MPI_Init(&argc, &argv);
	if(sc != MPI_SUCCESS)
	{
		printf("Failed to init MPI");
		MPI_Abort(MPI_COMM_WORLD, sc);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);

	//Parse the arguments
	numRows = strtoul(argv[1], NULL, 10);

	//Allocate the array
	mem = malloc(sizeof(double) * numRows * numRows);
	if (mem == NULL)
	{
		printf("Malloc Failed");
		return 1;
	}

	resultMem = malloc(sizeof(double) * numRows * numRows);
	if (resultMem == NULL)
	{
		printf("Malloc Failed");
		return 1;
	}

	array = malloc(sizeof(double*) * numRows);
	if(array == NULL)
	{
		printf("Malloc Failed");
		return 1;
	}

	resultArray = malloc(sizeof(double*) * numRows);
	if(resultArray == NULL)
	{
		printf("Malloc Failed");
		return 1;
	}

	for(int i = 0; i < numRows; ++i)
	{
		array[i] = &mem[i * numRows];
		resultArray[i] = &resultMem[i * numRows];
	}
	

	//Fill the array with 0s
	for(int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numRows; j++)
		{
			array[i][j] = 0;
			resultArray[i][j] = 0;
		}
	}

	//Fill the edges with 1s
	//Do this for both arrays so we can swap them later
	for(int i = 0; i < numRows; i++)
	{
		array[i][0] = 1;
		array[i][numRows-1] = 1;
		array[0][i] = 1;
		array[numRows-1][0] = 1;

		resultArray[i][0] = 1;
		resultArray[i][numRows-1] = 1;
		resultArray[0][i] = 1;
		resultArray[numRows-1][0] = 1;
	}

	int startRow = 0, endRow = 0;
	//How many rows are we actually operating on
	int numActiveRows = numRows - 2;
	calcStartEndRow(thisRank, numActiveRows, numProc, &startRow, &endRow);

	time_t startTime;
	startTime = time(NULL);


	double highestChange = 0.0;

	//Loop variables
	//declare them here so we dont have to allocate memory every iteration
	int i = 0, j = 0;
	double diff = 0.0;
	int stopFlag = 0, tempFlag = 0;
	double **temp;
	//Settles the cores rows
	//No functionization we die like men
	MPI_Request req;
	do
	{
		highestChange = 0.0;
		for (i = startRow; i <= endRow; i++)
		{
			for (j = 1; j <= numActiveRows; j++)
			{
				resultArray[i][j] = 
					(array[i - 1][j] +
					array[i + 1][j] +
					array[i][j - 1] +
					array[i][j + 1])/4.0;
				diff = fabs(resultArray[i][j] - array[i][j]);
				if(diff > highestChange)
				{
					highestChange = diff;
				}
			}
		}
		
		//The first processor doesnt need to send its top row 
		//anywhere since its the top row of the whole array
		if(thisRank != 0)
		{
			MPI_Isend(resultArray[startRow], numRows, MPI_DOUBLE, thisRank - 1, 
				SEND_TOP_ROW, MPI_COMM_WORLD, &req);
		}

		//likewise the last processor and the bottom row
		if(thisRank != numProc - 1)
		{
			MPI_Recv(resultArray[endRow + 1], numRows, MPI_DOUBLE, thisRank + 1, 
				SEND_TOP_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Isend(resultArray[endRow], numRows, MPI_DOUBLE, thisRank + 1, 
				SEND_BOTTOM_ROW, MPI_COMM_WORLD, &req);
		}

		if(thisRank != 0)
		{
			MPI_Recv(resultArray[startRow - 1], numRows, MPI_DOUBLE, 
				thisRank - 1, SEND_BOTTOM_ROW, MPI_COMM_WORLD, 
				MPI_STATUS_IGNORE);
		}

		//Use a flag to check whether we should finish
		tempFlag = highestChange <= precision;

		//If any processor isnt ready to stop, we all keep going
		MPI_Allreduce(&tempFlag, &stopFlag, 1, MPI_INT, 
			MPI_MIN, MPI_COMM_WORLD);

		temp = resultArray;
		resultArray = array;
		array = temp;
	}while(stopFlag == 0);

	time_t endTime;
	endTime = time(NULL);
	struct tm * timeinfo;
	timeinfo = localtime(&startTime);
	double timeTaken = difftime(endTime, startTime);

	if(thisRank == 0)
	{
		int procStart = 0, procEnd = 0;
		for (int i =1; i < numProc; ++i)
		{
			calcStartEndRow(i, numActiveRows, numProc, &procStart, &procEnd);
			for (int j = procStart; j <= procEnd; ++j)
			{
				MPI_Recv(array[j], numRows, MPI_DOUBLE, 
					i, SEND_ALL_ROWS, MPI_COMM_WORLD, 
					MPI_STATUS_IGNORE);
			}	
		}
		
	}

	
	if(thisRank == 0)
	{
		char string[100] = "";
		sprintf(string, "DistributedResults-%d-%d", numProc, numRows);


		FILE *f = fopen(string, "a+");
		fprintf(f, "Number of processors: %d\n", numProc);
		fprintf(f, "Size of array: %d\n", numRows);
		fprintf(f, "Time Taken: %lf\n", timeTaken);
		for (int j = 0; j < numRows; ++j)
		{
			for (int k = 0; k < numRows; ++k)
			{
				fprintf(f, "%lf ", array[j][k]);
			}
			fprintf(f, "\n");
		}

		fclose(f);
	}


	free(array); free(resultArray); free(mem); free(resultMem);

	MPI_Finalize();
    return 0;
}

