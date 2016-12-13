#include "mpi.h"

//Arguments
//Size - the width/height of the array
int main(int argc, char* argv[])
{
	int thisRank, numProc;
	double numRows;
	double **array, **resultArray;

	int precision = 0.001;

	//Do all the mpi stuff
	int sc = MPI_Init(&argc, &argv);
	if(sc != MPI_SUCCESS)
	{
		printf("Failed to init MPI");
		MPI_Abort(MPI_COMM_WORLD, sc);
	}

	MPI_Comm_Rank(MPI_COMM_WORLD, &thisRank);
	MPI_Comm_Size(MPI_COMM_WORLD, &numProc);

	//Parse the arguments
	numRows = strtoul(argv[1], NULL, 10);

	//Allocate the array
	array = malloc(sizeof(double*) * numRows);
	if (array == NULL)
	{
		printf("Malloc Failed");
		return 1;
	}

	resultArray = malloc(sizeof(double*) * numRows);
	if (resultArray == NULL)
	{
		printf("Malloc Failed");
		return 1;
	}

	for (int i = 0; i < numRows; i++)
	{
		array[i] = malloc(sizeof(double) * numRows);
		if (array[i] == NULL)
		{
			printf("Malloc Failed");
			return 1;
		}
		resultArray[i] = malloc(sizeof(double) * ARRAY_DIMENSIONS);
		if (resultArray[i] == NULL)
		{
			printf("Malloc Failed");
			return 1;
		}
	}

	//Fill the array with 0s
	for(double i = 0; i < numRows; i++)
	{
		for (double j = 0; j < count; j++)
		{
			array[i][j] = 0;
		}
	}

	//Fill the edges with 1s
	for(double i = 0; i < numRows; i++)
	{
		array[i][0] = 1;
		array[i][numRows-1] = 1;
		array[0][i] = 1;
		array[numRows-1][0] = 1;
	}

	int allocPerCore = numRows/numProc;
	int startRow, endRow;

	if (thisRank < numRows % numProc)
	{
		startRow = thisRank * (allocPerCore + 1);
		endRow = startRow + allocPerCore + 1;
	}
	else if (thisRank = numRows % numProc)
	{
		startRow = thisRank * (allocPerCore + 1);
		endRow = startRow + allocPerCore;
	}


	int highestChange;

	do
	{
		double **temp = resultArray;
		resultArray = array;
		array = temp;
		highestChange = 0;


	}while(highestChange > precision)

	MPI_Finalize();
    return 0;
}