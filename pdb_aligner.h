#ifndef PDBALIGNER
#define PDBALIGNER

#include <cmath>
#include <iostream>

using namespace std;

class Aligner
{
public:
	double solution[420][3] = {{0}};

	void transform( double translation[3], double rotation[3], double atoms[][3], int qtdAtoms );
	double calcRMSD( double reference[][3], int qtdAtoms );	
};

void Aligner::transform( double translation[3], double rotation[3], double atoms[][3], int qtdAtoms )
{
	//cout << "Initializing rotation..." << qtdAtoms << endl;
	double c[420][3] = {{0}};
	double d[420][3] = {{0}};

	double rotX[3][3] = { { 1.0, 0.0, 0.0 }, { 0.0, cos( rotation[0] ), -sin( rotation[0] ) }, { 0.0, sin( rotation[0] ), cos( rotation[0] ) } };
	double rotY[3][3] = { { cos( rotation[1] ), 0.0, sin( rotation[1] )}, { 0.0, 1.0, 0.0 }, { -sin( rotation[1] ), 0.0, cos( rotation[1] ) } };
    double rotZ[3][3] = { { cos( rotation[2] ), -sin( rotation[2] ), 0.0 }, { sin( rotation[2] ), cos( rotation[2] ), 0.0 }, { 0.0, 0.0, 1.0 } };

	for ( int i = 0; i < qtdAtoms; i++ )
	{
		// rotacao em X
		for ( int j = 0; j < 3; j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				c[i][j] += atoms[i][k] * rotX[k][j];
			}
		}

		// rotacao em Y
		for ( int j = 0; j < 3; j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				d[i][j] += c[i][k] * rotY[k][j];
			}
		}

		// rotacao em Z
		for ( int j = 0; j < 3; j++ )
		{
			for ( int k = 0; k < 3; k++ )
			{
				solution[i][j] += d[i][k] * rotZ[k][j];
			}
		}
	}

	//cout << "Initializing translation..." << endl;
	for ( int i = 0; i < qtdAtoms; i++ )
	{
		for ( int j = 0; j < 3; j++ )
		{
			solution[i][j] += translation[j];
			//cout << solution[i][j] << " ";
		}

		//cout << endl;
	}

	//cout << "finished" << endl;
}

double Aligner::calcRMSD( double reference[][3], int qtdAtoms )
{
	double sum_distance = 0;

	for ( int i = 0; i < qtdAtoms; i++ )
	{
		sum_distance += pow( reference[i][0] - solution[i][0], 2 );
		sum_distance += pow( reference[i][1] - solution[i][1], 2 );
		sum_distance += pow( reference[i][2] - solution[i][2], 2 );
	}

	return sqrt( sum_distance/2.0 );
}

#endif