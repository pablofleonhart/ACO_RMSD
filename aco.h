#ifndef ACO
#define ACO

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include "pdb_aligner.h"

using namespace std;

class Aco;
class Region;

double minT = -100.0, maxT = 100.0, minR = 0, maxR = 2.0*M_PI;
int ants = 150, G = 50, L = 100;
int qtdRegions = 200;

class Region
{
public:
	double component[6];
	double fitness;
};

class Aco
{
public:
	double refPosAtoms[420][3], modPosAtoms[420][3];
	int qtdAtoms;
	vector<Region> regions;

	void generateInitialRegions();
	void calcFitness();
};

void Aco::generateInitialRegions()
{
	for ( int i = 0; i < qtdRegions; i++ )
	{
		Region r;
		double val = (double)rand() / RAND_MAX;
		r.component[0] = minT + ( val * ( maxT - minT ) );
		val = (double)rand() / RAND_MAX;
		r.component[1] = minT + ( val * ( maxT - minT ) );
		val = (double)rand() / RAND_MAX;
		r.component[2] = minT + ( val * ( maxT - minT ) );
		val = (double)rand() / RAND_MAX;
		r.component[3] = minR + ( val * ( maxR - minR ) );
		val = (double)rand() / RAND_MAX;
		r.component[4] = minR + ( val * ( maxR - minR ) );
		val = (double)rand() / RAND_MAX;
		r.component[5] = minR + ( val * ( maxR - minR ) );
	    
	    regions.push_back( r );
	}
}

void Aco::calcFitness()
{
	double best = 10000000;
	Region bRegion;
	vector<Region>::iterator it;

	for ( it = regions.begin(); it != regions.end(); ++it )
	{
		double translation[3] = { it->component[0], it->component[1], it->component[2] };
    	double rotation[3] = { it->component[3], it->component[4], it->component[5] };

		Aligner aligner;

		aligner.transform( translation, rotation, modPosAtoms, qtdAtoms );

		double fitness = aligner.calcRMSD( refPosAtoms, qtdAtoms );

		it->fitness = fitness;

		if ( fitness < best )
		{
			best = fitness;
			bRegion = *it;
		}
	}

	cout << "Best initial fitness: " << best << endl;

	for ( int i =0; i < 6; i++ )
	{
		cout << bRegion.component[i] << " ";
	}

	cout << endl;
}

#endif