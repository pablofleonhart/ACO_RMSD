#include <iostream>
#include "pdb_reader.h"
#include "aco.h"

using namespace std;

int main( int argc, char *argv[] )
{
    PdbReader refPdb, modPdb;

    refPdb.readFile( "reference.pdb" );
    modPdb.readFile( "1ACW-01.pdb" );

    modPdb.adjustAtoms( refPdb.atomsPairs );
    //cout << refPdb.atoms.size() << " " << modPdb.atoms.size() << endl;
    
    refPdb.calcBackbonePos();
    refPdb.calcCaPos();
    modPdb.calcBackbonePos();
    modPdb.calcCaPos();

    //cout << refPdb.aSize << " " << modPdb.aSize << endl;
    //cout << refPdb.bSize << " " << modPdb.bSize << endl;

    /*for ( int i = 0; i < refPdb.atoms.size(); i++ )
    {
        cout << refPdb.posAtoms[i][0] << " " << refPdb.posAtoms[i][1] << " " << refPdb.posAtoms[i][2] << endl;
    }

    for ( int i = 0; i < modPdb.atoms.size(); i++ )
    {
        cout << modPdb.posAtoms[i][0] << " " << modPdb.posAtoms[i][1] << " " << modPdb.posAtoms[i][2] << endl;
    }*/

    double solution[6] = { 0, 4, 0, 0.5, 1.3, 6.1 };  // vetor de solucao obtido do ACO

    double translation[3] = { solution[0], solution[1], solution[2] };
    double rotation[3] = { solution[3], solution[4], solution[5] };

    Aco aco;

    aco.qtdAtoms = refPdb.atoms.size();

    for ( int i = 0; i < refPdb.atoms.size(); i++ )
    {
    	aco.refPosAtoms[i][0] = refPdb.posAtoms[i][0];
    	aco.refPosAtoms[i][1] = refPdb.posAtoms[i][1];
    	aco.refPosAtoms[i][2] = refPdb.posAtoms[i][2];

    	aco.modPosAtoms[i][0] = modPdb.posAtoms[i][0];
    	aco.modPosAtoms[i][1] = modPdb.posAtoms[i][1];
    	aco.modPosAtoms[i][2] = modPdb.posAtoms[i][2];
    }

    aco.generateInitialRegions();
    aco.calcFitness();

	return 0;
}