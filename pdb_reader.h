#ifndef PDBREADER
#define PDBREADER

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

string::size_type sz;

string ALPHA_CARBON = "CA";
vector<string> BACKBONE_ATOMS = { "N", "CA", "C", "O" };
string ATOM_TAG = "ATOM";
string END_TAG = "TER";

class PdbReader
{
public:
    vector<string> atoms, aminoAcids;
    vector<pair<string, string>> atomsPairs;
    double posAtoms[420][3], backbone[420][3], alfa[420][3];
    int bSize, aSize;

    bool isNumber( string value );
    void readFile( string name );
    void adjustAtoms( vector<pair<string, string>> pairs );
    void calcBackbonePos();
    void calcCaPos();
};

bool PdbReader::isNumber( string value )
{
    return !value.empty() && 
    		find_if( value.begin(), value.end(), [](char c) { return ! isdigit(c); } ) == value.end();
}

void PdbReader::readFile( string name )
{
    string line = "", tag, atom, trash, aacid, x, y, z;
    int count = 0, stop = 0;

    ifstream file;
    file.open( name );

    while ( !file.eof() && stop == 0 )
    {
        file >> tag;

        if ( tag == "ATOM" )
        {
            file >> trash >> atom >> trash >> aacid;

            if ( !isNumber( aacid ) )
            {
                file >> aacid >> x >> y >> z >> trash >> trash >> trash;
            }

            else
            {
                file >> x >> y >> z >> trash >> trash;
            }

            atoms.push_back( atom );
            aminoAcids.push_back( aacid );
            atomsPairs.push_back( make_pair( atom, aacid ) );

            posAtoms[count][0] = stod( x );
            posAtoms[count][1] = stod( y );
            posAtoms[count][2] = stod( z );

            count++;
        }
        else if ( tag == "TER" )
        {
            stop = 1;
        }
    }

    file.close();
}

void PdbReader::adjustAtoms( vector<pair<string, string>> pairs )
{
    vector<string> auxAtoms, auxAacid;
    double auxPos[atoms.size()][3];
    vector<pair<string, string>>::iterator it;
    int count = 0;

    for ( int i = 0; i < pairs.size(); i++ )
    {
    	it = find( atomsPairs.begin(), atomsPairs.end(), pairs[i] );

    	if ( it != atomsPairs.end() )
    	{
    		int index = it - atomsPairs.begin();
    		auxAtoms.push_back( atoms[index] );
    		auxAacid.push_back( aminoAcids[index] );
    		auxPos[count][0] = posAtoms[index][0];
    		auxPos[count][1] = posAtoms[index][1];
    		auxPos[count][2] = posAtoms[index][2];
            count++;
    	}
    }

    atoms = auxAtoms;
    aminoAcids = auxAacid;

    for ( int i = 0; i < atoms.size(); i++ )
    {
    	posAtoms[i][0] = auxPos[i][0];
    	posAtoms[i][1] = auxPos[i][1];
    	posAtoms[i][2] = auxPos[i][2];
    }
}

void PdbReader::calcBackbonePos()
{
    int countB = 0;

    for ( int i = 0; i < atoms.size(); i++ )
    {
        if ( find( BACKBONE_ATOMS.begin(), BACKBONE_ATOMS.end(), atoms[i] ) != BACKBONE_ATOMS.end() )
        {
            backbone[countB][0] = posAtoms[i][0];
            backbone[countB][1] = posAtoms[i][1];
            backbone[countB][2] = posAtoms[i][2];
            countB++;
        }
    }

    bSize = countB;
}
     
void PdbReader::calcCaPos()
{
    int countA = 0;

    for ( int i = 0; i < atoms.size(); i++ )
    {
        if ( atoms[i] == ALPHA_CARBON )
        {
            alfa[countA][0] = posAtoms[i][0];
            alfa[countA][1] = posAtoms[i][1];
            alfa[countA][2] = posAtoms[i][2];
            countA++;
        }
    }

    aSize = countA;
}

#endif