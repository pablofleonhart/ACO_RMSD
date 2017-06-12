import copy

class PDBReader:
	ALPHA_CARBON = "CA"
	BACKBONE_ATOMS = ( "N", "CA", "C", "O" )
	ATOM_TAG = "ATOM"
	END_TAG = "TER"

	fileName = ""
	atoms = []
	aminoAcids = []
	posAtoms = []
	backbone = []
	alpha = []

	def __init__( self, fileName ):
		self.fileName = fileName

		self.readFile()

	def isNumber( self, value ):
		try:
			int( value )
			return True

		except ValueError:
			return False

	def readFile( self ):
		self.atoms = []
		self.aminoAcids = []
		self.posAtoms = []

		finish = False
		file = open( self.fileName, "r" )

		while not finish:
			line = file.readline()

			if not line:
				finish = True

			else:
				line = line.split()

				if line[0] == self.END_TAG:
					finish = True

				elif line[0] == self.ATOM_TAG:
					atom = line[2]

					if self.isNumber( line[4] ):
						aminoAcid = int( line[4] )

					elif self.isNumber( line[5] ):
						aminoAcid = int( line[5] )

					posInit = 0

					for i in xrange( len( line ) ):
						if "." in line[i]:
							posInit = i
							break

					pos = map( float, line[posInit:posInit+3] )
					self.atoms.append( atom )
					self.posAtoms.append( pos )
					self.aminoAcids.append( aminoAcid )

		file.close()

	def adjustAtoms( self, refAtoms, refAminoAcids ):
		auxAtoms = []
		auxPos = []
		auxAminoAcids = []

		for i in xrange( len( refAtoms ) ):  
			if ( refAtoms[i], refAminoAcids[i] ) in zip( self.atoms, self.aminoAcids ):
				index = zip( self.atoms, self.aminoAcids ).index( ( refAtoms[i], refAminoAcids[i] ) )
				auxAtoms.append( self.atoms.pop( index ) )
				auxPos.append( self.posAtoms.pop( index ) )
				auxAminoAcids.append( self.aminoAcids.pop( index ) )

		self.atoms = copy.deepcopy( auxAtoms )
		self.posAtoms = copy.deepcopy( auxPos )
		self.aminoAcids = copy.deepcopy( auxAminoAcids )

	def calcBackbonePos( self ):
		self.backbone = []
		for i in range( len( self.atoms ) ):
			if self.atoms[i] in self.BACKBONE_ATOMS:
				self.backbone.append( self.posAtoms[i] )

	def calcCaPos( self ):
		self.alpha = []
		for i in range( len( self.atoms ) ):
			if self.atoms[i] == self.ALPHA_CARBON:
				self.alpha.append( self.posAtoms[i] )