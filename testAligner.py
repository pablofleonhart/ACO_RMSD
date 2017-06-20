import math
import numpy as np
import copy
import sys

class tAligner( object ):
    def align( self,transformation, mob_atoms):
        translation = np.matrix([transformation[0:3]]*len(mob_atoms))
        rot = transformation[3:6]

        #print translation
        #print rot
        
        rotX = np.matrix([[1.0, 0.0, 0.0], [0.0, math.cos(rot[0]), -math.sin(rot[0])], [0.0, math.sin(rot[0]), math.cos(rot[0])]])
        rotY = np.matrix([[math.cos(rot[1]), 0.0, math.sin(rot[1])], [0.0, 1.0, 0.0], [-math.sin(rot[1]), 0.0, math.cos(rot[1])]])
        rotZ = np.matrix([[math.cos(rot[2]), -math.sin(rot[2]), 0.0], [math.sin(rot[2]), math.cos(rot[2]), 0.0], [0.0, 0.0, 1.0]])
        rotXYZ = rotZ * rotY * rotX

        #print 'BG:'
        #print rotXYZ.transpose()

        transformed_atoms = np.matrix(copy.deepcopy(mob_atoms))
        transformed_atoms = transformed_atoms + translation
        print transformed_atoms
        transformed_atoms = transformed_atoms * rotXYZ.transpose()
       
        transformed_atoms = np.matrix.tolist(transformed_atoms)
        return transformed_atoms

    def rmsd( self, transformation, ref_atoms, mob_atoms):
        trans_atoms = self.align(transformation, mob_atoms)
        distance_sum = 0.0
        for coord in zip(ref_atoms, trans_atoms):
            distance_sum += (coord[0][0] - coord[1][0])**2
            distance_sum += (coord[0][1] - coord[1][1])**2 
            distance_sum += (coord[0][2] - coord[1][2])**2  
        score = math.sqrt(distance_sum / float(len(ref_atoms)))
        return score