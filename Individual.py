
'''
    @author: Nono Saha Cyrille Merleau 
    @email: nonosaha@mis.mpg.de 

    This class encode the structure of an RNA individual s
'''

import RNA


class Individual(object) : 


    def __init__(self, RNA_seq, RNA_structure, mfe, fitness,delta_mfe): 

        self.RNA_seq = RNA_seq 
        self.RNA_structure = RNA_structure
        self.fitness = fitness
        self.mfe = mfe
        self.delta_mfe = delta_mfe

    def get_bp_position(self, structure) : 
        position = RNA.ptable(structure)
        position = list(position[1:])

        base_paire_pos = []
        for i in range(len(position)) :
            if position[i] != 0 :
                if (position[i]-1,i) in base_paire_pos : 
          	    continue; 
                else : 
            	    base_paire_pos.append((i,position[i]-1))

        return base_paire_pos
    
def get_bp_position( structure) : 
    position = RNA.ptable(structure)
    position = list(position[1:])

    base_paire_pos = []
    for i in range(len(position)) :
         if position[i] != 0 :
             if (position[i]-1,i) in base_paire_pos : 
          	    continue; 
             else : 
            	    base_paire_pos.append((i,position[i]-1))

    return base_paire_pos
    
