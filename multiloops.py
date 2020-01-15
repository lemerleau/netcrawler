import numpy 
import RNA
import pandas
import os


def getMultiCoord(sequence, target) : 
    cmd_input = str(sequence)+"\n"+str(target)+"\n"
    os.system("echo '"+cmd_input+"'| RNAeval -v |grep Multi | tr -d A-Z-a-z-'('-')' | cut -d ':' -f 1 > multi.loop")
    
    multi = []
    with open("multi.loop", "r") as loops : 
        data_r = loops.read().split("\n")
        for dt in data_r : 
            if dt != '' : 
                multi.append(tuple(numpy.array(map(int, dt.strip().split(',')))-1))
    return multi 

def get_bp_position(structure) : 
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


def main() : 
    df = pandas.read_csv("eterna_notfound.csv", sep=',')
    targets_with_multiloops = {}
    targets = [".((...(((...((....))..(((...)))..)))....))"]
    print len(targets[0])
    for target in targets : 
        p_table = get_bp_position(target)
        init_depth =len(target)
    
        nucleotides = ["A", "U" , "G", "C"]
        base_paire = ["GC","CG","AU","UA", "GU", "UG"]
       
        arn = numpy.random.choice(nucleotides,init_depth)
        for bp_cord in p_table : 
            bp = numpy.random.choice(base_paire,1)
            arn[bp_cord[0]] = bp[0][0]
            arn[bp_cord[1]] = bp[0][1]
        multiloop_coord = getMultiCoord("".join(arn),target)
        if len(multiloop_coord) : 
            targets_with_multiloops[target] = multiloop_coord
    
    print targets_with_multiloops, len(targets_with_multiloops)

        



if __name__=="__main__": 
    main()