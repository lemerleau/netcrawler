import pandas 
import Individual
import random
import numpy 
import RNA
import collections
import pp 
import os
import platform
import subprocess
import multiprocess
import time
import argparse
import itertools
import matplotlib.pyplot as plt 
import multiprocessing as mp 

class Landscape (object) : 

    def __init__ (self, target) : 
        self.target = target 
    
    def fitness (self, structure) : 
        return 1./(1+RNA.hamming_distance(self.target, structure))

    def ens_defect(self, sequence) : 
        p = subprocess.Popen("defect", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_input = sequence+"\n"+self.target+"\n"
        
        defect, error = p.communicate(cmd_input)
        defect = defect.split("\n")
        
        return 1/float(defect[-3])


def getHairepinCoord(sequence, target) : 
    cmd_input = str(sequence)+"\n"+str(target)+"\n"
    os.system("echo '"+cmd_input+"'| RNAeval -v |grep Hairpin | tr -d A-Z-a-z-'('-')' | cut -d ':' -f 1 > hairpin.loop")
    
    hairpins = []
    with open("hairpin.loop", "r") as loops : 
        data_r = loops.read().split("\n")
        for dt in data_r : 
            if dt != '' : 
                hairpins.append(tuple(numpy.array(map(int, dt.strip().split(',')))-1))
    return hairpins  

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

def getInteriorCoord(sequence, target) : 
    cmd_input = str(sequence)+"\n"+str(target)+"\n"
    os.system("echo '"+cmd_input+"' | RNAeval -v | grep Interior\ loop | tr -d A-Z-a-z-'('-')'|cut -d ':' -f 1 > interior.loop")
    
    with open("interior.loop", "r") as loops : 
        data_r = loops.read().split(";")
        list_ = []
        for dt in data_r : 
            list_+=filter(None,dt.strip().split("\n"))
    return numpy.array(list(set([tuple(map(int,string.strip().split(","))) for string in list_]))) -1

def boost_hairpins(sequence,target, coord) :
    
    if (coord[1]-coord[0]-1) > 100 : 

        with open("DB3/Hp/Hp"+str(coord[1]-coord[0]-1)+".txt", 'r') as f: 
            Hps = f.read()
            hps = Hps.split("\n")
            hp = "".join(numpy.random.choice(hps[1:],1)[0])
            stems = list(hp) 
    else : 
        stems = list(numpy.random.choice(["A"],len(range(coord[1]-coord[0]-1))))
    if coord[1]-coord[0]-1 >=4 : 
        stems[0] = "G"
        stems.insert(0,"G")
        stems.append("C")
        sequence[coord[0]:coord[1]+1] = stems
    else : 
        sequence[coord[0]+1:coord[1]] = stems
    
    
    return sequence

def boostMulti(sequence, coord, pos) : 
    #sequence = list(sequence) 
    if (coord[1]-coord[0]-1) <100 : 
        with open("DB3/Ml/Ml"+str(coord[1]-coord[0]+1)+".txt", 'r') as f: 
            Mls = f.read()
            mls = Mls.split("\n")
            ml = "".join(numpy.random.choice(mls[1:],1)[0])
            if ml =="" :
                print "Empty"
                ml = "".join(numpy.random.choice(mls[1:],1)[0])
            stems = list(ml) 
        sequence[coord[0]:coord[1]+1] = stems
    
    if coord[1]-1 not in pos["p_table"] : 
        sequence[coord[1]-1] = "G"
    
    return sequence

def boostInterior(sequence, coord, pos) :
    
    if [coord[0]+1,coord[1]-1] not in pos["interior"] : 
        sequence[coord[0]+1] = "G"
    
    if coord[1]+1 not in pos["p_table"]: 
        if coord[1]+1 < len(sequence) : 
            sequence[coord[1]+1] = "G"   
    return sequence 

def pphamming(listOfStructures, landscape) : 
    pool = multiprocess.Pool(multiprocess.cpu_count())
    dists = pool.map(landscape.fitness,listOfStructures)
    pool.close()
    return dists

def ppens_defect(listOfSequences, landscape) : 
    pool = multiprocess.Pool(multiprocess.cpu_count())
    dists = pool.map(landscape.ens_defect,listOfSequences)
    pool.close()
    return dists

#Mutation function 
def mutateOne(seq, mut_probs,target,pos, p_n, p_c, mut_bp=.1) :  # best = 0.5, 0.1, 1/6
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    nucleotides = ["A", "G","U","C"]
    p_table = pos["bp_pos"]
    RNA_seq = numpy.array(list(seq))
    r = numpy.random.rand(len(seq))
    mut_pos =RNA_seq[r<mut_probs] 
    choices = numpy.random.choice(nucleotides, len(mut_pos), p=p_n)
    RNA_seq[r<mut_probs] = choices 
    apply = []
    for bp_cord in p_table : 
        r = random.uniform(0,1)
        
        if r < mut_bp : 
            bp = numpy.random.choice(base_paire,1, p=p_c)
            
            RNA_seq[bp_cord[0]] = bp[0][0]
            RNA_seq[bp_cord[1]] = bp[0][1]
        
        if bp_cord in pos["hairpins"] : 
                RNA_seq = boost_hairpins(RNA_seq,target, bp_cord)
        """
        elif bp_cord in pos["interior"] : 
                RNA_seq = boostInterior(RNA_seq,target, bp_cord)
        
        elif bp_cord  in pos["multi"] : 
            RNA_seq = boostMulti(RNA_seq, bp_cord,pos)
     
        elif bp_cord in pos["interior"] : 
            RNA_seq = boostInterior(RNA_seq, bp_cord, pos)
        """

    return ''.join(RNA_seq)


def mutateAll(population, mut_prob, target,pos, p_n, p_c) : 
    mutated_pop = [] 
    for individual in population : 
        mutated_pop.append(mutateOne(individual,mut_prob,target,pos,p_n, p_c))
    return mutated_pop



def bt_save_population(prev_pop, population,gen, root_path) : 
    data = []
    prev_data = []
    for i in range(len(population)) : 
        data.append([population[i].RNA_seq, population[i].RNA_structure,population[i].mfe, population[i].fitness])
        prev_data.append([prev_pop[i].RNA_seq, prev_pop[i].RNA_structure,prev_pop[i].mfe, prev_pop[i].fitness])

    
    dataFrame = pandas.DataFrame(data)
    prev_dataFrame = pandas.DataFrame(prev_data)
    prev_dataFrame.to_csv(root_path+"/prev_gen"+str(gen)+".csv")
    dataFrame.to_csv(root_path+"/gen"+str(gen)+".csv")


def save_population(population,gen, root_path) : 
    population.to_csv(root_path+"/gen"+str(gen)+".csv")


def ppeval(listOfSeqs, target, task) : 
    with open("rnaeval_in"+str(task), "w") as file_ : 
        for seq in listOfSeqs : 
            file_.write(seq.strip()+"\n"+target.strip()+"\n")
        file_.close()
          
    os.system("RNAeval -j -P vrna185x.par --infile=rnaeval_in"+str(task)+"  |tr -d A-Z,'(',')'|cut -d ' ' -f 2- > result_"+str(task))
    with open("result_"+str(task), "r") as file_ : 
        eval_ = file_.read().split()
    os.remove("result_"+str(task))
    os.remove("rnaeval_in"+str(task))
    return list(numpy.array(eval_, dtype=float))

def ppfold(listOfSeqs,task) : 
    dataFrame= pandas.DataFrame(listOfSeqs)
    dataFrame.to_csv("sequences"+str(task),sep=" ",index=False, header=False)    
    os.system("RNAfold -j -P vrna185x.par --infile=sequences"+str(task)+" --noPS > rnafold_result"+str(task)) #To change the energy par just add -P vrna185x.par  to RNAfold
    os.system("cut -d ' ' -f 1 rnafold_result"+str(task)+" | tr -d A-Z > result_"+str(task))
    os.system("cat rnafold_result"+str(task)+"|tr -d A-Z,'(',')' | cut -d ' ' -f 2- >mfes"+str(task))
    with open("result_"+str(task), "r") as file_ : 
        strc_ = file_.read().split()
        
    with open("mfes"+str(task), "r") as file_ : 
        mfes= file_.read().split()
     
    os.remove("result_"+str(task))
    
    return strc_, list(numpy.array(mfes, dtype=float))

def eval_proportion_selection(population, size) : 
    
    
    #ens_defects = numpy.array(population["Defects"], dtype = float)

    evals = numpy.array(population["Evals"], dtype = float)
    mfes  = numpy.array(population["Mfes"], dtype = float)
    delta_mfe = evals - mfes 
    delta_mfe = delta_mfe / sum(delta_mfe)
    weights = (1-delta_mfe)**100
    
    
    #selected = numpy.random.choice(list(population["RNA_sequence"]),size=size,p=weights/sum(weights))
    
    choices = list(population["RNA_sequence"])
    choices = numpy.array(choices)

    #weights = list(population["Mfes"])
    weights = list(weights)
    return choices[weights.index(max(weights))] , max(weights)


def ensDefect_proportion_selection(population, size, target) : 
    
    ensDefect = [ens_defect(ind.RNA_seq, target) for ind in population]
    selected = numpy.random.choice(population,size=size,p=numpy.array(ensDefect)/sum(ensDefect))
    return selected



def reproduce(population, size) : 
    list_fitness = []
    #ref = numpy.random.choice(a=[".",len(population[0].RNA_seq)])
    for ind in population : 
        list_fitness.append(ind.fitness)
    list_fitness = sorted(list_fitness, reverse=True) 

    sorted_pop = [ ] 
    for fitness in list_fitness : 
        for ind in population : 
            if ind.fitness == fitness : # or self.landscape.hamming_distance(ind.RNA_structure, ref) == 0 : 
                sorted_pop.append(ind)
                
    return sorted_pop[:size]

def mutateOnePoint(seq): 
    nucleotides = ["A", "C", "G" , "U"]
    mutate_pos = numpy.random.randint(low=0, high=len(seq))
    RNA_seq = list(seq)
    RNA_seq[mutate_pos] = numpy.random.choice(nucleotides, 1)[0]
    return "".join(RNA_seq)

def compute_relative_degree(sequence, nb_mutant=100) :
    
    target = ppfold([sequence], 1)[0]
    
    mutants = [ mutateOnePoint(sequence) for i in range(nb_mutant)]
    strs,_ = ppfold(mutants, 1)
    
    return strs.count(target)/(3.*len(sequence))

def select_on_relative_degree(pop, size) : 
    #degrees = map(compute_relative_degree, pop)
    """
    pool = mp.Pool(mp.cpu_count())
    degrees= list(pool.map(compute_relative_degree, pop))
    pool.close()
    """
    degrees= [compute_relative_degree(seq) for seq in pop]
    pop = numpy.array(pop)

    return pop[degrees.index(min(degrees))] 



def simple_EA(landscape, number_of_generation, mut_probs, init_pop, selection_method, log_folder,pos,p_n, p_c) : 
    
    print (" Starting of evolution ")
    prev_population = init_pop.copy(deep=True) #Initialize the population of RNA
    
    root_path = "../Logs/Defect/Test100/"+str(log_folder)+"/"+selection_method
    try:
        os.makedirs(root_path)
    except OSError :
        print (" Can not initialize the log folder ")
    
    population_size =  len(init_pop)
    n = number_of_generation
    max_fitness = max(numpy.array(prev_population["Fitness"], dtype = float))
    
    while (n > 0) and (max_fitness < 1.):
        
        if (number_of_generation - n)%10 == 0 : 
            print ('Generation '+str(number_of_generation - n)), "Max fitness = ", max_fitness
        newgeneration = []
    
        selected_ind = eval_proportion_selection(prev_population,population_size)
        
        mutated = mutateAll(selected_ind,mut_probs,landscape.target,pos, p_n, p_c)
        newgeneration.append(mutated)

        structures, mfes = ppfold(mutated,log_folder)
        newgeneration.append(structures)
        newgeneration.append(mfes)
        newgeneration.append(pphamming(structures,landscape))
        evals = ppeval(mutated, landscape.target,log_folder)
        newgeneration.append(evals)

        #defects = ppens_defect(mutated, landscape)
        #newgeneration.append(defects)
        
        prev_population = pandas.DataFrame(numpy.array(newgeneration).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
        new_max = max(numpy.array(prev_population["Fitness"], dtype=float))
        if new_max > max_fitness : 
            max_fitness = new_max
        n -=1

    return prev_population


def netcrawler(landscape, number_of_generation, mut_probs, init_pop, selection_method, log_folder,pos,p_n, p_c, wind_size) : 
    
    root_path = "Logs/degree_analysis/"+str(selection_method)+'/'+str(log_folder)
    try:
        os.makedirs(root_path)
    except OSError as error:
        #print "ERROR", error 
        pass

    print (" Starting of evolution ")
    init_sequence = ''.join(numpy.random.choice(['A', 'C', 'G','U'],len(landscape.target))) #Initialize the population of RNA
    current_structure, current_mfe = RNA.fold(init_sequence)
    print "initial seq :", init_sequence
    best_pop = []
    best_pop.append(init_sequence)
    
    w_i = landscape.fitness(current_structure)
    n = 0
    v_i = 1.
    oup = 1.
    mutants = []
    mutants.append(init_sequence)
    eval_ = 1.
    epoche = 0
    evals = []
    epoches = []
    epoches.append([epoche, w_i, v_i, v_i/eval_, eval_])
    fitness_data = []
    fitness_data.append([0., w_i])
    t= 0.
    neutral_set = []
    neutral_set.append(init_sequence)
    mut_bp = 0.05
    while (n<number_of_generation) and (w_i < 1.):
        
        if (number_of_generation - n)%10 == 0 : 
            print ('Generation '+str(number_of_generation - n)), "Max fitness = ", w_i
        
        for i in range(wind_size) :
            
            mutant = mutateOne(init_sequence, mut_probs,landscape.target,pos, p_n, p_c, mut_bp)
            mutants.append(mutant)
            mutant_strc, mutant_mfe = RNA.fold(mutant)
            w_j = landscape.fitness(mutant_strc)
            eval_ +=1
            t = t+1

            if w_j > w_i : 
                epoche +=1
                init_sequence = mutant
                current_structure = mutant_strc
                current_mfe = mutant_mfe
                w_i = w_j
                if t%100 == 0 : 
                    fitness_data.append([t/100,w_i])
                 
                print "v_i = ",v_i, "w_j = ", w_j, "v_obs = ", v_i/eval_, " total of mutants = ", eval_
                epoches.append([epoche, w_j,v_i,  v_i/eval_, eval_])
                evals.append(eval_)
                
                eval_ = 1.
                v_i = 1.
                oup = 1.
                mutants = []
                mutants.append(init_sequence)

                neutral_set = []
                neutral_set.append(init_sequence)
                
        
            elif w_j == w_i : 
                #if not np.array_equal(init_seq, mutant) : 
                init_sequence = mutant
                v_i +=1.
                neutral_set.append(init_sequence)
                if t%100 == 0 : 
                    fitness_data.append([t/100,w_i])
               
            else : 
                if t%100 == 0 : 
                    fitness_data.append([t/100,w_i])
                continue
            best_pop.append(init_sequence)
            """
            if len(neutral_set) > 100 : 
                print "Select one neighbor"
                init_sequence = select_on_relative_degree(neutral_set,1)
                neutral_set=[]
                #mut_bp = mut_bp + 0.05
            """
            
            if len(mutants) > 100 : 

                evals = ppeval(mutants, landscape.target,log_folder)
                strcs, mfes = ppfold(mutants,log_folder)
                fitnesses = numpy.array(pphamming(strcs,landscape))

                pop = pandas.DataFrame(numpy.array([mutants, strcs, mfes, fitnesses, evals]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
                init_sequence,current_mfe = eval_proportion_selection(pop,1)
                mutants =[]
            
            
        n +=1
        if n<=0 or w_i ==1.  : 
            epoches.append([epoche, w_i,v_i,  v_i/eval_, eval_])
    #print len(neutral_set), neutral_set
    save_population(pandas.DataFrame(best_pop), "best", root_path)
    return init_sequence, current_structure, current_mfe, w_i , epoches, fitness_data


def run_netcrawler(number_of_generation,pop_size, mut_probs, log_folder,landscape, p_n, p_c) : 
    target = landscape.target 
    p_table = get_bp_position(target)
    init_depth =len(target)
   
    nucleotides = ["A", "U" , "G", "C"]
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    pop = []
    i = 0
    
    
    while i < pop_size : 
        
        if i < 4: 
            arn = numpy.random.choice(nucleotides[i:i+1],init_depth)
        else : 
            arn = numpy.random.choice(nucleotides,len(target))

        for bp_cord in p_table : 
            bp = numpy.random.choice(base_paire,1)
            arn[bp_cord[0]] = bp[0][0]
            arn[bp_cord[1]] = bp[0][1]
        pop.append(''.join(arn))
        i = len(pop)
    
    pos = {
        "p_table" : numpy.array(RNA.ptable(target)) - 1,
        "bp_pos"  : p_table,
        "hairpins": getHairepinCoord(pop[0],target),
        "interior": getInteriorCoord(pop[0],target), 
        "multi"   : getMultiCoord(pop[0], target)}
    print len(pos["multi"]), " multi loop(s)"
    print len(pos["interior"]), " Interior loop(s)"
    print len(pos["hairpins"]), " Hairpin loop(s)"
    print len(pos["bp_pos"])," loop(s) in total"
    #init_pop = pandas.DataFrame(numpy.array([pop, strcs, mfes, fitnesses, evals]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
    tic = time.time()
    init_sequence, current_struc, current_mfe, w_i,  epoches, fitness_data  = netcrawler(landscape,number_of_generation, mut_probs,pop, "net_craw",log_folder, pos, p_n, p_c,10)
    toc = time.time()
    
    #for ind in best_pop.values : 
    #   if float(ind[3]) >= 0.4 : 
    #        print ind
    print "Time = ", toc-tic, "For method : ", "EVAL"

    print init_sequence, current_struc, current_mfe, w_i
    print  epoches

    return epoches, fitness_data
         


def run(number_of_generation,pop_size, mut_probs, log_folder,landscape, p_n, p_c) : 
    target = landscape.target 
    p_table = get_bp_position(target)
    init_depth =len(target)
   
    nucleotides = ["A", "U" , "G", "C"]
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    pop = []
    i = 0
    
    while i < pop_size : 
        
        if i < 4: 
            arn = numpy.random.choice(nucleotides[i:i+1],init_depth)
        else : 
            arn = numpy.random.choice(nucleotides,len(target))

        for bp_cord in p_table : 
            bp = numpy.random.choice(base_paire,1)
            arn[bp_cord[0]] = bp[0][0]
            arn[bp_cord[1]] = bp[0][1]
        pop.append(''.join(arn))
        i = len(pop)
    evals = ppeval(pop, target,log_folder)
    strcs, mfes = ppfold(pop,log_folder)
    fitnesses = pphamming(strcs, landscape)
    #ens_defects = ppens_defect(pop, landscape)
    pos = {
        "p_table" : numpy.array(RNA.ptable(target)) - 1,
        "bp_pos"  : p_table,
        "hairpins": getHairepinCoord(pop[0],target),
        "interior": getInteriorCoord(pop[0],target), 
        "multi"   : getMultiCoord(pop[0], target)}
    print len(pos["multi"]), " multi loop(s)"
    print len(pos["interior"]), " Interior loop(s)"
    print len(pos["hairpins"]), " Hairpin loop(s)"
    print len(pos["bp_pos"])," loop(s) in total"
    init_pop = pandas.DataFrame(numpy.array([pop, strcs, mfes, fitnesses, evals]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
    tic = time.time()
    best_pop = simple_EA(landscape,number_of_generation, mut_probs,init_pop, "EVAL",log_folder, pos, p_n, p_c)
    toc = time.time()
    
    for ind in best_pop.values : 
        if float(ind[3]) >= 0.4 : 
            print ind
    print "Time = ", toc-tic, "For method : ", "EVAL"
    
        
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
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('--target', type=str)
    parser.add_argument('--job', type=int)
    parser.add_argument('-g', type=int, help="Number of generation")
    parser.add_argument('-n', type=int, help="Population Size")
    
    args = parser.parse_args()
    target = args.target
    landscape = Landscape(target)
    print "====================================================================================================="
    print "Solving for the target = ", target, "Length = ",len(target) 
    print "====================================================================================================="
    number_of_run = args.job
    init_depth =len(target)
    mut_prob = 1./init_depth
    number_of_generation = args.g
    pop_size = args.n
    p_n = [0.7,0.1,0.1,.1] #default = [0.25,0.25,0.25,.25] [0.25,0.65,0.05,.05] [0.7,0.1,0.1,.1] ["A", "G","U","C"]
    p_c = [0.4, 0.5, 0.1, 0.,0.,0.] #[0.2,0.2,0.1,0.1,0.2,0.2] #[0.4, 0.5, 0.1, 0.,0.,0.] ["GC","CG","AU","UA", "GU", "UG"]
    ppservers = ()

    mut_probs = numpy.array(RNA.ptable(target)[1:])
    mut_probs = mut_probs + mut_prob
    mut_probs[mut_probs>mut_prob] = 0


    job_server = pp.Server(4, ppservers=ppservers)
    
    print "Start running job", number_of_run
    run_netcrawler(number_of_generation, pop_size, mut_probs, 1, landscape, p_n, p_c)
    
   #obs = [(task , job_server.submit(run_netcrawler, (number_of_generation,pop_size, mut_probs, task,landscape, p_n, p_c), (ppeval,ppfold, netcrawler,  getInteriorCoord,boostInterior,getMultiCoord, boostMulti, pphamming,mutateAll,get_bp_position, eval_proportion_selection,getHairepinCoord, boost_hairpins, mutateOne,  save_population, simple_EA, ppens_defect),("numpy", "Individual", "RNA", "random","pandas","os", "time", "collections","subprocess","multiprocess"))) for task in range(number_of_run)]
    """
    netcrawler_data = []
    for task, job in jobs : 
        netcrawler_data.append(job()[1])

    netcrawler_data = numpy.array(netcrawler_data)
    print netcrawler_data.shape
    mean_fitness = numpy.mean(netcrawler_data, axis=0)
   
    
    plt.plot(mean_fitness[:,0], mean_fitness[:,1],)
    plt.xlabel("Generation")
    plt.ylabel("Max Fitness")
    plt.savefig("necrawler2.eps")
    
    plt.show() 
    """
    
	
if __name__ == "__main__":
    main()
