import pandas as pd 
import numpy as np
import RNA
#import networkx as nx 
import matplotlib.pyplot as plt
#import graph_tool.all as gt 
from scipy.stats import levy


def mutate(sequence, mu) : 
    indexes_to_mutate = list(np.where(np.random.binomial(2,mu, len(sequence))==1)[0])
    mutant = list(sequence)
    for i in indexes_to_mutate: 
        mutant[i] = np.random.choice(['A', 'C','G','U'], 1)[0]
    return "".join(mutant)

def fitness(sequence, target): 
    strc = RNA.fold(sequence)[0]
    return 1./(1.+RNA.hamming_distance(strc, target))

def evaluate(pop,target): 
    return [fitness(seq, target) for seq in pop]
    

def select(pop, weights, size) : 
    
    return np.random.choice(pop,size=size, p=np.array(weights)/sum(weights))

def levy_rdv(a,b,size) : 
    levy_dat = levy.rvs(size=size)
    levy_ab = (b-a)*((levy_dat - min(levy_dat))/(max(levy_dat)-min(levy_dat))) + a
    return np.array(levy_ab, dtype=int)

def levy_mutation(pop) : 
    mutants = []
    nb_points = levy_rdv(1,len(pop[0]),len(pop))
    for i in range(len(pop)) : 
        indexes_to_mutate = np.random.choice(range(len(pop[i])), nb_points[i], replace=False)
        mutant = np.array(list(pop[i]))
        mutant[indexes_to_mutate] = np.random.choice(['A', 'C','G','U'], nb_points[i])
        
        mutants.append("".join(mutant))
        
    return mutants


def rna_evo(target,init_pop, t, mu) : 
    t_n = 0 
    current_pop = np.copy(init_pop)
    fitnesses = evaluate(current_pop, target)
    while t_n < t and max(fitnesses) < 1. : 
          
        mutants = [mutate(seq, mu) for seq in current_pop]
        fitnesses = evaluate(mutants, target)
        selected = select(mutants,fitnesses, len(mutants))
        
        current_pop = np.copy(selected)
        t_n +=1 
        
        if t_n%10==0 :
            print ("Generation ", t_n , " with max fitness ", max(fitnesses))
          
    return current_pop[fitnesses.index(max(fitnesses))]

def rna_evo_with_levy(target,init_pop, t) : 
    t_n = 0 
    current_pop = np.copy(init_pop)
    fitnesses = evaluate(current_pop, target)
    while t_n < t and max(fitnesses) < 1.: 
          
        mutants = levy_mutation(current_pop)
        fitnesses = evaluate(mutants, target)
        selected = select(mutants,fitnesses, len(mutants))
        
        current_pop = np.copy(selected)
        t_n +=1 
        
        if t_n%10==0 :
            print ("Generation ", t_n , " with max fitness ", max(fitnesses))
          
    return current_pop[fitnesses.index(max(fitnesses))]



def main() : 
    target = '((....))..((....))'
    pop = []

    for i in range(100): 
        rna = np.random.choice(['A','C','G','U'],len(target))
        pop.append(''.join(rna))

    print(pop)
    

if __name__=="__main__" : 
    main()
