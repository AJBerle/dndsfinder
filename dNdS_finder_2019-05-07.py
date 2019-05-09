"""
This program was written by Amelia Berle as part of her independent study project in Bioinformatics at Lewis & Clark College
The goal of this program is to use genetic algorithms to find areas in DNA sequences which are experiencing selection.
dN/dS values were used as a measure of selection.

As of this writing the project is still in progress.  You can get results that match the sequence well, but it is not entirely common.
In concept, the idea seems to work, but there are a number of improvements which need to be made before this program can be put to practical use

Required packages: Bioython & Pyevolve
Written for Python 2.7
Last update: 2019-05-07
"""

import Bio.AlignIO
from pyevolve import G1DList
from pyevolve import GSimpleGA

def alignmentRead(file, file_type):
    """Reads in alignment file"""
    return (Bio.AlignIO.read(file, file_type))

def inputFile():
    """Asks user for alignment file to read in"""
    ###Version that allows the user to enter the file and filetype###
    #file = input("Which alignment file would you like to use? ")
    #file_type = input("what type is it? ")

    ###Version for if you don't want to repeatedly enter the file name###
    file = "testalignlong2.nex"
    file_type = "nexus"

    return (alignmentRead(file, file_type))

def dNdS(start, end):
    """Calculates the dN/dS score of a section of the mutations string"""

    syn = 0 #number of synonymous mutations
    non = 0 #number of nonsynonymous mutations
    ssites = 0 #number of synonymous sites
    nsites = 0 #number of nonsynonymous sites

    # counts the number of nonsyonymous changes, nonsynonymous sites, synonymous changes, and synonymous sites
    for i in range(end - start + 1):
        if (mutations[start + i] == 'N'):
            non += 1
            nsites += 3
        elif (mutations[start + i] == 'n'):
            nsites += 3
        elif (mutations[start + i] == 'S'):
            syn += 1
            nsites += 2
            ssites += 1
        elif (mutations[start + i] == 's'):
            nsites += 2
            ssites += 1

    #returns the dN/dS score.  The '+1's ensure that there is no division by zero.
    return ((float(non+1) / float(nsites+1)) / (float(syn+1) / float(ssites+1)))

def printStr(individual):
    """Prints the mutations string so that it nicely lines up with the result string"""
    print "  result string:",
    for i in range(len(individual)):
        print individual[i] + ",",

def eval_func(individual):
    """Scores individuals based on group size, number of zeros, and dNdS values"""

    ###Values for editing the fitness function###
    zerobonus = 0.75 #bonus for each zero in the individual
    shortcost = -1.5 #cost for having a short group

    fitness = 0 #The fitness of the individual
    count = len(individual)
    section = 'q' #q is just a marker to specify that a new group is starting
    start = 0 #location of the start of the group being analyzed
    end = 0 #location of the end of the group being analyzed
    groups = 0 #number of groups

    while (count > 0):
        if (section == 'q'):  # Starts a new section
            section = individual[len(individual) - count]
            start = len(individual) - count
        if (section == 0):  # benefit for having 0s (otherwise function would always favor bits turned on)
            fitness += zerobonus
        if (individual[len(individual) - count] == section):  # if this bit is part of the same group
            count -= 1
        else:
            end = len(individual) - count - 1
            if(end-start <= 3): #if the group is short
                fitness += shortcost
                groups += 1
            if (section == 1):
                # adds the difference between the dNdS values*100 and 1 for dNdS values greater than 1
                fitness += (dNdS(start, end) - 1)*10 #the *100 is so that minimal differences stand out more
                groups += 1
            elif (section == 2):
                # adds the difference between the dNdS values*100 and 1 for dNdS values less than 1
                fitness += (1 - dNdS(start, end))*10 #the *100 is so that minimal differences stand out more
                groups += 1
            else:
                fitness -= (1 - abs(dNdS(start, end)))*10 #punishes the individual for groups which are turned off, but have dNdS values away from 1
            section = 'q'  # Resets section to start a new one
        if (count == 0):  # for the final codon.  This could be shortened by combining this and the else statement above, but was not a high priority and ran out of time
            end = len(individual) - count - 1
            if(end-start <= 3): #if the group is short
                fitness += shortcost
                groups += 1
            if (section == 1):
                # adds the difference between the dNdS values*100 and 1 for dNdS values greater than 1
                fitness += (dNdS(start, end) - 1)*100 #the *100 is so that minimal differences stand out more
                groups += 1
            elif (section == 2):
                # adds the difference between the dNdS values*100 and 1 for dNdS values less than 1
                fitness += (1 - dNdS(start, end))*100 #the *100 is so that minimal differences stand out more
                groups += 1
            else:
                fitness -= (1 - abs(dNdS(start, end))) #punishes the individual for groups which are turned off, but have dNdS values away from 1
    if (groups > 0):
        fitness = fitness / groups
        if(fitness<0):
            fitness = 0 #pyevolve cannot handle negative fitnesses
        return fitness
    else:
        return 0

if __name__ == "__main__":
    ###Values for modifying the genetic algorithm###
    numgens = 10000 #number of generations you would like to run
    numindls = 100 #number of individuals per generation
    printfreq = 10 #after this many generations, print the stats for current generation

    # this is a list of the first two sites for the codons in which changing the 3rd site cannot change the amino acid (AKA, the 3rd site is synonymous)
    non2 = ["ct", "gt", "tc", "cc", "ac", "gc", "cg", "gg"]
    mutations = []
    alignment = inputFile()
    for i in range(len(alignment[1]) // 3):  # reads codons out
        for c in alignment:
            if(c.format("fasta")!=alignment[0].format("fasta")):
                if(c[3*i:3*i+3].seq==alignment[0,3*i:3*i+3].seq): #if codons the same
                    if(c[3*i:3*i+2].seq in non2): #If the third site is synonymous
                        mutations.append('s')
                    else: #If the third site is nonsynonymous
                        mutations.append('n')
                elif(c[3*i:3*i+3].seq.translate()==alignment[0, 3*i:3*i+3].seq.translate()): #Synonymous mutation
                    mutations.append('S')
                else: #Nonsynonymous mutation
                    mutations.append('N')
    genome = G1DList.G1DList(len(mutations)) #Genome instance
    genome.setParams(rangemin=0, rangemax=2) #populates individuals with 0s, 1s, or 2s
    genome.evaluator.set(eval_func) #Sets the fitness function we wrote as the fitness function for pyevolve
    ga = GSimpleGA.GSimpleGA(genome) #sets up the genetic algorithm
    ga.setGenerations(numgens)
    ga.setPopulationSize(numindls)
    ga.evolve(freq_stats=printfreq) #run GA
    print(ga.bestIndividual())
    printStr(mutations)