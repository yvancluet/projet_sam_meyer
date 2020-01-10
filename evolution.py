###############################################################################
#############################        IMPORTS    ###############################
###############################################################################

import TSC as simul
import numpy as np
import math
import pandas as pd
import random as rd
import sys
import os
import shutil
import matplotlib.pyplot as plt
import csv


#PRINCIPE
#a chaque iteration, des mutations aleatoires sont executees (indel ou inversion)
#durant ces mutations, les fichiers contenant les positions des genes, sections
#codantes et non codantes etc, (present dans 'evol') sont lus sous forme de 
#dataframes, qui seront modifies par la mutation
#les fichiers sont regeneres dans le dossier 'evol2'
#on calcule la fitness des nouvelles valeurs
#si elle est meilleure qu'avant la mutation, alors la version mutee devient la
#nouvelle version du genome (les fichiers d'evol2 sont copies dans evol)
#meme si la fitness est moins bonne, l'update peut quand meme se faire avec une probabilite

#ERREURS :
#les fonctions de mutation marchent bien MAIS:
#proba de garder une fitness inférieure systématiquement = 1
#premiere ligne de tous genes identiques est tej durant la simulation --> erreurs
#vérifier le calcul de la fitness



###############################################################################
#############################      VARIABLES    ###############################
###############################################################################

#params.ini a ete modifie, travaille sur le dossier "evol" qui est une copie de
#"tousgenesidentiques", afin de continuer de faire évoluer le meme genome.

output_dir="output/" 
input_file = "params.ini"
input_folder = "evol/"

#Tous les parametres sont arbitraires (a réévaluer)

#parametre d'evolution (pondere la probabilite de garder une fitness inférieure a la fitness actuelle
q = 0.00005
#probabilité d'indel (proba d'inversion = 1-pi)
pid = 0.9

#unité du programme (nb de nucleotides)
unit = 60

#temps entre chaque sauvegarde
savetime = 50

#total genome size
genome_size = 30000

#amount of iterations
n_iter = 1000

#gff head text
def gff_head() :
    global genome_size
    return "##gff-version 3\n"+ \
    "#!gff-spec-version 1.20\n"+ \
    "#!processor NCBI annotwriter\n"+ \
    "##sequence-region\ttousgenesidentiques\t1\t"+str(genome_size)+"\n" +\
    "tousgenesidentiques\tRefSeq\tregion\t1\t"+str(genome_size)+"\t.\t+\t.\tID=id0;Name=tousgenesidentiques"
#pb : genome_size n'est pas toujours a jour (surtout si on fait un seul indel)

###############################################################################
#############################      FUNCTIONS    ###############################
###############################################################################

#copies files from evol2 in evol   
def Evol2ToEvol(): 
    files = ['evol2/prot.dat', 'evol2/TSS.dat', 'evol2/TTS.dat', 'evol2/tousgenesidentiques.gff']
    for f in files:
        shutil.copy(f, 'evol/')

#used to create folder
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

#used to add line at beginning of line
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


#gets result of transcription
def get_expression():

    tr_info = pd.read_csv(output_dir+"all_tr_info.csv", sep='\t')
    target = pd.read_csv("environment.dat", sep = "\t", header = None)
    target.columns = ["TU", "tr_rate"]
    return tr_info, target

#calculate fitness
def fitness():
    
    #get observed expression and target expression rates
    tr_info, target = get_expression()
    
    #only keep columns of interest
    f_cible = target["tr_rate"]
    f_obs = tr_info["number of transcripts generated"].values
    
    #normalize observed expression rate
    f_obs = f_obs/f_obs.sum(axis=0,keepdims=1)

    #number of genes
    nb_genes = len(f_cible)

    #calculate fitness
    fitness = 0
    for i in range(nb_genes):
        fitness += math.fabs(math.log(f_obs[i]/f_cible[i]))
    fitness = math.exp(-fitness)
    
    return fitness

#select random position and check that it's not problematic (in a gene)
#works for now, check again after inversions etc
def select_mut_pos(tss, tts, prot, tousgenesid, safety_mes):
    pos = rd.randint(1, genome_size)
    trunc = pos%unit
    pos -= trunc
    #checks if position is correct
    correct=True
    for i in range(tousgenesid.shape[0]):
        start = tousgenesid.at[i, "start"]
        end = tousgenesid.at[i, "end"]
        if start > end:
            a = end
            end = start
            start = a
        if pos > start - safety_mes and pos < end + safety_mes:
            correct = False
            break
    if correct == False:
        pos = select_mut_pos(tss, tts, prot, tousgenesid, safety_mes)
    return pos

#read all files from "evol" folder
def read_from_folder():
    tousgenesidentiques = pd.read_csv(input_folder+"tousgenesidentiques.gff", sep="\t", skiprows = 4, header = 0,
                        names = ["seqname", "source", "feature", "start", "end",
                                 "score", "strand", "frame", "attribute"])
    prot = pd.read_csv(input_folder+"prot.dat", sep="\t", header = 0)
    TTS = pd.read_csv(input_folder+"TTS.dat", sep="\t", header = 0)
    TSS = pd.read_csv(input_folder+"TSS.dat", sep="\t", header = 0)

    return tousgenesidentiques, prot, TTS, TSS

#writes files for cell after mutation to 'evol2'
def write_to_folder(tousgenesidentiques, TSS, TTS, prot):
    global genome_size
    tousgenesidentiques.to_csv("evol2/tousgenesidentiques.gff", sep="\t", index = False, header = False)
    line_prepender("evol2/tousgenesidentiques.gff", gff_head())
    TSS.to_csv("evol2/TSS.dat", sep="\t", index = False)
    TTS.to_csv("evol2/TTS.dat", sep="\t", index = False)
    prot.to_csv("evol2/prot.dat", sep="\t", index = False)

    

def indel():
    global genome_size
    #chooses an insertion or a deletion
    mut_type = np.random.choice(["insertion", "deletion"], 1, p = [0.5, 0.5])
    print (mut_type)
    
    #reads files
    tousgenesidentiques, prot, TTS, TSS = read_from_folder()
    
    #pick a random position
    pos = select_mut_pos(TSS, TTS, prot, tousgenesidentiques, unit)
    print (pos)

    
    #define wether we will delete or add 60 nucleotides
    if mut_type == "insertion" : change = unit
    elif mut_type == "deletion" : change = -1*unit
    else :
        print ("error : incorrect mutation type")
        sys.exit()
    genome_size += change
    
    for index, row in tousgenesidentiques.iterrows():
        if min(row["start"], row["end"]) > pos:
            tousgenesidentiques.at[index, "start"] += change
            tousgenesidentiques.at[index, "end"] += change
    
    for index, row in TSS.iterrows():
        if row["TSS_pos"] > pos:
            TSS.at[index, "TSS_pos"] += change
            
    for index, row in TTS.iterrows():
        if row["TTS_pos"] > pos:
            TTS.at[index, "TTS_pos"] += change
            
    for index, row in prot.iterrows():
        if row["prot_pos"] > pos:
            prot.at[index, "prot_pos"] += change
    
    print ("genome size : ", genome_size)
    
    #write new files to folder
    write_to_folder(tousgenesidentiques, TSS, TTS, prot)

    #works for the most part
def inversion():
    print ("inversion")
    #explication position
    #1(barriere 1) -- 1001 (TSS) -- 2000 (TTS) -- 3000 (fin) -- 3001(barriere 2)
    
    #read all files from "evol" folder
    tousgenesidentiques, prot, TTS, TSS = read_from_folder()
    
    #pick start and end positions for inversions
    inv_start = select_mut_pos(TSS, TTS, prot, tousgenesidentiques, unit)
    inv_end = select_mut_pos(TSS, TTS, prot, tousgenesidentiques, unit)
    print (inv_start, inv_end)
    
    if inv_start > inv_end:
        a = inv_end
        inv_end = inv_start
        inv_start = a
    
    #do the coordinates inversions
    for index, row in tousgenesidentiques.iterrows():
        if min(row["start"], row["end"]) > inv_start and max(row["start"], row["end"]) < inv_end:
            tousgenesidentiques.at[index, "start"] = inv_end - (row["start"] - inv_start)
            tousgenesidentiques.at[index, "end"] = inv_start + (inv_end - row["end"])
            if row["strand"] == "+" : tousgenesidentiques.at[index, "strand"] = "-"
            elif row["strand"] == "-" : tousgenesidentiques.at[index, "strand"] = "+"
            else : 
                print("strand problem in tousgenesidentiques !")
                sys.exit()
    
    for index, row in TSS.iterrows():
        if row["TSS_pos"] > inv_start and row["TSS_pos"] < inv_end:
            TSS.at[index, "TSS_pos"] = inv_end - (row["TSS_pos"] - inv_start)
            if row["TUorient"] == "+" : TSS.at[index, "TUorient"] = "-"
            elif row["TUorient"] == "-" : TSS.at[index, "TUorient"] = "+"
            
    for index, row in TTS.iterrows():
        if row["TTS_pos"] > inv_start and row["TTS_pos"] < inv_end:
            TTS.at[index, "TTS_pos"] = inv_end - (row["TTS_pos"] - inv_start)
            if row["TUorient"] == "+" : TTS.at[index, "TUorient"] = "-"
            elif row["TUorient"] == "-" : TTS.at[index, "TUorient"] = "+"
            
    for index, row in prot.iterrows():
        if row["prot_pos"] > inv_start and row["prot_pos"] < inv_end:
            prot.at[index, "prot_pos"] = inv_end - (row["prot_pos"] - inv_start)
    
    print ("genome size : ", genome_size)
    #copy results to evol2    
    write_to_folder(tousgenesidentiques, TSS, TTS, prot)
    



###############################################################################
#############################      MAIN         ###############################
###############################################################################

def _main():
    fitness_list=[0]
    mut_list = []

    #iterate
    for i in (range (n_iter)):
        print("\n\nITERATION : ", i)
        #select random mutation
        mut_type = np.random.choice([True, False], 1, p = [pid, 1-pid])
        if mut_type :
            indel()
            mut_list.append("indel")
        else :
            inversion()
            mut_list.append("inv")
        
        #transcribe new genome and calculate fitness
        simul.start_transcribing("params_evol2.ini", output_dir)
        fit = fitness()
        
        #if fitness is better, keep new genome
        if fit > fitness_list[-1]:
            print ("new_fitness : ", fit)
            fitness_list.append(fit)
            Evol2ToEvol()
            
        #otherwise, only keep it with probabilty
        else:
            df = fit-fitness_list[-1]
            prob = math.exp(df/q)
            print ("prob", prob, "df : ", df, "q : ", q, "fitness_list : ", fitness_list) 
            if rd.uniform(0,1) < prob:
                fitness_list.append(fit)
                Evol2ToEvol()
                #copier evol2 dans evol
            else:
                fitness_list.append(fitness_list[-1])
    
        if i%savetime == 0:
            with open('fitness_q_'+str(q)+'_pid_'+str(pid)+'.csv','w') as f:
                writer = csv.writer(f);
                writer.writerows(zip(fitness_list, mut_list));
    
    #ca marche pas le plot 
    plt.plot(fitness_list)   
    plt.show()
    




if __name__ == "__main__":
    _main()
