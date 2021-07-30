# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 10:00:59 2021
@author: CaiGroup
"""

#You will need to have the packages gc, os, copy, random, logging, numpy, and pandas installed

#The GeneWalk package can be installed with: pip install genewalk

#You will only need to do this once, but if you want to force an update to GeneWalk, just re-install it

#Running GeneWalk for the first time will create a folder labled genewalk where all future results will be sent

#This is located in your main drive under the currently installing User

#All nessecary information to run the program "should" be in this file

#For more information about GeneWalk and its inputs, go to the GeneWalk site: https://churchman.med.harvard.edu/genewalk/#tutorial

#The GeneWalk site gives examples on running GeneWalk from the command line, but they should still be appliacble to using this script
    
    
##############Import Packages and functions##############
import gc
import os
import copy
import random
import logging
import numpy as np
import pandas as pd
import genewalk
from genewalk.nx_mg_assembler import load_network
from genewalk.gene_lists import read_gene_list
from genewalk.deepwalk import run_walks
from genewalk.null_distributions import get_rand_graph, get_null_distributions
from genewalk.perform_statistics import GeneWalk
from genewalk import logger as root_logger, default_logger_format, default_date_format
from genewalk.resources import ResourceManager
from genewalk.plot import GW_Plotter
from genewalk.cli import create_folder, save_pickle, load_pickle


##############Initiate Logger##############
logger = logging.getLogger('genewalk.cli')


##############Set base folder##############
default_base_folder = os.path.join(os.path.expanduser('~/'), 'genewalk')


##############Argumants##############

#Set values for GeneWalk
class args: 
    project = "project"
    
    genes = "genes.csv"
    
    id_type = "ensembl_id"

    stage = "all"
    
    base_folder = default_base_folder
    
    network_source = "pc"
    
    network_file = None
    
    nproc = 1
    
    nreps_graph = 3
    
    nreps_null = 3
    
    alpha_fdr = 1
    
    dim_rep = 8
    
    save_dw = False
    
    random_seed = None

##############User Inputs##############
print("")
print("Please add user inputs. If you want to keep the default values, press enter. Program will beep to let you know its done.")

proj = str(input("Folder name to save results to(If no folder with that name exists, a new folder will be created): "))
    
gen = str(input("File containing the gene Ids: "))

check_ids = False

while check_ids == False:
    
    ids = str(input("Format of gene ids(values are hgnc_symbol, hgnc_id, ensembl_id ,mgi_id ,rgd_id, entrez_human, entrez_mouse, custom): "))
    
    if ids not in ["hgnc_symbol", "hgnc_id", "ensembl_id" , "mgi_id", "rgd_id", "entrez_human", "entrez_mouse", "custom", '']:
        
        print("")
        
        print("Invalid Id format")
        
    else:
        
        check_ids = True

check_stag = False

while check_stag == False:
    
    stag = str(input("Stage of GeneWalk do you wnt to run(values are all, node_vectors, null_distribution, statistics, visual): "))
    
    if stag not in ["all", "node_vectors", "null_distribution", "statistics", "visual", '']:
        
        print("")
        
        print("Invalid stage selection")
        
    else:
        
        check_stag = True

check_soc = False

while check_soc == False:
    
    net_soc = str(input("Network type to use(values are pc, indra, edge_list, sif, sif_annot, sif_full): "))

    if net_soc == "indra":
        
        net_fil = str(input("Python pickle file contating list of INDRA statements: "))
        
        check_soc = True
    
    elif net_soc == any(["edge_list", "sif", "sif_annot", "sif_full"]):
        
        net_fil = str(input("File contating network: "))
        
        check_soc = True
        
    elif net_soc in ['', "pc"]:
        
        net_file = None
        
        check_soc = True
        
    else:
        
        print("")
        
        print("Invalid network source")
        
check_npc = False

while check_npc == False:
    
    npc = str(input("number of processors to use(recomended number is 4): "))
    
    if npc == '':
        
        check_npc = True
    
    else:
    
        try:
            npc = int(npc)
            
            check_npc = True
            
        except ValueError:
            
            print("")
        
            print("Please enter an integer value")
            

check_rep = False

while check_rep == False:
    
    rep = str(input("Number of repeats to run when calculating node vectors on GeneWalk graph: "))
    
    if rep == '':
        
        check_rep = True
        
    else:
        
        try:
            rep = int(rep)
            
            check_rep = True
            
        except ValueError:
            
            print("")
        
            print("Please enter an integer value")
            
        except npc == '':
            
            check_rep = True
    
check_neps = False

while check_neps == False:
    
    neps = str(input("Number of repeats to run when calculating node vectors on the random network graphs for constructing null distribution: "))
    
    if neps == '':
        
        check_neps = True
        
    else:
        
        try:
            neps = int(neps)
        
            check_neps = True
        
        except ValueError:
        
            print("")
            
            print("Please enter an integer value")

check_fdr = False

while check_fdr == False:
    
    fdr = str(input("false discovery rate to use when outputting final statistics table: "))
    
    if neps == '':
        
        check_fdr = True
        
    else:
        
        try:
            fdr = float(fdr)
        
            check_fdr = True
        
        except ValueError:
        
            print("")
            
            print("Please enter a number value")

check_dim = False

while check_dim == False:
    
    dim = str(input("Dimension of vector representations: "))
    
    if dim == '':
        
        check_dim = True
        
    else:
        
        try:
            dim = int(dim)
        
            check_dim = True
        
        except ValueError:
        
            print("")
            
            print("Please enter an integer value")
            
save = str(input("Do you want to save the full file for each repeat(Yes/No): "))

check_ran = False

while check_ran == False:

    ran = str(input("Optional random seed: "))
    
    if ran == '':
        
        check_ran = True
        
    else:
        
        try:
            
            ran = float(ran)
            
            check_ran = True
            
        except ValueError:
            
            print("")
            
            print("Not a viable random seed")
            
args.project = proj

args.genes = gen
        
if ids == '':

    pass

else:
    
    args.id_type = ids

if stag == '':

    pass

else:
    
    args.stage = stag
    
if net_soc == '':

    pass

else:
    
    args.network_source = net_soc
    
    args.network_file = net_file
    
if npc == '':
    
    pass

else:
    
    args.nproc = npc
    
if rep == '':
    
    pass

else:
    
    args.nreps_graph = rep
    
if neps == '':
    
    pass

else:
    
    args.nreps_null = neps
    
if fdr == '':
    
    pass

else:
    
    args.alpha_fdr = fdr
    
if dim == '':
    
    pass

else: args.dim_rep = dim
    
if save in ['', "No"]:
    
    pass

else:
    
    args.save_dw = save
    
if ran == '':
    
    pass

else: 
    args.random_seed = ran
        
  
##############Processing############## 
def main():   
    
    #Create folder for this instance of GeneWalk
    project_folder = genewalk.cli.create_folder(args.base_folder, args.project)
    
    # Add a logger specific to the project and processing stage
    log_file = os.path.join(project_folder, 'genewalk_%s.log' % args.stage)
    
    formatter = logging.Formatter(default_logger_format,
                                  datefmt = default_date_format)
    
    project_log_handler = logging.FileHandler(log_file)
    
    project_log_handler.setFormatter(formatter)
    
    root_logger.addHandler(project_log_handler)
    
    # Make sure a network file was provided for custom network sources
    if args.network_source in {'indra', 'sif', 'sif_annot', 'sif_full',
                               'edge_list'}:
        if not args.network_file:
            raise ValueError('The network_file argument must be provided'
                             ' when using network_source %s.' %
                             args.network_source)
            
    # Make sure SIF network is provided when gene ID type is set to "custom"
    if args.id_type == 'custom':
        if args.network_source not in {'sif_annot', 'sif_full'}:
            raise ValueError('When using id_type custom, the network_source'
                             ' has to be either sif_annot or sif_full, with '
                             'the network_file argument pointing to a custom '
                             'SIF file.')
    
    #Set random seed is supplied
    if args.random_seed:
        logger.info('Running with random seed %d' % args.random_seed)
        
        random.seed(a = int(args.random_seed))
    
    #Instantiate the resource manager
    rm = ResourceManager(base_folder = args.base_folder)
    
    #Random walks of supplied genes
    if args.stage in ('all', 'node_vectors'):
        
        #Read in lists of genes
        genes = read_gene_list(args.genes, args.id_type, rm)
        
        #pickle the genes
        save_pickle(genes, project_folder, 'genes')
        
        #load network for genes
        MG = load_network(args.network_source, args.network_file, genes,
                          resource_manager = rm)
        
        #pickle the network
        save_pickle(MG.graph, project_folder, 'multi_graph')
        
        for i in range(args.nreps_graph):
            logger.info('%s/%s' % (i + 1, args.nreps_graph))
            
            #Perform the Walks
            DW = run_walks(MG.graph, workers = args.nproc, size = args.dim_rep)
    
            #Pickle the node vectors and DW object if desired
            if args.save_dw:
                save_pickle(DW, project_folder, 'deepwalk_%d' % (i + 1))
                
            nv = copy.deepcopy(DW.model.wv)
            
            save_pickle(nv, project_folder,
                        'deepwalk_node_vectors_%d' % (i + 1))
    
            # Delete the DeepWalk object to clear memory
            del DW, nv
            gc.collect()
    
    #Random walks for null graph
    if args.stage in ('all', 'null_distribution'):
        
        #Get base for random null graph
        MG = load_pickle(project_folder, 'multi_graph')
        
        #pre-allocate for null distribution
        srd = []
        
        for i in range(args.nreps_null):
            logger.info('%s/%s' % (i + 1, args.nreps_null))
            
            #Form radom null graph
            RG = get_rand_graph(MG)
            
            #Perform the walks
            DW = run_walks(RG, workers=args.nproc, size = args.dim_rep)
    
            #Pickle the node vectors and DW object for null is desired
            if args.save_dw:
                save_pickle(DW, project_folder, 'deepwalk_rand_%d' % (i + 1))
            nv = copy.deepcopy(DW.model.wv)
            save_pickle(nv, project_folder, 'deepwalk_node_vectors_rand_%d'
                                            % (i + 1))
            # Delete the DeepWalk object to clear memory
            del DW
            gc.collect()
    
            # Calculate the null distributions
            srd += get_null_distributions(RG, nv)
            del nv
            gc.collect()
            
        #Gather null distribution
        srd = np.asarray(sorted(srd))
        
        #Pickle the similarity distances
        save_pickle(srd, project_folder, 'genewalk_rand_simdists')
    
    ##############Results##############
    
    #Get stattistics for walks
    if args.stage in ('all', 'statistics'):
        
        #Gene network
        MG = load_pickle(project_folder, 'multi_graph')
        
        #list of genes
        genes = load_pickle(project_folder, 'genes')
        
        #Node vectors of the network
        nvs = [load_pickle(project_folder,
                           'deepwalk_node_vectors_%d' % (i + 1))
               for i in range(args.nreps_graph)]
        
        #null distribution
        null_dist = load_pickle(project_folder, 'genewalk_rand_simdists')
        
        #Create GeneWalk object for computation
        GW = GeneWalk(MG, genes, nvs, null_dist, gene_id_type = args.id_type)
        
        #Tabulate for output
        df = GW.generate_output(alpha_fdr = args.alpha_fdr)
        
        #Name result file
        fname = os.path.join(project_folder, 'genewalk_results.csv')
        
        logger.info('Saving final results into %s' % fname)
        
        #Save to .csv
        df.to_csv(fname, index=False, float_format = '%.3e')
    
    #Create visuals for statsitsics
    if args.stage in ('all', 'visual'):
        fname = os.path.join(project_folder, 'genewalk_results.csv')
        
        dGW = pd.read_csv(fname)
        
        figure_folder = create_folder(project_folder, 'figures')
        
        create_folder(figure_folder, 'barplots')
        
        GWp = GW_Plotter(figure_folder, dGW, args.alpha_fdr)
        
        GWp.generate_plots()
            
if __name__ == '__main__':
    main()