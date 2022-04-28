#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Time    :   2022/04/15
@Author  :   Quincey Lyu 
'''
# Modified from https://github.com/neocaleb/TENET

from jpype import *
import numpy
import os
import datetime
import argparse
import scanpy as sc
# import click

# Set environment variables
os.environ['JAVA_HOME']="/TJPROJ6/SC/personal_dir/lvguangqi/software/jdk-18"

# Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument('--h5ad', help='A h5ad file contains results of PAGA')
group.add_argument('--filelist', nargs=3, help='Space seprated files: expression.csv trajectory.txt cell_select.txt')
parser.add_argument('--hist', type=int, default=1, help='History length, default=1')
parser.add_argument('--out', default='.', help='Output dir, default=current dir')

# Parse arguments
args = parser.parse_args()
h5ad = args.h5ad
if args.filelist:
    filelist = args.filelist
    expr = args.filelist[0]
    traj = args.filelist[1]
    cell = args.filelist[2]
hist = int(args.hist)
out = args.out

# Functions: pre-process input data

def process_h5ad(h5ad):
    adata=sc.read_h5ad(h5ad)
    data_matrix=adata.X.todense().A

    # Cell select
    branch = []
    fileOut = open(f'{out}/PAGAcell_select.txt', 'w')

    for i in range(len(adata.obs_names)):
        branch.append(i)
        fileOut.write('1\n')

    fileOut.close()

    # Trajectory
    trajectory1 = []
    fileOut = open(f'{out}/PAGApseudotime.txt', 'w')
    for pseudotime in adata.obs['dpt_pseudotime']:
        trajectory1.append(pseudotime)
        fileOut.write(f'{str(pseudotime)}\n')

    fileOut.close()

    trajectory1 = numpy.array(trajectory1)
    trajectory1 = trajectory1[branch == 1]
    trajectory1SortIndex = numpy.argsort(trajectory1)

    # Expression
    gene_name = list(adata.var_names)
    cell_gene_all = numpy.transpose(data_matrix)

    return branch, trajectory1SortIndex, cell_gene_all, gene_name


def process_file_list(expr, traj, cell):
    # Process cell select file
    branch = []
    with open(cell) as fileIn:
        for line in fileIn:
            branch.append(int(line.rstrip().split()[0]))

    branch = numpy.array(branch)

    # Process trajectory file
    trajectory1 = []

    with open(traj) as fileIn:
        for line in fileIn:
            trajectory1.append(float(line.rstrip().split()[0]))

    trajectory1 = numpy.array(trajectory1)
    trajectory1 = trajectory1[branch == 1]
    trajectory1SortIndex = numpy.argsort(trajectory1)

    # Process expression matrix file
    cell_gene_all = []
    with open(expr) as fileIn:
        gene_name = fileIn.readline().rstrip().split(',')[1:]

        for line in fileIn:
            temp = line.rstrip().split(',')
            cell_gene_all_temp = []

            for i in range(1, len(temp)):
                cell_gene_all_temp.append(float(temp[i]))
            
            cell_gene_all.append(cell_gene_all_temp)
        cell_gene_all = numpy.transpose(numpy.array(cell_gene_all))
    
    return branch, trajectory1SortIndex, cell_gene_all, gene_name


def main(branch, trajectory1SortIndex, cell_gene_all, gene_name, out):
    # Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=/TJPROJ6/SC/personal_dir/lvguangqi/software/TENET/pkgs/infodynamics.jar","-Xmx16G")

    # Main function
    # List pairs
    list_pairs=[]
    for geneIndex1 in range(len(cell_gene_all)):
        for geneIndex2 in range(len(cell_gene_all)):
            if geneIndex1<geneIndex2:
                list_pairs.append([geneIndex1,geneIndex2])
    list_pairs=numpy.array(list_pairs)

    TEresult=[None] * len(list_pairs)
    for num_pair in range(len(list_pairs)):
        expression_data = cell_gene_all[list_pairs[num_pair,0]][numpy.newaxis]
        expression_data = numpy.append(expression_data, cell_gene_all[list_pairs[num_pair,1]][numpy.newaxis], axis=0)
        expression_data1=[]

        for i in range(len(expression_data)):
            data_temp1 = numpy.array(expression_data[i])
            data_temp1 = data_temp1[branch == 1]
            data_temp1 = data_temp1[trajectory1SortIndex]
            data_temp1 = list(data_temp1)
            expression_data1.append(data_temp1)
    
        expression_data = expression_data1
        del expression_data1

        # Create a TE calculator and run it:
        teCalcClass = JPackage("infodynamics.measures.continuous.kernel").TransferEntropyCalculatorKernel
        teCalc = teCalcClass()
        teCalc.setProperty("NORMALISE", "true") # Normalise the individual variables
        teCalc.initialise(hist, 0.5) # Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
        resultTemp = []
        teCalc.setObservations(JArray(JDouble, 1)(expression_data[0]), JArray(JDouble, 1)(expression_data[1]))
        resultTemp.append(teCalc.computeAverageLocalOfObservations())
        teCalc.initialise(hist, 0.5) # Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
        teCalc.setObservations(JArray(JDouble, 1)(expression_data[1]), JArray(JDouble, 1)(expression_data[0]))
        resultTemp.append(teCalc.computeAverageLocalOfObservations())
        TEresult[num_pair] = numpy.ndarray.tolist(list_pairs[num_pair,:]) + resultTemp
        if (num_pair % int(len(list_pairs) / 2)) == 0:
            print(datetime.datetime.now())
        else:
            print('Warning: number of pairs is an odd number')

    TEmatrix=[]
    for i in range(len(gene_name)):
        TEmatrixTemp=[]
        for j in range(len(gene_name)):
            TEmatrixTemp.append(0)
        TEmatrix.append(TEmatrixTemp)

    for i in range(len(TEresult)):
        temp=list(TEresult[i])
        TEmatrix[int(temp[0])-1][int(temp[1])-1]=float(temp[2])
        TEmatrix[int(temp[1])-1][int(temp[0])-1]=float(temp[3])
    
    fileOut = open(f'{out}/TE_result_matrix.txt', 'w')
    fileOut.write("TE")
    for i in range(len(gene_name)):
        fileOut.write("\t"+gene_name[i])
    for i in range(len(gene_name)):
        fileOut.write("\n"+gene_name[i])
        for j in range(len(gene_name)):
            fileOut.write("\t"+str(TEmatrix[i][j]))
    fileOut.close()


if __name__ == "__main__":
    if args.filelist:
        process_file_list(expr ,traj, cell)
    else:
        process_h5ad(h5ad)

    main()
