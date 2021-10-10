#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    i = 0
    with open(fastq_file) as filin:
        for line in filin:
            if i%4 == 1:
                yield line[:-1]
            i += 1


def cut_kmer(read, kmer_size):
    for i in range(len(read) - kmer_size + 1):
        yield read[i: i + kmer_size]



def build_kmer_dict(fastq_file, kmer_size):
    k_mer_dict = {}
    generator_reads = read_fastq(fastq_file)
    for read in generator_reads:
        generator_kmers = cut_kmer(read, kmer_size)
        for k_mer in generator_kmers:
            if k_mer not in k_mer_dict.keys():
                k_mer_dict[k_mer] = 1
            else:
                k_mer_dict[k_mer] += 1
    return k_mer_dict


def build_graph(kmer_dict):
    G = nx.DiGraph()
    for key, value in kmer_dict.items():
        G.add_edge(key[:-1], key[1:] , weight = value)
    return G


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif not delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path[1:])
        elif not delete_sink_node and delete_entry_node:
            graph.remove_nodes_from(path[:-1])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph

def std(data):
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list,delete_entry_node=False, delete_sink_node=False):
    path_list_tmp = path_list.copy()
    if len(weight_avg_list) >= 2 :
        max_weight = max(weight_avg_list)
        path_list = [path_list[i] for i in range(len(weight_avg_list)) if weight_avg_list[i] == max_weight]
        path_length = [path_length[i] for i in range(len(weight_avg_list)) if weight_avg_list[i] == max_weight]
    if len(path_length) >= 2 :
        max_length = max(path_length)
        path_list = [path_list[i] for i in range(len(path_length)) if path_length[i] == max_length]
    print(path_list)   
    best_path = path_list[random.randint(0, len(path_list) - 1)]
    print(best_path)
    path_list = [path for path in path_list_tmp if path != best_path]
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph
    

def path_average_weight(graph, path):
    weight_total = 0
    for weight in graph.subgraph(path).edges(data=True):
        weight_total += weight[2]["weight"]
    return weight_total / (len(path) - 1)

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = list(nx.all_simple_paths(graph, ancestor_node,descendant_node))
    length_list = []
    weight_list = []
    for path in path_list:
        weight_list.append(path_average_weight(graph, path))
        length_list.append(len(path))
    graph = select_best_path(graph, path_list, length_list, weight_list)
    return graph

def simplify_bubbles(graph):
    bubble = False 
    for nd in graph.nodes:
        # if nd != None :
        liste_predecesseurs = list(graph.predecessors(nd))
        if len(liste_predecesseurs) > 1:
            for i, pred1 in enumerate(liste_predecesseurs):
                for pred2 in liste_predecesseurs[i+1:]:
                    noeud_anc = nx.lowest_common_ancestor(graph, pred1, pred2)
                    if noeud_anc!= None:
                        bubble = True
                        break
    # La simplification ayant pour conséquence de supprimer des noeuds du hash
    # Une approche récursive est nécessaire avec networkx
    if bubble:                
        graph = simplify_bubbles(solve_bubble(graph,noeud_anc, nd))                    
    return graph

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    l_start = []
    for node in graph:
        list_it = list(graph.predecessors(node))
        if len(list_it) == 0:
            l_start.append(node)
    return l_start

def get_sink_nodes(graph):
    l_end = []
    for node in graph:
        list_it = list(graph.successors(node))
        if len(list_it) == 0:
            l_end.append(node)
    return l_end

def get_contigs(graph, starting_nodes, ending_nodes):
    list_contig = []
    seq = ""
    for start in starting_nodes:
        for end in ending_nodes:
            if nx.has_path(graph, start, end):
                for path in nx.all_simple_paths(graph, start, end):
                    for elem in path:
                        seq += elem[0]
                    seq += elem[1:]
                    list_contig.append((seq,len(seq)))
                    seq = ""
    return list_contig

def save_contigs(contigs_list, output_file):
    i = 0
    with open(output_file, "w") as filout:
        for contig in contigs_list:
            seq, length = contig
            filout.write(">contig_{} len={}\n".format(i, length))
            filout.write(fill(seq) + "\n")
            i += 1


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    g = build_graph(kmer)
    s = get_starting_nodes(g)
    e = get_sink_nodes(g)
    save_contigs(get_contigs(g, s, e), "eee")
    graph_1 = nx.DiGraph()
    graph_1.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_1 = select_best_path(graph_1, [[1,2], [3,2]], [1, 1], [5, 10], delete_entry_node=True)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()

    
