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
import random
import statistics
import argparse
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt

random.seed(9001)

__author__ = "Abdelmajid Omarjee and Ali Hamraoui"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Omarjee and Hamraoui"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Omarjee and Hamraoui"
__email__ = "abdelomarjee@hotmail.com and ali.hamraoui@etu.u-paris.fr"
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
    """
    Lecture du fichier fastq
    returns : generateur de séquences
    """
    i = 0
    with open(fastq_file) as filin:
        for line in filin:
            if i%4 == 1:
                yield line[:-1]
            i += 1


def cut_kmer(read, kmer_size):
    """
    Découpade de séquences en petits Kmer
    Parametres : Kmer_size = la longueur de Kmers
    return: générateur de kmers
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i: i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """
    parameters:
        kmer_size : int, la taille des kmers
        fastq_file: str, le fichier fastq à lire
    return:
        k_mer_dict: dict, dictionnaire de k_mers avec leur occurence
    """
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
    """
    parameters:
        k_mer_dict: dict, dictionnaire de k_mers avec leur occurence
    return:
        graph: instance de la classe Digraph
    le graphe contient une arrête par k_mers pondèré par son nombre
    d'occurence
    """
    graph = nx.DiGraph()
    for key, value in kmer_dict.items():
        graph.add_edge(key[:-1], key[1:], weight=value)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
    parameters:
        graph: Digraph
        path_list : liste de chemins à supprimer
        delete_entry_note: Booléen
        delete_sink_node : Booléen
    return:
        graph: Digraph nétoyé des chemins de list_path
    """
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
    """Calcul de std
    Parameters : data, liste d'entiers
    Returns: l'ecartype d'une liste d'entiers
    """
    if len(data) > 1:
        return statistics.stdev(data)
    return 0.0


def select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=False, delete_sink_node=False):
    """La Selection de meilleur chemins parmi tous les chemins
    Parameters : graphe, paths : list de chemains,
            path_length: liste de tailles de chemins
            path_weight: liste de poids de chemins
            two boleen parameters
    Returns : graph avec just le meilleur chemin"""
    path_list_tmp = path_list.copy()
    if std(weight_avg_list) > 0.0001:
        max_weight = max(weight_avg_list)
        path_list = [path_list[i] for i in range(len(weight_avg_list)) if weight_avg_list[i] == max_weight]
        path_length = [path_length[i] for i in range(len(weight_avg_list)) if weight_avg_list[i] == max_weight]
    if std(path_length) > 0.0001:
        max_length = max(path_length)
        path_list = [path_list[i] for i in range(len(path_length)) if path_length[i] == max_length]
    best_path = path_list[random.randint(0, len(path_list) - 1)]
    path_list = [path for path in path_list_tmp if path != best_path]
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    """Caclul des poids moyens de tous les chemins
    Parameters : graphe, paths : list de paths
    Returns : la moyenne des poids de tous les chemins
    """
    weight_total = 0
    for weight in graph.subgraph(path).edges(data=True):
        weight_total += weight[2]["weight"]
    return weight_total / (len(path) - 1)


def solve_bubble(graph, ancestor_node, descendant_node):
    """supprime les bubbles dans graphe
    Parameters :graphe, ancestor noeuds, descendant noeuds
    Returns: graphe sans bubbles entre ces noeuds"""
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    length_list = []
    weight_list = []
    for path in path_list:
        weight_list.append(path_average_weight(graph, path))
        length_list.append(len(path))
    graph = select_best_path(graph, path_list, length_list, weight_list)
    return graph


def simplify_bubbles(graph):
    """supprime les bubbles dans le graphe
    Parameters : graphe
    Returns: graph sans bubbles"""
    bubble = False
    for node in graph.nodes:
        liste_predecesseurs = list(graph.predecessors(node))
        if len(liste_predecesseurs) > 1:
            for i, pred1 in enumerate(liste_predecesseurs):
                for pred2 in liste_predecesseurs[i+1:]:
                    noeud_anc = nx.lowest_common_ancestor(graph, pred1, pred2)
                    if noeud_anc is not None:
                        bubble = True
                        break
                if bubble:
                    break
            if bubble:
                break
    # La simplification ayant pour conséquence de supprimer des noeuds du hash
    # Une approche récursive est nécessaire avec networkx
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, noeud_anc, node))
    return graph


def solve_entry_tips(graph, starting_nodes):
    """supprime tous les noeuds d'entrée et garde un seul
    Parameters: graph, liste de noueds d'entrée
    Returns: graph avec un seul noued d'entrée"""
    entry = False
    for node in graph.nodes:
        length_list = []
        weight_list = []
        path_list = []
        liste_predecesseurs = list(graph.predecessors(node))
        if len(liste_predecesseurs) > 1:
            for start in starting_nodes:
                path_list += list(nx.all_simple_paths(graph, start, node))
            for path in path_list:
                weight_list.append(path_average_weight(graph, path))
                length_list.append(len(path))
            entry = True
            break
    if entry:
        graph = select_best_path(graph, path_list, length_list, weight_list, True)
        starting_nodes = get_starting_nodes(graph)
        graph = solve_entry_tips(graph, starting_nodes)
    return graph


def solve_out_tips(graph, ending_nodes):
    """supprime tous les noeuds de sortie et garde un seul
    Parameters: graph, liste de noueds de sortie
    Returns: graph avec un seul noued de sortie"""
    out = False
    for node in graph.nodes:
        length_list = []
        weight_list = []
        path_list = []
        liste_successors = list(graph.successors(node))
        if len(liste_successors) > 1:
            for end in ending_nodes:
                path_list += list(nx.all_simple_paths(graph, node, end))
            for path in path_list:
                weight_list.append(path_average_weight(graph, path))
                length_list.append(len(path))
            out = True
            break
    if out:
        graph = select_best_path(graph, path_list, length_list, weight_list, False, True)
        ending_nodes = get_sink_nodes(graph)
        graph = solve_out_tips(graph, ending_nodes)
    return graph


def get_starting_nodes(graph):
    """ extraire tous les noeud d'entrée dans un graphe
    Parameters : graph, Digraphe
    Returns : liste de noeuds d'entrée """
    l_start = []
    for node in graph:
        list_it = list(graph.predecessors(node))
        if len(list_it) == 0:
            l_start.append(node)
    return l_start


def get_sink_nodes(graph):
    """ extraire tous les noeud de sortie dans un graphe
    Parameters : graph, Digraphe
    Returns : liste de noeuds de sortie """
    l_end = []
    for node in graph:
        list_it = list(graph.successors(node))
        if len(list_it) == 0:
            l_end.append(node)
    return l_end


def get_contigs(graph, starting_nodes, ending_nodes):
    """extraire tous les contigs possible
    Paramters: graph, starting_nodes: list de noeuds d'entrées, ending_nodes: liste de noeuds de sortie
    Returns: liste de tuples avec les contigs et leurs longuers"""
    list_contig = []
    seq = ""
    for start in starting_nodes:
        for end in ending_nodes:
            if nx.has_path(graph, start, end):
                for path in nx.all_simple_paths(graph, start, end):
                    for elem in path:
                        seq += elem[0]
                    seq += elem[1:]
                    list_contig.append((seq, len(seq)))
                    seq = ""
    return list_contig


def save_contigs(contigs_list, output_file):
    """Sauver les contigs en format fasta
    Parameters: contigs_list :liste de tuples avec les contigs et leurs longuers,
        output_file: fichier de sortie
    Returns : fichier fasta avec les contigs dedans"""
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
    kmers = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmers)
    #print(get_starting_nodes(graph))
    #print(get_sink_nodes(graph))
    graph = simplify_bubbles(graph)
    # if args.graphimg_file:
    #      draw_graph(graph, args.graphimg_file)
    #print(get_starting_nodes(graph))
    #print(get_sink_nodes(graph))
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    #save contigs in out_file
    if args.output_file:
        starting_nodes = get_starting_nodes(graph)
        ending_nodes = get_sink_nodes(graph)
        contigs = get_contigs(graph, starting_nodes, ending_nodes)
        save_contigs(contigs, args.output_file)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #      draw_graph(graph, args.graphimg_file)
    # #Save the graph in file
    # if args.graphimg_file:
    #      save_graph(graph, args.graphimg_file)

if __name__ == '__main__':
    main()
