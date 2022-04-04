import pandas as pd
from collections import Counter
import collections
import numpy as np
from typing import List, Dict
import time

class Node:
    """define a node class. A node has at most 2 children (left and right), a name,
    and a data. A node with no children is called a leaf"""
    def __init__(self, name, left = None, right = None, rawdata = None, data = None, metadata = None):
        self.rawdata = rawdata
        self.data = data
        self.metadata = metadata
        self.name = name
        self.left = left
        self.right = right

def to_dictionary(node : Node) -> Dict:
    node_list = get_node_list(node)
    node_names = [node.name for node in node_list]
    node_dict = {name:node for (name,node) in zip(node_names, node_list)}
    return node_dict

def is_leaf(node: Node) -> bool:
    if(not node.left and not node.right):
        return True
    else:
        return False

def get_leaf_rawdata(root: Node) -> Dict:
    """returns a dictionary where keys are leaf names and values are leaf raw data"""
    # If node is null, return
    leafdict = {}
    if (not root):
        return leafdict
    
    # If node is leaf node,
    # key --> leaf name, value --> leaf data
    if (not root.left and not root.right):
        leafdict[root.name] = root.rawdata
        #leaflist.append(root.name)
        #print(root.data,
        #	end = " ")
        return leafdict

    # If left child exists,
    # check for leaf recursively
    if root.left:
        temp = get_leaf_rawdata(root.left)
        leafdict.update(temp)

    # If right child exists,
    # check for leaf recursively
    if root.right:
        temp = get_leaf_rawdata(root.right)
        leafdict.update(temp)
    
    return leafdict

def get_leaf_data(root: Node) -> Dict:
    """returns a dictionary where keys are leaf names and values are leaf data"""
    # If node is null, return
    leafdict = {}
    if (not root):
        return leafdict
    
    # If node is leaf node,
    # key --> leaf name, value --> leaf data
    if (not root.left and not root.right):
        leafdict[root.name] = root.data
        #leaflist.append(root.name)
        #print(root.data,
        #	end = " ")
        return leafdict

    # If left child exists,
    # check for leaf recursively
    if root.left:
        temp = get_leaf_data(root.left)
        leafdict.update(temp)

    # If right child exists,
    # check for leaf recursively
    if root.right:
        temp = get_leaf_data(root.right)
        leafdict.update(temp)
    
    return leafdict

def get_leaf_list(root: Node) -> List:
    """returns a list of node objects (leaves)"""
    leaflist = []
    # If node is null, return
    if (not root):
        return leaflist

    # If node is leaf node,
    # get the node
    if (not root.left and not root.right):
        leaflist.append(root)
        #leaflist.append(root.name)
        #print(root.data,
        #	end = " ")
        return leaflist

    # If left child exists,
    # check for leaf recursively
    if root.left:
        leaflist.extend(get_leaf_list(root.left))

    # If right child exists,
    # check for leaf recursively
    if root.right:
        leaflist.extend(get_leaf_list(root.right))

    return leaflist

def get_node_list(root: Node) -> List:
    "returns a list of node objects (in order traversal)"
    #if node is null, return
    nodelist = []
    if (not root):
        return nodelist

    #left -> root -> right
    else:
        nodelist.extend(get_node_list(root.left))
        nodelist.append(root)
        nodelist.extend(get_node_list(root.right))
        return nodelist

def get_inner_node_list(root: Node) -> List:
    """get inner nodes (nodes that are not leaves) of the binary tree"""
    nodelist = get_node_list(root)
    leaflist = get_leaf_list(root)
    innernodelist = [k for k in nodelist if k not in leaflist]
    return innernodelist

def df_same(df1: pd.DataFrame, ind1: int, df2: pd.DataFrame, ind2: int) -> List:
    """return indexes of rows of two columns where two elements of the columns at the
    row are equal"""
    colnames1 = list(df1.columns)
    colnames2 = list(df2.columns)
    newdf = pd.concat([df1, df2], axis=1) #concatanate dfs to subset them together
    subdf = newdf.loc[newdf[colnames1[ind1]] == newdf[colnames2[ind2]]]
    idx = list(subdf.index)
    return idx

def df_equalto(df: pd.DataFrame, idx: int, ele) -> List:
    """return indexes of rows of a column where the element of the column at the
    row is equal to given element"""
    colnames = list(df.columns)
    subdf = df.loc[df[colnames[idx]] == ele]
    idx = list(subdf.index)
    return idx    
    pass

def df_diff(df1: pd.DataFrame, ind1: int, df2: pd.DataFrame, ind2:int) -> List:
    """return indexes of rows of two columns where two elements of the columns at the
    row are not equal"""
    colnames1 = list(df1.columns)
    colnames2 = list(df2.columns)
    newdf = pd.concat([df1, df2], axis=1) #concatanate dfs to subset them together
    subdf = newdf.loc[newdf[colnames1[ind1]] != newdf[colnames2[ind2]]]
    idx = list(subdf.index)
    return idx

# def inner_check(datadict: Dict) -> List:
#     """this function will take a dictionary of dataframes. Columns of dataframes are
#     0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
#     1. residues are the same
#     2. status is "C" or conserved (frequency above 0.9).
#     We will iterate over keys 2 at a time and check for the given conditions.
#     Function returns the indexes of said residues."""
#     #get the keys of the dictionary
#     keys = sorted(datadict.keys())
#     #create a list to hold indexes that satisfy given conditions
#     listofidx = []
#     #if there is only single dataframe, return indexes of the df
#     if len(keys) == 1:
#         dfsingle = datadict[keys[0]]
#         idx = df_equalto(dfsingle, 3, "C") #list(dfsingle.index)
#         idx.sort()
#         return idx
#     else: #more than one dataframe
#         for i in range(len(keys) - 1): #iterate over two elements at a time
#             keypair = keys[i:i+2]
#             df1 = datadict[keypair[0]] #first df
#             df2 = datadict[keypair[1]] #second df
#             temp =[df_same(df1, 1, df2, 1), df_equalto(df1, 3, "C"), df_equalto(df2, 3, "C")]
#             idx = set(temp[0]).intersection(*temp)
#             listofidx.append(idx)
#         #return the indexes that satisfy the conditions for all dataframes at the same time
#         commonels = list(set.intersection(*[set(x) for x in listofidx]))
#         commonels.sort()
#         return commonels

# def divergence(datadict1: Dict, datadict2: Dict) -> List:
#     """this function will take two dictionaries of dataframes. Columns of dataframes are
#     0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
#     1. residues are not the same
#     2. status is "C" or conserved (frequency above 0.9).
#     We will iterate over keys 2 at a time and check for the given conditions.
#     Function returns the indexes of said residues."""
#     listofidx = []
#     for value1 in datadict1.values(): #iterate over two dicts
#         for value2 in datadict2.values():
#             temp = [df_diff(value1, 1, value2, 1), df_equalto(value1, 3, "C"), df_equalto(value2, 3, "C")]
#             idx = set(temp[0]).intersection(*temp)
#             listofidx.append(idx)
#     #return the indexes that satisfy the conditions for all dataframes at the same time
#     commonels = list(set.intersection(*[set(x) for x in listofidx]))
#     commonels.sort()
#     return commonels

# def emergence(datadict1: Dict, datadict2: Dict) -> List:
#     """this function will take two dictionaries of dataframes. Columns of dataframes are
#     0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
#     1. residues are not the same
#     2. status is "C" or conserved for first dict and "NC" or nonconserved
#     for second dict. (frequency above 0.9). We will iterate over keys 2 at a time and 
#     check for the given conditions. Function returns the indexes of said residues."""
#     listofidx = []
#     for value1 in datadict1.values():
#         for value2 in datadict2.values():
#             temp = [df_diff(value1, 1, value2, 1), df_equalto(value1, 3, "C"), df_equalto(value2, 3, "NC")]
#             idx = set(temp[0]).intersection(*temp)
#             listofidx.append(idx)
#     #return the indexes that satisfy the conditions for all dataframes at the same time
#     commonels = list(set.intersection(*[set(x) for x in listofidx]))
#     commonels.sort()
#     return commonels

# def outer_check(datadict1: Dict, datadict2: Dict) -> List:
#     """this function will take two dictionaries of dataframes. Columns of dataframes are
#     0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
#     1. residues are not the same
#     2. status is "C" or conserved for first dict (frequency above 0.9).
#     We will iterate over keys 2 at a time and 
#     check for the given conditions. Function returns the indexes of said residues."""
#     listofidx = []
#     for value1 in datadict1.values():
#         for value2 in datadict2.values():
#             temp = [df_diff(value1, 1, value2, 1), df_equalto(value1, 3, "C")]
#             idx = set(temp[0]).intersection(*temp)
#             listofidx.append(idx)
#     #return the indexes that satisfy the conditions for all dataframes at the same time
#     commonels = list(set.intersection(*[set(x) for x in listofidx]))
#     commonels.sort()
#     return commonels

# def chisquare(observed, expected):
#     """This function takes two lists and perform chi square test using lists,
#     for info on chi-square test, visit https://www.statisticshowto.com/probability-and-statistics/chi-square/.
#     Lists should have same length, categories should share indexes, degree of freedom is equal to the length of lists.
#     returns a tuple (first element chi value, second degree of freedom)."""
#     chival = 0
#     for i in range(len(observed)):
#         a = ((observed[i] - expected[i])**2)/expected[i]
#         chival += a 
#     return (chival, len(observed) - 1)

def read_fasta(fasta_path):
    """This function reads the fasta file from specified path
    and returns all the fastas as a list"""

    fasta_seqs = []

    with open(fasta_path, 'r') as file:
        counter = 0
        temp = ""
        for line in file:
            if line[0] == '>':
                counter = counter + 1
                if counter == 2:
                    fasta_seq = temp
                    if fasta_seq != "":
                        fasta_seqs.append(fasta_seq)
                    temp = ""
                    counter = 1				
                temp += line
            elif line[0] != ">":
                temp += line
        fasta_seq = temp
        if fasta_seq != "":
            fasta_seqs.append(fasta_seq)

    return fasta_seqs

def fasta_to_dataframe(fasta_seqs):
    """This function takes a list of fasta sequences and
    turn them into a pandas dataframe where columns represent sequences and rows
    contain amino acid sequences at that position. It must be that every sequence
    in the list is in equal length"""

    list_of_fastas = []
    for seq in fasta_seqs:
        split_seq = seq.split("\n")
        header = split_seq[0]
        body = "".join(map(str,split_seq[1:]))
        body = list(body)
        list_of_fastas.append({"header": header, "body": body})

    df = pd.DataFrame()
    for fasta in list_of_fastas:
        df[fasta["header"]] = fasta["body"]

    return df

def Most_Common(lst):
    """this function returns the most common element in a list"""

    #count every element with counter
    data = Counter(lst)
    return data.most_common(1)[0][0]

def consensus(df, limit):
    """This function takes a pandas dataframe that has the multiple sequence
    alignment fasta data. It defines consensus as any amino acid that appears 
    above given limit in the position.It does not take gaps	into consideration.
    IMPORTANT: ADD A CONDITION FOR TOO MANY GAPS. For now if gaps are above 50%"""

    #create an empty dataframe to add consensus sequence with the frequency of the residues 
    consensus_seq = pd.DataFrame(columns = ["residue", "frequency", "status"])
    #get the dimensions of original dataframe to iterate over rows
    dims = df.shape

    for i in range(dims[0]):
        #get the values of the row (each row has the residue at that position for every protein)
        x = df.iloc[i].astype(str).tolist()
        #count the gaps
        gapcount = x.count("-")
        #count how many proteins there are
        poslength = len(x)
        #if more than 50% is gaps, there is no consensus
        if (gapcount/poslength) >= 0.5:
            consensus_seq = consensus_seq.append({"residue" : "-", "frequency" : 0.0, "status" : "NC"}, ignore_index = True)
        else:
            #get rid of the gaps
            y = [ elem for elem in x if elem != "-"]
            #get most common residue, possible consensus
            mostcommonres = Most_Common(y)
            #count how many times most common residue occurs
            mostcommonrescount = y.count(mostcommonres)
            #we take residues above 50% frequency so we can later use them (I CHANGED THE LIMIT TO 0.34 FROM 0.5) -- reverted back the change 
            if (mostcommonrescount/len(y)) > (0.5):
                #if most common residue is below limit it is not conserved (NC)
                if (mostcommonrescount/len(y)) < (limit):
                    consensus_seq = consensus_seq.append({"residue" : mostcommonres, "frequency" : mostcommonrescount/len(y), "status" : "NC"}, ignore_index = True)
                #if most common residue equal or above limit it is conserved
                else:
                    consensus_seq = consensus_seq.append({"residue" : mostcommonres, "frequency" : mostcommonrescount/len(y), "status" : "C"}, ignore_index = True)
            #if residues is below 50% it is not conserved
            else:
                consensus_seq = consensus_seq.append({"residue" : "-", "frequency": 0.0, "status" : "NC"}, ignore_index = True)

    return consensus_seq

def new_consensus(df: pd.DataFrame, consensus_limit: float, gap_tolerance: float):
    """This function takes a pandas dataframe that has the multiple sequence
    alignment fasta data. It defines consensus as any amino acid that appears 
    above given limit in the position.It does not take gaps	into consideration (gap_tolerance
    is defined by the user).gap tolerance should be between 0 and 1 (ideally <=0.5).
    consensus limit should be between 0 and 1 (ideally >0.5)"""
    
    df_transposed = df.transpose() #transpose the table because pd.value_counts counts columns not rows
    count_table = df_transposed.apply(pd.value_counts) #count the occurance of amino acids ("-" gap)
    
    gap_limit = (df.shape[1])*gap_tolerance #how many gaps are tolerated
    filter = count_table.iloc[0] > gap_limit #which positions have gap more than allowed
    count_table.loc[:,filter] = np.NaN #replace positions above limit wiht NaN
    count_table.iloc[0] = np.NaN #get rid of counts of gaps since we don't take them into account
    total_counts = count_table.sum(axis = 0, skipna = True) #how many positions are left after we cull the gaps
    frequency_table = count_table.div(total_counts, axis = "columns") #frequency of each amino acid
    #start_time = time.time()
    top_aa = frequency_table.idxmax(axis = 0) #get the most frequent amino acid
    top_freq = frequency_table.max(axis = 0) #get the frequency of most frequent amino acid
    #print("My program took", time.time() - start_time, "seconds to run")
    consensus_table = pd.concat([top_aa, top_freq], axis=1) #combine the amimo acids and frequencies
    consensus_table.columns = ["residue", "frequency"] #rename the columns to match old consensus function
    
    #label the positions as consensus or not consensus with respect to given limit (consensus_limit)
    consensus_table["status"] = (consensus_table["frequency"] >= consensus_limit).map({True : "C", False : "NC"})
    #residues with frequency below gap tolerance are discarded
    consensus_table.loc[consensus_table["frequency"]<gap_tolerance, ["residue", "frequency"]] = "-", 0
    values = {"residue" : "-", "frequency": 0} 
    consensus_table.fillna(value=values, inplace=True) # replace NaN with "-" or 0 depending on the column
    
    return consensus_table

