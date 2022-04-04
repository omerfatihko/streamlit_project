import pandas as pd
import numpy as np
import streamlit as st
import time
# from os import listdir #X_OK, access, ctermid, path
# from os import path
from os.path import join #isfile
import my_functions as mf
from typing import List #Dict
start_time = time.time()

#ascii representation of class B1 family tree with node names
st.image("ascii_tree.png", use_column_width=True)

#required directories
fastadir = join("data", "fastas")
csvdir = join("data", "csvs")

#class B1 information table
classB1 = pd.read_csv(join(csvdir, "class_B1_all_consensus_seqs.csv"))

#node initiation
ghrhr = mf.Node("ghrhr")
vipr1 = mf.Node("vipr1")
vipr2 = mf.Node("vipr2")
n13 = mf.Node("n13", left = ghrhr, right = vipr1)
crfr2 = mf.Node("crfr2")
crfr1 = mf.Node("crfr1")
calcr = mf.Node("calcr")
calrl = mf.Node("calrl")
gipr = mf.Node("gipr")
gcgr = mf.Node("gcgr")
pacr = mf.Node("pacr")
n10 = mf.Node("n10", left = vipr2, right = n13)
pth1r = mf.Node("pth1r")
pth2r = mf.Node("pth2r")
n11 = mf.Node("n11", left = crfr2, right = crfr1)
n12 = mf.Node("n12", left = calcr, right = calrl)
glp1r = mf.Node("glp1r")
n6 = mf.Node("n6", left = gipr, right = gcgr)
sctr = mf.Node("sctr")
n7 = mf.Node("n7", left = pacr, right = n10)
n8 = mf.Node("n8", left = pth1r, right = pth2r)
n9 = mf.Node("n9", left = n11, right = n12)
glp2r = mf.Node("glp2r")
n3 = mf.Node("n3", left = glp1r, right = n6)
n4 = mf.Node("n4", left = sctr, right = n7)
n5 = mf.Node("n5", left = n8, right = n9)
n1 = mf.Node("n1", left = glp2r, right = n3)
n2 = mf.Node("n2", left = n4, right = n5)
n0 = mf.Node("n0", left = n1, right = n2)

#fill the leaves with fasta raw data
leaf_list = mf.get_leaf_list(n0) #returns the list of leaf node objects
for leaf in leaf_list:
    #read the raw data
    leaf_raw_file = leaf.name + "_raw_data.csv"
    leaf.rawdata = pd.read_csv(join(csvdir, leaf_raw_file))


node_list = mf.get_node_list(n0) #list of all nodes
node_names = [node.name for node in node_list] #list of all node names
node_dictionary = mf.to_dictionary(n0) #dictionary node_name:node_object
choices = ["NC", "LC", "RC", "C", "DI", "DC"] #list of all possible position statuses
#LC -> only conserved for left sister
#RC -> only conserved for right sister
#NC -> not conserved
#C  -> conserved
#DI -> conserved for both sisters but conservation difference is bigger than allowed
#DC -> conserved for both sisters but residues are different

#sidebar for options
st.sidebar.header("Options")
cons_slider = st.sidebar.slider("Conservation limit (%)", min_value=50, max_value=100, step=1, value=90,
help="Frequency threshold for a residue to be considered consensus (gaps are not considered)")
gap_slider = st.sidebar.slider("Gap tolerance (%)", min_value=0, max_value=50, step=1, value=50,
help="Positions that have gaps above given limit are directly considered non-conserved")
difference_slider = st.sidebar.slider("Difference tolerance (%)", min_value=0, max_value=50, step=1, value=5,
help="While considering consensus, if difference of conservation between sisters are above given threshold, positions is not considered consensus.")
node_selected = st.sidebar.multiselect("Select nodes to display", node_names, node_names)
attribute_selected = st.sidebar.multiselect("Select residue types to display", choices, choices[1:])
st.sidebar.markdown("""
- **LC:** only conserved for left sister
- **RC:** only conserved for right sister
- **NC:** not conserved
- **C:** conserved
- **DI:** conserved for both sisters but conservation difference is bigger than allowed
- **DC:** conserved for both sisters but residues are different""")
ancestral_dif_bool = st.sidebar.checkbox(label="Ancestral difference", value=True, 
help= "Display ancestral difference of a node (consensus positions of the current node that are different from each of its ancestor nodes")

#global variables
consensus_limit = cons_slider/100
gap_tolerance = gap_slider/100
difference_tolerance = difference_slider/100

#functions
@st.cache
def leaf_consensus(node: mf.Node) -> pd.DataFrame:
    """get consensus of leaf nodes"""
    node_raw_data = node.rawdata #get raw data
    node_consensus = mf.new_consensus(node_raw_data, consensus_limit, gap_tolerance) #get consensus
    return node_consensus

@st.cache
def inner_consensus(node: mf.Node) -> pd.DataFrame :
    """get consensus of inner nodes. If a position:
    1. conserved for left child but not right child it is labeled LC (left conserved)
    2. conserved for right child but not left child it is labeled RC (right conserved)
    3. not conserved for both children it is labeled NC (not conserved)
    4. conserved for both children, has the same residue, and conservation frequency
    difference is within given difference tolerance limit it is labeled C (conserved)
    5. conserved for both children, has the same residue, but conservation frequency
    difference is bigger than allowed it is labeled DI (difference intolerant)
    6. conserved for both children but has different residues it is labeled 
    DC (different conserved)"""
    left_sister = node.left #get left sister
    right_sister = node.right #get rigth sister
    left_raw_data = mf.get_leaf_rawdata(left_sister) #get raw data dict
    right_raw_data = mf.get_leaf_rawdata(right_sister)
    left_values = list(left_raw_data.values()) #dict to list
    right_values = list(right_raw_data.values())
    node_values = left_values + right_values #concat left and right raw data (main node whole data)
    left_data = pd.concat(left_values, axis=1) #concat raw datas
    right_data = pd.concat(right_values, axis=1)
    node_data = pd.concat(node_values, axis=1)
    node_cons = mf.new_consensus(node_data, consensus_limit, gap_tolerance) #get consensus
    left_cons = mf.new_consensus(left_data, consensus_limit, gap_tolerance)
    right_cons = mf.new_consensus(right_data, consensus_limit, gap_tolerance)
    result_df = pd.DataFrame()
    result_df["residue"] = node_cons["residue"] #residues and frequencies are from main node whole data
    result_df["frequency"] = node_cons["frequency"]
    left_cons = left_cons.add_suffix("_left") #add suffixes to distinguish
    right_cons = right_cons.add_suffix("_right")
    return_df = pd.concat([left_cons, right_cons, result_df], axis=1)
    conditions = [
        return_df["status_left"].eq("NC") & return_df["status_right"].eq("NC"),
        return_df["status_left"].eq("C") & return_df["status_right"].eq("NC"),
        return_df["status_left"].eq("NC") & return_df["status_right"].eq("C"),
        return_df["status_left"].eq("C") & return_df["status_right"].eq("C") & return_df["residue_left"].eq(return_df["residue_right"]) & (abs(return_df["frequency_left"] - return_df["frequency_right"]) <= difference_tolerance),
        return_df["status_left"].eq("C") & return_df["status_right"].eq("C") & return_df["residue_left"].eq(return_df["residue_right"]) & (abs(return_df["frequency_left"] - return_df["frequency_right"]) > difference_tolerance),
        return_df["status_left"].eq("C") & return_df["status_right"].eq("C") & return_df["residue_left"].ne(return_df["residue_right"])
    ]
    return_df["status"] = np.select(conditions, choices)    
    return return_df

@st.cache
def new_df_diff(df1: pd.DataFrame, df2: pd.DataFrame, index_list: List) -> List:
    """return indexes of rows where residues are different. Only rows that are conserved
    (indicated by index list) are considered"""
    sub_df1 = df1.loc[index_list] #subset conserved positions
    sub_df2 = df2.loc[index_list]
    #find the rows where residues are different
    sub_df1["diff"] = np.where(sub_df1["residue"] != sub_df2["residue"], "True", "False")
    #list the indexes
    diff_idxs = list(sub_df1.loc[sub_df1["diff"] == "True"].index)
    return diff_idxs

#fill the nodes with consensus data, conserved position lists
for node in node_list:
    if mf.is_leaf(node):
        node_data = leaf_consensus(node)
        node.data = node_data
        node.metadata = {"internal_check" : list(node_data.loc[node_data["status"] == "C"].index)}
    else:
        node_data = inner_consensus(node)
        node.data = node_data
        node.metadata = {"internal_check" : list(node_data.loc[node_data["status"] == "C"].index)}

#define ancestral trails
#ancestral trail of each receptor
glp2rlist = [n0, n1, glp2r]
glp1rlist = [n0, n1, n3, glp1r]
giprlist = [n0, n1, n3, n6, gipr]
gcgrlist = [n0, n1, n3, n6, gcgr]
sctrlist = [n0, n2, n4, sctr]
pacrlist = [n0, n2, n4, n7, pacr]
vipr2list = [n0, n2, n4, n7, n10, vipr2]
ghrhrlist = [n0, n2, n4, n7, n10, n13, ghrhr]
vipr1list = [n0, n2, n4, n7, n10, n13, vipr1]
pth1rlist = [n0, n2, n5, n8, pth1r]
pth2rlist = [n0, n2, n5, n8, pth2r]
crfr2list = [n0, n2, n5, n9, n11, crfr2]
crfr1list = [n0, n2, n5, n9, n11, crfr1]
calcrlist = [n0, n2, n5, n9, n12, calcr]
calrllist = [n0, n2, n5, n9, n12, calrl]

ancestrydict = {"glp2r" : glp2rlist, "glp1r" : glp1rlist, "gipr" : giprlist, "gcgr" : gcgrlist, "sctr" : sctrlist, "pacr" : pacrlist, "vipr2" : vipr2list, "ghrhr" : ghrhrlist, 
"vipr1" : vipr1list, "pth1r" : pth1rlist, "pth2r" : pth2rlist, "crfr2" : crfr2list, "crfr1" : crfr1list, "calcr" : calcrlist, "calrl" : calrllist}

#get the total ancestral difference of each node
#total ancestral difference: consensus positions of the current node that are different from each of its ancestors 
for key in ancestrydict: #we will traverse from root to the node (follow the ancestral trail)
    ancestrylist = ancestrydict[key] #get the ancestral trail
    seen = [] #nodes that are alraady seen while traversing the trail
    for i in range(len(ancestrylist)): #iterate over the trail
        current_node = ancestrylist[i] #current node that we want to compare with its ancestors
        if seen: #if there is any node in the seen list
            current_consensus = current_node.metadata["internal_check"] #get the consensus positions of the current node
            differences = [] #node vs ancestor list
            for j in seen: #iterate over seen nodes (ancestors)
                diff = new_df_diff(current_node.data, j.data, current_consensus) #get the node vs ancestor difference
                differences.append(diff) #add them to the differences list
            total_ancestral_diff = list(set.intersection(*[set(x) for x in differences])) #get the intersections of all differences
            total_ancestral_diff.sort()
            current_node.metadata["total_ancestral_diff"] = total_ancestral_diff #assign the total ancestral difference to the metadata            
        seen.append(current_node)

#display the selected nodes
for node in node_selected:
    node_obj = node_dictionary[node] #get the node object
    if node_obj.name in node_selected:
        st.write(node_obj.name.upper()) #write the name of the node
        df = node_obj.data
        st.dataframe(df.loc[df["status"].isin(attribute_selected)])
        metadata = node_obj.metadata
        if ancestral_dif_bool:
            key = "total_ancestral_diff"
            if key in metadata:
                adp = metadata["total_ancestral_diff"]
                st.write("Ancestral difference positions")
                st.dataframe(df.loc[adp])
        st.markdown("""---""")


st.write("My program took", time.time() - start_time, "seconds to run")