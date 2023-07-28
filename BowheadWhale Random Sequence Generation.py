#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:39:27 2022

@author: brandonslockey

This Software will only count perfect inverted repeats.

"""

import random
import pandas as pd
from striprtf.striprtf import rtf_to_text


def generate_random_sequence(length,name,n):
    """
    Parameters
    ----------
    length : integer
        Length of Sequence.
    name : string
        name of textfile.
    n : integer
        how many text files.
    """
    
        
    for x in range(n):
        my_list = ["a","g","c","t"]
        random_mtDNA_Sequence = random.choices(my_list, weights=[32.7,13.16,28.07,26.06], k=length)
    
        with open(name + str(x) + ".txt", 'w') as output:
            for row in random_mtDNA_Sequence:
                output.write(str(row) + '\n')
            
        with open(name + str(x) + ".txt", 'r') as file:
            data = rtf_to_text(file.read()) 
            data.replace(" ","")

generate_random_sequence(16389, "RandomSequence1_WeightedBowhead_Run",10)
# Count Normal Repeats 16389


# Open FASTA File


# Strip RTF since RTF contains lots of unnecessary properties text
with open("Bowheadwhale 1 FASTA.rtf", 'r') as file:
    data = rtf_to_text(file.read())
    
length = len(data)
Length_bin = 9
data = data
dictionary = {}

list1 = list(range(9,17))
for Length_bin in list1:
    for index in range(len(data[0:length])): #stop 20bp before the end... 
        string = data[index:index+Length_bin]
        if len(string) >= 9: # avoid matching truncation at the end
            count = data.count(string)
            list1 = [count, Length_bin]
            if count >1:
                dictionary[string] = list1
    Length_bin = Length_bin +1
 
    
df = pd.DataFrame.from_dict(dictionary)

df = (df.T)


df.to_excel('dict1.xlsx')

# Reverse Complementary - count inverted repeats

# with open("Bowheadwhale 1 FASTA.rtf", 'r') as file:
#     data = rtf_to_text(file.read())
    
# with open("Human mtDNA FASTA.rtf", 'r') as file:
#     data = rtf_to_text(file.read())   


# Test Sequence
# Test_Sequence = "atcgtacgataaa"
# data = Test_Sequence
    
dictionary_raw_repeats = {}
dictionary_counts = {}
n=10
run_file = "RandomSequence1_WeightedBowhead_Run"
for run_index in list(range(0,n)):
    
    
    with open(run_file + str(run_index) + ".txt", 'r') as file:
        data = rtf_to_text(file.read()) 
    
# with open("RandomSequence1_WeightedBowhead2.txt", 'r') as file:
#     data = rtf_to_text(file.read()) 

    reverse_complementary = data[::-1]
    rc_list = []
    rc_list[:0] = reverse_complementary
    
    
    for index,x in enumerate(rc_list):
        if x == "a":
            rc_list[index] = "t"
        if x == "t":
            rc_list[index] = "a"
        if x == "g":
            rc_list[index] = "c"
        if x == "c":
            rc_list[index] = "g"
    
    rc_string = ''.join(rc_list)
    
    length = len(data)
    min_Length = 4
    Max_Length = 25
    OG_data = data

    
    length_list_choice = list(range(min_Length,Max_Length))
    
    for Length_bin in length_list_choice:
        index_tracker_list = []
        data = OG_data + OG_data[0:Length_bin-1] #mtDNA is a circular molecule. this allows repeats at the circularised end to be counted.
        count_Variable = 0
        for index in range(len(data[0:length])): 
            string = data[index:index+Length_bin] 
            if len(string) >= Length_bin: # avoid matching truncation at the end
                count = rc_string.count(string)
                list1 = [count, Length_bin, run_index,string, index]
                if count >0: #FOR INVERTED REPEATS, greater than 0, NOT 1 which works for repeats.
                    dictionary_raw_repeats[str(run_index)+"_"+str(Length_bin)+"_"+str(index)] = list1
                    count_Variable = count_Variable + 1
                    index_tracker_list.append(index)
        how_many_1k = sum(map(lambda x : x<1000, index_tracker_list))
        how_many_2k = sum(map(lambda x : x>=1000 and x<2000, index_tracker_list))
        how_many_3k = sum(map(lambda x : x>=2000 and x<3000, index_tracker_list))
        how_many_4k = sum(map(lambda x : x>=3000 and x<4000, index_tracker_list))
        how_many_5k = sum(map(lambda x : x>=4000 and x<5000, index_tracker_list))
        how_many_6k = sum(map(lambda x : x>=5000 and x<6000, index_tracker_list))
        how_many_7k = sum(map(lambda x : x>=6000 and x<7000, index_tracker_list))
        how_many_8k = sum(map(lambda x : x>=7000 and x<8000, index_tracker_list))
        how_many_9k = sum(map(lambda x : x>=8000 and x<9000, index_tracker_list))
        how_many_10k = sum(map(lambda x : x>=9000 and x<10000, index_tracker_list))
        how_many_11k = sum(map(lambda x : x>=10000 and x<11000, index_tracker_list))
        how_many_12k = sum(map(lambda x : x>=11000 and x<12000, index_tracker_list))
        how_many_13k = sum(map(lambda x : x>=12000 and x<13000, index_tracker_list))
        how_many_14k = sum(map(lambda x : x>=13000 and x<14000, index_tracker_list))
        how_many_15k = sum(map(lambda x : x>=14000 and x<15000, index_tracker_list))
        how_many_16k = sum(map(lambda x : x>=15000 and x<16000, index_tracker_list))
        how_many_17k = sum(map(lambda x : x>=16000 and x<17000, index_tracker_list))
        dictionary_counts[str(run_index)+"_"+str(Length_bin)] = [count_Variable, Length_bin, run_index, how_many_1k, how_many_2k, how_many_3k, how_many_4k, how_many_5k, how_many_6k, how_many_7k, how_many_8k, how_many_9k, how_many_10k, how_many_11k, how_many_12k, how_many_13k, how_many_14k, how_many_15k, how_many_16k, how_many_17k]
        Length_bin = Length_bin +1
    print("We just finished Run Index: ", run_index)
     
        
df = pd.DataFrame.from_dict(dictionary_raw_repeats)

df = (df.T)


df.to_excel(run_file + 'dict_Inverted_repeats.xlsx', header=["Repeat_No","Repeat Length","Run_No","Sequence", "Start_Position"])   


df2 = pd.DataFrame.from_dict(dictionary_counts)

df2 = (df2.T)


df2.to_excel(run_file + 'dict_Inverted_Repeats_Tally.xlsx',header=["Total Repeat_No","Repeat Length", "Run_No", "how_many_1k", "how_many_2k", "how_many_3k", "how_many_4k", "how_many_5k", "how_many_6k", "how_many_7k", "how_many_8k", "how_many_9k", "how_many_10k", "how_many_11k", "how_many_12k", "how_many_13k", "how_many_14k", "how_many_15k", "how_many_16k", "how_many_17k"])   


####################################
######### TOM CODE !!!! #########
####################################
import plotly.graph_objs as go

def get_theta_from_human_mtdna_position(pos):
    # human mtDNA length is hardcoded so it can be used in map()
    if pos < 1:
        pos = 16568 + pos
    return ((pos - 1) / 16568) * 360 % 360


def make_poslist(starts, length):
    """Returns a flat list of points covered by repeats given their start and length"""
    x = []
    for i in starts:
        x += list(range(i, i+length))
    return x

 # change these to the real ones
human_repeat_positions = (314, 965, 3576, 15545, 16193)
whale_repeat_positions = tuple(start_Positions_List)

# create lists of positions at which there is a repeat
human_poslist = make_poslist(human_repeat_positions, 9)
whale_poslist = make_poslist(whale_repeat_positions, 9)

# create sequences which will be used to label the angular axis
ticks = list(range(0, 16568, 1000))
tickvals = list(map(get_theta_from_human_mtdna_position, ticks))
ticktext = [f"{n//1000}k" for n in ticks]

fig = go.Figure()

# add the line representing the DNA (nb. you can't draw outside the axis line hence not using that)
fig.add_trace(
    go.Scatterpolar(
        r=[1] * 361,
        theta=list(range(361)),
        mode='lines',
        line_color='black',
        hoverinfo='skip',
        showlegend=False
))

# add human repeats
fig.add_trace(
    go.Scatterpolar(
        r=[0.985] * len(human_poslist),
        theta=list(map(get_theta_from_human_mtdna_position, human_poslist)),
        mode='markers',
        name='Human repeats',
        hovertext=[f"Human\tPos: {p}" for p in human_poslist],
        hoverinfo='text'
))

# add whale repeats
fig.add_trace(
    go.Scatterpolar(
        r=[1.015] * len(whale_poslist),
        theta=list(map(get_theta_from_human_mtdna_position, whale_poslist)),
        mode='markers',
        name='Whale repeats',
        hovertext=[f"Whale\tPos: {p}" for p in whale_poslist],
        hoverinfo='text'
))

# format the axes
fig.update_layout(
    template=None,
    height=800,
    polar=dict(
        angularaxis=dict(
            direction='clockwise',
            showgrid=False,
            showline=False,
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext,
            ticks='',
        ),
        radialaxis=dict(
            showgrid=False,
            showline=False,
            showticklabels=False,
            range=(0, 1.02),
            ticks='',
        )
    )
)

fig.show()


