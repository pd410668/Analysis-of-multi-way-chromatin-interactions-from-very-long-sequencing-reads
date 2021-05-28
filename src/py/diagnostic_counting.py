#!/usr/bin/env python

import pandas as pd

# experiment on "hs_k562_I_1"

df = pd.read_csv("hs_k562_I_1.tsv", sep='\t')
experiment_name = "hs_k562_I_1"

print("with condition: compare only between one chromosome")
print(df["RvsR"].value_counts())

df_RvsR = {x: y for x, y in df.groupby("RvsR")}

R1vsR1 = df_RvsR["R1 vs R1"]
R1vsR2 = df_RvsR["R1 vs R2"]
# R2vsR1 = df_RvsR["R2 vs R1"] # lack in hs_k562_I_1
R2vsR2 = df_RvsR["R2 vs R2"]

print("R1 vs R1: max and min abs value of two positions")
print(R1vsR1["abs_pos"].max())
print(R1vsR1["abs_pos"].min())

print("R1 vs R2: max and min abs value of two positions")
print(R1vsR2["abs_pos"].max())
print(R1vsR2["abs_pos"].min())

# print("R2 vs R1: max and min abs value of two positions")
# print(R2vsR1["abs_pos"].max())
# print(R2vsR1["abs_pos"].min())

print("R2 vs R2: max and min abs value of two positions")
print(R2vsR2["abs_pos"].max())
print(R2vsR2["abs_pos"].min())

print(df["chr_R1"].value_counts())



"""
Output for "hs_k562_I_1"

without conditions:

R1 vs R2    1254288
R1 vs R1       1144
R2 vs R2        358

R1 vs R1: max and min abs value of two positions
243838403
13

R1 vs R2: max and min abs value of two positions
243838440
0

R2 vs R2: max and min abs value of two positions
219367672
0

# example atypical chromosomes found in "hs_k562_I_1"

chrM                       4794
chrUn_gl000232              495
chrUn_gl000220              127
chrUn_gl000234               77
chrUn_gl000231               46
chrUn_gl000224               37
chrUn_gl000219               20
chrUn_gl000240                7
chr4_gl000194_random          3
chr7_gl000195_random          3
chr4_gl000193_random          2
chr1_gl000192_random          1
chrUn_gl000241                1
chr19_gl000208_random         1
chrUn_gl000238                1
chrUn_gl000216                1
chrUn_gl000214                1
chr17_gl000205_random         1


with condition: 
# compare only between one chromosome
# absolute value of their position does not exceed 1,000 bp
# after removed atypical chromosomes 

R1 vs R2    690899
R1 vs R1       139
R2 vs R2        76

R1 vs R1: max and min abs value of two positions
926
13

R1 vs R2: max and min abs value of two positions
1000
0

R2 vs R2: max and min abs value of two positions
926
0

# example atypical chromosomes found in "hs_k562_I_1"

chrM                       4442
chrUn_gl000232              428
chrUn_gl000220              103
chrUn_gl000224               33
chrUn_gl000234               26
chrUn_gl000231                8
chrUn_gl000240                5
chr7_gl000195_random          3
chr4_gl000194_random          2
chr4_gl000193_random          1
chr19_gl000208_random         1

"""
