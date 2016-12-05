#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:52:37 2016

@author: Jiajia
"""
import pandas as pd
from bs4 import BeautifulSoup

# open html
soup = BeautifulSoup(open('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_fi_c20_std/cluster_1_AgriGO.html'),"lxml")

# find table
table = soup.find("table")

# parse table
datasets = []
for row in table.findAll("tr"):
    cells = row.findAll("td")
    if len(cells)==7:
        GO_term = cells[0].find(text=True)
        Ontology = cells[1].find(text=True)
        Description = cells[2].find(text=True)
        Number_input_list = cells[3].find(text=True)
        Number_BG = cells[4].find(text=True)
        p_value = cells[5].find(text=True)
        fdr_value = cells[6].find(text=True)
        new_row = [GO_term,Ontology,Description,Number_input_list,Number_BG,p_value,fdr_value]
    datasets.append(new_row)

# transform to dataframe and filter
go_result = pd.DataFrame(datasets)

go_result.to_csv('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_fi_c20_std/cluster_1_AgriGO.csv',encoding='utf-8',header=None,index=None)

