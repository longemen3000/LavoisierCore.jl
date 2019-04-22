# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:53:29 2019

@author: richar56
"""

import csv
import numpy as np
import pandas as pd

import os
from os import listdir
from os.path import isfile, join
import re
pattern = re.compile(r'\s+')


mypath = os.path.dirname(os.path.realpath("__file__"))
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

onlyfiles.remove("combine_the_data.py")

filename1 = "Yaws_ocr_cp_gas.csv"

def separate_first(alist):
    first = alist[0]
    firstchar = 0
    for char in first:
        if char.isdigit():
            firstchar += 1
        else:
            break
    realfirst = first[:firstchar]
    realsecond = first[firstchar:]
    blist = [realfirst, realsecond]
    for i in range(1, len(alist)):
        blist.append(alist[i])
    return blist

def combine_name(alist):
    if len(alist) < 4:
        return []
    name1 = alist[2]
    name2 = alist[3]
    
    if len(name1) == 1 and name1[0].isdigit():
        #we need to see how many other digits there might be
        if len(name2) == 1 and name2[0].isdigit():
            
    
    else:
        if all([ch.isalpha() for ch in name2]) and len(name2) > 1:
            #then name2 is a real part of the name, hopefully
            name = name1 + " " + name2
            del alist[2]
            del alist[2]
            alist.insert(2, name)
            return alist
        else: 
            return alist

def preprocessing(afilename):
    
    filename2 = afilename[:-4] + "_out.csv"   
    output = open(filename2, "w")
    
    with open(afilename) as file1:
        headflag = True
        mywriter = csv.writer(output, delimiter=",")
        for line in file1.readlines():
            global headflag
            sentence = re.sub(pattern, ',', line)
            
            s2 =sentence.replace(",.,",",")
            s3 = s2.replace("|","")
            s4 = s3.replace("[","")
            s5 = s4.replace(",_,",",")
            s6 = s5.replace(",,",",")
            s7 = s6.replace(",=,",",")
            s8 = s7.replace("\"","")
            
            preprocessed = s8.split(",")
            try:
                i = int(preprocessed[0])
            except ValueError:
                preprocessed = separate_first(preprocessed)
                
            for i, item in enumerate(preprocessed):
                if item == '':
                    preprocessed.remove('')
            preprocessed = combine_name(preprocessed)
            ischem = all([ch.isdigit() for ch in preprocessed[0]])
            if ischem:
                headflag=False
            if preprocessed != [] and ischem:
                mywriter.writerow(preprocessed)
            
    output.close()
    
    
for file in onlyfiles:
    preprocessing(file)