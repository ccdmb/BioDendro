
#Author: P Moolhuijzen
#Date  : 26 February 2018
#Inputs: <input1>.mgf and <input1>.csv text files
#Ouputs: <msms_nonredundant_list>.txt 

import os
import sys
import argparse
from collections import defaultdict
import xlsxwriter
import csv
import operator
import pandas as pd


def get_record(lines):
    """ Make a dictionary key pepmass-retention-time
    and value of ions from the trigger data (MGF file)
    """

    mydic = {}
    val = []
    p = ''
    for fl in lines:
        if 'TITLE' in fl:
            ft, f, spt, sp, sct, sc = fl.split(' ')
            fls = f.split('\\')
            fna = fls[-1]
            fid, suf = fna.split('.')
        elif 'PEPMASS' in fl:
            pt, p = fl.split('=')
            m, a = p.split(' ')
        elif 'RTINSECONDS' in fl:
            rt, rtn = fl.split('=')
            key = fid + "_" + m + "_" + rtn
        elif fl[0].isdigit():
            ion, area = fl.split(' ')
            flag = 1
            val.extend([ion])
        elif 'END' in fl:
            if flag == 1:
                mydic[key] = val
                val = []
                flag = 0
            else:
                val = []
                flag = 0

    return mydic


def get_csv_record(lines2, rec):
    """ Collects all the csv list (real samples list) matches to the trigger
    data within a mass of 0.002 and retention time of 5 secs return a
    dictionary key real sample id and value of triggers (id and ions list) 
    """

    nwdic = {}
    for cl in lines2:
        if 'Components' not in cl:
           mass = cl.split('_')
           mass[3] = mass[3].strip('m/z') # real sample mass
           rto = mass[4].strip('RT')
           rto = float(rto) * 60
           rto = int(rto) # real sample retention time 
           umass = float(mass[3]) + 0.002
           lmass = float(mass[3]) - 0.002
           urt = rto + 5
           lrt = rto - 5

           comp = []
           for key in rec:
               mk = key.split('_')
               mk[3] = float(mk[3]) # mass of trigger data
               mk[4] = int(mk[4]) # rentention time of trigger data
               if (lmass < mk[3] < umass) and (lrt < mk[4] < urt):
                   strg = key, rec[key] 
                   comp.append(strg)

           if len(comp) > 0:
               nwdic[cl] = comp

    return nwdic 


def remove_redundancy(ndic, negloss):
    """ Selects the closest trigger mass to the real sample mass 
    Prints the best trigger id and ion list 
    """

    tmpdic ={}
    for cpkey in ndic: # for each real sample id
        km = cpkey.split("_")
        kms = float(km[3].strip('m/z')) # real sample mass
        #dms = abs(kms - 0.002)
        dms = 0.002
        krt = km[4].strip('RT')
        krt = float(krt) * 60 # convert retention time to seconds
        krt = int(krt)   # real sample retention time
        #drt = abs(krt-5)
        drt = 5
        best = ndic[cpkey] # list of trigger matches (id and ions list)
        kval = []
        if len(ndic[cpkey]) == 1: # if only one trigger match
            for ion in ndic[cpkey]:
                one = ion[1] # get ion list
                bestion = ion[0].strip() # get trigger id (QE..)
                tm = bestion.split("_")
                tm[3] = tm[3].strip('m/z')
                for ion1 in one: # foreach ion
                    if negloss:
                        ion1 = round(float(ion1) - float(tm[3]), 5) # get negative loss
                    else: #options
                        ion1 = float(ion1.strip())
                    fele = bestion, ion1
                    tmpdic[fele] = cpkey
        else:
            for ls in ndic[cpkey]: # for each trigger match
                y = ls[0] # get trigger id
                vm = y.split('_')
                vms = float(vm[3]) # get trigger mass
                vrt = int(vm[4]) # get trigger retention time

                # find the closest trigger to the sample
                if (abs(kms - vms) < dms) and (abs(krt - vrt) < drt): 
                    dms = abs(kms - vms) # new distance of trigger mass to sample
                    drt = abs(krt - vrt) # new distance of trigger retention time to sample
                    best = ls # set best trigger match

            multi = best[1] # get ion list
            for bion in multi: # for each ion 
                bion1 = float(bion.strip())
                bestlab = best[0].strip() # get trigger label
                btm = bestlab.split("_") # get trigger id (QE..)
                btm[3] = float(btm[3].strip('m/z'))
                if negloss:
                    bion1 = round(bion1 - btm[3],5) # get negative loss
                else:
                    bion1 = bion1

                fele1 = bestlab, bion1
                tmpdic[fele1] = cpkey

    flip = {}
    for k, v in tmpdic.items():
        if k not in flip:
            flip[k] = [v]
        else:
            flip[k].append(v)

    alist = []
    for key, value in flip.items():
        #print(key[0],"\t",key[1],sep='')
        alist.append((key[0], key[1]))

    # Return the table, sorted by mz
    # Columns are (sample, mz)

    #Now remove redundancy and print best trigger ion list

    table = pd.DataFrame(alist, columns=['sample', 'mz'])
    table.sort_values(by='mz', inplace=True)
    return table

