#!/usr/bin/python3
#Author: P Moolhuijzen
#Date  : 26 February 2018
#Inputs: <input1>.mgf and <input1>.csv text files
#Ouputs: <msms_nonredundant_list>.txt 

import sys, argparse
from collections import defaultdict
import xlsxwriter, csv
import operator
import pandas as pd

#make a dictionary key pepmass-retention-time and value of ions from the trigger data (MGF file)
def get_record(lines):
   mydic = {}
   val = []
   p = ''
   for fl in lines:
       if 'TITLE' in fl:
           ft,f,spt,sp,sct,sc = fl.split(' ')
           fls = f.split('\\')
           fna = fls[- 1]
           fid,suf = fna.split('.')
       elif 'PEPMASS' in fl:
           pt,p = fl.split('=')
           m,a = p.split(' ')
       elif 'RTINSECONDS' in fl:
           rt,rtn = fl.split('=')
           key = fid + "_" + m + "_" + rtn
       elif fl[0].isdigit():
           ion,area = fl.split(' ')
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

# collects all the csv list (real samples list) matches to the trigger data within a mass of 0.002 and retention time of 5 secs
# return a dictionary key real sample id and value of triggers (id and ions list) 
def get_csv_record(lines2, rec):

   nwdic = {}
   for cl in lines2:
      if 'Components' in cl:
         pass
      else:
         mass = cl.split('_')
         mass[3] = mass[3].strip('m/z') # real sample mass
         rto = mass[4].strip('RT')
         rto = float(rto)*60
         rto = int(rto) # real sample retention time 
         umass = float(mass[3])+0.002
         lmass = float(mass[3])-0.002
         urt = rto+5
         lrt = rto-5
         comp = []
         for key in rec:
            mk = key.split('_')
            mk[3] = float(mk[3]) # mass of trigger data
            mk[4] = int(mk[4]) # rentention time of trigger data
            if lmass < mk[3] < umass and lrt < mk[4] < urt:
               strg =  key,rec[key] 
               comp.extend([strg])
         if comp:         
             nwdic[cl] = comp
   
   return nwdic 

# Selects the closest trigger mass to the real sample mass 
# Prints the best trigger id and ion list 
def remove_redundancy(ndic, negloss):
   tmpdic ={}
   for cpkey in ndic: # for each real sample id
       km = cpkey.split("_")
       kms = km[3].strip('m/z') # real sample mass
       #dms = abs(float(kms)-0.002)
       dms = float(0.002)
       krt = km[4].strip('RT')
       krt = float(krt)*60 # convert retention time to seconds
       krt = int(krt)   # real sample retention time
       #drt = abs(krt-5)
       drt = int(5)
       best = ndic[cpkey] # list of trigger matches (id and ions list)
       kval = []
       if len(ndic[cpkey]) == 1: # if only one trigger match
           for ion in ndic[cpkey]:
               one = ion[1] # get ion list
               bestion=ion[0].strip() # get trigger id (QE..)
               tm = bestion.split("_")
               tm[3] = tm[3].strip('m/z')
               for ion1 in one: # foreach ion
                   if negloss:
                      ion1 = round(float(ion1)-float(tm[3]),5) # get negative loss
                   else: #options
                      ion1 = ion1.strip()
                   ion1 = float(ion1)
                   fele = bestion,ion1
                   tmpdic[fele] = cpkey
       else:
           for ls in ndic[cpkey]: # for each trigger match
               y = ls[0] # get trigger id
               vm = y.split('_')
               vms = vm[3] # get trigger mass
               vrt = vm[4] # get trigger retention time
               if (abs(float(kms)-float(vms)) < float(dms)) and (abs(int(krt)-int(vrt)) < int(drt)): # find the closest trigger to the sample
                   dms = abs(float(kms)-float(vms)) # set the new distance of trigger mass to sample
                   drt = abs(int(krt)-int(vrt)) # set the new distance of trigger retention time to sample
                   best = ls # set best trigger match
           multi = best[1] # get ion list
           for bion in multi: # for each ion 
               bion1 = bion.strip()
               bestlab = best[0].strip() # get trigger label
               btm = bestlab.split("_") # get trigger id (QE..)
               btm[3] = btm[3].strip('m/z')
               if negloss:
                  bion1 = round(float(bion1)-float(btm[3]),5) # get negative loss
               else:
                  bion1 = float(bion1)
               fele1 = bestlab,bion1
               tmpdic[fele1] = cpkey
               
   flip = {}
   for k, v in tmpdic.items():
      if k not in flip:
         flip[k] = [v]
      else:
         flip[k].append(v)
   
   alist=[]
   for key, value in flip.items():
       #print(key[0],"\t",key[1],sep='')
       alist.append((key[0],key[1]))

   with open('out.csv', 'w', newline='') as csvfile:
       spamwriter = csv.writer(csvfile, delimiter='\t')
       spamwriter.writerow(['sample','mz']) # file header
       sorted_x = sorted(alist, key=lambda x: x[1])
       for k1, k2 in sorted_x:
          spamwriter.writerow([k1,k2])
  
def main(argv):
   inputfile1 = ''
   inputfile2 = ''

   parser = argparse.ArgumentParser()
   parser.add_argument("file1", help="Please use a MGF file (file1.mgf)", type=str)    
   parser.add_argument("file2", help="Please use a file of components listed (file2.txt)", type=str)

   parser.add_argument("-n","--negative", help="You are applying negative loss", action="store_true")
   args = parser.parse_args()

   print(args.file1, args.file2, args.negative)


   #Open the trigger data <file>.msg
   f_in = open(args.file1, 'r')
   lines = filter(None, (line.rstrip() for line in f_in)) 
   
   rec = get_record(lines)

   #Open the sample list <file>.csv
   f2_in = open(args.file2, 'r')
   lines2 = filter(None, (line2.rstrip() for line2 in f2_in))  
   ndic = get_csv_record(lines2, rec)

   negloss = args.negative
   #Now remove redundancy and print best trigger ion list
   remove_redundancy(ndic, negloss)

   # get an excel file too
   pd.read_csv("out.csv", delimiter="\t").to_excel("out.xlsx", index=False)

if __name__ == "__main__":
   main(sys.argv[1:])


