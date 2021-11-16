# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 17:53:50 2021

@author: Maria
"""

#Filter reads by length (min and max) and trim bases at extreme 5' and 3'
#Parameters: length_min; length_max; quality; trim5; trim3
#Open file and read one sequence
 

qual_sanger = {'!': 0, '#': 2, '"': 1, '%': 4, '$': 3,       # Table of quality values for illumina 1.8+ and sanger
               "'": 6, '&': 5, ')': 8, '(': 7, '+': 10,
               '*': 9, '-': 12, ',': 11, '/': 14, '.': 13,
               '1': 16, '0': 15, '3': 18, '2': 17, '5': 20,
               '4': 19, '7': 22, '6': 21, '9': 24, '8': 23,
               ';': 26, ':': 25, '=': 28, '<': 27, '?': 30,
               '>': 29, 'A': 32, '@': 31, 'C': 34, 'B': 33,
               'E': 36, 'D': 35, 'G': 38, 'F': 37, 'I': 40,
               'H': 39, 'J': 41, 'K':42, 'L':43, 'M':44,
               'N':45, 'O':46,}

def trans_sanger(value):                                    #Translate a quality symbol to value
    return qual_sanger[value]

def translate_qual_read(line4):                             #Translate a read to values serie
    return [trans_sanger(value) for value in line4]

def qual_read(line4):                                       # Average value of quality read (values serie)
    val_ave = 0.0
    values = translate_qual_read(line4)
    val_ave = float(sum(values))/len(values)
    return val_ave


count_badS = 0
count_badL = 0
bad_qual = 0
count_seq = 0
result = []
read_len = []
"""

""" 
with open(snakemake.output[0],'w+') as fileOUT1:
    fileOUT1.readline()

with open (snakemake.input[0]) as file1:
    line1 = file1.readline().replace('\n','')
    line2 = file1.readline().replace('\n','')
    line3 = file1.readline().replace('\n','')
    line4 = file1.readline().replace('\n','')
    assert len(line2) == len(line4), 'different length line2 and line4'
    while line1 !='':
        #Filter by length min and max [user parameters]
        if len(line2) < int(snakemake.config["Lmin"]):
            count_badS += 1
        elif len(line2) > int(snakemake.config["Lmax"]):
            count_badL += 1
        else:
    #Trim bases at 5' extreme [user parameter]

            if snakemake.config["t3"] == 'auto':
                assert snakemake.config["Q"] != None, "You have to define a min quality value (-Q)"
                while len(line2) > int(snakemake.config["Lmin"])-1 and qual_read(line4) < float(snakemake.config["Q"]):
                    line2 = line2[1:]
                    line4 = line4[1:]
            elif snakemake.config["t5"] == 'N':
                while len(line2) > int(snakemake.config["Lmin"])-1 and line2[0:1] == 'N':
                    line2 = line2[1:]
                    line4 = line4[1:]
            elif int(snakemake.config["t5"]) > 0:
                line2 = line2[int(snakemake.config["t5"]):]
                line4 = line4[int(snakemake.config["t5"]):]
            elif snakemake.config["t5"] == None:
                pass

#trim bases at 3' extreme [user parameter]


            if snakemake.config["t3"] == 'auto':
                assert snakemake.config["Q"] != None, "You have to define a min quality value (-Q)"
                while len(line2) > int(snakemake.config["Lmin"])-1 and qual_read(line4) < float(snakemake.config["Q"]):
                    line2 = line2[:-1]
                    line4 = line4[:-1]
            elif snakemake.config["t3"] == 'N':
                while len(line2) > int(snakemake.config["Lmin"])-1 and line2[-1:] == 'N':
                    line2 = line2[:-1]
                    line4 = line4[:-1]
            elif int(snakemake.config["t3"]) > 0:
                line2 = line2[:-int(snakemake.config["t3"])]
                line4 = line4[:-int(snakemake.config["t3"])]
            elif snakemake.config["t3"] == None:
                pass


#count read
            if len(line2) < int(snakemake.config["Lmin"]):
                bad_qual += 1
            else:
                count_seq += 1
                q_seq = qual_read(line4)
                result.append(str(q_seq)[0:5])
                read_len.append(len(line2))
                with open (snakemake.output[0],'a') as fileOUT1:
                    fileOUT1.writelines(line1+'\n')
                    fileOUT1.writelines(line2+'\n')
                    fileOUT1.writelines(line3+'\n')
                    fileOUT1.writelines(line4+'\n') 
        line1 = file1.readline().replace('\n','')
        line2 = file1.readline().replace('\n','')
        line3 = file1.readline().replace('\n','')
        line4 = file1.readline().replace('\n','')
        assert len(line2) == len(line4), 'different length line2 and line4'