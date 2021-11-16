# Command line: file1, file2, -N, Q= value, -R [file report], -L min length
# 4-06-19 Version .4 Modified to filter by len

#############################################
with open(snakemake.output[0],'w+') as fileOUT1:
            fileOUT1.readline()

with open(snakemake.output[1],'w+') as fileOUT2:
            fileOUT2.readline()


qual_illumina = {     # Table of quality values for illumina 1.8+ and IonTorrent
'A': 32, '@': 31, 'C': 34, 'B': 33,
'E': 36, 'D': 35, 'G': 38, 'F': 37,
'I': 40, 'H': 39, 'J': 41, '!': 0,
'#': 2, '"': 1, '%': 4, '$': 3,
"'": 6, '&': 5, ')': 8, '(': 7,
'+': 10, '*': 9, '-': 12, ',': 11,
'/': 14, '.': 13, '1': 16, '0': 15,
'3': 18, '2': 17, '5': 20, '4': 19,
'7': 22, '6': 21, '9': 24, '8': 23,
';': 26, ':': 25, '=': 28, '<': 27,
'?': 30, '>': 29,'K': 42, 'L': 43,
'M':44, 'N':45, 'O': 46,}

def translate_qual_value(value):    #Translate a quality symbol to value
    return qual_illumina[value]

def translate_qual_read(line4):     #Translate a read to values serie
    return [translate_qual_value(value) for value in line4]

def qual_read(line4):     # Average value of quality read (values serie)
    val_ave = 0.0
    values = translate_qual_read(line4)
    val_ave = float(sum(values))/len(values)
    return val_ave

list_report = []
valid_bases = set('ACGTacgt')
count_N = 0
count_Q = 0
count_seq = 0

def validate_bases(line):
    return set(line) <= valid_bases

with open(snakemake.input[0]) as file1, open(snakemake.input[1]) as file2:
    line1_1 = file1.readline().replace('\n','')
    line1_2 = file1.readline().replace('\n','')
    line1_3 = file1.readline().replace('\n','')
    line1_4 = file1.readline().replace('\n','')
    line2_1 = file2.readline().replace('\n','')
    line2_2 = file2.readline().replace('\n','')
    line2_3 = file2.readline().replace('\n','')
    line2_4 = file2.readline().replace('\n','')
    while line1_1 and line2_1 != '':                 # Filter for N bases, Q value and length
        assert line1_1.split(" ")[0] == line2_1.split(" ")[0], "Not paired-end reads"
        count_seq += 1
        if validate_bases(line1_2) is False and snakemake.config["N"] == True:
            list_report.append(line1_1)#[0:41])
            count_N += 1
        elif qual_read(line1_4) < int(snakemake.config["Q"]) or len(line1_4) < int(snakemake.config["min_read_length"]):
            list_report.append(line1_1)#[0:41])
            count_Q += 1
        elif validate_bases(line2_2) is False and snakemake.config["N"] == True:
            list_report.append(line2_1)#[0:41])
            count_N += 1
        elif qual_read(line2_4) < int(snakemake.config["Q"]) or len(line2_4) < int(snakemake.config["min_read_length"]):
            list_report.append(line2_1)#[0:41])
            count_Q += 1
        else:
            with open(snakemake.output[0],'a') as fileOUT1,\
                 open(snakemake.output[1],'a') as fileOUT2:
                fileOUT1.write(line1_1+'\n')
                fileOUT1.write(line1_2+'\n')
                fileOUT1.write(line1_3+'\n')
                fileOUT1.write(line1_4+'\n')
                fileOUT2.write(line2_1+'\n')
                fileOUT2.write(line2_2+'\n')
                fileOUT2.write(line2_3+'\n')
                fileOUT2.write(line2_4+'\n')
        line1_1 = file1.readline().replace('\n','')
        line1_2 = file1.readline().replace('\n','')
        line1_3 = file1.readline().replace('\n','')
        line1_4 = file1.readline().replace('\n','')
        line2_1 = file2.readline().replace('\n','')
        line2_2 = file2.readline().replace('\n','')
        line2_3 = file2.readline().replace('\n','')
        line2_4 = file2.readline().replace('\n','')



#######################################################################

if snakemake.config["N"] == False:
    count_N = 'N/A'