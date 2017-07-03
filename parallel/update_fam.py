#!/usr/bin/env python
#Andrew Rupert 4-6-17 
import fileinput
import re
######################
#bulk vcf (yale)
######################
# r_parse_base_id=re.compile(r'^(.*?)(-\d\d)?$')
# new_columns = []
path = "tmp_fam.txt"
# tmpFamFileWriter = open("tmp_fam.txt", "w+")
# for line in fileinput.input():
#     line=line.strip()
#     columns=line.split(' ')
#     new_columns=columns[:]
#     #print(columns)
#     #print(new_columns)
#     match = r_parse_base_id.match(columns[0])

#     family_id = match.group(1)
#     sex = '0'
#     mother_id = '0'
#     father_id = '0'
#     if match.group(2) == '-01':
#         sex = '2'
#     elif match.group(2) == '-02':
#         sex = '1'
#     elif match.group(2) is None:
#         mother_id = family_id + '-01'
#         father_id = family_id + '-02'
#     else:
#         raise ValueError('HEEELP')

#     new_columns[0] = family_id
#     new_columns[2] = mother_id
#     new_columns[3] = father_id
#     new_columns[4] = sex
#     tmpFamFileWriter.write(' '.join(new_columns))
#     tmpFamFileWriter.write("\n")
#     # print(' '.join(new_columns))#bulk yale vcf use > to write to file in
# tmpFamFileWriter.close()
######################

######################
#fastq to kin pipeline
######################
r_parse_base_id=re.compile(r'^(.*?)(-\d\d)?$')
new_columns = []
path = "tmp_fam.txt"
tmpFamFileWriter = open("tmp_fam.txt", "w+")
for i,line in enumerate(fileinput.input()):
    line=line.strip()
    columns=line.split(' ')
    new_columns=columns[:]
    # print(columns)
    #print(new_columns)
    _id = new_columns[1].split(".")
   # print(_id)
    _id = _id[1]
    match = r_parse_base_id.match(columns[0])

    family_id = match.group(1)
    sex = '0'
    mother_id = '0'
    father_id = '0'

    if(i == 0):
        new_columns[0] = family_id
        new_columns[1] = _id
        new_columns[2] = family_id + '-01'
        new_columns[3] = family_id + '-02'
        new_columns[4] = sex
    elif (i % 2) == 1:
        if "-01" in _id:
            sex = '2'
        elif "-02" in _id:
            sex = '1'
        new_columns[0] = family_id
        new_columns[1] = _id
        new_columns[4] = sex
    elif (i % 2 ) == 0:
        if "-01" in _id:
            sex = '2'
        elif "-02" in _id:
            sex = '1'
        new_columns[0] = family_id
        new_columns[1] = _id
        new_columns[4] = sex
    tmpFamFileWriter.write(' '.join(new_columns))
    tmpFamFileWriter.write("\n")
    #print(' '.join(new_columns))#bulk yale vcf use > to write to file in
tmpFamFileWriter.close()
######################
