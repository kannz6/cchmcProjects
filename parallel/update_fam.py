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
#fasq to kin pipeline
######################
with open(path,'r') as f:
    lines = f.readlines()
    #begin 4-10-2017 
    blind_id = lines[0].split()
    fam_id = blind_id[0]
    blind_id[1] = blind_id[0]

    #mother id fix
    blind_id_m = lines[1].split()
    blind_id_m[1] = fam_id + "-01"
    blind_id_m[0] = fam_id                 # Replace column one with the file name as well.
    blind_id_m[4] = "2"
    lines[1] = ' '.join(blind_id_m)
    blind_id[2] = blind_id_m[1]

    #father id fix
    blind_id_f = lines[2].split()
    blind_id_f[1] = fam_id + "-02"
    blind_id_f[0] = fam_id                # Replace column one with the file name as well.
    blind_id_f[4] = "1"
    lines[2] = ' '.join(blind_id_f)
    blind_id[3] = blind_id_f[1]

    lines[0] = ' '.join(blind_id)
    # end changes 4-10-17
with open(path,'w') as f:
    f.write('\n'.join(lines)+'\n')
f.close()