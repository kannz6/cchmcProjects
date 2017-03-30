#! /bin/bash

echo "wall = \"18:00\"" > fromLoniTmpSettings.txt
echo "memory = \"18000\"" >> fromLoniTmpSettings.txt
echo "" >> fromLoniTmpSettings.txt
echo "bwa = \"0.7.15\"" >> fromLoniTmpSettings.txt
echo "samtools = \"1.3\"" >> fromLoniTmpSettings.txt
echo "bcftools = \"1.3\"" >> fromLoniTmpSettings.txt
echo "plink = \"1.90b\"" >> fromLoniTmpSettings.txt
echo "king = \"2.0\"" >> fromLoniTmpSettings.txt
echo "" >> fromLoniTmpSettings.txt

# echo "wall = \"$3\"" > fromLoniTmpSettings.txt
# echo "memory = \"$4\"" >> fromLoniTmpSettings.txt
# echo "" >> fromLoniTmpSettings.txt
# echo "bwa = \"$5\"" >> fromLoniTmpSettings.txt
# echo "samtools = \"$6\"" >> fromLoniTmpSettings.txt
# echo "bcftools = \"$7\"" >> fromLoniTmpSettings.txt
# echo "plink = \"$8\"" >> fromLoniTmpSettings.txt
# echo "king = \"$9\"" >> fromLoniTmpSettings.txt
# echo "" >> fromLoniTmpSettings.txt

python dynamicSettings.py > temp.txt
# cp temp.txt versionSettings.py

