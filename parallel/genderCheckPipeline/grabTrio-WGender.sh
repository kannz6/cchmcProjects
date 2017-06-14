#! /bin/bash

cd /scratch/kannz6/loni/vcf_pipeline/genderCheckPipeline
mkdir $2

cp * $2
cd $2

echo "wall = \"$3\"" > fromLoniTmpSettings.txt
echo "memory = \"$4\"" >> fromLoniTmpSettings.txt
echo "" >> fromLoniTmpSettings.txt
echo "bwa = \"$5\"" >> fromLoniTmpSettings.txt
echo "samtools = \"$6\"" >> fromLoniTmpSettings.txt
echo "bcftools = \"$7\"" >> fromLoniTmpSettings.txt
echo "plink = \"$8\"" >> fromLoniTmpSettings.txt
echo "king = \"$9\"" >> fromLoniTmpSettings.txt
echo "" >> fromLoniTmpSettings.txt

python dynamicSettings.py > tempSettings.txt

echo "$1" > triosToVerify.txt

awk '{print $1}' triosToVerify.txt | tr ',' '\n' > temp.txt
cp temp.txt triosToVerify.txt
rm temp.txt

sh wrapTriosToPSQLWGender.awk
echo $(basename $2)  > out.txt

# jid=$(bjobs | grep  $(expr ".* $9" : '.* \(.*\)') | awk '{print $1}')
jid=$(bjobs | grep  $(expr ".* 1.4" : '.* \(.*\)') | awk '{print $1}' | tac | egrep -m 1 .)
echo "$jid" > loniJobId.txt
echo "Last parameter to be used to get loni job-id: $9" >> out.txt

#sleep 10
#childJobsStart1=$(bjobs | grep $(date | awk '{print substr($4, 0, 5)}') | awk '{print $1}')
#childJobsStart2=$(bjobs | grep $(date | awk '{print substr($4, 0, 3) substr($4,4)+1}') | awk '{print $1}')
#childJobsStart1=$(date | awk '{print substr($4, 0, 5)}')
#childJobsStart2=$(date | awk '{print substr($4, 0, 3) substr($4,4)+1}')

python mainWGender.py
# last=$(cat query_output.dat | awk '{if($1 ~ /-02/){print substr($1,0,7)}}' > last.txt && echo $(tac last.txt | egrep -m 1 .)"_output/"$(tac last.txt | egrep -m 1 .)".king.kin")
last=$(echo $(tac $2/triosStructure.txt | egrep -m 1 .)"_output/"$(tac $2/triosStructure.txt | egrep -m 1 .)".king.kin")
echo "$last" >> out.txt
while ( ! ( test -e $last ) ); do sleep 180; done;
