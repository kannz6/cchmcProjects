#! /bin/bash

# jobid=$(bjobs | grep $(date | awk '{print substr($4, 0, 5)}') | awk '{print $1}')
cd /scratch/kannz6/loni/vcf_pipeline/pipeline
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

sh wrapTriosToPSQL.awk
echo $(basename $2) > out.txt
dir="$(expr $(basename $2) : '.*-[0-9]\(.*\)')"
# jid="$(expr ".* $9" : '.* \(.*\)')"#works
# echo "$(expr "1.90b 2.0" : '.* \(.*\)')"#####todo use this to get job id of loni job
# test -z "$( bjobs <jobid> )" && echo "empty" || echo "job running..."
#while ( ! ( test -z <jobId> ) ); do sleep 3; done;
#while ( ( test -n 715969 ) ); do sleep 3; ; test -z "$( bjobs 715969 )" && echo "empty" > testJobID.txt || echo "job running..." > runningJob.txt done;
#echo "dir try regex: $dir" >> out.txt
# echo "jobid w/ regex: $jobid" >> out.txt
# jid=$(bjobs | grep  $(expr ".* $9" : '.* \(.*\)') | awk '{print $1}')
jid=$(bjobs | grep  $(expr ".* 1.4" : '.* \(.*\)') | awk '{print $1}' | tac | egrep -m 1 .)
echo "$jid" > loniJobId.txt
# echo "$(tac loniJobId.txt | egrep -m 1 .)" > loniJobId.txt
echo "Last parameter to be used to get loni job-id: $9" >> out.txt

python main.py

last=$(cat query_output.dat | awk '{if($1 ~ /-02/){print substr($1,0,7)}}' > last.txt && echo $(tac last.txt | egrep -m 1 .)"_output/"$(tac last.txt | egrep -m 1 .)".king.kin")
echo "$last" >> out.txt

while ( ! ( test -e $last ) ); do sleep 180; done;

