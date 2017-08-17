#! /bin/bash

# jobid=$(bjobs | grep $(date | awk '{print substr($4, 0, 5)}') | awk '{print $1}')
cd /scratch/kannz6/loni/vcf_pipeline/pipeline
! test -d $2 && mkdir $2

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
jid=$(bjobs | grep  $(expr ".* $9" : '.* \(.*\)') | awk '{print $1}' | tac | egrep -m 1 .)
echo "$jid" > loniJobId.txt
# echo "$(tac loniJobId.txt | egrep -m 1 .)" > loniJobId.txt
echo "Last parameter to be used to get loni job-id: $9" >> out.txt

python main.py &

# last=$(cat query_output.dat | awk '{if($1 ~ /-02/){print substr($1,0,7)}}' > last.txt && echo $(tac last.txt | egrep -m 1 .)"_output/"$(tac last.txt | egrep -m 1 .)".king.kin")
last=$(echo $(tac $2/triosStructure.txt | egrep -m 1 .)"_output/"$(tac $2/triosStructure.txt | egrep -m 1 .)".king.kin")
echo "$last" >> out.txt
# while ( ! ( test -e $last ) ); do sleep 180; done;#7-13-17, commented out to test line 52
# while ( ( ! ( test -e $last ) ) && ( test -n "$(bjobs | grep $jid | awk '{print $1}')" ) ); do sleep 180; done;#comment out 7-25-17
#testing this loni job exit 7-25-17
while ( ! test -e triosLength.txt )
do
	sleep 10
done

( test -e "triosLength.txt") && expecting=$( cat triosLength.txt | awk '{print $1}') || [ $( cat triosToVerify.txt | egrep -m 1 .) == "%" ] && scenario=1 || scenario=1
# ( test -e "triosLength.txt") && expecting=$( cat triosLength.txt | awk '{print $1}') || [ $( cat triosToVerify.txt | egrep -m 1 .) == "%" ] && scenario=3 || scenario=1
echo "$scenario" > scenario.txt
echo "expecting: $expecting trios to be validated" >> scenario.txt

# while true; do [ $(date | awk '{print substr($4,4,2)}') == "44" ] && break || echo "still clocking..." && sleep 1; done
# while true; do [ $(date | awk '{print substr($4,0,5)}') == $(echo $oneMinTest | awk '{print substr($4,0,5)}') ] && break || echo "still clocking..." && sleep 1; done
case $scenario in
	1)
	actual=$(awk 'BEGIN{print "0"}')
	errors=0
	if test -n "$expecting"
	then
		while ( test -n "$(bjobs | grep $jid | awk '{print $1}')" && ! ( test $expecting -eq $actual ) && ( test $errors -eq 0 ) )
		do 

			if test -z "$(find -type f -name "*trio-validation-complete.txt" | tac | egrep -m 1 . | awk '{print $1}')"
			then
				actual=$(awk 'BEGIN{print "0"}') 
			else
				actual=$(wc -l $(find -type f -name "*trio-validation-complete.txt") | tac | egrep -m 1 . | awk '{print $1}')
			fi

			sleep 9
			
			test -e jobStatus.txt && echo "expecting[$expecting] vs actual[$actual]" >> jobStatus.txt || echo "expecting[$expecting] vs actual[$actual]" > jobStatus.txt;
			
			if test "$actual" -eq "$expecting"
			then
				errors=0
				kinSink=$2"-kinship"
				inputSink=$2"-inputs"
				test -d $kinSink && ( mkdir $kinSink && cp 1-*/*.kin $kinSink ) || cp 1-*/*.kin $kinSink
				test -d $inputSink && ( mkdir $inputSink && cp lsf/in* $inputSink ) || cp lsf*/in* $inputSink

				mv $kinSink ../../kinshipFiles/
				mv $inputSink ../../inputFiles/

				bsub < ../../../scripts/kinFileParse.sh
				sleep 36
				break
			elif test -n "$(find -type f -name "*log" | tac | egrep -m 1 . )"
				then
					echo "Exiting because log file found!" > error_report.txt
					errors=1
					errorReport="Exiting because plink.log file was found! Run 'find -type f -name \"*log\" | tac | egrep -m 1 . )' to locate file."
					break
			elif test $(echo "$(bjobs | grep "$(echo $jid | awk '{print substr($1,4,3)}')" | awk '{print $1}')" > runningJobs.txt && wc -l runningJobs.txt | tac |awk '{print $1}') -eq 1
				then
					bkill $jid
					errors=1
					errorReport="Only the main lsf job was running. Check both lsf/output* for information, and check email for any issues found!"
					break
			else
				sleep 180
			fi
		done
	else
		while ( ( ! ( test -e "$last" ) ) && ( test -n "$(bjobs | grep $jid | awk '{print $1}')" ) )
	 	do 
			if test -n "$(find -type f -name "*log" | tac | egrep -m 1 . )"
				then
					echo "Exiting because log file found!" > error_report.txt
					errors=1
					errorReport="Exiting because plink.log file was found! Run 'find -type f -name \"*log\" | tac | egrep -m 1 . )' to locate file."
					break
			elif test $(echo "$(bjobs | grep "$(echo $jid | awk '{print substr($1,4,3)}')" | awk '{print $1}')" > runningJobs.txt && wc -l runningJobs.txt | tac |awk '{print $1}') -eq 1
				then
					bkill $jid
					errors=1
					errorReport="Only the main lsf job was running. Check both lsf/output* for information, and check email for any issues found!"
					break
			else
				sleep 180
			fi
		done
	fi
	;;
	2)
	( test -e "triosLength.txt") && expecting=$( cat triosLength.txt | awk '{print $1}') || expecting=$(wc -l triosToVerify.txt | awk '{print $1}')
	while ( test -n "$(bjobs | grep $jid | awk '{print $1}')" )
	do 
		actual=$(wc -l $(find -type f -name "*trio-validation-complete.txt") | tac | egrep -m 1 . | awk '{print $1}')
		if test $actual -eq $expecting
		then 
			break
		else 
			sleep 180
		fi
	done
	;;
	# https://stackoverflow.com/questions/27555727/timeouting-a-while-loop-in-linux-bash-one-line-shell-scripting
	# https://stackoverflow.com/questions/1250079/how-to-escape-single-quotes-within-single-quoted-strings
	# timeout 2m bash -c -- 'while true; do echo "hello world"; sleep 3; done' (local working example)
	#this is case 1, but rewritten with a timeout command
	3)
	actual=$(awk 'BEGIN{print "0"}')
	if test -n $expecting
	then
		while ( test -n "$(bjobs | grep $jid | awk '{print $1}')" && ! ( test $expecting -eq $actual ) )
		do 

			if test -z $(find -type f -name "*trio-validation-complete.txt" | tac | egrep -m 1 . | awk '{print $1}')
			then
				actual=$(awk 'BEGIN{print "0"}') 
			else
				actual=$(wc -l $(find -type f -name "*trio-validation-complete.txt") | tac | egrep -m 1 . | awk '{print $1}')
			fi

			if test $actual -eq $expecting
			then 
				break
			else
				if test -n $(find -type f -name "*log" | tac | egrep -m 1 . )
				then
					echo "Exiting because log file found!" > error_report.txt
					break
				else
					sleep 180
				fi
			fi
		done
	else
		while ( ( ! ( test -e $last ) ) && ( test -n "$(bjobs | grep $jid | awk '{print $1}')" ) )
	 	do 
	 		sleep 180
		done
	fi
	;;
esac

test -z "$errors" && ( echo "Finished grabTrio-dir script!" > finishedJobResult.txt ) || ( echo "Finished grabTrio-dir script with $errors errors" > finishedJobResult.txt && test -n "$errorReport" && ( echo "" >> finishedJobResult.txt && echo "Info: $errorReport" >> finishedJobResult.txt) )
# while ( ( test -n "$(find -type f -name "*.bam" | cat | egrep -m 1 .)" ) && ( test -n "$(bjobs | grep $jid | awk '{print $1}')" ) ); do sleep 180; done;
