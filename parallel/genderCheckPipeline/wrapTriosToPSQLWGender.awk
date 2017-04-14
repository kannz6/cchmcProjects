awk '{print "SELECT core_person.blinded_id AS blinded_id, bulkup_uploadeddatafile.location AS location, core_subject.gender AS gender"
	print "FROM pcgc.bulkup_uploadeddatafile"
print "JOIN pcgc.core_sample"
  print "\tON bulkup_uploadeddatafile.sample_id = core_sample.sample_id"
print "JOIN pcgc.core_person"
  print "\tON core_sample.person_id = core_person.id"
print "LEFT OUTER JOIN pcgc.core_subject"
  print "\t ON core_person.blinded_id = core_subject.blinded_id"
print "WHERE"
  print "\t(bulkup_uploadeddatafile.location ILIKE \47%p1.fastq.gz\47 OR bulkup_uploadeddatafile.location ILIKE \47%p2.fastq.gz\47)"
  print "\tAND bulkup_uploadeddatafile.location LIKE \47%_EX_%\47"
  print "\tAND bulkup_uploadeddatafile.quarantined = false"
  print "\tAND core_person.blinded_id !~ \47-.*-0[3-9]\47"}
{print "\tAND (";i=0;vals[0]=$0;while(getline){vals[++i]=$0}; for(n=0;n < i; n++){print "\tcore_person.blinded_id ILIKE \47"vals[n]"%\47 or";}  print "\tcore_person.blinded_id ILIKE \47"vals[i]"%\47\n\t)"; }
{print "ORDER BY core_person.blinded_id;"}' triosToVerify.txt > trios.sql