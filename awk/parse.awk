awk 'BEGIN{i=0;while(getline)
	{ 
		trioValues[i++] = $1;
	} 
	for (n = 0; n < i; n++)
	{
		print "awk \47/" trioValues[n] "." trioValues[n] "-(01|02)|" trioValues[n] "-01." trioValues[n] "-02/{print $1, \" \", $3, \" \",  $8}\47 king.kin0 >> parseResult.txt"
		
	}
}' trios/badTrios2_3_17.txt > parseResult.awk