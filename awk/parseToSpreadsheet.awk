awk 'BEGIN{print "proband_blindID		proband_exome_file(s)		mom_blindID		mom_exome_file(s)	dad_blindID	dad_exomefile(s)	trio_validated_via WES (yes/no)";}
NR==FNR{ids[i++]= $1; next} 
NF{ 
	for (i in ids){
		#if($0 ~ ids[i] && ids[i] !~ /(-01|-02)/) {print ids[i] "	" $0; break;}
		if($0 ~ ids[i])
			if(ids[i] !~ /(-01|-02)/) 
			{
				for(x in ids){
					if(ids[x] ~ ids[i] && ids[x] ~ (-01))
						print i ": " ids[i] "	" $0; break;
				}
			}
	}
}' king.kin0 test2.txt > result.txt