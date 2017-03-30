with open("fromLoniTmpSettings.txt", "r") as dynamicSettingsReader:
		tmpFileContent = dynamicSettingsReader.readlines()
dynamicSettingsReader.close()
tmpFileContent = [line.strip() for line in tmpFileContent]

tmpDict = {}
for i in tmpFileContent:
	if ( i != ""):
		i=i.split("=")
		# print "i: %s= %s"%(i[0], i[1])
		if( i[1] != "\"\""):
			tmpDict.update( {i[0] : i[1]} )

with open("versionSettings.py", "r") as versionSettingsR:
		versionFileContent = versionSettingsR.readlines()
versionSettingsR.close()

versionFileContent = [line.strip() for line in versionFileContent]
versionSettingsDict = {}

for s in versionFileContent:
	if ( s != ""):
		s = s.split("=")
		# print "s: %s=%s"%(s[0],s[1])
		if( tmpDict.has_key( s[0] ) and s[1] != tmpDict.get( s[0] ) ):
			print "[changed] %s: %s ---> %s: %s"%(s[0] , s[1] , s[0], tmpDict.get( s[0] ))
			versionSettingsDict.update( { s[0] : tmpDict.get( s[0] ) } )
		elif( len(s) != 1 ):
				versionSettingsDict.update( { s[0] : s[1] } )
				print "[original] %s: %s"%(s[0] , s[1])

newSettings = "#version settings\n\n"
with open("versionSettings.py", "w+") as versionSettingsW:
	versionSettingsW.write(newSettings)
	for k,v in versionSettingsDict.items():
		setting = "{0}={1}\n".format( k , v )
		versionSettingsW.write( setting )
versionSettingsW.close()



