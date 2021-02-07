#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
// This set of procedures will read PDB (protein data bank) files:
// https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)
// Contributors:
//	Jan Ilavsky, 2018
//	free to use or modify, no copyright.
//	Major thanks to John Weeks who wrote most of early code when I needed it
//	Thanks to tcaus, Igor Exchange user who described the pdb file structure:
//	reference: https://www.wavemetrics.com/forum/general/igor-does-not-recognize-decim…

// Read PDB file keep carbon backbone only, separate chains as individual waves, make biological assembly

Function IgorPDB()
	ReadPDBFile()
//	ChainsAsWaves()
End

Function ReadPDBFile()
	String filePath = FindPDBFile()
	if (strlen(filePath) == 0)
		return -1
	endif
	String fileName = ParseFilePath(0, filePath, ":", 1, 0)
	String NewFoldername = UniqueName(StringFromList(0,fileName,"."), 11, 0)
	// create location for data
	NewDataFolder/O/S $("root:pdb" + NewFoldername)
	
	String columnInfo = ""				// This only works for rows starting with "ATOM"
	columnInfo += "C=1,W=6,N=keyword;"	// "ATOM" name column								 0
	columnInfo += "C=1,W=6,N='_skip_';"   // Atom number (last column unspecified)	 7
	columnInfo += "C=1,W=5,N=altLoc;"   // Atom name; last character is "alt Loc"   13
	columnInfo += "C=1,W=3,N='_skip_';"   // Residue name; starts at column 18		18
	columnInfo += "C=1,W=2,N=chain;"   // Chain identifier in column 22 (21=?) 21
	columnInfo += "C=1,W=4,N='_skip_';"   // Residue sequence number					  23
	columnInfo += "C=1,W=4,N='_skip_';"   // Blank stuff									  27
	columnInfo += "C=1,W=8,N=xval;"		// Should start at column 31			   31
	columnInfo += "C=1,W=8,N=yval;"		//												 39
	columnInfo += "C=1,W=8,N=zval;"		//												 47
	columnInfo += "C=1,W=6,N='_skip_';"   // Occupancy
	columnInfo += "C=1,W=6,N='_skip_';"   // Temperature info
	columnInfo += "C=1,W=11,N='_skip_';"  // Element symbol   (col 77,78)
	columnInfo += "C=1,W=2,N='_skip_';"   // Charge (col 79,80)
	// Now load the data...
	LoadWave/A/O/F={14,8,0}/K=2/B=columnInfo filePath
	Wave/T xval, yval, zval
	Wave/T keyword, altLoc, chain
	
	// remove all non ATOM lines from the waves
	Variable count = numpnts(keyword) - 1
	
	do 
		if (strsearch(keyword[count], "ATOM", 0) <0 )
			DeletePoints count, 1, keyword, altLoc, chain, xval, yval, zval
		endif
		count -= 1
	while (count >= 0 )

	// keep alpha carbon backbone
	count = numpnts(altLoc) - 1
	
	do 
		if (strsearch(altLoc[count], "CA", 0) <0 )
			DeletePoints count, 1, keyword, altLoc, chain, xval, yval, zval
		endif
		count -= 1
	while (count >= 0 )
	
	// create location for centers of atoms in Gizmo and assign the positions. 
	Make/O/D/N=(numpnts(xval),3) centersWave = NaN
	centersWave[][0] = str2num(xval[p])
	centersWave[][1] = str2num(yval[p])
	centersWave[][2] = str2num(zval[p])
	// cleanup
	KillWaves keyword, altLoc, xval, yval, zval
	
	ChainsAsWaves()
	
	// read symmetry operators
	columnInfo = ""				// parse REMARK lines
	columnInfo += "C=1,W=6,N=keyword;"	// "REMARK" name column								 0
	columnInfo += "C=1,W=4,N='_skip_';"   // Remark number (last column unspecified)	 6
	columnInfo += "C=1,W=9,N=sym;"   // SMTRY or BIOMT   10
	columnInfo += "C=1,W=4,N=op;"   // Operator; starts at column 19
	columnInfo += "C=1,W=10,N=xrot;"		// Should start at column 23
	columnInfo += "C=1,W=10,N=yrot;"		//												 39
	columnInfo += "C=1,W=10,N=zrot;"		//												 47
	columnInfo += "C=1,W=18,N='_skip_';"   // Start at 53
	LoadWave/A/O/F={8,8,0}/K=2/B=columnInfo filePath
	Wave/T xrot, yrot, zrot
	Wave/T keyword, sym, op
	
	// get large symmetry matrix
	count = numpnts(keyword)-1
	
	do 
		if (strsearch(keyword[count], "REMARK", 0) < 0 || strsearch(sym[count], "BIOMT", 0) < 0)
			DeletePoints count, 1, keyword, sym, op, xrot, yrot, zrot
		endif
		count -= 1
	while (count >= 0 )
	KillWaves/z keyword, sym
	Make/O/N=(numpnts(op)) matIndex = str2num(op[p])
	Make/O/N=(numpnts(xrot),3) bigMat
	bigMat[][0] = str2num(xrot[p])
	bigMat[][1] = str2num(yrot[p])
	bigMat[][2] = str2num(zrot[p])
	KillWaves/z op,xrot,yrot,zrot
	
	MakeBiologicalAssembly()
	KillWaves/Z bigMat, matIndex
	
	MakeGizmo(GetDataFolder(1))
	
	SetDataFolder root:
End

STATIC Function/S FindPDBFile()
	Variable refNum
	String fileName
	String filters = "PDB File (*.pdb):.pdb;"
	filters += "All Files:.*;"
	Open/D/R/F=filters/M="Select pdb file" refNum
	fileName = S_fileName			// S_fileName is set by Open/D
	if (strlen(fileName) == 0)		// User cancelled?
		return ""
	endif
	
	return fileName
End

Function ChainsAsWaves()
	Wave/T/Z chain
	Wave/Z centersWave
	// find unique list of chains
	FindDuplicates/RT=uChains/FREE chain
	Variable nChains = numpnts(uChains)
	String newName
	
	Variable i
	
	for(i = 0; i < nChains; i += 1)
		newName = "chain_" + uChains[i]
		Duplicate/O centersWave, $newName
		Wave w = $newName
		w[][] = (cmpstr(chain[p],uChains[i]) == 0) ? w[p][q] : NaN
		MatrixOp/O $newName = zapnans(w)
		Redimension/N=(numpnts(w)/3,3) w
	endfor
	KillWaves/Z centersWave, chain
End

STATIC Function MakeBiologicalAssembly()
	WAVE/Z matIndex, bigMat
	Variable nOps = WaveMax(matIndex) // rotation matrix numbers are 1-based
	String wList = WaveList("chain_*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName, newName

	Variable i,j
	
	for(i = 1; i <= nOps; i += 1)
		Duplicate/O/FREE bigMat, rotMat
		rotMat[][] = (matIndex[p] == i) ? rotMat[p][q] : NaN
		MatrixOp/O rotMat = zapnans(rotMat)
		Redimension/N=(numpnts(rotMat)/3,3) rotMat
		
		for(j = 0; j < nWaves; j += 1)
			wName = StringFromList(j, wList)
			Wave w = $wName
			newName = wName + "_" + num2str(i)
			MatrixOp/O $newName = w x rotMat
		endfor
	endfor
	
	// kill originals
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		KillWaves/Z $wName
	endfor
End

Function MakeGizmo(folderName)
	String folderName
	String plotname = "p_" + folderName
	
	String wList = WaveList("chain_*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName, objName
	NewGizmo/K=1/W=(83,132,981,947)/N=$plotName
	Variable bigNum = 0
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w = $wName
		objName = ReplaceString("chain",wName,"path")
		AppendToGizmo/N=$plotName/D path=$wName,name=$objName
		bigNum = max(bigNum,WaveMax(w))
	endfor
	
	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo setDisplayList=2, object=axes0
	ShowTools
	ModifyGizmo showInfo
	ModifyGizmo infoWindow={904,365,1721,662}
	// scale axes to be isometric and 5% larger than max
	bigNum *= 1.05
	ModifyGizmo setOuterBox={-bigNum, bigNum, -bigNum, bigNum, -bigNum, bigNum}
end