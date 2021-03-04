#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include "IgorPDB"

Menu "Macros"
	Submenu "3IYV Figure"
		"Vertices", /Q, MakeBasicFigure(0)
		"TDs Below Hub", /Q, MakeBasicFigure(1)
		"Vertices Reverse", /Q, MakeBasicFigure(2)
		"TDs Below Hub Reverse", /Q, MakeBasicFigure(3)
		"TDs Below Hub Schematic", /Q, MakeBasicFigure(7)
		"Vertices Schematic Side", /Q, MakeBasicFigure(10)
		"TDs Below Hub Schematic Side", /Q, MakeBasicFigure(11)
	End
End

Function MakeBasicFigure(vIn)
	Variable vIn
	CleanSlate()
	IgorPDB() // find 3IYV.pdb
	RemoveChainsFromDisplay("J;K;L;M;N;O;P;Q;R;")
	DownSampleStructure(5,25)
	FlushColours()
	if((vIn & 2^0) != 0)
		ColourFeetBelow()
	else
		ColourTriskeliaTol()
	endif
	BeautifyGizmo()	
	if((vIn & 2^1) != 0)
		//ModifyGizmo opName=ortho0, operation=ortho,data={1,-1,1,-1,-0.4,-2}
		LookAtReverse()
	endif
	if((vIn & 2^2) != 0)
		// this is the view from the vesicle with the TDs organised to match the schematic diagram
		ModifyGizmo SETQUATERNION={0.271037,-0.923299,0.183335,0.201122}
	endif
	if((vIn & 2^3) != 0)
		// this is the view with the vertex at the top with the TDs organised to match the schematic diagram
		ModifyGizmo SETQUATERNION={0.034713,-0.71831,-0.60628,0.33949}
		ModifyGizmo opName=ortho0, operation=ortho,data={-0.5,0.5,0.2,1.2,-0.5,1}
		ModifyGizmo modifyObject=light0,objectType=light,property={ position,-0.7233,0.5064,0.4695,0.0000}
		ModifyGizmo modifyObject=light0,objectType=light,property={ direction,-0.7233,0.5064,0.4695}
	endif
End

Function FlushColours()
	GetGizmo objectNameList
	String allObj = S_ObjectNames
	Variable nObj = ItemsInList(allObj)
	String oName
	
	Variable i
	
	for(i = 0; i < nObj; i += 1)
		oName = StringFromList(i,allObj)
		if(strsearch(oName,"path",0) >= 0)
			ModifyGizmo ModifyObject=$oName,objectType=path,property={ pathColor,0.733333,0.733333,0.733333,1}
			ModifyGizmo ModifyObject=$oName,objectType=path,property={ pathColorType,1}
			ModifyGizmo ModifyObject=$oName,objectType=path,property={ drawEndcaps,0}
			ModifyGizmo ModifyObject=$oName,objectType=path,property={ calcNormals,1}
			ModifyGizmo ModifyObject=$oName,objectType=path,property={ numSegments,60}
		endif
	endfor
End

Function ColourTriskeliaSochaki()
	// highlighted triskelion
	ColourThesePaths(1,"ABC",132,187,135) // light green
	// facing highlighted thigh
	ColourThesePaths(1,"DEF",40,86,53) // dark green
	ColourThesePaths(2,"ABC",228,163,95) // orange
	ColourThesePaths(6,"ABC",241,226,109) // yellow
	// under the thigh interaction
	ColourThesePaths(5,"ABC",191,149,180) // pink
	ColourThesePaths(2,"DEF",23,24,51) // dark purple
	ColourThesePaths(6,"GHI",122,112,159) // medium purple
End

Function ColourTriskeliaTolBright()
	// highlighted triskelion
	ColourThesePaths(1,"ABC",187,187,187) // grey
	// facing highlighted thigh
	ColourThesePaths(1,"DEF",34,136,51) // green
	ColourThesePaths(2,"ABC",68,119,170) // blue
	ColourThesePaths(6,"ABC",238,102,119) // red
	// under the thigh interaction
	ColourThesePaths(5,"ABC",204,187,68) // yellow
	ColourThesePaths(2,"DEF",170,51,119) // purple
	ColourThesePaths(6,"GHI",102,204,238) // cyan
End

Function ColourTriskeliaTol()
	// highlighted triskelion
	ColourThesePaths(1,"ABC",51,34,136) // Indigo
	// facing highlighted thigh
	ColourThesePaths(1,"DEF",68,170,153) // teal
	ColourThesePaths(2,"ABC",136,204,238) // cyan
	ColourThesePaths(6,"ABC",204,102,119) // rose
	// under the thigh interaction
	ColourThesePaths(5,"ABC",221,204,119) // sand
	ColourThesePaths(2,"DEF",170,68,153) // purple
	ColourThesePaths(6,"GHI",153,153,51) // olive
End

Function ColourFeetBelow()
	ColourThesePaths(1,"ABC",51,34,136) // Indigo
// feet below green triskelion
	ColourThesePaths(1,"GHI",136,34,85) // purple
	ColourThesePaths(6,"DEF",136,34,85) // purple
	ColourThesePaths(4,"ABC",136,34,85) // purple
End


Function ColourThesePaths(sym,triskelion,rr,gg,bb)
	Variable sym
	String triskelion
	Variable rr,gg,bb
	String pList = ""
	pList += "path_" + triskelion[0] + "_" + num2str(sym) + ";"
	pList += "path_" + triskelion[1] + "_" + num2str(sym) + ";"
	pList += "path_" + triskelion[2] + "_" + num2str(sym) + ";"
	
	Variable i
	
	for(i = 0; i < ItemsInList(pList); i += 1)
		ModifyGizmo ModifyObject=$(StringFromList(i,pList)),objectType=path,property={ pathColor,rr/255,gg/255,bb/255,1}
	endfor
End

Function BeautifyGizmo()
	ModifyGizmo scalingMode=8
	ModifyGizmo setOuterBox={-355,355,-355,355,-355,355}
	ModifyGizmo scalingOption=0
	// add light and properties
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo modifyObject=light0,objectType=light,property={ position,0.2962,0.1710,0.9397,0.0000}
	ModifyGizmo modifyObject=light0,objectType=light,property={ direction,0.2962,0.1710,0.9397}
	ModifyGizmo modifyObject=light0,objectType=light,property={ specular,1.000000,1.000000,1.000000,1.000000}
	AppendToGizmo sphere={0.4,25,25},name=sphere0
	ModifyGizmo modifyObject=sphere0,objectType=Sphere,property={colorType,1}
	ModifyGizmo modifyObject=sphere0,objectType=Sphere,property={color,0.7333,0.7333,0.7333,1}
	ModifyGizmo modifyObject=sphere0,objectType=Sphere,property={useGlobalAttributes,0}
	ModifyGizmo modifyObject=sphere0,objectType=Sphere,property={normals,100000}
	AppendToGizmo attribute specular={0.8,0.799954,1.5259e-05,1,1032},name=specular0
	AppendToGizmo attribute shininess={5,5},name=shininess0
	AppendToGizmo attribute diffuse={1,0,0,1,1032},name=diffuse0
	// now enforce order
	ModifyGizmo clearDisplayList
	ModifyGizmo setDisplayList=0, object=light0
	ModifyGizmo setDisplayList=1, attribute=diffuse0
	ModifyGizmo setDisplayList=2, attribute=shininess0
	ModifyGizmo setDisplayList=3, object=sphere0
	GetGizmo objectNameList
	String oName
	Variable nObj = ItemsInList(S_ObjectNames)
	Variable i
	for(i = 0; i < nObj; i += 1)
		oName = StringFromList(i,S_ObjectNames)
		if(stringmatch(oName,"path_*"))
			ModifyGizmo setDisplayList=i+4, object=$oName
		endif
	endfor
	ModifyGizmo currentGroupObject=""
	ModifyGizmo SETQUATERNION={-0.118680,0.251696,0.635953,0.719812}
End

Function LookAtReverse()
	ModifyGizmo appendRotation={1,0,0,0}
	ModifyGizmo modifyObject=light0,objectType=light,property={ position,0.0000,0.0000,-1.0000,0.0000}
	ModifyGizmo modifyObject=light0,objectType=light,property={ direction,0.0000,0.0000,-1.0000}
	ModifyGizmo SETQUATERNION={0.6618148684501648,-0.6959876418113708,0.09318086504936218,0.2625178694725037}
	ModifyGizmo opName=ortho0, operation=ortho,data={-1,1,-1,1,0.4,2}
End