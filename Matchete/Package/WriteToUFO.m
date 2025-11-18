(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete `WriteToUFO`*)


(* ::Subtitle:: *)
(*Create UFO model for Matchete.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


PackageExport["WriteToUFO"]


PackageScope["WritePYParticles"]
PackageScope["WritePyParameters"]
PackageScope["WritePyVertices"]
PackageScope["WritePyLorentz"]
PackageScope["WritePyCouplings"]
PackageScope["WritePyCouplingOrders"]


PackageExport["TranslateToUFO"]
PackageExport["CreateNewDir"]
PackageExport["ParamsTableInput"]
PackageExport["GetDiracStructureSingleTerm"]
PackageScope["ExplicitEinsteinSums"]
PackageExport["ExplicitFlavor"]
PackageExport["DefaultParamcard"]
PackageExport["InputFile"]


PackageExport["GetCouplingOrders"]
PackageExport["DefineCouplingOrder"]
PackageExport["RemoveCouplingOrder"]
PackageExport["ResetCouplingOrders"]
PackageExport["ExtractCouplingPower"]
PackageScope["ExtractCouplingLabel"]
PackageExport["ClassifyCouplingOrder"]


PackageScope["HasMassQ"]
PackageScope["PowerOf"]
PackageExport["RenameDummyIndicesToUFO"]


(* ::Section:: *)
(*Usage messages*)


PerturbativeExps::usage = "PerturbativeExps is an option for the routine DefineCouplingOrder that specifies the maximum power of the interaction that could appear in a single amplitude."
DefineCouplingOrder::usage = "DefineCouplingOrder[name, couplings, expansionOrder, hierarchy, opts] defines a coupling order based on the given coupling(s). The expansion order and hierarchy should be stated for numerical reasons."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Main*)


(* ::Subsection:: *)
(*WriteToUFO*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


WriteToUFO::nofile = "File not found: `1`";
WriteToUFO::missingdata = "Definitions for the following parameters are missing: `1`. MATCHETE only checks for parameters of the lagrangian (couplings, etc.)";
WriteToUFO::missingpdg = "pdg code for particle `1` is missing. Please define fields with DefineField[..., UFO$Options -> <|\"pdg\" -> ...|>] and gauge fields with DefineGaugeGroup[..., UFO$Options -> <|\"pdg\" -> ...|>].";


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[WriteToUFO] = {
	InputFile -> None,
	Tami -> False
};


(* ::Subsubsection:: *)
(*WriteToUFO*)


WriteToUFO[lagrangian_, opt:OptionsPattern[]] := Module[
	{
	vertices,
	allFields,
	vertexData = {},
	dir,
	fields,
	fieldTypes,
	expandedVertex,
	indexTypes,
	lorentzStructure,
	groupFlavorStructure,
	structure,
	couplings,
	lorentzData,
	couplingData,
	parametersNeeded,
	parametersExternal,
	parametersFile = OptionValue[InputFile]
	},
	
	(*Import Parameters data*)
	If[
		parametersFile === None,
		parametersExternal = {"parameters" -> {}},
		If[
			FileExistsQ[parametersFile],
			parametersExternal = Import[parametersFile] /. ("lhacode" -> {x_}) :> ("lhacode" -> x),
			Message[WriteToUFO::nofile, parametersFile];
			Return[$Failed]
		]
	];
	
	(*Calculate Feynman rules*)
	If[
		OptionValue[Tami],
		vertices = FeynmanRuleFinderAll[lagrangian],
		vertices = FeynmanRules[lagrangian]
	];
	
	allFields = DeleteDuplicates @ (GetFieldsInTerm[vertices] //. {
			Bar[x_] :> x,
			OverBar[x_] :> x,
			Subscript[x_, _] :> x,
			Field[EXT[x_], ___] :> x
	});
	
	Do[
		fields = SortBy[DeleteDuplicates @ GetFieldsInTerm[vertex], Cases[#, _Integer, \[Infinity]] &];
		fieldTypes = StringTake[SymbolName[If[AtomQ[#], #, Head[#]]], 1] & /@ (fields /. Bar[x_] -> x)[[All, 2]];
		expandedVertex = Expand @ vertex;
		indexTypes = DeleteDuplicates @ Replace[Flatten[GetFields[#][Indices] & /@ Keys @ GetFields[]], s_Symbol[___] :> s, {1}];
		
		(*Extract lorentz/spinor structures*)
		If[
			Head[expandedVertex] === Plus
			,
			lorentzStructure = (Times @@ Join[GetDiracStructureSingleTerm[#], Cases[#, Metric[___]|Subsuperscript["p", ___]|LCTensor[___]]]&) /@ List @@ expandedVertex;
			expandedVertex = (List @@ expandedVertex) /. Field[___]|Metric[___]|Subsuperscript["p", ___]|LCTensor[___]|DiracProduct[___] :> 1
			,
			lorentzStructure = {Times @@ Join[GetDiracStructureSingleTerm[expandedVertex], Cases[expandedVertex, Metric[___]|Subsuperscript["p", ___]|LCTensor[___]]]};
			expandedVertex = {expandedVertex /. Field[___]|Metric[___]|Subsuperscript["p", ___]|LCTensor[___]|DiracProduct[___] :> 1}
		];
		
		(*Extract all remaining structures (gauge groups, flavor)*)
		groupFlavorStructure = Association[
			Table[
				ToString @ indexType -> Map[
					Function[
						term, Module[
							{factors},
							(* Break the product into its factors if it is a Times product *)
							factors = If[Head[term] === Times, List @@ term, {term}];
							(* Recombine only those factors that contain the current index type *)
							Times @@ Select[factors, Not[FreeQ[#, indexType]] &]
						]
					],
					expandedVertex
				],
				{indexType, indexTypes}
			]
		];
		
		structure = Append[Map[RenameDummyIndicesToUFO, groupFlavorStructure, {2}], "Lorentz" -> lorentzStructure];
		
		(*The remaining part of the Feynman rule is the coupling*)
		couplings = expandedVertex /. CG[___]|Delta[___] :> 1;
		AppendTo[structure, "Couplings" -> couplings];
		
		AppendTo[structure, "FieldsContent" -> StringJoin @ Sort[fieldTypes]];
		
		AppendTo[structure, "Fields" -> fields //. {
			(*Bar[x_] :> x,*)
			(*OverBar[x_] :> x,*)
			Subscript[x_, _] :> x,
			Field[EXT[x_], ___] :> x
		}];
		
		fieldTypeTranslation = <|"V" -> 3, "S" -> 1, "F" -> 2, "G" -> -1|>;
		AppendTo[structure, "FieldsSorted" -> Replace[fieldTypes, fieldTypeTranslation, {1}]];
		
		AppendTo[structure, "Feynmanrule" -> vertex];
		
		AppendTo[vertexData, structure]
		,
		{vertex, vertices}
	];
	
	vertexData = vertexData /. Subsuperscript["p", \[Mu]_, Subscript[_, n_]] :> Subsuperscript["p", \[Mu], n]; 
	
	(*Create a list of unique couplings*)
	couplingData = DeleteDuplicates @ Flatten[
		Table[
			<|"Coupling" -> \[ScriptL]|>,
			{ver, vertexData},(* loop over each association *)
			{\[ScriptL], ver["Couplings"]}(* loop over each element of its Lorentz list *)
		]
	];
	
	(*Create a list of unique Lorentz structures*)
	lorentzData = DeleteDuplicates @ Flatten[
		Table[
			<|"Lorentz" -> \[ScriptL], "FieldsContent" -> ver["FieldsContent"], "FieldsSorted" -> ver["FieldsSorted"]|>,
			{ver, vertexData},(* loop over each association *)
			{\[ScriptL], ver["Lorentz"]}(* loop over each element of its Lorentz list *)
		]
	];
	
	(*Add names  to Lorentz structures*)
	lorentzData = Flatten[
		KeyValueMap[
			Function[
				{fieldCounts, group},
				MapIndexed[
					Append[#, "Name" -> fieldCounts <> ToString[First[#2]]] &,
					group
				]
			],
			GroupBy[lorentzData, #["FieldsContent"]&]
		]
	];
	
	(*Replace each Lorentz structure in vertex data with its corresponding name*)
	lorentzRules = Association[
		Table[
			{l["Lorentz"], l["FieldsContent"]} -> l["Name"],
			{l, lorentzData}
		]
	];

	vertexData = Map[
		Function[v,
			Module[
				{fieldsContent = v["FieldsContent"], oldL = v["Lorentz"], newL},
				(* produce a new list of names, one per old lorentz structure *)
				newL = Lookup[
					lorentzRules,
					Thread[{oldL, ConstantArray[fieldsContent, Length@oldL]}],
					oldL  (* fallback: keep original if no rule found *)
				];
				(* return a new association which is like v but with "Lorentz"->newL *)
				ReplacePart[v, "Lorentz" -> newL]
			]
		],
		vertexData
	];
	
	(*Drop any broken symmetry*)
	brokenSymmetries = Select[
		Keys[First[vertexData]],
		symmetry |-> AllTrue[vertexData, MatchQ[# [symmetry], {1 ..}] &]
	];

	vertexData = KeyDrop[#, brokenSymmetries] & /@ vertexData;
	
	(*Edit coupling data*)
	couplingData = MapIndexed[
		Join[
			#,
			<|
				"Name" -> "GC_" <> ToString[First @ #2],
				"UFO" -> TranslateToUFO[#["Coupling"]],
				"CouplingOrder" -> StringJoin @ {
					"{", 
					StringRiffle[
						"'" <>
						ToString[#1, InputForm] <> 
						"':" <> 
						ToString[#2, InputForm] & @@@ Select[Normal[ClassifyCouplingOrder[#["Coupling"]]], #[[2]] =!= 0 &],
						", "
					],
					"}"
				}
			|>
		] &,
		couplingData
	];
	
	(*Replace each Coupling in the vertex data with its corresponding name*)
	couplingRules = Association[ 
		(#["Coupling"] -> #["Name"]) & /@ couplingData
	];
	
	vertexData = Map[MapAt[Lookup[couplingRules, #, #] &, "Couplings"] @ # &, vertexData];

	(*Translate from MATCHETE convention to UFO format*)
	lorentzData = Map[MapAt[TranslateToUFO, #, Key["Lorentz"]] &, lorentzData /. (index_Symbol /; StringMatchQ[SymbolName[index], "d$$"~~__] :> ToExpression @ StringDrop[SymbolName[index], 3])];
	vertexData = Map[MapAt[TranslateToUFO /@ # &, #, Key["SU3c"]] &, vertexData /. (index_Symbol /; StringMatchQ[SymbolName[index], "d$$"~~__] :> ToExpression @ StringDrop[SymbolName[index], 3])];
	
	(*Get all used parameters of the theory (First extract naked symbols and then Coupling objects with indices)*)
	parametersNeeded = Union[
		Flatten[
			TranslateToUFO /@ Cases[
				DeleteCases[Lookup[couplingData, "Coupling"], c_Coupling, Infinity],
				s_Symbol /; Context[s] === "Global`",
				Infinity
			]
		],
		TranslateToUFO /@ Cases[Lookup[couplingData, "Coupling"], c_Coupling, Infinity]
	];

(*	Print[vertexData//NiceForm];
	Print[couplingData//NiceForm];
	Print[lorentzData//NiceForm];*)

	
	(*Writing files*)
	PY$InputFiles = {"object_library.py", 
                 "__init__.py",
                 "function_library.py",
                 "write_param_card.py",
                 "propagators.py"
                };
	
	dir = CreateNewDir[];

	CopyFile[$MatchetePath <> "Package/StandardFilesUFO/" <> #, dir <> "/" <> #]& /@ PY$InputFiles; 
	
	SetDirectory[dir];
	
	(*Write particles.py*)
	WritePYParticles[allFields, lagrangian];
	
	(*Write parameters.py*)
	WritePyParameters[parametersNeeded, parametersExternal];
	
	(*Write vertices.py*)
	WritePyVertices[vertexData];
	
	(*Write lorentz.py*)
	WritePyLorentz[lorentzData];
	
	(*Write couplings.py*)
	WritePyCouplings[couplingData];
	
	(*Write coupling_orders.py*)
	WritePyCouplingOrders[GetCouplingOrders[]];
	
	(*Write*)
	ResetDirectory[];
	Print["Done!"];
];


(* ::Subsection:: *)
(*Writing *)


(* ::Subsubsection:: *)
(*Particles*)


WritePYParticles[fields_, lagrangian_] := Module[
	{
	outfile,
	massAssociation,
	spinAssociation,
	colorAssociation,
	listKeys,
	fieldAttributesList,
	outputString
	},
	
	outfile = "particles.py";
	
	(*Get all relevant data*)
	massAssociation = HasMassQ[lagrangian];
	spinAssociation = <|Scalar->"1", Fermion->"2", Vector->"3", Ghost->"-1"|>;
	colorAssociation = <|fund->"3", adj->"8"|>;
	listKeys = {"pdg_code", "name", "antiname", "spin", "color", "mass", "width", "texname", "antitexname", "charge"};
	fieldAttributesList = KeyValueMap[
		{
		If[KeyExistsQ[#2[UFO$Options], "pdg"], ToString[#2[UFO$Options]["pdg"]], Message[WriteToUFO::missingpdg, #1]; Abort[]]
		(*ToString@Input["PDG Monte Carlo code for "<>ToString@#1<>":"]*),
		If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]], 
		If[#2[SelfConjugate], If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]], If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]] <> "~"], 
		spinAssociation[#2[Type]], 
		If[FreeQ[#2[Indices], Global`SU3c], "1", colorAssociation[Cases[#2[Indices], Global`SU3c[x_] :> x][[1]]]], 
		If[massAssociation[#1], "Param.M"<>If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]], "Param.ZERO"], 
		If[massAssociation[#1], "Param.W"<>If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]], "Param.ZERO"], 
		If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]],
		If[#2[SelfConjugate] == True, If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]], If[KeyExistsQ[#2[UFO$Options], "name"], ToString @ #2[UFO$Options]["name"], TranslateToUFO[#1]] <> "~"],
		If[FreeQ[#2[Charges], Global`U1em], "0", ToString[First[Cases[#2[Charges], Global`U1em[x_] :> x]], InputForm]],
		#2[SelfConjugate]
		}&, KeyTake[GetFields[], fields]
	];
	
	(*Create the output string*)
	outputString = StringJoin @ Table[
		StringJoin[
			fieldAttributes[[2]],
			" = Particle( \n", 
			StringJoin[
				Table[
					StringJoin[
						"   ",
						listKeys[[i]],
						If[
							StringContainsQ[listKeys[[i]], "name"],
							" = '" <> fieldAttributes[[i]] <> "',\n",
							" = " <> fieldAttributes[[i]] <> ",\n"
						]
					],
					{i, Length[fieldAttributes] - 1}
				],
				")\n\n",
				If[
					!(Last @ fieldAttributes),
					StringJoin[
						fieldAttributes[[2]],
						"__tilde__ = ",
						fieldAttributes[[2]],
						".anti()\n\n"
					],
					""
				]
			]
		],
		{fieldAttributes, fieldAttributesList}
	];
	
	WriteString[outfile, "from __future__ import division\n"];
	WriteString[outfile, "from object_library import all_particles, Particle\n"];
	WriteString[outfile, "import parameters as Param\n\n"];
	WriteString[outfile, "import propagators as Prop\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
]


(* ::Subsubsection::Closed:: *)
(*Vertices*)


WritePyVertices[data_] := Module[
	{
	outfile = "vertices.py",
	headers,
	outputString
	},
	
	outputString = StringJoin @ MapIndexed[
		Module[
			{vertex = #1, i = First@#2},
			
			StringJoin[
				"V_",
				ToString[i],
				" = Vertex(\n   name = 'V_",
				ToString[i],
				"',\n   particles = ",
				"[ ",
				Riffle["P." <> If[KeyExistsQ[GetFields[# /. OverBar[x_] :> x][UFO$Options], "name"], ToString @ GetFields[# /. OverBar[x_] :> x][UFO$Options]["name"], TranslateToUFO[# /. OverBar[x_] :> x]] <> If[MatchQ[#, OverBar[___]], "__tilde__", ""] & /@ vertex["Fields"], ", "],
				" ],\n   color = [ ",
				Riffle["'" <> ToString@# <> "'" & /@ vertex["SU3c"], ", "],
				" ],\n   lorentz = [ ",
				Riffle["L." <> StringReplace[SpokenString[#], " " -> ""] & /@ vertex["Lorentz"], ", "],
				" ],\n   couplings = {",
				StringRiffle[
					MapIndexed[
						Function[
							{name, idx}, 
							"(" <> ToString[idx[[1]] - 1] <> "," <> ToString[idx[[1]] - 1] <> "):C." <> ToString[name]
						],
						vertex["Couplings"]
					],
					", "
				],
				"}\n)\n"
			]
		] &,
		data
	];
	
	WriteString[outfile, "from object_library import all_vertices, Vertex\n"];
	WriteString[outfile, "import particles as P\n"];
	WriteString[outfile, "import couplings as C\n"];
	WriteString[outfile, "import lorentz as L\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
]


(* ::Subsubsection::Closed:: *)
(*Parameters*)


WritePyParameters[defaultParams_List, defaultData_List] := Module[
	{
	outfile,
	params,
	outputString,
	headers
	},
	
	outfile = "parameters.py";
	
	(*Get all relevant parameters via UI (need to introduce standard parameters for SM and SMEFT)*)
	headers = {"name", "nature", "type", "value", "texname", "lhablock", "lhacode"};
	params = Select[ParamsTableInput[headers, defaultData, defaultParams], !AllTrue[#, # == "" &] &]; (*Remove empty sublists*)
	params = DeleteCases[#, "-" | ""] & /@ params; (*Delete all "-" to ensure correct formatting of internal and external degrees of freedom*)

	outputString = StringJoin @ Table[
		StringJoin[
			paramValues[[1]],
			" = Parameter( \n", 
			StringJoin[
				Table[
					StringJoin[
						"   ",
						headers[[i]],
						Switch[
							headers[[i]],
							"value",
							If[
								paramValues[[2]] === "external",
								" = " <> paramValues[[i]] <> ",\n",
								" = '" <> paramValues[[i]] <> "',\n"
							],
							"lhacode",
							" = [ " <> paramValues[[i]] <> " ],\n",
							_,
							" = '" <> paramValues[[i]] <> "',\n"
						]
					],
					{i, Length[paramValues]}
				],
				")\n"
			]
		],
		{paramValues, params}
	];
	
	WriteString[outfile, "from object_library import all_parameters, Parameter\n"];
	WriteString[outfile, "from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot\n\n"];
	WriteString[outfile, StringJoin[
		"ZERO = Parameter(\n",
		"   name = 'ZERO',\n",
		"   nature = 'internal',\n",
		"   type = 'real',\n",
		"   value = '0.0',\n",
		"   texname = '0'\n",
		")\n"
		]
	];
	WriteString[outfile, outputString];
	Close[outfile];
];


(* ::Subsubsection::Closed:: *)
(*Lorentz*)


WritePyLorentz[data_] := Module[
	{
	outfile = "lorentz.py",
	outputString
	},
	 
	outputString = StringJoin @ Table[
		StringJoin[
			structure["Name"],
			" = Lorentz(\n   name = '",
			structure["Name"],
			"',\n   spins = [ ",
			StringJoin[Riffle[ToString /@ structure["FieldsSorted"], ", "]],
			" ],\n   structure = '",
			structure["Lorentz"],
			"'\n)\n\n"
		],
		{structure, data}
	];
	
	WriteString[outfile, "from object_library import all_lorentz, Lorentz\n"];
	WriteString[outfile, "from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot\n"];
	WriteString[outfile, "try:\n   import form_factors as ForFac\nexcept ImportError:\n   pass\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
];


(* ::Subsubsection::Closed:: *)
(*Couplings*)


WritePyCouplings[data_] := Module[
	{
	outfile = "couplings.py",
	outputString
	},
	
	outputString = StringJoin @ Table[
		StringJoin[
			coupling["Name"],
			" = Coupling(\n   name = '",
			coupling["Name"],
			"',\n   value = '",
			coupling["UFO"],
			"',\n   order = ",
			coupling["CouplingOrder"],
			"\n)\n\n"
		],
		{coupling, data}
	];
	
	WriteString[outfile, "from object_library import all_couplings, Coupling\n"];
	WriteString[outfile, "from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot\n\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
]


(* ::Subsubsection::Closed:: *)
(*Coupling Orders*)


WritePyCouplingOrders[data_] := Module[
	{
	outfile = "coupling_orders.py",
	outputString
	},
	
	outputString = StringJoin @ Table[
		StringJoin[
			ToString@couplingOrder[Name],
			" = CouplingOrder(\n   name = '",
			ToString@couplingOrder[Name],
			"',\n   expansion_order = ",
			ToString@couplingOrder[ExpansionOrder],
			",\n   hierarchy = ",
			ToString@couplingOrder[Hierarchy],
			If[
				couplingOrder[PerturbativeExps] =!= None,
				",\n   perturbative_expansion = " <> ToString[couplingOrder[PerturbativeExps]],
				""
			],
			"\n)\n\n"
		],
		{couplingOrder, data}
	];
	
	WriteString[outfile, "from object_library import all_orders, CouplingOrder \n\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
]


(* ::Section:: *)
(*Help Functions*)


(* ::Subsection::Closed:: *)
(*Miscellaneous*)


CreateNewDir[] := Module[
	{
	parent, 
	newDirName, 
	newDirPath
	},
	
	(* Let the user choose the parent directory *)
	parent = SystemDialogInput["Directory"];
	If[parent === $Canceled, Abort[]];
	
	(* Prompt the user for the new directory name *)
	newDirName = InputString["Enter the name for the new directory:"];
	If[!StringQ[newDirName] || newDirName === "", 
    Return[Message[InputString::nval, newDirName]]
	];
	
	(* Create the new directory if it doesn't already exist *)
	newDirPath = FileNameJoin[{parent, newDirName}];
	If[DirectoryQ[newDirPath],
		Print["Directory already exists at: ", newDirPath];
		Abort[],
		CreateDirectory[newDirPath];
		Print["Directory created at: ", newDirPath]
	];
	newDirPath
];


ExtractCouplingPower[expr_, couplingName_] := Module[
	{
	dummy,
	expanded
	},
	
	dummy = Unique["dummy"];
	expanded = Expand[expr /. Coupling[couplingName, ___] -> dummy];
	Max[Exponent[#, dummy] & /@ List[expanded]]
];


(* Not very efficient *)
HasMassQ[lagrangian_] := Module[
	{
	fields = Keys @ GetFields[],
	terms
	},
	
	(* Search for terms with 2 fields*)
	terms = AssociationMap[
		symbol |-> Select[
			List @@ Expand[lagrangian],
			(
				(* complex fields *)
				Count[#, symbol, {0, \[Infinity]}] === 2 &&
				Count[#, _FieldStrength, {0, \[Infinity]}] === 0 &&
				Count[#, _Field, {0, \[Infinity]}] === 2 &&
				Count[#, Field[symbol, ___]^2, {0, \[Infinity]}] === 0
			)||(
				(* real fields *)
				Count[#, Field[symbol, ___]^2, {0, \[Infinity]}] === 1 &&
				Count[#, _FieldStrength, {0, \[Infinity]}] === 0 &&
				Count[#, Field[___], {0, \[Infinity]}] === 1
			)&
		], 
		fields
	];
	
	(* Classify if massless or not *)
	AssociationMap[
		Function[
			label,
			Module[
				{fieldType, numEntries},
				fieldType = GetFields[label][Type];
				numEntries = Length @ terms[label];
				If[
					fieldType === Vector,
					numEntries >= 1,
					numEntries >= 2
				]
			]
		],
		Keys[terms]
	]
]


ExtractCouplingLabel[expr_] := Module[
	{head},
	head = Head[expr];
	Which[
		head === Symbol,
		expr,
		head === Power || head === Times,
		ExtractCouplingLabel[First[expr]],
		head === Coupling,
		First[expr],
		True,
		Missing["UnknownForm",expr]
	]
]


PowerOf[expr_,sym_] := Module[
	{
	factors,
	exp
	},
	(*1) Turn a product into a list of factors*)
	factors = If[Head[expr] === Times, List @@ expr, {expr}];
	(*2) Sum up all explicit Powers of sym*)
	exp = Total @ Cases[factors, Power[sym, n_] :> n];
	(*3) Count any bare sym\[CloseCurlyQuote]s (each is+1)*)
	exp += Count[factors, sym];
	exp
]


(* ::Subsection::Closed:: *)
(*Table Input for Parameters*)


ParamsTableInput[headers_List, defaultData_List : {"parameters" -> {}}, defaultParamsSymbols_List : {}] := DialogInput[DynamicModule[
    {
    defaultParamsData,
    table,
    paramsWOData,
    optionsNature,
    optionsType,
    addRow
    },
    
    (* Extract data out of defaultData*)
    defaultParamsData = Lookup[defaultData, "parameters"];
    table = Table[
		ToString @ Lookup[Lookup[defaultParamsData, paramName], header, ""],
		{paramName, Keys[defaultParamsData]},
		{header, headers}
    ];
	
    (* Check if all MATCHETE parameters get a value from the default data if applicable *)
    paramsWOData = Complement[defaultParamsSymbols, Keys @ defaultParamsData];
    
    If[
		defaultParamsData =!= {} && paramsWOData =!= {},
		Message[WriteToUFO::missingdata, paramsWOData];
		Abort[]
    ];
	
	(* Initialize table from parameters out of the lagrangian *)
    table = Join[
		table,
		PadRight[List @ #, Length[headers], ""] & /@ paramsWOData
	];
	
    (* Dropdown options *)
    optionsNature = {"external", "internal"};
    optionsType = {"real", "complex"};

    (* Add Row *)
    addRow[] := (table = Append[table, ConstantArray["", Length[headers]]]);

    (* UI *)
    Column[
        {Pane[Dynamic[Grid[
            Prepend[
                Table[
                    With[
                        {ii = i, jj = j},
                        DynamicModule[{},
                            Switch[
                                jj,
                                2, (* Nature Column *)
                                PopupMenu[
                                    Dynamic[
                                        table[[ii, jj]],
                                        (table[[ii, jj]] = #;
                                         Which[
											# === "external", 
											table[[ii, 3]] = "real",
											# === "internal",
											table[[ii, 6]] = "-";
											table[[ii, 7]] = "-"]) & (* Fix type to real *)
                                    ],
                                    optionsNature
                                ],
                                3, (* Type Column *)
                                PopupMenu[
                                    Dynamic[
                                        table[[ii, jj]],
                                        (If[table[[ii, 2]] === "external", table[[ii, jj]] = "real", table[[ii, jj]] = #]) &
                                    ],
                                    optionsType,
                                    Enabled -> (table[[ii, 2]] =!= "external") (* Disable type selection if external *)
                                ],
                                6 | 7, (* lhablock and lhacode Columns *)
                                InputField[
                                    Dynamic[table[[ii, jj]]],
                                    String,
                                    Enabled -> (table[[ii, 2]] =!= "internal") (* Disable when internal *)
                                ],
                                _, (* Default case *)
                                InputField[Dynamic[table[[ii, jj]]], String]
                            ]
                        ]
                    ],
                    {i, Length[table]},
                    {j, Length[headers]}
                ],
                headers
            ],
            Frame -> All
        ]],
		{Automatic, 300},
		Scrollbars -> True
		],
        Row[
            {
            Button["Add Row", addRow[]],
            Spacer[10],
            Button["Done", DialogReturn[table]]
            }
        ]}
    ]
]];


(* ::Subsection::Closed:: *)
(*GetDiracStructure*)


(* Contraction fermion*fermion *)
GetDiracStructureSingleTerm[term_] /; Cases[term, DiracProduct[___], {0, Infinity}] == {} && !FreeQ[term, Fermion, {0, Infinity}] := Module[
	{
	spinorStructure,
	spinorFields
	},
	
	spinorStructure = Cases[term, NonCommutativeMultiply[___], {0, Infinity}];
	spinorFields = Cases[spinorStructure, Field[___], {0, Infinity}] /. Bar[x_] :> x;
	
	{Delta[Index[spinorFields[[1, 1, 1, 2]], Spinor], Index[spinorFields[[2, 1, 1, 2]], Spinor]]}
];

(* Contraction fermion*spinorstrcture*fermion *)
GetDiracStructureSingleTerm[term_] /; Cases[term, DiracProduct[___], {0, Infinity}] =!= {} := Module[
	{
	diracStructures,
	positionsDirac,
	diracStructureFields,
	listInd = {},
	dummyIndexNumber = 1,
	result,
	dummyrules
	},
	
	diracStructures = Cases[term, DiracProduct[___], {0, Infinity}];
	positionsDirac = Position[term, _?(MemberQ[diracStructures, #] &), Infinity];
	
	Do[
		diracStructureFields = Cases[Extract[term, Most @ dirac], Field[___], {0, Infinity}] /. Bar[x_] :> x;
		Switch[
			Length[Extract[term, dirac]],
			1,
			AppendTo[listInd, {Index[diracStructureFields[[2, 1, 1, 2]], Spinor], Index[diracStructureFields[[1, 1, 1, 2]], Spinor]}],
			2,
			listInd = Join[
				{{Index[diracStructureFields[[2, 1, 1, 2]], Spinor], Index[-dummyIndexNumber, Spinor]}, {Index[-dummyIndexNumber, Spinor], Index[diracStructureFields[[1, 1, 1, 2]], Spinor]}},
				listInd
			];
			dummyIndexNumber++,
			_,
			AppendTo[listInd, {Index[diracStructureFields[[2, 1, 1, 2]], Spinor], Index[-dummyIndexNumber, Spinor]}];
			Do[
				AppendTo[listInd, {Index[-dummyIndexNumber, Spinor],Index[-dummyIndexNumber-1, Spinor]}];
				dummyIndexNumber++,
				{i, Length[Extract[term, dirac]] - 2}
			];
			AppendTo[listInd, {Index[-dummyIndexNumber, Spinor], Index[diracStructureFields[[1, 1, 1, 2]], Spinor]}];
			dummyIndexNumber++
		],
		{dirac, positionsDirac}
	];
	result = MapThread[
		Function[
			{expr,idxList},
			If[
				Head[expr] === Symbol,
				expr[idxList],
				Head[expr][Sequence @@ (List @@ expr), idxList]
			]
		],
		{Flatten[List @@@ diracStructures], listInd}
	];
	
	(* Replace potential Lorentz dummies *)
	dummyrules = MapIndexed[
		#1 -> Index[-dummyIndexNumber + First[#2] - 1, Lorentz]&, 
		Cases[FindDummyIndices[result], Index[_, Lorentz]]
	];
	
	result /. dummyrules
];

(* Exception *)
GetDiracStructureSingleTerm[term_] /; Cases[term, DiracProduct[___], {0, Infinity}] == {} := {};


(* ::Subsection:: *)
(*UFO translation functions*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


TranslateToUFO::unidentified = "Structure cannot be identified: `1`"


(* ::Subsubsection:: *)
(*TranslateToUFO*)


(* Metrics *)
TranslateToUFO[expr_Metric] := Module[
	{inds},
	inds = expr /. Metric[Index[d_, Lorentz], Index[e_, Lorentz]] :> {d, e};
	"Metric(" <> StringJoin[Riffle[ToString /@ inds, ","]] <> ")"
];

(* Gamma matrices *)
TranslateToUFO[expr_GammaM] := Module[
	{inds}, 
	inds = expr /. GammaM[Index[d_, Lorentz], {Index[e_, Spinor], Index[f_, Spinor]}] :> {d, e, f};
	"Gamma(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Charge conjugation operator *)
TranslateToUFO[expr_GammaCC] := Module[
	{inds}, 
	inds = expr /. GammaCC[{Index[e_, Spinor], Index[f_, Spinor]}] :> {e, f};
	"Identity(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Projectors *)
TranslateToUFO[expr_Proj] /; MatchQ[expr, Proj[-1, _]] := Module[
	{inds}, 
	inds = expr /. Proj[_, {Index[e_, Spinor], Index[f_, Spinor]}] :> {e, f};
	"ProjM(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];
TranslateToUFO[expr_Proj] /; MatchQ[expr, Proj[1, _]] := Module[
	{inds}, 
	inds = expr /. Proj[_, {Index[e_, Spinor], Index[f_, Spinor]}] :> {e, f};
	"ProjP(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Momenta *)
TranslateToUFO[expr_] /; MatchQ[expr, Subsuperscript["p", ___]] := Module[
	{inds}, 
	inds = expr /. Subsuperscript["p", Index[e_, Lorentz], f_] :> {e, f};
	"P(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Deltas *)
TranslateToUFO[Delta[Index[a_, _], Index[b_, _]]] := "Identity(" <> ToString@a <> "," <> ToString@b <> ")";

(* Levi Civita Tensors*)
TranslateToUFO[LCTensor[Index[a_, _], Index[b_, _], Index[c_, _], Index[d_, _]]] := "Epsilon(" <> StringJoin[Riffle[ToString /@ {a, b, c, d}, ","]] <> ")";

(* Structure functions*)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[fStruct[_], ___]] := Module[
	{inds},
	inds = expr /. CG[fStruct[_], {Index[d_, _], Index[e_, _], Index[f_, _]}] :> {d, e, f};
	"f(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Generators of SU(n) *)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[gen[_], ___]] := Module[
	{inds},
	inds = expr /. CG[gen[_], {Index[d_, _], Index[e_, _], Bar[Index[f_, _]]}] :> {d, f, e};
	"T(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Couplings (placeholder)*)
TranslateToUFO[c_Coupling] := StringJoin[
	TranslateToUFO @ c[[1]],
	Sequence @@ (ToString /@ c[[2]])
];

(* Operations *)
TranslateToUFO[Times[a_, b__]] := StringJoin @ Riffle[TranslateToUFO /@ {a, b}, "*"];
TranslateToUFO[DiracProduct[a___]] := StringJoin @ Riffle[TranslateToUFO /@ {a}, "*"];
TranslateToUFO[Rational[a_, b_]] := StringJoin @ Riffle[TranslateToUFO /@ {a, b}, "/"];
TranslateToUFO[Power[a_, 1/2]] := "cmath.sqrt(" <> TranslateToUFO[a] <> ")";
TranslateToUFO[Power[a_, -1/2]] := "1/cmath.sqrt(" <> TranslateToUFO[a] <> ")";
TranslateToUFO[Power[a_, exponent_]] := TranslateToUFO[a] <> "**(" <> TranslateToUFO[exponent] <> ")";
TranslateToUFO[Bar[a_]] := "complexconjugate(" <> TranslateToUFO@a <> ")";

(* Standard expressions *)
TranslateToUFO[1] := "1";
TranslateToUFO[0] := "0";
TranslateToUFO[n_Integer] := ToString @ n <> ".";
TranslateToUFO[n_Symbol] := StringReplace[SpokenString[n], " " -> ""];
TranslateToUFO[z_Complex] := "complex("<>TranslateToUFO[Re @ z]<>","<>TranslateToUFO[Im @ z]<>")";

(* Error *)
TranslateToUFO[expr_] := Message[TranslateToUFO::unidentified, expr];


(* ::Subsection:: *)
(*ExplicitFlavorStructure*)


ExplicitEinsteinSums[sum_Plus]:=ExplicitEinsteinSums/@sum

ExplicitEinsteinSums[term:Except[_Plus]] := Module[
	{
	repeatedInds,
	flavors = GetFlavorIndices[]
	},
	
	(* determine repeated indices *)
	repeatedInds = Matchete`PackageScope`FindDummyIndices[term];
	(* filter for flavor indices *)
	repeatedInds = DeleteCases[repeatedInds, Index[_, Except[Alternatives @@ Keys @ flavors]]];
	(* determine dimensions of these indices *)
	repeatedInds = Table[
		{ind, flavors[Last[ind]][IndexDimension]}
		,
		{ind, repeatedInds}
	];
	(* perform explicit sums over flavor indices *)
	If[Length[repeatedInds] > 0, Sum[term, Evaluate[Sequence @@ repeatedInds]], term]
]

ExplicitFlavor[expr_] := Module[
	{rules},
	RemoveField[Global`\[GothicU]];
	DefineField[Global`\[GothicU],Fermion,
		Indices->{Global`SU3c[fund]},
		Charges->{Global`U1em[+2/3]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->2,"name"->"u"|>
	];
	DefineField[Global`\[GothicC],Fermion,
		Indices->{Global`SU3c[fund]},
		Charges->{Global`U1em[+2/3]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->4,"name"->"c"|>
	];

	DefineField[Global`\[GothicT],Fermion,
		Indices->{Global`SU3c[fund]},
		Charges->{Global`U1em[+2/3]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->6,"name"->"t"|>
	];
	DefineField[Global`\[GothicB],Fermion,
		Indices->{Global`SU3c[fund]},
		Charges->{Global`U1em[-1/3]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->5,"name"->"b"|>
	];

	DefineField[Global`\[GothicS],Fermion,
		Indices->{Global`SU3c[fund]},
		Charges->{Global`U1em[-1/3]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->3,"name"->"s"|>
	];
	RemoveField[Global`\[GothicD]];
	DefineField[Global`\[GothicD],Fermion,
		Indices->{Global`SU3c[fund]},
		Charges->{Global`U1em[-1/3]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->1,"name"->"d"|>
	];
	RemoveField[Global`\[GothicE]];
	DefineField[Global`\[GothicE],Fermion,
		Indices->{},
		Charges->{Global`U1em[-1]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->11,"name"->"e"|>
	];
	DefineField[Global`\[Mu],Fermion,
		Indices->{},
		Charges->{Global`U1em[-1]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->13|>
	];
	DefineField[Global`\[Tau],Fermion,
		Indices->{},
		Charges->{Global`U1em[-1]},
		Mass->Heavy,
		UFO$Options-><|"pdg"->15|>
	];
	DefineField[Global`\[Nu]e,Fermion,
		Indices->{},
		Mass->Heavy,
		Chiral->LeftHanded,
		UFO$Options-><|"pdg"->12|>
	];
	DefineField[Global`\[Nu]\[Mu],Fermion,
		Indices->{},
		Mass->Heavy,
		Chiral->LeftHanded,
		UFO$Options-><|"pdg"->14|>
	];
	DefineField[Global`\[Nu]\[Tau],Fermion,
		Indices->{},
		Mass->Heavy,
		Chiral->LeftHanded,
		UFO$Options-><|"pdg"->16|>
	];
	rules = {
	Field[Global`\[GothicE],_,{n_},lorentz_] :>CD[lorentz,{Global`\[GothicE],Global`\[Mu],Global`\[Tau]}[[n]][]],
	Field[Global`\[GothicU],_,{Index[a_, Global`SU3c[_]],n_Integer},lorentz_] :>CD[lorentz,{Global`\[GothicU],Global`\[GothicC],Global`\[GothicT]}[[n]][a]],
	Field[Global`\[GothicD],_,{Index[a_, Global`SU3c[_]],n_Integer},lorentz_] :>CD[lorentz,{Global`\[GothicD],Global`\[GothicS],Global`\[GothicB]}[[n]][a]],
	Field[Global`\[Nu],_,{n_},lorentz_] :>CD[lorentz,{Global`\[Nu]e,Global`\[Nu]\[Mu],Global`\[Nu]\[Tau]}[[n]][]]
	};

	ExplicitEinsteinSums[expr] //. rules
]


(* ::Subsection::Closed:: *)
(*DefaultParamcard*)


DefaultParamcard[L_] := Module[
	{
	params,
	paramAssoc,
	fullJSON
	},
	
	params = Union[
		Flatten[
			TranslateToUFO /@ Cases[
				DeleteCases[L, c_Coupling|f_Field|fs_FieldStrength|i_Index, Infinity],
				s_Symbol /; Context[s] === "Global`",
				Infinity
			]
		],
		TranslateToUFO /@ Cases[L, c_Coupling, Infinity]
	];
	
	paramAssoc = Association @@ Table[param -> <|"name" -> param, "nature" -> "internal", "type" -> "", "value" -> "", "texname" -> ""|>, {param, params}];
	fullJSON = <|"description" -> "Minimal JSON with internal parameters only", "modelname" -> "Standard Model", "parameters" -> paramAssoc|>;
	Export["model_parameters.json", fullJSON, "JSON", "Indentation" -> 2]
]


(* ::Subsection:: *)
(*Coupling Orders*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineCouplingOrder::name = "`1` is already defined";
DefineCouplingOrder::couplings = "`1` are not defined couplings.";


(* ::Subsubsection::Closed:: *)
(*Define/Remove/Get/ResetCouplingOrders*)


$CouplingOrderAssociation = <||>;
GetCouplingOrders[names___] := Return @ $CouplingOrderAssociation[names];


Options[DefineCouplingOrder] = {PerturbativeExps -> None};

DefineCouplingOrder[name_, couplings_, expansionOrder_, hierarchy_, opts:OptionsPattern[]] := Module[
	{
	couplingNames = ExtractCouplingLabel[#] & /@ If[Head @ couplings =!= List, List @ couplings, couplings],
	couplingCount
	},
	
	If[
		Defined[name],
		Message[DefineCouplingOrder::name, name];
		Abort[]
	];
	
	If[
		Head[couplingNames] =!= List,
		couplingNames = List @ couplingNames
	];
	
	If[
		!AllTrue[couplingNames, KeyExistsQ[GetCouplings[], #] & ],
		Message[DefineCouplingOrder::couplings, StringRiffle[ToString /@ couplingNames, ", "]];
		Abort[]
	];
	
	couplingCount = AssociationThread[couplingNames, MapThread[PowerOf, {If[Head @ couplings =!= List, List @ couplings, couplings] //. Coupling[label_, ___] :> label, couplingNames}]];
	
	(* Add to global coupling orders*)
	AppendTo[
		$CouplingOrderAssociation,
		name -> <|
			Name -> name,
			Couplings -> couplingCount,
			ExpansionOrder -> expansionOrder,
			Hierarchy -> hierarchy,
			PerturbativeExps -> OptionValue @ PerturbativeExps
		|>
	];
	
	name::usage = "Coupling order " <> ToString @ name <> " for the couplings " <> StringRiffle[ToString /@ couplingNames, ", "];
];

RemoveCouplingOrder[name_] := Module[
	{},
	
	If[
		KeyExistsQ[$CouplingOrderAssociation, name],
		KeyDropFrom[$CouplingOrderAssociation, name]
	];
	
	ClearAll @ name;
];

ResetCouplingOrders[] := (RemoveCouplingOrder /@ Keys[$CouplingOrderAssociation];);


(* ::Subsubsection::Closed:: *)
(*ClassifyCouplingOrder*)


ClassifyCouplingOrder[term_] := Module[
	{
	couplingOrders = GetCouplingOrders[]
	},
	
	AssociationMap[
		Function[
			order,
			Total[
				GetCouplingOrders[order][Couplings][#]ExtractCouplingPower[term, #] & /@ Keys[GetCouplingOrders[order][Couplings]]
			]
		],
		Keys @ couplingOrders]
]


(* ::Subsection::Closed:: *)
(*RenameDummyIndicesToUFO*)


RenameDummyIndicesToUFO[term_] := Module[
	{
	dummies = FindDummyIndices[term],
	replacements,
	counter = -1
	},
	
	replacements = Association[
		Table[
			sym /. Index[sym_, grp_] :> Index[sym, grp] -> Index[ToString[counter--], grp],
			{sym, dummies}
		]
	];
	
	term /. Index[sym_, grp_] :> Lookup[replacements, Index[sym, grp], Index[sym, grp]]
];
