(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Counterterms`*)


(* ::Subtitle:: *)
(*Tools needed to calculate finite and divergent counterterms*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsection:: *)
(*Exported*)


PackageExport["OnShell"]


(* ::Subsection:: *)
(*Internal*)


PackageScope["UVDivergentAction"]


PackageScope["EvanescentContribution"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection:: *)
(*Exported*)


OnShell::usage = "Option for UVDivergentAction. Can be set to: Automatic, True, or False";


(* ::Subsection:: *)
(*Internal*)


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Divergence*)


(* ::Subsection:: *)
(*UVDivergentAction*)


Options@ UVDivergentAction = { EFTOrder -> Automatic , Verbose -> True, OnShell -> Automatic };


UVDivergentAction::orderArgument = "The option EFTOrder-> `1`  is neither an integer nor Automatic and cannot be used as an order specification.";


UVDivergentAction[expr_,OptionsPattern[]]:=Module[{powerTraces,logTraces,isOnShell=True,L=expr,order,
								ptr,LDiv,field,i},

	(* determine order argument and cut the Lagrangian to the given order *)
	{order, L}= TruncateLagrangian[OptionValue@ EFTOrder, L];
	
	(*If[IntegerQ[OptionValue@ EFTOrder],
		order = OptionValue@ EFTOrder;
		L = SeriesEFT[L,EFTOrder->order]
	,
		If[SameQ[OptionValue@ EFTOrder,Automatic],
			order = Max[OperatorDimension/@TermsToList[Expand@L]];
			L = SeriesEFT[L, EFTOrder -> order ]
		,
			Message[UVDivergentAction::orderArgument,OptionValue@ EFTOrder];Return[$Aborted]
		]
	];*)
	
	(* check if this Lagrangian contains redundant operators - if so, this Lagrangian is not in the on-shell basis *)
	If[Length @ FieldsToShift @ InternalSimplify @ L > 0, isOnShell = False];
	If[OptionValue@OnShell===False, isOnShell = False];
	
	(* set the expanded Lagrangian and determine X terms *)
	SetCurrentLagrangian[L, 1, order, Mode -> Divergence];
	 
	$MonitorString="";
	
	(* determine the Log traces *)
	logTraces = OptionalMonitor[TrueQ@ OptionValue@ Verbose,
					Sum[
						If[Or[
							Length@ GetFieldsByProperty[Sequence@@ Normal@ $FieldTypes@ field, Charges-> {__}] > 0,
							Length@ GetFieldsByProperty[Sequence@@ Normal@ $FieldTypes@ field,
								Indices-> inds_/; IntersectingQ[GroupFromRep/@ inds, Keys@ $GaugeGroups]] > 0
						],
							$MonitorString = StringForm["Evaluating log-type supertrace: `1`", field/. fieldFormat];
							LogTypeSTr[field, order, Mode -> Divergence]
						,
							0
						]
					,{field,{hScalar,hFermion,hVector,hGhost,lScalar,lFermion,lVector,lGhost}}
					],
					$MonitorString					
				];
	
	(* list the possible power traces, including light fields *)
	ptr = ListPowerTypeTraces[ order, True ];
	
	powerTraces = OptionalMonitor[TrueQ@ OptionValue@ Verbose,
						Sum[
							$MonitorString = StringForm["Evaluating power-type supertrace: `1` \t (`2` / `3`)",
											ptr[[i]] /. fieldFormat, i, Length@ ptr];
							PowerTypeSTr[ptr[[i]],order,Mode->Divergence],
						{i, Length@ptr}]
					,
						$MonitorString];
	
	(* the Dirac algebra induces terms with \[ScriptD] so expand once more and drop vacuum terms *)
	(* the Collect calls before the Series calls below bring significant performance improvements *)
	(*LDiv = Normal@ Series[Collect[Select[ BetterExpand@ (logTraces +powerTraces ), !FreeQ[#,Field|FieldStrength]&],\[Epsilon]],{\[Epsilon],0,-1}];*)
	
	(* ExtractePoles can be significantly faster than Series *)
	LDiv = Select[BetterExpand[logTraces + powerTraces], !FreeQ[#,Field|FieldStrength]&];
	LDiv = ExtractPoles[LDiv];
	
	(*Return[LDiv];*)

	If[isOnShell  ,
		(*Normal@ Series[Collect[EOMSimplify[L + LDiv,DummyCoefficients->True],\[Epsilon]],{\[Epsilon],0,-1}],*)
		ExtractPoles[EOMSimplify[L + LDiv,DummyCoefficients->True]], (* faster than Series *)
		GreensSimplify@ LDiv
	]
]


TruncateLagrangian[eftOrderOption_, lag_]:= Module[{order}, 
	Switch[eftOrderOption
	,_Integer,
		{eftOrderOption, SeriesEFT[lag, EFTOrder-> eftOrderOption]}
	,Automatic,
		{order= Max[OperatorDimension/@ TermsToList@ lag], SeriesEFT[lag, EFTOrder-> order]}
	,_,
		Message[UVDivergentAction::orderArgument, eftOrderOption]; Abort[];
	]
];


(* function to extractt 1/\[Epsilon] poles in an efficient way *)
ExtractPoles[expr_]:=Module[{x},
	Coefficient[expr/.{\[Epsilon]->\[Epsilon]/x}, x]
]


(* ::Section:: *)
(*Evanescent contribution*)


(* ::Subsection:: *)
(*EvanescentContribution*)


Options@ EvanescentContribution = {EFTOrder-> Automatic, Verbose-> True};


EvanescentContribution::orderArgument = "The option EFTOrder-> `1`  is neither an integer nor Automatic and cannot be used as an order specification.";


EvanescentContribution[lagrangian_, OptionsPattern[]]? OptionsCheck:= Module[
		{lag, order, out, traceTypes, fields, i= 0},
		
	(* Determine order argument and cut the Lagrangian to the given order *)
	{order, lag}= TruncateLagrangian[OptionValue@ EFTOrder, lagrangian];
	
	(* Set the expanded Lagrangian and determine X terms *)
	SetCurrentLagrangian[lag, 1, order, Mode-> Evanescent, Verbose-> OptionValue@ Verbose];
	
	(* Loop-level evanescent insertions *)
	traceTypes= ListEvanescentTypeTraces[order];
	out= OptionalMonitor[OptionValue@ Verbose,
			Sum[i++;
				PowerTypeSTr[fields, order, Mode-> Evanescent]
			, {fields, traceTypes}]
		, StringForm["Evaluating power-type supertrace: `1` \t (`2` / `3`)",
			fields/. fieldFormat, i, Length@ traceTypes] ];
	
	out// ContractCGs// MatchReduce
]


(* ::Subsubsection::Closed:: *)
(*Listing evanescent trace types*)


(* ::Text:: *)
(*ListEvanescentTypeTraces provides the list of all power type traces beginning the reading from an evanescent operator, while being indifferent as to the nature of the remaining insertions*)


ListEvanescentTypeTraces@ {order_Integer}:= ListEvanescentTypeTraces@ order;


ListEvanescentTypeTraces@ order_Integer:= Module[{possibilities, seed, temp},
	(*Recursively determines candidates for traces with order \[LessEqual] order*)
	possibilities= Flatten[Last@ Reap[
		seed= Sow@ {{lScalar}, {lFermion}, {lVector}};
		While[(Length@ seed> 0),
			temp= Flatten/@ Tuples[{seed, {lScalar, lFermion, lVector}}];
			seed= Sow@ Select[temp, $EvaOrdMin@ #[[;;2]]+ Total[BlockMap[$XOrdMin, #[[2;;]], 2, 1]] <= order&];
		];
	], 2];

	(*Eliminates traces with order > order*)
	possibilities= Select[possibilities,
		If[Length@# === 1, $EvaOrdMin@ #[[{-1, 1}]], 
			$EvaOrdMin@ #[[;;2]]+ Total@ BlockMap[$XOrdMin, #[[2;;]], 2, 1]+ $XOrdMin@ #[[{-1, 1}]]] <= order&];

	(*Delete duplicates under cyclic permutations*)
	DeleteDuplicatesBy[possibilities, (First@ Sort@ NestList[RotateLeft, #, Length@# - 1] &)]
];
