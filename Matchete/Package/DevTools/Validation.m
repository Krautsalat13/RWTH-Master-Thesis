(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Validation`*)


(* ::Text:: *)
(*Contains the functions needed for testing and validating the package*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["ValidateCurrentVersion"]
PackageScope["UpdateValidation"]
PackageScope["ActivateValidationMode"]
PackageScope["RunUnitTests"]
PackageScope["SaveValidationResults"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


ActivateValidationMode::usage= "SaveForValidation[\"model name\", bool] specifies that the next time the Match routine is called, the results are saved to an internal database for comparison. The String argument \"model name\" should match the shorthand name used for the model. The boolean argument \"reset\" specifies whether the derived results should replace the previous results or be treated as the new reults.";
UpdateValidation::usage= "Saves matching results for all model files to an internal directory for later comparison.";
ValidateCurrentVersion::usage= "Compares the current version of the code for a preset list of models to results obtainde on an earlier version of the code. Use option \"Tests\"-> <file name(s)> to specify which test files to use (default `All`). Use option \"Detail\"-> True to see the individual results of all tests.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Validation*)


(* ::Subsubsection::Closed:: *)
(*Functionality for the validation of the code*)


(* If True, the results of Match that are saved replace the previous results used for comparisson *)
$UpdateValidationResults= False;


(* ::Text:: *)
(*Function to activate the validation mode*)


ActivateValidationMode[str_String, reset_:False]:=Module[{},
	$ValidationRun= True;
	$UpdateValidationResults= reset;
	$ValidationModelName= str;
];


(* ::Text:: *)
(*Function calling the validation module*)


Options@ ValidateCurrentVersion= {
		"Tests"   -> All,
		"Details" -> False,
		"Models"  -> {"VLF_toy_model", "Singlet_Scalar_Extension", "E_VLL", "S1S3LQs"}
	};


ValidateCurrentVersion[opt:OptionsPattern[]]:= Module[{},
	ResetAll[];
	CheckDocumentation[];
	RunUnitTests["Tests"->OptionValue["Tests"],"Details"->OptionValue["Details"]];
	(*RunUnitTests@ opt;*)
	(* choose models to validate *)
	BranchValidation`$UVmodels= OptionValue["Models"];
	Get@ FileNameJoin[{$MatchetePath, "Validation", "Validation.m"}];
]


(* ::Text:: *)
(*Function to update the matching results used for the validation*)


Options@ UpdateValidation= {
		"Models" -> {"VLF_toy_model", "Singlet_Scalar_Extension", "E_VLL", "S1S3LQs"}
	};


UpdateValidation[OptionsPattern[]]:= Module[{},
	(* choose models to update *)
	UpdateBranchValidation`$UVmodels= OptionValue["Models"];
	Get@ FileNameJoin[{$MatchetePath, "Validation", "UpdateValidation.m"}]
]


(* ::Subsubsection::Closed:: *)
(*Routine to save the matching result for later comparison*)


(* ::Text:: *)
(*Function that writes the matching results to an internal database*)


SaveValidationResults[strResults_, LagrangianEFT_, time_, lag_]:= Module[{traceResults=strResults, Loff, tGreensSimplify, Lon, tEOMSimplify, tMapEffectiveCouplings=None, matchingCond=None,\[ScriptCapitalL]SMEFT},
		(* transform results of individual traces to a rule *)
		traceResults= Apply[Rule, Partition[traceResults,2], {1}];
		(* simplify all traces *)
		traceResults= traceResults /. {(Rule[a_,b_]:>Rule[a,b//ContractCGs//MatchReduce//GreensSimplify])};
		traceResults= Association@@ traceResults;

		(* simplify Lagrangian and save the computation time *)
		{tGreensSimplify, Loff} = Timing@GreensSimplify[LagrangianEFT];

		{tEOMSimplify, Lon}     = Timing@EOMSimplify[LagrangianEFT,DummyCoefficients->False, ReductionIdentities->dDimensional];
		
		(* check MapEffectiveCouplings only for SMEFT like theories*)
		If[!StringMatchQ[$ValidationModelName,"VLF_toy_model"],
			(* need to change context so SMEFT definitions match model definitions *)
			Begin["Global`"];
			\[ScriptCapitalL]SMEFT=LoadModel["SMEFT"];
			End[];
			
			{tMapEffectiveCouplings, matchingCond} = Timing@MapEffectiveCouplings[
				ReplaceEffectiveCouplings@EOMSimplify[LagrangianEFT,ReductionIdentities->EvanescenceFree], 
				ReplaceEffectiveCouplings@EOMSimplify[\[ScriptCapitalL]SMEFT,ReductionIdentities->FourDimensional]
				,
				SortByEFTOrder               -> True,
				KeepTrivalReplacements       -> True,
				AppendEffectiveCouplingsDefs -> True,
				EOMSimplify                  -> False,
				ShiftRenCouplings            -> False,
				ReductionIdentities          -> dDimensional
			];
		];

		(* save all relevant information to an internal directory *)
		SaveForComparison[<|
			"Model"                        -> $ValidationModelName,
			"Version"                      -> Global`$MatcheteVersion,
			"Date"                         -> Today,
			"Time (Match)"                 -> time,
			"Time (GreensSimplify)"        -> tGreensSimplify,
			"Time (EOMSimplify)"           -> tEOMSimplify,
			"Time (MapEffectiveCouplings)" -> tMapEffectiveCouplings,
			"UV Lagrangian"                -> lag,
			"Off-shell EFT Lagrangian"     -> Loff,
			"On-shell EFT Lagrangian"      -> ReplaceEffectiveCouplings[Lon],
			"SuperTraces"                  -> traceResults,
			"Matching Conditions"          -> matchingCond
		|>];

		(* end the validation run mode *)
		$UpdateValidationResults= False;
		$ValidationRun= False;
]


SaveForComparison[expr_Association]:=Module[
	{aux= If[$UpdateValidationResults, "previous", "current"]},
	Export[FileNameJoin@{$MatchetePath, "Validation", "MatchingResults", aux, ToString[expr["Model"]]<>"-EFT.m"}, expr];
];


(* ::Subsubsection::Closed:: *)
(*Run unit tests*)


Options@ RunUnitTests= {
		"Tests"-> All,
		"Details"-> False
	};


RunUnitTests[OptionsPattern[]]:= Module[{results, summary, test, testDi,  testFiles, testNames, testTally},
	testDir= FileNameJoin@{$MatchetePath, "Validation", "Tests"};
	testFiles= FileNames[{testDir <> "/*.wl", testDir <> "/*.m"}];
	(*Filter files*)
	If[(testNames= OptionValue@ "Tests") =!= All,
		testFiles= Cases[testFiles, str_String/; StringContainsQ[str, testNames]];
	];
	Monitor[
		results= Association@@ Table[
			FileNameTake@ file-> TestReport@ file
		, {file, testFiles}];
	, StringForm["Evaluating tests \"`1`\"...", FileNameTake@ file]];
	
	Print["____________________________"];
	If[OptionValue@ "Details",
		KeyValueMap[PrintTestResults[#1,#2] &, results]
	,
		testTally= 0;
		summary= Table[
			summary= SummarizeTestReport@ results@ test;
			If[Length@ summary[[2]] === 0, 
				testTally += summary[[1]];
				Nothing
			,
				test-> summary[[2]]
			]
		, {test, Keys@ results}];
		
		If[Length@ summary === 0, 
			Print[Style["\[CheckmarkedBox] ", Bold, Darker@ Green], Style["All "<> ToString@ testTally <> " unit tests succeeded!", Bold]];
			Print["From the file(s) ", Sequence@@ Flatten@ Table[{Style[test, "Code"], ", "} ,{test, Keys@ results}]];
		,
			Print[Style["\[WarningSign] " <> ToString@ Length@ Flatten@ summary <> " test(s) failed!", Bold, Darker@ Red]];
			summary= Association@@ summary;
			KeyValueMap[Print[Style[#1, "Code"], ":\t",  TableForm@ #2]&, summary] 
		];
	];
	Print["____________________________"];
	
]


SummarizeTestReport[testReport_TestReportObject]:= Module[{failures},
	failures= Table[
			If[test@ "Outcome" === "Success", 
				Nothing
			,
				test
			]
		, {test, testReport["Results"]}];
	{Length@ testReport["Results"], failures}
];


GetTestResults[tr_TestReportObject]:= Module[{abbreviations, fields, results},
	fields= {"TestID", "AbsoluteTimeUsed", "MemoryUsed", "Outcome"};
	abbreviations= {"AbsoluteTimeUsed"-> "Time [s]"};
	results= ReplaceRepeated[Outer[#1[#2]&, Values[tr["TestResults"]], fields],
		{x_Quantity:> QuantityMagnitude[x], x_Real:> Round[x, 0.001]}];
	Join[{fields/. abbreviations}, results]
]


PrintTestResults[name_String,testReport_TestReportObject]:= Module[{list, indx, noTests, time}, 
	list= GetTestResults@ testReport;
	indx= MapIndexed[If[MemberQ[#1, "Failure"|"MessagesFailure"|"Error"], First[#2], Nothing]&, list];
	time= Round[QuantityMagnitude[testReport["TimeElapsed"]], 0.01];
	noTests= Length[testReport["TestResults"]];
	If[TrueQ@ testReport["AllTestsSucceeded"], 
		Print[Style[name, "Code"], ":\t", Style["All tests succeeded! \[CheckmarkedBox]", Bold, Darker@ Green]], 
		Print[Style[name, "Code"], ":\t", 
			Style[ToString@ testReport["TestsFailedCount"] <> " tests failed! \[WarningSign]", Darker@ Red, Bold]];
	];
	Print@ Grid[list, Alignment->Left, Dividers-> {None, {2-> True}}, Background-> {None, Thread[indx->Pink]}];
]


(* ::Subsubsection::Closed:: *)
(*Check documentation*)


CheckDocumentation[]:=Module[{deprecatedSymbs,docuFiles,exportedSymb,missingDocu},
	deprecatedSymbs= {"DefineGroup", "DefineGroupRepresentation"};
	docuFiles= FileNameTake/@ FileNames@FileNameJoin[{$MatchetePath,"DocumentationSource","Documentation","English","ReferencePages","Symbols","*.nb"}];
	docuFiles= FileBaseName/@ docuFiles;
	exportedSymb= Names["Matchete`*"];
	exportedSymb= exportedSymb /. {
		"\[ScriptD]"->"ScriptD",
		"\[Gamma]"->"Gamma",
		"\[Epsilon]"->"Epsilon",
		"\[Mu]bar2"->"Mubar2",
		"\[Sigma]"->"Sigma"
	};
	missingDocu= Complement[exportedSymb,docuFiles,deprecatedSymbs];
	If[Length[missingDocu]>0,
		Print@Style["Warning:",Orange,Bold];
		Print["Missing documentation notebooks for the following exported symbols:"];
		Do[
			Print@Style[ToString[doc],Bold],
			{doc,missingDocu}
		];
		Print["Please add the missing documentation notebooks to the directory "<>FileNameJoin[{$MatchetePath,"DocumentationSouce"}]<>", and then rebuilt the documentation."];
	];
]
