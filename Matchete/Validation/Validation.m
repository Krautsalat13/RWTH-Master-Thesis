(* ::Package:: *)

Begin["BranchValidation`"]


(* ::Title:: *)
(*Validation of Matching Results*)


(* ::Text:: *)
(*List of UV models for which a model file and a saved EFT Lagrangian result exist.*)


(* this is set in Matchete.m now *)
(*$UVmodels= {"VLF_toy_model", "Singlet_Scalar_Extension", "E_VLL", "S1S3LQs"};*)


$validatedModels=0;


Do[
	(* reset the model *)
	ResetAll[];
	
	(* delete temporary files *)	
	DeleteFile/@FileNames[FileNameJoin@{$MatchetePath, "Validation", "MatchingResults", "current","*"}];
	
	Print["Validating model: \"", Style[model, Bold],"\""];
	
	(* turn on validation mode *)
	Matchete`PackageScope`ActivateValidationMode[model, False];
	
	LUV::usage="UV Lagrangian";
	(* define and run the model *)
	Begin["Global`"];
		LUV=LoadModel[model];
	End[];
	
	Match[LUV,EFTOrder->6,LoopOrder->1];
	
	(* load current result *)
	LEFT$now= Get@FileNameJoin[{$MatchetePath, "Validation", "MatchingResults", "current", model<>"-EFT.m"}];
	
	(* load previous result *)
	LEFT$previous= Get@FileNameJoin[{$MatchetePath, "Validation", "MatchingResults", "previous", model<>"-EFT.m"}];
	
	Print["\t- Previous results determined with version ", LEFT$previous["Version"], " on ", LEFT$previous["Date"], "."];
	
	(* check that the UV Lagrangian actually agree *)
	If[GreensSimplify[LEFT$now["UV Lagrangian"]-LEFT$previous["UV Lagrangian"]]=!=0,
		Print[Style["The UV Lagrangian do not agree!", Red, Bold]];
		Echo[LEFT$now["UV Lagrangian"], "UV Lagrangian NOW:", Format[#,NiceForm]&];
		Echo[LEFT$previous["UV Lagrangian"], "UV Lagrangian BEFORE:", Format[#,NiceForm]&];
	];
	
	(* compare computation time *)
	Print["\t- Computation time in sec."];
	Print["\t\t\[Rule] ", Style["Match","Code"], ": \t\t\t\t\t  ", LEFT$now["Time (Match)"], " (now) vs. ", LEFT$previous["Time (Match)"], " (before)"];
	Print["\t\t\[Rule] ", Style["GreensSimplify","Code"], ": \t\t\t ", LEFT$now["Time (GreensSimplify)"], " (now) vs. ", LEFT$previous["Time (GreensSimplify)"], " (before)"];
	Print["\t\t\[Rule] ", Style["EOMSimplify","Code"], ": \t\t\t\t", LEFT$now["Time (EOMSimplify)"], " (now) vs. ", LEFT$previous["Time (EOMSimplify)"], " (before)"];
	Print["\t\t\[Rule] ", Style["MapEffectiveCouplings","Code"], ": \t  ", LEFT$now["Time (MapEffectiveCouplings)"], " (now) vs. ", LEFT$previous["Time (MapEffectiveCouplings)"], " (before)"];
	
	(* compare off-shell Lagrangian *)
	diff$off= GreensSimplify@Contract@EvaluateLoopFunctions[LEFT$now["Off-shell EFT Lagrangian"]-LEFT$previous["Off-shell EFT Lagrangian"]];
	diff$off= If[Head[diff$off]===Plus,
		FullSimplify[#,And@@DeleteDuplicates@Cases[#,Log[arg_]:>(arg>0),All]]&/@diff$off,
		FullSimplify[diff$off,And@@DeleteDuplicates@Cases[diff$off,Log[arg_]:>(arg>0),All]]
	];
	If[diff$off===0,
		Print["\t- ", Style["\[CheckmarkedBox]", Darker@Green], " Off-shell EFT Lagrangian for the model agrees with the previous result."]
	,
		Print["\t- ", Style["\[WarningSign]", Red]," Off-shell EFT Lagrangian does not agrees with the previous result."];
		Echo[diff$off, "Difference off-shell EFT Lagrangian: ", Iconize[#,Format[#,NiceForm]]&];
	];
	
	(* check individual loops *)
	If[Sort@Keys@LEFT$now["SuperTraces"]=!=Sort@Keys@LEFT$previous["SuperTraces"],
		Print@Style["\t- The number of SuperTraces differs.", Red, Bold];
	,
		STrCounter=0;
		Do[
			\[CapitalDelta]STr=GreensSimplify@Contract@EvaluateLoopFunctions[LEFT$now["SuperTraces"][trace]-LEFT$previous["SuperTraces"][trace]];
			\[CapitalDelta]STr= If[Head[\[CapitalDelta]STr]===Plus,
				FullSimplify[#,And@@DeleteDuplicates@Cases[#,Log[arg_]:>(arg>0),All]]&/@\[CapitalDelta]STr,
				FullSimplify[\[CapitalDelta]STr,And@@DeleteDuplicates@Cases[\[CapitalDelta]STr,Log[arg_]:>(arg>0),All]]
			];
			If[\[CapitalDelta]STr=!=0,
				Print["\t- ", Style["\[WarningSign]", Red]," The SuperTrace ", trace, "does not agree with the previous result."];
				Echo[\[CapitalDelta]STr, "\[CapitalDelta]STr("<>trace<>"): ", Iconize[#,Format[#,NiceForm]]&];
				STrCounter++;
			]
		,
			{trace, Keys@LEFT$now["SuperTraces"]}
		];
		Print["\t- ", If[STrCounter==0,Style["\[CheckmarkedBox]", Darker@Green],Style["\[WarningSign]", Red]]," Agreement found for ", Length@Keys@LEFT$now["SuperTraces"]-STrCounter, "/", Length@Keys@LEFT$now["SuperTraces"], " SuperTraces."];
	];
	
	(* compare on-shell Lagrangian *)
	diff$on= GreensSimplify@Contract@EvaluateLoopFunctions[LEFT$now["On-shell EFT Lagrangian"]-LEFT$previous["On-shell EFT Lagrangian"]];
	diff$on= If[Head[diff$on]===Plus,
		FullSimplify[#,And@@DeleteDuplicates@Cases[#,Log[arg_]:>(arg>0),All]]&/@diff$on,
		FullSimplify[diff$on,And@@DeleteDuplicates@Cases[diff$on,Log[arg_]:>(arg>0),All]]
	];
	If[diff$on===0,
		Print["\t- ", Style["\[CheckmarkedBox]", Darker@Green], " On-shell EFT Lagrangian agrees with the previous result."];
		$validatedModels++;
	,
		Print["\t- ", Style["\[WarningSign]", Red]," On-shell EFT Lagrangian does not agrees with the previous result."];
		Echo[diff$on, "Difference on-shell EFT Lagrangian: ", Iconize[#,Format[#,NiceForm]]&];
	];
	
	diff$shift= GreensSimplify[ReplaceEffectiveCouplings@LEFT$now["Off-shell EFT Lagrangian"]-ReplaceEffectiveCouplings@ShiftRenCouplings@ReplaceEffectiveCouplings@LEFT$now["Off-shell EFT Lagrangian"], ReductionIdentities->dDimensional];
	If[0 === diff$shift, 
		Print["\t- ", Style["\[CheckmarkedBox]", Darker@Green], " ReplaceEffectiveCouplings@ShiftRenCouplings[\!\(\*SubscriptBox[\(\[ScriptCapitalL]\), \(off - shell\)]\)] \[Equal] \!\(\*SubscriptBox[\(\[ScriptCapitalL]\), \(off - shell\)]\)."];
		,
		Print["\t- ", Style["\[WarningSign]", Red]," ReplaceEffectiveCouplings@ShiftRenCouplings[\!\(\*SubscriptBox[\(\[ScriptCapitalL]\), \(off - shell\)]\)] \[NotEqual] \!\(\*SubscriptBox[\(\[ScriptCapitalL]\), \(off - shell\)]\)."];
		Echo[diff$shift, "\[ScriptCapitalL] - ReplaceEffectiveCouplings[ShiftRenCouplings[\[ScriptCapitalL]]]: ", Iconize[#,Format[#,NiceForm]]&];
	];
	
	(* check MapEffectiveCouplings *)
	If[!StringMatchQ[model,"VLF_toy_model"],
		\[ScriptCapitalL]SMEFT = LoadModel["SMEFT"];
		mc$diff = GreensSimplify[ReplaceInLagrangian[\[ScriptCapitalL]SMEFT, LEFT$now["Matching Conditions"]]-
			ReplaceInLagrangian[\[ScriptCapitalL]SMEFT, LEFT$previous["Matching Conditions"]], ReductionIdentities->dDimensional];
		If[mc$diff===0,
			Print["\t- ", Style["\[CheckmarkedBox]", Darker@Green], " Warsaw basis matching conditions agrees with the previous result."];
			(*$validatedModels++;*)
		,
			Print["\t- ", Style["\[WarningSign]", Red]," Warsaw basis matching conditions do not agrees with the previous result."];
			Echo[mc$diff, "Difference of Warsaw basis EFT Lagrangians: ", Iconize[#,Format[#,NiceForm]]&];
		];
	];
	
	Print["__________"];
,
	{model, $UVmodels}
]


Print["Validation of implemented models finished: ", $validatedModels, "/", Length@$UVmodels, " successfully verfied."];


End[];
