(* ::Package:: *)

Begin["UpdateBranchValidation`"]


(* ::Title:: *)
(*Validation of Matching Results*)


(* ::Text:: *)
(*List of UV models for which a model file and a saved EFT Lagrangian result exist.*)


(* this is set in Matchete.m now *)
(*$UVmodels= {"VLF_toy_model", "Singlet_Scalar_Extension", "E_VLL", "S1S3LQs"};*)


If[ChoiceDialog["Do you really want to overwrite the matching results that are used for validation?", {"Yes"->True, "No"->False}],
	Do[
		(* reset the model *)
		ResetAll[];
		
		(* turn on validation mode *)
		Matchete`PackageScope`ActivateValidationMode[model, True];
		
		(* define and run the model *)
		Begin["Global`"];
			LUV=LoadModel[model];
		End[];
		
		Match[LUV,EFTOrder->6,LoopOrder->1];
	,
		{model, $UVmodels}
	]
]


End[];
