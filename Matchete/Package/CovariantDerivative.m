(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`CovariantDerivative`*)


(* ::Subtitle:: *)
(*Covariant derivative rules for Matchete paclet.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsection:: *)
(*Exported*)


PackageExport["CD"]


(* ::Subsection:: *)
(*Internal*)


PackageScope["OpenCD"]


PackageScope["TakeDev"]


PackageScope["CommuteCDs"]
PackageScope["GAction"]
PackageScope["SeparateGeneratorsFromFS"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection:: *)
(*Exported*)


CD::usage     = "CD[ind,expr] returns the covariant derivative(s) of a given expression expr. The argument ind can either be a single Lorentz index or a list of Lorentz indices.";


(* ::Subsection:: *)
(*Internal*)


OpenCD::usage = "OpenCD[indices] denotes open covariant derivatives that act on everything to the right of the given expression. The argument indices is must be a list of Lorentz indices corresponding to the indices of the covariant derivatives.";


TakeDev::usage = "TakeDev[{Index[\[Mu],Lorentz],Index[\[Nu],Lorentz],...}, expr] is an internal function encoding the covariant derivative rules.";


CommuteCDs::usage               = "CommuteCDs[field,n] commutes derivatives n and n+1 (from the left) of a Field object D_{\[Mu]1, ... \[Mu]n, \[Mu](n+1), ...} f.";
GAction::usage                  = "GAction[{Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]},field] denotes the action of a field strength tensor on a field or FS object.";
SeparateGeneratorsFromFS::usage = "Matchete uses implicit contraction of FS with generators at certain stages in the evaluation. This functions separates F_{mu nu}^a_b = F_{mu nu}^A t^{Aa}_b, with t^{Aa}_b being the appropriate generator.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Covariant derivative definition and rules*)


(* ::Subsection:: *)
(*Error messages*)


CD::indices = "The index or indices '`1`' are invalid. These must be either a symbol or a list of symbols.";


(* ::Subsection:: *)
(*CD definition*)


(* ::Text:: *)
(*Covariant derivative definition*)


CD[ind_, expr_]:= Module[{},
	(* Check that \[Mu] has the correct form *)
	If[
		!MatchQ[ind, _Symbol | {_Symbol..} | Index[_,Lorentz] | {Index[_,Lorentz]..} | {} | _Pattern],
		Message[CD::indices,ind];
		Abort[]
	];
	
	Switch[ind,
		Index[_,Lorentz],     TakeDev[{ind},expr],
		{Index[_,Lorentz]..}, TakeDev[ind,expr],
		_Symbol,              TakeDev[Index[#,Lorentz]&/@{ind},expr],
		{_Symbol..},          TakeDev[Index[#,Lorentz]&/@ind,expr],
		{},                   expr
	]
]


(* ::Subsection:: *)
(*CD Test function*)


(* ::Text:: *)
(*This tells when the covariant derivative is zero*)


DZeroQ[f_? DZeroQ[x__]]:=And@@ DZeroQ/@ {x};
DZeroQ[f_[x__]]:= False;
DZeroQ[_]= True;
DZeroQ[Field]= False;
DZeroQ[FieldStrength]= False;
DZeroQ[OpenCD]= False;
DZeroQ[WilsonTerm]= False;


(* ::Subsection:: *)
(*Applying covariant derivative to an expression*)


(* ::Text:: *)
(*Basic rules for the covariant derivative*)


TakeDev[{}, expr_]:= expr;
TakeDev[{\[Mu]_, \[Nu]__}, expr_]:= TakeDev[{\[Mu]}, TakeDev[{\[Nu]}, expr]];
TakeDev[\[Mu]_, _?DZeroQ]= 0;
TakeDev[\[Mu]_, expr:Alternatives[_Plus, _List]]:= TakeDev[\[Mu], #]&/@ expr;
TakeDev[\[Mu]_, expr:Alternatives[_Times, _NonCommutativeMultiply, _NCProduct, _FuncNCM]]:=
	Module[{n}, Sum[MapAt[TakeDev[\[Mu], #]&, expr, n], {n, Length@expr}] ];
TakeDev[{\[Mu]_}, Power[expr_, n_]]:= n TakeDev[{\[Mu]}, expr] Power[expr, n-1];


TakeDev[\[Mu]_,Field[label_,type_,indices_,\[Nu]_]]:= Field[label,type,indices,Join[\[Mu],\[Nu]]];
TakeDev[\[Mu]_,Bar@Field[label_,type_,indices_,\[Nu]_]]:= Bar@Field[label,type,indices,Join[\[Mu],\[Nu]]];
TakeDev[\[Mu]_,Transp@Field[label_,type_,indices_,\[Nu]_]]:= Transp@Field[label,type,indices,Join[\[Mu],\[Nu]]];
TakeDev[\[Mu]_,Transp@Bar@Field[label_,type_,indices_,\[Nu]_]]:= Transp@Bar@Field[label,type,indices,Join[\[Mu],\[Nu]]];
TakeDev[\[Mu]_,FieldStrength[label_,li_,indices_,\[Nu]_]]:= FieldStrength[label,li,indices,Join[\[Mu],\[Nu]]];
TakeDev[\[Mu]_,Bar@FieldStrength[label_,li_,indices_,\[Nu]_]]:= Bar@FieldStrength[label,li,indices,Join[\[Mu],\[Nu]]];


TakeDev[\[Mu]_,X_EoM]:= CD[\[Mu], OperatorToNormalForm@ X];
TakeDev[\[Mu]_, op_Operator]:= Operator@ CD[\[Mu], OperatorToNormalForm@ op];


TakeDev[\[Mu]_, WilsonTerm[a__, \[Nu]_]]:= WilsonTerm[a, Join[\[Mu], \[Nu]]]; 


(* derivatives acting to the right *)
TakeDev[\[Mu]_List,OpenCD[\[Nu]_List]]:= OpenCD[Join[\[Mu],\[Nu]]];


(* ::Subsection:: *)
(*Canonize derivatives *)


(* ::Subsubsection::Closed:: *)
(*Action of field-strength tensor on field*)


(* ::Text:: *)
(*The action of a field strength tensor on a field/FS object*)


GAction[lInds_List, field:(Field|FieldStrength)[lab_, _, inds_, devs_]]:= Module[{charge, charges, gaugeCharges, gaugeInds, group, 
		indc, out, ind, rep, j, A},	
	
	(*Abelian FS*)
	charges= GetFields[lab, Charges];
	gaugeCharges= Cases[charges, _? (MemberQ[Keys@$GaugeGroups,Head@ #] &)];
	out= Sum[
			{group, charge}= {Head@ charge, First@ charge};
			charge 
			(*$GaugeGroups[group, Coupling][] (*added*)*)
			FieldStrength[$GaugeGroups[group, Field], lInds, {}, {}]
		, {charge, gaugeCharges}]* field;
	
	(*Non-Abelian FS*)
	gaugeInds= Cases[inds, _? (MemberQ[Keys@$GaugeGroups,GroupFromInd@ #] &)];
	out+= Sum[
		group= GroupFromInd@ ind;
		rep= Last[ind/. Bar-> Identity];
		indc= If[Head@ ind === Bar, Index[j, rep], Bar@ Index[j, rep]];
		If[rep === group@ adj, 
			I CG[fStruct@ group, {ind, Index[A, group@ adj], indc}], 
			If[Head@ ind === Bar, (*Hermitian generators*)
				- CG[gen@ rep, {Index[A, group@ adj], indc, ind}],
				CG[gen@ rep, {Index[A, group@ adj], ind, indc}] 
			]]
		(*$GaugeGroups[group, Coupling][] (*added*)*)
		FieldStrength[$GaugeGroups[group, Field], lInds, {Index[A, group@ adj]}, {}] 
		(field/. ind-> Bar@ indc)
	, {ind, gaugeInds}];
	
	out 
];


(* ::Subsubsection::Closed:: *)
(*Commute derivatives *)


(* ::Text:: *)
(*Commuting derivatives n and n+1 (from the left) of a Field object  D_{\[Mu]1, ... \[Mu]n, \[Mu](n+1), ...} f.*)


CommuteCDs[Field[lab_, type_, inds_, devs_], n_]:= Module[{},
	Field[lab, type, inds, Join[devs[[;;n-1]], devs[[{n+1, n}]], devs[[n+2;;]] ]] + 
		-I TakeDev[devs[[;;n-1]], GAction[devs[[{n, n+1}]], 
			Field[lab, type, inds, devs[[n+2;;]] ]]]
 ]; 


CommuteCDs[FieldStrength[lab_, linds_, inds_, devs_], n_]:= Module[{},
	FieldStrength[lab, linds, inds, Join[devs[[;;n-1]], devs[[{n+1, n}]], devs[[n+2;;]] ]] + 
		-I TakeDev[devs[[;;n-1]], GAction[devs[[{n, n+1}]], 
			FieldStrength[lab, linds, inds, devs[[n+2;;]] ]]]
 ]; 


(* ::Subsection:: *)
(*Separating generators from FS*)


(* ::Text:: *)
(*Matchete uses implicit contraction of FS with generators at stages in the evaluation. This functions separates F_{mu nu}^a_b = F_{mu nu}^A t^{Aa}_b, with t^{Aa}_b being the appropriate generator*)


SeparateGeneratorsFromFS@ HoldPattern@ FieldStrength[pre__, {ind_, indc_}, devs_]:= Module[{gr, rep, conj, A},
	gr= GroupFromInd@ ind;
	rep= Last[ind/. Bar-> Identity];
	conj= Head@ ind === Bar;
	FieldStrength[pre, {Index[A, gr@ adj]}, devs] *
		If[rep === gr@ adj, 
			I CG[fStruct@ gr, {ind, Index[A, gr@ adj], indc}], 
			If[conj, 
				- CG[gen@ rep, {Index[A, gr@ adj], indc, ind}],
				CG[gen@ rep, {Index[A, gr@ adj], ind, indc}] 
			]]
];


SeparateGeneratorsFromFS@ fs:HoldPattern@ FieldStrength[pre__, {}|{_}, devs_]:= fs;


SeparateGeneratorsFromFS@ expr_:= expr/. fs_FieldStrength:> SeparateGeneratorsFromFS@ fs;
