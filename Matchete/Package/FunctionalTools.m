(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`FunctionalTools`*)


(* ::Subtitle:: *)
(*Sub package providing the functional tools used in functional one-loop matching.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["DeriveEOM"]
PackageScope["FluctuationOperator"]


PackageScope["VarD"]
PackageScope["FD"]
PackageScope["FuncD"]
PackageScope["FuncDSimplify"]
PackageScope["ShiftCD"]
PackageScope["ShiftFS"]


PackageScope["GaugeFixing"]


PackageScope["TransposeThisSpinChain"]


PackageScope["SubstituteFieldStrength"]
PackageScope["SubstituteHeavyFieldStrength"]


PackageScope["GaugeFieldQ"]


PackageScope["ExpandVectorFluctuations"]


PackageScope["ChangeFieldIndex"]


PackageScope["NCProduct"]


PackageScope["FieldTransformsUnderGaugeGroupQ"]


PackageScope["BackgroundField"]
PackageScope["BackgroundCD"]
PackageScope["BackgroundFS"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


FD::usage =
"FD[expr,field] calculates the (partial) functional derivative of expr with respect to field, where expr is a Lagrangian or any other expression depending on fields. The argument field must be of type Field[label, type, {indices}, {CDerivs}] and the diferent covariant derivatives of the field are treated as independent.
Example: FD[\[ScriptCapitalL], Field[\[Phi], Scalar, {}, {\[Mu],\[Nu]}]] = \!\(\*FractionBox[\(\[Delta]\[ScriptCapitalL]\), \(\[Delta] \((\*SubscriptBox[\(D\), \(\[Mu]\)] \*SubscriptBox[\(D\), \(\[Nu]\)] \[Phi])\)\)]\).";


FuncD::usage =
"FuncD[expr,field] calculates the (partial) functional derivative of expr with respect to field, where expr is a Lagrangian or any other expression depending on fields. Results are returned in a non simplified form. The argument field must be of type Field[label, type, {indices}, {CDerivs}] and the diferent covariant derivatives of the field are treated as independent.
Example: FuncD[\[ScriptCapitalL], Field[\[Phi], Scalar, {}, {\[Mu],\[Nu]}]] = \!\(\*FractionBox[\(\[Delta]\[ScriptCapitalL]\), \(\[Delta] \((\*SubscriptBox[\(D\), \(\[Mu]\)] \*SubscriptBox[\(D\), \(\[Nu]\)] \[Phi])\)\)]\).";


FuncDSimplify::usage =
"FuncDSimplify[expr] simplifies the expression expr that was returned by FuncD. It contracts all indices in Delta[_,_] and Metric[_,_] and transposes spin-chains if necessary.  ";


VarD::usage=
"VarD[expr,field1,field2(optional)] calculates the variational derivative of expr with respect to field2 and then field1, where expr is a Lagrangian or any other expression depending on fields. The arguments field1 and field2 must be of the form Field[label, type, {indices}, {}].";


Grassmann::usage =
"Grassmann is an option for FuncD and VarDraw that specifies whether minus signs from anti-commuting functional derivatives of fermionic fields should be considered. The default value is Grassmann->True.";


DeriveEOM::usage =
"DeriveEOM[\[ScriptCapitalL],\[Phi]] derives the equation of motion for the field \[Phi] from the Lagrangian density \[ScriptCapitalL], i.e. VarD[\[ScriptCapitalL],Bar[\[Phi]]] is computed.";


FluctuationOperator::usage =
"FluctuationOperator[\[ScriptCapitalL],\[Phi],\[Psi]] derives the fluctuation operator, i.e. the X-term, for the fields \[Phi] and \[Psi] from the Lagrangian density \[ScriptCapitalL], i.e. \!\(\*FractionBox[\(\*SuperscriptBox[\(\[Delta]\), \(2\)] \[ScriptCapitalL]\), \(\[Delta]\[Phi]\\\ \[Delta]\[Psi]\)]\) ~ VarD[\[ScriptCapitalL],\[Psi],\[Phi]] is computed.";


RecursiveExpandCD::usage =
"RecursiveExpandCD[indices,group,field] shifts all covariant derivatives by the gauge field fluctuation and expands the result up to leading order in the gauge field fluctuation.";


SubstituteFieldStrength::usage =
"SubstituteFieldStrength[expr] substitues all instances of FieldStrength[...] in expr with the corresponding Fields using \!\(\*SubscriptBox[SuperscriptBox[\(F\), \(a\)], \(\[Mu]\[Nu]\)]\)=(\!\(\*SubscriptBox[\(D\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(A\), \(\[Nu]\)]\)\!\(\*SuperscriptBox[\()\), \(a\)]\)-(\!\(\*SubscriptBox[\(D\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(A\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\()\), \(a\)]\)."


ChangeFieldIndex::usage="ChangeFieldIndex[field, group, indexLabel] returns the index of field that belongs to the given group and returns the field with that index replaced by the given indexLabel.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*FD = FuncDSimplify@FuncD (partial functional derivative)*)


(* ::Text:: *)
(*This function calls FuncD (which does the actual partial functional derivative) and applies some simplifications to it.*)


(* ::Subsection:: *)
(*Main entry point for computing partial functional derivatives *)


(* ::Text:: *)
(*Power counting parameter*)


$\[Epsilon]FD::usage="$\[Epsilon]FD is a power counting parameter for the expansion of gauge covariant derivatives.";


(* directly drop order $\[Epsilon]FD^3 terms *)
$\[Epsilon]FD /: Power[$\[Epsilon]FD,n_/;n>2] := 0


(* ::Text:: *)
(*Option test*)


Options[FD] = {Grassmann->True};


OptionTest[_, Grassmann] = BooleanQ;


OptionMessage[Grassmann, func:FD, val_] := Message[General::optexpectsval, Grassmann, func, val, "boolean (True or False)"];


FD[x_, y_, OptionsPattern[]] := Module[
	{
		result,
		arg = Expand[x]
	}
	,

	(* computation of the partial functional derivative and simplifications *)
	result = FuncDSimplify@FuncD[arg, y, Grassmann->OptionValue[Grassmann]];

	Return[result]
]


(* ::Text:: *)
(*Check arguments*)


FD::invalidargument =
"The second argument `1` of FuncD is invalid. It must be either of the following forms: Field[...], Bar[Field[...]].";


FD[_,y_,OptionsPattern[]] := Message[FD::invalidargument,y] /; !MatchQ[y, Field[___] | Bar@Field[___] | Transp@Field[___] | Transp@Bar@Field[___] | HoldComplete[___] | NonCommutativeMultiply[DiracProduct[GammaCC],_] | NonCommutativeMultiply[_,DiracProduct[GammaCC]]]


(* ::Subsection:: *)
(*Expansion of Field Strength tensors and Covariant Derivatives*)


(* ::Subsubsection::Closed:: *)
(*Expand fluctuations of vector fields*)


ExpandVectorFluctuations[lag_, f1_, f2__]:= ExpandVectorFluctuations[ExpandVectorFluctuations[lag, f2], f1]


(* This functions expands out the vector field fluctuations of the vector field y in Field Strengh tensors and Covariant Derivatives *)
ExpandVectorFluctuations[x_, y_]:=
ExpandVectorFluctuations[x, y]= Module[
	{
		gaugeGroup,
		yField = First@Cases[y,Field[___],All],
		yFieldLabel,
		arg = Expand[x]
	}
	,
	yFieldLabel = First[yField];

	(* Shift covariant derivatives and field-strength tensors in the case of gauge fields *)
	If[GaugeFieldQ@yField,

		(* shift explicit gauge fields *)
		arg = Expand[arg /. Field[yFieldLabel,type_Vector,rest___]:>(Field[BackgroundField[yFieldLabel],type,rest]+$\[Epsilon]FD*Field[yFieldLabel,type,rest])];

		(* remove powers for pattern matching *)
		arg = RemovePower[arg];

		(* determine gauge group *)
		gaugeGroup = First@First@FirstPosition[$GaugeGroups, KeyValuePattern[Field->First[yField]]];

		(* shift covariant derivatives iff field is charged under gauge group *)
		arg = arg/.Field[l_,t_,ind_List,deriv_List] :> ShiftCD[deriv,gaugeGroup,Field[l,t,ind,{}],$\[Epsilon]FD] /; (Length[deriv]>0 && FieldTransformsUnderGaugeGroupQ[Field[l,t,ind,{}], gaugeGroup]);
		arg = RemovePower@ReleaseHold[arg];

		(* shift covariant derivatives acting on field-strength tensors *)  (*crosscheck*)
		arg = arg/.FieldStrength[l_,lorentz_,ind_List,deriv_List] :> ShiftCD[deriv,gaugeGroup,FieldStrength[l,lorentz,ind,{}],$\[Epsilon]FD] /; (Length[deriv]>0 && FieldTransformsUnderGaugeGroupQ[FieldStrength[l,lorentz,ind,{}], gaugeGroup]);

		(* shift the field-strength tensors *)
		arg = arg /. FieldStrength[arguments___] :> ShiftFS[FieldStrength[arguments],$\[Epsilon]FD] /; MatchQ[First[{arguments}],$GaugeGroups[gaugeGroup][Field]];
		arg = RemovePower@ReleaseHold[arg];

	];

	(* Shift only the field-strength tensors in the case of vectors that are not gauge *)
	If[VectorFieldQ@yField && !GaugeFieldQ@yField,

		(* shift explicit vector fields *)
		arg = Expand[arg /. Field[yFieldLabel,type_Vector,rest___]:>(Field[BackgroundField[yFieldLabel],type,rest]+$\[Epsilon]FD*Field[yFieldLabel,type,rest])];

		(* remove powers for pattern matching *)
		arg = RemovePower[arg];

		(* shift the field-strength tensors *)
		arg = arg /. FieldStrength[First@yField,args___] :> ShiftFS[FieldStrength[First@yField,args],$\[Epsilon]FD];
	];

	(*Release BackgroundCD*)
	arg = arg/. BackgroundCD[ind__,exp_]:>CD[ind,exp] /.BackgroundFS[all___]:>FieldStrength[all];

	(* expand to quadratic order in the fluctuation *)
	arg = Normal@Series[Expand[ReleaseHold@arg],{$\[Epsilon]FD,0,2}];

	RelabelIndices[arg, Unique->True]
]


(* ::Subsubsection::Closed:: *)
(*Gauge fixing*)


(*ExpandVectorFluctuations[,_]:=0;*)
(*ExpandVectorFluctuations[GaugeFixing[exp_],_]:= exp*$\[Epsilon]FD^2
ExpandVectorFluctuations[Plus[a___,GaugeFixing[exp_],b___],f__]:= Plus[ExpandVectorFluctuations[If[a===Blank,0,Plus@a],f],exp*$\[Epsilon]FD^2,ExpandVectorFluctuations[If[b===Blank,0,Plus@b],f]]*)

(*Change CD to BackgroundCD and FS to BackgroundFS*)
GaugeFixing[exp_]:= exp/.Field[label_,type_,ind_, derivs:{\[Mu]_,\[Nu]___}]:>BackgroundCD[derivs, Field[label, type, ind,{}]]/.FieldStrength[all___]:>BackgroundFS[all];


(* ::Subsubsection::Closed:: *)
(*Check if a given field transforms under a given gauge group*)


(* returns true if the given gauge group is abelian *)
AbelianQ[gaugeGroup_] := GetGaugeGroups[gaugeGroup][Abelian];


(* checks whether the given field transforms under the given gauge group *)
FieldTransformsUnderGaugeGroupQ[Field[l_,t_,ind_List,deriv_List], gaugeGroup_] := If[
	AbelianQ[gaugeGroup],
	(* for abelian groups *)
	Return[
		MemberQ[GetFields[l][Charges][[;;,0]], gaugeGroup]
	]
	,
	(* for non-abelian groups *)
	Return[
		MemberQ[Map[GroupFromRep, GetFields[l][Indices]], gaugeGroup]
	]
]


(* checks whether the given field-strength tensor transforms under the given gauge group *)
FieldTransformsUnderGaugeGroupQ[FieldStrength[l_,lorentz_,ind_List,derive_List], gaugeGroup_] := Module[
	{groups}
	,
	If[AbelianQ[gaugeGroup],
		Return[False] (* FS can never be charged under abelian groups *)
	];

	groups= GetGaugeGroups[];
	If[MemberQ[groups, KeyValuePattern[Field->l]],
		If[MatchQ[gaugeGroup, First@ First@ First@ Position[groups, KeyValuePattern[Field->l]]],
			Return[True] (* If FS belongs go same group as gaugeGroup *)
		]
	];
	Return[False]
]


(* ::Subsubsection::Closed:: *)
(*Auxiliary function to check for gauge fields and heavy vectors*)


(* ::Text:: *)
(*gauge fields*)


GaugeFieldQ@ label_Symbol:= MemberQ[$GaugeGroups, KeyValuePattern[Field-> label] ];


GaugeFieldQ@ Field[label_, _Vector, __]:= GaugeFieldQ@ label;


GaugeFieldQ@ _Field= False;


(* ::Text:: *)
(*heavy vectors*)


VectorFieldQ@ Field[_, type_, __]:= MatchQ[type, Vector[_]]


(* ::Subsubsection::Closed:: *)
(*Shifting all covariant derivatives containing gauge field fluctuations*)


(* Expand the quantum in fluctuations of the gauge field associated to the given group in covariant derivatives *)
ShiftCD[indices_List, group_, field_, globalCounter_:1] := Module[{result,counter},

	(* if there is no CD just return the field*)
	If[Length[indices]==0, Return[field]];

	(* else: shift all the CDs *)
	If[AbelianQ[group],
		result = Expand@RecursiveExpandAbelianCD[indices, group, field, counter];
		,
		result = Expand@RecursiveExpandCD[indices, group, field, counter];
	];

	(* fluctuations should not be truncated to leading order when multiple derivatives are present *)
	(*
	(* truncate at leading order in the fluctuation *)
	result = Normal@Series[result,{counter,0,1}];
	*)

	(* remove counting parameter *)
	result = Expand[result/.counter->globalCounter];

	Return[result]
]


(* ::Text:: *)
(*Non-Abelian groups*)


RecursiveExpandCD[indices_List, group_, initField_, \[Lambda]_:1] := Module[
	{
		\[Mu]1,
		more\[Mu],
		indexAdj=Unique[],
		newRepIndex=Unique[],
		generatorIndices,
		field
	}
	,
	(* replace the group index in field *)
	{field,generatorIndices}=ChangeFieldIndex[initField, group, newRepIndex];

	(* split the Lorentz indices *)
	{\[Mu]1,more\[Mu]} = TakeDrop[indices,1];
	(* recursive formula for building up the fluctuation expansion in covariant derivatives *)
	If[Length[more\[Mu]]>=1,
		Return@Plus[
			Times[-1,
				\[Lambda],I,
				(*$GaugeGroups[group][Coupling][],*)
				$GaugeGroups[group][Field][\[Mu]1[[1,1]],indexAdj],
				FieldGenerators[GetFieldLabel[field], group, {indexAdj, First[generatorIndices], Last[generatorIndices]}],
				RecursiveExpandCD[more\[Mu],group,field,\[Lambda]]
			]
			,
			CD[\[Mu]1,RecursiveExpandCD[more\[Mu],group,initField,\[Lambda]]]
		]
		,
		If[Length[more\[Mu]]==0,
			Return@Plus[
				Times[-1,
					\[Lambda],I,
					(*$GaugeGroups[group][Coupling][],*)
					$GaugeGroups[group][Field][\[Mu]1[[1,1]],indexAdj],
					FieldGenerators[GetFieldLabel[field], group, {indexAdj, First[generatorIndices], Last[generatorIndices]}],
					field
				]
				,
				CD[\[Mu]1,initField]
			]
		]
	]
]


(* ::Text:: *)
(*Abelian groups*)


RecursiveExpandAbelianCD[indices_List, group_, initField_, \[Lambda]_:1] := Module[
	{
		\[Mu]1,
		more\[Mu],
		field= initField
	}
	,
	(* split the Lorentz indices *)
	{\[Mu]1,more\[Mu]} = TakeDrop[indices,1];
	(* recursive formula for building up the covariant derivative expansion *)
	If[Length[more\[Mu]]>=1,
		Return@Plus[
			Times[-1, (*crosschecked*)
				\[Lambda],I,
				(*$GaugeGroups[group][Coupling][],*)
				$GaugeGroups[group][Field][\[Mu]1[[1,1]]],
				FieldGenerators[GetFieldLabel[field], group],
				RecursiveExpandAbelianCD[more\[Mu],group,field,\[Lambda]]
			]
			,
			CD[\[Mu]1, RecursiveExpandAbelianCD[more\[Mu],group,initField,\[Lambda]]]
		]
		,
		If[Length[more\[Mu]]==0,
			Return@Plus[
				Times[-1, (*crosschecked*)
					\[Lambda],I,
					(*$GaugeGroups[group][Coupling][],*)
					$GaugeGroups[group][Field][\[Mu]1[[1,1]]],
					FieldGenerators[GetFieldLabel[field], group],
					field
				]
				,
				CD[\[Mu]1, initField]
			]
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*Auxiliary definitions*)


GetFieldLabel[field:(_Field | _FieldStrength)]:=Module[
	{
		label = First[field]
	}
	,
	(* treat Background vector fields *)
	If[Head[label]===BackgroundField,
		label = First[label]
	];
	(* treat fields expanded in the EFT power counting *)
	If[Head[label]===List,
		label = First[label]
	];
	Return[label]
]


FD::multiplegroupindices="The field `1` has several indices belonging to the group `2`.";


(* Function finding the index of a field that belongs to a specific group, returning that index and the field with the index replaced by a given label *)
ChangeFieldIndex[field_, group_, label_]:=Module[
	{
		(* find all indices of the field *)
		indices = Cases[field,_Index,All],
		groupIndex={}
	}
	,
	(* Find all indices belongin to the given group *)
	Do[
		If[MatchQ[GroupFromRep@Last[ind],group],
			AppendTo[groupIndex,ind]
		]
		,
		{ind,indices}
	];

	(* check that there is only one index *)
	If[Length[groupIndex]>1,
		Message[FD::multiplegroupindices,field,group];
		Abort[]
		,
		groupIndex = First[groupIndex]
	];

	(* return the field with replaced index and the original and new index labels *)
	Return[{
		field/.groupIndex->Index[label,Last[groupIndex]],
		{First[groupIndex],label}
	}]
]


(* ::Subsubsection::Closed:: *)
(*Shifting the field strength by their background configuration*)


(* ::Text:: *)
(*counter is a power-counting variable, default=1 -> no power counting*)


ShiftFS[FieldStrength[label_, {\[Mu]_Index,\[Nu]_Index}, {}, lorentzIndices_List], counter_:1] :=
	RelabelIndices[Plus[
		FieldStrength[label, {\[Mu],\[Nu]}, {}, lorentzIndices],
		counter*Field[label, Vector[\[Nu]], {}, Append[lorentzIndices,\[Mu]]],
		-counter*Field[label, Vector[\[Mu]], {}, Append[lorentzIndices,\[Nu]]]
	],Unique->True]


(* ::Text:: *)
(*Nb. this is not a conclusive way of telling if a field strength tensor is for a gauge field or not*)


ShiftFS[FieldStrength[label_, {\[Mu]_Index,\[Nu]_Index}, gaugeIndices:{Index[a_,group_[adj]]}, lorentzIndices_List], counter_:1] :=
	Module[{b,c},
	RelabelIndices[Plus[
		FieldStrength[label, {\[Mu],\[Nu]}, gaugeIndices, lorentzIndices],
		counter*Field[label, Vector[\[Nu]], gaugeIndices, Append[lorentzIndices,\[Mu]]],
		-counter*Field[label, Vector[\[Mu]], gaugeIndices, Append[lorentzIndices,\[Nu]]],
		counter^2(**$GaugeGroups[group, Coupling][]*)* CG[fStruct[group],{Index[a,group[adj]],Index[b,group[adj]],Index[c,group[adj]]}]*
		CD[lorentzIndices,Field[label, Vector[\[Mu]], {Index[b,group[adj]]},{}]*Field[label, Vector[\[Nu]], {Index[c,group[adj]]}, {}]]
	],Unique->True]
	]


ShiftFS[FieldStrength[label_, {\[Mu]_Index,\[Nu]_Index}, gaugeIndices:{Index[a_,group_[adj]],Index[b_,group_[adj]]}, lorentzIndices_List], counter_:1] :=
	Module[{c},
	 RelabelIndices[I CG[fStruct[group],{Index[a,group[adj]],Index[c,group[adj]],Index[b,group[adj]]}]*
	ShiftFS[FieldStrength[label, {\[Mu],\[Nu]}, {Index[c,group[adj]]}, lorentzIndices], counter],Unique->True]
	]


(* ::Section:: *)
(*FuncD \[LongDash] internal functionality*)


Options[FuncD] = {Grassmann->True};


(* ::Subsection:: *)
(*Differentiation rules*)


(* ::Subsubsection::Closed:: *)
(*Sum rule*)


FuncD[Plus[x1_,x2__],y_,OptionsPattern[]] :=
	Plus[
		FuncD[x1,y,Grassmann->OptionValue[Grassmann]],
		FuncD[Plus[x2],y,Grassmann->OptionValue[Grassmann]]
	]


(* ::Subsubsection::Closed:: *)
(*Product rules*)


(* ::Text:: *)
(*Times[]*)


FuncD[Times[x1_,x2__],y_,OptionsPattern[]] :=
	Plus[
		Times[FuncD[x1,y,Grassmann->OptionValue[Grassmann]],x2],
		Times[x1,FuncD[Times[x2],y,Grassmann->OptionValue[Grassmann]]]
	]


(* ::Text:: *)
(*NonCommutativeMultiply[]*)


FuncD[NonCommutativeMultiply[x1_,x2__],y_,OptionsPattern[]] := Module[
	{
		fieldType,
		sign = 1
	}
	,
	(* when Grassmann minus sign should be considered *)
	If[OptionValue[Grassmann],
		(*extract field type*)
		fieldType = Cases[y,Field[_,typ_,___]:>typ,All];
		If[MatchQ[fieldType,{Fermion|Ghost}],
			(*get the number of anti-commuting fields in x1*)
			sign = (-1)^(Plus@@Cases[x1,Field[_,Fermion|Ghost,___]->1,All]);
		];
	];

	(*apply the product rule*)
	Plus[
		NonCommutativeMultiply[FuncD[x1,y,Grassmann->OptionValue[Grassmann]],x2],
		sign * NonCommutativeMultiply[x1,FuncD[NonCommutativeMultiply[x2],y,Grassmann->OptionValue[Grassmann]]]
	]
] /; FreeQ[x1,Plus]



FuncD[NonCommutativeMultiply[x:(Field[___] | Bar@Field[___] | Transp@Field[___] | Transp@Bar@Field[___])], y_, OptionsPattern[]] := NonCommutativeMultiply@FuncD[x, y, Grassmann->OptionValue[Grassmann]]


(* ::Text:: *)
(*NCProduct[]*)


FuncD[NCProduct[x1_,x2__],y_,OptionsPattern[]] := Module[
	{
		fieldType,
		sign = 1
	}
	,
	(*when Grassmann minus sign should be considered*)
	If[OptionValue[Grassmann],
		(*extract field type*)
		fieldType = Cases[y,Field[_,typ_,___]:>typ,All];
		If[MatchQ[fieldType,{Fermion|Ghost}],
			(*Get the number of anti-commuting fields in x1*)
			sign = (-1)^(Plus@@Cases[x1,Field[_,Fermion|Ghost,___]->1,All]);
		];
	];

	(*apply the product rule*)
	Plus[
		NCProduct[FuncD[x1,y,Grassmann->OptionValue[Grassmann]],x2],
		sign * NCProduct[x1,FuncD[NCProduct[x2],y,Grassmann->OptionValue[Grassmann]]]
	]
] /; FreeQ[x1,Plus]


FuncD[NCProduct[x_], y_, OptionsPattern[]]:= NCProduct[FuncD[x, y, Grassmann->OptionValue[Grassmann]]]


(* ::Subsubsection::Closed:: *)
(*Power rule*)


FuncD[Power[x_,pow_Integer],y_,OptionsPattern[]]:=
	pow*Power[x,pow-1]*FuncD[x,y,Grassmann->OptionValue[Grassmann]] /; pow>=2


(* ::Subsection:: *)
(*Evaluation of derivatives*)


(* ::Subsubsection::Closed:: *)
(*Vanishing derivatives*)


FuncD[x_,_,OptionsPattern[]] := 0 /; FreeQ[x, Field] (* FieldStrength does not contain vector field fluctuations after the shift in FD *)


FuncD[Bar@Field[___], f:Field[___],OptionsPattern[]] := 0
FuncD[Field[___], Bar[f:Field[___]],OptionsPattern[]] := 0


FuncD[Field[_,_,_,\[Mu]_List], f:Field[_,_,_,\[Nu]_List], OptionsPattern[]] := 0 /; (Length[\[Mu]]!=Length[\[Nu]])
FuncD[Bar@Field[_,_,_,\[Mu]_List], Bar[f:Field[_,_,_,\[Nu]_List]], OptionsPattern[]] := 0 /; (Length[\[Mu]]!=Length[\[Nu]])


FuncD[Field[label1_,__], f:Field[label2_,__], OptionsPattern[]] := 0 /; (label1=!=label2)
FuncD[Bar@Field[label1_,__], Bar[f:Field[label2_,__]], OptionsPattern[]] := 0 /; (label1=!=label2)


(* ::Subsubsection::Closed:: *)
(*Non-vanishing derivatives*)


FuncD::repeatedindices = "The indices `1` in argument x and the indices `2` in argument y of FuncD[x,y] must not have identical labels.";


FuncD::invalidtype = "Invalid combination of field types found: `1` and `2`.";


FuncD::indexnumber = "Unmatched number of group indices found.";


FuncD::unmatchedindex = "Unmatched types of group indices found.";


FuncD[Bar@Field[label_,rest1__], Bar@Field[label_,rest2__], OptionsPattern[]] := FuncD[Field[label,rest1], Field[label,rest2], Grassmann->OptionValue[Grassmann]]


FuncD[Field[label_,type1_,index1_,derivs1_], Field[label_,type2_,index2_,derivs2_], OptionsPattern[]] :=
	Module[
		{ result = 1 }
		,
		(*check that the group indices match up*)
		If[Length[index1]!=Length[index2],
			Message[FuncD::indexnumber];
			Abort[]
			,
			If[(index1/.Index[_,rep_]:>rep)=!=(index2/.Index[_,rep_]:>rep),
				Message[FuncD::unmatchedindex];
				Abort[]
			]
		];

		(*check that there are no duplicate index labels *)
		If[Or@@(!FreeQ[index1,#]&/@index2),
			Message[FuncD::repeatedindices, index1, index2];
			Abort[]
		];
		If[Or@@(!FreeQ[derivs1,#]&/@derivs2),
			Message[FuncD::repeatedindices, derivs1, derivs2];
			Abort[]
		];

		(*functionality depending on field type*)
		Switch[{type1,type2},
			(*fermions*)
			{Fermion,Fermion},
				result *= 1, (*This used to be the identity in Dirac space*)
			(*vectors*)
			{Vector[Index[_,Lorentz]],Vector[Index[_,Lorentz]]},
				If[type1===type2,
					Message[FuncD::repeatedindices, First@type1, First@type2]; Abort[],
					result *= Metric[First@type1, First@type2]
				],
			(*scalars*)
			{Scalar,Scalar},
				result *= 1,
			(*ghosts*)   (*to be crosschecked*)
			{Ghost,Ghost},
				result *= 1,
			(*other*)
			{_,_},
				Message[FuncD::invalidtype, type1, type2];
				Abort[]
		];

		(*include Delta[_,_] for all group indices*)
		Do[
			result *= Delta@@pair
			,
			{pair,Transpose[{index1,index2}]}
		];

		(*include Metric[_,_] for all Lorentz indices of the covariant derivatives*)
		Do[
			result *= Metric@@pair
			,
			{pair,Transpose[{derivs1,derivs2}]}
		];

		Return[result]
	] /; (Length[derivs1]==Length[derivs2])


(* ::Subsection:: *)
(*FuncD with transposed and charge-conjugated fields *)


(* ::Subsubsection::Closed:: *)
(*FuncD w.r.t. Field[...]*)


FuncD[x_, Transp[y_], OptionsPattern[]] := (TransposeThisSpinChain ** FuncD[x, y, Grassmann->OptionValue[Grassmann]])/; FreeQ[x,Plus]


FuncD[Transp[x:(Field[___] | Bar@Field[___])], y:Except[_NonCommutativeMultiply] , OptionsPattern[]] := (TransposeThisSpinChain ** FuncD[x, y, Grassmann->OptionValue[Grassmann]])


(* ::Text:: *)
(*If y is a NonCommutativeMultiply[...] then y is the charge conjugate of some field for which the rules below should be used instead, since transposition generates minus signs due to the presence of GammaCC.*)


(* ::Subsubsection::Closed:: *)
(*FuncD  w.r.t. CC@Field[...]*)


FuncD[
	Bar[x_Field],
	NonCommutativeMultiply[DiracProduct[GammaCC], Transp@Bar[y_Field]],
	OptionsPattern[]
] := NonCommutativeMultiply[
	DiracProduct[GammaCC],
	TransposeThisSpinChain,
	FuncD[Bar[x], Bar[y], Grassmann->OptionValue[Grassmann]]
]


FuncD[
	Transp@Bar[x_Field],
	NonCommutativeMultiply[DiracProduct[GammaCC], Transp@Bar[y_Field]],
	OptionsPattern[]
] := (-1)*NonCommutativeMultiply[
	DiracProduct[GammaCC],
	FuncD[Bar[x], Bar[y], Grassmann->OptionValue[Grassmann]]
]


FuncD[
	x:(Field[___] | Transp@Field[___]),
	NonCommutativeMultiply[DiracProduct[GammaCC], Transp@Bar[y_Field]],
	OptionsPattern[]
] := 0


(* ::Subsubsection::Closed:: *)
(*FuncD  w.r.t. Bar@CC@Field[...]*)


FuncD[
	x_Field,
	NonCommutativeMultiply[Transp[y_Field], DiracProduct[GammaCC]],
	OptionsPattern[]
] := NonCommutativeMultiply[
	DiracProduct[GammaCC],
	TransposeThisSpinChain,
	FuncD[x, y, Grassmann->OptionValue[Grassmann]]
]


FuncD[
	Transp[x_Field],
	NonCommutativeMultiply[Transp[y_Field], DiracProduct[GammaCC]],
	OptionsPattern[]
] := (-1) * NonCommutativeMultiply[
	DiracProduct[GammaCC],
	FuncD[x, y, Grassmann->OptionValue[Grassmann]]
]


FuncD[
	x:(Bar@Field[___] | Transp@Bar@Field[___]),
	NonCommutativeMultiply[Transp[y_Field], DiracProduct[GammaCC]],
	OptionsPattern[]
] := 0


(* ::Section:: *)
(*FuncDSimplify*)


Transp[TransposeThisSpinChain] := TransposeThisSpinChain


(* ::Text:: *)
(*Simplifies the expression returned by FuncD by transposing spin-chains where necessary*)


FuncDSimplify[expr_] := Module[{solution=expr},
	(* avoid unnecessary transpositions *)
	solution=solution//.HoldPattern@NonCommutativeMultiply[A___,TransposeThisSpinChain,B___,TransposeThisSpinChain,C___]:>NonCommutativeMultiply[A,B,C]/;FreeQ[B,TransposeThisSpinChain];
	solution=solution//.HoldPattern@NonCommutativeMultiply[A___,TransposeThisSpinChain,TransposeThisSpinChain,C___]:>NonCommutativeMultiply[A,C];
	(* transpose all necessary spin-chains *)
	solution=solution/.HoldPattern@NonCommutativeMultiply[A___,TransposeThisSpinChain,B___]:>Transp@NonCommutativeMultiply[A,B];
	Contract[solution]
]


(* ::Section:: *)
(*Variational derivative \[LongDash] VarD*)


Options[VarDraw] = {Grassmann->True, EFTOrder -> 6};


(* ::Subsection:: *)
(*Check Arguments*)


OptionMessage[Grassmann, func:VarDraw, val_] := Message[General::optexpectsval, Grassmann, func, val, "boolean (True or False)"];


VarDraw::invalidargument =
"The second argument `1` of VarDraw is invalid. It must be either of the following forms:
Field[_,_,_,{}], Bar[Field[_,_,_,{}]], Transp[Field[_,_,_,{}]], Transp[Bar[Field[_,_,_,{}]]], CConj[Field[_,_,_,{}]], Bar[CConj[Field[_,_,_,{}]]].";


SecondArgVarD = Alternatives[
		Field[_,_,_,{}],
		Bar@Field[_,_,_,{}],
		Transp@Field[_,_,_,{}],
		Transp@Bar@Field[_,_,_,{}],
		NonCommutativeMultiply[DiracProduct[GammaCC], Transp@Bar@Field[_,_,_,{}]],
		NonCommutativeMultiply[Transp@Field[___], DiracProduct[GammaCC]]
	];


VarDraw[
	_,
	y : Except@SecondArgVarD,
	OptionsPattern[]
] := Message[VarDraw::invalidargument,y]


(* ::Subsection:: *)
(*The total functional derivative*)


(* ::Subsubsection::Closed:: *)
(*Single application of variational derivative*)


CovariantLoop[lagrangian_, fields_List, opts:OptionsPattern[]]? OptionsCheck:=
CovariantLoop[lagrangian, fields, opts]= Module[{},0];


VarDraw[L_, x : SecondArgVarD, opts:OptionsPattern[]]:=
VarDraw[L, x, opts]= Module[
  	{
   		n,
   		arg,
   		result,
   		indices = {}
   	},

  	If[Head@(Cases[x, Field[___], All][[1,2]])===Vector,arg=Expand[L/$\[Epsilon]FD];, arg=L;];

  	(*Add FD w.r.t. Field[...] to result*)
  	result = FD[arg, x, Grassmann -> OptionValue[Grassmann]];

  	(*Add FD w.r.t. CD[Field[...]] to result, included number of CD[] depends on mass dimension*)
  	Do[
   		(*create unique lorentz indices for covariant derivatives*)
   		indices = Table[Unique[],{i, n}];
   		result += ((-1)^n) * FuncDSimplify@CD[
       			Reverse@indices, (* indices need to be reversed due to IBP *)
       			FD[arg, CD[indices, x], Grassmann -> OptionValue[Grassmann]]
       		];
   		,
   		{n, OptionValue[EFTOrder] - 2}
   	];

  	result
  ]


(* ::Subsubsection::Closed:: *)
(*Variational derivatives*)


Options[VarD] = {EFTOrder -> 6};
Options[VarD1] = {EFTOrder -> 6};


(* ::Text:: *)
(*Master function removing all terms which trivially vanish under differentiation*)


VarD[lag_, fs__, opt:OptionsPattern[]]:= Module[{gr, terms, labels},
	labels= Cases[{fs}, Field[l_, __]-> l, All];
	terms= BetterExpand@ lag;
	terms= If[Head@ terms === Plus, List@@ terms, {terms}];

	(* Remove terms which trivially vanish under differentiation *)
	Do[
		terms= If[MemberQ[$GaugeGroups, KeyValuePattern[Field-> lab]],
			(* gauge field *)
			gr= First@ GetGaugeGroupByProperty[Field -> lab];
			Select[terms, NotTrivialWRTGaugeFieldQ[gr]]
		,
			(* matter field *)
			Cases[terms, a_/;! FreeQ[a, Field[lab, __]] ]
		];
	, {lab, labels}];

	terms= ExpandVectorFluctuations[Total@ terms, Sequence@@ DeleteDuplicatesBy[{fs},
		FirstCase[{#}, Field[l_, __] :> l, Nothing, All]& ] ];
	VarD1[terms, fs, opt]
];


(* ::Text:: *)
(*Checks if the term is trivially zero under differentiation w.r.t. the field gauge group gr*)


NotTrivialWRTGaugeFieldQ[gr_][term_] := !FreeQ[term, f : (Field[__, {__}] | _FieldStrength) /; FieldTransformsUnderGaugeGroupQ[f, gr]] ||
	!FreeQ[term, (Field|FieldStrength)[$GaugeGroups[gr, Field], __]];


(* ::Text:: *)
(*Variational derivative w.r.t. 1 field*)


VarD1[lag_, f1 : SecondArgVarD, opt:OptionsPattern[]] :=
	RelabelIndices@VarDraw[lag, f1,opt]/.$\[Epsilon]FD->0 /. BackgroundField[l_]:>l


(* ::Text:: *)
(*Variational derivative w.r.t. 2 fields*)


(* THIS VERSION DOES NOT WORK: the result must we wrapped inside a FuncNCM with the OpenCD last *)
(*VarD1[lag_, f1 : SecondArgVarD, f2 : SecondArgVarD, opt:OptionsPattern[]] :=
	RelabelIndices@VarDraw[Expand[VarDraw[lag, f2,opt] * FuncNCM[OpenCD[{}]]], f1, opt]/.
		OpenCD[{}]-> 1/. $\[Epsilon]FD-> 0/. BackgroundField[l_]:>l*)


VarD1[lag_, f1 : SecondArgVarD, f2 : SecondArgVarD, opt:OptionsPattern[]] := Module[{res},
	res = RelabelIndices@VarDraw[RelabelIndices[VarDraw[lag, f2,opt] * FuncNCM[OpenCD[{}]],Unique->True], f1, opt]/.
		OpenCD[{}]-> 1/. $\[Epsilon]FD-> 0/. BackgroundField[l_]:>l;
	res = res/.Times->FuncNCM;
	res = res/.{FuncNCM[a___,ocd_OpenCD,b___]:>FuncNCM[a,b,ocd]}
]


(* ::Section:: *)
(*Equations of Motion & fluctuation operator*)


(* ::Subsection:: *)
(*EOM*)


Options[DeriveEOM] = {EFTOrder->6};


(* ::Subsubsection::Closed:: *)
(*EoM for one field*)


DeriveEOM[lagrangian_, field:(Field[label_,___,{}] | Bar@Field[label_,___,{}]), OptionsPattern[]] := Module[
	{
		eom,
		deriveField = Bar[field],
		fieldProps = GetFieldsUpdated@ label
	}
	,

	(* check for massive chiral fermions *)
	If[(fieldProps[Type] === Fermion) && (fieldProps[Chiral] =!= False) && (fieldProps[Mass] =!= 0) && !fieldProps[SelfConjugate],
		deriveField = CConj@deriveField
	];

	(* derive EOM *)
	eom = VarD[lagrangian, deriveField, EFTOrder->OptionValue[EFTOrder]];

	Return[eom]
]


(* ::Subsubsection::Closed:: *)
(*EoM for list of fields*)


DeriveEOM[lagrangian_, fields_List, OptionsPattern[]] := (DeriveEOM[lagrangian,#,EFTOrder->OptionValue[EFTOrder]]&/@fields)


(* ::Subsubsection::Closed:: *)
(*Errors*)


DeriveEOM::invalarg = "The argument `1` is invalid; only Field[...] or Bar[Field[...]] are allowed.";


DeriveEOM[_,arg:Except[List[___]|Field[___]|Bar@Field[___]], OptionsPattern[]] := Message[DeriveEOM::invalarg,arg]


(* ::Subsection:: *)
(*Fluctuation operator*)


Options[FluctuationOperator] = {EFTOrder -> 6};


(* ::Text:: *)
(*Acceptable input for functional derivatives *)


fdPattern= Alternatives[
		Field[_,_,_,{}],
		Bar@Field[_,_,_,{}],
		Transp@Field[_,_,_,{}],
		Transp@Bar@Field[_,_,_,{}],
		NonCommutativeMultiply[DiracProduct[GammaCC], Transp@Bar@Field[_,_,_,{}]],
		NonCommutativeMultiply[Transp@Field[___], DiracProduct[GammaCC]]
	];


FluctuationOperator::invldField= "`1` is not a valid field variable.";


FluctuationOperator[lag_, field1_, field2_, OptionsPattern[]] := Module[
		{f1, f2, eom, xTerm, f, sign, lagrangian= NCProduct[lag]},
	(*Remove projection operators from functional derivative field (chiral fields always have explicit projectors in the Lagrangian).*)
	{f1, f2}= {field1, field2}/. {DiracProduct[a___, (Transp@ _Proj| _Proj)]:> DiracProduct@ a}/.
		NonCommutativeMultiply@ x_-> x;

	Do[
		If[! MatchQ[f, fdPattern],
			Message[FluctuationOperator::invldField, f];
			Abort[];
		];
		,
		{f, {f1, f2}}
	];

	(*
	if field2 is fermionic or ghost add a minus sign:
	since VarD[\[ScriptCapitalL],field2] has an odd number of fermions and the fluctuation of field2 has to be anticommuted through the entire expression.
	 *)
	If[FreeQ[field2, Field[_, Fermion | Ghost, ___], All],
		sign= 1,
		sign= -1
	];

	(* apply functional derivatives *)
	xTerm= sign * VarD[lagrangian,f1,f2, EFTOrder->OptionValue[EFTOrder]];

	xTerm= SortNCProduct[xTerm]/.NCProduct->NonCommutativeMultiply;

	Return[xTerm]
]


(* ::Section:: *)
(*Auxiliary NCProduct*)


(* ::Subsection:: *)
(*Properties of NCProduct*)


NCProduct[a___, NCProduct[b__], c___]:= NCProduct[a, b, c]
NCProduct[] = 1;
NCProduct@NCProduct[x___]:=NCProduct[x];


(* ::Subsubsection::Closed:: *)
(*General properties*)


(* ::Text:: *)
(*NCProduct contractions*)


(* Extracting commuting objects *)
NCProduct[a___, b_?CommutativeQ, c___]/; FreeQ[b,_NonCommutativeMultiply]:= b * NCProduct[a, c]
NCProduct[a___, b_?CommutativeQ * x_, c___]/; FreeQ[b,_NonCommutativeMultiply]:= b * NCProduct[a, x, c]
(* Distributivity *)
NCProduct[a___, b_Plus, c___]:= NCProduct[a, #, c] & /@ b


NCProduct[a___, Times[b_NonCommutativeMultiply, c___], d___]:= NCProduct[a,b,Times[c],d]


SortNCProduct[expr_]:= Module[
	{
		rules={
			NCProduct[a___,b_NonCommutativeMultiply,c___]/;ClosedSpinChainQ[b] :> NCProduct[a,c,b],
			NCProduct[a___,x_NonCommutativeMultiply,y_NonCommutativeMultiply,b___] /;(ROpenSpinChainQ[x]&&LOpenSpinChainQ[y]) :> (-1)^(GrassmannCount[x]*GrassmannCount[y]) * NCProduct[a,y,x,b]
		}
	}
	,
	expr/.rules
]


(* counts the grassmann parameters in arg *)
GrassmannCount[arg_]:= Count[arg, Fermion|Ghost, All]
