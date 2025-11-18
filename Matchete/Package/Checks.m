(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Checks`*)


(* ::Subtitle:: *)
(*Various modules for checking if inputs are correct.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsection::Closed:: *)
(*Exported*)


PackageExport["CheckLagrangian"]


PackageExport["Hermiticity"]
PackageExport["ContractedIndices"]
PackageExport["ClosedSpinChains"]
PackageExport["CanonicallyNormalized"]
PackageExport["HeavyMassBasis"]
PackageExport["FreeOfHeavyTadpoles"]
PackageExport["ChargeNeutral"]
PackageExport["FreeOfGaugeFields"]
PackageExport["UndefinedObject"]
PackageExport["GaugeAnomalies"]
PackageExport["DetailedOutput"]


(* ::Subsection::Closed:: *)
(*Internal*)


PackageScope["LagrangianLikeCheck"]


PackageScope["HermitianQ"]
PackageScope["UncontractedIndices"]
PackageScope["OpenSpinChains"]
PackageScope["ComplexSpinChains"]
PackageScope["HeavyMassBasisQ"]
PackageScope["HeavyTadpoleQ"]
PackageScope["HeavyTadpoles"]
PackageScope["ChargeNeutralQ"]
PackageScope["KineticCanonicalQ"]
PackageScope["GaugeVectorOccurenceQ"]
PackageScope["GaugeSingletQ"]
PackageScope["GaugeAnomaliesQ"]


PackageScope["GetCharge"]
PackageScope["GaugeAnomalyContribution"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection::Closed:: *)
(*Exported*)


CheckLagrangian::usage = "CheckLagrangian[Lagrangian, options] performs a series of checks on the input Lagrangian. 
The options {Hermiticity, ContractedIndices, ClosedSpinChains, CanonicallyNormalized, HeavyMassBasis, ChargeNeutral, FreeOfHeavyTadpoles, FreeOfGaugeFields, UndefinedObject, GaugeAnomalies} take values True/False and determine which of the checks are performed. 
If the option DetailedOutput->True, an association with the result of each test is returned.";


Hermiticity::usage    = ContractedIndices::usage = ClosedSpinChains::usage = CanonicallyNormalized::usage = HeavyMassBasis::usage = UndefinedObject::usage = ChargeNeutral::usage = FreeOfGaugeFields::usage = FreeOfHeavyTadpoles::usage = GaugeAnomalies::usage = "Option for the CheckLagrangian function. It takes the values True/False to determine whether the check is performed.";
DetailedOutput::usage = "Option for the CheckLagrangian function. It takes the values True/False. An association with the result of each test is returned when set to True.";


(* ::Subsection::Closed:: *)
(*Internal*)


HermitianQ::usage               = "HermitianQ[expr] returns True if the expression is hermitian, False otherwise.";
UncontractedIndices::usage      = "UncontractedIndices[expr] returns the list of open indices and mismatch between bared and not bared indices.";
OpenSpinChains::usage           = "OpenSpinChains[expr] returns the list of open spin chains in expr."
ComplexSpinChains::usage        = "ComplexSpinChains[expr] returns the list of spin chains in expr with more than 2 fermions."
HeavyMassBasisQ::usage          = "HeavyMassBasisQ[Lag] returns True if the Lagrangian is in the heavy fields mass basis, False otherwise.";
HeavyTadpoleQ::usage            = "HeavyTadpoleQ[operator] identifies terms with only one heavy field."
HeavyTadpoles::usage            = "HeavyTadpoles[Lag] returns the list of heavy tadpole terms with mass dimension less than 1. "
ChargeNeutralQ::usage           = "ChargeNeutralQ[term1 (+ term2 + ...)] returns True if the term(s) are neutral under all Abelian gauge groups defined, False otherwise.";				 
KineticCanonicalQ::usage        = "KineticCanonicalQ[Lag] returns True if the Lagrangian's kinetic terms are canonically normalized and diagonal, False otherwise.";
GaugeVectorOccurenceQ::usage    = "GaugeVectorOccurenceQ[expr] returns True if the expression contains gauge vectors outside field strength tensors or covariant derivatives, False otherwise";
GaugeSingletQ::usage            = "GaugeSingletQ[field] returns True if the field is a gauge singlet."
GaugeAnomaliesQ::usage          = "GaugeAnomaliesQ[Lag] returns True if the Lagrangian is gauge anomaly free, False otherwise.";


GetCharge::usage                = "GetCharge[term, Agroup, ChargeNeutral -> True/False] returns the total charge of the term under the abelian group Agroup. It issues a warning if not neutral when ChargeNeutral -> True.";
GaugeAnomalyContribution::usage = "GaugeAnomalyContribution[field/fieldList] returns the contributions to the gauge anomaly triangles from the field 'field' or the fields in the 'fieldList' list.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Modules*)


(* ::Subsection::Closed:: *)
(*Lagrange-like expression*)


(* ::Text:: *)
(*Throws an error if the overall form of an expression can not conceivably correspond to a Lagrangian or an operator.  *)


LagrangianLikeCheck::nofield= "The following term(s) in the expression does not contain fields or field strength tensors and the expression is not compatible with a Lagrangian or an operator: `1`";


LagrangianLikeCheck@ 0:= Null;
LagrangianLikeCheck@ expr_:= Module[{temp= TermsToList@ expr},
	If[(Or@@ FreeQ[Field|FieldStrength]/@ temp) && temp=!={0} ,
		Message[LagrangianLikeCheck::nofield, Format[Cases[temp,a_/;FreeQ[a,_Field|_FieldStrength,All],1],NiceForm]];
		Abort[]; 
	];
	Null
];


(* ::Subsection:: *)
(*Hermiticity*)


(* ::Text:: *)
(*Applies a FullSimplify to each term in a sum, while assuming that the arguments of all Log are real positive numbers, as long as the arguments do not contain explicit complex numbers.*)


SimplifyLogs[expr_]:=If[Head[expr]===Plus,
	FullSimplify[#,And@@DeleteDuplicates@Cases[#,Log[arg_/;FreeQ[arg,_Complex,All]]:>(arg>0),All]]&/@ CollectOperators[expr]
	,
	FullSimplify[expr,And@@DeleteDuplicates@Cases[expr,Log[arg_/;FreeQ[arg,_Complex,All]]:>(arg>0),All]]
]
SimplifyLogs[0]=0


(* ::Text:: *)
(*If "manifest" -> True return 2 booleans for the Lagrangian being Hermitian, and it being manifestly Hermitian.*)


HermitianQ[expr_,OptionsPattern[{"manifest"->False}]]:= If[OptionValue["manifest"],
	Module[{tmp=GreensSimplify[expr - Bar[expr]]},
		{SimplifyLogs[tmp/.lf_LF:>EvaluateLoopFunctions[lf]] === 0 , tmp === 0}
	]
	,
	SimplifyLogs[GreensSimplify[expr - Bar[expr]]/.lf_LF:>EvaluateLoopFunctions[lf]] === 0
];


(* ::Subsection::Closed:: *)
(*Contracted indices*)


UncontractedIndices@ 0= {};


UncontractedIndices[Lag_]:= Module[{indOpen=Flatten[FindOpenIndices/@ TermsToList@ Lag], indDummyBar},
	If[indOpen=!={}, Return[indOpen]];
	indDummyBar=Flatten[Complement[Bar/@FindDummyIndices[#],
			Flatten[(List@@#)//.{Power[x_,n_]/;Element[n,NonNegativeIntegers]:>ConstantArray[x,n],FieldStrength[_,ind_,indG_,indL_]:>Flatten@Join[ind,indG,indL], _?NumericQ->Nothing,_Proj->Nothing,NonCommutativeMultiply->List,DiracProduct->List,GammaM[ind__]:>List@ind,Field[_,Vector[indV_],ind_,indL_]:>Join[List@indV,ind,indL],Field[_,_,ind_,indL_]:>Join[ind,indL], Coupling[_,ind_,_]:>ind, CG[_,ind_]:>ind}]] 
		& /@TermsToList@ Lag];
	If[indDummyBar=!={}, Return[Bar/@indDummyBar]];
	{}
]


(* ::Subsection::Closed:: *)
(*Open spin chains*)


(* ::Text:: *)
(*Identifies and returns open and complex spin chains in an expression*)


OpenSpinChains[expr_]:= Select[Flatten@Join[Cases[{expr}, _NonCommutativeMultiply,Infinity],Cases[{{expr}/._NonCommutativeMultiply->0},Field[_,Fermion, ___] ,All]], OpenSpinChainQ@# &]


ComplexSpinChains[expr_]:= Select[Cases[{expr}, _NonCommutativeMultiply,Infinity], Count[#,Fermion, Infinity]>2 &]


(* ::Text:: *)
(*Identifies inconsistent spin chains *)


InconsistentSpinChains@ expr_:= Cases[{expr},
		NonCommutativeMultiply[Transp@ _Field, DiracProduct[Except[GammaCC], ___], ___]|
		NonCommutativeMultiply[___, DiracProduct[Except[GammaCC], ___], Transp@ Bar@ _Field]
	, Infinity]


(* ::Subsection::Closed:: *)
(*Canonically normalized kinetic term*)


KineticCanonicalQ[L_]:=Module[{fields, LShouldBe,LIs},
	fields=OccuringFields[L];
	LShouldBe=KineticTerms[Sum[FreeLag[f],{f, fields}]];
	LIs = KineticTerms[L];
	Return[GreensSimplify[LShouldBe -LIs] === 0]
]


(* ::Subsection:: *)
(*Heavy mass basis*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


HeavyMassBasisQ::nomasscoupling  = "There is no Coupling in the coefficient of the heavy mass terms `1`."
HeavyMassBasisQ::mixedfields     = "There should not be any mixed field mass terms with coupling with EFT order 0."


(* ::Subsubsection::Closed:: *)
(*HeavyMassBasisQ [NEW]*)


(*check if heavy mass is diagonal in the input*)
HeavyMassBasisQ[Lag_]:=Module[{terms,coupling,fields,indicestypes},
	terms = Collect[Contract@ Operator@ IsolateMassTerms[Lag, Heavy->True], _Operator] //TermsToList;
	Catch[
	If[terms=={0},Throw[True]];
	(*for each term...*)
	(fields = Cases[#,_Field,Infinity];
	(*...check if it is twice the same field*) 
	If[!MatchQ[fields, {Field[f_,__],Field[f_,__]}], Message[HeavyMassBasisQ::mixedfields]; Throw[False]];
	(*...check there is at least one Coupling[...] in the coefficient*) 
	If[FreeQ[#, _Coupling, All], Message[HeavyMassBasisQ::nomasscoupling,#]; Throw[False]]
	) & /@ terms;
	Throw[True]]
];


(* ::Subsubsection::Closed:: *)
(*HeavyMassTermQ [OLD]*)


(*PackageScope["HeavyMassTermQ"]*)


(*HeavyMassTermQ::usage           = "HeavyMassTermQ[operator] returns True if the operator contains two heavy fields and nothing else."*)


(*(*terms with exactly two fields, no derivatives and at least one heavy field count as heavy mass terms*)
HeavyMassTermQ[(c_:1)q_Operator]:=
		MatchQ[OperatorType[q,CountEoMDerivatives->True], {{Field[f1_,_,__],Field[f2_,_,__]},0}/; (GetFields[f1][Heavy]||GetFields[f2][Heavy])]*)


(* ::Subsubsection::Closed:: *)
(*HeavyMassBasisQ [OLD]*)


(*HeavyMassBasisQ::masscoupling    = "There should only be one coupling in the heavy mass terms."
HeavyMassBasisQ::masslabel       = "The coupling label `1` and the mass label of the heavy field `2` should be the same. To find the defined mass label, use GetFields[`2`]."*)


(*(*check if heavy mass is diagonal in the input*)
HeavyMassBasisQ[Lag_]:=Module[{terms,coupling,fields,indicestypes},
	terms = (Operator[Lag]//.o_Operator /;(!HeavyMassTermQ[o]):>0) //TermsToList;
	terms = Select[terms, OperatorDimension[#/._Operator:>1]<1 &];
	Catch[
	If[terms=={0},Throw[True]];
	(*for each term...*)
	(fields = Cases[#,_Field,Infinity];
	(*...check if it is twice the same field*) 
	If[!MatchQ[fields, {Field[f_,__],Field[f_,__]}],Message[HeavyMassBasisQ::mixedfields];Throw[False]];
	coupling=Abs[#/.{_Operator:>If[GetFields[fields[[1,1]]][SelfConjugate]===True || GetFields[fields[[1,1]]][Chiral]=!=False,2,1]}]//.{Abs[x_]:>x,Power[m_,_]:>m};
	(*...check if there is only one coupling*)
	If[Head[Echo@coupling]=!=Coupling, Message[HeavyMassBasisQ::masscoupling];Throw[False]];
	(*...check if the coupling has the same label as the field mass*)
	If[GetFields[First@First@fields][Mass]=!=First@coupling,Message[HeavyMassBasisQ::masslabel,First@coupling, First@First@fields];Throw[False]];
	) & /@ terms;
	Throw[True]]
];*)


(* ::Subsection::Closed:: *)
(*Heavy tadpoles*)


HeavyTadpoleQ[(c_:1)op_Operator, cond_]:=
		MatchQ[OperatorType[op,CountEoMDerivatives->True], {{Field[f_,_,__]},0}/; (cond[f])]


HeavyTadpoles[Lag_]:=Module[{tadterms, massterms=IsolateMassTerms[Lag,Heavy->True], heavyQ},
	heavyQ = !FreeQ[massterms,#,All]&; (*heavy masss mixing was checked before*)
	tadterms = (Operator[Lag]//.o_Operator /;(!HeavyTadpoleQ[o, heavyQ]):>0) //TermsToList;
	tadterms = Select[tadterms, OperatorDimension[#/._Operator:>1]<1 &];
	tadterms
];


(* ::Subsection::Closed:: *)
(*Charge Neutral*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


GetCharge::ChargeNeutral = "The term `1` is not invariant under `2` gauge group."


(* ::Subsubsection::Closed:: *)
(*GetCharge*)


GetCharge[field_Symbol,Agroup_]:=Cases[GetFields[field][Charges],_Agroup]/.{Agroup[x_]}:>x/.{}->0;
GetCharge[Field[f_,___],Agroup_]:=GetCharge[f,Agroup]
GetCharge[Bar@Field[f_,___],Agroup_]:=-GetCharge[f,Agroup]
GetCharge[list_List,Agroup_]:=GetCharge[#,Agroup]& /@list


GetCharge[term:Except[_Symbol|_Field|_Bar|_Power|_List],Agroup_,OptionsPattern[{ChargeNeutral->False}]]:=Module[{tot},
			tot=Total[GetCharge[#,Agroup]& /@ Cases[RemovePower@term,Bar@Field@___, Infinity]] //Expand;
			tot=1/2 Total[GetCharge[#,Agroup]& /@ Cases[RemovePower@term,_Field, Infinity]]+  tot //Expand;
			If[ tot =!= 0 && OptionValue@ChargeNeutral == True, Message[GetCharge::ChargeNeutral,Format[term,NiceForm], Agroup]];
			Return[tot]
			]


(* ::Subsubsection::Closed:: *)
(*ChargeNeutralQ*)


ChargeNeutralQ[Lag_]:= Total@(Abs/@Flatten@Outer[GetCharge, TermsToList@ Lag,Keys@Select[$GaugeGroups,#[Abelian] === True&],{ChargeNeutral->True}]) == 0


(* ::Subsection::Closed:: *)
(*Isolated gauge field*)


GaugeVectorOccurenceQ[expr_]:= !FreeQ[expr, Field[Alternatives@@(Query[Transpose]@GetGaugeGroups[])[Field],___], Infinity];


(* ::Subsection:: *)
(*Gauge anomalies*)


(* ::Subsubsection::Closed:: *)
(*GaugeSingletQ*)


GaugeSingletQ::NotAField = "'`1`' is not a previously defined field.";


GaugeSingletQ[field_]:=Module[{},
(* Check that the argument is a previously defined field *)
	If[!MemberQ[Keys@$FieldAssociation,field],
		Message[GaugeSingletQ::NotAField,field];
		Abort[]
	];

	Join[$FieldAssociation[field][Charges],Intersection[Keys@$GaugeGroups,GroupFromRep/@$FieldAssociation[field][Indices]]]==={}
]


(* ::Subsubsection::Closed:: *)
(*Function that compute gauge anomalies*)


GaugeAnomalyContribution[fieldList_List]:=Merge[GaugeAnomalyContribution[#]&/@fieldList,Total];


GaugeAnomalyContribution[field_Symbol]:=Module[
	{
	i,j,k,
	A=Global`d$$1,B=Global`d$$2,C=Global`d$$3,
	ind=$FieldAssociation[field][Indices],
	gaugeInd=Select[$FieldAssociation[field][Indices],MemberQ[Keys@$GaugeGroups,GroupFromRep@#]&],
	charges=$FieldAssociation[field][Charges],
	chiralFactor=If[$FieldAssociation[field][Chiral]===LeftHanded, 1, -1],
	nonAbelianAnomalies,
	nonAbelianAnomalyFactors,
	mixedAnomalies,
	mixedAnomalyFactors,
	AbelianAnomalies,
	AbelianAnomalyFactors,
	gravityAnomalies,
	gravityAnomalyFactors
	},
	
	If[$FieldAssociation[field][Chiral]===False||GaugeSingletQ[field]===True,Return[<||>]];

	(* Contruct all possible non-Abelian triangles and compute the corresponding anomalies, given by Tr[{T[rep],T[rep]}T[rep]]*)
	nonAbelianAnomalies=Table[{i,i,i},{i,gaugeInd}];
	nonAbelianAnomalyFactors=Head/@#->If[Length[ind]>1,Times@@DimRep/@DeleteCases[ind,#[[1]]],1]chiralFactor ContractCGs[(CG[gen[#[[1]]],{A,i,j}]CG[gen[#[[3]]],{B,j,k}]+CG[gen[#[[2]]],{B,i,j}]CG[gen[#[[1]]],{A,j,k}])CG[gen[#[[3]]],{C,k,i}]]&/@nonAbelianAnomalies;

	(* Contruct all possible mixed non-Abelian/Abelian triangles and compute the corresponding anomalies, given by Tr[T[rep] T[rep]] U1[Charge] *)
	mixedAnomalies=Flatten[Table[{i,i,j},{i,gaugeInd},{j,charges}],1];
	mixedAnomalyFactors=Head/@#->If[Length[ind]>1,Times@@DimRep/@DeleteCases[ind,#[[1]]],1]chiralFactor ContractCGs[CG[gen[#[[1]]],{A,i,j}]CG[gen[#[[2]]],{B,j,i}]] #[[3]][[1]]&/@mixedAnomalies;

	(* Contruct all possible Abelian triangles and compute the corresponding anomalies, given by U1[Charge]^3 *)
	(*AbelianAnomalies=If[Length[charges]===1,
		Table[{i,i,i},{i,charges}],
		Flatten[DeleteCases[#/.{},_?NumericQ]]&/@(Expand[(Plus@@charges)^3]/.{Times->List,Plus->List,Power->Apply[List]@*ConstantArray})
	];*)
	AbelianAnomalies=Switch[Length[charges],
		0,
			{}, (* Exception for purely non-Abelian theories *)
		1,
			Table[{i,i,i},{i,charges}],
		_,
			Flatten[DeleteCases[#/.{},_?NumericQ]]&/@(Expand[(Plus@@charges)^3]/.{Times->List,Plus->List,Power->Apply[List]@*ConstantArray})
	];
	AbelianAnomalyFactors=Head/@#->If[Length[ind]>0,Times@@DimRep/@ind,1]chiralFactor #[[1]][[1]] #[[2]][[1]]#[[3]][[1]]&/@AbelianAnomalies;

	(* Contruct all possible mixed Gravity/Abelian triangles and compute the corresponding anomalies, given by U1[Charge] *)
	gravityAnomalies=Table[{"Gravity"[],"Gravity"[],i},{i,charges}];
	gravityAnomalyFactors=Head/@#->If[Length[ind]>0,Times@@DimRep/@ind,1]chiralFactor #[[3]][[1]]&/@gravityAnomalies;

	Association@Join [nonAbelianAnomalyFactors,mixedAnomalyFactors,AbelianAnomalyFactors,gravityAnomalyFactors]
]


(* ::Subsubsection::Closed:: *)
(*Function that tests if the Lagrangian has anomalies*)


GaugeAnomaliesQ[lagrangian_]:=Total@Values@GaugeAnomalyContribution@OccuringFields@lagrangian===0;


(* ::Section:: *)
(*Check Lagrangian function*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


CheckLagrangian::Hermiticity          = "The input Lagrangian is not hermitian. The non-Hermitian part \[ScriptCapitalL]-\[ScriptCapitalL]\[ConjugateTranspose] is `1`";
CheckLagrangian::ContractedIndices    = "The input Lagrangian contains open indices or a mismatch of bared and not bared indices.";
CheckLagrangian::ClosedSpinChains     = "The input Lagrangian contains open, complex, or inconsistent (using transposed spinors with no corresponding C-matrix) spin chains."
CheckLagrangian::CanonicallyNormalized = "The input Lagrangian kinetic part is not canonically normalized."
CheckLagrangian::HeavyMassBasis       = "The input Lagrangian is not in the mass basis for heavy fields."
CheckLagrangian::FreeOfHeavyTadpoles  = "The input Lagrangian contains heavy tadpoles that should be removed."
CheckLagrangian::ChargeNeutral        = "The input Lagrangian is not fully neutral."
CheckLagrangian::GaugeVectorOccurence = "The Lagrangian cannot contain gauge vector boson outside of covariant derivatives or field strength tensors."
CheckLagrangian::UndefinedObject      = "The object `1` is undefined."
CheckLagrangian::GaugeAnomalies       = "The input Lagrangian is not gauge anomaly free."


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[CheckLagrangian]={
		CanonicallyNormalized-> True,
		ChargeNeutral-> True,
		ContractedIndices-> True, 
		ClosedSpinChains-> True,
		DetailedOutput-> False,
		FreeOfGaugeFields-> True,
		GaugeAnomalies-> True,
		HeavyMassBasis-> True,
		FreeOfHeavyTadpoles-> True,
		Hermiticity-> True, 
		UndefinedObject-> True
	};


(* ::Subsubsection::Closed:: *)
(*CheckLagrangian*)


(* ::Text:: *)
(*Function for checking if the Lagrangian is valid for use in the matching routines *)


CheckLagrangian[Lagrangian_,opt:OptionsPattern[]]? OptionsCheck:=CheckLagrangian[Lagrangian,opt]=
	Module[{Lag=RelabelIndices@HcExpand@Lagrangian,DetOutput=<||>,OSpinChains,HeavTadpoles,UncIndices,ExtraHeads, mHermiticity=False,mContractedIndices=False,mClosedSpinChains=False, mCanonicallyNormalized=False,mMassBasis=False,mHeavyTadpoles=False,mChargeNeutral=False,mFreeOfGaugeFields=False,mUndefinedObject=False,mGaugeAnomalies=False},
	(*Checks that Lag is a series of terms with fields or FS tensors*)
	
	LagrangianLikeCheck@ Lag;
	
	(*check if all spin chains are closed and contain only 2 fermions*)
	If[OptionValue@ ClosedSpinChains,
		OSpinChains=Join[OpenSpinChains[Lag],ComplexSpinChains[Lag], InconsistentSpinChains[Lag]];
		If[ (mClosedSpinChains = OSpinChains =!= {}) , 
			Message[CheckLagrangian::ClosedSpinChains];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"Open/Complex/InconsistentSpinChains"->(Format[#,NiceForm]&/@OSpinChains)]]
		];
	
	(*check if L is hermitian *)
	(*HermitianQ cannot evaluate improperly contracted spinchains*)	
	If[OptionValue@Hermiticity && !mClosedSpinChains,
		If[(mHermiticity = !HermitianQ[Lag]), 
			Message[CheckLagrangian::Hermiticity,Format[(Lag - Bar[Lag] //GreensSimplify),NiceForm]]; 
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"Hermiticity"->!mHermiticity]]
		];
		
	(*check if all indices are contracted *)	
	If[OptionValue@ ContractedIndices,
		UncIndices = UncontractedIndices[Lag];
		If[ (mContractedIndices = UncIndices =!= {}), 
			Message[CheckLagrangian::ContractedIndices,UncIndices]; 
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"UncontractedIndices"->UncIndices]]
		];
		
	(*check if kinetaic part is canonically normalized*) 
	If[OptionValue@ CanonicallyNormalized,
		If[ (mCanonicallyNormalized = !KineticCanonicalQ[Lag]) , 
			Message[CheckLagrangian::CanonicallyNormalized]; 
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"CanonicallyNormalized"->!mCanonicallyNormalized]]
		];
		
	(*check if the heavy mass is diagonal*) 
	If[OptionValue@ HeavyMassBasis,
		If[ (mMassBasis= !HeavyMassBasisQ[Lag]), 
			Message[CheckLagrangian::HeavyMassBasis]; 
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"HeavyMassBasis"->!mMassBasis]]
		];
	
	If[OptionValue@ FreeOfHeavyTadpoles,
		HeavTadpoles=HeavyTadpoles[Lag];
		If[ (mHeavyTadpoles= HeavTadpoles =!= {}), 
			Message[CheckLagrangian::FreeOfHeavyTadpoles]; 
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"HeavyTadpoles"->(Format[#,NiceForm]&/@HeavTadpoles)]]
		];
		
	(*check if L is neutral under all charges*) 
	If[OptionValue@ ChargeNeutral,
		If[ (mChargeNeutral = !ChargeNeutralQ[Lag]) , 
			Message[CheckLagrangian::ChargeNeutral]; 
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"ChargeNeutral"->!mChargeNeutral]]
		];
		
	(*check if gauge fields occur outside of field strength*)
	If[OptionValue@ FreeOfGaugeFields,
		If[ (mFreeOfGaugeFields = GaugeVectorOccurenceQ[Lag]) , 
			Message[CheckLagrangian::GaugeVectorOccurence]; 
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"FreeOfGaugeFields"->!mFreeOfGaugeFields]]
		];
		
	(*check if all objects are defined*)	
	If[OptionValue@ UndefinedObject,
		ExtraHeads= DeleteDuplicates@ DeleteCases[Flatten[{Lag}//.
				{Times->List,
				Plus->List, 
				NonCommutativeMultiply->List, 
				Power[a_,___]:> a, Bar[a_]:> a , Transp[a_]:> a, Log[a_]:>a, 
				hbar->1, \[Mu]bar2->1, ev-> 1, \[Epsilon]->1, \[Pi]->1, \[ScriptD]->1}], 
			_Field | _CG | _Coupling | _?NumberQ | _DiracProduct  | _FieldStrength | _LCTensor |  _Delta];
		If[(mUndefinedObject)= (ExtraHeads=!={}),
			Message[CheckLagrangian::UndefinedObject,#]& /@ ExtraHeads;
		];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"AllObjectsDefined"->!mUndefinedObject]]
	];
	
    (*check if the Lagrangian is gauge anomaly free *)	
    If[OptionValue@ GaugeAnomalies,
		If[(mGaugeAnomalies = !GaugeAnomaliesQ@Lag), 
			Message[CheckLagrangian::GaugeAnomalies]; 
		];
		If[OptionValue@DetailedOutput, 
			AppendTo[DetOutput,"GaugeAnomalies"-> Format[Normal@GaugeAnomalyContribution@OccuringFields@Lag, NiceForm]]
		];
	];
	
	(* output the test results combined, we can output more details if wanted *)
	If[OptionValue@DetailedOutput,
		Return[DetOutput],
		Return[Nor[mHermiticity,mContractedIndices,mClosedSpinChains,mCanonicallyNormalized,mMassBasis,mChargeNeutral,mHeavyTadpoles,mFreeOfGaugeFields,mUndefinedObject,mGaugeAnomalies]]
	]
]
