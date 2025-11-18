(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`FieldRedef`*)


(* ::Subtitle:: *)
(*Paclet for performing field redefinition and obtaining on-shell EFT basis *)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["EOMSimplify"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["EoMSimplificationStep"]
PackageScope["ReduceField"]
PackageScope["FieldsToShift"]
PackageScope["CoefficientOperator"]
PackageScope["KineticTerms"]
PackageScope["KineticTermQ"]
PackageScope["OccuringFields"]
PackageScope["WriteEFTLambda"]
PackageScope["ShiftOrder"]
PackageScope["CDExp"]
PackageScope["OperatorType"]
PackageScope["GetMinOrder"]
PackageScope["GetMaxOrder"]
PackageScope["$EOMSimplifyEchoIntermediateLagrangians"]


PackageScope["RenormalizeMatterFields"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


EOMSimplify::usage =                    "EOMSimplify[L] takes a Lagrangian and removes redundant operators by the means of field redefinitions. The option 'EFTOrder' allows to define a maximal order in EFT counting. The standard value (Automatic) determines it from the input Lagrangian. The option 'DummyCoefficients' internally replaces matching coefficients with symbols, potentially improving performance when simplifying very large Lagrangians (standard is False).";
PrintEffectiveCouplings::usage =        "PrintEffectiveCouplings[L] prints a list of effective couplings defined by EOMSimplify appearing in the input Lagrangian L, along with their definitions. With the option Rules\[Rule] True, the function will return the replacement list for these couplings, instead.";
ReplaceEffectiveCouplings::usage =      "ReplaceEffectiveCouplings[L] resubstitutes effective couplings defined by EOMSimplify or by IntroduceEffective couplings.";
Rules::usage =                          "Rules is an option to specify whether output should be a list of replacement rules instead of printed.";
IntroduceEffectiveMasses::usage =       "IntroduceEffectiveMasses[L] replaces heavy mass terms with newly defined, generic mass couplings to compress the structure of the Lagrangian. Use ReplaceEffectiveCouplings to resubstitute them or PrintEffectiveCouplings to show a list of the couplings and their definitions.";


(* ::Chapter:: *)
(*Private:*)


(* ::Subsection:: *)
(*Operator classification*)


(* ::Subsubsection::Closed:: *)
(*Operator type*)


(* ::Text:: *)
(*Classifies an operator*)


Options[OperatorType]={CountEoMDerivatives -> False};

OperatorType[(c_:1)o_Operator, OptionsPattern[]]/;FreeQ[c, Operator]:=Module[
	{to=List@@ o,labels={}, derivs=0, tf,cEoMDer}
,
	cEoMDer=OptionValue[CountEoMDerivatives];

	If[cEoMDer,
		to=to//.{
			EoM[Field[l_,Scalar,i_,_]]:>  (derivs+=2;labels=Join[labels, {Field[l,Scalar,{},{}]}];tf[l,Scalar,i]),
			EoM[Field[l_,Fermion,i_,_]]:>  (derivs+=1;labels=Join[labels, {Field[l,Fermion,{},{}]}];tf[l,Fermion,i]),
			EoM[Bar@Field[l_,Fermion,i_,_]]:>  (derivs+=1;labels=Join[labels, {Field[l,Fermion,{},{}]}];Bar@tf[l,Fermion,i]),
			EoM[Transp@Field[l_,Fermion,i_,_]]:>  (derivs+=1;labels=Join[labels, {Field[l,Fermion,{},{}]}];Bar@tf[l,Fermion,i]),
			EoM[Field[l_,Vector,i_,_]]:>  (derivs+=1;labels=Join[labels, {Field[l,Vector,{},{}]}];tf[l,Vector,i])
			}
	];

	to=to//.{
		Field[l_, t_, i_, \[Mu]_/;Length[\[Mu]]>0]:> (derivs+=Length[\[Mu]];Field[l,t,i,{}]),
		FieldStrength[l_, mn_, i_, \[Mu]_/;Length[\[Mu]]>0]:> (derivs+=Length[\[Mu]];FieldStrength[l,mn,i,{}]),
		Field[l_,t_,i_,{}]:> (labels=Join[labels, {Field[l,t,{},{}]}];tf[l,t,i]),
		FieldStrength[l_,mn_,i_,{}]:> (labels=Join[labels, {FieldStrength[l,mn,{},{}]}];tf[l,mn,i])
	};
	{Sort[labels],derivs}
]


(* ::Subsubsection::Closed:: *)
(*Find kinetic terms*)


KineticTermQ[L_Plus]:=False
KineticTermQ[(c_:1)q_Operator]:=Module[{ot1,ot2},
	ot1=OperatorType[q,CountEoMDerivatives->True];
	ot2=OperatorType[q,CountEoMDerivatives->False];
	MatchQ[ot1, {{Field[_,Scalar,__],Field[_,Scalar,__]},2}]||
	(MatchQ[ot1, {{Field[_,Fermion,__],Field[_,Fermion,__]},1}]&&MatchQ[ot2, {{Field[_,Fermion,__],Field[_,Fermion,__]},0}])||MatchQ[ot1,{{FieldStrength[__],FieldStrength[__]},0}]
]


KineticTerms[L_]:=Operator[L]//.o_Operator /;(!KineticTermQ[o]):>0


OccuringFields[L_]:=Module[{f={},rField,rFS},
	L/. {Field[k_, A__]:>(AppendTo[f, k];rField[k,A]),FieldStrength[k_, A__]:>(AppendTo[f, k];rFS[k,A])};
	DeleteDuplicates[f]
]


(* ::Section:: *)
(*Field redefinitions*)


(* ::Subsubsection::Closed:: *)
(*Intermediate functions for power counting*)


(* make the EFT counting explicit so we can expand with the traditional series function *)
WriteEFTLambda[a_+b_, \[Lambda]_]:=WriteEFTLambda[a,\[Lambda]]+WriteEFTLambda[b,\[Lambda]]
WriteEFTLambda[a_*b_, \[Lambda]_]:=WriteEFTLambda[a,\[Lambda]]*WriteEFTLambda[b,\[Lambda]]
WriteEFTLambda[a_^b_, \[Lambda]_]:=WriteEFTLambda[a,\[Lambda]]^b
WriteEFTLambda[x_/;!(Head[x]===Plus ||Head[x]===Times), \[Lambda]_]:= x \[Lambda]^OperatorDimension[x]


(* splits Lagrangian by powercounting, returning a list of two Lagrangians {L0, L1}, where L0 contains all terms with of order splitOrder or below and L1 the rest *)
SplitLagrangianByPower[L_, splitOrder_]:=Module[{below=0,above=0,LExp},
	Do[
		If[OperatorDimension[term]<=splitOrder,
			below+=term,
			above+=term],
		{term,List@@(BetterExpand[L]+Nothing)}];
	Return[{below,above}]
]


(* further separates the Lagrangian into terms that will be unaffected by field redefinitions of f *)
SeparateByFieldContent[{L0_,L1_}, f_Symbol]:=Module[{l0,ol0=0,ol1=L1},
	l0 = List@@(L0+Nothing);
	(*l0 = TermsToList[L0];*)
	If[Length[l0]===0,Return[{L0,L1}]];
	Do[
		If[FreeQ[l0,f],
			ol1+=term,
			ol0+=term
		]
	,{term,l0}];

	{ol0,ol1}
]


DropFreeQ[L_,f_Symbol]:=Total@DeleteCases[List@@(Expand@L+Nothing), x_/;FreeQ[x,f]]
DropFreeQ[L_,f_List]:=Total@DeleteCases[List@@(Expand@L+Nothing), x_/;FreeQ[x,Alternatives@@f]]


GetMaxOrder[L_]:=Max[OperatorDimension/@(List@@(Expand@(L+Nothing)))]
GetMinOrder[L_]:=Min[OperatorDimension/@(List@@(Expand@(L+Nothing)))]


(* ::Subsection:: *)
(*Aux modules needed for the vector fields*)


(* ::Subsubsection::Closed:: *)
(*Shorthand for collecting by operators (without simplification)*)


CBOp[ex_]:=Collect[ex,_Operator]


(* ::Subsubsection::Closed:: *)
(*Expanding out the covariant derivatives*)


CDGaugeTerm[mu_,Field[l_,t_,i_,d_],gr_]:=Module[{a,i2,in,j,gf},
	(* get the index of the field we're contracting *)
	i2=First[Cases[i/.Index[index_,repr_]:> in[index, GroupFromRep@repr], in[_,gr]]/.in[f_,gr]:> f];
	-I ($GaugeGroups[gr][Coupling][])RelabelIndices[$GaugeGroups[gr][Field][mu,a]FieldGenerators[l,gr,{a,i2,j}](Field[l,t,i,d]/.i2 -> j),Unique->True]
]

CDGaugeTerm[mu_,FieldStrength[l_,t_,i_,d_],gr_]:=Module[{a,i2,in,j,gf},
	(* get the index of the field we're contracting *)
	i2=First[Cases[i/.Index[index_,repr_]:> in[index, GroupFromRep@repr], in[_,gr]]/.in[f_,gr]:> f];
	-I ($GaugeGroups[gr][Coupling][])RelabelIndices[$GaugeGroups[gr][Field][mu,a]FieldGenerators[l,gr,{a,i2,j}](FieldStrength[l,t,i,d]/.i2 -> j),Unique->True]
]

CDGaugeTermsAbel[mu_,FieldStrength[l_,t_,i_,d_]]:=Module[{grps,qs,As,gs},
	(* all the groups a field is charged under *)
	grps=Head/@($FieldAssociation[l][Charges]);
	(* all the charges corresponding to those *)
	qs=#[[1]]&/@($FieldAssociation[l][Charges]);
	(* the gauge fields and couplings *)
	As=$GaugeGroups[#][Field][mu]&/@grps;
	gs=$GaugeGroups[#][Coupling][]&/@grps;
	-I Total[qs As gs]FieldStrength[l,t,i,d]
]

CDGaugeTermsAbel[mu_,Field[l_,t_,i_,d_]]:=Module[{grps,qs,As,gs},
	(* all the groups a field is charged under *)
	grps=Head/@($FieldAssociation[l][Charges]);
	(* all the charges corresponding to those *)
	qs=#[[1]]&/@($FieldAssociation[l][Charges]);
	(* the gauge fields and couplings *)
	As=$GaugeGroups[#][Field][mu]&/@grps;
	gs=$GaugeGroups[#][Coupling][]&/@grps;
	-I Total[qs As gs]Field[l,t,i,d]
]


CDExp[ind_List/;Length[ind]>1,expr_]:=Module[{},
	CDExp[ind[[;;-2]],CDExp[{Last@ind},expr]]
]

CDExp[{Index[f_,Lorentz]}, x_]:=CDExp[{f},x]

CDExp[{mu_},a_+b_]:=CDExp[{mu},a]+CDExp[{mu},b]

CDExp[{mu_Symbol},a:_Times|_NonCommutativeMultiply]:=Module[{t},
	(* the argument is a composite object - first let the covariant derivative act *)
	t= CD[mu,a];

	(* now we replace the derivative wrt mu with CDExp again *)
	t/.{Field[l_,t_,i_, {Index[mu,Lorentz],v___}]:> CDExp[{mu},Field[l,t,i,{v}]]}]

CDExp[{mu_Symbol},Bar[x_]]:=Bar@CDExp[{mu},x]

CDExp[{mu_Symbol},Transp[x_]]:=Transp@CDExp[{mu},x]


CDExp[{mu_Symbol},Field[l_,t_,ind_,d_]]:=Module[{groups},
	groups=DeleteCases[DeleteDuplicates[GroupFromRep/@(ind//.Index[_,g_]:> g)],None];
	CD[mu,Field[l,t,ind,d]]+Total[CDGaugeTerm[mu,Field[l,t,ind,d],#]&/@groups]+CDGaugeTermsAbel[mu,Field[l,t,ind,d]]
]

CDExp[{mu_Symbol},FieldStrength[l_,li_,ind_,d_]]:=Module[{groups},
	groups=DeleteCases[DeleteDuplicates[GroupFromRep/@(ind//.Index[_,g_]:> g)],None];
	CD[mu,FieldStrength[l,li,ind,d]]+Total[CDGaugeTerm[mu,FieldStrength[l,li,ind,d],#]&/@groups]+CDGaugeTermsAbel[mu,FieldStrength[l,li,ind,d]]
]



(* ::Subsubsection::Closed:: *)
(*VectorFieldReplacement*)


(* out of a pair of V[\[Mu]] and \[Delta]V[\[Mu]], this function creates the appropriate replacement rules for fields and field-strengths *)
VectorFieldReplacement::indexError="The expression '`1`' does not have the same open indices as '`2`'.";
VectorFieldReplacement[{ls:Field[_,_Vector,__],rs_}]:=Module[{oinds,pinds,inds,lsp,rsp,rslist,rslistJ,rslistV,outV,outFS,nu,pnu,mu,pmu,jmu,jnu,vectorlabel,vectormu,vectorindices,pvectorindices,drvs},
(*If[!(Sort[FindOpenIndices[ls]]===Sort[FindOpenIndices[rs]]),Message[FieldReplacement::indexError, ls,rs];Abort[]];*)

(* set up patterns and their non-pattern counterparts for the replacement rules *)
oinds=FindOpenIndices[ls];
inds=Symbol["ix"<>ToString[#]]&/@Range@Length@oinds;
pinds=Pattern[#, _]&/@inds;

rsp=Inactive[CD][drvs,OperatorToNormalForm[rs,Unique->True]/.Thread[oinds->inds]];
lsp=(ls/.Field[l_,t_,i_,_List]:>Field[l,t,i,Pattern[drvs,_]])/.Thread[oinds -> pinds];

(* the replacement rule for the vector field itself *)
With[{trhs = rsp},
	outV=lsp :> RelabelIndices[trhs,Unique->True]/.Inactive[CD]->CD
];

(* two indices for the FS object *)
mu = Unique["a"];
pmu=Pattern[Evaluate@mu,_];
nu = Unique["b"];
pnu=Pattern[Evaluate@nu,_];

(* label of the field we are replacing *)
vectorlabel = (ls/.Field[a_,__]:> a);
vectormu=(ls/.Field[_,Vector@Index[a_,Lorentz],___]:> a);
(* non-Lorentz indices of the field *)
vectorindices=(ls/.Field[_,_,i_,___]:> i)/.Index[in_,_]:> in;
pvectorindices=Pattern[#,_]&/@vectorindices;

(* on the rs, drop all things that are just constant*vector field *)
rslist=List@@(rs+Nothing);

(* turn the vectors in this list into FS objects *)
rslistV = (
	Total[
		Cases[rslist, (c_:1)(x:Operator@Field[_, _Vector,___])/;FreeQ[c, Field|FieldStrength]]
		] /.
	Index[vectormu,Lorentz]:> mu) /. Field[lab_,Vector[mu],in_,_] :> FieldStrength[lab,{Index[mu,Lorentz],Index[nu,Lorentz]},in,{}];
(* currents that we will act on with the CD *)
rslistJ=OperatorToNormalForm[RelabelIndices[Total[DeleteCases[rslist, (c_:1)(x:Operator[Field[_, _Vector,___]])/;FreeQ[c, Field|FieldStrength]]],Unique->True], Unique->True];


jmu= rslistJ/.Index[vectormu,Lorentz]:> Index[mu,Lorentz];
jnu= rslistJ/.Index[vectormu,Lorentz]:> Index[nu,Lorentz];
outFS=OperatorToNormalForm[rslistV+CD[mu,jnu]-CD[nu,jmu], Unique->True];

With[{trhs = outFS},
	outFS=FieldStrength[vectorlabel,{Index[pmu,Lorentz],Index[pnu,Lorentz]},pvectorindices,Pattern[drvs,_]]:> RelabelIndices[CD[drvs,trhs],Unique->True]/.Inactive[CD]-> CD
];

{outFS,outV}
]


(* ::Subsubsection::Closed:: *)
(*ShiftVectorFields*)


(*ShiftVectorFields::usage="ShiftVectorFields[L, V, dV] redefines (non-gauge) vector fields V inside L by dV.";*)
ShiftVectorFields::shiftmismatch="The fields list '`1`' and the list of shifts '`2`' are not of the same length.";
ShiftVectorFields::gaugeinvariance= "The vector field shifts contain at least one gauge term, which is currently impossible.";

ShiftVectorFields[expr_,fields_List,shift_List]:=Module[
	{pairs,pairsA,pairsV,
	A, V, dAA, dAV, dVA, dVV,temp,
	SplitVectorShift,GaugeContributionQ,
	chargeShiftMatrix
	},

	(* make sure shapes are compatible *)
	If[Length@fields =!= Length@shift,
	Message[ShiftVectorFields::shiftmismatch,fields,shift];Return[expr]];

	(* split fields and their shifts into gauge and non-gauge *)
	pairs=Transpose@{fields,Operator/@(Expand@shift)};
	pairsA=Cases[pairs,{x_,y_}/;GaugeFieldQ@x];
	pairsV=Cases[pairs,{x_,y_}/;!GaugeFieldQ@x];

	If[pairsV==={}, pairsV={{},{}},pairsV=Transpose@pairsV];
	If[pairsA==={}, pairsA={{},{}},pairsA=Transpose@pairsA];

	V=pairsV[[1]];
	A=pairsA[[1]];

	(* is a term a shift by a gauge field? *)
	GaugeContributionQ[x_]:=MatchQ[x, (c_:1)*Operator[Field[X_,FT_Vector,ARGS__]]/;GaugeFieldQ[Field[X, FT,ARGS]]];

	(* splitting shifts into gauge and nongauge. And yes, this is a module inside a module. *)
	SplitVectorShift[x_]:=Module[{shiftlist},
		shiftlist=List@@(x+Nothing);
		{Total@Cases[shiftlist,k_?GaugeContributionQ],Total@Cases[shiftlist,Except[k_?GaugeContributionQ]]}
	];


	temp=(SplitVectorShift/@pairsA[[2]]);
	If[temp=={}, temp={{},{}},temp=Transpose@temp];
	{dAA,dAV}=temp;

	temp=(SplitVectorShift/@pairsV[[2]]);
	If[temp=={}, temp={{},{}},temp=Transpose@temp];
	{dVA,dVV}=temp;

	(* make sure any doubly occuring indices are replaced by unique dummies *)
	dVV=RelabelIndices[dVV, Unique->True];
	dVA=RelabelIndices[dVA, Unique->True];
	dAA=RelabelIndices[dAA, Unique->True];
	dAV=RelabelIndices[dAV, Unique->True];

	(* check that there is no shift of a proca field involving a gauge field *)
	If[Length@Cases[Flatten[dVA], x_/;x=!=0] >0,
		Message[ShiftVectorFields::gaugeinvariance];
		Return[expr];
	];

	(* now we need to remove EoM operators, so the shifts trigger on them, but still need to deactivate Mathematica's x*x -> x\.b2 rule *)
	temp=RemovePower@(RelabelIndices[OperatorToNormalForm[expr, Unique ->True, CanonizeKinetic->True],Unique->True])(*/.(a1_:1) Power[b1_/;(!FreeQ[b1, Field|FieldStrength]), k1_Integer?Positive]:>a1 Inactive[Times]@@ConstantArray[b1,k1]*);

	(* set up replacement rules for the non-gauge vectors, meaning V -> dVV *)
	temp=temp/.Flatten[VectorFieldReplacement/@({V,dVV}\[Transpose])];

	(* now we expand the covariant derivatives everywhere, so we have gauge vectors explicit and can replace them by proca fields A -> dAV *)
	temp=temp/.{
		Field[l_,t_,i_,d_List/;Length[d]>0]:> (CDExp[d,Field[l,t,i,{}]]),
		FieldStrength[l_,t_,i_,d_List/;Length[d]>0]:> (CDExp[d,FieldStrength[l,t,i,{}]])
		};
	(* now we replace the gauge fields *)
	temp=temp/.Flatten[VectorFieldReplacement/@({A,dAV+dAA}\[Transpose])];
	(*temp=temp/.(Map[Inactive@RelabelIndices[#,Unique->True]&,Association@Flatten[VectorFieldReplacement/@({A,dAV+dAA}\[Transpose])],1] /.Association->List );*)

	temp=temp//Activate;
	(* now we remove the gauge fields explicitly showing since they are a relic of our method *)
	temp=temp/.Field[Alternatives@@(GetGaugeGroups[#][Field]&/@Keys@GetGaugeGroups[]), ___]:>0;
	temp

]


(* ::Subsection:: *)
(*Field redefinitions*)


(* ::Subsubsection::Closed:: *)
(*HBarExpand*)


(* quickly handles large expressions expanding them step by step - which is somehow faster than Mathematica's in-house solution "Expand" *)
HBarExpand[expr_Plus]:=HBarExpand/@ expr


HBarExpand[expr_Times]:=Module[{reexp, t=1,inc},
	reexp =List@@(HBarExpand /@  expr);
	For[inc=1,inc<=Length@reexp,inc++,
		t =Expand[ t * reexp[[inc]]]
	];
	t
]


HBarExpand[Except[_Times|_Plus,x_]]:=x


(* ::Subsubsection::Closed:: *)
(*CoefficientOperator*)


CoefficientOperator[0,_] = 0;
CoefficientOperator[a_Plus,b_]:=(CoefficientOperator[#,b]& /@ a)
CoefficientOperator[a_ n_^(m_/;m<0),b_]:=(CoefficientOperator[a,b]n^m)
CoefficientOperator[a_,b_]/;(Head[a]=!= Plus&&(!FreeQ[a, Plus])):=(CoefficientOperator[Expand[a],b])

CoefficientOperator[(cin_:1)oin_Operator, b_]:=Module[{c,o,internalize,bPattern, rule={}, repcounter=0, closedChains},

	c=cin;
	o=oin;

	If[(LOpenSpinChainQ@ OperatorToNormalForm[o] && LOpenSpinChainQ@ b)||(ROpenSpinChainQ@ OperatorToNormalForm[o] && ROpenSpinChainQ@ b),
		o = Operator[OperatorToNormalForm[o] //.x_?ClosedSpinChainQ:>(c=c*x;1)]
	];

	internalize[m_]:=Symbol[SymbolName[m]<>"int"];
	internalize[Index[m_, t_]]:=Index[internalize[m],t];

	bPattern=b/. Index[mu_, t_]:> Index[Pattern[Evaluate@internalize[mu], _],t];

	AppendTo[rule,RuleDelayed[Condition[bPattern,repcounter++==0],Evaluate[(Times@@Table[Delta[k,internalize[k]],{k,FindOpenIndices[b]}])]]];

	Operator[ ContractDelta[OperatorToNormalForm[If[FreeQ[o,bPattern],0,c*ReplaceAll[o,rule]], Unique -> True]]]
]


(* ::Subsubsection::Closed:: *)
(*ReduceField*)


NonGaugeIndicesQ[f_]:=Module[{ind},
	ind = GetFields[f][Indices];

	If[Length@Cases[GroupFromRep/@ind,None]>0,
		True
	,
		False
	]
]


ReduceField::unimplemented = "The field type of '`1`' is currently not supported. Returning the input Lagrangian.";


(* redefines the field f such that reducible operators proportional to EoM[f] are removed from L *)
Options[ReduceField]={ShiftOrder -> All, ResultOrder->All};

ReduceField[L_,f_List/;Length[f]>1,opt:OptionsPattern[]]:=ReduceField[ReduceField[L, First@f, opt], f[[2;;]] , opt]

ReduceField[L_,{f_Symbol},opt:OptionsPattern[]]:=ReduceField[L,f,opt]

ReduceField[iL_, f_Symbol, OptionsPattern[]]:=Module[{real, type,oShift,oRes, res},
	type = GetFields[][f][Type];
	real = GetFields[][f][SelfConjugate];

	(* handle arguments *)
	oShift=If[MatchQ[OptionValue[ShiftOrder],_Integer],OptionValue[ShiftOrder],GetMaxOrder[iL]];
	oRes=If[MatchQ[OptionValue[ResultOrder],_Integer],OptionValue[ResultOrder],GetMaxOrder[iL]];

	$MonitorString2="-> Reducing "<> ToString@f<> " in terms of dimension "<>ToString@OptionValue@ShiftOrder <>".";

	Switch[{type,real},
		{Scalar,True},
			(* call the workaround solution if there are non-gauge indices on the real scalar *)
			If[Length@Cases[GroupFromRep/@(GetFields[f][Indices]), None|(x_/;MemberQ[Keys@GetGlobalGroups[],x])]>0,
				         res= ReduceRealScalarTemp[iL,f,oShift,oRes]
			,
				         res= ReduceRealScalar[iL,f,oShift,oRes]
			],
		{Scalar,False},  res= ReduceComplexScalar[iL,f, oShift, oRes],
		{Fermion,False}, res= ReduceDiracFermion[iL,f, oShift, oRes],
		{Fermion,True},  res= ReduceMajoranaFermion[iL,f, oShift,oRes],
		{Vector,True},   res= ReduceRealVector[iL, f, oShift, oRes],
		_, Message[ReduceField::unimplemented, f]; res= iL
	];

	CollectOperators[res, NormalForm->False]
]


(* ::Subsubsection::Closed:: *)
(*ReduceRealScalar*)


ReduceRealScalar[iL_, f_Symbol, oShift_,oRes_]:=
Module[
	{LNoShift=0,L,Ltemp,LR,inds,pinds, rule0,rule1, shift,shift0,shift1, field, fieldPattern0,fieldPattern1, drvs,Lread, a,b,c,ri, LAMBDA}
,
	inds=Symbol["i"<>ToString[#]]&/@Range@Length@GetFields[][f][Indices];
	pinds=Pattern[#, _]&/@inds;

	(* prepare the Lagrangian for reading off shift: expand and remove the FreeLag of the field we are reducing *)
	L=SeriesEFT[DropFreeQ[iL,f]-InternalSimplify@FreeLag[f],EFTOrder->oShift];
	(* need to relabel indices so the dummy indices don't lead to problems when substituting *)
	shift= RelabelIndices[OperatorToNormalForm[CoefficientOperator[L, EoM[f[Sequence@@inds]]], Unique -> True],Unique->True];
	shift0 = Coefficient[shift, hbar, 0];
	shift1 = Coefficient[shift, hbar, 1];

	(* set up the redefinition *)
	field=First@Cases[{f[Sequence@@inds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[a,b,c,drvs];
	fieldPattern0=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[{a, oShift-4,0},b,c,Pattern[drvs,_]];
	fieldPattern1=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[{a, oShift-4,1},b,c,Pattern[drvs,_]];

	With[{temp0=shift0,temp1=shift1,tfield=field,tdrvs=drvs},
		rule0 = (fieldPattern0 :> RelabelIndices[CD[tdrvs, temp0],Unique->True]);
		rule1 = (fieldPattern1 :> RelabelIndices[CD[tdrvs, temp1],Unique->True])
	];

	(* split the Lagrangian into terms affected and unaffected by the shift *)
	{LR,LNoShift}=SeparateByFieldContent[SplitLagrangianByPower[iL,4+oRes-oShift],f];

	(* insert the field expansion, return to NormalForm since we're inserting fields into EoM objects as well *)
	LR = SeriesEFT[BetterExpand @ OperatorToNormalForm[ LR , Unique->True, CanonizeKinetic->True] /. Field[f,args__]:>Field[f,args] + Field[{f, oShift-4,0},args] + hbar Field[{f, oShift-4,1},args],EFTOrder->oRes];

	(* we need to prepare the Lagrangian: powers need to be removed but we cannot use operators, because we need EoM-type objects to stay unchanged *)
	LR = RemovePower@(RelabelIndices[LR,Unique->True])(*/.(a1_:1) Power[b1_/;(!FreeQ[b1, Field|FieldStrength]), k1_Integer?Positive]:>a1 Inactive[Times]@@ConstantArray[b1,k1]*);

	(* plug in the rule, reactivate the powers and IBPSimplify *)
	LR = InternalSimplify @ Activate[LR/.{rule0,rule1}];

	LR+LNoShift
]


(* ::Subsubsection::Closed:: *)
(*ReduceRealScalarTemp*)


(* ::Text:: *)
(*Special implementation for real scalars in which the redundant operator is dropped without needing it to cancel against field redefinitions. Useful for scalars with flavor indices, where GreensSimplify fails to identify the result of the shift with the redundant operators and gets stuck.*)


ReduceRealScalarTemp[iL_, f_Symbol, oShift_,oRes_]:=
Module[
	{LNoShift=0,L,Ltemp,LR,inds,pinds, rule0,rule1, shift,shift0,shift1, field, fieldPattern0,fieldPattern1, drvs,Lread, a,b,c,ri, LAMBDA, Lkin, LDelta, Lred,\[Xi],Lkin\[Chi], LDelta\[Chi], Lred\[Chi]}
,
	inds=Symbol["i"<>ToString[#]]&/@Range@Length@GetFields[][f][Indices];
	pinds=Pattern[#, _]&/@inds;

	(* prepare the Lagrangian for reading off shift: expand and remove the FreeLag of the field we are reducing *)
	L=SeriesEFT[DropFreeQ[iL,f]-InternalSimplify@FreeLag[f],EFTOrder->oShift];
	(* need to relabel indices so the dummy indices don't lead to problems when substituting *)
	shift= RelabelIndices[OperatorToNormalForm[CoefficientOperator[L, EoM[f[Sequence@@inds]]], Unique -> True],Unique->True];
	shift0 = Coefficient[shift, hbar, 0];
	shift1 = Coefficient[shift, hbar, 1];

	(* set up the left-hand side of the field-redefinition *)
	field=First@Cases[{f[Sequence@@inds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[a,b,c,drvs];
	fieldPattern0=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[{a, oShift-4,0},b,c,Pattern[drvs,_]];
	fieldPattern1=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[{a, oShift-4,1},b,c,Pattern[drvs,_]];

	(* set up the right-hand side of the field-redefinition *)
	With[{temp0=shift0,temp1=shift1,tdrvs=drvs},
		rule0 = (fieldPattern0 :> RelabelIndices[CD[tdrvs, temp0],Unique->True]);
		rule1 = (fieldPattern1 :> RelabelIndices[CD[tdrvs, temp1],Unique->True])
	];

	(* split the Lagrangian into terms affected and unaffected by the shift *)
	(*{LR,LNoShift} = SeparateByFieldContent[SplitLagrangianByPower[iL,4+oRes-oShift],f];*)

	(* split Lagrangian into the kinetic term, the redundant term and the residual *)
	Lkin = OperatorToNormalForm[ KineticTerms@ FreeLag[f], Unique->True, CanonizeKinetic->True];
	Lred = shift CD[{\[Xi],\[Xi]},f[Sequence@@inds]];
	LDelta = InternalSimplify[ iL -  Lred - Lkin ];

	(* output the kinetic term and its shift but without the linear term that would only cancel the redundant operator *)
	Lkin\[Chi]=SeriesEFT[InternalSimplify[Lkin+1/2 CD[\[Xi],RelabelIndices[shift,Unique->True]]CD[\[Xi],RelabelIndices[shift,Unique->True]]], EFTOrder->oRes];

	(* next is the term that comes from redefining the redundant operator *)
	Lred\[Chi] = RelabelIndices@ SeriesEFT[
								((shift /. Field[f,args__]:> Field[f,args] + Field[{f, oShift-4,0},args] + hbar Field[{f, oShift-4,1},args])-shift)CD[{\[Xi],\[Xi]},f[Sequence@@inds]]
								+(shift /. Field[f,args__]:>Field[f,args] + Field[{f, oShift-4,0},args] + hbar Field[{f, oShift-4,1},args])CD[{\[Xi],\[Xi]},shift]
							,EFTOrder->oRes];
	Lred\[Chi] = RemovePower@(RelabelIndices[Lred\[Chi],Unique->True])(*/.(a1_:1) Power[b1_/;(!FreeQ[b1, Field|FieldStrength]), k1_Integer?Positive]:>a1 Inactive[Times]@@ConstantArray[b1,k1]*);
	Lred\[Chi] = InternalSimplify @ Activate[Lred\[Chi]/.{rule0,rule1}];

	(* now we run the full replacement routine on all the other terms in the Lagrangian *)
	{LR,LNoShift} = SeparateByFieldContent[SplitLagrangianByPower[LDelta,4+oRes-oShift],f];

	(* insert the field expansion, return to NormalForm since we're inserting fields into EoM objects as well *)
	LR = SeriesEFT[BetterExpand @ OperatorToNormalForm[ LR, Unique->True, CanonizeKinetic->True ] /. Field[f,args__]:>Field[f,args] + Field[{f, oShift-4,0},args] + hbar Field[{f, oShift-4,1},args],EFTOrder->oRes];

	(* we need to prepare the Lagrangian: powers need to be removed but we cannot use operators, because we need EoM-type objects to stay unchanged *)
	LR = RemovePower@(RelabelIndices[LR,Unique->True])(*/.(a1_:1) Power[b1_/;(!FreeQ[b1, Field|FieldStrength]), k1_Integer?Positive]:>a1 Inactive[Times]@@ConstantArray[b1,k1]*);

	(* plug in the rule, reactivate the powers and IBPSimplify *)
	LR = InternalSimplify @ Activate[LR/.{rule0,rule1}];

	Total[{LR,Lred\[Chi],Lkin\[Chi],LNoShift}]
]


(* ::Subsubsection::Closed:: *)
(*ReduceComplexScalar*)


ReduceComplexScalar[iL_, f_Symbol, oShift_,oRes_]:=Module[
	{L,LNoShift,Ltemp,LR,inds,pinds,  chi1, chi2, shift,
	shift0,rule0,shift1,rule1,LR0,LR1,
	field, fieldPattern0,fieldPattern1, drvs,Lread,a,b,c,ri},

	inds=Symbol["i"<>ToString[#]]&/@Range@Length@GetFields[][f][Indices];
	pinds=Pattern[#, _]&/@inds;

	(* expand Lagrangian to ShiftOrder, keep only terms depending on the field f and IBPSimplify them *)
	L=SeriesEFT[DropFreeQ[iL,f]-InternalSimplify@FreeLag[f],EFTOrder->oShift];

	(* factor out the EoM[f] term *)
	chi1=RelabelIndices[OperatorToNormalForm[CoefficientOperator[L, EoM[f[Sequence@@inds]]],Unique -> True],Unique->True];

	(* ...and subtract the terms we found, to now get the EoM[Bar@f] coefficients *)
	chi2=RelabelIndices[CoefficientOperator[RelabelIndices@(L - Operator[chi1 EoM[f[Sequence@@inds]]]),EoM[Bar@f[Sequence@@inds]]],Unique->True];


	shift=1/2 Bar@OperatorToNormalForm[(chi1+Bar@chi2),Unique -> True];
	shift0 = Coefficient[shift, hbar, 0];
	shift1 = Coefficient[shift, hbar, 1];

	field=First@Cases[{f[Sequence@@inds]},_Field,Infinity] /. Field[a_,b_,c_,{}] :> Field[a,b,c,drvs];
	fieldPattern0=First@Cases[{f[Sequence@@pinds]},_Field,Infinity] /. Field[a_,b_,c_,{}] :> Field[{a,oShift-4,0},b,c,Pattern[drvs,_]];
	fieldPattern1=First@Cases[{f[Sequence@@pinds]},_Field,Infinity] /. Field[a_,b_,c_,{}] :> Field[{a,oShift-4,1},b,c,Pattern[drvs,_]];

	With[{temp0=shift0,temp1=shift1,tfield=field,tdrvs=drvs},
		rule0 = (fieldPattern0 :> RelabelIndices[CD[tdrvs, temp0],Unique->True]);
		rule1 = (fieldPattern1 :> RelabelIndices[CD[tdrvs, temp1],Unique->True])
	];

	(* split the Lagrangian into terms affected and unaffected by the shift *)
	{LR,LNoShift}=SeparateByFieldContent[SplitLagrangianByPower[iL,4+oRes-oShift],f];

	(* we now expand the field and then the Lagrangian *)
	LR = SeriesEFT[BetterExpand @ OperatorToNormalForm[ LR, Unique->True, CanonizeKinetic->True] /. Field[f, args__] :> Field[f, args ] + Field[{f,oShift-4,0},args] + hbar Field[{f,oShift-4,1},args],EFTOrder->oRes];


	(* we need to prepare the Lagrangian: powers need to be removed but we cannot use operators, because we need EoM-type objects to stay unchanged *)
	LR = RemovePower@(RelabelIndices[LR,Unique->True])(*/.(a1_:1) Power[b1_/;(!FreeQ[b1, Field|FieldStrength]), k1_Integer?Positive]:>a1 Inactive[Times]@@ConstantArray[b1,k1]*);

	LR = Activate[LR /. {rule0,rule1}];

	LR=InternalSimplify @ LR;

	LR+LNoShift
]


(* ::Subsubsection::Closed:: *)
(*ReduceMajoranaFermion*)


ReduceMajoranaFermion[iL_, f_Symbol, oShift_, oRes_]:=Module[{start=Now,LNoShift=0,L,LInt,Ltemp,LR,inds,pinds, rule0,rule1, chi1, chi2, shift,shift0,shift1, field, fieldPattern0,fieldPattern1, drvs,Lread, a,b,c,ri, LAMBDA},

	inds=Symbol["i"<>ToString[#]]&/@Range@Length@GetFields[][f][Indices];
	pinds=Pattern[#, _]&/@inds;

	(* expand Lagrangian to ShiftOrder, keep only terms depending on the field f and IBPSimplify them *)
	L=SeriesEFT[DropFreeQ[iL,f]-InternalSimplify@FreeLag[f],EFTOrder->oShift];

	(* factor out the EoM[f] term *)
	chi1=RelabelIndices[OperatorToNormalForm[CoefficientOperator[L, EoM[First@Cases[{f[Sequence@@inds]},_Field,Infinity]]],Unique->True],Unique->True];

	(* ...and subtract the terms we found, to now get the EoM[Bar@f] coefficients *)
	chi2=RelabelIndices[OperatorToNormalForm[CoefficientOperator[RelabelIndices[L-Operator[chi1**EoM[First@Cases[{f[Sequence@@inds]},_Field,Infinity]]]],EoM[First@Cases[{Transp@f[Sequence@@inds]},Transp[_Field],Infinity]]],Unique->True],Unique->True];

	shift=RelabelIndices[I*CC**(Transp@chi1-chi2)];
	shift0 = Coefficient[shift, hbar, 0];
	shift1 = Coefficient[shift, hbar, 1];

	field=First@Cases[{f[Sequence@@inds]},_Field,Infinity]/. Field[a_,b_,c_,{}]-> Field[a,b,c,drvs];
	fieldPattern0=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]-> Field[{a,oShift-4,0},b,c,Pattern[drvs,_]];
	fieldPattern1=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]-> Field[{a,oShift-4,1},b,c,Pattern[drvs,_]];

	With[{temp0=shift0,temp1=shift1,tfield=field,tdrvs=drvs},
		rule0 = (fieldPattern0 :> RelabelIndices[CD[tdrvs, temp0],Unique->True]);
		rule1 = (fieldPattern1 :> RelabelIndices[CD[tdrvs, temp1],Unique->True])
	];

	(* the shift should be homogeneous in power-counting, so we split the Lagrangian into two pieces *)
	{LR,LNoShift}=SplitLagrangianByPower[iL,4+oRes-oShift];

	LR = SeriesEFT[ OperatorToNormalForm[ LR , Unique->True, CanonizeKinetic->True] /. Field[f, args__] :> Field[f, args] + Field[{f, oShift-4,0}, args]+hbar Field[{f, oShift-4,1}, args] , EFTOrder->oRes];

	(* we need to prepare the Lagrangian: powers need to be removed but we cannot use operators, because we need EoM-type objects to stay unchanged *)
	LR = RemovePower@(RelabelIndices[LR,Unique->True])(*/.(a1_:1) Power[b1_/;(!FreeQ[b1, Field|FieldStrength]), k1_Integer?Positive]:>a1 Inactive[Times]@@ConstantArray[b1,k1]*);

	LR=Activate[LR/.{rule0,rule1}];

	LR=InternalSimplify @ LR;

	LR+LNoShift
]


(* ::Subsubsection::Closed:: *)
(*Anomaly Factors*)


AnomalyFactor[field_]:=Module[{grAbel,grNAbel,dims,charges,tempInd,tI,tB,tJ,tA,t\[Mu],t\[Nu],t\[Rho],t\[Sigma]},
	grAbel = Head/@ GetFields[field, Charges];
	grNAbel = DeleteDuplicates@ DeleteCases[GroupFromRep/@GetFields[field ,Indices],None];

	(* dimensionality of each represenation including a 1 for the abelian *)
	dims = Join[ConstantArray[1,Length@grAbel],Table[Delta[#,#]&@Index[tempInd,FirstCase[GetFields[field, Indices], _gr]],{gr,grNAbel}]];
	(* anomaly trace factors for each group *)
	charges = Join[
		(GetGaugeGroups[#,Coupling]^2 FS[GetGaugeGroups[#,Field],t\[Mu],t\[Nu]]FS[GetGaugeGroups[#,Field],t\[Rho],t\[Sigma]]FirstCase[GetFields[field,Charges], #[q_]:>q]^2)&/@ grAbel,
		(GetGaugeGroups[#,Coupling]^2 FS[GetGaugeGroups[#,Field],t\[Mu],t\[Nu],tA]FS[GetGaugeGroups[#,Field],t\[Rho],t\[Sigma],tB]FieldGenerators[field,#,{tA,tI,tJ}]FieldGenerators[field,#,{tB,tJ,tI}])&/@grNAbel
	];

	Return[-2 I hbar RelabelIndices[Sum[Part[charges,inc]Times@@Drop[dims,{inc}],{inc,Length@dims}],Unique->True]]
]


(* ::Subsubsection::Closed:: *)
(*ReduceDiracFermion*)


(* check if a field transforms under any of the gauge groups *)
FieldChargedQ[f_Symbol] := (Length@Join[DeleteCases[GroupFromRep/@GetFields[f,Indices],None],GetFields[f,Charges]] > 0)


ReduceDiracFermion[iL_, f_Symbol, oShift_, oRes_]:=Module[
	{L,LInt,Ltemp,LR,inds,pinds, rule0,rule1, chi1, chi2, shift,shift0,shift1, field, fieldPattern0,fieldPattern1, drvs,Lread, a,b,c,ri, LAMBDA,
	LNoShift,g5mark,g5coeff, LAnomaly=0,abelgroups,nabelgroups},

	inds=Symbol["i"<>ToString[#]]&/@Range@Length@GetFields[][f][Indices];
	pinds=Pattern[#, _]&/@inds;

	(* expand Lagrangian to ShiftOrder, keep only terms depending on the field f and IBPSimplify them *)
	L=SeriesEFT[DropFreeQ[iL,f]-InternalSimplify@FreeLag[f],EFTOrder->oShift];

		(* factor out the EoM[f] term *)
	chi1=RelabelIndices[OperatorToNormalForm[CoefficientOperator[L , EoM[First@Cases[{f[Sequence@@inds]},_Field,Infinity]]], Unique->True],Unique->True];
	(* ...and subtract the terms we found, to now get the EoM[Bar@f] coefficients *)
	chi2=RelabelIndices[
			CoefficientOperator[
				RelabelIndices[L-Operator[chi1**EoM[First@Cases[{f[Sequence@@inds]},_Field,Infinity]]]],
				EoM[First@Cases[{Bar@f[Sequence@@inds]},Bar[_Field],Infinity]]],
			Unique->True];

	shift=-I/2*Bar@OperatorToNormalForm[( chi1+ Bar@chi2),Unique->True];
	shift0 = Coefficient[shift, hbar, 0];
	shift1 = Coefficient[shift, hbar, 1];

	field=First@Cases[{f[Sequence@@inds]},_Field,Infinity]/. Field[a_,b_,c_,{}]-> Field[a,b,c,drvs];
	fieldPattern0=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[{a,oShift-4,0},b,c,Pattern[drvs,_]];
	fieldPattern1=First@Cases[{f[Sequence@@pinds]},_Field,Infinity]/. Field[a_,b_,c_,{}]:> Field[{a,oShift-4,1},b,c,Pattern[drvs,_]];

	With[{temp0=shift0,temp1=shift1,tfield=field,tdrvs=drvs},
		rule0 = (fieldPattern0 :> RelabelIndices[CD[tdrvs, temp0],Unique->True]);
		rule1 = (fieldPattern1 :> RelabelIndices[CD[tdrvs, temp1],Unique->True])
	];

	If[FieldChargedQ@ f && False,
		(* check for a gamma5 in a shift *)
		g5coeff= RelabelIndices@OperatorToNormalForm[Coefficient[
			CoefficientOperator[
				Operator[shift0],
				First@Cases[{f[Sequence@@inds]},_Field,Infinity]
			  ]/.NonCommutativeMultiply[DiracProduct[Proj[chi_]]]:>chi/2*g5mark, g5mark
			],Unique->True, CanonizeKinetic->True];
		(* generate the chiral Jacobian term *)
		If[!(g5coeff===0),
			(* groups the field transforms under *)
			LAnomaly = g5coeff AnomalyFactor[f]
		];
	];

	{LR,LNoShift}=SeparateByFieldContent[SplitLagrangianByPower[iL,4+oRes-oShift],f];

	LR = SeriesEFT[BetterExpand @ OperatorToNormalForm[ LR ,Unique->True, CanonizeKinetic->True]/. Field[f, args__] :> Field[f,args] + Field[{f, oShift-4,0}, args] +hbar Field[{f, oShift-4, 1}, args] , EFTOrder-> oRes];

	(* we need to prepare the Lagrangian: powers need to be removed but we cannot use operators, because we need EoM-type objects to stay unchanged *)
	LR=RemovePower@(RelabelIndices[LR,Unique->True])(*/.(a1_:1) Power[b1_/;(!FreeQ[b1, Field|FieldStrength]), k1_Integer?Positive]:>a1 Inactive[Times]@@ConstantArray[b1,k1]*);

	LR= Activate[LR/.{rule0,rule1}];

	LR= InternalSimplify @  LR;

	LR+LNoShift
]


(* ::Subsubsection::Closed:: *)
(*ReduceRealVector*)


(* extracts a list of all vector fields that mix kinetically with the field f *)
GetMixingVectorFields[L_,f_]:=Module[{kinT},
	(*get all kinetic terms involving f*)
	kinT=Total@DeleteCases[List@@Expand[OperatorToNormalForm[KineticTerms@L,Unique->True, CanonizeKinetic->True] + Nothing], k_/;FreeQ[k,f]];
	DeleteDuplicates@Join[{f},DeleteDuplicatesBy[Cases[kinT,FieldStrength[___],Infinity],First]/.FieldStrength[k_,__]:>k]
]


GetZMatrix::NotDiagonalizable="The kinetic terms are not canoncially normalized at leading order and cannot be diagonalized.";

GetZMatrix[L_,f_List]:=Module[{Z,emptyList,\[Lambda],\[Mu],\[Nu],GetMaxMassDimension},
	Z=(#+Transpose@#)&@UpperTriangularize@Table[-2Coefficient[KineticTerms@L,Operator[RelabelIndices[FS[i,\[Mu],\[Nu]]FS[j,\[Mu],\[Nu]]]]],{i,f},{j,f}];

	(* take the matrix at leading order, make sure it's identity *)
	emptyList=DeleteCases[Flatten[IdentityMatrix[Length@f]-(Map[SeriesEFT[# , EFTOrder->0]&,Z,{2}]/.hbar -> 0)],0];
	If[Length@emptyList>0,
		Message[GetZMatrix::NotDiagonalizable];
		Return[IdentityMatrix[Length@f]]
	];


	(* the kinetic matrix can be inverted by perturbative methods. check if mixing occurs only at one-loop order *)

	If[Length@DeleteCases[Flatten[IdentityMatrix[Length@f]-(Z/.hbar -> 0)],0]>0,
		(* no, we have tree-level subleading power mixing *)
		Return[(BiPerturbativeInvertRegular[Map[WriteEFTLambda[#,\[Lambda]]&,Z,{2}],{hbar,0,1},{\[Lambda],0,GetMaxOrder[L]-4}]/.\[Lambda]->1)]
	,
		(* kinetic mixing only at loop-level, inversion is trivial *)
		Return[(2IdentityMatrix[Length@f]-Z)]
	]


]


ReduceRealVector[iL_, f_Symbol, oShift_, oRes_]:=Module[
	{LNoShift,L,Ltemp,LR,inds,pinds, rule, shift, field, fieldPattern, drvs,Lread, LAMBDA, mixingFields,Z},

	inds=Symbol["i"<>ToString[#]]&/@Range@(Length@GetFields[][f][Indices]+1);
	pinds=Pattern[#, _]&/@inds;

	L=iL;

	(* find if the the field f mixes with any other field *)
	mixingFields=GetMixingVectorFields[L,f];

	If[mixingFields==={f},
		(* f does not mix with any other field -> straightforward redefinition *)
		shift= -RelabelIndices[OperatorToNormalForm[CoefficientOperator[(SeriesEFT[DropFreeQ[L,f],EFTOrder->oShift]), EoM[f[Sequence@@inds]]],Unique->True, CanonizeKinetic->True],Unique->True];

		field=First@Cases[{f[Sequence@@inds]},_Field,Infinity];

		If[OptionValue[ShiftOrder]===All,
				{LR,LNoShift}={L,0}
			,

				{LR,LNoShift}=SplitLagrangianByPower[L,4+oRes-oShift];
		];
		LR=(*InternalSimplify@*)RelabelIndices@Activate@ShiftVectorFields[OperatorToNormalForm[LR, Unique->True, CanonizeKinetic->True], {field},{field+shift}]
	,

		(* there is mixing, so we need to work a bit harder, first extract the rotation matrix *)
		Z = GetZMatrix[L,mixingFields];

		(* the shift is now vector-valued, with the standard shift in the 1 component and then rotated by Z *)
		shift=Z . (-(SparseArray[{1->1},{Length@mixingFields}])RelabelIndices[OperatorToNormalForm[CoefficientOperator[Operator@(SeriesEFT[DropFreeQ[L,mixingFields],EFTOrder->oShift]), EoM[f[Sequence@@inds]]], Unique->True],Unique->True]);
		field=(First@Cases[{#[Sequence@@inds]},_Field,Infinity])&/@mixingFields;


		{LR,LNoShift}=SplitLagrangianByPower[L,4+oRes-oShift];


		LR=RelabelIndices@OperatorToNormalForm[Activate@ShiftVectorFields[OperatorToNormalForm[LR, Unique->True, CanonizeKinetic->True], field,field+shift], Unique->True, CanonizeKinetic->True]
	];

	LR=InternalSimplify@SeriesEFT[HBarExpand @ LR, EFTOrder->oRes];
	LR+LNoShift
]


(* ::Subsection:: *)
(*EOMSimplify*)


(* ::Subsubsection::Closed:: *)
(*Check for a standard kinetic term present*)


EOMInvalidQ[L_]:=Module[{fieldsInL, freeLag, kinLag},
	fieldsInL = DeleteDuplicates[Join[
		Cases[L, _Field, Infinity]/.Field[f_,__]:>f,
		Cases[L, _FieldStrength, Infinity]/.FieldStrength[f_,__]:>f
		]];
	(* free Lagrangian of all fields appearing in L *)
	freeLag=Sum[FreeLag[field],{field,fieldsInL}];
	(* subtract freeLag from the kinetic terms in L at tree-level and at dimension 4 *)
	kinLag = SeriesEFT[InternalSimplify@KineticTerms[freeLag - L]/.hbar -> 0, EFTOrder->4];
	Return[!(kinLag === 0)]
]


ContainsHeavyEOMQ[L_]:=Module[{EOMFields},
	(* occurences of EOM objects and which field labels they are *)
	EOMFields = DeleteDuplicates@Cases[L,EoM[Field[f_,___]]|EoM[Bar@Field[f_,___]]|EoM[Transp@Field[f_,___]]->f,Infinity];
	(* check if they are registered as heavy *)
	EOMFields = (Heavy//$FieldAssociation@#)&/@EOMFields;
	(* return true if one of them is true *)
	Return[Or@@EOMFields]
]


(* ::Subsubsection::Closed:: *)
(*FieldsToShift*)


(* return the label of a field that needs to be redefined along with the lowest order in EFT counting at which the shift has to be performed *)
FieldsToShift[L0_]:=Module[{freeL,L,fields, list},
	(* all fields appearing in the input Lagrangian *)
	fields=DeleteDuplicates@Cases[L0, HoldPattern@Field[a_,__]:>a,Infinity];
	(* the free Lagrangian of all of these fields - so we get non-standard kinetic terms *)
	freeL = InternalSimplify@Total[FreeLag[#]&/@fields];
	L = Collect[L0 - freeL,_Operator];

	(* vector fields get a special treatment, so their kinetic terms have to be removed completely *)
	L = L - SelectOperatorClass[L,{},4];

	(* now extract all EoM terms *)
	list = Cases[TermsToList@L, x_/;MemberQ[x, EoM[___],{0,Infinity}]];
	list = SortBy[{(First@Cases[#, _EoM,{0,Infinity}]&/@list)/.EoM[Field[f_,___]]|EoM[Bar@Field[f_,___]]|EoM[Transp@Field[f_,___]]:>f,(OperatorDimension@OperatorToNormalForm@#)&/@list}\[Transpose],{First,Last}];
	(* ignore EOM operators of heavy fields *)
	list = Select[list,!(Heavy//$FieldAssociation@First[#])&];
	Return[SortBy[First/@Gather[SortBy[list,{First,Last}],First[#1]===First[#2]&],Last]]
]



(* ::Subsubsection::Closed:: *)
(*Main module*)


EOMSimplify::EffectiveCoupling = "The Lagrangian contains terms of lower power than dimension 4. Defining effective couplings and assuming these terms to be dimension 4. Use 'PrintEffectiveCouplings' and 'ReplaceEffectiveCouplings' to recover explicit expressions.";
EOMSimplify::InvalidLagrangian = "The Lagrangian contains no valid kinetic terms for the occuring Kinetic terms and thus cannot be simplified.";
EOMSimplify::HeavyFieldEOM = "The Lagrangian contains redundant operators of fields that are flagged to be integrated out. Operators of this kind are not reduced by EOMSimplify.";


Options[EOMSimplify] = {
		DummyCoefficients -> False,
		EFTOrder->All,
		EffectiveCouplingSymbol -> "C",
		ReductionIdentities->EvanescenceFree,
		Verbose -> True,
		"New" -> True
	};
Options[EoMSimplificationStep] = {EFTOrder -> Automatic};


EOMSimplify[Lagrangian_, OptionsPattern[]] (*? OptionsCheck*):=Module[{out,maxOrder,operatorList,L3,L4,La,L,consts,redID=OptionValue@ReductionIdentities,redIDIntermediate},
	$MonitorString1="Preparing field redefinitions...";
	$MonitorString2="Determining order of input Lagrangian...";
	$auxTime=AbsoluteTime[];

	If[redID === EvanescenceFree,
		AddToBibliography["EvanescentTreatment", "Simplified expression to evanescent-free scheme (with EOMSimplify)"]
	];

	(* if the user chooses to keep evanescent operators explicit, do intermediate steps without introducing them, otherwise we keep expanding and reintroducing them in every step *)
	If[redID === Evanescent, redIDIntermediate = dDimensional, redIDIntermediate = redID];

	(*subtract off constants*)
	OptionalMonitor[OptionValue@Verbose,
		If[FreeQ[Lagrangian,_Field,All]&&FreeQ[Lagrangian,_FieldStrength,All], Return[Lagrangian]];
		L= HcExpand@ BetterExpand@ Lagrangian;
		{consts, L}= SeparateOutConstants@ L;
	,"Checking for constant terms \[Ellipsis]"];

	If[EOMInvalidQ[L],Message[EOMSimplify::InvalidLagrangian]; Return[L]];

	operatorList=OperatorDimension/@TermsToList[L];

	If[Head@OptionValue[EFTOrder]===Integer, maxOrder = OptionValue@EFTOrder, maxOrder=Max[operatorList]];
	If[Min[operatorList]<4,
		(* there are superleading terms in the Lagrangian, redefine them *)
		OptionalMonitor[OptionValue@Verbose,
			La = SeriesEFT[SubstituteCoefficients[InternalSimplify[L], EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol]], EFTOrder -> maxOrder];
		,"Redefining superleading terms \[Ellipsis]"];
		If[TrueQ@OptionValue@Verbose, Echo[EOMSimplify::EffectiveCoupling]]; (* maybe we could also print the operators here ? *)
		,
		La = L;
	];

	If[OptionValue@DummyCoefficients === True,
		OptionalMonitor[OptionValue@Verbose,
			La = InternalSimplify[ IntroduceDummyCoefficients @ La, ReductionIdentities->redIDIntermediate ];
		,"Introducing dummy coefficients and simplifying the Lagrangian \[Ellipsis]"]
		,
		OptionalMonitor[OptionValue@Verbose,
			La = InternalSimplify[ La, ReductionIdentities->redIDIntermediate ];
		,"Simplifying the Lagrangian \[Ellipsis]"]
	];

	If[ContainsHeavyEOMQ[La], Message[EOMSimplify::HeavyFieldEOM]];
	If[OptionValue@ "New",
		out= PerformSystematicFieldRedefs[La, maxOrder, OptionValue@ Verbose];
	,
		out = OptionalMonitor[OptionValue@Verbose,FixedPoint[EoMSimplificationStep[#,EFTOrder->maxOrder]&, La],$MonitorString1<>"\n"<>$MonitorString2];
	];
	$MonitorString1="";
	$MonitorString2="";

	If[OptionValue@DummyCoefficients === True,
		OptionalMonitor[OptionValue@Verbose,
			out = RelabelIndices[ out , Unique -> True ] //. $IntCouplingRules;
			(* fast expansion *)
			out = Collect[out,hbar];
			out = FastExpand@out;
		,"Resubstituting dummy coefficients \[Ellipsis]"];

		(* due to some issue it is necessary to truncate the EFT series again if the option DummyCoefficients -> True is used *)
		out = SeriesEFT[out,EFTOrder->maxOrder]
	];

	OptionalMonitor[OptionValue@Verbose,
		out = consts + GreensSimplify[out, ReductionIdentities->redID];
	,"Simplifying the results \[Ellipsis]"];

	out
]


(* ::Subsubsection::Closed:: *)
(*Individual step of the main module*)


(* debug flag that shows the intermediate Lagrangians *)
$EOMSimplifyEchoIntermediateLagrangians = False;


EoMSimplificationStep[L_,OptionsPattern[]]:=Module[{task, fields,order,temp,maxorder},
	$MonitorString1="Determining redundant operators...";
	$MonitorString2="";
	maxorder=OptionValue[EFTOrder];
	(* look at all fields that are involved in redundant operators *)
	task=GatherBy[FieldsToShift[L],Last];
	(* if nothing to simplify -> break out *)
	If[Length[task]==0, Return[L]];

	(* the task list is sorted by mass dimension -> start at the lowest order *)
	task=Transpose@First@task;
	fields=First@task;
	order=First@Last@task;
	$MonitorString1="Field redefinitions needed for "<> ToString@fields<> " at dimension "<> ToString@order<>".";
	(* call ReduceField with the current list of fields *)
	temp=ReduceField[L,fields, ShiftOrder-> order,ResultOrder->maxorder];

	(*LCurrent=temp;*)
	$MonitorString2="Simplifying Lagrangian...";
	temp = InternalSimplify @ temp;

	If[TrueQ[$EOMSimplifyEchoIntermediateLagrangians], Echo[Iconize@temp,"Intermediate Lagrangian result"]];

	Return[ temp ]
]


(* ::Subsection:: *)
(*New simplification *)


(* ::Subsubsection::Closed:: *)
(*PerformSystematicFieldRedefs*)


(* ::Text:: *)
(*Performs all field redefinitions to remove all EOM terms from the Lagrangian *)


PerformSystematicFieldRedefs[lag_, maxOrder_, verbose_? BooleanQ] := Module[{devs, gaugeNormalizations, kinMix, ord, out},
	OptionalMonitor[verbose,
		out= RenormalizeMatterFields@ lag;
		(*The gauge field normalization is compensated for in future shifts of the gauge fields*)
		{gaugeNormalizations, kinMix}= GaugeFieldNormalization[out, maxOrder];
	, "Renormalizing fields \[Ellipsis]"];
	OptionalMonitor[verbose, 
		Do[
			out= ShiftLagrangian[out, gaugeNormalizations, kinMix, devs, ord, maxOrder];
		, {ord, 5, maxOrder}, {devs, Reverse@ Range[ord- 2]}];
	, StringForm["Redefining fields at dimension `1` with `2` derivatives \[Ellipsis]", ord, devs]];
	out
]


(* ::Subsubsection::Closed:: *)
(*Kinetic normalization*)


$emptyKinMix= <|FieldMap-> <||>, Fields-> {}, Zinv-> {{}}|>;


(* ::Text:: *)
(*Renormalizes the kinetic terms for the matter fields. *)


RenormalizeMatterFields@ lag_:= Module[{fields, jacobianShift, out, replacementRules, terms},
	terms= SelectOperatorDevsAndDim[lag, _? Positive, 4];
	(* all fields appearing in the input Lagrangian *)
	fields= DeleteDuplicates@ Cases[terms, Field[lab_, __]-> lab, Infinity];
	(* the free Lagrangian of all of these fields - so we get non-standard kinetic terms *)
	terms= terms- InternalSimplify[FreeLag@@ fields]// BetterExpand;

	(*Determine shifts*)
	fields= DeleteDuplicates@ Cases[terms, EoM@ Alternatives[
		Field[lab_, __], _@ Field[lab_, __], _@ _@ Field[lab_, __]]:> lab, All];
	If[Length@ fields === 0, Return@ lag];

	{replacementRules, jacobianShift}= DetermineShifts[terms, fields, <||>, $emptyKinMix];
	(*TODO: implement check that the shifts are equal to \[Delta] at leading order (Z= 1 + O(hbar)) *)

	(* insert the field expansion, return to NormalForm since we're inserting fields into EoM objects as well *)
	out= BetterExpand@ OperatorToNormalForm @ lag/. Field[f_/; MemberQ[fields, f], args__]:>
			Field[f, args]+ hbar Field[{f, 0, 1}, args];

	(* we need to prepare the Lagrangian: powers need to be removed but we cannot use operators, because we need EoM-type objects to stay unchanged *)
	out= (RelabelIndices[out, Unique-> True])/.
		Power[f_/; (!FreeQ[f, Field|FieldStrength]), k_Integer? Positive]:>
		Inactive[Times]@@ ConstantArray[f, k];

	(* plug in the rule, reactivate the powers and IBPSimplify *)
	(* Might be better to simplify elsewhere *)
	InternalSimplify[Activate[out/. replacementRules]+ jacobianShift]
	(*Collect[Operator@ Activate[out/. replacementRules], _Operator]*)
]


(* ::Text:: *)
(*Determines the inverse of the the kinetic normalization for the gauge groups. *)
(*gaugeNormalization gives the Z^-1 for the fields w/o kinetic mixing*)
(*kinMix gives the Z^-1 matrix for the Abelian fields if there is more than one (potential mixing) *)


GaugeFieldNormalization::mixnormalization= "The kinetic mixing is not canonically normalized a LO: `2`";


GaugeFieldNormalization[lag_, ord_]:= Module[
		{abelianFields, fields, gaugeNormalizations, kinMix, mu1, mu2, terms, z, z1},
	terms= SelectOperatorDevsAndDim[lag/. _Field-> 0, _? Positive, 4];
	(* all fields appearing in the input Lagrangian *)
	fields= Intersection[DeleteDuplicates@ Cases[terms, FieldStrength[lab_, __]-> lab, Infinity],
		List@@ Query[All, Key@ Field]@ $GaugeGroups];
	If[Length@ fields === 0, Return@ {<||>, $emptyKinMix}; ];
	
	abelianFields= Intersection[fields, 
		List@@ Query[Key/@ GetGaugeGroupByProperty[Group-> U1], Key@ Field]@ $GaugeGroups];
	
	(*Check if there might be kinetic mixing*)
	kinMix= If[Length@ abelianFields > 1, 
			fields= Complement[fields, abelianFields];
			z= Table[
					If[f1 =!= f2, 1 / 2, 1]* GetOperatorCoefficient[terms, 
						-1/4 FS[f1, mu1, mu2]FS[f2, mu1, mu2]]
				, {f1, abelianFields}, {f2, abelianFields}];
			{z, z1}= Transpose[CoefficientList[z, hbar, 2], {2, 3, 1}];
			If[!DiagonalMatrixQ@ z,
				Message[GaugeFieldNormalization::mixnormalization, z];
				Abort[];
			];
			(*Normalization of hbar is taking Z^{-1} at one-loop order*)
			z= DiagonalMatrix@ Power[Diagonal@ z, -1];
			z= z - hbar z . z1 . z;
			<|
				FieldMap-> Association@@ MapIndexed[#1-> First@ #2&, abelianFields],
				Fields-> abelianFields,
				Zinv-> z
			|>
		,
			$emptyKinMix
		];
	
	(*Non-mixing fields*)
	(*Read-off the normalization of the kinetic terms*)
	gaugeNormalizations= Association@@ Table[
			(*Remove canonical gauge noramlization from gauge field*)
			z= GetOperatorCoefficient[terms, FreeLag@ f/. _Coupling-> 1];
			{z, z1}= CoefficientList[z, hbar, 2];
			(*Normalization of hbar is taking Z^{-1} at one-loop order*)
			z= z^(-1) - hbar z1 z^(-2);
			f-> z
		, {f, fields}];

	{gaugeNormalizations, kinMix}
]


(* ::Subsubsection::Closed:: *)
(*ShiftLagrangian*)


(* ::Text:: *)
(*Shift the fields in the Lagrangian to remove EOM terms at a given order and number of derivatives*)


ShiftLagrangian[lag_, gaugeNormalization_, kinMix_, devs_Integer, dim_Integer, maxOrder_Integer]:= Module[
		{eomTerms, fields, jacobianShift, replacementRules, shiftOrd, termsToShift, termsToKeep},
	eomTerms= (*InternalSimplify@ *)SelectOperatorDevsAndDim[lag, devs, dim];
	shiftOrd= dim- 4;

	(*Determine shifts*)
	fields= DeleteDuplicates@ Cases[eomTerms, EoM@ Alternatives[
		Field[lab_, __], _@ Field[lab_, __], _@ _@ Field[lab_, __]]:> lab, All];
	If[Length@ fields === 0, Return@ lag];

	{replacementRules, jacobianShift}= DetermineShifts[eomTerms, fields, gaugeNormalization, kinMix];

	(* Insert the field expansion, return to NormalForm since we're inserting fields into EoM objects as well *)
	(* Shift all kinetic-mixed fields if at least one appears with an EOM *)
	If[ContainsAny[fields, kinMix@ Fields], fields= Union[fields, kinMix@ Fields]]; 
	{termsToShift, termsToKeep}= SplitLagByDims[lag, maxOrder- shiftOrd];
	termsToShift= DummyShiftTerms[termsToShift, fields, shiftOrd, maxOrder];

	(* Prepare the Lagrangian: powers need to be removed but we cannot use operators, because we need EoM-type objects to stay unchanged *)
	termsToShift= (*(RelabelIndices[termsToShift, Unique-> True])*)termsToShift/.
		Power[f_/; (!FreeQ[f, Field|FieldStrength]), k_Integer?Positive]:>
		Inactive[Times]@@ ConstantArray[f, k];

	(* Plug in the rule, reactivate the powers and IBPSimplify *)
	termsToShift= InternalSimplify[Activate[termsToShift/. replacementRules]+ jacobianShift];

	termsToShift + termsToKeep//BetterExpand
]


(* ::Text:: *)
(*Determine the field shifts needed to cancel all occurrences of EoM in the Lagrangian terms*)


DetermineShifts::notimpl = "Shift has not been implemented for type `1`";


DetermineShifts[lagTerms_, fields_, gaugeNormalization_, kinMix_]:= Module[
		{jacobianShift= 0, real, replacementRules, type, mixedFields, unmixedFields, temp},
	
	mixedFields= Intersection[fields, kinMix@ Fields];	
	unmixedFields= Complement[fields, mixedFields];

	replacementRules= Join@@ Table[
			{real, type}= Lookup[GetFields@ f, {SelfConjugate, Type}];
			Switch[type
			,Scalar, ScalarShift[lagTerms, f, real]
			,Fermion, 
				(*format: {field replacement rule, jacobian shift (from chiral fermions)}*)
				temp= FermionShift[lagTerms, f, real];
				jacobianShift+= Last@ temp;
				First@ temp
			,Vector,
				VectorShift[lagTerms, f, gaugeNormalization@ f, real]
			]
		, {f, unmixedFields}];
	
	If[Length@ mixedFields > 0, 
		replacementRules= Join[replacementRules, KinMixingShift[lagTerms, mixedFields, kinMix]];
	];
	
	{replacementRules, jacobianShift}
]


(* ::Text:: *)
(*Perform a dummy shift of all the fields that are redefined *)


DummyShiftTerms[lagTerms_, fields_, shiftOrd_, maxOrd_]:= Module[{gaugeFields, matterFields, out},
	gaugeFields= Intersection[fields, List@@Query[All, Key@ Field]@ $GaugeGroups];
	matterFields= Complement[fields, gaugeFields];

	(* Insert the field expansion, return to NormalForm since we're inserting fields into EoM objects as well *)
	out= OperatorToNormalForm @lagTerms;

	Do[
		(*All gauge fields are from the shift and should eventually be substituted*)
		out= DummyGaugeShift[out, f]/. Field[f, args__]:> Field[{f, shiftOrd, 0}, args]+ hbar Field[{f, shiftOrd, 1}, args];
		out= TruncateShift[out, f, maxOrd];
	, {f, gaugeFields}];

	out= SeriesEFT[out/. Field[f_/; MemberQ[matterFields, f], args__]:>
			Field[f, args]+ Field[{f, shiftOrd, 0}, args]+ hbar Field[{f, shiftOrd, 1}, args],
		EFTOrder-> maxOrd];

	out
]


TruncateShift[lag_, fieldLab_, maxOrd_]:=
	Plus@@ (If[FreeQ[#, Field[{fieldLab, __}, __]], #, TruncateOperator[#, maxOrd]]&)/@ TermsToList@ lag;


(* ::Text:: *)
(*Function to ensure that terms with EOMs of multiple distinct fields  get shifted the right amount *)


AdjustEOMShifts[shift_, field_Symbol]:= Module[{terms},
	terms= TermsToList@ shift;
	(*Include only the shift for the field of the leading canonical (lexiographic) order*)
	terms= If[FreeQ[Operator@ #, EoM@ Alternatives[Field[lab_, __], _@ Field[lab_, __], _@ _@ Field[lab_, __]]/;
		Order[field, lab] === -1], #, 0]&/@ terms;
	(*terms= Power[Length@ DeleteDuplicates@ Cases[Operator@ #,
		EoM@ Alternatives[Field[lab_, __], _@ Field[lab_, __], _@ _@ Field[lab_, __]]/; lab=!= field:> lab,
		All]+ 1, -1]* #&/@ terms;*)
	Plus@@ terms
]


(* ::Subsubsection::Closed:: *)
(*Scalar shift*)


ScalarShift[lagTerms_, f_Symbol, real_? BooleanQ]:= Module[{chi1, chi2, devPat, fieldPattern, fieldShift, inds, pInds, rules, terms},
	inds= Symbol["i" <> ToString@ #]&/@ Range@ Length@ GetFields[f, Indices];
	pInds= Pattern[#, _]&/@inds;

	(*Keep only relevant terms*)
	terms= Plus@@ Select[TermsToList@ lagTerms, Not@* FreeQ[f]];

	fieldShift= If[real,
			RelabelIndices[OperatorToNormalForm[CoefficientOperator[terms,
				EoM@ f[Sequence@@ inds]], Unique-> True], Unique->True]
		,
			(* factor out the EoM[f] term *)
			chi1= RelabelIndices[OperatorToNormalForm[CoefficientOperator[terms,
				EoM@ f[Sequence@@ inds]], Unique-> True], Unique-> True];

			(* ...and subtract the terms we found, to now get the EoM[Bar@f] coefficients *)
			chi2= RelabelIndices[CoefficientOperator[
				RelabelIndices[terms - Operator[chi1 EoM@f[Sequence@@ inds] ] ],
				EoM@ Bar@ f[Sequence@@ inds] ], Unique->True];
			1/2 Bar@ OperatorToNormalForm[(chi1 + Bar@chi2), Unique-> True]
		];
	fieldShift= AdjustEOMShifts[fieldShift, f];

	rules= Table[
		fieldPattern= First@ Cases[{f[Sequence@@ pInds]},
			Field[lab_, type_, is_, {}]:> Field[{lab, _, n}, type, is, devPat_], Infinity];
		With[{shift= Coefficient[fieldShift, hbar, n], tDevs= devPat},
			fieldPattern:> RelabelIndices[CD[tDevs, shift], Unique-> True]
		]
	, {n, 0, 1}]
]


(* ::Subsubsection::Closed:: *)
(*Fermion shift*)


FermionShift[lagTerms_, f_Symbol, real_? BooleanQ]:= Module[
		{chi1, chi2, devPat, eomField, fieldPattern, fieldShift, inds, jacobianShift, pInds, rules, terms},
	inds= Symbol["i" <> ToString@ #]&/@ Range@ Length@ GetFields[f, Indices];
	pInds= Pattern[#, _]&/@inds;

	(*Keep only relevant terms*)
	terms= Plus@@ Select[TermsToList@ lagTerms, Not@* FreeQ[f]];

	eomField= EoM@ First@ Cases[{f[Sequence@@inds]}, _Field, Infinity];

	fieldShift= If[real,
			(* factor out the EoM[f] term *)
			chi1=RelabelIndices[OperatorToNormalForm[CoefficientOperator[terms, eomField],
				Unique-> True], Unique-> True];
			(* ...and subtract the terms we found, to now get the EoM[Bar@f] coefficients *)
			chi2=RelabelIndices[OperatorToNormalForm[CoefficientOperator[RelabelIndices[
				terms-Operator[chi1** eomField]], Transp/@ eomField], Unique->True], Unique->True];
			RelabelIndices[I*CC**(Transp@chi1-chi2)]
		,
			(* factor out the EoM[f] term *)
			chi1= RelabelIndices[OperatorToNormalForm[CoefficientOperator[terms, eomField],
				Unique-> True], Unique-> True];
			(* ...and subtract the terms we found, to now get the EoM[Bar@f] coefficients *)
			chi2=RelabelIndices[CoefficientOperator[RelabelIndices[
				terms- Operator[chi1** eomField]], Bar/@ eomField], Unique-> True];
			-I/2* Bar@ OperatorToNormalForm[chi1+ Bar@ chi2, Unique-> True]
		];
	fieldShift= AdjustEOMShifts[fieldShift, f];

	(*TODO: include chiral Jacobian*)
	jacobianShift= 0;

	rules= Table[
		fieldPattern= First@ Cases[{f[Sequence@@ pInds]},
			Field[lab_, type_, is_, {}]:> Field[{lab, _, n}, type, is, devPat_], Infinity];
		With[{shift= Coefficient[fieldShift, hbar, n], tDevs= devPat},
			fieldPattern:> RelabelIndices[CD[tDevs, shift], Unique-> True]
		]
	, {n, 0, 1}];

	{rules, jacobianShift}
]


(* ::Subsubsection::Closed:: *)
(*Vector shift*)


VectorShift::cmplx= "Complex vectors are not supported currently."


(* ::Text:: *)
(*The normalization is there to compensate for the normalization of the gauge kinetic term. *)


VectorShift[lagTerms_, f_Symbol, normalization_, real_? BooleanQ]:= Module[{devPat, fieldPattern, fieldShift, inds, pInds, rules, terms},
	inds= Symbol["i" <> ToString@ #]&/@ Range[Length@ GetFields[f, Indices]+ 1];
	pInds= Pattern[#, _]&/@inds;

	terms= Plus@@ Select[TermsToList@ lagTerms, Not@* FreeQ[f]];

	fieldShift= If[real,
			-normalization* RelabelIndices[OperatorToNormalForm[
				CoefficientOperator[terms, EoM@ f[Sequence@@inds] ],Unique->True, CanonizeKinetic->True]
				, Unique-> True]
		,
			Message[VectorShift::cmplx]; Abort[];
		];
	fieldShift= AdjustEOMShifts[fieldShift, f];

	rules= Table[
		fieldPattern= First@ Cases[{f[Sequence@@ pInds]},
			Field[lab_, type_, is_, {}]:> Field[{lab, _, n}, type, is, devPat_], Infinity];
		With[{shift= Coefficient[fieldShift, hbar, n], tDevs= devPat},
			fieldPattern:> RelabelIndices[CD[tDevs, shift], Unique-> True]
		]
	, {n, 0, 1}]
]


(* ::Text:: *)
(*Kinetic mixing assumes (currently) mixing between real U(1) fields *)


KinMixingShift[lagTerms_, eomFields_List, kinMix_Association]:= Module[
		{devPat, fieldPattern, fieldShifts, inds, pInds, rules, fShift, terms},
	inds= Symbol["i" <> ToString@ #]&/@ Range[1];
	pInds= Pattern[#, _]&/@inds;

	terms= Plus@@ Select[TermsToList@ lagTerms, Not@* FreeQ[Alternatives@@ eomFields]];
	fieldShifts= Sum[
			fShift= RelabelIndices[OperatorToNormalForm[
						CoefficientOperator[terms, EoM@ f[Sequence@@ inds] ],Unique->True, CanonizeKinetic->True]
						, Unique-> True];
			fShift= AdjustEOMShifts[fShift, f];
			fShift* UnitVector[Length@ kinMix@ Fields, kinMix[FieldMap, f] ]
		, {f, eomFields}];
	fieldShifts= -kinMix@ Zinv . fieldShifts;
	
	rules= Flatten@ Table[
		fieldPattern= First@ Cases[{f[Sequence@@ pInds]},
			Field[lab_, type_, is_, {}]:> Field[{lab, _, n}, type, is, devPat_], Infinity];
		With[{shift= Coefficient[fieldShifts[[kinMix[FieldMap, f]]], hbar, n], tDevs= devPat},
			fieldPattern:> RelabelIndices[CD[tDevs, shift], Unique-> True]
		]
	, {f, kinMix@ Fields}, {n, 0, 1}]
]


(* ::Subsubsection::Closed:: *)
(*Utility*)


(* ::Text:: *)
(*Selects all operators with the required number of derivatives and operator dimension (in internal Operator form)*)


SelectOperatorDevsAndDim[lag_, devs_, dim_]:= Module[{gaugeFields},
	(*lag/. op_Operator/; !MatchQ[OpDevsAndDim@op, {devs, dim}]-> 0 *) (*coefficient is part of the counting*)
	(*Plus@@ Select[TermsToList@ lag, MatchQ[{devs, dim}]@* OpDevsAndDim]*)
	gaugeFields= List@@ Query[All, Key@ Field]@ $GaugeGroups;
	Plus@@ Select[Select[TermsToList@lag, MatchQ[dim]@* OperatorDimension],
		MatchQ[devs]@ OpDevs[#, gaugeFields]&]
];


OpDevs[op_, gaugeFields_]:= Block[{devs= 0},
	devs+= Plus@@ Cases[op, EoM@
		Alternatives[Field[_, type_, __], _@ Field[_, type_, __], _@ _@ Field[_, type_, __]]:>
			EOMDevs@ type, All];
	devs+= Plus@@ Cases[op, FieldStrength[lab_, __, linds_]:>
		If[MemberQ[gaugeFields, lab], 2, 1]+ Length@ linds, All];
	devs+= Plus@@ Cases[op, Field[___, linds_]:> Length@ linds, All];
	devs
];


(*OpDevsAndDim@ op_:= Block[{devs= 0, gaugeFields, dim},
	gaugeFields= List@@ Query[All, Key@ Field]@ $GaugeGroups;
	dim= OperatorDimension@ op;

	devs+= Plus@@ Cases[op, EoM@
		Alternatives[Field[_, type_, __], _@ Field[_, type_, __], _@ _@ Field[_, type_, __]]:>
			EOMDevs@ type, All];
	devs+= Plus@@ Cases[op, FieldStrength[lab_, __, linds_]:>
		If[MemberQ[gaugeFields, lab], 2, 1]+ Length@ linds, All];
	devs+= Plus@@ Cases[op, Field[___, linds_]:> Length@ linds, All];
	{devs, dim}
];*)


(* ::Text:: *)
(*Separates the Lagrangian into two parts depending on the dimension *)


SplitLagByDims[lag_, dim_]:= Block[{},
	Plus@@ (If[OperatorDimension@ # > dim, {0, #}, {#, 0}]&)/@ TermsToList@ lag
];


(* ::Subsubsection::Closed:: *)
(*Gauge group shifts*)


DummyGaugeShift[expr_, gaugeField_Symbol]:= Module[{group, out},
	group= First@ GetGaugeGroupByProperty[Field-> gaugeField];

	(*If group is Abelian include its own field strength even though it doesn't transform*)
	out= If[$GaugeGroups[group, Group] === U1,
			expr/. f:FieldStrength[gaugeField, __]:> f+ ShiftAbelianFS@ f
		,
			expr
		];

	(*Shift of CDs*)
	out= out/. (f:(Field[__, {__}]| _FieldStrength))/; FieldTransformsUnderGaugeGroupQ[f, group]:>
		ShiftCDnew[f, group, gaugeField];

	out
];


ShiftAbelianFS[FieldStrength[fsLab_, {l1_, l2_}, {}, devs_]]:= Block[{},
	Field[fsLab, Vector@ l2, {}, Append[devs, l1]]- Field[fsLab, Vector@ l1, {}, Append[devs, l2]]
];


(* ::Text:: *)
(*Applies iteratively*)


ShiftCDnew[fieldOrFS_[fieldLab_, typeOrLors_, inds_, devs:{__}], group_, gaugeField_]:= Module[
		{charge, next, newAdj, newFieldInds, newInd, origInd},
	If[$GaugeGroups[group, Abelian],
		charge= FirstCase[GetFields[fieldLab, Charges], group@ q_-> q];
		next= ShiftCDnew[fieldOrFS[fieldLab, typeOrLors, inds, Rest@ devs], group, gaugeField];
		-I* charge* Field[gaugeField, Vector@ First@ devs, {}, {}]* next + CD[First@ devs, next]
	,
		{newFieldInds, origInd, newInd}= ChangeGroupIndex[inds, group, Unique[a]];
		newAdj= Index[Unique@ A, group@ adj];
		CD[First@ devs, ShiftCDnew[fieldOrFS[fieldLab, typeOrLors, inds, Rest@ devs], group, gaugeField]] +
			Field[gaugeField, Vector@ First@ devs, {newAdj}, {}] If[Head@ origInd === Bar,
					+I* Bar@ CG[gen@ origInd[[1, 2]], {newAdj, Bar@ origInd, newInd}],
					-I* CG[gen@ origInd[[2]], {newAdj, origInd, Bar@ newInd}]
				]* ShiftCDnew[fieldOrFS[fieldLab, typeOrLors, newFieldInds, Rest@ devs], group, gaugeField]
	]
];
ShiftCDnew[f:Field[_, _, _, {}], _, _]:= f;
ShiftCDnew[f:FieldStrength[_, _, _, {}], group_, gaugeField_]:= f+ ShiftCoreFS[f, group, gaugeField];


ShiftCoreFS[FieldStrength[fsLab_, {l1_, l2_}, grInds_, {}], group_, gaugeField_]:= Module[
		{charge, newAdj, newFsInds, newInd, origInd},
	If[fsLab === $GaugeGroups[group, Field], Return[
		(* Only non-Abelian gauge groups get to this point *)
		Field[fsLab, Vector@ l2, grInds, {l1}]-
		Field[fsLab, Vector@ l1, grInds, {l2}] +
		Module[{b, c}, CG[fStruct@ group, {grInds[[1, 1]], b, c}]*
			fsLab[l1, b]fsLab[l2, c] ]
	]; ];
	If[$GaugeGroups[group, Abelian],
		charge= FirstCase[GetFields[fsLab, Charges], group@ q_-> q];
		-I* charge* CD[devs, Field[gaugeField, Vector@ l1, {}, {}] Field[fsLab, Vector@ l2, grInds, {}]-
			Field[gaugeField, Vector@ l2, {}, {}] Field[fsLab, Vector@ l1, grInds, {}]]
	,
		{newFsInds, origInd, newInd}= ChangeGroupIndex[grInds, group, Unique[a]];
		newAdj= Index[Unique@ A, group@ adj];
		If[Head@ origInd === Bar,
			+I* Bar@ CG[gen@ origInd[[1, 2]], {newAdj, Bar@ origInd, newInd}],
			-I* CG[gen@ origInd[[2]], {newAdj, origInd, Bar@ newInd}] ]*
			CD[devs,
				Field[gaugeField, Vector@ l1, {newAdj}, {}]* Field[fsLab, Vector@ l2, newFsInds, {}]-
				Field[gaugeField, Vector@ l2, {newAdj}, {}]* Field[fsLab, Vector@ l1, newFsInds, {}]
			]
	]
];


(* Function finding the index in a list that belongs to a specific group, returning the list with the index replaced by a given label, the old index, and the new *)
ChangeGroupIndex[indices_, group_, newLab_]:=Module[{groupIndex, newInd},
	(* Find the index belongin to the given group *)
	groupIndex= FirstCase[indices, (Index[_, rep_]| Bar@ Index[_, rep_])/; GroupFromRep@ rep === group, {}];

	(* Return the field with replaced index and the original and new index labels *)
	newInd= If[Head@ groupIndex === Bar,
			Bar@ Index[newLab, groupIndex[[1, 2]]]
		,
			Index[newLab, Last[groupIndex]]
		];
	{indices/. groupIndex-> newInd, groupIndex, newInd}
];


(* ::Subsection:: *)
(*Handling off-diagonal kinetic terms*)


(* ::Subsubsection::Closed:: *)
(*Modules used in all cases*)


Options[BiPerturbativeInvert]={Triangular -> True};

BiPerturbativeInvert[M_,{\[Lambda]1_,l1_,n1_},{\[Lambda]2_,l2_,n2_},OptionsPattern[]]:=Module[{val},
	val = OptionValue[Triangular];
	If[val===True,
		BiPerturbativeInvertTriangular[M,{\[Lambda]1,l1,n1},{\[Lambda]2,l2,n2}]
	,
		BiPerturbativeInvertRegular[M,{\[Lambda]1,l1,n1},{\[Lambda]2,l2,n2}]
	]
]

BiPerturbativeInvertTriangular[M_,{\[Lambda]1_,l1_,n1_},{\[Lambda]2_,l2_,n2_}]:=Module[{R,r,m,eqs,var},
m=Length[M];
R= (IdentityMatrix[m]+UpperTriangularize@Table[Sum[(1-KroneckerDelta[k1,0]KroneckerDelta[k2,0])\[Lambda]1^k1 \[Lambda]2^k2 r[k1,k2,i,j],{k1,l1,n1},{k2,l2,n2}],{i,m},{j,m}]);
eqs=DeleteCases[Thread[0==Flatten[CoefficientList[Normal@Series[ComplexExpand[R\[ConjugateTranspose], r[_,_,_,_],TargetFunctions->Conjugate] . M . R-IdentityMatrix[m],{\[Lambda]1,0,n1},{\[Lambda]2,0,n2}],{\[Lambda]1,\[Lambda]2}]]],True]/.r[a_,b_,i_,j_]\[Conjugate]:> r[a,b,j,i];

var=Cases[eqs,r[__],{0,Infinity}];

Quiet[R/.(Last@Solve[Reduce@eqs,var]),{Solve::incnst}]
]

BiPerturbativeInvertRegular[M_,{\[Lambda]1_,l1_,n1_},{\[Lambda]2_,l2_,n2_}]:=Module[{R,r,m,eqs,var},

m=Length[M];
R= (IdentityMatrix[m]+Table[Sum[(1-KroneckerDelta[k1,0]KroneckerDelta[k2,0])\[Lambda]1^k1 \[Lambda]2^k2 r[k1,k2,i,j],{k1,l1,n1},{k2,l2,n2}],{i,m},{j,m}]);
eqs=DeleteCases[Thread[0==Flatten[CoefficientList[Normal@Series[R . M -IdentityMatrix[m],{\[Lambda]1,0,n1},{\[Lambda]2,0,n2}],{\[Lambda]1,\[Lambda]2}]]],True]/.r[a_,b_,i_,j_]\[Conjugate]:> r[a,b,j,i];

var=Cases[eqs,r[__],{0,Infinity}];

Quiet[R/.(Last@Solve[Reduce@eqs,var]),{Solve::incnst}]
]
