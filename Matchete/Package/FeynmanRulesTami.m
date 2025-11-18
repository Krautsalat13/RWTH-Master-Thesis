(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Feynman Rules Finder`*)


(* ::Subtitle:: *)
(*Feynman Rule finder for for Matchete packet.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


PackageImport["Matchete`"]


(* ::Subsection:: *)
(*Exported*)


PackageExport["PrepareLagrangian"]
PackageExport["DerivativesToMomenta"]
PackageExport["FeynmanRule"]
PackageExport["EXT"]
PackageExport["ExpandFieldStrengthTensors"]
PackageExport["ExpandLagrangian"]
PackageExport["FeynmanRuleFinderAll"]
PackageExport["FeynmanRules"]

(*for debugging*)
PackageExport["splitFS"]
PackageExport["position"]
PackageExport["FindCountofFields"]
PackageExport["Feynman"]
PackageExport["AssignExternalFields"]
PackageExport["Feynman"]
PackageExport["ScalarOneTerm"]
PackageExport["VectorOneTerm"]
PackageExport["FermionOneTerm"]
PackageExport["SubstituteFieldIndices"]
PackageExport["CountFieldOccurrences"]
PackageExport["GatherFieldsByType"]
PackageExport["ComputeGrassmannSign"]
PackageExport["ExpandCovariantDerivatives"]
PackageExport["ReplaceCouplingConvention"]
PackageExport["RenameDummyIndices"]
PackageExport["ReplaceFieldIndices"]


(* ::Subsection:: *)
(*Package Scope*)


PackageScope["RemoveMomentumLabels"]
PackageScope["ExtractExternalLabel"]
PackageScope["SplitSingleDerivative"]
PackageScope["CovariantDerivativeSplit"]
PackageScope["charges"]
PackageScope["NonCharges"]
PackageScope["splitFS"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection:: *)
(*Exported*)


FeynmanRules::usage  = "FeynmanRule[L_,input_], calculates Feynman Rules for Lagrangian L with optional Correlator input_";
ExpandLagrangian::usage = "ExpandLagrangian[L] fully expands sums and products in L while preserving non-commutative field order when needed.";
ExpandCovariantDerivatives::usage = "ExpandCovariantDerivatives[term] expands covariant derivatives acting on fields, splitting gauge and non-gauge contributions.";
DerivativesToMomenta::usage = "DerivativesToMomenta[L] converts partial-derivative operators to momentum insertions consistent with the chosen sign conventions.";
AssignExternalFields::usage = "AssignExternalFields[L, ext] labels external legs with unique symbols for later momentum/Grassmann bookkeeping; ext may be Automatic.";
ComputeGrassmannSign::usage = "ComputeGrassmannSign[order, term] computes the sign needed to reorder Grassmann-odd fields to the canonical order.";
GatherFieldsByType::usage = "GatherFieldsByType[L, type] filters interactions in L that contain fields of the given Type (Scalar|Fermion|Vector|Ghost).";
ExpandFieldStrengthTensors::usage = "ExpandFieldStrengthTensors[L] expands F_{\\[Mu]\\[Nu]} terms into their derivative and commutator pieces using the model's gauge structure.";
Momenta::usage = "Momenta[solution] extracts momentum replacement rules generated during derivative-to-momentum mapping for a single vertex.";
ReplaceCouplingConvention::usage = "ReplaceCouplingConvention[expr] rewrites couplings and field-strength normalizations to the library's internal conventions.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Calculate combinatorial factors*)


(* ::Subsection:: *)
(*Functions that are unify all three cases*)


Feynman[L_,pD_]:= Module[{ExpandLagrangianAndDerivatives,listofoperators,feynmanrulesD,sumornot},
ExpandLagrangianAndDerivatives = Expand[ExpandCovariantDerivatives[Expand[L]]];
sumornot[A_]:=If[Head[A]===Plus,List@@A,{A}];
feynmanrulesD = I*Total[(VectorOneTerm /@ sumornot@Total@(FermionOneTerm /@ sumornot @ Total @ (ScalarOneTerm/@ sumornot[ExpandLagrangianAndDerivatives])))];
If[pD == D, Return[feynmanrulesD]];
DerivativesToMomenta[feynmanrulesD]
]


(* ::Subsection:: *)
(*Scalar*)


(* ::Subsubsection:: *)
(*ScalarFeynman*)


ScalarFeynman[L_]:=
 Module[{listofoperators},

listofoperators = List @@ L ;
ScalarOneTerm /@listofoperators
]


(* ::Subsubsection:: *)
(*ScalarOneTerm*)


ScalarOneTerm[Lterm_]:= Module[{listoffieldsbarred,listoffieldsjoined,listoffields, listfieldn, emptyL,Lint, nlist,fieldname},

Lint ={Lterm/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}]}; 

listoffieldsjoined = GatherFieldsByType[Lint,Scalar];
emptyL = Lint/.Bar[Field[f_,Scalar,type2_,indices_]]:>Bar[Field[f,Scalar,{},{}]]/.Field[f_,Scalar,type2_,indices_]:>Field[f,Scalar,{},{}]/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];

nlist = CountFieldOccurrences[#,emptyL] &/@ listoffieldsjoined;
listfieldn = Transpose[{listoffieldsjoined,nlist}];

Total[{Activate[Fold[factorScalar,Lint, listfieldn]]}, Infinity]
]


(* ::Subsubsection:: *)
(*FactorScalar*)


factorScalarDirect[L_, {field_,n_}]:=factorScalar[ L/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}],{field,n}];
factorScalar[L_, {field_,n_}]:= Module[{list,sign,newList,signList,LastList,external,counter, EXT},
If[n== 0,Return[L]];

list = Permutations[Range[n]];

LastList= Table[SubstituteFieldIndices[field, Flatten[L][[i]], #]&/@ list,{i,Length[Flatten[L]]}];
Flatten[LastList]

]


(* ::Subsection:: *)
(*Fermion*)


(* ::Subsubsection:: *)
(*FermionFeynman*)


FermionFeynman[L_]:= Module[{listofoperators},

listofoperators = List @@ L ;
FermionOneTerm /@listofoperators
]


(* ::Subsubsection:: *)
(*FermionOneTerm*)


FermionOneTerm[Lterm_]:= Module[{listoffieldsjoined,listfieldn, emptyL,Lint, nlist,fieldname,vzg,Lintold,reswithnovz},

Lintold ={Lterm/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}]}; 
Lint = Lintold/.Times:>Inactive[Times];
listoffieldsjoined = GatherFieldsByType[Lint,Fermion];
emptyL = Lint/.Bar[Field[f_,Fermion,type2_,indices_]]:>Bar[Field[f,Fermion,{},{}]]/.Field[f_,Fermion,type2_,indices_]:>Field[f,Fermion,{},{}]/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];


nlist = CountFieldOccurrences[#,emptyL] &/@ listoffieldsjoined;
listfieldn = Transpose[{listoffieldsjoined,nlist}];
vzg = newGrassmanncomparision[Lint,listfieldn];
reswithnovz =Total[Activate[{Fold[factorFermion,Lint, listfieldn]}], Infinity];
vzg *reswithnovz
]


(* ::Subsubsection:: *)
(*FactorFermion*)


factorFermionDirect[L_, {field_,n_}]:=factorFermion[ L/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}],{field,n}];


factorFermion[L_, {field_,n_}]:= Module[{list,LastList,levilist,vz},
If[n== 0,Return[L]];
list = Permutations[Range[n]];

LastList= Table[SubstituteFieldIndices[field, Flatten[L][[i]], #]&/@ list,{i,Length[Flatten[L]]}];
vz = Extract[LeviCivitaTensor[n],#]&/@ list;

levilist = Map[#*vz&,LastList];
levilist
]


(* ::Subsubsection:: *)
(*Incorporate change difference due to Grassmann type variables*)


GrassMann[Lint_]:=Module[{pos,bar},
pos = Position[Lint,Field[a_?(Head[#]=!=EXT&),Fermion|Ghost,___]];
bar = Position[Head[Extract[Lint,Most[#]]]&/@ Position[Lint,Field[a_?(Head[#]=!=EXT&),Fermion|Ghost,___]],Bar];
Do[pos[[i]]=Most[pos[[i]]],{i,Flatten[bar]}];
pos
]
newGrassmanncomparision[interactions_,listfieldn_]:=Module[{vorzeichen,onelesstoworry},
vorzeichen = {};
onelesstoworry[intold_,fieldlist_]:= Module[{field,int,grassmannsinL,positions,rule,Vorz,intnew},
field = fieldlist[[1]];
int = intold/. Field[x_,type_,_,_]:> Field[x,type,{},{}]/;!MatchQ[x,EXT[_]] /.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];
grassmannsinL= GrassMann[int];
positions = Position[int,field];
rule=Field[x_,type_,a_,b_]:>Field[EXT[x],type,a,b];
Vorz =Times@@( Power[-1,#+1] &/@(Flatten[Position[grassmannsinL,#]&/@ positions]));
AppendTo[vorzeichen,Vorz];
intnew =ReplacePart[int,#:> (Extract[int, #]/.rule)&/@positions];
intnew
];
Fold[onelesstoworry,interactions,listfieldn];
Times@@vorzeichen
]


(* ::Subsection:: *)
(*Vector*)


(* ::Subsubsection:: *)
(*VectorFeynman*)


VectorFeynman[L_]:=
 Module[{listofoperators},

listofoperators = List @@ L ;
VectorOneTerm /@listofoperators
]


(* ::Subsubsection:: *)
(*VectorOneTerm*)


VectorOneTerm[Lterm_]:= Module[{listoffieldsbarred,listoffieldsjoined,listoffields, listfieldn, emptyL,CountFieldOccurrences,Lint, nlist,fieldname},

Lint ={Lterm/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}]}; 
listoffieldsbarred=DeleteDuplicates[Cases[Lint,Bar[Field[f_,Vector[_],type2_,indices_]]:>Bar[Field[f,Vector,{},{}]],\[Infinity]]];
listoffields=DeleteDuplicates[Cases[Lint,Field[f_,Vector[_],type2_,indices_]:>Field[f,Vector,{},{}],\[Infinity]]];
listoffieldsjoined = Join[listoffieldsbarred, listoffields];

emptyL = Lint/.Bar[Field[f_,Vector[lorentz_],type2_,indices_]]:>Bar[Field[f,Vector,{},{}]]/.Field[f_,Vector[lorentz_],type2_,indices_]:>Field[f,Vector,{},{}]/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];

CountFieldOccurrences[field_]:=Module[{n},
fieldname =field/. Field[x_,___]:>x;
Which[
GetFields[fieldname][SelfConjugate], (*Is field real*)
n = Count[emptyL,field, Infinity],
MatchQ[field,Bar[_]], (*Is field Bared of complex*)
n = Count[emptyL, field, Infinity],
True, (*Field is nonbared but complex*)
n =  Count[emptyL,field, Infinity] - Count[emptyL,Bar[field], Infinity]]
];

nlist = CountFieldOccurrences/@ listoffieldsjoined;
listfieldn = Transpose[{listoffieldsjoined,nlist}];

Total[{Activate[Fold[factorVector,Lint, listfieldn]]}, Infinity]
]




(* ::Subsubsection:: *)
(*FactorVector*)


factorVectorDirect[L_, {field_,n_}]:=factorVector[ L/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}],{field,n}];
factorVector[L_, {field_,n_}]:= Module[{list,sign,newList,signList,LastList,external,SubstituteFieldIndices,counter},
If[n== 0,Return[L]];

list = Permutations[Range[n]];

SubstituteFieldIndices[f_,expression_,list_] := Module[{fieldname},
counter=1;
fieldname =f/. Field[x_,___]:>x;
Which[
GetFields[fieldname][SelfConjugate], (*Is field real*)
 expression /. Field[fieldname,Vector[lorentz___], flavor_, index_]:>Module[{result = Field[EXT[Subscript[fieldname,list[[counter]]]], Vector[lorentz],flavor,index]},counter++;result],

MatchQ[field,Bar[_]], (*Is field Bared of complex*)
expression/. Bar[Field[fieldname,Vector[lorentz___], flavor_, index_]]:>Module[{result = Bar[Field[EXT[Subscript[OverBar[fieldname],list[[counter]]]], Vector[lorentz],flavor,index]]}, counter++;result],

True, (*Field is nonbared but complex*)
expression/. Field[fieldname,Vector[lorentz___], flavor_, index_]:>Module[{result = Field[EXT[Subscript[fieldname,list[[counter]]]], Vector[lorentz],flavor,index]}, counter++;result]
]];

LastList= Table[SubstituteFieldIndices[field, Flatten[L][[i]], #]&/@ list,{i,Length[Flatten[L]]}];
Flatten[LastList]

]



(* ::Subsection:: *)
(*SubstituteFieldIndices, CountFieldOccurrences and ListofFieldsjoined*)


(* ::Subsubsection:: *)
(*SubstituteFieldIndices Function*)


SubstituteFieldIndices[f_,expression_,list_] := Module[{fieldname,type,counter},
counter=1;
fieldname =f/. Field[x_,___]:>x;
type =f/. Field[_,y_ , ___]:>y;

Which[
GetFields[fieldname][SelfConjugate], (*Is field real*)
 expression /. Field[fieldname,type, flavor_, index_]:>Module[{result = Field[EXT[Subscript[fieldname,list[[counter]]]], type ,flavor,index]},counter++;result],

MatchQ[f,Bar[_]], (*Is field Bared of complex*)
expression/. Bar[Field[fieldname,type, flavor_, index_]]:>Module[{result = Bar[Field[EXT[Subscript[OverBar[fieldname],list[[counter]]]], type,flavor,index]]}, counter++;result],

True, (*Field is nonbared but complex*)
expression/. Field[fieldname,type, flavor_, index_]:>Module[{result = Field[EXT[Subscript[fieldname,list[[counter]]]], type,flavor,index]}, counter++;result]
]];


(* ::Subsubsection:: *)
(*Count occurence of given Field*)


CountFieldOccurrences[field_,empty_]:=Module[{n,fieldname},
fieldname =field/. Field[x_,___]:>x;
Which[
GetFields[fieldname][SelfConjugate], (*Is field real*)
n = Count[empty,field, Infinity],
MatchQ[field,Bar[_]], (*Is field Bared of complex*)
n = Count[empty, field, Infinity],
True, (*Field is nonbared but complex*)
n =  Count[empty,field, Infinity] - Count[empty,Bar[field], Infinity]]
];


(* ::Subsubsection:: *)
(*Gather Fields by Type for Lagrangian*)


GatherFieldsByType[L_,typeoffield_]:=Module[{listoffieldsbarred,listoffields},
	listoffieldsbarred=DeleteDuplicates[Cases[L,Bar[Field[f_,typeoffield,type2_,indices_]]:>Bar[Field[f,typeoffield,{},{}]],\[Infinity]]];
	listoffields=DeleteDuplicates[Cases[L,Field[f_,typeoffield,type2_,indices_]:>Field[f,typeoffield,{},{}],\[Infinity]]];
	Join[listoffieldsbarred, listoffields]
];


(* ::Section:: *)
(*Expand Covariant Derivatives and Field Strength Tensors*)


(* ::Subsection:: *)
(*Derivatives into momenta*)


DerivativesToMomenta[L_]:= Module[{Lsub,Lsubbar},
	Lsub = L /. Field[f_,type_,flavor_,index___] :>  
	Field[f,type,flavor,{}]*Subscript["p",{Flatten[{RemoveMomentumLabels[index]}],ExtractExternalLabel[f]}]/. Subscript["p",{{},_}]:> 1;

	Lsubbar = Lsub/. Bar[Subscript["p",{index_List,f_}]]:>  
Module[{table},table=Table[-I*Subsuperscript["p",Bar[RemoveMomentumLabels[index[[i]]]],ExtractExternalLabel[f]],{i,Length[index]}]; Times@@table];
	Lsubbar/. Subscript["p",{index_List,f_}]:>
		Module[{table},table=Table[-I*Subsuperscript["p",index[[i]],f],{i,Length[index]}]; Times@@table]
]

RemoveMomentumLabels[expr_]:=expr/. Index[s_Symbol,Lorentz]:>Index[Symbol[StringDrop[SymbolName[s],-1]],Lorentz]/;StringEndsQ[SymbolName[s],"p"]

ExtractExternalLabel[expr_]:=expr/. EXT[x_]:>x


(* ::Subsection:: *)
(*Separate the Covariant Derivative into its pieces*)


ExpandCovariantDerivatives[termold_]:=Module[{term,n},
term = termold /.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];
n = Max[Cases[term,Field[_,_,_,indices_]:>Length[indices],{0,Infinity}]];
If[n == -\[Infinity], n = 0];
Nest[SplitSingleDerivative,term,n]]

SplitSingleDerivative[term_]:=term/. Field[f___,derivative_]:>Module[{DD},
DD=SelectFirst[derivative,!MatchQ[#,Index[s_Symbol,Lorentz]/;StringEndsQ[SymbolName[s],"p"]]&,None];
If[DD===None,Field[f,derivative],CovariantDerivativeSplit[Field[f,derivative],{DD}]]];

CovariantDerivativeSplit[field_, derivative_]:= Module[{lorentz,lorentzrest,partiallabel,derivativedirect,sign,fieldname,type,index,derivatives,grouppart,chargepart,partialdiff,lorentzderivedhere},
	If[Length[derivative]== 0, Return[field]];
	derivativedirect = derivative[[1]];
	If[MatchQ[derivativedirect,s_Symbol/;StringEndsQ[SymbolName[s],"p"]],Return[field], Nothing];
	sign =  If[!FreeQ[field,Bar],-1,1];
	{fieldname,type,index,lorentz} = Cases[field, Field[x_,type_,index_,lorentz_]:>{x,type, index,lorentz} ,{0,Infinity}][[1]];
	grouppart = NonCharges[fieldname,type,index,derivativedirect];
	chargepart = charges[fieldname,type,index,derivativedirect];
	partiallabel=derivativedirect/. Index[s_Symbol,type_]:>Index[Symbol[SymbolName[s]<>"p"],type];
	partialdiff = Field[fieldname,type,index,{partiallabel}];
	lorentzrest = DeleteCases[lorentz,derivativedirect,1,1];
	CD[lorentzrest,partialdiff - sign*grouppart - sign*chargepart]
]

charges[field_,type_,groupindex_,derivative_]:= Module[{CouplingChargeGroupCorrect,ChargeFieldswithlabels,lorentz,lorentzderivedhere,lorentzrestderivatives,charge,chargeGroup,CouplingChargeGroup,ChargeGroupField,chargedecomp},
If[GetFields[field][Charges] ==Missing["KeyAbsent",field][Charges], Return[0]];

charge = Last /@ GetFields[field][Charges];
chargeGroup = Head /@GetFields[field][Charges];

ChargeGroupField = GetGaugeGroups[#][Field]&/@ chargeGroup;
ChargeFieldswithlabels = MapThread[Field[#1,Vector[derivative],{},{}]&,{ChargeGroupField}];



chargedecomp = I*Field[field,type,groupindex,{}]*ChargeFieldswithlabels*charge;
Total[chargedecomp]
]

NonCharges[field_,type_,groupindex_,derivative_]:= Module[{CouplingGroupCorrect,isitgroup,lorentz,lorentzderivedhere,lorentzrestderivatives,repandflavors,repgroup,rep,Groups,CouplingGroup,indices,dummy,generators,GaugeFields,GaugeFieldstransformation,GaugeFieldindex,GaugeFieldswithlabels,replace,Fieldindex,Fieldindices,Fieldwithlabels,groupdecomp},
repandflavors =GetFields[field][Indices];
If[GetFields[field][Indices] ==Missing["KeyAbsent",field][Indices], Return[0]];
isitgroup[X_]:=Module[{assoc,group},
assoc = GetGaugeGroups[];
group = SelectFirst[Keys[assoc],MemberQ[assoc[#][Representations],X]&];
If[ToString[group] ==ToString[Missing["NotFound"]], {}, {X,group}]
];

repgroup = DeleteCases[isitgroup/@ repandflavors,{}];
If[repgroup == {}, Return[0]];
{rep,Groups} = Transpose[repgroup];

indices = First/@ Select[groupindex,MatchQ[#,Index[_,Alternatives@@rep]]&];

dummy= Table[Table[Unique["z"],Length[rep]],2];

generators =MapThread[CG[gen[#1],{#2,#3,#4}]&,{rep,dummy[[1]],indices,dummy[[2]]}] ;

GaugeFields = GetGaugeGroups[#][Field]&/@ Groups;
GaugeFieldstransformation = GetFields[#][Indices] &/@GaugeFields;
GaugeFieldindex = MapThread[Index[#1,#2]&,{dummy[[1]],Flatten[GaugeFieldstransformation]}];

GaugeFieldswithlabels = MapThread[Field[#1,Vector[derivative],{#2},{}]&,{GaugeFields,GaugeFieldindex}] ;

replace[indicelist_,replacementpart_]:=Replace[indicelist,Index[_,grp_]/;grp===replacementpart[[2]]:>replacementpart,{1}];
Fieldindex = MapThread[Index[#1,#2]&,{dummy[[2]],rep}];
Fieldindices = replace[groupindex,#] &/@ Fieldindex;

Fieldwithlabels = MapThread[Field[field,type,#,{}]&,{Fieldindices}] ;



groupdecomp = I*generators*GaugeFieldswithlabels*Fieldwithlabels;
Total[groupdecomp]

]



(* ::Subsection:: *)
(*Fieldstrength Tensors - Convert Kinetic Term into vector fields*)


ExpandFieldStrengthTensors[L_]:= Module[{Lsub,Lnopowers},
Lnopowers = L /.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];
Lsub = Lnopowers/. FieldStrength[f_,lorentz_,flavor_,index_] :> splitFS[FieldStrength[f,lorentz,flavor,index]];
Activate[Lsub]
]

splitFS[fieldstrength_]:= Module[{labels,name,vectorlabel,grouplabel,GaugeFields, IsNameGaugeField},
labels =Cases[fieldstrength, FieldStrength[x_,lorentz_,index_,derivatives_]:>{x,lorentz,index,derivatives} ,{0,Infinity}][[1]];
name = labels[[1]];
vectorlabel =labels[[2]];
grouplabel = labels[[3]];

GaugeFields = Key[Field]/@ Values[GetGaugeGroups[]];
IsNameGaugeField = MemberQ[GaugeFields,name];

If[IsNameGaugeField == True, FSGaugeField[name,vectorlabel,grouplabel],FSVectorField[name,vectorlabel,grouplabel]]
]

FSGaugeField[name_,vectorlabel_,grouplabel_]:=Module[{CouplingGroupCorrect,Abelianpart,repandflavors,isitgroup,assoc,group,repgroup,rep,Groups,CouplingGroup,dummy,structureConst,dummy1,dummy2,vector1,vector2,nonAbelianpart,FSsplit},
Abelianpart = -Field[name,Vector[vectorlabel[[1]]],grouplabel,{vectorlabel[[2]]/. Index[s_Symbol, type_] :> Index[Symbol[SymbolName[s] <> "p"], type]}] + Field[name,Vector[vectorlabel[[2]]],grouplabel,{vectorlabel[[1]]/. Index[s_Symbol, type_] :> Index[Symbol[SymbolName[s] <> "p"], type]}];

repandflavors =GetFields[name][Indices];
If[Length[repandflavors]==0, Return[Abelianpart]];
isitgroup[X_]:=Module[{},
assoc = GetGaugeGroups[];
group = SelectFirst[Keys[assoc],MemberQ[assoc[#][Representations],X]&];
If[ToString[group] ==ToString[Missing["NotFound"]], {}, {X,group}]
];
repgroup = DeleteCases[isitgroup/@ repandflavors,{}];
{rep,Groups} = Transpose[repgroup];

dummy= Table[Table[Unique["z"],Length[rep]],2];

structureConst =MapThread[CG[fStruct[#1],{#2,Index[#3,#5],Index[#4,#5]}]&,{Groups,grouplabel,dummy[[1]],dummy[[2]],rep}] ;

dummy1 = Table[ReplacePart[grouplabel,i->Index[dummy[[1]][[i]],rep[[i]]]],{i,Length[grouplabel]} ];
dummy2=  Table[ReplacePart[grouplabel,i->Index[dummy[[2]][[i]],rep[[i]]]],{i,Length[grouplabel]} ];

vector1 =Field[name,Vector[vectorlabel[[1]]],#,{}] &/@dummy1;
vector2 =Field[name,Vector[vectorlabel[[2]]],#,{}] &/@dummy2;

nonAbelianpart = Total[structureConst*vector1*vector2];
FSsplit =Abelianpart  + nonAbelianpart
]

FSVectorField[name_,vectorlabel_,grouplabel_]:=Module[{},
-Field[name,Vector[vectorlabel[[1]]],grouplabel,{vectorlabel[[2]]}] + Field[name,Vector[vectorlabel[[2]]],grouplabel,{vectorlabel[[1]]}]
]


(* ::Subsection:: *)
(*Function that hopefully works*)


Momenta[solution_]:=Module[{pairs,findPartner,inlist,newindex},
	pairs=Cases[solution,Metric[arg1_,arg2_]:>{arg1,arg2},{1}];
	findPartner[element_,pairs_]:= If[pairs[[1]]=== element, pairs[[2]],pairs[[1]]];
	solution/.Subsuperscript["p",index_,label_]:> Module[{},
	inlist = SelectFirst[pairs,MemberQ[#,index]&];
	If[inlist==Missing["NotFound"], Return[Subsuperscript["p",index,label]]];
	newindex = findPartner[index, inlist];
	Subsuperscript["p",newindex,label]]/. Metric@@inlist :> 1
];


position[term_,tolook_]:=Module[{init,inlist,inTermAndInput,lengthofterms,lengthofinput,Positionsincorrespondenc},
init[t_,l_]:= Transpose[Position[t,l]][[1]];
inlist =init[term,#]&/@ tolook;

inTermAndInput = Apply[Intersection,init[term,#]&/@ tolook];
lengthofterms = Length /@term[[Apply[Intersection,init[term,#]&/@ tolook]]];
lengthofinput = Length[tolook];
Positionsincorrespondenc = Flatten[Position[lengthofterms,lengthofinput]];
inTermAndInput[[Positionsincorrespondenc]] 
];


(* ::Section:: *)
(*Feynman Rules*)


(* ::Subsection:: *)
(*Prepare Full Lagrangian into an accessible term*)


PrepareLagrangian[L_]:= Expand[ExpandCovariantDerivatives[Expand[ExpandFieldStrengthTensors[L/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}]]/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}]]]] 


ExpandLagrangian[L_]:= Expand[Activate[PrepareLagrangian[L]]];


(* ::Subsection::Closed:: *)
(*Find the Count of each Field in Term*)


FindCountofFields[termold_]:=Module[{resultArray,groupedFieldData,term,fieldPositions,fieldData},

term = termold /.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];
fieldPositions=Position[term,Field[_,___],{0,Infinity}];


fieldData=Table[Module[{currentField,fieldLabel,isInsideBar},
currentField=Extract[term,pos];

fieldLabel=currentField[[1]];

isInsideBar=(Length[pos]>1&&Head[Extract[term,Most[pos]]]===Bar);
{fieldLabel,isInsideBar}],{pos,fieldPositions}];


groupedFieldData=GatherBy[fieldData,Identity];


resultArray=Map[{First[First[#]],Length[#],Last[First[#]]}&,groupedFieldData]
]


(* ::Subsection::Closed:: *)
(*Change the code intern labels into the ones given externally*)


AssignExternalFields[newinput_,ext_]:=Module[{deltacreater,lin,lext, clebschgordens,listwithgroup,listwithoutcharge,fieldname,typeext,flavorext,lorentz,listwithchargeandnogroup,Z},
deltacreater[aa_,bb_]:= CG[del[aa],bb];
{fieldname,typeext,flavorext,lorentz} = Cases[newinput, Field[x_,type_,index_,lorentz_]:>{x,type, index,lorentz} ,{0,Infinity}][[1]];
Z[x_] := x;
If[Head[newinput]==Bar,Z = Bar];

ext/. Z[Field[f_,typein_,flavorint_,index___]]:>  Module[{},

listwithchargeandnogroup =GatherBy[Join[flavorext,flavorint],#[[2]]&];
listwithgroup = Map[{#[[1,2]],#}&,listwithchargeandnogroup];
clebschgordens = Apply[Times,Apply[deltacreater,#] &/@ listwithgroup];
If[!(Head[typeext] ===Vector), Z[Field[EXT[f],typein,flavorext,index]]*clebschgordens,
lin = typein[[1]];
lext = typeext[[1]];
Z[Field[EXT[f],typeext,flavorext,index]]*clebschgordens*Metric[lin,lext]
]
]
];


(* ::Subsection:: *)
(*Feynman Rules Finder*)


FeynmanRule[L_,input_]:=Module[{cleaninput,preresult2,preresult1,occur,termsinL,fieldOccurrences,groupedFields,fieldinput,pos,interactions,GeneralRules,listinput,listinputbar,listofnamesoffields,countfieldsinput,occurencelist,result,signsfrompermutations,shapedLagrangian,interactionswithsign,filteredList,listWithSubscripts},
shapedLagrangian = ExpandLagrangian[L]/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}] ;
termsinL = If[Head[shapedLagrangian]===Plus,List@@shapedLagrangian,{shapedLagrangian}];

occur[term_]:=Module[{dummy},
  dummy=Cases[term, Bar[Field[x_,___]]:>{x,True},{0,Infinity}];
If[Length[dummy]==1, Return[dummy[[1]]]];
Cases[term, Field[x_,___]:>{x,False} ,{0,Infinity}][[1]]
];


fieldOccurrences= occur/@ input;
groupedFields=GatherBy[fieldOccurrences,Identity];
fieldinput=Map[{First[First[#]],Length[#],Last[First[#]]}&,groupedFields];

If[Length[termsinL]==1, If[Sort[fieldinput]===Sort[FindCountofFields[termsinL]],pos = {1},Return[Print["Not matching"]]], pos = Position[Sort /@ FindCountofFields/@termsinL,Sort[fieldinput]]];



listinput =Cases[input, Field[x_,type_,index_,lorentz_]:>x ,{0,Infinity}];
listinputbar=Table[Boole[!FreeQ[input[[i]], Bar]],{i,Length[input]}];
listofnamesoffields = Table[If[listinputbar[[i]]==1, OverBar[listinput[[i]]],listinput[[i]]],{i,Length[listinputbar]}];
countfieldsinput=Module[{counts=<||>},Map[Function[x,counts[x]=Lookup[counts,x,0]+1],listofnamesoffields]];
occurencelist =MapThread[Subscript,{listofnamesoffields,countfieldsinput}];
listWithSubscripts = Table[Subscript[listofnamesoffields[[i]], i], {i, Length[listofnamesoffields]}];


cleaninput[term_]:=Module[{dummy},
  dummy=Cases[term, Bar[Field[x___]]:>Bar[Field[x]],{0,Infinity}];
If[Length[dummy]==1, Return[dummy[[1]]]];
Cases[term, Field[x___]:> Field[x],{0,Infinity}][[1]]
];



Print["There are ",Length[pos], " Interaction terms"];
interactions =termsinL[[Flatten[pos]]];

typeslist =GetFields[#][Type]&/@listinput;
positionsofferms=Position[typeslist,x_/;MemberQ[{Fermion,Ghost},x]];filteredList =Extract[occurencelist,positionsofferms];



signsfrompermutations =ComputeGrassmannSign[filteredList,# ]&/@interactions;
interactionswithsign = interactions*(Power[-1,#]&/@signsfrompermutations);
(*Print[interactions //RelabelIndices //NiceForm];*)
GeneralRules = DerivativesToMomenta[Feynman[Total[interactionswithsign],D]];
(*Print["Solution done, fitting to external legs"];*)


Dictionaryofextfields=AssociationThread[occurencelist,cleaninput/@input];

preresult1 = RelabelIndices[GeneralRules,Unique -> False] /. Bar[Field[EXT[l1_],l2_,l3_,l4_]] :> Module[{},
AssignExternalFields[Dictionaryofextfields[l1],Bar[Field[l1,l2,l3,l4]]]];

preresult2 = preresult1/. Field[EXT[l1_],l2_,l3_,l4_] :> Module[{}, 
If[MatchQ[l1,Subscript[OverBar[_],_]],Field[EXT[l1],l2,l3,l4],AssignExternalFields[Dictionaryofextfields[l1],Field[l1,l2,l3,l4]]]
];


result = If[Head[preresult2]===Plus,List@@preresult2,{preresult2}];

Simplify[Contract[Total[result]]]/. Normal[AssociationThread[occurencelist,listWithSubscripts]]

]


(* ::Subsection::Closed:: *)
(*ComputeGrassmannSign*)


ComputeGrassmannSign[inputlist_,Lterm_]:= Module[{lfj,emptyL,nlist,listfieldn,calclistpre,annotateWithSubscript,calclist,countAdjacentSwaps},
lfj= GatherFieldsByType[Lterm,Fermion];
If[Length[lfj]==0, Return[0]];
emptyL = Lterm/.Bar[Field[f_,Fermion,type2_,indices_]]:>Bar[Field[f,Fermion,{},{}]]/.Field[f_,Fermion,type2_,indices_]:>Field[f,Fermion,{},{}]/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}];
nlist = CountFieldOccurrences[#,emptyL] &/@ lfj;
listfieldn = Transpose[{lfj,nlist}];
calclistpre = listfieldn/.Bar[Field[x_,Fermion,f__]] :> OverBar[Field[x,Fermion,f]] /.Field[x_,___]:> x;


calclist = Flatten[Table[Subscript[#[[1]],k],{k,#[[2]],1,-1}]&/@calclistpre];
countAdjacentSwaps[list1_,list2_]:=Module[{n=Length[list1],swaps=0,arr,temp}, arr=Ordering[Ordering[list2/. Thread[list1->Range[n]]]];
Do[Do[If[arr[[j]]>arr[[j+1]],temp=arr[[j]];
arr[[j]]=arr[[j+1]];
arr[[j+1]]=temp;
swaps++],{j,n-i}],{i,n}];
swaps];

countAdjacentSwaps[inputlist,calclist]

]


(* ::Subsection:: *)
(*Find all Rules*)


FeynmanRuleFinderAll[L_]:=Module[{SubstituteFieldIndicesules,expandedL,termsinLold,termsinL,numberofFieldsinL,gathersimilarterms,Lnew,numberofFieldsinLnew,outputlabels,preresult,numberofFieldsinLreallynew},
SubstituteFieldIndicesules = ReplaceCouplingConvention[L];
expandedL = ExpandLagrangian[L]/.Power[x__,k_Integer?Positive]:>Inactive[Times]@@Table[x,{k}]/.SubstituteFieldIndicesules;
termsinLold = List@@expandedL;
termsinL = Select[termsinLold,(Count[#,Field[___],{0,Infinity}]>=3)&];

 numberofFieldsinL = FindCountofFields/@termsinL;
gathersimilarterms = GatherBy[Range[Length[numberofFieldsinL]],numberofFieldsinL[[#]]&];
Lnew= Total[termsinL[[#]]]&/@gathersimilarterms;
numberofFieldsinLnew = numberofFieldsinL[[First[#]]]&/@ gathersimilarterms;
numberofFieldsinLreallynew=CustomSort/@numberofFieldsinLnew;
outputlabels = Listmaker[#]&/@numberofFieldsinLreallynew;
preresult = Block[{Print=(Null&)},MapThread[FeynmanRule, {Lnew,outputlabels}]];
ReplaceFieldIndices/@(RenameDummyIndices/@preresult) //RelabelIndices//Simplify 


]


generateFields[{name_,n_,wrap_},Counter_]:=Module[{fields},fields=Table[Field[name,PoincareStructure[name,i+Counter],GroupStructure[name,i+Counter],{}],{i,1,n}];
If[wrap,Bar/@(fields/. Subscript[label_,k_]:> Subscript[OverBar[label],k]),fields]]

PoincareStructure[Name_,number_]:= Module[{},
type =  GetFields[Name][Type];
If[type === Vector, 
Vector[Index[Unique["z"], Lorentz]]
,type ]
]

GroupStructure[Name_,number_]:= Module[{},
type =  GetFields[Name][Indices];
indices = Index[Unique["z"],#]&/@ type
]

Listmaker[inputList_]:=Module[{counter=0,result={}},Do[With[{fields=generateFields[elem,counter]},AppendTo[result,fields];
counter+=elem[[2]];],{elem,inputList}];
Flatten[result]]


ReplaceFieldIndices[expr_]:=Module[{rules,assoc},rules=Cases[expr,Field[EXT[Subscript[F_,n_]],args___]:>Module[{zlist},zlist=Cases[{args},Index[z_Symbol,_],\[Infinity]][[All,1]];
Thread[zlist->n]],\[Infinity]];
rules=Flatten[rules];
assoc=Association[rules];
expr/. Index[z_,t_]:>If[KeyExistsQ[assoc,z],Index[Symbol["d$$"<>ToString[assoc[z]]],t],Index[z,t]]]

RenameDummyIndices[expr_]:=Module[{rules,fromList,toList},
fromList=DeleteDuplicates[Cases[expr,s_Symbol/;StringStartsQ[SymbolName[s],"d$$"],\[Infinity]]];
toList=Table[Unique["z"],{Length[fromList]}];

rules=Thread[fromList->toList];

expr/. rules]


ReplaceCouplingConvention[expr_]:= Module[{GaugeFields,real,complex,rescaleRulesreal,rescaleRulescomplex,rules},
GaugeFields = Keys@Select[GetFields[],(#[Type]===Vector)&];
real = Select[List@@expr,Module[{fs,hasField},fs=Cases[#,FieldStrength[x_,___]^2:>x,\[Infinity]];
hasField=!FreeQ[#,Field];
Length[fs]==1&&MemberQ[GaugeFields,First[fs]]&&!hasField]&];
complex = Select[List@@expr,MatchQ[#,_*Bar[FieldStrength[x_,___]]*FieldStrength[y_,___]/;x===y&&MemberQ[GaugeFields,x]]&];
rescaleRulesreal=Cases[real,-1/4*Power[coupling_,-2]*Power[FieldStrength[X_,___],2]:>({Field[X,args___]:>coupling*Field[X,args]})];
rescaleRulescomplex=Cases[complex,-1/2*Power[coupling_,-2]*Bar[FieldStrength[X_,___]]*FieldStrength[X_,___]:>({Field[X,args___]:>coupling*Field[X,args]})];
rules = Flatten[Join[rescaleRulesreal,rescaleRulescomplex]]
]


fieldOrder[Fermion]:=0;
fieldOrder[Vector]:=1;
fieldOrder[Scalar]:=2;
fieldOrder[Ghost]:=3;

CustomSort[list_List]:=SortBy[list,{fieldOrder[GetFields[#[[1]]][Type]]&,#[[1]]&,-Boole[#[[3]]]&}];


(* ::Subsection:: *)
(*Function that automatically checks input*)


FeynmanRules[x_,a_:Automatic]:=If[a===Automatic,FeynmanRuleFinderAll[x],FeynmanRule[x,a]]
