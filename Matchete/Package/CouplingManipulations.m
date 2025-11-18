(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`CouplingManipulations`*)


(* ::Subtitle:: *)
(*Paclet for manipulations of coupling coefficients*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["PrintEffectiveCouplings"]
PackageExport["ReplaceEffectiveCouplings"]
PackageExport["Rules"]
PackageExport["DummyCoefficients"]
PackageExport["EffectiveCouplingSymbol"]
PackageExport["IntroduceEffectiveCouplings"]


PackageExport["ShiftRenCouplings"]
PackageExport["MapEffectiveCouplings"]
PackageExport["SortByEFTOrder"]
PackageExport["KeepTrivalReplacements"]
PackageExport["AppendEffectiveCouplingsDefs"]


PackageExport["ReplaceInLagrangian"]


PackageScope["ExportSMEFT"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["IntroduceEffectiveMasses"]
PackageScope["ResetTempCouplings"]
PackageScope["ResetInternalCouplings"]


PackageScope["OverrideDuplicateCouplingCheck"]
PackageScope["SubstituteCoefficients"]


PackageScope["IntroduceDummyCoefficients"]
PackageScope["$IntCouplingRules"]


PackageScope["CouplingProduct"]
PackageScope["MassFunction"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


IntroduceEffectiveCouplings::usage =    "IntroduceEffectiveCouplings[L] replaces operator coefficients with newly defined, generic couplings to compress the structure of the Lagrangian. Use ReplaceEffectiveCouplings to resubstitute them or PrintEffectiveCouplings to show a list of the couplings and their definitions.";
PrintEffectiveCouplings::usage =        "PrintEffectiveCouplings[L] prints a list of effective couplings defined by EOMSimplify appearing in the input Lagrangian L, along with their definitions. With the option Rules\[Rule] True, the function will return the replacement list for these couplings, instead.";
ReplaceEffectiveCouplings::usage =      "ReplaceEffectiveCouplings[L] resubstitutes effective couplings defined by EOMSimplify or by IntroduceEffective couplings. It can be called with a list of labels as second argument, allowing to specify which couplings to replace.";
Rules::usage =                          "Rules is an option to specify whether output should be a list of replacement rules instead of printed.";
EffectiveCouplingSymbol::usage =        "EffectiveCouplingSymbol an option of IntroduceEffectiveCouplings and is the symbol used for automatically defined effective couplings. The standard option is EffectiveCouplingSymbol -> \"C\"";
OverrideDuplicateCouplingCheck::usage = "OverrideDuplicateCouplingCheck is an option for IntroduceEffectiveCouplings (standard is False). If set to true, the function will introduce effective couplings even if an identical one has been defined previously.";
ShiftRenCouplings::usage              = "ShiftRenCouplings[lag] introduces effective couplings for all renormalizable operators, and expresses the matching correction to the coefficients of the higher-dimensional operators in terms of these effective couplings, effectively shifting the matching corrections from the renormalizable couplings to the higher-dimensional coefficients. Only tree-level generated coefficients are shifted. Gauge couplings (and thus also gauge boson kinetic terms) are never shifted. To shift gauge couplings see MapEffectiveCouplings.";
IntroduceEffectiveMasses::usage =       "IntroduceEffectiveMasses[L] replaces mass terms with newly defined, generic mass couplings to compress the structure of the Lagrangian. Use ReplaceEffectiveCouplings to resubstitute them or PrintEffectiveCouplings to show a list of the couplings and their definitions.";
MapEffectiveCouplings::usage =          "MapEffectiveCouplings[L0,LTarget] attempts to relate the two input Lagrangians, expressing the effective couplings in LTarget through the matching coefficients in L0. The option SortByEFTOrder (default is True) specifies whether the output is sorted by power-counting.";
SortByEFTOrder::usage =                 "SortByEFTOrder is an option for MapEffectiveCouplings, specifying whether the output rules are sorted by power-counting. The default is True, but setting it to False can lead to minor speed-up for very large Lagrangians.";


ReplaceInLagrangian::usage= "ReplaceInLagrangian[lag, rule(s)] performs replacements in a Lagrangian expression while respecting Einstein summation convention, i.e. w/ proper relabeling of repeated indices to avoid collisions.";


ExportSMEFT::usage= "ExportToSmelli[<filename>, <SMEFT lagrangian>, Options] creates a file with the Warsaw basis SMEFT coefficients at the current working directory. The json file option contains all the functional dependence of the Warsaw basis coefficients on the UV parameters, in such a way that it can be imported in jelli.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Operator Coefficients*)


(* ::Subsection:: *)
(*Handling redefined effective couplings*)


(* ::Text:: *)
(*Generating new effective couplings parametrizing the coefficients in Lagrangians *)


(* ::Subsubsection::Closed:: *)
(*Data structures holding the information*)


ResetTempCouplings[]:= (
		$TempCouplings={};
		$TempCouplingRules={};
	);
ResetTempCouplings[];


ResetInternalCouplings[]:= (
		$InternalCouplings={};
		$IntCouplingRules={};
	);
ResetInternalCouplings[];


(* ::Subsubsection::Closed:: *)
(*Replace couplings  helpers*)


(* check if two operators are the same (or the hc of each other) *)
SameOperatorQ[(c1_:1) o1_Operator,(c2_:1) o2_Operator] := Module[{baro1,o1r,o2r},
	o1r = c1 RelabelIndices @ o1;
	o2r = c2 RelabelIndices @ o2;
	baro1 = Operator @ Bar @ RelabelIndices @ OperatorToNormalForm @ (c1 o1);
	Return[SameQ[o1r,o2r] || SameQ[baro1, o2r]]
]


(* check if an operator is hermitian *)
BarredOpQ[x_Operator] := Expand[RelabelIndices @ x - Operator @ Bar @ RelabelIndices @ OperatorToNormalForm @ x]===0
MHermitianQ[x_]:= GreensSimplify[x - Bar @ x] === 0


(* hacky check to see if a rule has already been defined - SameQ doesn't work reliably on RuleDelayed objects *)
SameRuleQ[r_, q_]:=Module[{rx=r,qx=q,tlbl,rd},
	rx = rx /.RuleDelayed[Coupling[a_,f___],b_]:> rd[Coupling[tlbl,f],RelabelIndices@b];
	qx = qx /.RuleDelayed[Coupling[a_,f___],b_]:> rd[Coupling[tlbl,f],RelabelIndices@b];
	SameQ[rx, qx]
]


(* build a coupling symbol and appropriate NiceForm *)
CouplingNameFromOperator[op:(_AtomicOp|_CompOp), couplingstring_]:=Module[
	{fieldtally,opdevcount,baseString,incr=1,prettylabel,temp},
		(* fields appearing - we want FS objects to appear here as well hence convert to Operator *)
		fieldtally = Tally@Cases[(*Operator@RelabelIndices@AtomicToNormalForm*)AtomicToOperatorForm@op,Field[f_,__]|FieldStrength[f_,__]:> f,Infinity];
		(* number of derivatives minus (twice) number of FS objects *)
		opdevcount = op[[Switch[Head@ op, AtomicOp, 1, CompOp, 2],2]] - 
					2*Length@Cases[(*Operator@RelabelIndices@AtomicToNormalForm*)AtomicToOperatorForm@op,_FieldStrength,Infinity] - 
					2*Length@Cases[(*Operator@RelabelIndices@AtomicToNormalForm*)AtomicToOperatorForm@op,EoM[Field[_,_Vector,___]],Infinity];

		baseString = couplingstring <> StringJoin[ 
					ToString/@ Flatten[ 
								If[MatchQ[#,{_,k_/;k>1}], #, {#[[1]]}]& /@ Join[fieldtally, {If[opdevcount>0,{"D",opdevcount},Nothing]}] ] ];

		If[Defined[Symbol@baseString],
			While[Defined[Symbol[baseString<>ToString[incr]]], incr++];
			temp = Symbol[baseString<>ToString[incr]];
			prettylabel = StandardForm[
							Subsuperscript[StandardForm[couplingstring],
									  Row@ Join[If[#[[2]]>1,Superscript[StandardForm[#[[1]]],#[[2]]],StandardForm[#[[1]]]]&/@fieldtally,
									            {If[#>0,If[#>1,Superscript[StandardForm[D],#],StandardForm[D]],Nothing]&@opdevcount}],Row[{"(",incr,")"}]]]
		,
			temp = Symbol@baseString;
			prettylabel = StandardForm[
							Subscript[StandardForm[couplingstring],
									  Row@ Join[If[#[[2]]>1,Superscript[StandardForm[#[[1]]],#[[2]]],StandardForm[#[[1]]]]&/@fieldtally,
									            {If[#>0,If[#>1,Superscript[StandardForm[D],#],StandardForm[D]],Nothing]&@opdevcount}]]]
		];
	
	{temp, prettylabel}
]


(* ::Subsubsection::Closed:: *)
(*ToEffectiveCoupling (for operators)*)


ToEffectiveCoupling::IndexProblem="The expression `1` contains a mixture of diagonal and non-diagonal indices and cannot be treated.";

(* option Superleading specifies if this coupling acts as a replacement of superleading terms or just as an internal coupling
   option Internal specifies if this coupling is saved to the temporary internal list or the public one *)
Options @ ToEffectiveCoupling = { Superleading -> True, Internal -> True, EffectiveCouplingSymbol -> "C", HermitianTerm -> False, OverrideDuplicateCouplingCheck -> False };

(* exp is only the coupling, this expression returns only the effective coupling as well, not the product with the operator *)
ToEffectiveCoupling[exp_, operator_Operator, OptionsPattern[] ] :=Module[
		{superleading,internal,power,diagIndex,nDiagIndex,index,indexType,\[Lambda]EFT, outCoupling,temp,inds, pinds, ruleLHS,ruleRHS,hermite, derivativePhase = 1, preRHS,
		preExisting, candidates, candidateLabels, hermitianTerm, pOrder, permConversion, permConversionInv, symmetries = {},
		tbdassociation, couplingsSameProperty, ignoreduplicates, isdiagonal = False, opFlavProp, openIndices,
		baseString, incr=1, opdevcount, fieldtally, prettylabel, couplingstring = "C",
		IndexNotRelabeled, indCounter1=1}
	,
		(* get options *)
		hermitianTerm = TrueQ @ OptionValue @ HermitianTerm;
		superleading = TrueQ @ OptionValue @ Superleading;
		internal = TrueQ @ OptionValue @ Internal;
		ignoreduplicates = TrueQ @ OptionValue @ OverrideDuplicateCouplingCheck;

		If[Head @ OptionValue[ EffectiveCouplingSymbol ] === String, couplingstring = OptionValue[ EffectiveCouplingSymbol ] ];

		(* build the coupling name from the operator field content and number of derivatives *)
		opdevcount = Length@Flatten@Cases[operator, Field[__,k_/;Length[k]>0]|FieldStrength[__,k_/;Length[k]>0] :> k,Infinity] +
					2 Length@Flatten@Cases[operator, EoM[Field[_,Scalar,__]]|EoM[Bar @ Field[_,Scalar,__]],Infinity] +
					Length@Flatten@Cases[operator, EoM[Field[_,Except@ Scalar,__]]|EoM[Bar @ Field[_,Except@ Scalar,__]],Infinity];
		fieldtally = Tally@Cases[operator,Field[f_,__]|FieldStrength[f_,__]:> f,Infinity];
		baseString = couplingstring <> StringJoin[ ToString/@ Flatten[ If[MatchQ[#,{_,k_/;k>1}], #, {#[[1]]}]& /@ Join[fieldtally, {If[opdevcount>0,{"D",opdevcount},Nothing]}] ] ];
		If[Defined[Symbol@baseString],
			While[Defined[Symbol[baseString<>ToString[incr]]], incr++];
			temp = Symbol[baseString<>ToString[incr]];
			prettylabel = StandardForm[
							Subsuperscript[StandardForm[couplingstring],
									  Row@ Join[If[#[[2]]>1,Superscript[StandardForm[#[[1]]],#[[2]]],StandardForm[#[[1]]]]&/@fieldtally,
									            {If[#>0,If[#>1,Superscript[StandardForm[D],#],StandardForm[D]],Nothing]&@opdevcount}],Row[{"(",incr,")"}]]]
		,
			temp = Symbol@baseString;
			prettylabel = StandardForm[
							Subscript[StandardForm[couplingstring],
									  Row@ Join[If[#[[2]]>1,Superscript[StandardForm[#[[1]]],#[[2]]],StandardForm[#[[1]]]]&/@fieldtally,
									            {If[#>0,If[#>1,Superscript[StandardForm[D],#],StandardForm[D]],Nothing]&@opdevcount}]]]
		];

		(* determine the open indices that the coupling needs to have *)
		nDiagIndex = FindOpenIndices @ First[List @@ (OperatorToNormalForm[operator, Unique->True]+ Nothing)];
		(* check the prefactor for diagonal indices *)
		diagIndex = FindDiagonalIndices @ exp;
		(* drop diagonal indices that are not contracted into to operator as they should be summed over only inside the coefficient *)
		diagIndex = Intersection[diagIndex, DeleteDuplicates@Cases[operator, _Index, All]];

		index = Join[nDiagIndex,diagIndex];
		indexType = Last /@ index;
		isdiagonal = Join[ ConstantArray[False, Length@nDiagIndex] , ConstantArray[True, Length@diagIndex] ];

		(* determine the phase coming from derivatives *)
		derivativePhase = (operator /. {Field[A__, ders_List] :> I^Length@ders Field[A, ders],EoM[Field[f_,Fermion, A___,ders_List]] :> I EoM @ Field[f, Fermion, A, ders]})/operator;

		(* power-counting of the object depending on context *)
		If[superleading,
			power = 4 - OperatorDimension @ operator
		,
			power = GetMinOrder @ exp
		];

		(* the coupling as it appears in the Lagrangian term *)
		outCoupling = derivativePhase Coupling[temp, index, power];

		(* we need this below inside and out of the next If clause *)
		opFlavProp = LookupOperatorFlavorProperties @ operator;
		openIndices = FindOpenIndices @ First[List @@ (OperatorToNormalForm[operator,Unique->True]+ Nothing)];
		(* if this has indices, determine the conversion to the IBPSimplify database and its inverse *)
		If[Length @ openIndices > 0,
			With[{operatorPatternIndices = Last @ MatchOperatorPatterns @ OperatorToNormalForm[ operator, Unique->True]},
				permConversion = FindPermutationOrder[operatorPatternIndices, openIndices];
				permConversionInv = FindPermutationOrder[ openIndices, operatorPatternIndices ];
			];
		];

		(*
		(* is the operator hermitian? *)
		If[ superleading || internal || !hermitianTerm ,
			(* for the automated routes or if this is not overridden by HermitianTerm -> True, just determine it from the operator *)
			hermite = MHermitianQ[derivativePhase NormalForm @ operator]
			(*hermite = {2,1}*)
		,
			(* called by IntroduceEffectiveCouplings with HermitianTerm -> True, need to determine symmetries *)
			If[Length @ openIndices > 0,
				(* take from IBPSimplify database and convert to the ordering found here *)
				pOrder = permConversionInv[[ opFlavProp[ConjugateIndexExchange] ]][[ permConversion ]];

				(* if the list is empty, this operator is simply hermitian without further constraints *)
				If[Length @ pOrder > 0,
					hermite = pOrder
				,
					hermite=True
				]
			,
				(* actually no open indices so just hermitian right away *)
				hermite=True
			];
		];
		*)

		(* is the operator hermitian? *)
		If[ superleading || internal || !hermitianTerm ,
			(* for the automated routes or if this is not overridden by HermitianTerm -> True, just determine it from the operator *)
			If[Length @ openIndices > 0,
				(* take from IBPSimplify database and convert to the ordering found here *)
				If[opFlavProp[ConjugateIndexExchange]===False,
					hermite = False
					,
					pOrder = permConversionInv[[ opFlavProp[ConjugateIndexExchange] ]][[ permConversion ]];
					(* if the list is empty, this operator is not hermitian *)
					If[Length @ pOrder > 0,
						hermite = pOrder
					,
						hermite = False
					]
				]
			,
				(* check hermiticity if coupling/operator do not contain open indices *)
				hermite = MHermitianQ[derivativePhase OperatorToNormalForm[operator,Unique->True]]
			]
		,
			(* called by IntroduceEffectiveCouplings with HermitianTerm -> True, need to determine symmetries *)
			If[Length @ openIndices > 0,
				(* take from IBPSimplify database and convert to the ordering found here *)
				pOrder = permConversionInv[[ opFlavProp[ConjugateIndexExchange] ]][[ permConversion ]];

				(* if the list is empty, this operator is simply hermitian without further constraints *)
				If[Length @ pOrder > 0,
					hermite = pOrder
				,
					hermite=True
				]
			,
				(* actually no open indices so just hermitian right away *)
				hermite=True
			];
		];

		symmetries = Association[Range[Length@index]->1];
		If[Length @ openIndices > 0,
			Module[{symrules1,symrules2, ineqp, symrulesrhs, replacedexp},
				(* we begin by finding the symmetries that we know from the operator *)
				symrules1 = (permConversionInv[[ #[[1]] ]][[ permConversion ]] -> #[[2]])&/@(opFlavProp[Symmetries]);
				symrulesrhs = (#[[1]])&/@symrules1;
				(* now we take the inequivalent permutations, delete the ones that are already in the symmetries list *)
				ineqp = Select[ ( permConversionInv[[ # ]][[ permConversion ]] )&/@ opFlavProp[InequivalentPermutations], !MemberQ[symrulesrhs,#]& ];
				(* iterate over all of them and check if they are actual symmetries *)
				symrules2 = Table[
					replacedexp = exp/.((#[[1]]->#[[2]])&/@DeleteCases[Transpose[{openIndices , openIndices[[cand]]}], {here_, here_} ]);
					If[RelabelIndices @ Contract[(replacedexp - exp)(OperatorToNormalForm[operator,Unique->True])] === 0,cand -> 1,
						If[RelabelIndices @ Contract[(replacedexp + exp)(OperatorToNormalForm[operator,Unique->True])] === 0, cand -> -1, Nothing]],
				{cand, ineqp}];
				symmetries = Association@@Join[symrules1,symrules2]
			];
		];

		(* indices for the replacement rules *)
		indCounter1 = 1;
		(*
		inds = Symbol["i" <> ToString[#]]& /@ Range@Length@index;
		pinds = Pattern[#, _]& /@ inds;
		*)
		inds = Index[Symbol["i" <> ToString[indCounter1++]],#]&/@ (Last/@index);
		pinds = IndexToPattern[inds];

		(* setting up the rules *)
		ruleRHS = exp /. Thread[index->inds];
		(* we need to temporarly rename the open diagonal indices, otherwise they get relabeled below *)
		ruleRHS = ruleRHS /. (dInd:(Alternatives@@(inds[[Count[isdiagonal,False]+1;;]])):> (dInd/.Index->IndexNotRelabeled));
		ruleLHS = Coupling[temp, pinds, power];

		preRHS = With[{RHS = RelabelIndices[ruleRHS, Unique->True] derivativePhase^-1}, ruleLHS :> (RelabelIndices[RHS,Unique -> True]/.IndexNotRelabeled->Index)];

		If[!ignoreduplicates,
			(* check for any couplings with the identical properties to the one we are about to define *)
			tbdassociation = CouplingAssociationEntry[power, indexType, hermite, symmetries, isdiagonal];
			couplingsSameProperty = Cases[{Keys@$CouplingAssociation,Values@$CouplingAssociation}\[Transpose],{val_, tbdassociation} :> val]
		];

		If[superleading == True ||(superleading == False && internal == False),
			(* this branch is for couplings visible to the user *)

			If[!ignoreduplicates,
				(* check we do not define duplicate couplings *)
				candidates =  Cases[$TempCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ]];
				candidateLabels = Cases[$TempCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ] :> f ];
				preExisting = (SameRuleQ[preRHS,#]& /@ candidates)
			];

			If[Or@@preExisting && !ignoreduplicates,
				(* the coupling has already be defined so figure out which one it is and use it instead *)
				outCoupling = derivativePhase Coupling[candidateLabels[[FirstPosition[preExisting,True][[1]]]], index, power]
			,
				(* does not exist so define it *)
				DefineCoupling[temp, Indices -> indexType, EFTOrder -> power, SelfConjugate -> hermite, Symmetries -> SymmetryOverride[symmetries], DiagonalCoupling -> isdiagonal];
				AppendTo[$TempCouplings,temp];
				AppendTo[$TempCouplingRules, preRHS];
				(* define the NiceForm of this new coupling *)
				AppendTo[LabelsNiceForm[Coupling],temp->prettylabel];
				(*Format[Coupling[temp, in_,_], NiceForm] := UpDownIndices[prettylabel,in];
				If[hermite===False, Format[Bar @ Coupling[temp, in_,_], NiceForm] := UpDownIndices[OverBar @ prettylabel, Bar /@ in]]*)
			]

		,
			(* this branch is for couplings invisible to the user *)

			If[!ignoreduplicates,
				(* check we do not define duplicate couplings *)
				candidates = Cases[$IntCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] ]];
				candidateLabels = Cases[$IntCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ] :> f ];
				preExisting = (SameRuleQ[preRHS,#]&/@ candidates)
			];

			If[Or@@preExisting && !ignoreduplicates,
				(* the coupling has already be defined so figure out which one it is and use it instead *)
				outCoupling = derivativePhase Coupling[candidateLabels[[FirstPosition[preExisting,True][[1]]]], index, power]
			,
				(* does not exist so define it *)
				DefineCoupling[temp, Indices -> indexType, EFTOrder -> power, SelfConjugate -> hermite, DiagonalCoupling -> isdiagonal];
				AppendTo[$InternalCouplings,temp];
				AppendTo[$IntCouplingRules, preRHS];
				(* define the NiceForm of this new coupling *)
				AppendTo[LabelsNiceForm[Coupling],temp->prettylabel];
				(*Format[Coupling[temp, in_,_], NiceForm] := UpDownIndices[prettylabel,in];
				If[hermite===False, Format[Bar @ Coupling[temp, in_,_], NiceForm] := UpDownIndices[OverBar @ prettylabel, Bar /@ in]]*)
			]

		];
		Return[outCoupling]
	]


(* ::Subsubsection::Closed:: *)
(*ToEffectiveCoupling (for atomics)*)


(* exp is only the coupling, this expression returns only the effective coupling as well, not the product with the operator *)
ToEffectiveCoupling[exp_, operator:(_CompOp|_AtomicOp), OptionsPattern[] ] :=Module[
		{superleading,internal,power,diagIndex,nDiagIndex,index,indexType,\[Lambda]EFT, outCoupling,temp,inds, pinds, ruleLHS,ruleRHS,hermite, derivativePhase = 1, preRHS,
		preExisting, candidates, candidateLabels, hermitianTerm, pOrder, permConversion, permConversionInv, symmetries = {},
		tbdassociation, couplingsSameProperty, ignoreduplicates, isdiagonal = False, opFlavProp, openIndices,
		baseString, prettylabel, couplingstring = "C",opdevcount,
		IndexNotRelabeled, indCounter1=1, preFactor}
	,
		(* get options *)
		hermitianTerm = TrueQ @ OptionValue @ HermitianTerm;
		superleading = TrueQ @ OptionValue @ Superleading;
		internal = TrueQ @ OptionValue @ Internal;
		ignoreduplicates = TrueQ @ OptionValue @ OverrideDuplicateCouplingCheck;

		(* if an object is already just a coupling (times a number), don't replace it - as long as OverrideDuplicateCouplingCheck is switched off *)
		If[(MatchQ[exp, (c_:1) x_Coupling /;NumericQ[c]]|| MatchQ[exp, (c_:1) HoldPattern[Bar[x_Coupling]] /;NumericQ[c]]) && !ignoreduplicates, Return[exp]];

		If[Head @ OptionValue[ EffectiveCouplingSymbol ] === String, couplingstring = OptionValue[ EffectiveCouplingSymbol ] ];

		{temp, prettylabel} = CouplingNameFromOperator[operator, couplingstring];

		(* determine the open indices that the coupling needs to have *)
		nDiagIndex = Last@ operator;
		nDiagIndex = First/@ DeleteCases[Tally[nDiagIndex],{_,Except[1]}]; (* drop indices repeated within operator *)
		(* check the original coupling for diagonal indices and delete duplicates *)
		diagIndex = DeleteDuplicates@ FindDiagonalIndices @ exp;
		(* drop diagonal indices that are not contracted into to operator as they should be summed over only inside the coefficient *)
		diagIndex = Intersection[diagIndex, DeleteDuplicates@Cases[operator, _Index, All]];

		index = Join[nDiagIndex,diagIndex];
		indexType = Last /@ index;
		isdiagonal = Join[ ConstantArray[False, Length@nDiagIndex] , ConstantArray[True, Length@diagIndex] ];

		(* determine the phase coming from derivatives *)
		opdevcount = operator[[Switch[Head@ operator, AtomicOp, 1, CompOp, 2],2]] - 2*Length@Cases[(*Operator@RelabelIndices@AtomicToNormalForm*)AtomicToOperatorForm@operator,_FieldStrength,Infinity];
		derivativePhase = I^opdevcount;

		(* power-counting of the object depending on context *)
		If[superleading,
			power = 4 - OperatorDimension @ AtomicToNormalForm @ operator
		,
			power = GetMinOrder @ exp
		];

		(* the coupling as it appears in the Lagrangian term *)
		outCoupling = derivativePhase Coupling[temp, index, power];

		(* we need this below inside and out of the next If clause *)
		opFlavProp = LookupOperatorFlavorProperties @ operator;
		(*openIndices = operator[[2]];*)
		openIndices = Last@ operator;
		(* if this has indices, determine the conversion to the IBPSimplify database and its inverse *)
		If[Length @ openIndices > 0,
			With[{operatorPatternIndices = Last @ operator},
				permConversion = FindPermutationOrder[operatorPatternIndices, openIndices];
				permConversionInv = FindPermutationOrder[ openIndices, operatorPatternIndices ];
			];
		];

		(* is the operator hermitian? *)
		If[ superleading || internal || !hermitianTerm ,
			(* for the automated routes or if this is not overridden by HermitianTerm -> True, just determine it from the operator *)
			If[Length @ openIndices > 0 && MatchQ[opFlavProp[ConjugateIndexExchange], _List],
				(* take from IBPSimplify database and convert to the ordering found here *)
				If[opFlavProp[ConjugateIndexExchange]===False,
					hermite = False
					,
					pOrder = permConversionInv[[ opFlavProp[ConjugateIndexExchange] ]][[ permConversion ]];
					(* if the list is empty, this operator is not hermitian *)
					If[Length @ pOrder > 0,
						hermite = pOrder
					,
						hermite = False
					]
				]
			,
				(* check hermiticity if coupling/operator do not contain open indices *)
				hermite = MHermitianQ[derivativePhase AtomicToNormalForm@ operator]
			]
		,
			(* called by IntroduceEffectiveCouplings with HermitianTerm -> True, need to determine symmetries *)
			If[Length @ openIndices > 0 && MatchQ[opFlavProp[ConjugateIndexExchange], _List],
				(* take from IBPSimplify database and convert to the ordering found here *)
				pOrder = permConversionInv[[ opFlavProp[ConjugateIndexExchange] ]][[ permConversion ]];

				(* if the list is empty, this operator is simply hermitian without further constraints *)
				If[Length @ pOrder > 0,
					hermite = pOrder
				,
					hermite=True
				]
			,
				(* actually no open indices so just hermitian right away *)
				hermite=True
			];
		];

		symmetries = Association[Range[Length@index]->1];
		If[Length @ openIndices > 0,
			Module[{symrules1,symrules2, ineqp, symrulesrhs, replacedexp},
				(* we begin by finding the symmetries that we know from the operator *)
				symrules1 = (permConversionInv[[ #[[1]] ]][[ permConversion ]] -> #[[2]])&/@(opFlavProp[Symmetries]);
				symrulesrhs = (#[[1]])&/@symrules1;
				(* now we take the inequivalent permutations, delete the ones that are already in the symmetries list *)
				ineqp = Select[ ( permConversionInv[[ # ]][[ permConversion ]] )&/@ opFlavProp[InequivalentPermutations], !MemberQ[symrulesrhs,#]& ];
				(* iterate over all of them and check if they are actual symmetries *)
				symrules2 = Table[
					replacedexp = exp/.((#[[1]]->#[[2]])&/@DeleteCases[Transpose[{openIndices , openIndices[[cand]]}], {here_, here_} ]);
					If[RelabelIndices @ Contract[(replacedexp - exp)((*RelabelIndices@*)AtomicToNormalForm@operator)] === 0,cand -> 1,
						If[RelabelIndices @ Contract[(replacedexp + exp)((*RelabelIndices@*)AtomicToNormalForm@operator)] === 0, cand -> -1, Nothing]],
				{cand, ineqp}];
				symmetries = Association@@Join[symrules1,symrules2]
			];
		];


		(* indices for the replacement rules *)
		indCounter1 = 1;
		inds = Index[Symbol["i" <> ToString[indCounter1++]],#]&/@ (Last/@index);
		pinds = IndexToPattern[inds];

		(* setting up the rules *)
		ruleRHS = exp /. Thread[index->inds];
		
		(* check if there is a tree-level contribution in exp and try to mimic index ordering from there *)
		
		(* if there is a renormalizable tree-level term with a single coupling, take its indices *)
		If[OperatorDimension@AtomicToNormalForm[operator]<=4,
			Module[{treeCouplingIndices,auxTerm},
				(* get tree-level terms with a single coupling *)
				auxTerm=DeleteCases[
					TermsToList@BetterExpand[ruleRHS/.hbar->0], (* w/o TermsToList the level spec does not properly work *)
					x_/;Count[x,Coupling[_,{__},_],All]=!=1,
					1
				];
				auxTerm=Plus@@auxTerm;
				(* get indices of that coupling *)
				treeCouplingIndices=FirstCase[
					auxTerm,
					Coupling[_,indAUX:{__},_]:>indAUX,
					{},
					All
				];
				(* check if its indices match the ones of the effective coupling *)
				If[MatchQ[inds,{OrderlessPatternSequence@@treeCouplingIndices}],
					(* determine the pre-factor for that tree-level coupling *)
					preFactor = auxTerm /. _Coupling->1;
					(* catch case of no pre-factors, i.e., only coupling *)
					If[preFactor===0,preFactor=1];
					(* mimick this prefactor for effective coupling and flip its indices according to the index ordering of the tree-coupling *)
					ruleRHS = ruleRHS/preFactor/.Thread[treeCouplingIndices->inds];
					(* do inverse for the output coupling *)
					outCoupling = preFactor * outCoupling /. Thread[index->inds] /. Thread[inds->treeCouplingIndices] /. Thread[inds->index];
				];
			];
		];
		
		(* we need to temporarly rename the open diagonal indices, otherwise they get relabeled below *)
		ruleRHS = ruleRHS /. (dInd:(Alternatives@@(inds[[Count[isdiagonal,False]+1;;]])):> (dInd/.Index->IndexNotRelabeled));
		ruleLHS = Coupling[temp, pinds, power];

		preRHS = With[{RHS = RelabelIndices[ruleRHS, Unique->True] derivativePhase^-1}, ruleLHS :> (RelabelIndices[RHS,Unique -> True]/.IndexNotRelabeled->Index)];

		If[!ignoreduplicates,
			(* check for any couplings with the identical properties to the one we are about to define *)
			tbdassociation = CouplingAssociationEntry[power, indexType, hermite, symmetries, isdiagonal];
			couplingsSameProperty = Cases[{Keys@$CouplingAssociation,Values@$CouplingAssociation}\[Transpose],{val_, tbdassociation} :> val]
		];

		If[superleading == True ||(superleading == False && internal == False),
			(* this branch is for couplings visible to the user *)

			If[!ignoreduplicates,
				(* check we do not define duplicate couplings *)
				candidates =  Cases[$TempCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ]];
				candidateLabels = Cases[$TempCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ] :> f ];
				preExisting = (SameRuleQ[preRHS,#]& /@ candidates)
			];

			If[Or@@preExisting && !ignoreduplicates,
				(* the coupling has already be defined so figure out which one it is and use it instead *)
				outCoupling = derivativePhase Coupling[candidateLabels[[FirstPosition[preExisting,True][[1]]]], index, power]
			,
				(* does not exist so define it *)
				DefineCoupling[temp, Indices -> indexType, EFTOrder -> power, SelfConjugate -> hermite, Symmetries -> SymmetryOverride[symmetries], DiagonalCoupling -> isdiagonal];
				AppendTo[$TempCouplings,temp];
				AppendTo[$TempCouplingRules, preRHS];
				(* define the NiceForm of this new coupling *)
				AppendTo[LabelsNiceForm[Coupling],temp->prettylabel];
				(*Format[Coupling[temp, in_,_], NiceForm] := UpDownIndices[prettylabel,in];
				If[hermite===False, Format[Bar @ Coupling[temp, in_,_], NiceForm] := UpDownIndices[OverBar @ prettylabel, Bar /@ in]]*)
			]
			
		,
			(* this branch is for couplings invisible to the user *)

			If[!ignoreduplicates,
				(* check we do not define duplicate couplings *)
				candidates = Cases[$IntCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] ]];
				candidateLabels = Cases[$IntCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ] :> f ];
				preExisting = (SameRuleQ[preRHS,#]&/@ candidates)
			];

			If[Or@@preExisting && !ignoreduplicates,
				(* the coupling has already be defined so figure out which one it is and use it instead *)
				outCoupling = derivativePhase Coupling[candidateLabels[[FirstPosition[preExisting,True][[1]]]], index, power]
			,
				(* does not exist so define it *)
				DefineCoupling[temp, Indices -> indexType, EFTOrder -> power, SelfConjugate -> hermite, DiagonalCoupling -> isdiagonal];
				AppendTo[$InternalCouplings,temp];
				AppendTo[$IntCouplingRules, preRHS];
				(* define the NiceForm of this new coupling *)
				AppendTo[LabelsNiceForm[Coupling],temp->prettylabel];
				(*Format[Coupling[temp, in_,_], NiceForm] := UpDownIndices[prettylabel,in];
				If[hermite===False, Format[Bar @ Coupling[temp, in_,_], NiceForm] := UpDownIndices[OverBar @ prettylabel, Bar /@ in]]*)
			]

		];
		Return[outCoupling]
	]


(* ::Subsubsection::Closed:: *)
(*ToEffectiveCoupling (for gauge kinetic terms)*)


ToEffectiveCoupling::invalidGaugeKineticTerm = "The operator `1` is not a gauge-kinetic (or kinetic mixing) term.";
ToEffectiveCoupling::noSuperleading = "Gauge-kinetic terms cannot be defined as superleading.";


ToEffectiveCoupling[exp_, operator:AtomicOp[{{},4},_,{}], OptionsPattern[] ]:=Module[
	{
		hermitianTerm, superleading, internal, ignoreduplicates,
		indepfields, kineticMixing = False, 
		ceff, power, outCoupling,
		couplingstring = "C", temp, prettylabel,
		ruleLHS, predefinedRule, tbdassociation, couplingsSameProperty, candidates, candidateLabels, preExisting
	},

	(* get options *)
	hermitianTerm = TrueQ @ OptionValue @ HermitianTerm;
	superleading = TrueQ @ OptionValue @ Superleading;
	internal = TrueQ @ OptionValue @ Internal;
	ignoreduplicates = TrueQ @ OptionValue @ OverrideDuplicateCouplingCheck;
	
	
	(* TEMP: ignore loop function replacement *)
	If[internal, Print["[ToEffectiveCoupling] Dummy coefficient functionality for gauge-kinetic terms inactive!"]; Return@ exp ];
	
	(* throw various warnings for blind options *)
	If[superleading, Message[ToEffectiveCoupling::noSuperleading]];
	
	(* kinetic mixing terms get a separate treatment *)
	If[(indepfields = Length@DeleteDuplicates[First/@Cases[RelabelIndices@AtomicToNormalForm@ operator,_FieldStrength,Infinity]])==2, 
		kineticMixing= True];

	(*  if this is already in standard form, we do not change the term here *)
	If[
		(!kineticMixing && (MatchQ[exp, -1/(4*_Coupling^2)]|| MatchQ[exp, -1/(4 HoldPattern[Bar[_Coupling]]^2) ]) && !ignoreduplicates)
		||
		(kineticMixing&&!ignoreduplicates&&(MatchQ[exp, (c_:1) x_Coupling /;NumericQ[c]]|| MatchQ[exp, (c_:1) HoldPattern[Bar[x_Coupling]] /;NumericQ[c]]))
	, 
		Return[exp]
	];
	
	If[kineticMixing,
		(* kinetic mixing just takes the prefactor and calls it a coupling *)
		ceff = exp;
	,
		(* match the expression to the form -1/4g^2 F^2 (or -1/2g^2 F.G) *)
		(* Simplify the Sqrt[c^2] from the series of the gauge couplings *)
		ceff = Assuming[0 < hbar < 1, Normal@ Series[Sqrt[(-1/4)/exp],{hbar,0,1}]]/. 
			Power[Power[c:Coupling[_, {}, 0], 2], Rational[n_Integer, 2]]:> c^n
	];

	If[Head @ OptionValue[ EffectiveCouplingSymbol ] === String, couplingstring = OptionValue[ EffectiveCouplingSymbol ] ];
	{temp, prettylabel} = CouplingNameFromOperator[operator, couplingstring];
	
	power = GetMinOrder@ exp;
	
	outCoupling = Coupling[temp, {}, power];
	ruleLHS = outCoupling;
	predefinedRule = With[{RHS=RelabelIndices[ceff,Unique->True]}, ruleLHS :> RHS];
	
	If[!ignoreduplicates,
			(* check for any couplings with the identical properties to the one we are about to define *)
			tbdassociation = CouplingAssociationEntry[power, {}, True, <|{}->1|>, {}];
			couplingsSameProperty = Cases[{Keys@$CouplingAssociation,Values@$CouplingAssociation}\[Transpose],{val_, tbdassociation} :> val]
		];
	
	If[internal,
		(* this branch is for couplings invisible to the user - stale branch! *)
		0
	,
		(* this branch is for couplings visible to the user *)
		If[!ignoreduplicates,
				candidates =  Cases[$TempCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ]];
				candidateLabels =  Cases[$TempCouplingRules, HoldPattern[ Coupling[f_/;MemberQ[couplingsSameProperty,f],__] :> _ ] :> f ];
				preExisting =  (SameRuleQ[predefinedRule,#]& /@ candidates)
		];
		
		If[Or @@ preExisting && !ignoreduplicates,
			(* this coupling exists already - reuse it *)
			outCoupling = Coupling[candidateLabels[[FirstPosition[preExisting,True][[1]]]], {}, power]
		,
			(* coupling is new, or check is disabled *)
			DefineCoupling[temp, EFTOrder -> power, SelfConjugate -> True];
			AppendTo[$TempCouplings,temp];
			AppendTo[$TempCouplingRules, predefinedRule];
			AppendTo[LabelsNiceForm[Coupling],temp -> prettylabel]
		]
	];
	
	If[kineticMixing,
		Return@ outCoupling
	,
		Return[ - 1/(4 outCoupling^2)]
	]
]


(* ::Subsubsection::Closed:: *)
(*More helpers*)


Options @ SubstituteCoefficients = {EffectiveCouplingSymbol -> "c"};

SubstituteCoefficients[L_, OptionsPattern[]]:=Module[{op3,tempCplPre,newCpls,LOut,res},
	(* list of all operators that appear at superleading power *)
	op3 =  DeleteDuplicates [ SortBy[DeleteDuplicates@Cases[SeriesEFT[L, EFTOrder -> 3],_Operator,Infinity],LeafCount], SameOperatorQ];

	(* replace the coefficients of all operators appearing in the above list by a newly defined effective coupling *)
	(*tempCplPre = $TempCouplings;*)
	tempCplPre = DeleteCases[$TempCouplings, x_/;FreeQ[L,x,All]];
	
	(*
	LOut=Collect[L, _Operator]/.{
		x_Operator*c_ /;(MemberQ[op3, x] && !BarredOpQ[x]) :> Operator@PlusHc@NormalForm[x ToEffectiveCoupling[c, x, EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol]]]-Operator@Bar@NormalForm[x c],
		x_Operator * c_ /;(MemberQ[op3, x] && BarredOpQ[x]) :> x ToEffectiveCoupling[c, x, EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol]]};
	*)
	(* ensure that flavor Deltas are factored out of the operators for proper relabeling of diagonal flavor indices *)
	LOut=CollectOperators[L, NormalForm->False]/.{
		x_Operator*c_ /;(MemberQ[op3, x] && !BarredOpQ[x]) :> Operator@PlusHc@OperatorToNormalForm[x ToEffectiveCoupling[c, x, EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol]],Unique->True]-Operator@Bar@OperatorToNormalForm[x c,Unique->True],
		x_Operator * c_ /;(MemberQ[op3, x] && BarredOpQ[x]) :> x ToEffectiveCoupling[c, x, EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol]]};

	(*newCpls=TreeReplacement /@ DeleteCases[$TempCouplings, x_/;MemberQ[tempCplPre,x]];*)
	newCpls=TreeReplacement /@ DeleteCases[$TempCouplings, x_/;(FreeQ[LOut,x,All] || MemberQ[tempCplPre,x])];

	(* this can be very slow for very comlex models *)
	(*LOut /. newCpls*)

	(* The Collect below can significantly improve performance for more complicated models *)
	With[{aux= newCpls[[;;,1]]},
		res= Collect[L, aux] /. newCpls
	];

	res
]


(* ::Text:: *)
(*Old version (issues with flavor indices & delta functions)*)


(*
(* derive replacement rules to eliminate tree-level objects from matching coefficients once the above replacements have been made *)
TreeReplacement[c_] := Module[{finalRep,inds,pinds,fullRHS,treeRHS,deltaRHS, \[Delta]RHS ,variable,rawRule,finalLHS,TP},
	inds=Table[Unique["i"],{incr, Length@GetCouplings[c][Indices]}];
	pinds=Pattern[#, _]&/@inds;
	fullRHS = c[Sequence@@inds]/.$TempCouplingRules;
	treeRHS = fullRHS /. hbar -> 0;
	Echo[Format[(c[Sequence@@inds]//RelabelIndices)==(fullRHS//RelabelIndices),NiceForm]];
	If[treeRHS===0, Return[Nothing]];
	deltaRHS = fullRHS - treeRHS;
	variable = First@Cases[treeRHS, _Coupling, All];
	rawRule = (Last@Quiet@Solve[c[Sequence@@inds] == treeRHS + \[Delta]RHS, variable]/.\[Delta]RHS -> deltaRHS)[[1]];
	(* hacky way of generating patterns *)
	finalLHS = (rawRule[[1]] /. Index[f_,g_] :> Index[TP[f,_],g])/.TP -> Pattern;
	EchoFunction[Format[#,NiceForm]&][With[{theRHS = rawRule[[2]]}, finalLHS :> RelabelIndices[theRHS,Unique->True]]]
]
*)


TreeReplacement::notree="No tree-level coupling found that could be solved for."
TreeReplacement::nosolution="Could not solve for the coupling `1` at tree level. Possibly, this coupling appears with different flavor index contractions at tree-level: `2`"


(* derive replacement rules to eliminate tree-level objects from matching coefficients once the above replacements have been made *)
TreeReplacement[c_] := Module[{finalRep,inds,pinds,fullRHS,treeRHS,deltaRHS,\[Delta]RHS,variable,var,rawRule,finalLHS,TP,deltas,MyRelabelIndices,lhs,rhs,pow,dropGaugeCouplings},
	inds    = Table[Unique["i"],{incr, Length@GetCouplings[c][Indices]}];
	pinds   = Pattern[#, _]&/@inds;
	fullRHS = c[Sequence@@inds]/.$TempCouplingRules;
	treeRHS = fullRHS /. hbar -> 0;
	If[treeRHS===0, Return[Nothing]];
	deltaRHS = fullRHS - treeRHS;

	(*Echo[Format[(c[Sequence@@inds]//RelabelIndices)==(fullRHS//RelabelIndices),NiceForm]];*)

	(* find all flavor Deltas in treeRHS *)
	treeRHS = ContractDelta[treeRHS];
	deltas  = Times@@Cases[treeRHS,_Delta];
	(* multiply everything with these Deltas, such that treeRHS no longer contains Deltas *)
	treeRHS = ContractDelta[treeRHS * deltas];
	lhs     = ContractDelta[c[Sequence@@inds] * deltas];
	(* set set placeholder variable \[Delta]RHS as complex *)
	NonBarableQ[\[Delta]RHS]=False;
	(* construct RHS *)
	rhs     = treeRHS + \[Delta]RHS * deltas;

	(* find variable to solve for *)
	(* 
	Take all couplings on the RHS, except for:
		- terms with inverse powers of couplings/masses
	*)
	variable = Cases[treeRHS/.Times[___,Power[_Coupling,_?Negative]]->0/.Power[_Coupling,_?Negative]->0, _Coupling, All];
	
	(* filter for variables with same property as the effective couplings *)
	variable = First[
		Select[variable, MatchQ[KeyDrop[GetCouplings[First@#],EFTOrder],KeyDrop[GetCouplings[c],EFTOrder]]&]
		,Message[TreeReplacement::notree]; Abort[]
	];
	
	(* check that this variable is not appearing with a different flavor index compination *)
	If[!FreeQ[Hold[rhs]/.variable->0, First[variable], All], 
		Message[TreeReplacement::nosolution, Format[variable,NiceForm], Format[treeRHS,NiceForm]]; Abort[]
	];
	
	(* catch case where powers of couplings appear *)
	If[FreeQ[treeRHS/.{Power[variable,_]|Bar@Power[variable,_]->0}, variable, All] && !FreeQ[treeRHS,Power[variable,_]|Bar@Power[variable,_],All],
		variable = First@SortBy[Cases[treeRHS,Power[variable,_],All],Last];
	];
	
	(* solve for variable *)
	NonBarableQ[var]=False; (* make auxiliary variable var complex *)
	With[{lhsLOC=lhs, rhsLOC = rhs /. variable->var},
		rawRule = First[Last@Solve[lhsLOC == rhsLOC, var, InverseFunctions->True]/.HoldPattern[InverseFunction[Bar,1,1]]->Bar];
	];
	rawRule=rawRule/.var->variable;
	
	(* prepare RHS *)
	rawRule[[2]] = MyRelabelIndices[ContractDelta@rawRule[[2]],Unique->True] /. \[Delta]RHS -> deltaRHS;
	
	(* remove powers from rule if it was solved for Power@Coupling *)
	If[Head[variable]===Power,
		pow = Last@variable;
		rawRule[[1]] = Power[rawRule[[1]],1/pow]//.Power[Power[coup_Coupling,a_],b_]:>Power[coup,a*b];
		rawRule[[2]] = Power[rawRule[[2]],1/pow];
	];
	
	(* hacky way of generating patterns *)
	finalLHS = IndexToPattern[rawRule[[1]]];

	With[{theRHS = rawRule[[2]]}, (finalLHS :> theRHS)/.MyRelabelIndices->RelabelIndices]
]


(* ::Subsubsection::Closed:: *)
(*User output for the coupling*)


(* ::Text:: *)
(*Function for printing or outputting the rules of relevant replacement couplings for a given Lagrangian *)


Options @ PrintEffectiveCouplings= { Rules -> False}


PrintEffectiveCouplings[L_, OptionsPattern[]]? OptionsCheck:= Module[{rules={},couplings,inds,$ind},
	couplings = Cases[DeleteDuplicates[Cases[L, _Coupling, Infinity]/.Coupling[a_,__]:> a], x_/;MemberQ[ $TempCouplings,x]];
	If[Length@couplings === 0, Print["The supplied Lagrangian does not contain effective couplings."];Return[]];

	(* return replacement rules if Rules-> True*)
	If[OptionValue@ Rules,
		Return@ DeleteCases[$TempCouplingRules, _? (FreeQ[First@ #, Alternatives@@ couplings]&) ];
	];
	Do[
		inds=With[{l=Alphabet[], n=Length@(GetCouplings[term][Indices])},
			(Symbol/@ (Flatten@Table[StringJoin[#,If[i=!=0,ToString[i],""]]&/@l,{i,0,(Ceiling[n/Length[l]]-1)}])[[;;n]])
		];
		Print[NiceForm[term[Sequence@@inds]], " = ", NiceForm[ Collect[RelabelIndices[(term[Sequence@@inds]/.$TempCouplingRules)/.Index[i:Alternatives@@inds,rep_]:>$ind[i,rep]]/.$ind->Index,hbar,Simplify]]]
	,
		{term,couplings}];
]


(* replace all effective couplings in the given expression *)
ReplaceEffectiveCouplings[inputL_]:=Module[{L=HcExpand[inputL],LagTree,LagLoop,TempCouplingRules,TempCouplingRulesTree},
	(* split Lagrangian and replavement rules into tree- and loop-level part for better performance *)
	LagTree               = L/.hbar->0;
	LagLoop               = L-LagTree;
	TempCouplingRules     = $TempCouplingRules;
	TempCouplingRulesTree = TempCouplingRules/.Times[hbar,___]->0;
	
	(* perform replacements and expand result to drop remaining hbar^2 conributions *)
	Normal@Series[FastExpand[(LagTree//.TempCouplingRules)],{hbar,0,1}]+FastExpand[(LagLoop//.TempCouplingRulesTree)]
]


(* replace only the effective couplings with the labels given in the second argument *)
ReplaceEffectiveCouplings[inputL_, labels_List]:=Module[{L=HcExpand[inputL],LagTree,LagLoop,TempCouplingRules,TempCouplingRulesTree},
	(* split Lagrangian and replavement rules into tree- and loop-level part for better performance *)
	LagTree               = L/.hbar->0;
	LagLoop               = L-LagTree;
	TempCouplingRules     = Cases[$TempCouplingRules, HoldPattern[Coupling[x_/;MemberQ[labels,x],__]:>_ ]];
	TempCouplingRulesTree = TempCouplingRules/.Times[hbar,___]->0;
	
	(* perform replacements and expand result to drop remaining hbar^2 conributions *)
	Normal@Series[FastExpand[(LagTree//.TempCouplingRules)],{hbar,0,1}]+FastExpand[(LagLoop//.TempCouplingRulesTree)]
]


(* ::Text:: *)
(*Old version w/o splitting tree and loop parts for the replacements*)


(*ReplaceEffectiveCouplings[L_]:=ExpandReplacedEffectiveCouplings[
	L //. $TempCouplingRules
]*)


(*(* replace only the effective couplings with the labels given in the second argument *)
ReplaceEffectiveCouplings[L_, labels_List]:=ExpandReplacedEffectiveCouplings[
	L //. Cases[$TempCouplingRules, HoldPattern[Coupling[x_/;MemberQ[labels,x],__]:>_ ]]
]*)


(*(* significant simplififcation of one-loop reults by Collecting hbar *)
(* FastExpand is much faster than BetterExpand for products of large sums *)
ExpandReplacedEffectiveCouplings[arg_]:=FastExpand@Collect[arg,hbar]*)


(* ::Text:: *)
(*Very old version w/o expanding the result*)


(*
ReplaceEffectiveCouplings[L_]:=L /. $TempCouplingRules
*)


(*
(* replace only the effective couplings with the labels given in the second argument *)
ReplaceEffectiveCouplings[L_, labels_List]:=L //. Cases[$TempCouplingRules, HoldPattern[Coupling[x_/;MemberQ[labels,x],__]:>_ ]]
*)


(* ::Subsection:: *)
(*Replacing coefficients with effective couplings*)


(* ::Subsubsection::Closed:: *)
(*Setting up effective couplings for loop functions*)


IntroduceDummyCoefficients[L_]:=Module[{tL, L0, tLHcTerms, tLHTerms},
	(* extract all one-loop terms of our Lagrangian *)
	L0 = GreensSimplify @ L;
	(*tL = List@@(Nothing + (HcSimplify@Coefficient[L0, hbar]));*)
	tL = TermsToList[HcSimplify@Coefficient[L0, hbar]];

	(*
	tLHcTerms = Collect[Operator[(Total@Cases[tL, _HcTerms])/.HcTerms -> Identity],_Operator];
	tLHTerms = Collect[Operator@Total@(tL /. _HcTerms -> 0),_Operator];
	*)
	(* ensure that flavor Deltas are factored out of the operators for proper relabeling of diagonal flavor indices *)
	tLHcTerms = CollectOperators[(Total@Cases[tL, _HcTerms])/.HcTerms -> Identity, NormalForm->False];
	tLHTerms =  CollectOperators[ Total@(tL /. _HcTerms -> 0), NormalForm->False];

	tLHcTerms = tLHcTerms /. x_Operator * c_ :> x ToEffectiveCoupling[c, x, Superleading -> False];
	tLHTerms = tLHTerms /. x_Operator * c_ :> x ToEffectiveCoupling[c, x, Superleading -> False];
	
	RelabelIndices[hbar(Operator@PlusHc[OperatorToNormalForm[tLHcTerms,Unique->True]]+1/2 Operator@PlusHc[OperatorToNormalForm[tLHTerms,Unique->True]])+Operator@(L0/.hbar->0)]
]


(* ::Subsubsection::Closed:: *)
(*Effective couplings for the final Lagrangian*)


Options @ IntroduceEffectiveCouplings = {EffectiveCouplingSymbol -> "C", OverrideDuplicateCouplingCheck -> False};


IntroduceEffectiveCouplings[L_ , OptionsPattern[]] := Module[{L0, LHcTerms, LHTerms, noCouplingTerms},
	(* subtract gauge-kinetic terms, those get a special treatment *)
	L0 = TermsToList@ L;
	
	(* determine which terms have no couplings at all and subtract these from the Lagrangian *)
	noCouplingTerms = Total @ Select[L0, (FreeQ[#,hbar] && FreeQ[#,Coupling]&)];
	
	L0 = HcSimplifyInternal[ Total @ Select[ L0, !(FreeQ[#,hbar] && FreeQ[#,Coupling])&], InternalOpRepresentation-> True];
	L0 = If[Head@L0===Plus, List@@L0,{L0}];
	
	LHcTerms = Total @ Cases[L0, _HcTerms] /. HcTerms -> Identity;
	LHTerms  = Total[ L0 /. _HcTerms -> 0 ];

	LHcTerms = LHcTerms /. (x:(_AtomicOp|_CompOp)) * c_ :> x ToEffectiveCoupling[c, x, Superleading -> False, Internal -> False, HermitianTerm -> False , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ];
	LHTerms  = LHTerms  /. (x:(_AtomicOp|_CompOp)) * c_ :> x ToEffectiveCoupling[c, x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ];

	HcTerms[ RelabelIndices@AtomicToNormalForm @ LHcTerms ] + RelabelIndices@AtomicToNormalForm @ LHTerms + noCouplingTerms
]


(* ::Subsubsection::Closed:: *)
(*Effective masses*)


Options @ IntroduceEffectiveMasses = {EffectiveCouplingSymbol -> "\[ScriptCapitalM]", OverrideDuplicateCouplingCheck -> False};


IntroduceEffectiveMasses[L_ , OptionsPattern[]] := Module[{L0, LHcTerms, LHTerms, noCouplingTerms, Lmass, field, prop, newMasses},
	L0 = Expand @ L;

	(* Extrac mass terms *)
	Lmass = IsolateMassTerms[L0, Heavy -> True];
	L0 = L0 - Lmass;

	Lmass = List @@ (Nothing + HcSimplify @ Lmass);
	Lmass = Contract/@Lmass;

	LHcTerms = Collect[Operator[(Total @ Cases[Lmass, _HcTerms])/.HcTerms -> Identity], _Operator];
	LHTerms  = Collect[Operator @ Total @ (Lmass /. _HcTerms -> 0), _Operator];

	Switch[Head[LHTerms],
		Plus, LHTerms = List@@ LHTerms,
		Times, LHTerms = {LHTerms}
	];

	Switch[Head[LHcTerms],
		Plus, LHcTerms = List@@ LHcTerms,
		Times, LHcTerms = {LHcTerms}
	];

	newMasses = Last@ Reap[
	LHcTerms = Plus@@ Table[
		(*find fields*)
		field = DeleteDuplicates[Cases[term, Field[lab_,___]:>lab, All]];
		If[Length@field > 1,
			term (*ignore mass mixing terms*)
			,
			field = First@ field;
			prop = GetFields[field];
			(* only chiral fermions have non-hermitian mass terms *)
			If[prop[Type]===Fermion && (!prop[SelfConjugate]) && MatchQ[prop[Chiral], LeftHanded|RightHanded],
				If[MatchQ[term,-(1/2)*_Coupling*_Operator],
					term,
					term /. x_Operator * c_ :> -(1/2)*x Sow@ToEffectiveCoupling[-2*c, x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ]
				],
				term
			]
		]
		,
		{term, LHcTerms}
	];

	LHTerms = Plus@@ Table[
		(*find fields*)
		field = DeleteDuplicates[Cases[term, Field[lab_,___]:>lab, All]];
		If[Length@field > 1,
			term (*ignore mass mixing terms*)
			,
			field = First@ field;
			prop = GetFields[field];
			Switch[prop[Type],
				Scalar,
					If[prop[SelfConjugate],
						If[MatchQ[term,-(1/2)*(_Coupling)^2*_Operator],
							term,
							term /. x_Operator * c_ :> -(1/2)*x Sow@(ToEffectiveCoupling[Sqrt[-2*c], x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ]^2)
						], (*real*)
						If[MatchQ[term,-(_Coupling)^2*_Operator],
							term,
							term /. x_Operator * c_ :> -x Sow@(ToEffectiveCoupling[Sqrt[-c], x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ]^2)
						] (*complex*)
					],
				Vector,
					If[prop[SelfConjugate],
						If[MatchQ[term,1/2*(_Coupling)^2*_Operator],
							term,
							term /. x_Operator * c_ :> 1/2*x Sow@(ToEffectiveCoupling[Sqrt[2*c], x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ]^2)
						], (*real*)
						If[MatchQ[term,(_Coupling)^2*_Operator],
							term,
							term /. x_Operator * c_ :> x Sow@(ToEffectiveCoupling[Sqrt[c], x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ]^2)
						] (*complex*)
					],
				Fermion,
					If[prop[SelfConjugate],
						If[MatchQ[term,-(1/2)*_Coupling*_Operator],
							term,
							term /. x_Operator * c_ :> -(1/2)*x Sow@ToEffectiveCoupling[-2*c, x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ]
						], (*Majorana*)
						If[MatchQ[prop[Chiral],False],
							If[MatchQ[term,-_Coupling*_Operator],
								term,
								term /. x_Operator * c_ :> -x Sow@ToEffectiveCoupling[-c, x, Superleading -> False, Internal -> False, HermitianTerm -> True , EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol] , OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ]
							], (*Dirac*)
							term
						]
					],
				_ , term
			]
		]
		,
		{term, LHTerms}
	];
	];

	If[Length@newMasses>0,
		Print["The following effective mass terms have been defined to simpify the Lagrangian:"];
		PrintEffectiveCouplings[newMasses];
	];

	HcTerms[ RelabelIndices @ CollectOperators @ OperatorToNormalForm[ LHcTerms,Unique->True]] +  RelabelIndices @ CollectOperators[ ( PlusHc @ OperatorToNormalForm[ LHTerms,Unique->True]) / 2 ] + L0
]


(* ::Subsection:: *)
(*Shifting corrections to renormalizable operators*)


(* ::Subsubsection::Closed:: *)
(*Main module*)


(* ::Text:: *)
(*Function to redefine renormalizable couplings of the EFT, thereby absorbing threshold corrections into the UV parameter. E.g.*)
(*	y_IR = y_UV + \[Delta]y -> y_UV = y_IR - \[Delta]y allows us to write all matching corrections in terms of y_IR.*)


Options[ShiftRenCouplings] = {EFTOrder->All};


ShiftRenCouplings[Lag_,OptionsPattern[]]:=Module[
		{OldCp,NewCp,LEFTRen,RenOpList,LEFTCorrections,CouplingCorrections,RenCouplings,RepRules,MaxOrder,
			LagTree,LagLoop,RepRulesTree,dropGaugeCouplings,Lagrangian},
	(* 
	Take all renormalizable (d<=4) terms and drop the ones containing:
		- loop corrections: 
			only tree-level couplings get shifted;
	*)
	(*dropGaugeCouplings = (Power[Coupling[#,___],n_/;n>1]->0)&/@Table[GetGaugeGroups[group][Coupling],{group,Keys@GetGaugeGroups[]}];*)
	(* the rule below ensures that all positive powers of gauge couplings are ignored for the redefinition *)
	dropGaugeCouplings = (Except[Power[Coupling[#,___],n_/;n<0], _[___,Coupling[#,___],___]]->0)&/@Table[GetGaugeGroups[group][Coupling],{group,Keys@GetGaugeGroups[]}];
	LEFTRen            = CollectOperators[SeriesEFT[Lag/.hbar->0, EFTOrder->4],NormalForm->False]/.dropGaugeCouplings;
	(*TruncateOperator admits also higher-EFT order  corrections to marginal operators*)
	LEFTCorrections    = CollectOperators[Operator[Lag]/.op_Operator:> TruncateOperator[op,4],NormalForm->False];
	
	(* Renormalizable operators not present in LEFTRen are set to zero *)
	RenOpList       = Cases[LEFTRen, _Operator,All];
	LEFTCorrections = CollectOperators[LEFTCorrections /. x_Operator?(!MemberQ[RenOpList,#]&) :> 0];
	
	(* Define effective couplings that absorb the corrections after the shift *)
	OldCp           = $TempCouplings;
	LEFTCorrections = IntroduceEffectiveCouplings[LEFTCorrections(*,OverrideDuplicateCouplingCheck->True*),EffectiveCouplingSymbol->"c"];
	NewCp           = Complement[$TempCouplings,OldCp];	
	
	(* Get rules to substitute original couplings for the new effective couplings [at loop and tree level] *)
	RepRules     = TreeReplacement/@NewCp;
	RepRulesTree = RepRules //.Plus[Times[hbar,___],rest___]:>Plus[rest] //.Times[hbar,___]->0 //.hbar->0;
	
	(* Separate Lagrangian into tree- and loop-level parts *)
	With[{patterns = Alternatives@@Cases[RepRules[[;;,1]],Coupling[l_,___]:>Coupling[l,___],All]},
	LagTree = Collect[Lag/.hbar->0, patterns, BetterExpand];
	LagLoop = Collect[Lag-LagTree, patterns, BetterExpand];
	];
	
	(* substitute in replacments *)
	LagTree = LagTree/. (*Account for inverse gauge couplings: consider implementing manual expansion if performance issues arise*)
		g:Power[Coupling[_, {}, 0], -2]:> Normal@ Series[g/. RepRules, {hbar, 0, 1}]/. 
		RepRules/.RepRulesTree;
	LagLoop = LagLoop/.RepRulesTree;
	
	(* expand results *)
	LagTree = FastExpand[LagTree];
	LagLoop = FastExpand[LagLoop];
	
	(* Determine desired EFT order for truncation *)
	If[Head@OptionValue[EFTOrder]===Integer, MaxOrder = OptionValue@EFTOrder, MaxOrder=Max[OperatorDimension/@ TermsToList[Lag]]];
	(* EFT series truncation *)
	LagTree = SeriesEFT[LagTree, EFTOrder->MaxOrder];
	LagLoop = SeriesEFT[LagLoop, EFTOrder->MaxOrder];
	
	(* the RHS of loop replacements might again contain the params, thus a 2nd tree replacement is required *)
	CollectOperators[LagTree+LagLoop]
]


(* ::Text:: *)
(*Old version*)


(* ::Subsection:: *)
(*Mapping effective Lagrangians*)


(* ::Text:: *)
(*Determine the values of the coefficients in a target EFT, such that it reproduces the input Lagrangian (possibly in different bases)*)


(* ::Subsubsection::Closed:: *)
(*Auxiliary functions*)


(* ::Text:: *)
(*Function extracting the coefficient equalities from the difference of two Lagrangians *)


CoefficientEqualities@ lagDiff_:= Module[{diff, term},
	diff= HcSimplifyInternal[lagDiff, InternalOpRepresentation-> True]/. HcTerms-> Identity;
	diff= If[Head@ diff === Plus, List@@ diff, {diff}];
	Table[
			(*Check if gauge term*)
			If[FreeQ[term, AtomicOp[{{}, 4}, _, _]],
				term/. _AtomicOp|_CompOp-> 1
			,
				term= term/. _AtomicOp-> 1// Expand;
				term* Times@@ DeleteDuplicates@ Cases[term, Power[g:Coupling[_, {}, 0], -2]:> g^2, All]//Expand
			]
		, {term, diff}]
];


(* ::Text:: *)
(*Function to create the replacement rules*)


OpenIndexToPattern[rule_]:=Module[{lhs=rule[[1]],rhs=rule[[2]],index,inds,pinds,coeff},
	 (* looking for open indices and creating the pattern *)
	index = FindOpenIndices@lhs/.Index[a_,_]:>a;
	inds = Symbol["i" <> ToString[#]]& /@ Range@Length@index;
	pinds = Pattern[#, _]& /@ inds;

	coeff=lhs/(lhs/.Times[pre___,fact_?NumericQ,post___]:>pre post);

	(* setting up the rule *)
	If[Head[lhs]===Bar,
		(Bar[lhs/coeff]/. Thread[index->pinds])-> (RelabelIndices[Bar[rhs/coeff],Unique->True]/. Thread[index->inds])
		,
		(lhs/coeff/. Thread[index->pinds])-> (RelabelIndices[rhs/coeff,Unique->True]/. Thread[index->inds])
	]
]


(* remove rules of the form x[a_,b_] -> x[a,b] from the output *)
TrivialRuleQ[HoldPattern[lhs_ -> Except[_Coupling]]]=False;

TrivialRuleQ[HoldPattern[lhs_ -> rhs_ ]]:=Module[{lhindices,pattern},
(* do the coupling labels fit? *)
If[lhs[[1]]===rhs[[1]],
	(* transform index patterns to indices and check if they are the same as the rhs *)
	lhindices = (lhs[[2]] /.Pattern-> pattern)/.pattern[a_,Blank[]]:>a;
	If[ lhindices === rhs[[2]],True, False]
,
False]
]


(* ::Subsubsection::Closed:: *)
(*Main exported function*)


MapEffectiveCouplings::SolutionNotFound = "No solution was found: either the two Lagragians are inequivalent or Matchete simplification routines are unable to relate them.";
MapEffectiveCouplings::GaugeFieldMismatch = "The gauge fields appearing in the input and target Lagrangian do not match. Input field strengths: `1`. Target field strengths: `2`.";


MapEffectiveCouplings::NormMismatch =
"The normalization of the kinetic term of the gauge field `1` is different in both Lagrangians.
Redefining the gauge field \!\(\*SubscriptBox[\(`1`\), \(\[Mu]\)]\) \[Rule] (1+\[Delta]) \!\(\*SubscriptBox[\(`1`\), \(\[Mu]\)]\) and gauge coupling `2` \[Rule] (1-\[Delta]) `2` in the input Lagrangian to match their definitions to the ones of the target Lagrangian.
This leaves covariant derivatives invariant but transforms field-strength tensors by \!\(\*SubscriptBox[\(`1`\), \(\[Mu]\[Nu]\)]\) \[Rule] (1+\[Delta]) \!\(\*SubscriptBox[\(`1`\), \(\[Mu]\[Nu]\)]\).
The shift was determined as: \[Delta] = `3`.
The relation of the target Lagrangian gauge coupling to the input Lagrangian gauge coulping (\!\(\*SuperscriptBox[\(`2`\), \(target\)]\) \[Rule] (1+\[Delta]) \!\(\*SuperscriptBox[\(`2`\), \(input\)]\)) is added to the matching conditions in the output."


MapEffectiveCouplings::const = "The following constant terms were detected and removed in the Lagrangians:
Input Lagrangian: const. = `1`,
Target Lagrangian: const. = `2`."


Options@ MapEffectiveCouplings = {
		EOMSimplify -> True, 
		ReductionIdentities -> EvanescenceFree,
		ShiftRenCouplings -> False, 
		SortByEFTOrder -> True, 
		AppendEffectiveCouplingsDefs -> False, 
		KeepTrivalReplacements -> True, 
		Verbose -> True
	};


MapEffectiveCouplings[inputLag_,targetLag_,opts:OptionsPattern[]]:=Module[
		{lagInput,lagTarget,constsInput,constsTarget,
		lagInputTree,lagTargetTree,ruleTree},
	
	(* drop constant terms and throw a warning *)
	{constsInput, lagInput}= SeparateOutConstants@ inputLag;
	{constsTarget, lagTarget}= SeparateOutConstants@ targetLag;
	
	If[constsInput =!= 0 || constsTarget =!= 0,
		Echo@StringForm[MapEffectiveCouplings::const, Format[constsInput,NiceForm], Format[constsTarget,NiceForm]];
	];

	(* simplify input & target Lagrangian *)
	If[OptionValue@ EOMSimplify,
		lagTarget = OptionalMonitor[OptionValue@ Verbose, 
			EOMSimplify[lagTarget, ReductionIdentities-> OptionValue@ ReductionIdentities], 
			"Simplifying target Lagrangian \[Ellipsis]"];
		lagInput  = OptionalMonitor[OptionValue@ Verbose, 
			EOMSimplify[lagInput, ReductionIdentities-> OptionValue@ ReductionIdentities], 
			"Simplifying input Lagrangian \[Ellipsis]"];
	];
	
	(* ensure that all flavor deltas have been extracted from operator *)
	(*lagInput = CollectOperatorsWithGenricFlavorStructure[lagInput];*)

	If[OptionValue@ ShiftRenCouplings,
		(* target Lagrangian should not be shifted *)
		lagInput  = OptionalMonitor[OptionValue@ Verbose, 
			ShiftRenCouplings@ lagInput, "Shifting renomarlizable couplings of the input Lagrangian \[Ellipsis]"];
	];
	
	If[OptionValue@ ReductionIdentities === EvanescenceFree,
		lagInputTree= lagInput/. hbar->0;
		lagTargetTree= lagTarget/. hbar->0;
		
		ruleTree= MapEffectiveCouplingsInternal[lagInputTree,lagTargetTree,opts];
		
		lagInput= lagInput - (hbar Coefficient[lagTarget, hbar]/. ruleTree);
		lagTarget= lagTargetTree;
	];
	
	MapEffectiveCouplingsInternal[lagInput,lagTarget,opts]
]


(* ::Subsubsection::Closed:: *)
(*Internal function doing most of the job*)


Options@ MapEffectiveCouplingsInternal = Options@ MapEffectiveCouplings;


MapEffectiveCouplingsInternal[Lag_,TargetLag_,OptionsPattern[]]:=Module[
		{oldcp,newcp,lagInput,LDiff,TargetCouplings,terms,sol,output,orderAssociation,
			mislabeledTerms,indexReplacement, var={}, opt$CollectCoefficients},
		
	(* MapEffectiveCouplings is not compatible with flavor indices that are contracted within operators -> turn of this functionality *)
	opt$CollectCoefficients = Options[CollectCoefficients]; (* possibly add a CheckAbort below to make sure the settings are always reverted *)
	SetOptions[CollectCoefficients, "ContractOverallDelta" -> False];

	(* save all temporary couplings, must be done after ShiftRenCouplings *)
	(*oldcp=$TempCouplings;*)
	oldcp = DeleteCases[$TempCouplings, x_/;FreeQ[Lag,x,All]];

	(* use dummy coefficients for input Lagrangian *)
	lagInput = OptionalMonitor[OptionValue@Verbose,
		IntroduceEffectiveCouplings[Lag,OverrideDuplicateCouplingCheck->True,EffectiveCouplingSymbol->"TEMP"]
	,"Introducing dummy coefficients to simplify the input Lagrangian \[Ellipsis]"];
	
	(* determine newly introduced effective couplings *)
	(*newcp=Complement[$TempCouplings,oldcp];*)
	newcp = Complement[
		DeleteCases[$TempCouplings, x_/;FreeQ[lagInput,x,All]],
		oldcp
	];

	(* compute the difference of the two Lagrangians *)
	LDiff= OptionalMonitor[OptionValue@Verbose, CollectOperators[lagInput-TargetLag], "Computing difference of both Lagrangians \[Ellipsis]"];

	(* determine appropriate variables, and solve for them *)
	OptionalMonitor[OptionValue@Verbose,
		TargetCouplings= DeleteDuplicates@Cases[TargetLag, Coupling[l_,___]:>l,Infinity];
		terms= CoefficientEqualities@LDiff;
		
		(* Here goes the annoying process of selecting the right variables for the system to solve *)
		var=Module[{tmpvar,varSingle={},varMultiple={}},
			Do[
				(* For each term, couplings of the target Lagrangian are stored in tmpvar *)
				tmpvar={};
				terms[[i]]/. coup_Coupling?(MemberQ[TargetCouplings,#[[1]]]&):>(AppendTo[tmpvar, coup];);
				(* If there is a single target Lagrangian coupling in the term, this is cleary a variable, otherwise we need to do extra work *)
				If[Length@tmpvar===1,AppendTo[varSingle,tmpvar[[1]]],AppendTo[varMultiple,tmpvar]]
			,{i,Length@terms}];
			
			(* Here we deal with multiple variable terms *)
			Do[
				(* For those terms with more than one target Lagrangian coupling, I exclude those couplings already appearing in varSingle *)
				tmpvar=Select[varMultiple[[i]],!MemberQ[varSingle,#]&];
				(* If more than one coupling remains after this, I make sure that the one chosen is either new or a different flavor combination of a coupling already in varSingle *)
				If[tmpvar=!={},
					Module[{filterVar=Select[tmpvar,!MemberQ[varSingle/.Coupling[name_,___]:>name,#/.Coupling[name_,___]:>name]&]},If[filterVar=!={},tmpvar=filterVar]];
					AppendTo[varSingle,First@tmpvar];
				]
			,{i,Length@varMultiple}];
			varSingle
		];
			
		(* get all couplings with mislabeled flavor indices *)
		mislabeledTerms= Cases[
			DeleteDuplicates@Cases[
				terms,
				Coupling[Alternatives@@(var[[;;,1]]),___],
				All
			],
			x_/;!MatchQ[x,Alternatives@@var]
		];
		
		Do[
			(* determine required index replacement *)
			indexReplacement= MapThread[
				(#1->#2)&,
				{coeff[[2]], FirstCase[var,Coupling[First@coeff,inds_,___]:>inds]}
			];
			
			(* apply replacements to the required terms *)
			terms= Table[
				If[FreeQ[term,coeff,All],
					term
					,
					(* if term contains coupling with wrong flavor indices, check if further couplings are present, if so relabel them in the variable list var, then relabel the entire term *)
					var= var/. {y:Alternatives@@Cases[term,(x:Coupling[Except[coeff],___]),All]:>(y/.indexReplacement)};
					term/.indexReplacement
				]
				,{term,terms}
			]
			,{coeff,mislabeledTerms}
		];
		sol=Quiet@Solve[#==0&/@terms,var]/.HoldPattern[InverseFunction[Bar,1,1]]->Bar;
	,"Solving the matching conditions for the coefficients of the target Lagrangian \[Ellipsis]"];
	
	(* replace temporary couplings in solutions and expand the result *)
	If[sol==={},
		$TempCouplingRules = DeleteCases[$TempCouplingRules, HoldPattern[Coupling[f_/;MemberQ[newcp,f],__] :> _]];
		$TempCouplings = DeleteCases[$TempCouplings, Alternatives@@newcp];
		RemoveCoupling/@newcp;
		Message[MapEffectiveCouplings::SolutionNotFound];
		Abort[]
	,
		OptionalMonitor[OptionValue@Verbose,
			output=ReplaceEffectiveCouplings[ OpenIndexToPattern/@Last@sol, newcp];
			output=MapAt[RelabelIndices,output,{All,2}]; (* canonically relabel indices on rhs *)
			$TempCouplingRules = DeleteCases[$TempCouplingRules, HoldPattern[Coupling[f_/;MemberQ[newcp,f],__] :> _]];
			$TempCouplings = DeleteCases[$TempCouplings, Alternatives@@newcp];
			RemoveCoupling/@newcp;
		,"Simplifying the results \[Ellipsis]"];
	];
	
	(* This might brake at subleading power in the EFT expansion... *)
	Module[{maxOrder,dimWC},
		(* maximum EFT order of the target Lagrangian *)
		maxOrder = Max[OperatorDimension/@TermsToList[TargetLag]]; (* we should add a check that this agrees with the same number obtained for the input Lagrangian *)
		(* Loop through matching conditions*)
		output = Table[
			(* determine EFTOrder of all terms containing this coupling and take the lowest occuring order *)
			dimWC = Min[OperatorDimension/@Cases[TermsToList[TargetLag],x_/;!FreeQ[x,First@cond,All]]];
			(* truncate RHS at EFTOrder given my max order we work at, minus lowest order of terms containing this coupling, plus the EFTOrder of the coupling itself *)
			First@cond -> SeriesEFT[Last@cond,EFTOrder->(maxOrder-dimWC+Last@First@cond)]
			,
			{cond,output}
		];
	];

	(* order output by EFT order *)
	If[OptionValue@ SortByEFTOrder,
		OptionalMonitor[OptionValue@Verbose,
			orderAssociation = Merge[Table[Association@DeleteDuplicates@Cases[SeriesEFT[TargetLag,EFTOrder->{o}], Coupling[l_,___]:>(l->o),Infinity],{o,GetMinOrder@TargetLag,GetMaxOrder@TargetLag}],First];
			output = SortBy[output, orderAssociation@ #[[1,1]]&];
		,"Sorting results by EFT order \[Ellipsis]"]
	];

	(* delete the rules that, after replacing couplings, become trivial *)
	If[!OptionValue@KeepTrivalReplacements,
		output = Select[output,!TrivialRuleQ[#]&]
	];

	(* delete duplicate solutions *)
	output= DeleteDuplicates[output];

	(* check if still multiple rules for the same coupling have been found *)
	If[!DuplicateFreeQ[output[[;;,1,1]]],
		Module[{tmp},
			tmp=FilterRules[output,pat_/;!FreeQ[pat,Alternatives@@(DeleteCases[Tally[output[[;;,1,1]]],{_,1}][[;;,1]]),All]];
			Print["The system of equations was overdetermined. The following solutions might be inconsistent:"];
			Print@Format[tmp,NiceForm];
			Print["Only the first solution is taken for the output"];
			output = DeleteCases[output, Alternatives@@tmp[[2;;]]];
		]
	];

	(* If set to true, add replacements for effective couplings in Lagrangian *)
	If[OptionValue@AppendEffectiveCouplingsDefs === True,
		(* check for remaining temporary couplings that where already present in the input, add their definitions to the output *)
		If[!FreeQ[output[[;;,2]],Alternatives@@$TempCouplings,All],
			output= Join[
				output,
				(PrintEffectiveCouplings[output[[;;,2]], Rules->True]/.(RelabelIndices[arg_,___]:>RelabelIndices[arg])/.RuleDelayed->Rule)/.Rule[a_,b_]:>Rule[a, RelabelIndices[b]]
			];
		];
	];
		
	(* go back to original options for CollectCoefficients *)
	Options[CollectCoefficients] = opt$CollectCoefficients;
	
	Return@ output
]


(* ::Subsection:: *)
(*Replacing with dummy indices*)


(* ::Subsubsection::Closed:: *)
(*Replacing symbols in a Lagrangian*)


(* ::Text:: *)
(*Function for replacing symbols in a Lagrangian expressions in a way so as to maintain integrity of the dummy indices after substitution*)


Options@ ReplaceInLagrangian= {Simplify-> True};


OptionTest[ReplaceInLagrangian, Simplify]= BooleanQ;
OptionMessage[Simplify, ReplaceInLagrangian, val_]:= Message[General::optexpectsval, BooleanQ, ReplaceInLagrangian, val, "Boolean"];


ReplaceInLagrangian[expr_, rule:(_Rule|_RuleDelayed), opts:OptionsPattern[]]:= 
	ReplaceInLagrangian[expr, {rule}, opts];


ReplaceInLagrangian[expr_, rules:{(_Rule|_RuleDelayed)...}, OptionsPattern[]]? OptionsCheck:= Module[
		{newRules, out},
	(*Current implementation always evaluates RHS in order to ensure initial uniqueness of dummy indices. 
	If we wish to support true RuleDelayed behavior a different approach is needed*)
	(*newRules= rules/. (Rule|RuleDelayed)[lhs_, rhs_]:> With[{newrhs= RelabelIndices[rhs, Unique-> True]},
		RuleDelayed[lhs, RelabelIndices[newrhs, Unique-> True]] ];*)
	(*Damnit! Stupid Mathematica! Insertion directly into a new delayed rule changes the pattern names on LHS for
	no f***ing reason! Hence, the hack with first a List and a Hold. Let's transition to Symbolica. I've had it!*)
	newRules= rules/. (Rule|RuleDelayed)[lhs_, rhs_]:> With[{newrhs= RelabelIndices[rhs, Unique-> True]},
		RuleDelayed@@ List[lhs, Hold@ RelabelIndices[newrhs, Unique-> True]] ]//ReleaseHold;
	out= expr/. x:Power[_, _Integer? Positive]:> PseudoTimes@ x/. newRules// ReleasePseudoTimes;
	
	(* hbar expansion to account for possible powers in hbar in the denominator *)
	out=Normal@Series[out,{hbar,0,1}];
	
	If[OptionValue@ Simplify, 
		GreensSimplify@ out
	,
		out
	]	
]


(* ::Section:: *)
(*SMEFT export*)


(* ::Subsection:: *)
(*Main export function*)


(* ::Text:: *)
(*Function for exporting a dimension-6 SMEFT *)


(* ::Text:: *)
(*Known limitations: *)
(*	SM fields and gauge groups need to be defined as in the "SM" or "SMEFT" model files. No flavor merging is currently possible.*)


(* ::Subsubsection::Closed:: *)
(*ExportSMEFT*)


ExportSMEFT::filename= "The filename must be a valid string.";
ExportSMEFT::filetype= "The filetype must be \"JSON\".";


Options@ ExportSMEFT= {
	FileType-> "JSON", (*More options could be added in the future*)
	ParameterNames-> Automatic,
	Verbose-> True
};


OptionTest[ExportSMEFT, ParameterNames]= 
	# === Automatic || MatchQ[#, {Rule[_Symbol, _String? PrintableASCIIQ]...}] &;


OptionMessage[ExportSMEFT, ParameterNames, val_]:= Message[General::optexpectsval, ParameterNames, ExportSMEFT, val, 
	"value Automatic or a list of rules from coupling labels to printable ASCII strings."]


ExportSMEFT[fileName_, lagEFT_, opts: OptionsPattern[]]:= Module[{coefficientList, newFileName, options},
	(*Check filename*)
	If[!MatchQ[fileName, _String],
		Message[ExportToSmelli::filename]; Abort[]; ];

	coefficientList= MapToWarsaw[lagEFT];

	Switch[OptionValue@ FileType
	,"JSON",
		newFileName= StringReplace[fileName<> ".json",
			RegularExpression["(\\.(?i)json\\.json)"]:> ".json"];
		options= FilterRules[{opts}, Options@ ExportAsJSON];
		ExportAsJSON[newFileName, coefficientList, Sequence@@ options];
	,_,
		Message[ExportToSmelli::filetype]; Abort[];
	];
]


(* ::Subsubsection::Closed:: *)
(*Obtain Warsaw basis coefficients*)


(* ::Text:: *)
(*Maps the EFT Lagrangian to the Warsaw basis provided by the SMEFT.m model file*)


MapToWarsaw::warsawReps= "Not recognized input.";


MapToWarsaw[lagEFT_]:= Module[{lagWarsaw, warsawCoefReplacements},
	lagWarsaw= LoadModel["SMEFT"];
	(*To be removed!*)
	Echo@ "The Warsaw Lagrangian currently simplifies with 4D IDENTITIES!";
	lagWarsaw= GreensSimplify[lagWarsaw, ReductionIdentities-> FourDimensional];

	(*Get solutions for the Warsaw basis coefficients*)
	warsawCoefReplacements= QuietEcho@ MapEffectiveCouplings[lagEFT, lagWarsaw];
	If[!MatchQ[warsawCoefReplacements, {Rule[_Coupling, _]..}],
		Message[MapToWarsaw::warsawReps]; Abort[]; ];
	warsawCoefReplacements[[;;, 2]]= ReplaceEffectiveCouplings/@ warsawCoefReplacements[[;;, 2]];
	
	(*We need to express the EFT coefficients in terms of the SMEFT renormalizable parameters rather than 
	the corresponding UV parameters*)
	warsawCoefReplacements
]


(* ::Subsubsection::Closed:: *)
(*JSON export*)


(* ::Text:: *)
(*For creating a jelli compatible json file with the SMEFT matching coefficients as contained in the Warsaw coefficient replacement list.*)


Options@ ExportAsJSON= {
	ParameterNames-> Automatic,
	Verbose-> True
};


ExportAsJSON::parnames= "Option \"ParameterNames\" must be a list of replacement rules for Symbols to valid python variable name (as Strings).";


ExportAsJSON[fileName_, warsawCoefReplacements_, OptionsPattern[]]:= Module[
		{ReplaceSymb, conjDependence, conjAssociation, coefDict, dim, einsums, flavorIndices, heavyMasses, internalFlavorAlphabet, lagWarsaw, massFunctions, massFunctionIndices, MF, nameReplacements, out, parameters, symb, temp, n},
	(*Change to first replace the specified symbols according to Option
		+ identify SM parameters*)
	If[OptionValue@ ParameterNames === Automatic,
		nameReplacements= {};
	,
		If[!MatchQ[OptionValue@ ParameterNames, {Rule[_Symbol, _String? PrintableASCIIQ]...}],
			Message[ExportAsJSON::parnames]; Abort[]];
		nameReplacements= OptionValue@ ParameterNames;
	];
	ReplaceSymb@ symb_:= symb/. nameReplacements/. s_Symbol:> SymbolToASCII@ s;

	(*Identify parameters*)
	parameters= DeleteDuplicates@ Cases[warsawCoefReplacements[[;;, 2]],
		Coupling[label_, __]-> label, All];
	parameters= Query[Key/@ parameters,
		{Key@ Indices, Key@ SelfConjugate}]@ $CouplingAssociation;
	parameters= KeyValueMap[Function[{key, val}, Module[{name= ReplaceSymb@ key},
		MapAt[$FlavorIndices[#, IndexDimension]&, If[val@ SelfConjugate === False,
			{{name, name <> "_conj" },
				{{name, val@ Indices}, {name <> "_conj", val@ Indices}}},
			{{}, {{name, val@ Indices}} }
		], {2, All, 2, All} ] ]
		], parameters];

	conjDependence= Cases[parameters[[;;, 1]], {_, _}];
	conjAssociation= Association@@ Rule@@@ conjDependence;
	parameters= Flatten[parameters[[;;, 2]], 1];

	(*Determine flavor indices*)
	flavorIndices= {};
	Do[
		If[!MatchQ[dim, _Symbol], Continue[];];
		If[FreeQ[warsawCoefReplacements[[;;, 2]], dim] && FreeQ[parameters, dim],
			Continue[]; ];
		AppendTo[flavorIndices, dim];
	, {dim, Values@ $FlavorIndices}];

	(*Determine heavy masses*)
	heavyMasses= List@@ Query[Key/@ GetFieldsByProperty[Heavy-> True], Key@ Mass]@ $FieldAssociation;

	(*Repalcement rules to smelli dictionary*)
	coefDict= {#/. Coupling[c_, __]-> c/. $smelliCoefNames,
		#/. warsawCoefReplacements// LogsToLF}&/@ $warsawCoefficientsWithSmelliFlavors;

	(*Match Warsaw basis conventions*)
	(*Flip sign of gauge couplings (from the sing in covariant derivatives)*)
	coefDict= coefDict/. c:Coupling[lab_, __]/;
		MemberQ[Values@ Query[All, Key@ Coupling]@ $GaugeGroups, lab]-> -c;

	(*Fix mass functions in coefficients*)
	OptionalMonitor[OptionValue@ Verbose, 
		Do[
			symb= coefDict[[n, 1]];
			coefDict[[n, 2]]= CanonizeAndCollectCoefficients@ coefDict[[n, 2]];
		, {n, Length@ coefDict}];
	, StringForm["Gathering mass functions of `1` (`2`/64)", symb, n] ];
	
	(*Determine 'mass_functions'*)
	OptionalMonitor[OptionValue@ Verbose, 
		massFunctions= DeleteDuplicates@ Cases[coefDict, _MassFunction, All];
		massFunctions= MapIndexed[#1-> "mass_function_" <> ToString@@ #2&, massFunctions];
		coefDict= coefDict/. massFunctions;
		massFunctions= List@@@ Reverse/@ massFunctions/. MassFunction-> Identity;
		massFunctions[[;;, 2]]= CollectMassPowers/@ massFunctions[[;;, 2]];
		
		(*Setup association with the dummy indices of each mass function*)
		massFunctionIndices= Association@@ 
			(#1-> DeleteDuplicates@ Cases[#2, Index[ind_, _]-> ind, All] &)@@@ massFunctions;
		
		(*Mass functions to python form*)
		(*Replace LF to prevent automatic collapsing of powers when indices are removed from mass labels*)
		massFunctions[[;;, 2]]= massFunctions[[;;, 2]]/. LF-> LFtemp/. 
			Coupling[lab_, _, _]:> ReplaceSymb@ lab; 
		massFunctions[[;;, 2]]= StringDelete[#, "\""]&/@ ToString/@ FortranForm/@ (massFunctions[[;;, 2]]/.
		(head:LFtemp|MassPower)[masses_, powers_]:> StringJoin[
			Switch[head, LFtemp, "LF", MassPower, "MassPower"], "(", 
			StringRiffle[masses, {"[", ", ", "]"}], ", ",
			StringRiffle[powers, {"[", ", ", "]"}], ")"]);
	, "Collect mass functions"];
	
	(*Gather Einsums*)
	OptionalMonitor[OptionValue@ Verbose, 
		Do[
			symb= coefDict[[n, 1]];
			temp= BetterExpand@ coefDict[[n, 2]]/. CouplingProduct-> Identity;
			temp= temp/. 
				Bar@ Coupling[lab_, inds_, _]:> Einsum[{inds[[;;, 1]]}, {Conj@ lab}]/.
				Coupling[lab_, inds_, _]:> Einsum[{inds[[;;, 1]]}, {lab}]/.
				mf_String:> Einsum[{massFunctionIndices@ mf}, {mf}]/. 
				d_Delta:> Einsum[{List@@ d[[;;, 1]]}, {"kroneckerdelta"}]/.
				_FlavorSum:> 1;
			temp= temp//. 
				Einsum[aInds_, aTens_]* Einsum[bInds_, bTens_]:> Einsum[Join[aInds, bInds], Join[aTens, bTens]]/.
				Power[Einsum[aInds_, aTens_], 2]:> Einsum[{aInds, aInds}, {aTens, aTens}];
			coefDict[[n, 2]]= Collect[temp/. e_Einsum:> CanonizeEinsum@ e,
				{_Einsum}, Simplify];
		, {n, Length@ coefDict}];
	, StringForm["Gathering Einsums of `1` (`2`/64)", symb, n] ];
	
	(*Determine 'einsums'*)
	OptionalMonitor[OptionValue@ Verbose, 
		einsums= MapIndexed[#1-> "einsum_" <> ToString@@ #2&, DeleteDuplicates@ Cases[coefDict, _Einsum, All]];
		coefDict= coefDict/. einsums;
		einsums= List@@@ Reverse/@ einsums;
		
		(*Einsums to python form*)
		(*specify the external coefficients (which are then not summed)*)
		einsums[[;;, 2]]= MapAt[Replace[#, s_Symbol:> ReplaceSymb@ s, -1]&,
			MapAt[(# <> "->" <> StringJoin@@ Sort@ DeleteDuplicates@ StringCases[#, RegularExpression["i|j|k|l"]]&)@* 
				(StringRiffle[#, ",", ""]&), einsums[[;;, 2]],
			{All, 1}], {All, 2}]/. Conj@ s_:> conjAssociation@ s;
		einsums[[;;, 2]]= einsums[[;;, 2]]/. Einsum[inds_, labs_]:> "einsum(\""<> inds<> "\", "<>
			StringRiffle[labs, ", "]<> ")";
	, "Collect Einsums"];

	(*Final symbol replacement*)
	coefDict[[;;, 2]]= Replace[coefDict[[;;, 2]]/. hbar-> 1/(16 \[Pi]^2)/. Bar-> Conj,
		s_Symbol? (Not@* NumericQ):> ReplaceSymb@ s, -1]/.
		Conj@ s_:> conjAssociation@ s;
	(*Python form*)
	coefDict[[;;, 2]]= StringDelete[ToString@ FortranForm@ #, "\""]&/@ coefDict[[;;, 2]];

	(*Construct JSON output*)
	out= {
		{"parameters", Map[If[Head@ # === Symbol, ToString@ #, #]&, parameters, {-1}]},
		{"conjugate_parameters", conjDependence},
		{"flavor_indices", ToString/@ flavorIndices},
		{"masses", ReplaceSymb/@ heavyMasses},
		{"einsums", einsums},
		{"mass_functions", massFunctions},
		{"wilson_coefficients", coefDict}
		};

	Export[fileName, out, "JSON"];
]


(* ::Subsection:: *)
(*Utility functions*)


(* ::Subsubsection::Closed:: *)
(*Loop functions instead of logs*)


(* ::Text:: *)
(*Insert loop function in place of the logs*)


LogsToLF@ expr_:= expr/. Log[m_^-2* \[Mu]bar2]-> LF[{m}, {1, 0}]/ m^2 - 1;


Einsum[{{}}, {lab_}]= lab;


internalFlavorAlphabet= CharacterRange["A", "W"];
CanonizeEinsum@ Einsum[inds_, labs_]:= Module[{indReplacement, out, ord= Ordering@ labs},
	out= Einsum[indReplacement= inds[[ord]], labs[[ord]]];
	(*i, j, k, l are always external flavor indices*)
	indReplacement= DeleteCases[DeleteDuplicates@ Flatten@ indReplacement, i | j | k | l];
	(*indReplacement= Cases[Tally@ Flatten@ indReplacement, {i_, 2}-> i];*)
	indReplacement= Thread@ Rule[indReplacement, internalFlavorAlphabet[[;;Length@ indReplacement]]];
	out/. indReplacement/. ext:Alternatives[i, j, k, l]:> ToString@ ext
]


(* ::Subsubsection::Closed:: *)
(*Canonize loop functions in coefficients*)


(* ::Text:: *)
(*Collect a coefficient by non-mass couplings, then simplify the mass function  *)


CanonizeAndCollectCoefficients@ coef_:= Module[{heavyMasses, out},
	heavyMasses= List@@ Query[Key/@GetFieldsByProperty[Heavy-> True], Key@ Mass]@ $FieldAssociation;
	(*Couplings  (not including heavy masses) are gathered in each term*)
	out= BetterExpand@ coef/. c:Coupling[lab_/; !MemberQ[heavyMasses, lab], __]:> CouplingProduct[c]/.
		f:(_FlavorSum|_Delta):> CouplingProduct@ f/.
		Bar@ CouplingProduct@ c_:> CouplingProduct@ Bar@ c/.
		Power[CouplingProduct@ c_, n_]:> CouplingProduct@ Power[c, n]//.
		CouplingProduct@ x_ CouplingProduct@ y_:> CouplingProduct[x * y];
	(*Terms are grouped by what couplings they contain.*)
	(*The function of the heavy masses appearing with each coupling combination is simplified*)
	out= Collect[out, {hbar, _CouplingProduct}, 
		ExtractCommonFactorFromMFs@* DecomposeInMFs@* SimplifyMassFunction];
	out// Expand 
]


(* ::Text:: *)
(*Split terms into Mass functions grouped by what masses appear in each term*)


DecomposeInMFs@ expr_:= Module[{terms},
	terms= TermsToList@ expr;
	terms= GatherBy[terms, Sort@ DeleteDuplicates@ Cases[#, _Coupling, All]&];
	Plus@@ (MassFunction@* Plus)@@@ terms
]


(* ::Text:: *)
(*Extract the common factors appearing with the LFs in the mass function*)


ExtractCommonFactorFromMFs@ expr_:= expr/. mf_MassFunction:> ExtractCommonFactorFromMFs@ mf;


(*Maybe also extract common mass (without indices) in addition to numerical factors*)
ExtractCommonFactorFromMFs@ MassFunction@ a_?NumberQ:= a;
ExtractCommonFactorFromMFs@ mf_MassFunction:= Module[{out= Simplify/@ mf},
	out= out/. MassFunction@ Times[a_?NumberQ, rest__]:> a MassFunction@ Times@ rest//.
		MassFunction@ Times[a:(Coupling[_, {}, _] | Power[Coupling[_, {}, _], _]), rest__]:> a MassFunction@ Times@ rest;
	out= out/. MassFunction[a:(Coupling[_, {}, _] | Power[Coupling[_, {}, _], _])]-> a;
	out/. m_MassFunction:> ExtractCommonRationalFromMF@ m
];


ExtractCommonRationalFromMF@ MassFunction@ expr_:= Module[{out= Expand@ expr, nums},
	nums= Cases[out, Times[a_?NumberQ, _LF | _Coupling | Power[_Coupling, _], ___]-> a, All];
	If[Length@ nums === 0, Return@ MassFunction@ expr];
	nums= GCD@@ nums;
	out= out/ nums// Expand;
	(*Obtain a canonical sign*)
	{nums, out}*= Order[out, -out];
	
	nums* MassFunction[out]
];


(* ::Subsubsection::Closed:: *)
(*Collect mass powers in mass functions*)


(* ::Text:: *)
(*Gather up loose powers of masses MassPower function:*)
(*MassPower[{m1, m2, ...}, {a, b, ...}] = m1^a * m2^ b *...*)


CollectMassPowers@ expr_:= Module[{len, masses, out},
	masses= Sort@ DeleteDuplicates@ Cases[expr, _Coupling, All];
	len= Length@ masses;
	out= ReplaceShieldSubexpressions[expr/. Power[c_Coupling, n_]:> MassPower[n UnitVector[len, First@ FirstPosition[masses, c]]], 
		c_Coupling:> MassPower[UnitVector[len, First@ FirstPosition[masses, c]]], _LF];
	out= out//. MassPower[inds1_]* MassPower[inds2_]:> MassPower[inds1+ inds2];
	out/. MassPower[inds_]:> MassPower[masses, inds]
];


(* ::Subsubsection::Closed:: *)
(*Symbol to ASCII*)


SymbolToASCII@ symb_:= StringReplace[
	StringReplace[ToString[symb, CharacterEncoding-> "ASCII"], {
		RegularExpression["\\\\\[(Capital)([[:upper:]])(.*?)\\]"]:>
			"_" <> ToUpperCase@"$2" <> "$3" <> "_",
		RegularExpression["\\\\\[([[:upper:]])(.*?)\\]"]:>
			"_" <> ToLowerCase@"$1" <> "$2" <> "_",
		RegularExpression["([[:alpha:]])(\\d)"]:>
			"$1" <> "_" <> "$2" 
	}],
	{
		RegularExpression["^_"]:> "",
		RegularExpression["_$"]:> "",
		RegularExpression["_+"]:> "_"
	}
]
(*SymbolToASCII@ symb_:= StringReplace[Echo@ ToString[symb, CharacterEncoding-> "ASCII"], {
	RegularExpression["(^|.)\\\\\[(Capital)([[:upper:]])(.*?)]($|.)"]:>
		If["$1" === "", "", "$1_"] <> ToUpperCase@"$3" <> "$4" <> If["$5" === "", "", "_$5"],
	RegularExpression["(^|.)\\\\\[([[:upper:]])(.*?)]($|.)"]:>
		If["$1" === "", "", "$1_"] <> ToLowerCase@"$2" <> "$3" <> If["$4" === "", "", "_$4"]
}]*)


(* ::Subsubsection::Closed:: *)
(*Smelli coefficient information*)


(* ::Text:: *)
(*These lists must correspond to the couplings in "SMEFT.m" in the model files *)


$warsawCoefficientsWithSmelliFlavors= {
	Coupling[cllHH,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cdB,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cdd,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cdG,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cdH,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cduq,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cduu,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cdW,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[ceB,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[ced,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cee,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[ceH,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[ceu,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[ceW,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cG,{},0],
	Coupling[cGt,{},0],
	Coupling[cH,{},0],
	Coupling[cHB,{},0],
	Coupling[cHBox,{},0],
	Coupling[cHBt,{},0],
	Coupling[cHd,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHD,{},0],
	Coupling[cHe,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHG,{},0],
	Coupling[cHGt,{},0],
	Coupling[cHl1,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHl3,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHq1,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHq3,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHu,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHud,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cHW,{},0],
	Coupling[cHWB,{},0],
	Coupling[cHWt,{},0],
	Coupling[cHWtB,{},0],
	Coupling[cld,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cle,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cledq,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[clequ1,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[clequ3,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cll,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[clq1,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[clq3,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[clu,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqd1,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqd8,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqe,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqq1,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqq3,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqqq,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqqu,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqu1,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cqu8,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cquqd1,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cquqd8,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cuB,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cud1,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cud8,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cuG,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cuH,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cuu,{Index[i,Flavor],Index[j,Flavor],Index[k,Flavor],Index[l,Flavor]},0],
	Coupling[cuW,{Index[i,Flavor],Index[j,Flavor]},0],
	Coupling[cW,{},0],
	Coupling[cWt,{},0]
};


$smelliCoefNames= {
	cllHH -> "C_llphiphi",
	cdB -> "C_dB",
	cdd -> "C_dd",
	cdG -> "C_dG",
	cdH -> "C_dphi",
	cduq -> "C_duql",
	cduu -> "C_duue",
	cdW -> "C_dW",
	ceB -> "C_eB",
	ced -> "C_ed",
	cee -> "C_ee",
	ceH -> "C_ephi",
	ceu -> "C_eu",
	ceW -> "C_eW",
	cG -> "C_G",
	cGt -> "C_Gtilde",
	cH -> "C_phi",
	cHB -> "C_phiB",
	cHBox -> "C_phiBox",
	cHBt -> "C_phiBtilde",
	cHd -> "C_phid",
	cHD -> "C_phiD",
	cHe -> "C_phie",
	cHG -> "C_phiG",
	cHGt -> "C_phiGtilde",
	cHl1 -> "C_phil1",
	cHl3 -> "C_phil3",
	cHq1 -> "C_phiq1",
	cHq3 -> "C_phiq3",
	cHu -> "C_phiu",
	cHud -> "C_phiud",
	cHW -> "C_phiW",
	cHWB -> "C_phiWB",
	cHWt -> "C_phiWtild",
	cHWtB -> "C_phiWtildeB",
	cld -> "C_ld",
	cle -> "C_le",
	cledq -> "C_ledq",
	clequ1 -> "C_lequ1",
	clequ3 -> "C_lequ3",
	cll -> "C_ll",
	clq1 -> "C_lq1",
	clq3 -> "C_lq3",
	clu -> "C_lu",
	cqd1 -> "C_qd1",
	cqd8 -> "C_qd8",
	cqe -> "C_qe",
	cqq1 -> "C_qq1",
	cqq3 -> "C_qq3",
	cqqq -> "C_qqql",
	cqqu -> "C_qque",
	cqu1 -> "C_qu1",
	cqu8 -> "C_qu8",
	cquqd1 -> "C_quqd1",
	cquqd8 -> "C_quqd8",
	cuB -> "C_uB",
	cud1 -> "C_ud1",
	cud8 -> "C_ud8",
	cuG -> "C_uG",
	cuH -> "C_uphi",
	cuu -> "C_uu",
	cuW -> "C_uW",
	cW -> "C_W",
	cWt -> "C_Wtilde"
};
