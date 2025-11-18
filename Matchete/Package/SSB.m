(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`SSB`*)


(* ::Subtitle:: *)
(*Procedures for spontaneous symmetry breaking*)


(* ::Chapter:: *)
(*Public:*)


(* ::Subsection:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["SetSymmetryBreaking"]
PackageExport["RepresentationDecomposition"]
PackageExport["FieldDecomposition"]
PackageExport["CGDecomposition"]
PackageExport["BrokenPhase"]


PackageExport["Singlet"]
PackageExport["ComplexSinglet"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


(* ::Subsection:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Setting up the breaking pattern*)


(* ::Subsubsection::Closed:: *)
(*Index projectors*)


(* ::Text:: *)
(*The BrokenInds facilitate contraction between  the right broken components from various fields and CGs *)


ResetBrokenInds[]:= Module[{},
	Clear@ BrokenInds;
	BrokenInds/: Bar@ b_BrokenInds:= b;
	BrokenInds[Bar@ ind_, m_]:= BrokenInds[ind, m];
]


(* ::Subsection:: *)
(*Breaking pattern*)


(* ::Subsubsection::Closed:: *)
(*Groups*)


$symmetryBreaking= <||>;


SetSymmetryBreaking[group_, stabilityGroup_]:= Module[{},
	(*Check inputs are list of defined gauge groups*)
	ResetBrokenInds[];
	$symmetryBreaking= <|
			CGSubstitutions-> <||>,
			FieldSubstitutions-> <||>,
			OriginalGroup-> group,
			Representations-> <||>,
			StabilityGroup-> stabilityGroup
		|>;	
]


(* ::Subsubsection::Closed:: *)
(*Representations*)


RepresentationDecomposition[rep_, decomposition_]:= RepresentationDecomposition[rep, {decomposition}];
RepresentationDecomposition[rep_, decomposition_List]:= Module[{fs, count, conj, pairs},
	(*Check representations and dimensions*)
	$symmetryBreaking[Representations, rep]= decomposition;
	
	Switch[$Representations[rep, Reality]
		,1,
			(*Complex representations embedded into real representations of the original symmetry are 
			assumed to form conjugate pairs with the neighboring rep*)
			conj= False; count= 0;
			pairs= Table[
					r= r/. Bar-> Identity;
					fs= If[!MatchQ[r, Singlet|ComplexSinglet], 
							$Representations[r, Reality] === +1
						,
							False
						];
					count++;
					If[r === Singlet || fs,
						{count, count}-> 1
					,
						If[conj,
							conj= False;
							{count, count- 1}-> 1
						,
							conj= True;
							{count, count+ 1}-> 1
						]
					]
				, {r, decomposition}];
			pairs= Normal@ SparseArray@ pairs;
			With[{array= pairs},
				BrokenInds/: BrokenInds[Index[a_, rep], m_]BrokenInds[Index[a_, rep], n_]:= 
					array[[n, m]];
				BrokenInds/: Power[BrokenInds[Index[a_, rep], m_], 2]:= array[[m, m]]; 
			];
		,_,
			BrokenInds/: BrokenInds[Index[a_, rep], m_]BrokenInds[Index[a_, rep], n_]:= 
				KroneckerDelta[n, m];
			BrokenInds/: Power[BrokenInds[Index[a_, rep], m_], 2]:= 1;
	];
]


(* ::Subsection:: *)
(*Decomposition of objects*)


(* ::Subsubsection::Closed:: *)
(*Field decomposition*)


FieldDecomposition[NonCommutativeMultiply[DiracProduct@ _Proj, field_Field], decomposition_]:=
	FieldDecomposition[field, decomposition];
FieldDecomposition[field_Field, decomposition_]:= FieldDecomposition[field, {decomposition}];


FieldDecomposition[field_Field, decomposition_List]:= Module[{brokenInds, fieldInds, indReplace, depth, rhs, \[Mu], lhs, repInds},
	(*Checks open indices match in all components*)
	
	(*Determins the broken indices in field*)
	brokenInds= Cases[field[[3]], _Index? (MemberQ[GroupFromInd@ #]@ $symmetryBreaking@ OriginalGroup &), All];
	depth= Length@ brokenInds;
	
	fieldInds= Cases[field, _Index, All];
	indReplace= fieldInds/. ind:Index[lab_, rep_]:> 
		With[{temp= lab}, Rule[ind, Index[Pattern[temp, Blank[]], rep] ] ];
	fieldInds= fieldInds/. Index[l_, _]:> l;
	
	lhs= field/. indReplace;
	lhs[[-1]]= \[Mu]___;
	
	(*Check dimensions of input array*)
	rhs= MapIndexed[(#1 Times@@ Thread@ BrokenInds[brokenInds, #2]&), decomposition, {depth}];
	rhs= Plus@@ Flatten[rhs, depth];
	rhs= rhs/. Field[rest__, {}]-> Field[rest, \[Mu]];
	
	(*Make repeated indices in the replacement unique*)
	repInds= Complement[Cases[decomposition, Index[l_, r_]:> l, All], fieldInds];
	
	$symmetryBreaking[FieldSubstitutions, First@ field]= With[
		{temp= rhs, sub= (Index[#, rep_]-> Hold@ Index[Unique@ #, rep]&)/@ repInds},
			RuleDelayed[lhs, temp/. ReleaseHold@ sub]
		];
]


DecomposeFields@ expr_:= expr/. List@@ $symmetryBreaking@ FieldSubstitutions;


(* ::Subsubsection::Closed:: *)
(*CG decomposition*)


CGDecomposition[cg_, decomposition_]:= Module[{brokenInds, indReplace, depth, rhs, lhs},
	(*Check that cg is of a broken group*)
	(*Broken indices in field*)
	brokenInds= Cases[cg[[-1]], _Index, All];
	indReplace= brokenInds/. ind:Index[lab_, rep_]:> 
		With[{temp= lab}, Rule[ind, Index[Pattern[temp, Blank[]], rep] ] ];
	depth= Length@ brokenInds;
	
	lhs= cg/. indReplace;
	
	(*Check dimensions of input array*)
	rhs= MapIndexed[(#1 Times@@ Thread@ BrokenInds[brokenInds, #2]&), decomposition, {depth}];
	rhs= Plus@@ Flatten[rhs, depth];
	
	(*Make repeated indices in the replacement unique*)
	
	$symmetryBreaking[CGSubstitutions, First@ cg]= With[{temp= rhs},
			(*Make conjugate replacement pattern if CG has unique conjugate:*)
			If[$CGproperties[First@ cg, UniqueConj],
				{lhs:> temp, Bar@ lhs:> Bar@ temp}
			,
				{lhs:> temp}
			]
		];
] 


DecomposeCGs@ expr_:= expr/. Flatten[List@@ $symmetryBreaking@ CGSubstitutions];


(* ::Subsubsection::Closed:: *)
(*Decomposition of the structure constants*)


(* ::Text:: *)
(*Use the known decomposition of the generators of a representation to determine the decomposition of the structure constants*)


AddFStructDecomposition[group_, rep_]:= Module[{decomposition, dynk, a, b, c, X, Y, Z},
	(*Check group, representation, and generator decomposition*)
	
	dynk= DynkinIndex[$GaugeGroups[group, Group], 
		$Representations[rep, DynkinCoefficients] ];
	decomposition= -I/ dynk ContractCGs@ Expand@ DecomposeCGs[
		CG[gen@ rep, {X, a, b}] (CG[gen@ rep, {Y, b, c}] CG[gen@ rep, {Z, c, a}]- 
			CG[gen@ rep, {Z, b, c}] CG[gen@ rep, {Y, c, a}])]// Expand;
		
	$symmetryBreaking[CGSubstitutions, fStruct@ group]= With[
		{temp= decomposition, inds= Index[#, group@ adj]&/@ {X_, Y_, Z_}},
			CG[fStruct@ group, inds]:> temp
		];
	
	(*Check normalization of (f^ABC)^2 before/after symmetry breaking*)

]


(* ::Subsubsection::Closed:: *)
(*Decomposition of the Kronecker deltas*)


(* ::Text:: *)
(*Use the knowledge of a representation decomposition to determine the decomposition of a delta*)


AddDeltaDecomposition[brokenRep_]:= Module[{bar, bi1, bi2, decomposition, group, realBrokenRep, fs, repDecomposition, r, a, b},
	(*Check group and representation decomposition*)
	
	group= GroupFromRep@ brokenRep;
	realBrokenRep= $Representations[brokenRep, Reality]  === 1;
	repDecomposition= $symmetryBreaking[Representations, brokenRep];
	{bi1, bi2}= {Index[a, brokenRep], Bar@ Index[b, brokenRep]};
	
	conj= False; count= 0;
	decomposition= Sum[
			count++; 
			bar= Head@ r === Bar;
			If[bar, r= First@ r];
			
			fs= If[MatchQ[r, Singlet|ComplexSinglet], False, 
				$Representations[r, Reality]  === 1];
			If[realBrokenRep && !(r === Singlet || fs),
				conj= !conj;
				If[conj,
					BrokenInds[bi1, count] BrokenInds[bi2, count+ 1]
				,
					BrokenInds[bi1, count] BrokenInds[bi2, count- 1]
				]
			,
				BrokenInds[bi1, count] BrokenInds[bi2, count]
			] *
			If[MatchQ[r, Singlet|ComplexSinglet], 1,
				If[bar, Bar, Identity][CG[del@ r, {a, b}]]]
		, {r, repDecomposition}];
	
	(*Check normalization of \[Delta]_ab^2 before/after symmetry breaking*)
	
	$symmetryBreaking[CGSubstitutions, del@ brokenRep]= With[
		{dec= decomposition, inds= {Index[a_, brokenRep], Bar@ Index[b_, brokenRep]}},
			CG[del@ brokenRep, inds]:> dec
		];
]


(* ::Subsection:: *)
(*Derivatives*)


(* ::Subsubsection::Closed:: *)
(*Heavy vector action*)


(* ::Text:: *)
(*The action of a heavy vector on a field *)


VAction[\[Mu]_, f:(Field[lab_, _, inds_, _]| FieldStrength[lab, _, inds_, _])]:= Module[{ind, charge, newInd, newLab},
	Sum[
		ContractDelta@ If[MemberQ[$symmetryBreaking@ OriginalGroup, GroupFromInd@ ind],
			newInd= If[Head@ ind === Bar, Bar, Identity]@ Index[newLab, Last[ind/. Bar-> Identity]];
			- I VActionOnIndex[\[Mu], ind, newInd] BrokenPhaseCD[f/. ind-> newInd] 
		,
			0
		]
	, {ind, inds}] + Sum[
		If[MemberQ[$symmetryBreaking@ OriginalGroup, Head@ charge],
			- I VActionOnCharge[\[Mu], charge]
		,
			0
		]
	,{charge, GetFields[lab, Charges]}] * BrokenPhaseCD@ f
]


VActionOnIndex[\[Mu]_, brokenInd_, newInd_]:= Module[{generator, group, l, adjInd},
	group= GroupFromInd@ brokenInd;
	adjInd= Index[l, group@ adj];
	generator= If[Head@ brokenInd === Bar,
			-Bar@ CG[gen@ brokenInd[[1, 2]], {adjInd, Bar@ brokenInd, newInd}]
		,
			CG[gen@ brokenInd[[2]], {adjInd, brokenInd, Bar@ newInd}]
	];
	ContractDelta@ Expand@ LightGaugeFieldsToZero@ DecomposeCGs@ DecomposeFields[
		(*$GaugeGroups[group, Coupling][]*) $GaugeGroups[group, Field][\[Mu], l] generator
	]
]


VActionOnCharge[\[Mu]_, brokenGroup_@ charge_]:= LightGaugeFieldsToZero@ DecomposeFields[
		charge (*$GaugeGroups[brokenGroup, Coupling][]*) $GaugeGroups[brokenGroup, Field][\[Mu]]
	]


LightGaugeFieldsToZero@ expr_:= expr/.
	(Field[#, __]:> 0& )/@ List@@ Query[Key/@ $symmetryBreaking@ StabilityGroup, Key@ Field]@ $GaugeGroups;


(* ::Subsubsection::Closed:: *)
(*Decomposing derivatives*)


(* ::Text:: *)
(*We identify the action of a heavy vector in the decomposition [D_mu, f] = [d_mu, f] -i [V_mu, f] of the full covariant derivative on a field.   *)


BrokenPhaseCD[Field[lab__, inds_, {\[Mu]_, \[Nu]___}]]:= Module[{f= Field[lab, inds, {\[Nu]}]},
	CD[\[Mu], BrokenPhaseCD@ f] + VAction[\[Mu], f]
]
BrokenPhaseCD[f:Field[___, {}]]:= DecomposeFields@ f;


BrokenPhaseCD[FieldStrength[rest__, inds_, {\[Mu]_, \[Nu]___}]]:= Module[{f= FieldStrength[rest, inds, {\[Nu]}]},
	CD[\[Mu], BrokenPhaseCD@ f] + VAction[\[Mu], f] 
]
BrokenPhaseCD[FieldStrength[lab_, {\[Mu]_, \[Nu]_}, inds_, {}]]:= Module[{FSs, $a, $b},
	(*Include FSs of the light fields*)
	FSs= DecomposeFields@ Field[lab, Vector@ \[Mu], inds, {}]/. 
		Field[l_, _, gInds_, {}]:> FieldStrength[l, {\[Mu], \[Nu]}, gInds, {}];
	FSs+ LightGaugeFieldsToZero@ VAction[\[Mu], Field[lab, Vector@ \[Nu], inds, {}]]
]


(* ::Subsection:: *)
(*Group algebra of the broken phase*)


(* ::Text:: *)
(*To compute how the structure constant of the product gauge group decomposes in the broken phase*)


(* ::Subsubsection::Closed:: *)
(*Determine heavy vectors *)


(* ::Text:: *)
(*The heavy vectors of the breaking pattern are*)


HeavyVectors[]:= Module[{fields},
	fields= DeleteDuplicates@ Cases[List@@ $symmetryBreaking[FieldSubstitutions][[;;, 2]], Field[lab_, __]-> lab, All];
	Intersection[GetFieldsByProperty[Type-> Vector, Heavy-> True], fields]
];


HeavyVectorsDOFs[]:= Module[{vecs},
	vecs= HeavyVectors[];
	Flatten[If[GetFields[#, SelfConjugate], {#}, {#, Conj@ #}]&/@ vecs]
]


(* ::Text:: *)
(*The light vectors are the gauge fields themself *)


LightVectors[]:= Module[{fields},
	fields= DeleteDuplicates@ Cases[List@@ $symmetryBreaking[FieldSubstitutions][[;;,2]], Field[lab_, __]-> lab, All];
	Intersection[GetFieldsByProperty[Type-> Vector, Mass-> 0], fields]
];


(* ::Subsubsection::Closed:: *)
(*Structure constants of 3 heavy vectors  *)


HeavyIndexMetric[]:= Module[{conj, count, vec, dofs},
	dofs= HeavyVectorsDOFs[]/. Conj-> Identity;
	conj= False; count= 0;
	Normal@ SparseArray@ Table[
			count++; 
			If[GetFields[vec, SelfConjugate], 
				{count, count}-> 1
			,
				conj= !conj;
				If[conj,
					{count, count+ 1}-> 1
				,
					{count, count- 1}-> 1
				]
			]
		, {vec, dofs}] 	
]


GaugeIndexDecompositionHeavy[]:= Module[{decomposition, gaugeFields, f, mu, inds, out, vec, heavies},
	gaugeFields= Intersection[Keys@ $symmetryBreaking[FieldSubstitutions], 
		GetFieldsByProperty[Type-> Vector, Heavy-> False] ];
	heavies= HeavyVectorsDOFs[];
	
	out= A|-> Evaluate[Table[
		inds= Sequence@@ ConstantArray[A, Length@ GetFields[f, Indices]];
		First@ GetGaugeGroupByProperty[Field-> f]-> 
			Table[
				f[mu, inds]/. $symmetryBreaking[FieldSubstitutions, f]/. 
					If[Head@ vec === Conj, {Bar@ Field[First@ vec, __]-> 1}, 
						{Bar@ Field[vec, __]-> 0, Field[vec, __]-> 1}]/. 
					_Field-> 0
			, {vec, heavies}]
	, {f, gaugeFields}]];
	Apply[Association]@* out
]


(* ::Text:: *)
(*Returns the structure constants restricted to the space and basis of the heavy vectors, f^{i}_{jk}*)


StructureConstant3H[]:= Module[{cg, ha, hb, hc, hi, gr, groups},
	hi= GaugeIndexDecompositionHeavy[];
	groups= Keys@ hi@ a;
	(*Tensor multiplying the heavy vector metric to raise the first index.*)
	{a, b, c}|-> Evaluate[TensorContract[TensorProduct[HeavyIndexMetric[], Sum[
		If[!$GaugeGroups[gr, Abelian],
			cg= $GaugeGroups[gr, Coupling] DecomposeCGs@ CG[fStruct@ gr, {a, b, c}];
			Table[
				cg ha hb hc// Expand
			, {ha, hi[a][gr]},{hb, hi[b][gr]},{hc, hi[c][gr]}]
		,
			0
		]
	, {gr, groups}] ], {2, 3}]]
]


(* ::Subsubsection::Closed:: *)
(*Structure constants with mixed heavy and light vectors  *)


LightIndexMetric[]:= Module[{vec},
	DiagonalMatrix@ Table[
			$GaugeGroups[First@ GetGaugeGroupByProperty[Field-> vec], Coupling]^2
		, {vec, LightVectors[]}] 	
]


GaugeIndexDecompositionLight[]:= Module[{decomposition, unbrokenGaugeFields, f, mu, inds, out, vec, lights},
	unbrokenGaugeFields= Intersection[Keys@ $symmetryBreaking[FieldSubstitutions], 
		GetFieldsByProperty[Type-> Vector, Heavy-> False] ];
	lights= LightVectors[];
	
	out= A|-> Evaluate[Table[
		inds= Sequence@@ ConstantArray[A, Length@ GetFields[f, Indices]];
		First@ GetGaugeGroupByProperty[Field-> f]-> 
			Table[
				f[mu, inds]/. $symmetryBreaking[FieldSubstitutions, f]/. 
					Field[vec, __]-> 1/. _Field-> 0
			, {vec, lights}]
	, {f, unbrokenGaugeFields}]];
	Apply[Association]@* out
]


(* ::Text:: *)
(*Returns the structure constants restricted to the space and basis of the heavy vectors, f^{i}_{\alpha j}*)


StructureConstant2H1L[]:= Module[{cg, ha, hb, hc, heavyIndex, lightIndex, lMet12, gr, groups},
	heavyIndex= GaugeIndexDecompositionHeavy[];
	groups= Keys@ heavyIndex@ a;
	lightIndex= GaugeIndexDecompositionLight[];
	lMet12= DiagonalMatrix@ Table[
			$GaugeGroups[First@ GetGaugeGroupByProperty[Field-> vec], Coupling]^(-1)
		, {vec, LightVectors[]}];
	(*Tensor multiplying the heavy vector metric to raise the first index.*)
	{a, b, c}|-> Evaluate[Expand@ TensorContract[TensorProduct[HeavyIndexMetric[], lMet12, Sum[
		If[!$GaugeGroups[gr, Abelian],
			cg= $GaugeGroups[gr, Coupling] DecomposeCGs@ CG[fStruct@ gr, {a, b, c}];
			Table[
				cg ha lb hc// Expand
			, {ha, heavyIndex[a][gr]},{lb, lightIndex[b][gr]},{hc, heavyIndex[c][gr]}]
		,
			0
		]
	, {gr, groups}] ], {{2, 5}, {4, 6}}]]
]


(* ::Subsection:: *)
(*Scalar fields and VEVs*)


(* ::Subsubsection::Closed:: *)
(*Determine scalar multiplets receiving a VEV *)


(* ::Text:: *)
(*A list of all the scalar fields receiving a VEV *)


SymmetryBreakingScalars[]:= Module[{scalars, f, vev},
	scalars= Intersection[GetFieldsByProperty[Type-> Scalar], Keys@$symmetryBreaking[FieldSubstitutions]];
	(*Select scalars with a non-zero VEV*)
	Table[
		vev= $symmetryBreaking[FieldSubstitutions, f][[2]]/. _Field-> 0;
		If[vev === 0, Nothing, f]
	, {f, scalars}]
];


SymmetryBreakingScalarsDOFs[]:= Module[{scals},
	scals= SymmetryBreakingScalars[];
	Flatten[If[GetFields[#, SelfConjugate], {#}, {#, Conj@ #}]&/@ scals]
]


(* ::Text:: *)
(*VEVs of all the scalars developing a VEV*)


ScalarVEVs[]:= Module[{scalars, f, inds, out, vev},
	scalars= SymmetryBreakingScalars[];
	out= a|-> Evaluate[Table[
		inds= Sequence@@ Array[a, Length@ GetFields[f, Indices]];
		f-> DecomposeFields[f@ inds]/. _Field-> 0
	, {f, scalars}]];
	Apply[Association]@* out
];


(* ::Text:: *)
(*All the scalar degrees of freedom embedded in those fields*)


ScalarFields[]:= Module[{scalars, f, inds, out, vev, fields},
	scalars= SymmetryBreakingScalars[];
	out= a|-> Evaluate[Table[
		inds= Sequence@@ Array[a, Length@ GetFields[f, Indices]];
		(*Subtract of the vev*)
		fields= DecomposeFields[f@ inds];
		f-> Expand[fields- (fields/. _Field-> 0)]
	, {f, scalars}]];
	Apply[Association]@* out
];


(* ::Subsubsection::Closed:: *)
(*Broken generators *)


(* ::Text:: *)
(*Returns the broken generator x^a_{ib} acting on the representation of a scalar field. *)


BrokenGeneratorsOnField[field_]:= Module[{charges, rep, reps, indNo, indI, indJ, n, k},
	(*Adjust the index representations and the charges depending on conjugation of field*)
	{reps, charges}= If[Head@ field === Conj, 
		{Bar/@ GetFields[First@ field, Indices], MapAt[Minus, GetFields[First@ field, Charges], {All, 1}]}
	,
		{GetFields[field, Indices], GetFields[field, Charges]}
	];
	indNo= Length@ reps;
	
	{A, i, j}|-> Evaluate[
		{indI, indJ}= {Array[i, indNo], Array[j, indNo]};
		(*Non-Abelian*)
		Expand@ Sum[
			Product[
				If[k === n,
					BrokenGeneratorOnRep[reps[[k]]][A, indI[[k]], indJ[[k]]]
				,
					DecomposeCGs@ If[Head@ reps[[k]] === Bar, 
						rep= First@ reps[[k]];
						Bar
					, 
						rep = reps[[k]];
						Identity
					]@ CG[del@ rep, {Index[indI[[k]], rep], Bar@ Index[indJ[[k]], rep]}]
				]
			, {k, indNo}]
		, {n, indNo}] + 
		(*Abelian*)
		Expand[
			Sum[
					BrokenGeneratorOnCharge[ch]
				, {ch, charges}] *
			Product[
					DecomposeCGs@ If[Head@ reps[[k]] === Bar, 
						rep= First@ reps[[k]];
						Bar
					, 
						rep = reps[[k]];
						Identity
					]@ CG[del@ rep, {Index[indI[[k]], rep], Bar@ Index[indJ[[k]], rep]}]
				, {k, indNo}]
		]
	]
]


BrokenGeneratorOnCharge@ group_@ q_:= Module[{A},
	q $GaugeGroups[group, Coupling] GaugeIndexDecompositionHeavy[][A][group]
]


BrokenGeneratorOnRep[representation_]:= Module[{rep, group, conj},
	rep= If[(conj= Head@ representation === Bar), First@ representation, representation];
	group= GroupFromRep@ rep;
	{A, i, j}|-> Evaluate@ Table[
		hA $GaugeGroups[group, Coupling] DecomposeCGs@ If[conj, Minus@* Bar, Identity]@ CG[gen@ rep, {A, i, j}]//Expand
	, {hA, GaugeIndexDecompositionHeavy[][A][group]}]
]


(* ::Subsubsection::Closed:: *)
(*Decay constants *)


(* ::Text:: *)
(*The decay constant defined by f_i^a = - i x^a_{ib} v^b*)


DecayConstantMatrix[field_]:= Module[{b, conj, scalar},
	scalar= If[(conj= Head@ field === Conj), First, Identity]@ field;
	{i, a}|-> Evaluate[
		Contract/@ (-I BrokenGeneratorsOnField[field][i, a, b] If[conj, Bar, Identity]@ ScalarVEVs[][b][scalar])
	]
]


(* ::Section:: *)
(*Prepare Lagrangian for the broken phase*)


(* ::Subsection:: *)
(*Going to the broken phase*)


(* ::Subsubsection::Closed:: *)
(*Broken phase*)


BrokenPhase@ expr_:= Module[{out},
	out= PseudoTimes@ RelabelIndices[expr, Unique-> True];
	out= DecomposeCGs@ out/. f:(_Field| _FieldStrength):> BrokenPhaseCD@ f/. Bar@ d_Delta-> d;
	out= Contract@ ContractCGs@ ContractDelta@ Expand@ ReleasePseudoTimes@ out;
	(*Remove vacuum*)
	out= out- (out/. (_Field| _FieldStrength)-> 0);
	CollectOperators@ out
]


(* ::Subsubsection::Closed:: *)
(*TrigonometricSimplify*)


$trigIdentities= {};


CoSinePair[c_, s_]:= (
	$trigIdentities= $trigIdentities~ Join~ {c^2+s^2-> 1, -c^2-s^2-> -1, c^4-s^4-> c^2-s^2, s^4-c^4-> s^2-c^2}; 
);


TrigSimplify@ expr_:= expr/. $trigIdentities;


LagrangianSimplify[expr_, rules_]:= Module[{},
	TrigSimplify@ CollectOperators[expr//. rules]
]


(* ::Subsection:: *)
(*Prepare Lagrangian *)


(* ::Text:: *)
(*Functions to prepare the broken phase Lagrangian for further computations*)


(* ::Subsubsection::Closed:: *)
(*Prepare *)


PrepareLagrangianForMatching::incKinTerm= "The Lagrangian cannot be canonically normalized. The kinetic term involving the fields `1` is normalized to `2`."; 


(* ::Text:: *)
(*Takes a Lagrangian in the broken phase and determine the conditions on the coefficients required to put it on a sensible form. *)
(*It returns the Lagrangian with these conditions applied, and the conditions are kept in $SImplificationConditions*)


$SimplificationConditions= <||>;


PrepareLagrangianForMatching@ lag_:= Module[{out, terms, conditions},
	$SimplificationConditions= <||>;
	(*out= SubstituteCoefficients@ lag;*)
	out= IntroduceEffectiveCouplings@ lag;
	terms= InternalSimplify[out, InternalOpRepresentation-> True];
	(*Tadpole coeficients*)
	$SimplificationConditions@ TadpoleConditions= ExtractTadpoleConditions@ terms;
	(*Off-diagonal masses*)
	$SimplificationConditions@ MassConditions= ExtractOffDiagonalMasses[out, terms];
	(*Coefficients of non-GB mixing terms*)
	$SimplificationConditions@ ScalarVectorMixingConditions= ExtractScalarVectorMixing[out, terms];
	(*Kinetic terms*)
	$SimplificationConditions@ KineticTermConditions= ExtractKineticTermConditions[out];
	
	(*Constrain Lagrangian*)
	conditions= Flatten[List@@ List@@@ $SimplificationConditions];
	out= out/.conditions// CollectOperators;
	
	(*Put simplification conditions in terms of the original couplings*)
	$SimplificationConditions= Query[All,All,ReplaceEffectiveCouplings]@ $SimplificationConditions;
	
	out
]


(* ::Subsubsection::Closed:: *)
(*Name coefficients *)


(* ::Text:: *)
(*To gather the coefficients of Lagrangian terms in single coefficients*)


(*$couplingNo= 1;
$coefSubstitutions= {};*)


(*SubstituteCoefficients@ expr_:= Module[{out, real, complex},
	out= InternalSimplify[expr, InternalOpRepresentation-> True];
	If[out === 0, Return@ out;];
	out= If[Head@ out =!= Plus, {out}, List@@ out];
	
	(*Include test to ensure that expression is explicitly Hermitian save for the kinetic terms*)
	
	real= Cases[out, (Times[op:(_AtomicOp|_CompOp), rest__]|op:(_AtomicOp|_CompOp)/; RealOpQ@ op):> 
		{op, Times@ rest}];
	complex= Cases[out, (Times[op:(_AtomicOp|_CompOp), rest__]|op:(_AtomicOp|_CompOp)/; !RealOpQ@ op):> 
		{op, Times@ rest}];
	complex= DeleteDuplicatesBy[complex, (Sort@ {#, MapAt[OpClassConjugate, #, 1]}&@* First@* First)];
	(*Update to factor out expected numeric factor*)
	$coefSubstitutions= Reap[
		real= MapAt[MakeCouplingPattern, real, {All, 2}];
		complex= MapAt[MakeCouplingPattern[#, False]&, complex, {All, 2}];
	][[2, 1]];
	
	out= Plus@@ Times@@@ real+
		Plus@@ (#[[1]]* #[[2]] + ConjugateAtomicOperator@ #[[1]]* Bar@ #[[2]] &)/@ complex;
	
	AtomicToNormalForm@ out
]*)


(* ::Text:: *)
(*Test if an AtomicOp is real or complex (exception if the operator is a kinetic term)*)


(*MakeCouplingPattern[expr_, real_:True]:= Module[{coup},
	(*Test if a new coupling should be introduced*)
	If[MatchQ[expr, numPattern| Times[numPattern, Alternatives[_Coupling, Bar@ _Coupling]] ],
		Return@ expr;
	];
	
	(*Include indices (external and internal)? *)
	If[!FreeQ[expr, _Index], Return@ expr;];
	
	coup= Symbol["cof"<> ToString@ $couplingNo++];
	(*Determine and include EFT order*)
	DefineCoupling[coup, SelfConjugate-> real];
	Sow[coup[]-> expr];
	coup[]
];
numPattern= Alternatives[_Rational, _Integer, _Complex];*)


(* ::Subsubsection::Closed:: *)
(*Find tadpole coefficients*)


(* ::Text:: *)
(*Finds the coefficients of all tadpole terms and set them to zero*)


ExtractTadpoleConditions@ expr_:= Module[{terms},
	(*terms= expr/. (AtomicOp| CompOp)[{Except@{{_}, 0}, _}, _]-> 0;*)
	terms= expr/. {AtomicOp[Except@ {{_}, 0}, __]-> 0, CompOp[_, Except@ {{_}, 0}, __]-> 0};
	If[terms === 0, Return@ <||>;];
	terms= If[Head@ terms === Plus, List@@ terms, {terms}];
	terms= terms/. Times[AtomicOp[{{f_}, 0}, __], rest___]:> Rule[f, Times@ rest];
	Rule[#, 0]&/@ Association@@ terms
]


(* ::Subsubsection::Closed:: *)
(*Find off-diagonal masses*)


(* ::Text:: *)
(*Finds all instances of off-diagonal masses and determines the coefficients required to eliminate them*)


ExtractOffDiagonalMasses[lag_, terms_]:= Module[{expr, fields, masses, i, j},
	(*Group fields by quantum numbers*)
	fields= DeleteDuplicates@ Cases[lag, Field[f_, __]:>f, All];
	fields= Query[Key/@ fields, Key/@ {Type, Indices, Charges}]@ GetFields[];
	fields= List@@(Keys/@ GroupBy[MapAt[Sort, fields, {All,Key@Indices}], Identity]);
	
	expr= AtomicToNormalForm[terms/. 
		{AtomicOp[Except@ {{_, _}, 0}, __]-> 0, CompOp[_, Except@ {{_, _}, 0}, __]-> 0}];
	masses= Association@@ (#-> Simplify@ MassMatrix[#, expr, Global`d$$1, Global`d$$2]&)/@ fields;
	
	Rule[#, 0]&/@ Association@@ Flatten@ KeyValueMap[
		Table[
			If[#2[[i,j]]=!= 0, {#1[[i]], #1[[j]]}-> #2[[i,j]], Nothing]
		, {j, Length@ #1}, {i, j- 1}]&
	, masses]
]


(* ::Text:: *)
(*Mass matrix*)


MassMatrix[fields_, lag_, i_, j_]:= Module[{deltas, f1, f2, fs, fBars, gaugeInds, fieldType, indNo},
	(*For projecting over open gauge indices*)
	gaugeInds= DeleteCases[GetFields[First@ fields, Indices], 
		_? (!MemberQ[Keys@ $GaugeGroups, GroupFromRep@ #]&)];
	deltas= Times@@(Delta[Index[i, #],Index[j, #]]&/@ gaugeInds);
	deltas= deltas/ Power[deltas, 2];
	
	(*Fields used in taking functional derivatives*)
	fieldType= GetFields[First@ fields, Type];
	{fBars, fs}= Transpose@ Table[
			indNo= Length@ GetFields[f1, Indices]+ If[fieldType === Vector, 1, 0];
			{Bar[f1@@ ConstantArray[i, indNo]], f1@@ ConstantArray[j, indNo]}
		, {f1, fields}];
	
	Map[Contract, Table[VarD[RelabelIndices[lag, Unique->True], f1, f2]/. (_Field| _FieldStrength)-> 0
		, {f1, fBars}, {f2, fs}] deltas If[fieldType === Vector, Metric[i, j]/\[ScriptD], 1], {2}]
]


(* ::Subsubsection::Closed:: *)
(*Kinetic mixing of non-GB scalars with vectors*)


(* ::Text:: *)
(*Finds all cases of kinetic mixing between heavy vectors and non-GB scalars and determine what coefficients are put to zero*)


ExtractScalarVectorMixing[lag_, expr_]:= Module[{vectors, terms, nonGBscalars},
	vectors= DeleteDuplicates@ Cases[lag, Field[f_, _Vector, __]:> f, All];
	(*Keep all scalar-vector mixing terms*)
	(*terms= expr/. (AtomicOp| CompOp)[{Except@ {{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, _}, _]-> 0;*)
	terms= expr/. {AtomicOp[Except@ {{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, __]-> 0,
		CompOp[_, Except@ {{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, __]-> 0};
	If[terms === 0, Return@ <||>;];
	(*Remove the Goldstone terms*)
	nonGBscalars= Alternatives@@ GetFieldsByProperty[Type-> Scalar, GoldstoneBoson-> False];
	terms= If[Head@ terms === Plus, List@@ terms, {terms}];

	terms= terms/. Times[AtomicOp[{f_, _}, __], rest___]:> 
		If[FreeQ[f, nonGBscalars], Nothing, Rule[f, Times@ rest]];
	Rule[#, 0]&/@ Association@@ terms	
]


(* ::Subsubsection::Closed:: *)
(*Find non-canonical kinetic terms*)


(* ::Text:: *)
(*Finds all kinetic terms and determine the coefficients required for canonical normalization and elimination of all off-diagonal terms*)


ExtractKineticTermConditions[lag_]:= Module[{inconsistencies, terms, var},
	(*Might have to use simplified Lagrangian*)
	terms= Collect[Operator@ lag, _Operator];
	If[terms === 0, Return@ <||>;];
	terms= If[Head@ terms === Plus, List@@ terms, {terms}];
	
	(*Finds all instances of kinetic terms*)
	terms= Cases[terms, op_Operator|Times[op_Operator, rest__]/; KineticOpQ@ op || OtherScalarKinQ@ op:> {op, Times@ rest}];
	(*Set coefficients equal to their canonical values*)
	terms= Rule[Cases[#1, Field[f_, __]|FieldStrength[f_, __]:> f, All], 
		#2 == KineticTermNormalization@ #1]&@@@ terms;
	
	(*Check for inconsistencies*)
	inconsistencies= DeleteCases[terms, Except@ Rule[_, False]];
	If[Length@ inconsistencies> 0,
		KeyValueMap[(Message[PrepareLagrangianForMatching::incKinTerm, #1, First@ #2]&), inconsistencies];
		Abort[];
	];
	
	(*Solve non-trivial conditions*)
	terms= DeleteCases[terms, HoldPattern@ Rule[_, True]];
	terms[[;;, 2]]= Table[
			SmartSolve[cond]
		, {cond, terms[[;;, 2]]}];
	
	Association@@ terms
]


(* ::Text:: *)
(*Additional kinetic forms, which are not standard for the IBPSimplifications*)


OtherScalarKinQ= MatchQ[Alternatives[
		HoldPattern@ Operator[Bar@Field[_, Scalar, _, {\[Mu]_}], Field[_, Scalar, _, {\[Mu]_}]],
		HoldPattern@ Operator[Field[_, Scalar, _, {\[Mu]_}], Field[_, Scalar, _, {\[Mu]_}]],
		HoldPattern@ Operator[Field[_, Scalar, _, {}], EoM@ Bar@ Field[_, Scalar, _, {}]]
	] ];  


(* ::Text:: *)
(*Standard normalization for kinetic terms *)


KineticTermNormalization@ op_:= op/.{
		(*Scalars*)
		HoldPattern@ Operator[Bar@Field[f_, Scalar, _, {\[Mu]_}], Field[f_, Scalar, _, {\[Mu]_}]]-> 1,
		HoldPattern@ Operator[Field[f_, Scalar, _, {\[Mu]_}], Field[f_, Scalar, _, {\[Mu]_}]]-> 1/2,
		HoldPattern@ Operator[Bar@Field[f_, Scalar, _, {}], EoM@ Field[f_, Scalar, _, {}]]-> -1/2,
		HoldPattern@ Operator[Field[f_, Scalar, _, {}], EoM@ Field[f_, Scalar, _, {}]]-> -1/2,
		HoldPattern@ Operator[Field[f_, Scalar, _, {}], EoM@ Bar@ Field[f_, Scalar, _, {}]]-> -1/2,
		(*Fermions*)
		HoldPattern@ Operator[Bar@ Field[f_, Fermion, _, {}]** EoM@ Field[f_, Fermion, _, {}] ]-> 1,
		HoldPattern@ Operator[Bar@ Field[f_, Fermion, _, {}]** DiracProduct@ _Proj**
			EoM@ Field[f_, Fermion, _, {}] ]-> 1,
		HoldPattern@ Operator[Transp@ Field[_, Fermion, _, {}]** DiracProduct[GammaCC]** 
			EoM@ Field[_, Fermion, _, {}] ]-> 1/2,
		HoldPattern@ Operator[Transp@ Field[_, Fermion, _, {}]** DiracProduct[GammaCC, _Proj]**
			EoM@ Field[_, Fermion, _, {}] ]-> 1,
		(*Vectors*)
		HoldPattern@ Operator[FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}],
			FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}]] -> -1/4,
		HoldPattern@ Operator[Bar@ FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}],
			FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}]] -> -1/2,
		(*Off-diagonals*)
		_Operator-> 0
		}; 


(* ::Text:: *)
(*A custom solve that attempts to find a good variable to solve in, by looking for an equality of the form *)
(*	0 = a + b x^n *)


SmartSolve@ equality:Equal[lhs_, rhs_]:= Module[
		{candidates, eq, expr, powers, var, vars, temp, tempExpr, terms},
	expr= lhs- rhs;
	terms= TermsToList@ expr;
	(*Consider adjusting for the posibility of coliding symbol and coupling names*)
	vars= DeleteDuplicates@ Cases[equality, _Coupling, All];
	tempExpr= equality/. {_Coupling-> 1, _Delta-> 1};
	
	(*TEMPORARY*)
	If[Length@ vars === 0, 
		Return@ Nothing
	];
	
	(*Check if expression is of the form a+ b var^n = 0 for some var*)
	candidates= Table[Catch[
			If[!PolynomialQ[expr, var], Throw[Nothing]];
			(*What powers of the variable is present in the terms*)
			powers= DeleteDuplicates@ Cases[terms, Power[var, n_]-> n, All];
			Switch[Length@ powers
			,0,
				{var, 1}
			,1,
				If[FreeQ[expr/. Power[var, _]-> 1, var],
					{var, First@ powers}
				,
					Nothing
				]
			,_,
				Nothing
			]
		], {var, vars}];
	
	(*If not candidate found*)
	If[Length@ candidates === 0, 
		(*Default variable: Improve this*)
		Return@ First@ Solve[equality, First@ var]
	];
	
	(*Solve in the smalest power of the symbols*)
	var= First@ MinimalBy[candidates, Last];
	If[Last@ var > 1,
		eq= equality/. Power@@ var-> temp;
		First@ Solve[eq, temp]/. temp:> Power@@ var
	,
		First@ Solve[equality, First@ var]
	]
]


(* ::Subsection:: *)
(*Construct gauge-fixing*)


(* ::Subsubsection::Closed:: *)
(*Goldstone mixing terms *)


GBMixingTerms@ expr_:= Module[{vectors, terms},
	vectors= DeleteDuplicates@ Cases[expr, Field[f_, _Vector, __]:> f, All];
	terms= InternalSimplify[expr, InternalOpRepresentation-> True];
	terms= terms/. (AtomicOp| CompOp)[{Except@{{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, _}, _]-> 0;
	AtomicToNormalForm@ terms
]


(* ::Subsubsection::Closed:: *)
(*Find Goldstone bosons and decay constants*)


(* ::Text:: *)
(*We determine the GBs corresponding to the heavy vectors in the broken phase along with their associated decay constants.*)


DetermineGBs[lag_]:= Module[{terms, vectors},
	terms= GBMixingTerms@ lag;
	vectors= DeleteDuplicates@ Cases[terms, Field[f_, _Vector, __]:> f, All];
	Association@@ (#-> GBAndDecayConstant[#, terms]&)/@ vectors
];


GBAndDecayConstant::mult= "Multiple potential GB candidates for `1`";


GBAndDecayConstant[vector_, lag_]:= Module[{a, decayConst, gb, inds, term},
	inds= ConstantArray[a, Length@ GetFields[vector, Indices] +1];
	term= VarD[RelabelIndices[lag, Unique->True], vector@@ inds];
	term= Collect[term, _Field];
	
	(* Check if multiple potential GB candidates to the Vector *)
	If[Count[term, _Field, All] > 1,
		Message[GBAndDecayConstant::mult, vector];
		Abort[];
	];
	
	(*Determine GB*)
	gb= FirstCase[term, Bar@ Field[f_, __]-> Conj@ f, 
		FirstCase[term, Field[f_, __]-> f, All], All];

	(*Determine decay constant*)
	decayConst= term/. _Field-> 1;

	{gb, decayConst}
];


(* ::Subsubsection::Closed:: *)
(*Construct gauge-fixing terms*)


(* ::Text:: *)
(*Construct a gauge-fixing term for the partially gauge-fixed G/H coset*)


GaugeFixingTerms[lag_, gaugeParam_]:= Module[{gbAndDecayConst},
	(*Determine GBs and decay constants*)
	gbAndDecayConst= DetermineGBs[lag];

	(*Construct gauge-fixing terms*)
	-1/gaugeParam Plus@@ KeyValueMap[ConstructGaugeFixingTerm[gaugeParam], gbAndDecayConst]// RelabelIndices
]


(* ::Text:: *)
(*Construct gauge-fixing term for a single heavy vector*)


ConstructGaugeFixingTerm[gaugeParam_][vec_, {gb_, fDecay_}]:= Module[{scalar, a, inds, mu1, mu2},
	inds= Sequence@@ ConstantArray[a, Length@ GetFields[vec, Indices]];
	scalar= If[Head@ gb === Conj, Bar@ First[gb][inds], gb@ inds];

	If[GetFields[vec, SelfConjugate], 
		1/2 (CD[mu1, vec[mu1, inds]]- gaugeParam fDecay scalar) (CD[mu2, vec[mu2, inds]] - gaugeParam fDecay scalar)
	,
		(CD[mu1, vec[mu1, inds]] - gaugeParam fDecay Bar@ scalar) (Bar@ CD[mu2, vec[mu2, inds]] - gaugeParam fDecay scalar)
	]
]


(* ::Subsubsection::Closed:: *)
(*Initialize ghost fields*)


(* ::Text:: *)
(*Function for initializing the ghost fields associated with heavy vector fields. Returns the label of the ghost and anti ghost field*)


InitializeGhostFields@ vector_:= Module[{antiGhostLabel, ghostLabel},
	ghostLabel= Symbol["$c"<> ToString@ vector];
	antiGhostLabel= Symbol["$c"<> ToString@ vector <> "bar"];
	DefineField[ghostLabel, Ghost, 
		Charges-> $FieldAssociation[vector, Charges],
		Indices-> $FieldAssociation[vector, Indices],
		Mass-> {Heavy, $FieldAssociation[vector, Mass]},
		SelfConjugate-> $FieldAssociation[vector, SelfConjugate] ]; 
	DefineField[antiGhostLabel, Ghost, 
		Charges-> MapAt[Minus, $FieldAssociation[vector, Charges], {All, 1}],
		Indices-> Bar/@ $FieldAssociation[vector, Indices],
		Mass-> {Heavy, $FieldAssociation[vector, Mass]},
		SelfConjugate-> $FieldAssociation[vector, SelfConjugate] ]; 
	{antiGhostLabel, ghostLabel}
]


(* ::Subsubsection::Closed:: *)
(*Construct Ghost terms*)


ConstructGhostTerms[lag_, gaugeParam_]:= Module[{fTensor3H, fTensor2H1L, gbs, ghosts, heavyMetric, inds, lightMetric, vec, mu, c, cbar, V, i, j, k, l, alpha},
	gbs= DetermineGBs@ lag; 
	ghosts= AssociationMap[InitializeGhostFields, Keys@ gbs];
	(*anti-ghosts*)
	cbar= a|-> Evaluate[Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ GetFields[vec, Indices]];
			If[$FieldAssociation[vec, SelfConjugate], 
				{ghosts[vec][[1]][inds]}, {ghosts[vec][[1]][inds], Bar@ ghosts[vec][[1]][inds]}
			]
		, {vec, HeavyVectors[]}] ];
	(*ghosts*)
	c= a|-> Evaluate[Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ GetFields[vec, Indices]];
			If[$FieldAssociation[vec, SelfConjugate], 
				{ghosts[vec][[2]][inds]}, {ghosts[vec][[2]][inds], Bar@ ghosts[vec][[2]][inds]}
			]
		, {vec, HeavyVectors[]}] ];
	(*heavy vectors*)
	V= a|-> Evaluate[Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ GetFields[vec, Indices]];
			If[$FieldAssociation[vec, SelfConjugate], 
				{vec[mu, inds]}, {vec[mu, inds], Bar@ vec[mu, inds]}
			]
		, {vec, HeavyVectors[]}] ];
		
	mass2Matrix= DiagonalMatrix@ Flatten@ Table[
			If[$FieldAssociation[vec, SelfConjugate], {gbs[vec][[2]]^2}, {gbs[vec][[2]]^2, gbs[vec][[2]]^2}]
		, {vec, HeavyVectors[]}];
	fTensor3H= StructureConstant3H[];
	fTensor2H1L= StructureConstant2H1L[];
	heavyMetric= HeavyIndexMetric[];
	lightMetric= LightIndexMetric[];
	
	
	Inner[NonCommutativeMultiply, cbar[i],
		-CD[{mu, mu}, c[i]] -
		mass2Matrix . c[i] -
		2 TensorContract[TensorProduct[fTensor3H[i, k, j], CD[mu, V@ k]], {2, 4}] . c[j] -
		TensorContract[TensorProduct[fTensor3H[i, k, j], V@ k], {2, 4}] . CD[mu, c[j]] +
		TensorContract[TensorProduct[fTensor2H1L[i, alpha, k], fTensor2H1L[l, alpha, j], lightMetric, 
			V@ k, heavyMetric . V@ l], {{2, 7}, {5, 8}, {3, 9}, {4, 10}}] . c[j]  +
		2 I gaugeParam heavyMetric . Fxphi[i, j] . c[j] -
		I gaugeParam heavyMetric . Transpose[Fxphi[j, i]] . c[j]
	]// Contract// RelabelIndices
]


(* ::Text:: *)
(*Auxiliary function for the coupling of scalars with ghosts: f_i^a x^a_{jb} \[CurlyPhi]^b*)


Fxphi[i_, j_]:= Module[{a, b, fields},
	fields= SymmetryBreakingScalars[][[;;1]];
	Map[ContractCGs@* Contract@* Expand, Sum[
		If[GetFields[f, SelfConjugate],
			Outer[Times, DecayConstantMatrix[f][i, a], 
				BrokenGeneratorsOnField[f][j, a, b] ScalarFields[][b][f]]
		,
			Outer[Times, DecayConstantMatrix[Conj@f][i, a], 
				BrokenGeneratorsOnField[f][j, a, b] ScalarFields[][b][f]]+
			Outer[Times, DecayConstantMatrix[f][i, a], 
				BrokenGeneratorsOnField[Conj@ f][j, a, b] Bar@ ScalarFields[][b][f]]/. Bar@ d_Delta-> d
		]
	, {f, fields}], {2}]
	
]
