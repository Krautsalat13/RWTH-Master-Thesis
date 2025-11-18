(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`LoopIntegration`*)


(* ::Subtitle:: *)
(*Procedures for evaluating one-loop integrals and CDE*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["\[Mu]bar2"]
PackageExport["hbar"]
PackageExport["\[Epsilon]"]
PackageExport["LF"]


PackageExport["EvaluateLoopFunctions"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["Prop"]
PackageScope["InvProp"]
PackageScope["LoopMom"]
PackageScope["MomDev"]


PackageScope["LoopIntegrate"]
PackageScope["LFFull"]
PackageScope["LFFull2LF"]
PackageScope["CanonizeLoopFunctions"]
PackageScope["ReduceLoopFunctions"]
PackageScope["SimplifyMassFunction"]


PackageScope["GatherLoopMomenta"]
PackageScope["SymmetricLorentzInds"]
PackageScope["EvaluateSymmetricLorentzInds"]


PackageScope["EpsExpand"]
PackageScope["PerformMomDerivatives"]


PackageScope["FuncNCM"]


PackageScope["MomNCM"]
PackageScope["MatrixNCM"]
PackageScope["MomAntiCommutator"]
PackageScope["ReleaseMomNCM"]
PackageScope["ExpandMatrixNCM"]


PackageScope["ev"]


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


\[Mu]bar2::usage=
"\!\(\*SuperscriptBox[\(\[Mu]\), \(2\)]\) is the matching scale squared.";
hbar::usage=
"\[HBar] is the power-counting parameter for the loop-order in an expression. It is understood as a \!\(\*FractionBox[\(1\), \(16 \*SuperscriptBox[\(\[Pi]\), \(2\)]\)]\) factor in the final result.";
\[Epsilon]::usage=
"\[Epsilon] is the parameter defined by expanding the space-time dimension \[ScriptD] around four dimensions, i.e. \[ScriptD]=4-2\[Epsilon].";
LF::usage=
"LF[{\!\(\*SubscriptBox[\(m\), \(1\)]\),...,\!\(\*SubscriptBox[\(m\), \(n\)]\)},{\!\(\*SubscriptBox[\(i\), \(1\)]\),..., \!\(\*SubscriptBox[\(i\), \(n\)]\), \!\(\*SubscriptBox[\(i\), \(n + 1\)]\)}] is a placeholder for the finite piece of the loop integration performed over \!\(\*FractionBox[\(1\), \(\(\*SuperscriptBox[\((\*SuperscriptBox[\(k\), \(2\)] - \*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)])\), SubscriptBox[\(i\), \(1\)]] ... \) \*SuperscriptBox[\((\*SuperscriptBox[\(k\), \(2\)] - \*SuperscriptBox[SubscriptBox[\(m\), \(n\)], \(2\)])\), SubscriptBox[\(i\), \(n\)]] \*SuperscriptBox[\((\*SuperscriptBox[\(k\), \(2\)])\), SubscriptBox[\(i\), \(n + 1\)]]\)]\) where k is the loop momentum.";


EvaluateLoopFunctions::usage=
"EvaluateLoopFunctions[expr] replaces the placeholders LF in expr with the appropriate finite pieces. ";


(* ::Subsubsection::Closed:: *)
(*Internal*)


InvProp::usage=
"InvProp[m] is the inverse propagator: \!\(\*SuperscriptBox[\(k\), \(2\)]\)-\!\(\*SuperscriptBox[\(m\), \(2\)]\), where k is the loop momentum.";
LoopMom::usage=
"LoopMom, or LoopMom[Index[\[Mu],Lorentz]], is the momentum running inside the loop.";
MomDev::usage=
"MomDev[\[Mu]] is an object representing a derivative with respect to the loop momentum with Lorentz index \[Mu], acting to the right.";


LoopIntegrate::usage=
"LoopIntegrate[expr, LogTerm ->True/False] performs the scalar integral over the loop momentum propagators in the expression, written in the form: InvProp[m1\!\(\*SuperscriptBox[\(]\), \(-i1\)]\)...InvProp[mn\!\(\*SuperscriptBox[\(]\), \(-in\)]\)InvProp[0\!\(\*SuperscriptBox[\(]\), \(j\)]\);
By default, LogTerm is False. If turned on (LogTerm-> True), there is an additional factor corresponding to the integral over \[Xi] from the integral expression of the log."
LFFull::usage="LFFull[{\!\(\*SubscriptBox[\(m\), \(1\)]\),...,\!\(\*SubscriptBox[\(m\), \(n\)]\)},{\!\(\*SubscriptBox[\(i\), \(1\)]\),..., \!\(\*SubscriptBox[\(i\), \(n\)]\), \!\(\*SubscriptBox[\(i\), \(n + 1\)]\)}] is a placeholder for the loop integration (including divergent pieces) performed over \!\(\*FractionBox[\(1\), \(\(\*SuperscriptBox[\((\*SuperscriptBox[\(k\), \(2\)] - \*SuperscriptBox[SubscriptBox[\(m\), \(1\)], \(2\)])\), SubscriptBox[\(i\), \(1\)]] ... \) \*SuperscriptBox[\((\*SuperscriptBox[\(k\), \(2\)] - \*SuperscriptBox[SubscriptBox[\(m\), \(n\)], \(2\)])\), SubscriptBox[\(i\), \(n\)]] \*SuperscriptBox[\((\*SuperscriptBox[\(k\), \(2\)])\), SubscriptBox[\(i\), \(n + 1\)]]\)]\) where k is the loop momentum."
LFFull2LF::usage="LFFull2LF takes an expression with LFFull, and separates the loop function into divergent and finite pieces, the latter given in terms of LF."


ExtractMomenta::usage=
"Extract slashed LoopMom from DiracProduct by replace product of momentum vectors with symmetrized product of metrics,
e.g. LoopMom[Index[mu, Lorentz]] LoopMom[Index[nu, Lorentz]] -> Prop[0]^-1 Metric[Index[mu, Lorentz], Index[nu, Lorentz]] /\[ScriptD], etc.";
EpsExpand::usage=
"EpsExpand[expr,Order->n] expands an expression to nth order (by default 0) in \[Epsilon], the regulator from dimensional regularization, after replacing the space-time dimension \[ScriptD] by 4-2\[Epsilon]."
PerformMomDerivatives::usage=
"Perform all momentum derivatives in a loop integral expression using IbP to reduce the number of terms. ";


TermsToList::usage=
"Transform a sum of terms into a list, or convert a single term into a list."


MomNCM::usage=
"Internal NCM head for use in the CDE of the SuperTraces."
MatrixNCM::usage=
"Temporary NCM for use in X-term multiplication in SuperTraces (non-distributive)."
MomAntiCommutator::usage=
"Anti commutator used in the CDE expansion."
ReleaseMomNCM::usage=
"Drop the MomNCM head after tracing."
ExpandMatrixNCM::usage=
"Expand out the MatrixNCM after X-term substitutions."


ev::usage=
"Dummy label for evanescent operators."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Momentum integrals*)


(* ::Subsection:: *)
(*Properties*)


(* ::Text:: *)
(*Proporties of loop labels 'hbar' and 'ev,'  used to account for ordering of terms *)


(* ::Subsubsection::Closed:: *)
(*hbar *)


(* ::Text:: *)
(*Dummy label counting loop order*)


Bar@hbar ^:=hbar
Transp@hbar ^:=hbar


Bar@\[Mu]bar2 ^:=\[Mu]bar2
Transp@\[Mu]bar2 ^:=\[Mu]bar2


(* ::Text:: *)
(*Limit to one - loop order*)


hbar^x_/;x>=2 ^:= 0


(* ::Subsubsection::Closed:: *)
(*ev*)


(* ::Text:: *)
(*Dummy label for evanescent operators*)


ev^x_/;x>=2 ^:= ev;


(* ::Subsection:: *)
(*Epsilon expansion*)


(* ::Text:: *)
(*Function for expanding an expression around \[Epsilon] = 0, keeping poles and terms up to power 'order.'*)


Options@ EpsExpand= {Order-> 0};


EpsExpand[expr_, OptionsPattern[]]:= Module[{term},
	Sum[
		Normal@ Series[term/. SymGammaFactor:> EvaluateGammaFactor/. \[ScriptD]-> 4- 2\[Epsilon], {\[Epsilon], 0, OptionValue@Order}]
	, {term, TermsToList@ expr}]
]


(* ::Subsection:: *)
(*Scalar integrals*)


(* ::Subsubsection::Closed:: *)
(*General loop integration function*)


(* ::Text:: *)
(*Performs the 1-loop integral over an expression. Assumes that all loop-momenta occurs in 'Prop' objects. *)


LoopIntegrate::onescale="There is more than one scale in the LogTerm integration."


Options@ LoopIntegrate= {LogTerm-> False};


(* 
TO DO:
apply HcSimplify first, then integralType & ToLoopFunctions, and eventually use PlusHc, to ensure manifest hermiticity 
*)
LoopIntegrate[expr_, opt:OptionsPattern[]]:=Module[{out},
	out= Expand[expr* integralType[<||>, 0]];
	If[OptionValue@ LogTerm,
		If[!AllTrue[Cases[out, integralType[x_, _]:> Length@ x, Infinity], # < 2&],
			Message[LoopIntegrate::onescale];
			Abort[]
		];
		out/. integralType[<|m_-> \[Beta]_|>, \[Alpha]_]:> integralType[{m, \[Beta]}, \[Alpha], LogTerm-> True]/.
			x_integralType:> SingleScaleIntegral@@ x
	,
		out/. x_integralType:> ToLoopFunctions@@ x
	]
];


(* ::Subsubsection::Closed:: *)
(*Auxiliary functions*)


(*GatherPropagatorsInTerm@ term_:= Module[{props},
	props= If[Head@ term === Times, List@@ term, {term}];
	props= Proppagators@@ Cases[temp, Prop@ m_|Power[Prop@ m_, p_.]:> {m, Times@ p}];
	props* (term/. _Prop-> 1)
]*)


(* ::Text:: *)
(*Identifying propagators in the expressions:*)
(*mprops is an association of mass propagators and their corresponding power, for example a denominator 1/((p^2-m1^2)^n1(p^2-m1^2)^n2) reads <|m1->n1,m2->n2|> in mprops.*)
(*\[Alpha] is the power of (1/p^2)^\[Alpha].*)


integralType/: integralType[mprops_, \[Alpha]_] Prop@ 0:= integralType[mprops, \[Alpha]+ 1];
integralType/: integralType[mprops_, \[Alpha]_] Prop@ m_:= integralType[Merge[{mprops, <|m-> 1|>}, Total], \[Alpha]];
integralType/: integralType[mprops_, \[Alpha]_] Power[Prop@ 0, n_]:= integralType[mprops, \[Alpha]+ n];
integralType/: integralType[mprops_, \[Alpha]_] Power[Prop@ m_, n_]:= integralType[Merge[{mprops, <|m-> n|>}, Total], \[Alpha]];


(* ::Subsubsection::Closed:: *)
(*Single scale (and Log term) integrals*)


(* ::Text:: *)
(*For the log term integral (LogTerm->True) there is an additional term due to the \[Xi] integral.*)
(*SingleScaleIntegral can either evaluate the pole (Pole ->True) or the finite part of the integral (Pole ->False). In the latter case, expansion with respect to \[Epsilon] is already performed.*)


LoopIntegration::error="The LogTerm option might have been turned on by mistake."


SingleScaleIntegral[{mass_, \[Beta]_Integer}, \[Alpha]_Integer, OptionsPattern[{LogTerm->False, Pole->False}]]:=
Module[{preFact},
	If[OptionValue[LogTerm],
		If[\[Beta] <= 0 || \[Alpha] >= 2 || \[Alpha]+ \[Beta] <= 2,
			Message[LoopIntegration::error];
			Abort[];
		];
		Return[(-1)^(\[Alpha]+\[Beta]) I mass^(2(2- \[Alpha]- \[Beta]))/ Gamma@ \[Beta] Gamma[2- \[Alpha]] Gamma[\[Alpha]+ \[Beta]- 2] (1+ (Log[\[Mu]bar2/ mass^2]+ 1+ PolyGamma[\[Alpha]+\[Beta]-2]- PolyGamma[2- \[Alpha]])\[Epsilon] )]
	];

	If[\[Beta] <= 0, Return@ 0;];
	preFact= (-1)^(\[Alpha]+\[Beta]) I mass^(2(2- \[Alpha]- \[Beta]))/ Gamma@ \[Beta];

	If[\[Alpha] >= 2, (* Expand around IR pole*)
		Return[-preFact Gamma[\[Alpha]+ \[Beta]- 2] (-1)^\[Alpha]/ (\[Alpha]-2)! *
			If[OptionValue[Pole],
				1/\[Epsilon]
			,
				(-HarmonicNumber[\[Alpha]-2]+ Log[\[Mu]bar2/ mass^2]+ 1+ PolyGamma[\[Alpha]+\[Beta]-2]+ EulerGamma)
			]
		];
	];
	If[\[Alpha]+ \[Beta] <= 2, (* Expand around UV pole*)
		Return[preFact Gamma[2- \[Alpha]] (-1)^(\[Alpha]+\[Beta])/(2-\[Alpha]-\[Beta])! *
			If[OptionValue[Pole],
				1/\[Epsilon]
			,
				(+HarmonicNumber[2-\[Alpha]-\[Beta]]+ Log[\[Mu]bar2/ mass^2]+ 1- PolyGamma[2- \[Alpha]]- EulerGamma)
			]
		];
	];
	If[OptionValue[Pole],
		0
	,
		preFact Gamma[2- \[Alpha]]Gamma[\[Alpha]+ \[Beta]- 2]
	]
];


(* ::Subsubsection::Closed:: *)
(*Multiscale integrals*)


MultiScaleIntegral[propPowers_Association, \[Alpha]_, OptionsPattern[{Pole->False}]]:= Module[{mi, pi, mj, temp, isFinite, firstMass, out},
	isFinite = Plus@@ propPowers + \[Alpha]> 2 && \[Alpha]<2;
	If[OptionValue[Pole],If[isFinite, Return@ 0;]];
	out= Sum[temp= KeyDrop[propPowers, mi];
		SingleScaleIntegral[{mi, propPowers@ mi - pi}, \[Alpha], Pole -> OptionValue[Pole]] /pi! D[Product[Power[M2@ mi -M2@ mj, -propPowers@ mj],
			{mj, Keys@ temp}], {M2@ mi, pi}]
	,{mi, Keys@ propPowers}, {pi, 0, propPowers@ mi -1}]/. M2@ mass_:> mass^2;
	If[isFinite, firstMass = First@Sort@Keys@propPowers; out /. \[Mu]bar2 -> firstMass^2, out]
];


(* ::Subsubsection::Closed:: *)
(*Loop functions *)


(* ::Text:: *)
(*Removes duplicate masses in LF and sum over its powers.*)


LF[m_List?(!DuplicateFreeQ@# &),ind_List]:= Module[{dupm,posm,newind,newm},
	dupm=Select[Tally@m,#[[2]]>1 &][[1,1]];
	posm=First/@Position[m,dupm];
	newind=DeleteCases[(ReplacePart[Append[Association@@(#->xxx&/@posm[[2;;]]),First@posm->Total@(Part[ind,#]&/@posm)]][ind]),xxx];
	newm=DeleteCases[ReplacePart[Association@@(#->xxx&/@posm[[2;;]])][m],xxx];
	LF[newm,newind]
]


(* ::Text:: *)
(*Canonically order the loop function*)


CanonizeLoopFunctions@ expr_:= expr/. lf:LF[{_,__}, _]:> CanonizeLoopFunctions@ lf;


CanonizeLoopFunctions@ LF[masses_, subscripts_]:= Module[{props},
	props= Transpose@ {masses, Most@ subscripts};
	props= Transpose@ SortBy[props, {Minus@* Last, First}];
	LF[First@ props, Append[Last@ props, Last@ subscripts]]
];


(* ::Text:: *)
(*Simplification rules for the loop functions*)


(* Canonical order *) (* this ordering is not canonical once flavor indices are considered on masses *)
LFFull[m_?(!OrderedQ[#]&), ord_]:= LFFull[Sort[m],Append[ord[[;;-2]][[Ordering[m]]],ord[[-1]]]]

(* Remove zero entries (other than the last one) *)
LFFull[m_,ord_?(MemberQ[#[[;;-2]],0]&)]:=LFFull[Delete[m,Position[ord[[;;-2]],0]],Append[DeleteCases[ord[[;;-2]],0],ord[[-1]]]]

(* Set to zero scaleless integrals *)
LFFull[{},ord_]:=0


(* ::Text:: *)
(*The  reduction  identities  below  are  problematic  as  they  can  significantly  complicate  coefficients .*)
(* It  is  also  not  guaranteed  that  a  term  and  its  hermitian  conjugate  are  reduced  in  the  same  way .*)


(*
(* IBP rule to make the first power one *)
LFFull[m_,ord_?(#[[1]]>1 &)]:=(\[ScriptD]/2-ord[[-1]]-1)/(ord[[1]]-1)LFFull[m,ord+UnitVector[Length@ord,Length@ord]-UnitVector[Length@ord,1]]-Sum[ord[[i]]/(ord[[1]]-1)LFFull[m,ord+UnitVector[Length@ord,i]-UnitVector[Length@ord,1]],{i,2,Length@ord-1}]

(* Mass relation to reduce massless propagators (only used after first power is one) *)
LFFull[m_,ord_?(#[[1]]==1 && #[[-1]]>0 &)]:= - 1/m[[1]]^2 LFFull[m[[2;;]],ord[[2;;]]] + 1/m[[1]]^2 LFFull[m,ord - UnitVector[Length@ord,Length@ord]]

(* Mass relation to reduce positive powers of mass terms multiplying loop functions *)
LFFull/: Power[mi_, 2] LFFull[m:{m1___,mi_,mn___},ord_]:= LFFull[{m1,mi,mn},ord-UnitVector[Length@ord,Length@ord]] - LFFull[{m1,mi,mn},ord-UnitVector[Length@ord,Position[m,mi][[1,1]]]]
*)


(* ::Text:: *)
(*A factor "i" is extracted from the definition of the loop function, so that it is real*)


ToLoopFunctions[propPowers_Association, \[Alpha]_]:= Module[{props,full,pole},
	props= SortBy[KeyValueMap[{#1, #2}&, propPowers], (-#[[2]]&)];
	(*full= I LFFullEntry[props[[;;, 1]], props[[;;, 2]] ~ Join ~ {\[Alpha]}];*)
	full= I LFFull[props[[;;, 1]], props[[;;, 2]] ~ Join ~ {\[Alpha]}];
	LFFull2LF[full]
]


LFFull2LF[full_]:= Normal@Series[Expand@full/. {\[ScriptD]->4-2\[Epsilon], LFFull[args__]:> LF[args] + EvaluateLoopFunctions[LF[args], Pole -> True]},{\[Epsilon],0,0}]


EvaluateLoopFunctions[LF[denoms_, powers_], opt:OptionsPattern[{Pole->False}]]:= 
EvaluateLoopFunctions[LF[denoms, powers], opt]= Module[{association},
	association= Association@@ (#[[1]]->#[[2]]&)/@ Transpose@ {denoms, powers[[;;-2]]};
	-I MultiScaleIntegral[association, Last@ powers, Pole-> OptionValue[Pole]]//Simplify
]

EvaluateLoopFunctions@ expr_:= Module[{res,summedInd},
	(* substitute loop functions *)
	res= expr/. lf_LF:> EvaluateLoopFunctions@ lf;
	(* expand *)
	res= TermsToList@BetterExpand[res];
	(* remove unnecessary FlavorSums *)
	res=Sum[
		If[FreeQ[term,_FlavorSum,All],
			term
			,
			summedInd= FirstCase[term,FlavorSum[ind_]:>ind,Nothing,All];
			If[FreeQ[term/._FlavorSum->1,summedInd,All],
				term/._FlavorSum->DimRep[Last@summedInd],
				term
			]
		]
		,
		{term,res}
	];
	res
]


(* ::Subsection:: *)
(*Numerator momentum structures*)


(* ::Text:: *)
(*Bring momenta to propagator (with zero mass) form.*)


LoopMom/: Power[_LoopMom, 2] := Power[Prop@ 0, -1];


(* ::Subsubsection::Closed:: *)
(*Extract momenta*)


GatherLoopMomenta@ expr_:= CollectMomenta@ ExtractMomenta@ expr


(* ::Text:: *)
(*Extract LoopMom-slashes from DiracProduct*)


ExtractMomenta@ expr_:= expr//. {
		DiracProduct[a___, GammaM@ LoopMom, b___] :>
			Module[{nu}, DiracProduct[a, GammaM@ Index[nu, Lorentz], b] LoopMom@ Index[nu, Lorentz]],
		DiracProduct[a___, Transp@ GammaM@ LoopMom, b___] :>
			Module[{nu}, DiracProduct[a, Transp@ GammaM@ Index[nu, Lorentz], b] LoopMom@ Index[nu, Lorentz]]
	};


(* ::Text:: *)
(*Replaces product of momentum vectors with symmetrized product of metrics, e.g. LoopMom[Index[mu, Lorentz]] LoopMom[Index[nu, Lorentz]] -> Prop[0]^-1 Metric[Index[mu, Lorentz], Index[nu, Lorentz]] /\[ScriptD], etc. *)


CollectMomenta@ expr_:= Module[{out= BetterExpand@ expr},
	If[Head@ out === Plus, Return[CollectMomenta/@ out]; ];
	LoopMoms@@ Cases[out, LoopMom@ l_:> l, All]* (out/. _LoopMom-> 1)
]


LoopMoms[]:= 1;
LoopMoms@ inds__:= If[EvenQ@ Length@ {inds}, Power[Prop@ 0, -Length@ {inds}/ 2] SymmetricLorentzInds@ inds, 0];


(* ::Subsubsection::Closed:: *)
(*Evaluating the symmetric tensor *)


(* ::Text:: *)
(*The symmetric Lorentz tensor coming from the product of loop momenta*)


EvaluateSymmetricLorentzInds@ expr_:=
	expr/. SymmetricLorentzInds-> SymmetricLorentzIndsReplacement// ContractMetric;


SetAttributes[SymmetricLorentzInds, Orderless]


SymmetricLorentzInds[a_, a_, rest___]:= SymmetricLorentzInds@ rest;


(*This method only works for non-repeating indices due to behavior of Permutations*)
SymmetricLorentzIndsReplacement@ lorentzInds___:= Module[{inds, n, symTensor},
	inds= List@ lorentzInds;
	n= Length@ inds/ 2;
	symTensor= If[OddQ@ Length@ inds, 0,
		Plus@@ Times@@@ Apply[Metric,
			DeleteDuplicatesBy[Partition[#, 2]&/@ Permutations@ inds, (Sort[Sort/@ #] &)], {2}]];
	(*SymGammaFactor[n] = Gamma[d/2]/ (Gamma[d/2 +n] 2^n)*)
	symTensor SymGammaFactor@ n
];


SymGammaFactor@ 0= 1;


EvaluateGammaFactor@ n_Integer:= EvaluateGammaFactor@ n= Normal@ Series[Gamma[2- \[Epsilon]]/(2^n Gamma[2- \[Epsilon]+ n]), {\[Epsilon], 0, 1}];


(* ::Subsection:: *)
(*Loop function reduction*)


(* ::Subsubsection::Closed:: *)
(*Simplification of loop function sum*)


(* ::Text:: *)
(*Function to simplify a sum of loop functions \[Dash] a mass function*)
(*The strategy is to reduce the propagator power associated with a given mass to 1, then try this for every possible mass to obtain the best result.*)


SimplifyMassFunction@ expr_Plus:= Module[{lfs, out= List@@ expr, keep, simp, temp, types, target},
	out= LFtoTempLF@ out;
	types= SortBy[FirstCase[#, TempLF[masses_, _]:> masses, {}, All]&/@ out// DeleteDuplicates, 
		Minus@* Length];
	Do[
		If[Length@ masses< 2, Continue[];];
		{keep, simp}= SelectAndDeleteCases[out, a_/; FreeQ[a, TempLF[masses, _]]];
		If[(target= Length@ simp) < 2, Continue[];];
		Do[
			temp= Plus@@ simp/. lf_TempLF:> PermuteTempLF[lf, {i, 1}];
			temp= temp// SimpTempLFRules// TermsToList;
			If[Length@ temp> target, Continue[]; ];
			target= Length@ temp;
			out= Join[keep, temp];
		, {i, Length@ masses}];
	, {masses, types}];
	Plus@@ out/. lf:TempLF[masses_, inds_]:> PermuteTempLF[lf, Ordering@ masses]/. TempLF-> LF
];
SimplifyMassFunction@ expr_:= expr;


(* ::Text:: *)
(*Replace LF with a temporary loop function*)


LFtoTempLF@ expr_:= expr/. LF-> TempLF /. lf:TempLF[masses_, inds_]:> PermuteTempLF[lf, Ordering@ masses];


PermuteTempLF[TempLF[masses_, inds_], order_]:= Module[{newMasses= masses, newInds= inds, old= Sort@ order},
	newMasses[[old]]= masses[[order]];
	newInds[[old]]= newInds[[order]];
	TempLF[newMasses, newInds]
]


(* ::Text:: *)
(*The IBP reduction rules are based on full loop functions with divergent pieces. By contrast LF is only the finite part. The pole part difference between LF and TempLF never multiplies O(\[Epsilon]) pieces before simplification and can be ignored. The O(\[Epsilon]) pieces picked up by IBP can combine with poles in TempLF to give finite contributions, which must be kept. After truncation of the O(\[Epsilon]) pieces the poles in TempLF can once again be ignored, when we convert back to LF. We know that the original LF expression was finite.*)


SimpTempLFRules@ expr_:= Block[{out},
	out= expr//. TempLF[masses_, ord:{first_/; first>1, mid__, last_}]:>
	(1 -\[Epsilon] -last)/(first -1) TempLF[masses, {first -1 , mid, last+1}] -
		Sum[{mid}[[i]]/(first-1) TempLF[masses, {first- 1, mid, last}+ UnitVector[Length@ ord, i+1]], {i, Length@ {mid}}];
	TempLFFiniteExtraction@ out
]


TempLFFiniteExtraction@ expr_:= Normal@ Series[
		Expand@ expr/. lf_TempLF:> lf + EvaluateLoopFunctions[LF@@ lf, Pole -> True]
	, {\[Epsilon], 0, 0}]/. Power[\[Epsilon], -1]-> 0;


(* ::Subsubsection::Closed:: *)
(*IBP reduction (not in use)*)


(* ::Text:: *)
(*For IBP reduction of the  loop functions *)


ReduceLoopFunctions@ expr_:= expr/. lf_LF:> ReduceLoopFunction@ lf;


ReduceLoopFunction@ lf_LF:= Module[{out},
	(*Possible divergent pieces would cancel in the end anyway*)
	out= LFFull@@ lf;
	out= out//. lfFullReductionRules;
	LFFull2LF@ out/. Power[\[Epsilon], -1]-> 1
]


lfFullReductionRules= {
	(* IBP rule to make the first power one *)
	LFFull[masses_, ord:{first_/; first>1, mid___, last_}]:>
		(2 -\[Epsilon] -last -1)/(first-1) LFFull[masses, {first- 1, mid, last+ 1}]
			-Sum[ord[[i]]/(first -1)LFFull[masses, {first- 1, mid, last}+ UnitVector[Length@ ord, i]], {i, 2, Length@ ord-1}],
	(* Mass relation to reduce massless propagators (only used after first power is one) *)
	LFFull[masses_,ord:{1, mid___, last_/; last>0}]:> 
		- 1/masses[[1]]^2 LFFull[masses[[2;;]], ord[[2;;]]] + 1/masses[[1]]^2 LFFull[masses, {1, mid, last -1}],
	(*Mass relation to increase massless propagators (only used after first power is one)*)
	LFFull[masses_, ord:{1, mid___, last_/; last<0}]:> 
		LFFull[masses[[2 ;;]], {mid, last + 1}] + masses[[1]]^2 LFFull[masses, {1, mid, last + 1}]
}


(* ::Section:: *)
(*Momentum structures for CDE*)


(* ::Subsection:: *)
(*Functional Non-commutative product *)


(* ::Text:: *)
(*A non-commutative product used in the evaluation on functional loops*)


(* ::Subsubsection::Closed:: *)
(*Commutative check*)


(* ::Text:: *)
(*Default assumption: everything is commutative*)


FuncCommuteQ[f_?FuncCommuteQ[x___]] := And@@ FuncCommuteQ/@ {x};
FuncCommuteQ[f_[x___]] := False;
FuncCommuteQ[_] := True;


(* ::Text:: *)
(*Define the non-commutative objects*)


FuncCommuteQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold] = False;


(* ::Text:: *)
(*Commutative is a head that can be used to temporarily treat an object as commutative*)


FuncCommuteQ@ Commutative@ _ ^= True;


(* ::Subsubsection::Closed:: *)
(*Non-commutative symbols used in Matchete declared here*)


(* ::Text:: *)
(*Define non-commutative objects*)


((FuncCommuteQ@ # ^= False) &) /@ {
	Field,
	FieldStrength,
	DiracProduct,
	OpenCD,
	Xop, Xterm,
	WilsonLine, WilsonTerm
	};


(* ::Subsubsection::Closed:: *)
(*Properties of FuncNCM*)


(* ::Text:: *)
(*Flatness*)


FuncNCM[a___, FuncNCM[b__], c___]:= FuncNCM[a, b, c];
FuncNCM[] = 1;
(*FuncNCM@ x_FuncNCM:= x;*)


(* ::Text:: *)
(*Extracting commuting objects*)


FuncNCM[a___, b_?FuncCommuteQ, c___]:= b * FuncNCM[a, c];
FuncNCM[a___, b_?FuncCommuteQ * x_, c___]:= b * FuncNCM[a, x, c];


(* ::Text:: *)
(*Distributivity*)


FuncNCM[a___, b_Plus, c___]:= FuncNCM[a, #, c] & /@ b


(* ::Subsection:: *)
(*Non-commutative products*)


(*To be deleted*)


(* ::Text:: *)
(*A specific non-commutative product for use in the CDE*)


(* ::Subsubsection::Closed:: *)
(*Commutative check*)


(* ::Text:: *)
(*Default assumption: everything is commutative*)


MomCommuteQ[f_?MomCommuteQ[x___]] := And@@ MomCommuteQ/@ {x};
MomCommuteQ[f_[x___]] := False;
MomCommuteQ[_] := True;


(* ::Text:: *)
(*Define the non-commutative objects*)


MomCommuteQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold] = False;


(* ::Text:: *)
(*Commutative is a head that can be used to temporarily treat an object as commutative*)


MomCommuteQ@ Commutative@ _ ^= True;


(* ::Subsubsection::Closed:: *)
(*Non-commutative symbols used in Matchete declared here*)


(* ::Text:: *)
(*Define non-commutative objects*)


((MomCommuteQ@ # ^= False) &) /@ {
	Field,
	FieldStrength,
	DiracProduct,
	MomDev,
	LoopMom,
	InvProp,
	Gop,
	Xop,
	HoldPart,
	MassOp
	};


(* ::Subsubsection::Closed:: *)
(*Properties of MomNCM*)


(* ::Text:: *)
(*Flatness*)


MomNCM[a___, MomNCM[b__], c___]:= MomNCM[a, b, c];
MomNCM[] = 1;
MomNCM@ x_MomNCM:= x;


(* ::Text:: *)
(*Extracting commuting objects*)


MomNCM[a___, b_?MomCommuteQ, c___]:= b * MomNCM[a, c];
MomNCM[a___, b_?MomCommuteQ * x_, c___]:= b * MomNCM[a, x, c];


(* ::Text:: *)
(*Distributivity*)


MomNCM[a___, b_Plus, c___]:= MomNCM[a, #, c] & /@ b


(* ::Subsubsection::Closed:: *)
(*Properties of Matrix NCM*)


(* ::Text:: *)
(*Flatness*)


MatrixNCM[a___, MatrixNCM[b__], c___]:= MatrixNCM[a, b, c];
MatrixNCM[] = 1;


(* ::Text:: *)
(*Working on lists and matrices*)


MatrixNCM@ expr_List:= expr;


MatrixNCM[a___, x_List, y_List, b___]:= MatrixNCM[a, Inner[MatrixNCM, x, y, Plus], b];
MatrixNCM[a___, x_List, y: Except[_List], b___]:= MatrixNCM[a, MatrixNCM[#, y]&/@ x, b];
MatrixNCM[a___, x: Except[_List], y_List, b___]:= MatrixNCM[a, MatrixNCM[x, #]&/@ y, b];


(* ::Text:: *)
(*Extracting commuting objects (no deltas on account of summation over heavy mass indices)*)


MatrixNCM[a___, (x:Except[_Delta|_Times|_Coupling]?MomCommuteQ), c___]:= x * MatrixNCM[a, c];


(* ::Subsubsection::Closed:: *)
(*Utility functions*)


MomAntiCommutator[x_, y_] := MomNCM[x, y]+ MomNCM[y, x];


(* ::Text:: *)
(*Drop the MomNCM after tracing*)


ReleaseMomNCM@ expr_:= expr/. MomNCM-> NonCommutativeMultiply;


(* ::Text:: *)
(*Expand out the MatrixNCM after X-term substitutions*)


ExpandMatrixNCM@ expr_:= expr//. {
		MatrixNCM[a___, (x:Except[_Delta|_Times]?MomCommuteQ) b_, c___]:> x * MatrixNCM[a, b, c],
		MatrixNCM[a___, (x:(InvProp[0]|InvProp[0]^_))* b_, c___]:> x * MatrixNCM[a, b, c],
		MatrixNCM[a___, b_Plus, c___]:> (MatrixNCM[a, #, c] & /@ b)
	}


(* ::Subsection:: *)
(*Momentum derivatives*)


(*To be deleted*)


(* ::Subsubsection::Closed:: *)
(*Auxiliary functions*)


LastPosition[x_, y_]:= Module[{temp},
	temp = FirstPosition[Reverse@ x, y, None, {1}];
	If[temp=== None, temp, Length@ x- First@ temp +1]
];


(* ::Text:: *)
(*Test if expression is independent of loop momentum (Subscript[LoopMom, \[Mu]]) *)


MomFreeQ[f_?MomFreeQ[x___]]:= And@@ MomFreeQ/@ {x};
MomFreeQ[f_[x___]]:= False;
MomFreeQ[_]:= True;
MomFreeQ[LoopMom| InvProp]:= False;


(* ::Subsubsection::Closed:: *)
(*Partial Momentum Derivative*)


(* ::Text:: *)
(*Rules for applying a partial momentum derivative to an individual term*)


MomDerivative[_?MomFreeQ, _]:= 0;
MomDerivative[x_Plus, \[Mu]_]:= MomDerivative[#, \[Mu]]&/@ x;
MomDerivative[Times[x_, y__], \[Mu]_]:= MomDerivative[x, \[Mu]] y + x MomDerivative[Times@ y, \[Mu]];
MomDerivative[Power[x_, n_], \[Mu]_]:= n Power[x, n-1] MomDerivative[x, \[Mu]];


MomDerivative[LoopMom@ Index[\[Mu]_,Lorentz], Index[\[Nu]_,Lorentz]]:= Metric[Index[\[Mu],Lorentz], Index[\[Nu],Lorentz]];
MomDerivative[_InvProp, Index[\[Mu]_,Lorentz]]:= 2 LoopMom@ Index[\[Mu],Lorentz];
MomDerivative[\[Gamma]@ LoopMom, Index[\[Mu]_,Lorentz]]:= \[Gamma]@ Index[\[Mu],Lorentz];


(* ::Text:: *)
(*Act with momentum derivatives to the right in an NCM product, and eventually terminate them.*)


TakeMomDerivatives@ expr_:= expr//. {
	MomNCM[a___, MomDev@ \[Mu]_, x:Except[_MomDev], b___]:>
		MomNCM[a, x, MomDev@ \[Mu], b] + MomNCM[a, MomDerivative[x, \[Mu]], b],
	MomNCM[___, _MomDev]:> 0
};


(* ::Text:: *)
(*Perform all momentum derivatives in a loop integral expression using IbP to reduce the number of terms.  *)


PerformMomDerivatives@ expr_:= expr/. x_MomNCM:> PerformMomDerivatives@ x;
PerformMomDerivatives@ expr_MomNCM:= Module[{pos},
	pos= LastPosition[expr, _MomDev];
	If[pos === None, Return@ expr];
	(*Decides whether to do IbP or not*)
	PerformMomDerivatives@ If[Count[expr[[;; pos- 1]], Except[_?MomFreeQ], {1}]- Count[expr[[pos+ 1;; ]], Except[_?MomFreeQ], {1}] > 0,
		MomNCM[expr[[;; pos- 1]], (TakeMomDerivatives@ expr[[pos;;]])]
	,
		pos= First@ FirstPosition[expr, _MomDev, None, {1}];
		-MomNCM[TakeMomDerivatives@ MomNCM[expr[[pos]], expr[[;;pos-1]]], expr[[pos+1;;]] ]
	]
];
