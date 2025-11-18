(* ::Package:: *)

(* ::Chapter:: *)
(*Simplification tests*)


(* ::Section:: *)
(*Setup*)


If[!MemberQ[$ContextPath, "Matchete`PackageScope`"], PrependTo[$ContextPath,"Matchete`PackageScope`"];]


LoadModel@ "SM";


DefineField[\[Psi], Fermion];
DefineField[\[Nu], Fermion, SelfConjugate-> True];
DefineField[\[Phi], Scalar, SelfConjugate-> True, Indices-> Flavor];
DefineCoupling[cA, SelfConjugate-> True, Indices-> {Flavor, Flavor}, Symmetries-> {AntisymmetricIndices[1, 2]}]
DefineCoupling[cS, SelfConjugate-> True, Indices-> {Flavor, Flavor}, Symmetries-> {SymmetricIndices[1, 2]}]
DefineCoupling[M, SelfConjugate-> True, Indices-> Flavor, DiagonalCoupling-> True]
DefineCoupling[m, SelfConjugate-> True]
DefineField[#, Scalar, Indices-> SU2L@ fund]&/@ {H1, H2, H3, H4};
DefineGlobalGroup[SU2g, SU@2];
DefineField[l2, Fermion, Indices-> SU2g@ fund, Chiral-> LeftHanded, Mass-> 0];
DefineField[e2, Fermion, Chiral-> RightHanded, Mass-> 0];


(* ::Section:: *)
(*Tests*)


(* ::Subsubsection::Closed:: *)
(*Reducible4dOpQ*)


VerificationTest[Module[{op, i, p, s, r, t, \[Mu], \[Nu], \[Rho], \[Sigma]},
	op= Bar@ l[i, p]** \[Gamma][\[Mu], \[Nu]]** e[r] Bar@ e[s]** \[Gamma][\[Rho], \[Sigma]]** l[i, t] LCTensor[\[Mu], \[Nu], \[Rho], \[Sigma]];
	Matchete`Simplifications`PackagePrivate`Reducible4dOpQ@ Operator@ op
],
	True
, TestID-> "4DReducibleOp1"]


VerificationTest[Module[{op, i, p, s, r, t, \[Mu], \[Nu], \[Rho]},
	op= Bar@ e[p]** \[Gamma][\[Mu], \[Nu]]** \[Gamma][\[Rho]]** e[r] Bar@ e[s]** \[Gamma][\[Mu], \[Nu], \[Rho]]** e[t];
	Matchete`Simplifications`PackagePrivate`Reducible4dOpQ@ Operator@ op
],
	True
, TestID-> "4DReducibleOp2"]


VerificationTest[Module[{op, i, p, s, r, t, \[Mu]},
	op= Bar@ e[p]** \[Gamma][\[Mu]]** e[r] Bar@ e[s]** \[Gamma][\[Mu]]** e[t];
	Matchete`Simplifications`PackagePrivate`Reducible4dOpQ@ Operator@ op
],
	False
, TestID-> "4DReducibleOp3"]


(* ::Subsection:: *)
(*Standard terms*)


(* ::Subsubsection::Closed:: *)
(*Canonical terms*)


VerificationTest[
	- Bar@ H@ i CD[{\[Mu],\[Mu]}, H@ i]// GreensSimplify
,
	Bar@ CD[d$$1, H@ d$$1] CD[d$$1, H@ d$$1]
, TestID-> "Higgs kinetic term"]


VerificationTest[
	- I CD[\[Mu], Bar@ l[i, p]] **\[Gamma]@\[Mu] ** l[i,p] //GreensSimplify //RelabelIndices
,
	I Bar@ l[d$$1, d$$1] ** \[Gamma]@ d$$1 ** CD[d$$1, l[d$$1, d$$1]]
, TestID-> "Lepton kinetic term"]


(* ::Subsubsection::Closed:: *)
(*Fermion mass terms*)


VerificationTest[
	Bar@ \[Psi][]** \[Psi][] + Bar@ \[Psi][]** PL** \[Psi][] + Bar@ \[Psi][]** PR** \[Psi][]// GreensSimplify
,
	2 Bar@ \[Psi][]** \[Psi][]
, TestID-> "VL fermion mass term"]


VerificationTest[
	Bar@ \[Nu][]** \[Nu][] + Bar@ \[Nu][]** PL** \[Nu][] + Bar@ \[Nu][]** PR** \[Nu][]// GreensSimplify
,
	2 Bar@ \[Nu][]** \[Nu][]
, TestID-> "Majorana fermion mass term"]


(* ::Subsection:: *)
(*Identities *)


(* ::Subsubsection::Closed:: *)
(*Group relations *)


(* ::Text:: *)
(*SU(2) identities with GroupFierz*)


VerificationTest[
	Block[{op1, op2, op3},
		op1=(CD[\[Nu],Bar@H@i]H@i Bar@H@j CD[\[Mu],H@k] - Bar@H@i CD[\[Nu],H@i] CD[\[Mu],Bar@H@j] H@k)CG[gen@SU2L@fund,{a,j,k}]FS[W,\[Mu],\[Nu],a];
		op2= CD[\[Mu],Bar@H@j] CD[\[Nu],H@k] Bar@H@i H@i CG[gen@SU2L@fund,{a,j,k}]FS[W,\[Mu],\[Nu],a];
		op3= CD[\[Mu],Bar@H@i] CD[\[Nu],H@i] Bar@H@j H@k CG[gen@SU2L@fund,{a,j,k}]FS[W,\[Mu],\[Nu],a];
		op1+op2+op3//GreensSimplify
	]
,
	0
, TestID-> "IdentitiesGroupFierz"]


(* ::Text:: *)
(*Schouten identities*)


VerificationTest[Module[{i, j, k, l},
	H1@ i H2@ j H3@ k H4@ l (
		Bar@ eps[SU2L][i, j] Bar@ eps[SU2L][k, l]+ 
		Bar@ eps[SU2L][j, k] Bar@ eps[SU2L][i, l]+ 
		Bar@ eps[SU2L][k, i] Bar@ eps[SU2L][j, l]
		)// GreensSimplify
],
	0
, TestID-> "IdentitiesGroupSchouten"]


(* ::Text:: *)
(*Symmetries*)


VerificationTest[
	Bar@ CG[eps@ SU2L, {i, j}] H@i H@ j// GreensSimplify
,
	0
, TestID-> "Higgs with eps[SU2L] 0"]


VerificationTest[
	Bar@ H@ i CD[\[Mu], H@ i]Bar@ H@ j CD[\[Nu], H@ j]LCTensor[\[Mu], \[Nu], \[Rho], \[Sigma]]FS[B, \[Rho], \[Sigma]]// GreensSimplify
,
	0
, TestID-> "Levi-Civita with derivatives 0"]


VerificationTest[
	Bar@CConj@\[Psi][]**\[Gamma][\[Mu],\[Nu]]**\[Psi][] Bar@\[Psi][]**\[Gamma][\[Mu],\[Nu]]**CConj@\[Psi][]// GreensSimplify
,
	0
, TestID-> "Spinor line transposition 0"]


VerificationTest[
	Bar@H[i]Bar@H[j]Bar@H[k](CD[\[Mu],Bar@e[p]]**\[Gamma][\[Mu]]**CConj[l[m,r]])CG[eps[SU2L],{i,k}]CG[eps[SU2L],{m,j}]// GreensSimplify
,
	0
, TestID-> "Levi-Civita with derivatives and derivatives/EoM in spin chain"]


VerificationTest[Module[{\[Mu], \[Nu], \[Rho], a},
	FS[B, \[Mu], \[Nu]]FS[B, \[Rho], \[Nu]]FS[B, \[Mu], \[Rho]]+ FS[B, \[Mu], \[Nu]]FS[W, \[Rho], \[Nu], a]FS[W, \[Mu], \[Rho], a]// GreensSimplify
],
	0
, TestID-> "Symmetry-vanishing triple gauge operators"]


(* ::Subsection:: *)
(*Evanescent reduction*)


(* ::Subsubsection::Closed:: *)
(*Global group*)


VerificationTest[Module[{lag, a, \[Mu]},
	lag= Bar@ l2[a]** e2[] Bar@ e2[]** l2[a] + 1/2 Bar@ l2@ a** \[Gamma]@\[Mu]** l2@ a Bar@ e2[]** \[Gamma][\[Mu]]** e2[];
	(*Simplifies to purely evanescent operators*)
	GreensSimplify[lag, ReductionIdentities-> FourDimensional]
],
	0
, TestID-> "Fierzing with global indices"]


(* ::Subsubsection::Closed:: *)
(*Gamma + LC reductions*)


VerificationTest[Module[{op},
	op= (Bar@ e[d$$1]** \[Gamma][d$$1, d$$3]** l[d$$1, d$$2]) (Bar@ l[d$$1, d$$3]** \[Gamma][d$$2, d$$4]** e[d$$4])
		LCTensor[d$$1, d$$2, d$$3, d$$4];
	(*Simplifies to purely evanescent operators*)
	op= GreensSimplify[op, ReductionIdentities-> Evanescent];
	{Count[op, _EvaOp, All], op/. _EvaOp-> 0}
],
	{1, 0}
, TestID-> "Gamma reduction + LCTensor"]


VerificationTest[Module[{op, i, p, \[Mu], \[Nu], \[Rho], \[Sigma]},
	op= Bar@ H@ i H@ i Bar@ e@ p** \[Gamma][\[Nu], \[Rho], \[Sigma]]** CD[\[Mu], e@ p] LCTensor[\[Mu], \[Nu], \[Rho], \[Sigma]];
	op= GreensSimplify[op, ReductionIdentities-> Evanescent];
	FreeQ[op, EvaOp[_, _, {Index[p_, Flavor], Index[p_, Flavor]}]]
],
	False
, TestID-> "Gamma reduction + LCTensor (flavor indices)"]


(* ::Subsubsection::Closed:: *)
(*Schouten identity*)


(* ::Text:: *)
(*Vanishing of BWW~ in SMEFT*)


VerificationTest[Module[{op, \[Mu], \[Nu], \[Rho], \[Sigma], \[Alpha], a},
	op= FS[B, \[Mu], \[Nu]]FS[W, \[Nu], \[Rho], a]FS[W, \[Sigma], \[Alpha], a]LCTensor[\[Mu], \[Rho], \[Sigma], \[Alpha]];
	op= GreensSimplify[op, ReductionIdentities-> Evanescent];
	MatchQ[op, _EvaOp|Times[_EvaOp, _?NumericQ]]
],
	True
, TestID-> "LorentzSchouten simplifcication"]


(* ::Subsection:: *)
(*Compound operators*)


(* ::Subsubsection::Closed:: *)
(*Spin chain transposition*)


(* ::Text:: *)
(*Identifies symmetry in spin-chain transposition and identifies the flavor symmetry of Weinberg + dim-7 operators*)


VerificationTest[Module[{op5, op7, lag},
		DefineCoupling[Cll, Indices-> {Flavor, Flavor}, Symmetries-> {SymmetricIndices[1, 2]}];
		DefineCoupling[CllB, Indices-> {Flavor, Flavor}, Symmetries-> {AntisymmetricIndices[1, 2]}];
		op5= Bar@ CConj@ l[i, p]**l[j, r]H@ k H@ l Bar@CG[eps@ SU2L, {i, k}]Bar@ CG[eps@ SU2L, {j, l}];
		op7= FS[B, \[Mu], \[Nu]]Bar@ CConj@ l[i, p]** \[Sigma][\[Mu], \[Nu]]** l[j,r]H@ k H@ l Bar@ CG[eps@ SU2L, {i, k}]Bar@ CG[eps@ SU2L, {j, l}];
		lag= -Cll[p, r] op5- CllB[p, r] op7// PlusHc;
		lag= InternalSimplify[lag, InternalOpRepresentation-> True];
		{
			(*No composite operators*)
			FreeQ[lag, CompOp]
		, 
			(*The two starting operators + conjugates*)
			Sort@ Cases[lag, AtomicOp[type_, id_, _]-> {type, id}, Infinity]
		,
			(*op5 is symmetric. op7 is antisymmetric*)
			{$operators[{{H, H, l, l}, 0}, 1]@ Symmetries, $operators[{{H, H, l, l}, 2}, 1]@ Symmetries}
		}
	]
,
	{
		True
	,
		{
			{{{H, H, l, l}, 0}, 1}, 
			{{{H, H, l, l}, 2}, 1},
			{{{Conj[H], Conj[H], Conj[l], Conj[l]}, 0}, 1},
			{{{Conj[H], Conj[H], Conj[l], Conj[l]}, 2}, 1}
		}
	,
		{{{1, 2}-> 1, {2, 1}-> 1}, {{1, 2}-> 1, {2, 1}-> -1}}
	}
, TestID-> "Flavor symmetries of \[CapitalDelta]L=2 operators"]


(* ::Subsubsection::Closed:: *)
(*Symmetrization of four-fermion operators*)


VerificationTest[Module[{expr, \[Mu], p, r},
		expr= Bar@ e[p]** \[Gamma][\[Mu]]** e[r] Bar@ e[r]** \[Gamma][\[Mu]]** e[p] - 
			Bar@ e[p]** \[Gamma][\[Mu]]** e[p] Bar@ e[r]** \[Gamma][\[Mu]]** e[r];
		expr= GreensSimplify[expr, ReductionIdentities-> Evanescent];
		{Count[expr, _EvaOp, All]> 0, expr/. _EvaOp-> 0}
	]
,
	{True, 0}
, TestID-> "Four-e operator symmetries"]


(* ::Subsection:: *)
(*Coefficient simplification*)


(* ::Subsubsection::Closed:: *)
(*CoefficientPatterns*)


(* ::Text:: *)
(*Symmetrization of coefficients*)


VerificationTest[Module[{coef, i, j, k, l, openInds},
	coef= 2cS[i, k]cS[j, k]+ 2 cA[i, k]cA[j, k]+ cS[i, l]cS[l, j]+ 
		cA[i, l]cA[l, j] +cA[i, j] +cS[i, j]+ 5 cS[i, k]cA[k, l] cS[l, j];
	openInds= {Index[i, Flavor], Index[j, Flavor]};
	{ (*Symmetric coefficient*)
		Matchete`Simplifications`PackagePrivate`CanonizeCouplingContractions[coef, openInds, {{1, 2}-> 1, {2, 1}-> 1}]- 
		(3 cS[i, k]cS[j, k] + cA[i, k]cA[j, k]+ cS[i, j])
	, (*Antisymmetric coefficient*)
		Matchete`Simplifications`PackagePrivate`CanonizeCouplingContractions[coef, openInds, {{1, 2}-> 1, {2, 1}-> -1}]-
		(5 cS[i, k]cA[k, l] cS[l, j]+ cA[i, j])
	}
],
	{0, 0}
, TestID-> "CanonizeCouplingContractions: couplings"]


(* ::Text:: *)
(*Tests with constants and Deltas:*)
(*Nb. if the coefficient does not meaningfully contract into the operator (a mismatch of open indices) it should not be removed by anti-symmetrization.*)


VerificationTest[Module[{coef, i, j, k, l, openInds},
	openInds= {Index[i, Flavor], Index[j, Flavor]};
	coef= 1+ 2 Delta@@ openInds;
	{ (*Symmetric coefficient*)
		Matchete`Simplifications`PackagePrivate`CanonizeCouplingContractions[coef, openInds, {{1, 2}-> 1, {2, 1}-> 1}]- 
		(1+ 2 Delta@@ openInds)		
	, (*Antisymmetric coefficient*)
		Matchete`Simplifications`PackagePrivate`CanonizeCouplingContractions[coef, openInds, {{1, 2}-> 1, {2, 1}-> -1}]-
		1
	}
],
	{0, 0}
, TestID-> "CanonizeCouplingContractions: deltas and constants"]


(* ::Text:: *)
(*Tests Diagonal indices *)


VerificationTest[Module[{coef, i, j, openInds},
	openInds= {Index[i, Flavor], Index[j, Flavor]};
	coef= M@ i (1 + M@ j) Delta@@ openInds;
	{ (*Symmetric coefficient*)
		Matchete`Simplifications`PackagePrivate`CanonizeCouplingContractions[coef, openInds, {{1, 2}-> 1, {2, 1}-> 1}]- 
		Expand[M@ i (1 + M@ j) Delta@@ openInds]
	, (*Antisymmetric coefficient*)
		Matchete`Simplifications`PackagePrivate`CanonizeCouplingContractions[coef, openInds, {{1, 2}-> 1, {2, 1}-> -1}]
	}
],
	{0, 0}
, TestID-> "CanonizeCouplingContractions: Diagonal indices"]


(* ::Subsubsection::Closed:: *)
(*Coupling symmetries*)


VerificationTest[
	cA[i, j] \[Phi][i] \[Phi][j]//GreensSimplify
,
	0
, TestID-> "Coefficient vs operator symmetry"]


VerificationTest[Module[{s,t,a,b,op},
	op= LF[{M[s], m[]}, {2, 1, 0}] (Bar@ u[a, t]** \[Gamma][\[Mu]]** u[a, s])(Bar@ u[b, s]** \[Gamma][\[Mu]]** u[b, t]);
	op- Bar@ op// GreensSimplify
],
	0
, TestID-> "Term Hermitianity with flavor deltas"]


(* ::Section:: *)
(*Cleanup*)


ResetAll[]
