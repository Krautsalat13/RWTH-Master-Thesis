(* ::Package:: *)

(* ::Title:: *)
(*CovariantDerivative test unit*)


If[!MemberQ[$ContextPath, "Matchete`PackageScope`"], PrependTo[$ContextPath,"Matchete`PackageScope`"];]


(* ::Section:: *)
(*Covariant derivatives acting on fields *)


DefineGaugeGroup[SU3c, SU@ 3, g3, G];
DefineGaugeGroup[U1Y, U1, g1, B];
DefineFlavorIndex[Flavor, 3];


DefineField[S, Scalar, Indices-> {SU3c@ adj}, SelfConjugate-> True, Mass-> Heavy];
DefineField[H, Scalar, Charges-> {U1Y[1/2]}];
DefineField[T, Scalar, Indices-> {Bar@ SU3c@ fund}];
DefineField[F, Fermion, Indices-> {SU3c@ fund, Flavor}, Charges-> {U1Y[-2]}];


DefineCoupling[y]


(* ::Subsection:: *)
(*Taking derivatives of expressions*)


VerificationTest[
	TakeDev[{}, 3 y[] S[A] CG[gen@ SU3c@ fund, {A, i ,j}] Bar@ F[i, a]** PL** F[j, a]]
,
	3 y[] CG[gen@ SU3c@ fund, {A, i ,j}] S[A] Bar@ F[i, a]** PL** F[j, a]
,
	TestID-> "Take 0 derivatives of a compund expression"
]


VerificationTest[
	Expand@ TakeDev[{Index[\[Mu], Lorentz]}, 3 y[] S[A] CG[gen@ SU3c@ fund, {A, i ,j}] Bar@ F[i, a]** PL** F[j, a]]
,
	3 y[] CG[gen@ SU3c@ fund, {A, i ,j}] (CD[\[Mu], S[A]] Bar@ F[i, a]** PL** F[j, a]+ 
		S[A] CD[\[Mu], Bar@ F[i, a]]** PL** F[j, a] + S[A] Bar@ F[i, a]** PL** CD[\[Mu], F[j, a]] )// Expand
,
	TestID-> "Take 1 derivative of a compund expression"
]


VerificationTest[
	Expand@ TakeDev[{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, S[A] F[i, a]]
,
	CD[{\[Mu], \[Nu]}, S@ A] F[i, a] + CD[\[Mu], S[A]] CD[\[Nu], F[i, a]] + 
		CD[\[Nu], S[A]] CD[\[Mu], F[i, a]] +S[A] CD[{\[Mu], \[Nu]}, F[i, a]]
,
	TestID-> "Take 2 derivatives of a compund expression"
]


(* ::Subsection:: *)
(*Commute derivatives*)


(* ::Subsubsection::Closed:: *)
(*GAction*)


VerificationTest[
	GAction[{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, H[]]
,
	1/2 FS[B, \[Mu], \[Nu]] H[]
,
	TestID-> "GAction on charged scalar"
]


VerificationTest[
	GAction[{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, S@ A]// RelabelIndices
,
	I CG[fStruct@ SU3c, {A, B, C}] FS[G, \[Mu],  \[Nu], B] S@ C// RelabelIndices
,
	TestID-> "GAction on colored scalar"
]


VerificationTest[
	GAction[{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, T@ a]// RelabelIndices
,
	-FS[G, \[Mu], \[Nu], A] Bar@ CG[gen@ SU3c@ fund, {A, a, b}] T@ b// RelabelIndices
,
	TestID-> "GAction on anti-fundamental scalar"
]


VerificationTest[
	GAction[{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, F[a, i]]// RelabelIndices
,
	FS[G, \[Mu], \[Nu], A] CG[gen@ SU3c@ fund, {A, a, b}]F[b, i] +
	(-2) FS[B, \[Mu], \[Nu]] F[a ,i] //RelabelIndices
,
	TestID-> "GAction on colored and charged fermion"
]


(* ::Subsubsection::Closed:: *)
(*CommuteCDs*)


VerificationTest[
	CommuteCDs[CD[{\[Mu], \[Nu]}, H[]], 1]
,
	CD[{\[Nu], \[Mu]}, H[]] -1/2 I FS[B, \[Mu], \[Nu]] H[]
,
	TestID-> "Commute derivatives on charged scalar"
]


VerificationTest[
	CommuteCDs[CD[{\[Mu], \[Nu], \[Rho]}, H[]], 2] // Expand
,
	CD[{\[Mu], \[Rho], \[Nu]}, H[]] -1/2 I CD[\[Mu], FS[B, \[Nu], \[Rho]] H[]]// Expand
,
	TestID-> "Commute inner derivatives on charged scalar"
]


VerificationTest[
	CommuteCDs[CD[{\[Mu], \[Nu]}, F[a, i]], 1]// RelabelIndices
,
	CD[{\[Nu], \[Mu]}, F[a, i]] - I GAction[{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, F[a, i]]// RelabelIndices
,
	TestID-> "Commute derivatives on charged and colored Fermion"
]


(* ::Subsection:: *)
(*Remove field definitions*)


ResetAll[];


(* ::Section:: *)
(*Loading the package*)


(* ::Text:: *)
(*Un-comment and run first to see the tests directly in this file. *)
(*Comment when finished and add file name to the list testfiles in  the script "RunTests.wls" .*)


(*	PrependTo[$Path, ParentDirectory[NotebookDirectory[]]];
	<< Matchete`*)
