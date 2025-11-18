(* ::Package:: *)

(* ::Title:: *)
(*NCM test unit*)


(* ::Section:: *)
(*Tests*)


If[!MemberQ[$ContextPath, "Matchete`PackageScope`"], PrependTo[$ContextPath,"Matchete`PackageScope`"];]


ResetAll[];
DefineField[\[Psi]M,Fermion,SelfConjugate->True];
DefineField[\[Psi]M2,Fermion,SelfConjugate->True];
DefineField[\[Psi]D,Fermion,SelfConjugate->False];
DefineField[\[Psi]L,Fermion,  Chiral->LeftHanded];
DefineField[\[Psi]R,Fermion,  Chiral->RightHanded];
DefineField[\[Phi],Scalar];


(* ::Subsection::Closed:: *)
(*Contractions*)


VerificationTest[
	Bar@\[Psi]L[]**\[Phi][]**\[Psi]R[],
	\[Phi][]Bar@\[Psi]L[]**\[Psi]R[],
	TestID->"Extract scalar from spin chain"
]


VerificationTest[
	Bar@\[Psi]L[]**\[Gamma][\[Mu]]** (I CD[\[Mu],\[Phi][]])**\[Psi]L[],
	I CD[\[Mu],\[Phi][]]Bar@\[Psi]L[]**\[Gamma][\[Mu]]**\[Psi]L[],
	TestID->"Extract derivative scalar from spin chain"
]


VerificationTest[
	Bar@\[Psi]D[]**(\[Gamma][\[Mu]] Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[])**\[Psi]D[],
	Bar@\[Psi]D[]**\[Gamma][\[Mu]]**\[Psi]D[] Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[],
	TestID->"Contraction of spinors 1"
]


VerificationTest[
	Bar@\[Psi]D[]**(\[Gamma][\[Mu]]**\[Gamma][5] Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[])**\[Psi]D[],
	Bar@\[Psi]D[]**\[Gamma][\[Mu]]**\[Gamma][5]**\[Psi]D[] Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[],
	TestID->"Contraction of spinors 2"
]


VerificationTest[
	Bar@\[Psi]L[]**(Coupling[m,0]Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[])**\[Psi]R[],
	Coupling[m,0]Bar@\[Psi]L[]**\[Psi]R[] Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[],
	TestID->"Contraction of spinors 3"
]


(*VerificationTest[
	\[Psi]R[]**(Bar@\[Psi]M[]**\[Gamma][5])**(Bar@\[Psi]D[]**\[Gamma][\[Mu]]) **\[Psi]D[]**\[Psi]M[]**Bar@\[Psi]L[],
	-Bar@\[Psi]L[]**\[Psi]R[]Bar@\[Psi]M[]**\[Gamma][5]**\[Psi]M[] Bar@\[Psi]D[]**\[Gamma][\[Mu]]**\[Psi]D[],
	TestID->"Contraction of spinors 4"
]*)
(*AET: I do not think this is the behavior we want*)


(*VerificationTest[
	(\[Gamma][\[Mu]]**\[Psi]D[]) **NonCommutativeMultiply@Bar@\[Psi]D[] //Simplify,
	-Bar@\[Psi]D[]**\[Gamma][\[Mu]]**\[Psi]D[] //Simplify,
	TestID->"Trace of spin chain"
]*)
(*I have no idea why this test does not want to run in the script but works here...*)


(*Bar@\[Psi]L[]**Bar@\[Psi]D[] **\[Psi]D[]**\[Psi]R[]  //TraditionalForm*) (*Even single fields need NCM warpping*)


(* ::Subsection::Closed:: *)
(*Bar & Transpose*)


VerificationTest[
	Bar@ (Bar@CConj@\[Psi]M[]**\[Psi]D[]),
	Bar@\[Psi]D[]**\[Psi]M[],
	TestID->"Bar of Majorana with Dirac fermions 1"
]


VerificationTest[
	Bar@ (Bar@\[Psi]D[]**\[Psi]M[]),
	Bar@CConj@\[Psi]M[]**\[Psi]D[],
	TestID->"Bar of Majorana with Dirac fermions 2"
]


VerificationTest[
	Bar@ (Bar@\[Psi]M[]**\[Psi]M[]),
	Bar@\[Psi]M[]**\[Psi]M[] ,
	TestID->"Bar of Majorana mass"
]


(* ::Subsection:: *)
(*Spin Chain*)


(* ::Subsubsection:: *)
(*LOpen*)


VerificationTest[
	LOpenSpinChainQ[Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[]],
	False ,
	TestID->"LOpenSpinChainQ 1"
]


VerificationTest[
	LOpenSpinChainQ[\[Gamma][\[Mu]]**\[Psi]M[]],
	True ,
	TestID->"LOpenSpinChainQ 2"
]


VerificationTest[
	LOpenSpinChainQ[Bar@\[Psi]M[]**\[Gamma][\[Mu]]],
	False ,
	TestID->"LOpenSpinChainQ 3"
]


VerificationTest[
	LOpenSpinChainQ[\[Gamma][\[Mu]]],
	True ,
	TestID->"LOpenSpinChainQ 4"
]


VerificationTest[
	LOpenSpinChainQ[\[Psi]M[]],
	True ,
	TestID->"LOpenSpinChainQ 5"
]


VerificationTest[
	LOpenSpinChainQ[Transp@\[Psi]M[]],
	False ,
	TestID->"LOpenSpinChainQ 6"
]


VerificationTest[
	LOpenSpinChainQ[Bar@\[Psi]M[]],
	False ,
	TestID->"LOpenSpinChainQ 7"
]


VerificationTest[
	LOpenSpinChainQ[Transp@Bar@\[Psi]M[]],
	True ,
	TestID->"LOpenSpinChainQ 8"
]


VerificationTest[
	LOpenSpinChainQ[Bar@CD[\[Mu],\[Psi]M[]]],
	False ,
	TestID->"LOpenSpinChainQ 9"
]


(* ::Subsubsection::Closed:: *)
(*ROpen*)


VerificationTest[
	ROpenSpinChainQ[Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[]],
	False ,
	TestID->"ROpenSpinChainQ 1"
]


VerificationTest[
	ROpenSpinChainQ[\[Gamma][\[Mu]]**\[Psi]M[]],
	False ,
	TestID->"ROpenSpinChainQ 2"
]


VerificationTest[
	ROpenSpinChainQ[Bar@\[Psi]M[]**\[Gamma][\[Mu]]],
	True ,
	TestID->"ROpenSpinChainQ 3"
]


VerificationTest[
	ROpenSpinChainQ[\[Gamma][\[Mu]]**\[Gamma][5]],
	True ,
	TestID->"ROpenSpinChainQ 4"
]


VerificationTest[
	ROpenSpinChainQ[\[Psi]M[]],
	False ,
	TestID->"ROpenSpinChainQ 5"
]


VerificationTest[
	ROpenSpinChainQ[Transp@\[Psi]M[]],
	True ,
	TestID->"ROpenSpinChainQ 6"
]


VerificationTest[
	ROpenSpinChainQ[Bar@\[Psi]M[]],
	True ,
	TestID->"ROpenSpinChainQ 7"
]


VerificationTest[
	ROpenSpinChainQ[-Transp@Bar@\[Psi]M[]],
	False ,
	TestID->"ROpenSpinChainQ 8"
]


VerificationTest[
	ROpenSpinChainQ[Bar@CD[\[Mu],\[Psi]M[]]],
	True ,
	TestID->"ROpenSpinChainQ 9"
]


(* ::Subsubsection::Closed:: *)
(*Closed*)


VerificationTest[
	Bar@\[Psi]M[]**\[Gamma][\[Mu]]**\[Psi]M[]// ClosedSpinChainQ,
	True, 
	TestID->"ClosedSpinChainQ 1"
]


VerificationTest[
	Bar@\[Psi]M[]**\[Gamma][\[Mu]]// ClosedSpinChainQ,
	False, 
	TestID->"ClosedSpinChainQ 2"
]


VerificationTest[
	\[Psi]D[]// ClosedSpinChainQ,
	False, 
	TestID->"ClosedSpinChainQ 3"
]


VerificationTest[
	\[Phi][]// ClosedSpinChainQ,
	False, 
	TestID->"ClosedSpinChainQ 4"
]


(* ::Subsection::Closed:: *)
(*Transpose canonical ordering*)


VerificationTest[
	Bar@\[Psi]L[]**\[Psi]R[],
	Transp[Bar@\[Psi]L[]**\[Psi]R[]],
	TestID->"Canonical ordering 1"
]


VerificationTest[
	Bar@CConj@\[Psi]M[]**\[Psi]D[],
	Transp[Bar@CConj@\[Psi]M[]**\[Psi]D[]],
	TestID->"Canonical ordering 2"
]


VerificationTest[
	Bar@CConj@\[Psi]M[]**\[Psi]M2[],
	Transp[Bar@CConj@\[Psi]M[]**\[Psi]M2[]],
	TestID->"Canonical ordering 3"
]


VerificationTest[
	Bar@CConj@\[Psi]D[]**\[Psi]R[],
	Transp[Bar@CConj@\[Psi]D[]**\[Psi]R[]],
	TestID->"Canonical ordering 4"
]


VerificationTest[
	Bar@\[Psi]D[]**Transp@Bar@\[Psi]R[],
	Transp[Bar@\[Psi]D[]**Transp@Bar@\[Psi]R[]],
	TestID->"Canonical ordering 5"
]


(* ::Subsection:: *)
(*Actively acting on spin chains*)


(* ::Subsubsection:: *)
(*CanonizeSpinorLines*)


VerificationTest[
	Bar@\[Psi]L[]**\[Gamma][\[Mu],\[Nu]]**\[Psi]R[]// CanonizeSpinorLines,
	Bar@\[Psi]L[]**\[Gamma][\[Mu],\[Nu]]**\[Psi]R[],
	TestID->"CanonizeSpinorLines: nothing 1"
]


VerificationTest[
	\[Psi]M[]**Bar@\[Psi]R[]// CanonizeSpinorLines,
	\[Psi]M[]**Bar@\[Psi]R[],
	TestID->"CanonizeSpinorLines: nothing 2"
]


VerificationTest[
	CConj@\[Psi]D[]**Bar@CConj@ \[Psi]R[]// CanonizeSpinorLines,
	CConj@\[Psi]D[]**Bar@CConj@ \[Psi]R[],
	TestID->"CanonizeSpinorLines: nothing 3"
]


VerificationTest[
	\[Gamma][\[Mu], \[Nu]]** PL** \[Psi]M[]// CanonizeSpinorLines,
	\[Gamma][\[Mu], \[Nu]]** PL** \[Psi]M[],
	TestID->"CanonizeSpinorLines: nothing 4"
]


VerificationTest[
	CheckAbort[\[Psi]L[]**\[Gamma][\[Mu],\[Nu]]**\[Psi]R[]// CanonizeSpinorLines,True],
	True,
	Matchete`PackageScope`CanonizeSpinorLines::lopens,
	TestID->"CanonizeSpinorLines: abort 1"
]


VerificationTest[
	CheckAbort[Bar@\[Psi]L[]**\[Gamma][\[Mu],\[Nu]]**Bar@\[Psi]R[]// CanonizeSpinorLines, True],
	True,
	Matchete`PackageScope`CanonizeSpinorLines::ropens,
	TestID->"CanonizeSpinorLines: abort 2"
]


VerificationTest[
	Bar@ CConj@\[Psi]M[] **Bar@\[Psi]L[]**\[Gamma][\[Mu],\[Nu]]**\[Psi]R[]// CanonizeSpinorLines,
	Bar@\[Psi]L[]**\[Gamma][\[Mu],\[Nu]]**\[Psi]R[] Bar@ CConj@\[Psi]M[],
	TestID->"CanonizeSpinorLines: split chain 1"
]


VerificationTest[
	Bar@CConj@\[Psi]L[] **Bar@\[Psi]D[]**\[Gamma][\[Mu]]**PL**\[Psi]D[]** \[Gamma]@\[Mu]**CConj@\[Psi]L[]// CanonizeSpinorLines,
	-Bar@\[Psi]D[]**\[Gamma][\[Mu]]**PL**\[Psi]D[] Bar@\[Psi]L[]**\[Gamma]@\[Mu]**\[Psi]L[],
	TestID->"CanonizeSpinorLines: split chain 2"
]


VerificationTest[
	Bar@CConj@\[Psi]L[] ** \[Gamma]@\[Mu]**CConj@\[Psi]L[] **Bar@\[Psi]D[] ** \[Psi]D[]// CanonizeSpinorLines,
	-Bar@\[Psi]D[] ** \[Psi]D[] Bar@\[Psi]L[]**\[Gamma]@\[Mu]**\[Psi]L[],
	TestID->"CanonizeSpinorLines: split chain 3"
]


(* ::Subsubsection:: *)
(*FermionTrace*)


(*VerificationTest[
	CheckAbort[\[Psi]D[]**Bar@ \[Psi]D[] NonCommutativeMultiply@ \[Gamma]@\[Mu]// FermionTrace, True],
	True,
	Matchete`PackageScope`FermionTrace::multopens,
	TestID->"FermionTrace: abort 1"
]*)


(*VerificationTest[
	CheckAbort[\[Psi]D[]**Bar@\[Psi]D[]**\[Gamma]@\[Mu]**\[Psi]L[]// FermionTrace, True],
	True,
	Matchete`PackageScope`FermionTrace::cantclose,
	TestID->"FermionTrace: abort 2"
]*)


(*VerificationTest[
	CConj@ \[Psi]D[]**Bar@ CConj@ \[Psi]D[] ** PR// FermionTrace,
	-Bar@ \[Psi]D[]** PR** \[Psi]D[],
	TestID->"FermionTrace: single line 1"
]*)


(*VerificationTest[
	\[Gamma]@\[Mu]**\[Gamma]@\[Nu]// FermionTrace,
	4 Metric[Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]],
	TestID->"FermionTrace: single line 2"
]*)


(*VerificationTest[
	Bar@ CConj@ \[Psi]D[] ** PR** CConj@ \[Psi]D[] \[Psi]M[]**Bar@\[Psi]M[]// FermionTrace,
	-Bar@ \[Psi]D[]** PR** \[Psi]D[] Bar@ \[Psi]M[]** \[Psi]M[],
	TestID->"FermionTrace: multiple lines 1"
]*)


(*VerificationTest[
	Bar@ \[Psi]D[] ** PR** \[Psi]D[] + \[Psi]M[]**Bar@\[Psi]M[]// FermionTrace,
	4 Bar@ \[Psi]D[]** PR** \[Psi]D[] - Bar@ \[Psi]M[]** \[Psi]M[],
	TestID->"FermionTrace: multiple lines 2"
]*)


(* ::Subsection:: *)
(*Examples from the CDE*)


(*VerificationTest[
	Matchete`SuperTrace`PackagePrivate`GopFerm[f,{i,j},2]//RelabelIndices
	,
	1/2 I DiracProduct[GammaM[Index[d$$1,Lorentz]]]**
		Matchete`PackageScope`Gop[f,{i,j},{Index[d$$1,Lorentz],Index[d$$2,Lorentz]},{}]**
		Matchete`PackageScope`MomDev[Index[d$$2,Lorentz]]
	,
	TestID->"CDE_ferm_1"
]*)


(*VerificationTest[
	Matchete`SuperTrace`PackagePrivate`GopFerm[f,{i,j},3]//RelabelIndices
	,
	1/3 DiracProduct[GammaM[Index[d$$1,Lorentz]]]**
		Matchete`PackageScope`Gop[f,{i,j},{Index[d$$1,Lorentz],Index[d$$2,Lorentz]},{Index[d$$3,Lorentz]}]**
		Matchete`PackageScope`MomDev[Index[d$$3,Lorentz]]**
		Matchete`PackageScope`MomDev[Index[d$$2,Lorentz]]
	,
	TestID->"CDE_ferm_2"
]*)


(* ::Section:: *)
(*Cleanup after tests *)


ResetAll[]


(* ::Section:: *)
(*Loading the package*)


(* ::Text:: *)
(*Un-comment and run first to see the tests directly in this file. *)
(*Comment when finished and add file name to the list testfiles in  the script "RunTests.wls" .*)


(*	PrependTo[$Path, ParentDirectory[NotebookDirectory[]]];
	<< Matchete`*)
