(* ::Package:: *)

(* ::Title:: *)
(*Flavor index test unit*)


(* ::Section:: *)
(*Tests*)


If[!MemberQ[$ContextPath, "Matchete`PackageScope`"], PrependTo[$ContextPath,"Matchete`PackageScope`"];]


DefineFlavorIndex[Flavor,nf];


DefineCoupling[M,SelfConjugate->True,Indices->{Flavor},DiagonalCoupling->{True}]


DefineCoupling[Y,Indices->{Flavor,Flavor}]


(* ::Subsection:: *)
(*Delta contractions for flavor indices*)


VerificationTest[
	ContractDelta[
		Delta[Index[p,Flavor],Index[r,Flavor]]Delta[Index[r,Flavor],Index[s,Flavor]]
	],
	Delta[Index[p,Flavor],Index[s,Flavor]],
	TestID->"FlavorDelta1"
]


VerificationTest[
	ContractDelta[
		Delta[Index[p,Flavor],Index[r,Flavor]]Delta[Index[r,Flavor],Index[p,Flavor]]
	],
	nf,
	TestID->"FlavorDelta2"
]


VerificationTest[
	ContractDelta[
		1/M[p]^2 Delta[Index[p,Flavor],Index[r,Flavor]]Delta[Index[r,Flavor],Index[s,Flavor]]
	],
	Delta[Index[p,Flavor],Index[s,Flavor]]/M[p]^2,
	TestID->"FlavorDelta3"
]


VerificationTest[
	ContractDelta[
		1/M[p]^2 Delta[Index[p,Flavor],Index[r,Flavor]]Delta[Index[r,Flavor],Index[p,Flavor]]
	],
	FlavorSum[Index[p,Flavor]]/M[p]^2,
	TestID->"FlavorDelta4"
]


(* ::Subsection:: *)
(*RelabelIndices*)


(* ::Text:: *)
(*Un-comment to see MessageFailure. (commented for commit...)*)


VerificationTest[
	RelabelIndices[
		FlavorSum[Index[p,Flavor]]1/M[p]
	],
	FlavorSum[Index[d$$1,Flavor]] * 1/M[d$$1],
	TestID->"RelabelFlavorIndices1"
]


VerificationTest[
	RelabelIndices[
		1/M[p] Y[p,r]Bar@Y[s,r]
	],
	1/M[p] Y[p,d$$1]Bar@Y[s,d$$1],
	TestID->"RelabelFlavorIndices2"
]


(* ::Subsection:: *)
(*Evaluate Loop Functions*)


VerificationTest[
	EvaluateLoopFunctions@ContractDelta[
		Delta[Index[p,Flavor],Index[r,Flavor]]Delta[Index[r,Flavor],Index[p,Flavor]]LF[{M[p]},{1,1}]
	],
	nf+FlavorSum[Index[p,Flavor]] Log[\[Mu]bar2/M[p]^2],
	TestID->"FlavorSum"
]


(* ::Subsection:: *)
(*Contractions for FlavorSum*)


VerificationTest[
	FlavorSum[i]/.$DropDiagonalCouplings[[1]],
	1,
	TestID->"DropDiagonalCouplings"
](* check that "_FlavorSum->1" is first argument of DropDiagonalCouplings *)


VerificationTest[
	Quiet@ContractDelta[2(1+FlavorSum[Index[p,Flavor]]Delta[Index[p,Flavor],Index[r,Flavor]]),Expand->False],
	2 (1+FlavorSum[Index[r,Flavor]]),
	TestID->"Contract FlavorSum"
]


(* ::Section:: *)
(*Info for evaluating test:*)


(* ::Text:: *)
(*1. Each test can be evaluated individually directly in this file by first loading the package locally.*)
(**)
(*2. All tests contained in "TestUnit_Template.wl" can be evaluated at the same time from any other notebook. (See an example in "EvaluateTemplateTests.nb")*)
(**)
(*3. All tests in "TestUnit_Template.wl" and other .wl unit tests specified in the list testfiles in "RunTests.wls" script will be run before any commit to the local master branch and push to the online master branch. *)
(*If a test unit contains at least one failed test, the commit/push will be interrupted and cancelled (the testing procedure will not proceed to the other test units) . *)
(*A minimal error message with the failed tests ID will be displayed in the terminal, more extensive information can be found from a notebook using point 2, or directly in the unit test using point 1. *)
