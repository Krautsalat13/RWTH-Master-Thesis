(* ::Package:: *)

(* ::Title:: *)
(*Definitions test unit*)


(* ::Section:: *)
(*Tests*)


ResetAll[]


(* ::Subsection:: *)
(*Input definitions*)


(* ::Subsubsection::Closed:: *)
(*Field definition*)


VerificationTest[
	DefineGaugeGroup[SU2,SU@2,g2,W];
	DefineField[\[Psi],Fermion,Indices->{SU2[fund]}];\[Psi][a],
	Field[\[Psi],Fermion,{Index[a,SU2[fund]]},{}],
	TestID->"Fields definition"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Field deletion*)


VerificationTest[
	DefineField[\[Psi],Fermion];RemoveField[\[Psi]]; Head[\[Psi][]],
	\[Psi],
	TestID->"Field deletion"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Masses and indices (heavy fields)*)


VerificationTest[
	DefineFlavorIndex[Flavor,3];
	DefineField[\[CapitalPhi],Scalar,Mass->{Heavy,M\[CapitalPsi],{Flavor}},Indices->{Flavor}];
	M\[CapitalPsi][f]
	,
	Coupling[M\[CapitalPsi],{Index[f,Flavor]},0],
	TestID->"Masses and indices (heavy fields)"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Masses and indices (light fields)*)


VerificationTest[
	DefineFlavorIndex[Flavor,3];
	DefineField[\[Psi],Fermion,Mass->{Light,m\[Psi],{Flavor}},Indices->{Flavor}];
	m\[Psi][f1,f2]
	,
	Coupling[m\[Psi],{Index[f1,Flavor],Index[f2,Flavor]},1],
	TestID->"Masses and indices (light fields)"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Projectors onto chiral fields*)


VerificationTest[
	DefineField[\[Psi],Fermion,Chiral->LeftHanded];
	PR**\[Psi][]
	,
	0,
	TestID->"Projectors onto chiral fields"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Fields strength tensors are created *)


VerificationTest[
	DefineGaugeGroup[SU2,SU@2,g2,W];
	DefineGaugeGroup[SU3,SU@3,g3,G];
	DefineField[V,Vector,Indices->{SU3[adj],SU2[adj]},SelfConjugate->True];
	FS[V,\[Mu],\[Nu],a,i]
	,
	FieldStrength[V,{Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]},{Index[a,SU3[adj]],Index[i,SU2[adj]]},{}],
	TestID->"Fields strength tensors are created"
]
ResetAll[];


(* ::Subsection:: *)
(*Bar operation*)


(* ::Subsubsection::Closed:: *)
(*Onto fields*)


VerificationTest[
	DefineGaugeGroup[SU2,SU@2,g2,W];
	DefineGaugeGroup[SU3,SU@3,g3,G];
	DefineField[V,Vector,Indices->{SU3[adj],SU2[adj]},SelfConjugate->False];
	Bar@V[\[Mu],a,i]
	,
	Bar@Field[V,Vector[Index[\[Mu],Lorentz]],{Index[a,SU3[adj]],Index[i,SU2[adj]]},{}],
	TestID->"Bar onto fields"
]


(* ::Subsubsection::Closed:: *)
(*Onto fields strength tensors*)


VerificationTest[
	Bar@FS[V,\[Mu],\[Nu],a,i]
	,
	Bar@FieldStrength[V,{Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]},{Index[a,SU3[adj]],Index[i,SU2[adj]]},{}],
	TestID->"Bar onto fields strength tensors"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Onto selfconjugate fields*)


VerificationTest[
	DefineGaugeGroup[SU2,SU@2,g2,W];
	DefineGaugeGroup[SU3,SU@3,g3,G];
	DefineField[V,Vector,Indices->{SU3[adj],SU2[adj]},SelfConjugate->True];
	Bar@V[\[Mu],a,i]
	,
	Field[V,Vector[Index[\[Mu],Lorentz]],{Index[a,SU3[adj]],Index[i,SU2[adj]]},{}],
	TestID->"Bar onto selfconjugate fields"
]


(* ::Subsubsection::Closed:: *)
(*Onto fields strength tensors of a selfconjugate field*)


VerificationTest[
	Bar@FS[V,\[Mu],\[Nu],a,i]
	,
	FieldStrength[V,{Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]},{Index[a,SU3[adj]],Index[i,SU2[adj]]},{}],
	TestID->"Bar onto fields strength tensors of a selfconjugate field"
]


(* ::Subsubsection::Closed:: *)
(*Bar acting on gauge field implicitly contracted with generator *)


VerificationTest[
	Bar@ FieldStrength[W,{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, 
		{Index[i, SU3@ fund], Bar@ Index[j, SU3@ fund]}, {}]
,
	FieldStrength[W,{Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, 
		{Bar@ Index[i, SU3@ fund], Index[j, SU3@ fund]}, {}]
,
	TestID-> "Bar acting on gauge field implicitly contracted with generator"
]
ResetAll[];


(* ::Subsection:: *)
(*Free Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*Scalar*)


VerificationTest[
	DefineGaugeGroup[SU2,SU@2,g2,W];
	DefineField[\[Phi],Scalar,Mass->{Light,m},Indices->{SU2[fund]}];
	FreeLag[\[Phi]]
	,
	CD[d$$1,Bar[\[Phi][d$$1]]]CD[d$$1,\[Phi][d$$1]] - m[]^2 Bar[\[Phi][d$$1]]\[Phi][d$$1],
	TestID->"Free Lagrangian: complex scalar"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Fermions*)


VerificationTest[
	DefineFlavorIndex[Flavor,3];
	DefineGaugeGroup[SU2,SU@2,g2,W];
	DefineGaugeGroup[SU3,SU@3,g3,G];
	DefineField[\[CapitalPsi],Fermion,Mass->{Heavy,M,{Flavor}},Indices->{SU2[fund],SU3[fund],Flavor}];
	FreeLag[\[CapitalPsi]]
	,
	I Bar[\[CapitalPsi][d$$1,d$$1,d$$1]]**\[Gamma][d$$1]**CD[d$$1,\[CapitalPsi][d$$1,d$$1,d$$1]]-M[d$$1] Bar[\[CapitalPsi][d$$1,d$$1,d$$1]]**\[CapitalPsi][d$$1,d$$1,d$$1],
	TestID->"Free Lagrangian: fermions"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Chiral fermions*)


VerificationTest[
	DefineFlavorIndex[Flavor,3];
	DefineField[\[Psi],Fermion,Mass->{Light,m,{Flavor}},Indices->{Flavor},Chiral->LeftHanded];
	FreeLag[\[Psi]]//Expand
	,
	I Bar[\[Psi][d$$1]]**\[Gamma][d$$1]**CD[d$$1,\[Psi][d$$1]]-1/2 m[d$$1,d$$2] Bar[CConj[\[Psi][d$$1]]]**\[Psi][d$$2]- 1/2 Bar@m[d$$1,d$$2] Bar[\[Psi][d$$1]]** CConj@\[Psi]@d$$2,
	TestID->"Free Lagrangian: chiral fermions"
]
ResetAll[];


VerificationTest[
	DefineFlavorIndex[Flavor, 3];
	DefineGaugeGroup[SU2, SU@ 2, g, W];
	DefineField[\[Psi], Fermion, Mass-> {Light, m, {Flavor}}, Indices-> {SU2@ fund, Flavor}, Chiral-> LeftHanded];
	
	I Bar@ \[Psi][a, i]** \[Gamma]@\[Mu]**CD[\[Mu],\[Psi][a, i]] - PlusHc[1/2 m[i, j] Bar@ CG[eps@ SU2, {a, b}] Bar@ CConj@ \[Psi][a, i]** \[Psi][b, j]]- 
		FreeLag@ \[Psi]// GreensSimplify
,
	0
,
	TestID->"Free Lagrangian: chiral fermions in pseudoreal rep"
]
ResetAll[];


(* ::Subsubsection::Closed:: *)
(*Vector*)


VerificationTest[
	DefineGaugeGroup[SU2,SU@2,g2,W];
	DefineGaugeGroup[SU3,SU@3,g3,G];
	DefineField[V,Vector,Indices->{SU3[adj],SU2[adj]},SelfConjugate->True];
	FreeLag[V]//Simplify
	,
	-1/4 FS[V,d$$1,d$$2,d$$1,d$$1]^2 + 1/2 MV[]^2 V[d$$1,d$$1,d$$1]^2 //Simplify,
	TestID->"Free Lagrangian: real vector"
]
ResetAll[];


(* ::Section:: *)
(*Loading the package*)


(* ::Text:: *)
(*Un-comment and run first to see the tests directly in this file. *)
(*Comment when finished and add file name to the list testfiles in  the script "RunTests.wls" .*)


(*	PrependTo[$Path, ParentDirectory[NotebookDirectory[]]];
	<< Matchete`*)
