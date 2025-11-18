(* ::Package:: *)

(* ::Title:: *)
(*Dirac Algebra test unit*)


(* ::Section:: *)
(*Tests*)


If[!MemberQ[$ContextPath, "Matchete`PackageScope`"], PrependTo[$ContextPath,"Matchete`PackageScope`"];]


(* ::Subsubsection::Closed:: *)
(*Input shortcuts*)


VerificationTest[
	\[Gamma][\[Mu]],
	DiracProduct[GammaM[Index[\[Mu],Lorentz]]],
	TestID->"\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)"
]


VerificationTest[
	\[Gamma][\[Mu],\[Nu],\[Rho]],
	DiracProduct[GammaM[Index[\[Mu],Lorentz],Index[\[Nu],Lorentz],Index[\[Rho],Lorentz]]],
	TestID->"\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\[Nu]\[Rho]\)]\)"
]


VerificationTest[
	\[Gamma][Index[\[Mu],Lorentz]],
	DiracProduct[GammaM[Index[\[Mu],Lorentz]]],
	TestID->"\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) with Index[_,Lorentz]"
]


VerificationTest[
	\[Gamma][Index[\[Mu],Lorentz],Index[\[Nu],Lorentz],Index[\[Rho],Lorentz]],
	DiracProduct[GammaM[Index[\[Mu],Lorentz],Index[\[Nu],Lorentz],Index[\[Rho],Lorentz]]],
	TestID->"\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\[Nu]\[Rho]\)]\) with Index[_,Lorentz]"
]


VerificationTest[
	1**\[Sigma][\[Mu],\[Nu]],
	I/2 (\[Gamma][\[Mu]]**\[Gamma][\[Nu]]-\[Gamma][\[Nu]]**\[Gamma][\[Mu]])  //RefineDiracProducts,
	TestID->"\[Sigma][\[Mu],\[Nu]]"
]


VerificationTest[
	1**\[Sigma][Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]],
	I/2 (\[Gamma][\[Mu]]**\[Gamma][\[Nu]]-\[Gamma][\[Nu]]**\[Gamma][\[Mu]])  //RefineDiracProducts,
	TestID->"\[Sigma][\[Mu],\[Nu]] with Index[_,Lorentz]"
]


(* ::Subsubsection::Closed:: *)
(*Dirac product*)


VerificationTest[
	DiracProduct[Gamma5,Gamma5],
	1,
	TestID->"\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\) contraction"
]


VerificationTest[
	DiracProduct[Gamma5, Gamma[Index[\[Mu],Lorentz]]],
	-DiracProduct[Gamma[Index[\[Mu],Lorentz]],Gamma5],
	TestID->"\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\) \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) anticommutation"
]


VerificationTest[
	DiracProduct[Proj[+1], Proj[-1]],
	0,
	TestID->"Orthogonal projectors"
]


VerificationTest[
	DiracProduct[Proj[-1], Gamma5,GammaM[Index[\[Mu],Lorentz]], Proj[+1]],
	-DiracProduct[GammaM[Index[\[Mu],Lorentz]],Proj[1]],
	TestID->"Projectors, \!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\) and \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)"
]


(* ::Subsubsection::Closed:: *)
(*Transp*)


VerificationTest[
	DiracProduct[Transp@Gamma5,Transp@Gamma5],
	1,
	TestID->"Subscript[\[Gamma], 5]^T contraction"
]


VerificationTest[
	DiracProduct[Transp@Proj[+1], Transp@Proj[-1]],
	0,
	TestID->"Orthogonal transposed projectors"
]


VerificationTest[
	DiracProduct[Transp@Gamma5, Transp@Gamma[Index[\[Mu],Lorentz]]],
	-DiracProduct[Transp@Gamma[Index[\[Mu],Lorentz]],Transp@Gamma5],
	TestID->"Subscript[\[Gamma], 5]^T Subscript[\[Gamma], \[Mu]]^T anticommutation"
]


VerificationTest[
	DiracProduct[Transp@Proj[-1], Gamma5,GammaM[Index[\[Mu],Lorentz]], Transp@Proj[+1]],
	-DiracProduct[GammaM[Index[\[Mu],Lorentz]],Transp@Proj[1]],
	TestID->"Tranposed projectors, \!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\) and \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)"
]


(* ::Subsection:: *)
(*Dirac manipulations*)


(* ::Subsubsection::Closed:: *)
(*Dirac Trace*)


VerificationTest[
	DiracTrace[c*DiracProduct[GammaM[Index[\[Mu],Lorentz]], GammaM[Index[\[Nu],Lorentz]],GammaM[Index[\[Rho],Lorentz]],GammaM[Index[\[Sigma],Lorentz]],Proj[+1]]] //Simplify,
	c(2 Metric[Index[\[Mu],Lorentz],Index[\[Sigma],Lorentz]] Metric[Index[\[Nu],Lorentz],Index[\[Rho],Lorentz]]-2 Metric[Index[\[Mu],Lorentz],Index[\[Rho],Lorentz]] Metric[Index[\[Nu],Lorentz],Index[\[Sigma],Lorentz]]+2 Metric[Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]] Metric[Index[\[Rho],Lorentz],Index[\[Sigma],Lorentz]]-2 I LCTensor[Index[\[Mu],Lorentz],Index[\[Nu],Lorentz],Index[\[Rho],Lorentz],Index[\[Sigma],Lorentz]]) //Simplify,
	TestID->"Dirac Trace of 4 \[Gamma] and a projector"
]


VerificationTest[
	DiracTrace[DiracProduct[GammaM[Index[\[Alpha],Lorentz]],GammaM[Index[\[Beta],Lorentz]],GammaM[Index[\[Mu],Lorentz],Index[\[Nu],Lorentz]]]] //Simplify,
	4 Metric[Index[\[Alpha],Lorentz],Index[\[Nu],Lorentz]] Metric[Index[\[Beta],Lorentz],Index[\[Mu],Lorentz]]-4 Metric[Index[\[Alpha],Lorentz],Index[\[Mu],Lorentz]] Metric[Index[\[Beta],Lorentz],Index[\[Nu],Lorentz]] //Simplify, 
	TestID->"Dirac Trace with \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \(\(\*SubscriptBox[\(\[Mu]\), \(1\)] ... \) \*SubscriptBox[\(\[Mu]\), \(n\)]\)]\)"
]


(* ::Subsubsection::Closed:: *)
(*Refine Dirac Product*)


VerificationTest[
	RefineDiracProducts[DiracProduct[GammaM[Index[\[Mu],Lorentz]], GammaM[Index[\[Alpha],Lorentz]],GammaM[Index[\[Mu],Lorentz]]] ] //ContractMetric //Simplify,
	(2-\[ScriptD]) DiracProduct[GammaM[Index[\[Alpha],Lorentz]]]//Simplify,
	TestID->"Contract two \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)"
]


VerificationTest[
	RefineDiracProducts[DiracProduct[Transp@GammaM[Index[\[Mu],Lorentz]], Transp@GammaM[Index[\[Alpha],Lorentz]],Transp@GammaM[Index[\[Mu],Lorentz]]] ] //ContractMetric //Simplify,
	(2-\[ScriptD]) DiracProduct[Transp@GammaM[Index[\[Alpha],Lorentz]]]//Simplify,
	TestID->"Contract two transposed \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)"
]


VerificationTest[
	RefineDiracProducts[\[Gamma]@\[Mu] ** \[Gamma]@\[Nu]] - (Metric[Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]] + NonCommutativeMultiply@ \[Gamma][\[Mu], \[Nu]]),
	0,
	TestID->"RefineDiracProducts of \[Gamma][\[Mu]]**\[Gamma][\[Nu]]"
]


VerificationTest[
	RefineDiracProducts@ Transp@DiracProduct[GammaM@Index[\[Sigma],Lorentz],GammaM[Index[\[Nu],Lorentz],Index[\[Rho],Lorentz]]]-
	(DiracProduct[Transp[GammaM[Index[\[Nu], Lorentz], Index[\[Rho], Lorentz], Index[\[Sigma], Lorentz]]]] + 
 DiracProduct[Transp[GammaM[Index[\[Rho], Lorentz]]]]*Metric[Index[\[Nu], Lorentz], Index[\[Sigma], Lorentz]] - 
 DiracProduct[Transp[GammaM[Index[\[Nu], Lorentz]]]]*Metric[Index[\[Rho], Lorentz], Index[\[Sigma], Lorentz]]),	
	0,
	TestID->"RefineDiracProducts of \!\(\*TemplateBox[{TagBox[SuperscriptBox[SubscriptBox[\"\[CapitalGamma]\", StyleBox[TemplateBox[{\"\[Nu]\", \"\[Rho]\"}, \"RowDefault\"], StripOnInput -> False, FontSize -> 12]], \"T\"], DisplayForm], TagBox[SuperscriptBox[SubscriptBox[\"\[Gamma]\", StyleBox[\"\[Sigma]\", StripOnInput -> False, FontSize -> 12]], \"T\"], DisplayForm]},\n\"RowDefault\"]\)"
]


VerificationTest[
	RefineDiracProducts[\[Gamma][\[Sigma],\[Nu],\[Mu]]],
	-\[Gamma][\[Mu],\[Nu],\[Sigma]],
	TestID->"RefineDiracProducts of \[Gamma][\[Sigma],\[Nu],\[Mu]]"
]


VerificationTest[
	\[Gamma][\[Mu],\[Nu],\[Rho]] ** \[Gamma][\[Mu],\[Nu],\[Rho]]// RefineDiracProducts
,
	-2 \[ScriptD]+ 3 \[ScriptD]^2- \[ScriptD]^3
,
	TestID-> "RefineDiracProduct \[Gamma][\[Mu],\[Nu],\[Rho]] ** \[Gamma][\[Mu],\[Nu],\[Rho]]"
]


(* ::Subsubsection::Closed:: *)
(*CollectGammaMatrices (internal function)*)


VerificationTest[
	CollectGammaMatrices[CC**\[Gamma]@\[Mu]**PL ** \[Gamma]@\[Nu]**\[Gamma]@5],
	CC ** \[Gamma]@\[Mu] ** \[Gamma]@\[Nu] ** PR,
	TestID->"CollectGammaMatrices on the left"
]


(* ::Subsection:: *)
(*Operations*)


(* ::Subsubsection::Closed:: *)
(*Fermion tests*)


VerificationTest[
	LOpenSpinChainQ[\[Gamma][\[Mu]]**PL**Field[cL,Fermion,{},{}]], 
	True,
	TestID->"LOpenSpinChainQ True"
]


VerificationTest[
	LOpenSpinChainQ[Bar@Field[cL,Fermion,{},{}]**\[Gamma][\[Mu]]**PL], 
	False,
	TestID->"LOpenSpinChainQ False"
]


VerificationTest[
	ROpenSpinChainQ[\[Gamma][\[Mu]]**PL**Field[cL,Fermion,{},{}]],
	False,
	TestID->"ROpenSpinChainQ False"
]


VerificationTest[
	ROpenSpinChainQ[Bar@Field[cL,Fermion,{},{}]**\[Gamma][\[Mu]]**PL],
	True,
	TestID->"ROpenSpinChainQ True"
]


VerificationTest[
	ClosedSpinChainQ[Bar@Field[cL,Fermion,{},{}]**\[Gamma][\[Mu]]**PL**Field[cL,Fermion,{},{}]],
	True,
	TestID->"ClosedSpinChainQ True"
]


VerificationTest[
	ClosedSpinChainQ[Field[c1,Fermion,{},{}]**(Bar@Field[c2,Fermion,{},{}])**\[Gamma][\[Mu]]**PL**Field[c3,Fermion,{},{}]],
	False,
	TestID->"ClosedSpinChainQ True"
]


VerificationTest[
	MajoranaQ[Field[c1,Fermion,{},{}]**(Bar@Field[c2,Fermion,{},{}])**\[Gamma][\[Mu]]**PL**Field[c3,Fermion,{},{}]],
	False,
	TestID->"MajoranaQ False"
]


DefineField[c4,Fermion, SelfConjugate->True];


VerificationTest[
	MajoranaQ[c4[]**(Bar@Field[c2,Fermion,{},{}])**\[Gamma][\[Mu]]**PL**Field[c3,Fermion,{},{}]],
	True,
	TestID->"MajoranaQ True"
]


(* ::Subsubsection::Closed:: *)
(*CC*)


VerificationTest[
	CC**CC,
	-1,
	TestID->"Product of two C matrices"
]


VerificationTest[
	\[Gamma][5]**Transp@CC,
	-CC**Transp@\[Gamma][5],
	TestID->"Product of \!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\) and \!\(\*SuperscriptBox[\(C\), \(T\)]\)"
]


VerificationTest[
	Transp@PR**CC,
	CC**PR,
	TestID->"Product of \!\(\*SuperscriptBox[SubscriptBox[\(P\), \(R\)], \(T\)]\) and C"
]


VerificationTest[
	\[Gamma][\[Mu]]**\[Gamma][\[Nu]]**CC,
	CC**Transp@\[Gamma][\[Mu]]**Transp@\[Gamma][\[Nu]],
	TestID->"Product of \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\) and C"
]


VerificationTest[
	\[Gamma][\[Mu],\[Nu]]**CC,
	-CC**Transp@\[Gamma][\[Mu],\[Nu]],
	TestID->"Product of \!\(\*SubscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\) and C"
]


(* ::Subsubsection::Closed:: *)
(*Transp*)


VerificationTest[
	Transp@Transp@\[Gamma][\[Mu]],
	\[Gamma][\[Mu]],
	TestID->"Twice transposed"
]


VerificationTest[
	Transp[],
	1,
	TestID->"Empty transpose"
]


VerificationTest[
	Transp[-1/2],
	-1/2,
	TestID->"Tranpose of a number"
]


VerificationTest[
	Transp[Field[\[Phi],Vector[\[Mu]],{},{}]],
	Field[\[Phi],Vector[\[Mu]],{},{}],
	TestID->"Tranpose of a vector"
]


VerificationTest[
	Transp[Bar@Field[\[Phi],Scalar,{},{}]],
	Bar@Field[\[Phi],Scalar,{},{}],
	TestID->"Tranpose of a scalar bar"
]


VerificationTest[
	Head[Transp[Field[\[Psi],Fermion,{},{}]]]===Transp,
	True,
	TestID->"Tranpose of a fermion "
]


VerificationTest[
	Head[Transp[Bar@Field[\[Psi],Fermion,{},{}]]]===Transp,
	True,
	TestID->"Tranpose of a fermion bar"
]


VerificationTest[
	Transp@(Bar@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][\[Mu]]** Field[\[Psi]2,Fermion,{},{}]**Transp@Field[\[Psi]3,Fermion,{},{}]),
	-Field[\[Psi]3,Fermion,{},{}]**Transp@Field[\[Psi]2,Fermion,{},{}]**Transp@\[Gamma][\[Mu]]**Bar@Transp@Field[\[Psi]1,Fermion,{},{}],
	TestID->"Tranpose on multiple fermions in NCM"
]


VerificationTest[
	CanonizeSpinorLines@Transp@(Bar@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][\[Mu]]** Field[\[Psi]2,Fermion,{},{}]**Transp@Field[\[Psi]3,Fermion,{},{}]),
	Transp[(Bar@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][\[Mu]]** Field[\[Psi]2,Fermion,{},{}])**Transp@Field[\[Psi]3,Fermion,{},{}]],
	TestID->"Tranpose on multiple fermions in NCM"
]


(* ::Subsubsection::Closed:: *)
(*Bar On Fermion lines *)


VerificationTest[
	Bar@Bar@Field[\[Psi]3,Fermion,{},{}],
	Field[\[Psi]3,Fermion,{},{}],
	TestID->"Bar of bar"
]


VerificationTest[
	Bar[Bar@Field[\[Psi]1,Fermion,{},{}]**Field[\[Psi]2,Fermion,{},{}]],
	Bar@Field[\[Psi]2,Fermion,{},{}]**Field[\[Psi]1,Fermion,{},{}],
	TestID->"Bar of scalar current"
]


VerificationTest[
	Bar[Bar@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][5]**Field[\[Psi]2,Fermion,{},{}]],
	-Bar@Field[\[Psi]2,Fermion,{},{}]**\[Gamma][5]**Field[\[Psi]1,Fermion,{},{}],
	TestID->"Bar of pseudo-scalar current"
]


VerificationTest[
	Bar[Bar@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][\[Mu]]**Field[\[Psi]2,Fermion,{},{}]],
	Bar@Field[\[Psi]2,Fermion,{},{}]**\[Gamma][\[Mu]]**Field[\[Psi]1,Fermion,{},{}],
	TestID->"Bar of vector current"
]


VerificationTest[
	Bar[Bar@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][\[Mu]]**\[Gamma][5]**Field[\[Psi]2,Fermion,{},{}]],
	Bar@Field[\[Psi]2,Fermion,{},{}]**\[Gamma][\[Mu]]**\[Gamma][5]**Field[\[Psi]1,Fermion,{},{}],
	TestID->"Bar of axial vector current"
]


VerificationTest[
	Bar[Bar@Field[\[Psi]1,Fermion,{},{}]**\[Sigma][\[Mu],\[Nu]]**Field[\[Psi]2,Fermion,{},{}]]//RefineDiracProducts,
	Bar@Field[\[Psi]2,Fermion,{},{}]**\[Sigma][\[Mu],\[Nu]]**Field[\[Psi]1,Fermion,{},{}],
	TestID->"Bar of tensor current"
]


(* ::Subsubsection::Closed:: *)
(*CC in closed fermion line (CConj)*)


VerificationTest[
	Transp[Bar@CConj@Field[\[Psi]1,Fermion,{},{}]**Field[\[Psi]2,Fermion,{},{}]],
	Bar@CConj@Field[\[Psi]2,Fermion,{},{}]**Field[\[Psi]1,Fermion,{},{}],
	TestID->"Transpose scalar current with CC"
]


VerificationTest[
	Transp[Bar@CConj@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][5]**Field[\[Psi]2,Fermion,{},{}] ],
	Bar@CConj@Field[\[Psi]2,Fermion,{},{}]**\[Gamma][5]**Field[\[Psi]1,Fermion,{},{}] ,
	TestID->"Transpose pseudo-scalar current with CC"
]


VerificationTest[
	Transp[Bar@CConj@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][\[Mu]]**Field[\[Psi]2,Fermion,{},{}]],
	-Bar@CConj@Field[\[Psi]2,Fermion,{},{}]**\[Gamma][\[Mu]]**Field[\[Psi]1,Fermion,{},{}] ,
	TestID->"Transpose vector current with CC"
]


VerificationTest[
	Transp[Bar@CConj@Field[\[Psi]1,Fermion,{},{}]**\[Gamma][\[Mu]]**\[Gamma][5]**Field[\[Psi]2,Fermion,{},{}]],
	Bar@CConj@Field[\[Psi]2,Fermion,{},{}]**\[Gamma][\[Mu]]**\[Gamma][5]**Field[\[Psi]1,Fermion,{},{}] ,
	TestID->"Transpose axial-vector current with CC"
]


VerificationTest[
	Transp[Bar@CConj@Field[\[Psi]1,Fermion,{},{}]**\[Sigma][\[Mu],\[Nu]]**Field[\[Psi]2,Fermion,{},{}]],
	-Bar@CConj@Field[\[Psi]2,Fermion,{},{}]**\[Sigma][\[Mu],\[Nu]]**Field[\[Psi]1,Fermion,{},{}] ,
	TestID->"Transpose tensor current with CC"
]


(* ::Subsection:: *)
(*Four-dimensional identities*)


(* ::Subsubsection::Closed:: *)
(*Fierz in 4D*)


LoadModel["SM"];


VerificationTest[
	Fierz[(Bar@ l[i,p]**e[t])( Bar@q[a,i,r]**d[a,s]),Evanescent->False]//Expand, 
	-(1/2)(Bar@ l[i,p]** d[a,s])(Bar@ q[a,i,r] ** e[t])-1/8 (Bar@ l[i,p]**\[Sigma][\[Mu],\[Nu]]** d[a,s])(Bar@ q[a,i,r] **\[Sigma][\[Mu],\[Nu]]** e[t]) //RelabelIndices,
	TestID->"Fierz: PR x PR"
]


VerificationTest[
	Fierz[(Bar@ l[i,p]** e[r] )( Bar@ e[s]** l[i,t]),Evanescent->False], 
	-(1/2)(Bar@ l[i,p] ** \[Gamma][\[Mu]] ** l[i,t])(Bar@ e[s] ** \[Gamma][\[Mu]] ** e[r] ) //RelabelIndices,
	TestID->"Fierz: PR x PL"
]


VerificationTest[
	Fierz[(Bar@ l[i,p]**  \[Gamma][\[Mu]]**l[j,r] )(Bar@ l[k,s]**  \[Gamma][\[Mu]]**l[m,s] ),Evanescent->False], 
	(Bar@ l[i,p]**\[Gamma][\[Mu]]**l[m,s] )(Bar@ l[k,s]**  \[Gamma][\[Mu]]**l[j,r] ) //RelabelIndices,
	TestID->"Fierz: \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)PL x \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)PL"
]


VerificationTest[
	Fierz[(Bar@ l[i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r])( Bar@d[\[Alpha],s]**\[Gamma][\[Nu]]**e[t]),Evanescent->False], 
	- 2(Bar@ l[i,p]** e[t])(Bar@ d[\[Alpha],s] **q[\[Alpha],i,r]) //RelabelIndices,
	TestID->"Fierz: \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)PL x \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)PR"
]


VerificationTest[
	Fierz[-(Bar@ l[i,p]**\[Sigma][\[Mu],\[Nu]]**e[t])( Bar@q[\[Alpha],i,r]**\[Sigma][\[Mu],\[Nu]]**d[\[Alpha],s]),Evanescent->False], 
	 6(Bar@ l[i,p]** d[\[Alpha],s])(Bar@ q[\[Alpha],i,r] **e[t]) - 1/2 (Bar@ l[i,p]**\[Sigma][\[Mu],\[Nu]]** d[\[Alpha],s])(Bar@ q[\[Alpha],i,r]**\[Sigma][\[Mu],\[Nu]] **e[t]) //RelabelIndices,
	TestID->"Fierz: \!\(\*SubscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)PR x \!\(\*SubscriptBox[\(\[Sigma]\), \(\[Mu]\[VeryThinSpace]\[Nu]\)]\)PR"
]


VerificationTest[
	Fierz[(Bar@CConj@ e[p]**  e[r] )(Bar@ e[s]**  CConj@e[t] ),Evanescent->False], 
	1/2 (Bar@ e[s]**\[Gamma][\[Mu]]** e[r] )(Bar@ e[t]**\[Gamma][\[Mu]]** e[p] )//RelabelIndices,
	TestID-> "Fierz: C PL x C PR"
]


VerificationTest[
	Fierz[(Bar@ l[i,p]**CConj@q[\[Alpha],i,r])( Bar@CConj@d[\[Alpha],s]**e[t]),Evanescent->False]//Expand, 
	-(1/2)(Bar@ l[i,p]** e[t])(Bar@ q[a,i,r] ** d[a,s])+1/8 (Bar@ l[i,p]**\[Sigma][\[Mu],\[Nu]]** e[t])(Bar@ q[a,i,r] **\[Sigma][\[Mu],\[Nu]]** d[a,s]) //RelabelIndices,
	TestID-> "Fierz: C PR x C PR"
]


DefineField[\[Psi],Fermion];


VerificationTest[
	-Fierz[(Bar@ \[Psi][]**\[Psi][])( Bar@\[Psi][]**\[Psi][]),Evanescent->False]//ProjExpand, 
	1/4 (Bar@\[Psi][]**\[Psi][])^2+1/4 (Bar@\[Psi][]**\[Gamma][5]**\[Psi][])^2+1/4 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Psi][])^2+1/8 (Bar@\[Psi][]**\[Sigma][\[Mu],\[Nu]]**\[Psi][])^2-1/4 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Gamma][5]**\[Psi][])^2//RelabelIndices,
	TestID->"Fierz: S x S"
]


VerificationTest[
	-Fierz[(Bar@ \[Psi][]**\[Gamma][\[Mu]]**\[Psi][])( Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Psi][]),Evanescent->False]//ProjExpand, 
	(Bar@\[Psi][]**\[Psi][])^2- (Bar@\[Psi][]**\[Gamma][5]**\[Psi][])^2-1/2 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Psi][])^2-1/2 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Gamma][5]**\[Psi][])^2//RelabelIndices,
	TestID->"Fierz: V x V"
]


VerificationTest[
	Fierz[-(Bar@ \[Psi][]**\[Sigma][\[Mu],\[Nu]]**\[Psi][])( Bar@\[Psi][]**\[Sigma][\[Mu],\[Nu]]**\[Psi][]),Evanescent->False]//ProjExpand, 
	3 (Bar@\[Psi][]**\[Psi][])^2+3 (Bar@\[Psi][]**\[Gamma][5]**\[Psi][])^2-1/2 (Bar@\[Psi][]**\[Sigma][\[Mu],\[Nu]]**\[Psi][])^2//RelabelIndices,
	TestID->"Fierz: T x T"
]


VerificationTest[
	-Fierz[(Bar@ \[Psi][]**\[Gamma][\[Mu]]**\[Gamma][5]**\[Psi][])( Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Gamma][5]**\[Psi][]),Evanescent->False]//ProjExpand, 
	-(Bar@\[Psi][]**\[Psi][])^2+ (Bar@\[Psi][]**\[Gamma][5]**\[Psi][])^2-1/2 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Psi][])^2-1/2 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Gamma][5]**\[Psi][])^2//RelabelIndices,
	TestID->"Fierz: A x A"
]


VerificationTest[
	-Fierz[(Bar@ \[Psi][]**\[Gamma][5]**\[Psi][])( Bar@\[Psi][]**\[Gamma][5]**\[Psi][]),Evanescent->False]//ProjExpand, 
	1/4 (Bar@\[Psi][]**\[Psi][])^2+1/4 (Bar@\[Psi][]**\[Gamma][5]**\[Psi][])^2-1/4 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Psi][])^2+1/8 (Bar@\[Psi][]**\[Sigma][\[Mu],\[Nu]]**\[Psi][])^2+1/4 (Bar@\[Psi][]**\[Gamma][\[Mu]]**\[Gamma][5]**\[Psi][])^2//RelabelIndices,
	TestID->"Fierz: P x P"
]


(* ::Subsubsection::Closed:: *)
(*Gamma Reduction*)


VerificationTest[
	EpsExpand[GammaReduction[(Bar@ e[p]**\[Gamma][\[Mu]]**\[Gamma][\[Nu]]**l[i,r])( Bar@ e[s]**\[Gamma][\[Nu]]**\[Gamma][\[Mu]]**l[i,t]),Evanescent->False],Order->1], 
	(4-2\[Epsilon]) (Bar@ e[p]**l[i,r])( Bar@ e[s]**l[i,t])  + (Bar@ e[p]**\[Sigma][\[Mu],\[Nu]]**l[i,r])( Bar@ e[s]**\[Sigma][\[Mu],\[Nu]]**l[i,t])   //RelabelIndices,
	TestID->"GammaReduction: \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\) PL x \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\) PL"
]


VerificationTest[
	EpsExpand[GammaReduction[(Bar@ e[p]**\[Gamma][\[Mu]]**\[Gamma][\[Nu]]**l[i,r])(Bar@q[\[Alpha],i,s]**\[Gamma][\[Nu]]**\[Gamma][\[Mu]]**d[\[Alpha],t]) ,Evanescent->False],Order->1], 
	4(1-2\[Epsilon]) (Bar@ e[p]**l[i,r])(Bar@q[\[Alpha],i,s]**d[\[Alpha],t])   //RelabelIndices,
	TestID->"GammaReduction: \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\) PL x \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\) \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\) PR"
]


VerificationTest[
	EpsExpand[GammaReduction[(Bar@ l[i,p]**\[Gamma][\[Mu]]**\[Gamma][\[Nu]]**\[Gamma][\[Lambda]]**l[i,r])(Bar@q[\[Alpha],j,s]**\[Gamma][\[Lambda]]**\[Gamma][\[Nu]]**\[Gamma][\[Mu]]**q[\[Alpha],j,t]) ,Evanescent->False],Order->1], 
	4(1-2\[Epsilon]) (Bar@ l[i,p]**\[Gamma][\[Mu]]**l[i,r])(Bar@q[\[Alpha],j,s]**\[Gamma][\[Mu]]**q[\[Alpha],j,t])  //RelabelIndices,
	TestID->"GammaReduction: \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Lambda]\)]\)PL x \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Lambda]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)PL"
]


VerificationTest[
	EpsExpand[GammaReduction[(Bar@ l[i,p]**\[Gamma][\[Mu]]**\[Gamma][\[Nu]]**\[Gamma][\[Lambda]]**l[i,r])( Bar@ e[s]**\[Gamma][\[Lambda]]**\[Gamma][\[Nu]]**\[Gamma][\[Mu]]** e[t]),Evanescent->False],Order->1]//Expand, 
	16(1-\[Epsilon]) (Bar@ l[i,p]**\[Gamma][\[Mu]]**l[i,r])( Bar@ e[s]**\[Gamma][\[Mu]]** e[t]) //RelabelIndices,
	TestID->"GammaReduction: \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Lambda]\)]\)PL x \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Lambda]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)PR"
]


VerificationTest[
	Expand@EpsExpand[GammaReduction[(Bar@ e[p]**\[Gamma][\[Mu]]**\[Gamma][\[Nu]]**\[Sigma][\[Lambda],\[Rho]]**l[i,r])(Bar@d[\[Alpha],s]**\[Sigma][\[Lambda],\[Rho]]**\[Gamma][\[Nu]]**\[Gamma][\[Mu]]**q[\[Alpha],i,t]),Evanescent->False ],Order->1], 
	16(3-5\[Epsilon]) (Bar@ e[p]**l[i,r])(Bar@d[\[Alpha],s]**q[\[Alpha],i,t]) + 2(6-7\[Epsilon])  (Bar@ e[p]**\[Sigma][\[Mu],\[Nu]]**l[i,r])(Bar@d[\[Alpha],s]**\[Sigma][\[Mu],\[Nu]]**q[\[Alpha],i,t]) //RelabelIndices,
	TestID->"GammaReduction: \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\[Sigma]\), \(\[Lambda]\[Rho]\)]\)PL x \!\(\*SubscriptBox[\(\[Sigma]\), \(\[Lambda]\[Rho]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)PL"
]


ResetAll[];


(* ::Subsubsection::Closed:: *)
(*LC2Gamma5*)


VerificationTest[
	LCTensor[\[Mu], \[Nu], \[Rho], \[Sigma]] \[Gamma][\[Sigma]]// LC2Gamma5
,
	-I \[Gamma][\[Mu], \[Nu], \[Rho]]** \[Gamma][5]
,
	TestID-> "LCTensor contracted on \[Gamma]_{\[Mu]}"
]


VerificationTest[
	LCTensor[\[Mu], \[Nu], \[Rho], \[Sigma]] \[Gamma][\[Rho], \[Sigma]]// LC2Gamma5
,
	-2I \[Gamma][\[Mu], \[Nu]]** \[Gamma][5]
,
	TestID-> "LCTensor contracted on \[Gamma]_{\[Mu]\[Nu]}"
]


VerificationTest[
	LCTensor[\[Mu], \[Nu], \[Rho], \[Sigma]] \[Gamma][\[Nu], \[Rho], \[Sigma]]// LC2Gamma5
,
	6I \[Gamma][\[Mu]]** \[Gamma][5]
,
	TestID-> "LCTensor contracted on \[Gamma]_{\[Mu]\[Nu]\[Rho]}"
]


VerificationTest[
	LCTensor[\[Mu], \[Nu], \[Rho], \[Sigma]] \[Gamma][\[Mu], \[Nu], \[Rho], \[Sigma]]// LC2Gamma5
,
	24 I NonCommutativeMultiply@ \[Gamma][5]
,
	TestID-> "LCTensor contracted on \[Gamma]_{\[Mu]\[Nu]\[Rho]\[Sigma]}"
]
