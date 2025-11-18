(* ::Package:: *)

(* ::Title:: *)
(*LoopIntegration test unit*)


(* ::Section:: *)
(*Tests*)


If[!MemberQ[$ContextPath, "Matchete`PackageScope`"], PrependTo[$ContextPath,"Matchete`PackageScope`"];]


(* ::Subsubsection:: *)
(*LoopIntegrate*)


VerificationTest[
	LoopIntegrate[Prop[m]^2 Prop[0]^2],
	-(I/(m^4 \[Epsilon])) + I LF[{m},{2,2}],
	TestID->"Integration of \!\(\*FractionBox[\(1\), \(\*SuperscriptBox[\(p\), \(2\)] \*SuperscriptBox[\((\*SuperscriptBox[\(p\), \(2\)] - \*SuperscriptBox[\(m\), \(2\)])\), \(2\)]\)]\)"
]


VerificationTest[
	LoopIntegrate[Prop[m]^2 Prop[0]^2] //EvaluateLoopFunctions //Expand,
	-(I/m^4)(1/\[Epsilon]+Log[\[Mu]bar2/m^2]+2) //Expand,
	TestID->"Evaluation of loop function from \!\(\*FractionBox[\(1\), \(\*SuperscriptBox[\(p\), \(2\)] \*SuperscriptBox[\((\*SuperscriptBox[\(p\), \(2\)] - \*SuperscriptBox[\(m\), \(2\)])\), \(2\)]\)]\)"
]


VerificationTest[
	LoopIntegrate[Prop[m]^2 Prop[0]^-2] ,
	3m^4 I/\[Epsilon] + I LF[{m},{2,-2}],
	TestID->"Integration of \!\(\*FractionBox[SuperscriptBox[\(p\), \(2\)], SuperscriptBox[\((\*SuperscriptBox[\(p\), \(2\)] - \*SuperscriptBox[\(m\), \(2\)])\), \(2\)]]\)"
]


VerificationTest[
	LoopIntegrate[Prop[m]^2 Prop[0]^-2] //EvaluateLoopFunctions //Expand,
	 I m^4  (3/\[Epsilon] + 3 Log[\[Mu]bar2/m^2] + 2) //Expand,
	TestID->"Evaluation of loop function from \!\(\*FractionBox[SuperscriptBox[\(p\), \(2\)], SuperscriptBox[\((\*SuperscriptBox[\(p\), \(2\)] - \*SuperscriptBox[\(m\), \(2\)])\), \(2\)]]\)"
]


(* ::Subsubsection:: *)
(*SimplifyMassFunction*)


VerificationTest[
	Block[{lf= LF[{a, b}, {2, 1, 1}], rulesApplied},
		rulesApplied= Matchete`LoopIntegration`PackagePrivate`LFtoTempLF@
			Matchete`LoopIntegration`PackagePrivate`SimpTempLFRules@ lf/.
			Matchete`LoopIntegration`PackagePrivate`TempLF-> LF;
		lf - rulesApplied// EvaluateLoopFunctions// Simplify 
	]
,
	0
,
	TestID-> "SimpTempLFRules: pole treatment"
]


VerificationTest[
	LF[{M1, M3}, {1, 1, 1}]+ LF[{M1, M3}, {2, 1, 0}]- LF[{M3, M1}, {2, 1, 0}]// SimplifyMassFunction
,
	2 LF[{M1, M3}, {2, 1, 0}]
,
	TestID-> "SimplifyMassFunction: 2 mass simple sum 1"
]


VerificationTest[
	LF[{M1, M3}, {1, 1, 1}]- LF[{M1, M3}, {2, 1, 0}]+ LF[{M3, M1}, {2, 1, 0}]// SimplifyMassFunction
,
	2 LF[{M1, M3}, {1, 2, 0}]
,
	TestID-> "SimplifyMassFunction: 2 mass simple sum 2"
]


VerificationTest[
	LF[{Md, Mq}, {2, 2, 0}]+ LF[{Md, Mq}, {3, 1, 0}]- LF[{Md, Mq}, {3, 2, -1}]-
	LF[{Md, Mq}, {4, 1, -1}]+ LF[{Mq, Md}, {3, 1, 0}]- LF[{Mq, Md}, {3, 2, -1}]- 
	LF[{Mq, Md}, {4, 1, -1}] // SimplifyMassFunction
,
	(*The result of %// EvaluateLoopFunctions// Simplify*)
	1/6 Mq^-2 Md^-2
,
	TestID-> "SimplifyMassFunction: 2 mass complicated sum full reduction"
]


VerificationTest[
	Block[{lf8= -LF[{Coupling[Mq, {}, 0], Coupling[Mu, {}, 0]}, {2, 2, 0}]-
		LF[{Coupling[Mq, {}, 0], Coupling[Mu, {}, 0]}, {3, 1, 0}]+ 
		LF[{Coupling[Mq, {}, 0], Coupling[Mu, {}, 0]}, {3, 2, -1}]+ 
		LF[{Coupling[Mq, {}, 0], Coupling[Mu, {}, 0]}, {4, 1, -1}]+
		LF[{Coupling[Mu, {}, 0], Coupling[Mq, {}, 0]}, {3, 1, 0}]+
		LF[{Coupling[Mu, {}, 0], Coupling[Mq, {}, 0]}, {3, 2, -1}]-
		3 LF[{Coupling[Mu, {}, 0], Coupling[Mq, {}, 0]}, {4, 1, -1}]+
		2 LF[{Coupling[Mu, {}, 0], Coupling[Mq, {}, 0]}, {5, 1, -2}]},
		{SimplifyMassFunction@ lf8 -lf8// EvaluateLoopFunctions// Simplify, Length@ SimplifyMassFunction@ lf8}
	]
,
	{0, 2}
,
	TestID-> "SimplifyMassFunction: 2 mass complicated sum partial reduction"
]


(* ::Section:: *)
(*Loading the package*)


(* ::Text:: *)
(*Un-comment and run first to see the tests directly in this file. *)
(*Comment when finished and add file name to the list testfiles in  the script "RunTests.wls" .*)


(*	PrependTo[$Path, ParentDirectory[NotebookDirectory[]]];
	<< Matchete`*)
