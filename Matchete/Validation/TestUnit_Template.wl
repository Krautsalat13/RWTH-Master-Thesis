(* ::Package:: *)

(* ::Title:: *)
(*Template file for unit test*)


(* ::Section:: *)
(*Loading the Matchete` package*)


(* ::Text:: *)
(*Un-comment to run the tests directly in this file. *)
(*However, it must commented out to work with the pre-commit automatic check. *)


(*
	PrependTo[$Path,ParentDirectory[NotebookDirectory[]]];
	<< Matchete`
*)


(* ::Section:: *)
(*Add all test like this:*)


(* ::Subsubsection:: *)
(*test 1*)


VerificationTest[
	True==False,
	False,
	TestID->"some test"
]


(* ::Subsubsection:: *)
(*test 2*)


(* ::Text:: *)
(*Un-comment to see Failed test. (commented for commit...)*)


(*VerificationTest[
	True==False,
	True,
	TestID->"some other test"
]*)


(* ::Subsubsection:: *)
(*test 3*)


(* ::Text:: *)
(*Un-comment to see MessageFailure. (commented for commit...)*)


(*VerificationTest[
	1/0,
	ComplexInfinity,
	TestID->"infinity test"
]*)


(* ::Subsubsection:: *)
(*test 4*)


VerificationTest[
	Metric[Index[a,Lorentz],Index[b,Lorentz]]Metric[Index[b,Lorentz],Index[c,Lorentz]],
	Metric[Index[a,Lorentz],Index[c,Lorentz]],
	TestID->"metric"
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
