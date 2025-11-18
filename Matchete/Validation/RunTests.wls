#!/usr/bin/env wolframscript

Print["  Running \"Matchete\" tests..."(*,"\n"*)];

(* Path of the directory from which the script is started (git root directory by default) . *)

PacletDirectoryLoad@Directory[];

Block[{Print},
  Get["Matchete`"];
]

#Print["  Using Mathematica: ",$VersionNumber];

Module[
    {testfiles,report,time,results,failIdx,commitpass=0, j},

    testfiles = {
		"GroupMagic.wl",
		"CovariantDerivative.wl",
		"DiracAlgebra.wl",
		"Definitions.wl",
		"LoopIntegration.wl",
		"NCM.wl"};

    Do[
    	report=TestReport[
		FileNameJoin[{Directory[],"Validation/Tests", testfiles[[j]]}]
    	];
    	time=Round[
      		QuantityMagnitude@report["TimeElapsed"],
    		0.001
    	];
    	results=report["TestResults"];
    	failIdx=report["TestsFailedIndices"];

    	Print["\n","  ", testfiles[[j]]," : ",Length[results]," tests run in ",time," seconds."];

    	If[
        		TrueQ@report["AllTestsSucceeded"]
        		,
        		Print["     All tests succeeded!"];
        		,
        		Print["     ",Length[failIdx]," tests failed!"];
        		Do[
            		Print["        ID: ",i," / ",results[i]["TestID"]," -> ",results[i]["Outcome"]],
            		{i,failIdx}
       		 ];
		Print["   "]
        		Quit[1] (* exit code if any test fails *)
   	],
	{j,Length[testfiles]}
    ]
    Print["  "]
    Quit[commitpass] (* exit code with 0 if all test pass *)
]
