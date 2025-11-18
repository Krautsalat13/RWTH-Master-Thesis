(* ::Package:: *)

(* ::Title:: *)
(*Initialization [Matchete`]*)


(* Generate warning for Mathematica versions earlier than 11.0.0 *)
If[$VersionNumber < 11.0, 
CellPrint[{
TextCell["Matchete` was developed for Mathematica 11.0.0 and later. 
Your current Mathematica version [" <> ToString@$Version <> "] might not be compatible. 
In case you are experiencing problems with Matchete`, please update your Mathematica version.",
"Text",Background->LightRed]
}]
]


(* Load the package *)
Module[{logoMatchete, packageDirectory},
	(* Avoid double loading the package *)
	If[MemberQ[$Packages,"Matchete`"],
		Print[Style["The package Matchete` is already loaded. Please restart the kernel to reload it.",RGBColor[.8,.4706,0.2573]]]
	];
	
	(* Set directory of Matchete package *)
	packageDirectory= DirectoryName[$InputFileName, 2]<>"/Package";
	logoMatchete= Import[FileNameJoin[{DirectoryName[$InputFileName],"LogoMatchete.eps"}]];
	
	(*Matchete version*)
	$MatcheteVersion= StringTrim[First@ StringSplit[
		Import@ FileNameJoin[{DirectoryName@ packageDirectory, "version"}], "\n"] ];
	
	(*The whole package contents are loaded here*)
	Check[
		(*Ensure Bar updates as first thing*)
		Get[FileNameJoin[{packageDirectory, "NCM.m"}]];
		,
		Print[Style["Loading failed!",RGBColor[.6,.0706,0.1373]]];
		Abort[]
	];

	(*Print the logo*)
	CellPrint@ ExpressionCell[
		Row[{
			Style[Graphics@ logoMatchete, Magnification-> 0.5],
			Row[{TextCell["v"<>$MatcheteVersion, Bold, 
				FontColor-> RGBColor[0, 0.3168,0.5928], 
				FontSize-> 15,
				FontFamily-> "Source Code Pro"]}, BaselinePosition->Scaled@.5]
		}, Spacer[7]]
	, CellMargins->{{70,5},{5,5}} ];

	Print[
		"by Javier Fuentes-Mart\[IAcute]n, Matthias K\[ODoubleDot]nig, Julie Pag\[EGrave]s, Anders Eller Thomsen, and Felix Wilsch \n",
		"Reference: ", Hyperlink["arXiv:2212.04510", "https://arxiv.org/abs/2212.04510"],"\n",
		"Website: ",Hyperlink["https://gitlab.com/matchete/matchete","https://gitlab.com/matchete/matchete"]
	];
];


Matchete`Matchete`PackagePrivate`AddAutoCompletion["LoadModel"][
	FileBaseName/@ FileNameTake/@ FileNames["*.m", FileNameJoin[{$MatchetePath, "Models"}]]];


(* Protect the symbols used in the package *)
(*
SetAttributes[
  Evaluate @ Flatten[Names /@ {"Matchete`*"}],
  {Protected, ReadProtected}
]
*)
