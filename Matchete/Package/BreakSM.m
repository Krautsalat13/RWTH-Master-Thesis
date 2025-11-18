(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`SFV Theory`*)


(* ::Subtitle:: *)
(*Feynman Rule finder for for Matchete packet.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(*PackageImport["GroupMagic`"]*)
PackageImport["Matchete`"]


(* ::Subsection:: *)
(*Exported*)


PackageExport["BreakingSM"]
PackageExport["BreakingSMEFT"]
PackageExport["DiagonaliseBrokenLagrangian"]


(* ::Subsection:: *)
(*Package Scope*)


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection:: *)
(*Exported*)


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Break the SM*)


Begin["Global`"]
BreakingSM[\[ScriptCapitalL]SM_]:= Module[{\[ScriptCapitalL]SMbroken,\[ScriptCapitalL]SMbrokenSimpl,preResult},
(*\[ScriptCapitalL]SM=LoadModel["SM", ModelParameters->{
	Global`W->Global`WW,
	Global`B->Global`BB,
	Global`q->Global`qq,
	Global`u->Global`uu,
	Global`d->Global`dd,
	Global`l->Global`ll,
	Global`e->Global`ee
}];*)

DefineGaugeGroup[Global`U1em,U1,Global`ge,Global`\[ScriptCapitalA],NiceForm->{"\[ScriptE]","\[ScriptCapitalA]"}, UFO$Options-><|"pdg"->22,"name"->"A"|>];
DefineCoupling[Global`gZ,EFTOrder->0,Indices->{},SelfConjugate->True,Symmetries->{},DiagonalCoupling->{},NiceForm->"\!\(\*SubscriptBox[\(g\), \(Z\)]\)"];
DefineCoupling[Global`CKM, Indices -> {Global`Flavor,Global`Flavor},NiceForm->"\!\(\*SubscriptBox[\(V\), \(CKM\)]\)"];

DefineCoupling[Global`v, EFTOrder->0, Indices -> {}, DiagonalCoupling->{},SelfConjugate->True,Symmetries->{}];

DefineCoupling[Global`me, Indices -> Global`Flavor, DiagonalCoupling -> True];
DefineCoupling[Global`mu, Indices -> Global`Flavor, DiagonalCoupling -> True];
DefineCoupling[Global`md, Indices -> Global`Flavor, DiagonalCoupling -> True];


DefineField[Global`\[ScriptCapitalZ],Vector,
	SelfConjugate->True,
	Mass->Heavy,
	UFO$Options -> <|"pdg" -> 23,"name"->"Z"|>
];

DefineField[Global`\[ScriptCapitalW],Vector,
	Charges->{Global`U1em[1]},
	Mass->Heavy,
	UFO$Options -> <|"pdg" -> 24,"name"->"W"|>
];


DefineField[Global`h,Scalar,
	SelfConjugate->True,
	Mass->Heavy,
	UFO$Options -> <|"pdg" -> 25|>
];

DefineField[Global`\[GothicU],Fermion,
	Indices->{Global`SU3c[fund],Global`Flavor},
	Charges->{Global`U1em[+2/3]},
	Mass->Heavy
];

DefineField[Global`\[GothicD],Fermion,
	Indices->{Global`SU3c[fund],Global`Flavor},
	Charges->{Global`U1em[-1/3]},
	Mass->Heavy
];

DefineField[Global`\[GothicE],Fermion,
	Indices->{Global`Flavor},
	Charges->{Global`U1em[-1]},
	Mass->Heavy
];

DefineField[Global`\[Nu],Fermion,
	Indices->{Global`Flavor},
	Mass->Heavy,
	Chiral->LeftHanded
];


SetSymmetryBreaking[
	{Global`SU2L,Global`U1Y},
	{Global`U1em}
];

RepresentationDecomposition[Global`SU2L[adj],  {ComplexSinglet,ComplexSinglet,Singlet}];
RepresentationDecomposition[Global`SU2L[fund],{Singlet,Singlet}];

FieldDecomposition[Global`W[Global`\[Mu],Global`a],{Bar@Global`\[ScriptCapitalW][Global`\[Mu]],Global`\[ScriptCapitalW][Global`\[Mu]],(Global`c\[Theta]*Global`\[ScriptCapitalZ][Global`\[Mu]]/Coupling[Global`gZ,{},0]+Global`s\[Theta]*Global`\[ScriptCapitalA][Global`\[Mu]]/Coupling[Global`ge,{},0])*Coupling[Global`gL,{},0]}];
FieldDecomposition[Global`B[Global`\[Mu]],{(-Global`s\[Theta]*Global`\[ScriptCapitalZ][Global`\[Mu]]/Coupling[Global`gZ,{},0]+Global`c\[Theta]*Global`\[ScriptCapitalA][Global`\[Mu]]/Coupling[Global`ge,{},0])*Coupling[Global`gY,{},0]}];



FieldDecomposition[Global`H[Global`i],{0,(Global`v+Global`h[])/Sqrt[2]}];




FieldDecomposition[Global`H[Global`i],{0,(Coupling[Global`v,{},0]+Global`h[])/Sqrt[2]}];


FieldDecomposition[Global`q[Global`a,Global`i,Global`p],{Global`CKM[Global`p, Global`r]*Global`\[GothicU][Global`a,Global`r], Global`\[GothicD][Global`a,Global`p]}];
FieldDecomposition[Global`u[Global`a,Global`p],{Global`\[GothicU][Global`a,Global`p]}];
FieldDecomposition[Global`d[Global`a,Global`p],{Global`\[GothicD][Global`a,Global`p]}];
FieldDecomposition[Global`l[Global`i,Global`p],{Global`\[Nu][Global`p],Global`\[GothicE][Global`p]}];
FieldDecomposition[Global`e[Global`p],{Global`\[GothicE][Global`p]}];

CGDecomposition[
	CG[gen@Global`SU2L@fund,{Global`J,Global`i,Global`j}],
	Normal@SparseArray[{
		(* indices below correspond to {J,i,j} above *)
		{1,1,2}->+(1/Sqrt[2])(*T+*),
		{2,2,1}->+(1/Sqrt[2])(*T-*),
		{3,1,1}->+(1/2),{3,2,2}->-(1/2)(*T3*)
	}]
];

Matchete`SSB`PackagePrivate`AddFStructDecomposition[Global`SU2L,Global`SU2L@fund];

CGDecomposition[
	CG[eps@Global`SU2L,{i,j}],
	Normal@SparseArray[{
		{1,2}->+1,
		{2,1}->-1
	}]
];

\[ScriptCapitalL]SMbroken=GreensSimplify[BrokenPhase[\[ScriptCapitalL]SM]];
\[ScriptCapitalL]SMbrokenSimpl=Matchete`SSB`PackagePrivate`PrepareLagrangianForMatching[\[ScriptCapitalL]SMbroken];
preResult = \[ScriptCapitalL]SMbrokenSimpl //ReplaceEffectiveCouplings;
preResult /.Global`c\[Theta]^2 ->1- Global`s\[Theta]^2 //Expand
]
End[]




(* ::Section:: *)
(*Break the SMEFT*)


Begin["Global`"]
BreakingSMEFT[LSMEFT_]:= Module[{LSMEFTcanonical,Rot,m,Rot2,Convention,LSMEFTcanonicalbroken,\[ScriptCapitalL]SMEFTbrokenSimpl,\[ScriptCapitalL]SMEFTbrokenSimplEff,\[ScriptCapitalL]SMEFTbrokenSimplEfforder,LSMEFTfinal},

DefineGaugeGroup[Global`U1em,U1,Global`ge,Global`\[ScriptCapitalA],NiceForm->{"\[ScriptE]","\[ScriptCapitalA]"}, UFO$Options-><|"pdg"->22,"name"->"A"|>];
DefineCoupling[Global`gZ,EFTOrder->0,Indices->{},SelfConjugate->True,Symmetries->{},DiagonalCoupling->{},NiceForm->"\!\(\*SubscriptBox[\(g\), \(Z\)]\)"];
DefineCoupling[Global`CKM, Indices -> {Global`Flavor,Global`Flavor},NiceForm->"\!\(\*SubscriptBox[\(V\), \(CKM\)]\)"];


DefineCoupling[Global`me, Indices -> Global`Flavor, DiagonalCoupling -> True];
DefineCoupling[Global`mu, Indices -> Global`Flavor, DiagonalCoupling -> True];
DefineCoupling[Global`md, Indices -> Global`Flavor, DiagonalCoupling -> True];


DefineCoupling[Global`v, EFTOrder->0, Indices -> {}, DiagonalCoupling->{},SelfConjugate->True,Symmetries->{}];

DefineField[Global`\[ScriptCapitalZ],Vector,
	SelfConjugate->True,
	Mass->Heavy,
	UFO$Options -> <|"pdg" -> 23,"name"->"Z"|>
];

DefineField[Global`\[ScriptCapitalW],Vector,
	Charges->{Global`U1em[1]},
	Mass->Heavy,
	UFO$Options -> <|"pdg" -> 24,"name"->"W"|>
];


DefineField[Global`h,Scalar,
	SelfConjugate->True,
	Mass->Heavy,
	UFO$Options -> <|"pdg" -> 25|>
];

DefineField[Global`\[GothicU],Fermion,
	Indices->{Global`SU3c[fund],Global`Flavor},
	Charges->{Global`U1em[+2/3]},
	Mass->Heavy
];

DefineField[Global`\[GothicD],Fermion,
	Indices->{Global`SU3c[fund],Global`Flavor},
	Charges->{Global`U1em[-1/3]},
	Mass->Heavy
];

DefineField[Global`\[GothicE],Fermion,
	Indices->{Global`Flavor},
	Charges->{Global`U1em[-1]},
	Mass->Heavy
];

DefineField[Global`\[Nu],Fermion,
	Indices->{Global`Flavor},
	Mass->Heavy,
	Chiral->LeftHanded
];


SetSymmetryBreaking[
	{Global`SU2L,Global`U1Y},
	{Global`U1em}
];

RepresentationDecomposition[Global`SU2L[adj],  {ComplexSinglet,ComplexSinglet,Singlet}];
RepresentationDecomposition[Global`SU2L[fund],{Singlet,Singlet}];



FieldDecomposition[Global`H[Global`i],{0,(Global`vT+(1+Global`cHkin)*Global`h[])/Sqrt[2]}];


LSMEFTcanonical = LSMEFT /.Coupling[Global`gs,{},0]:>(1-Coupling[Global`cHG,{},0]*Global`vT^2)* Coupling[Global`gs,{},0]/.Coupling[Global`gY,{},0]:>(1-Coupling[Global`cHB,{},0]*Global`vT^2)* Coupling[Global`gY,{},0]/.Coupling[Global`gL,{},0]:>(1-Coupling[Global`cHW,{},0]*Global`vT^2)* Coupling[Global`gL,{},0];



Rot ={{Global`c\[Theta],Global`s\[Theta]},{-Global`s\[Theta],Global`c\[Theta]}};
m = -ExplicitFlavor1/2*Global`vT^2*Coupling[Global`cHWB,{},0];
Rot2 = {{1,m},{m,1}};
Convention ={{1/Coupling[Global`gL,{},0] ,0},{0,1/Coupling[Global`gY,{},0]}};
{Global`W3Field,Global`BField} =Inverse[Convention] . Rot2 . Rot . {Global`\[ScriptCapitalZ][Global`\[Mu]]/Coupling[Global`gZ,{},0],Global`\[ScriptCapitalA][Global`\[Mu]]/Coupling[Global`ge,{},0]};

FieldDecomposition[Global`W[Global`\[Mu],Global`a],{Global`\[ScriptCapitalW][Global`\[Mu]],Bar@Global`\[ScriptCapitalW][Global`\[Mu]],Global`W3Field}];
FieldDecomposition[Global`B[Global`\[Mu]],{Global`BField}];


FieldDecomposition[Global`q[Global`a,Global`i,Global`p],{Global`CKM[Global`p, Global`r]*Global`\[GothicU][Global`a,Global`r], Global`\[GothicD][Global`a,Global`p]}];
FieldDecomposition[Global`u[Global`a,Global`p],{Global`\[GothicU][Global`a,Global`p]}];
FieldDecomposition[Global`d[Global`a,Global`p],{Global`\[GothicD][Global`a,Global`p]}];
FieldDecomposition[Global`l[Global`i,Global`p],{Global`\[Nu][Global`p],Global`\[GothicE][Global`p]}];
FieldDecomposition[Global`e[Global`p],{Global`\[GothicE][Global`p]}];

CGDecomposition[
	CG[gen@Global`SU2L@fund,{Global`J,Global`i,Global`j}],
	Normal@SparseArray[{
		(* indices below correspond to {J,i,j} above *)
		{1,1,2}->+(1/Sqrt[2])(*T+*),
		{2,2,1}->+(1/Sqrt[2])(*T-*),
		{3,1,1}->+(1/2),{3,2,2}->-(1/2)(*T3*)
	}]
];

Matchete`SSB`PackagePrivate`AddFStructDecomposition[Global`SU2L,Global`SU2L@fund];

CGDecomposition[
	CG[eps@Global`SU2L,{i,j}],
	Normal@SparseArray[{
		{1,2}->+1,
		{2,1}->-1
	}]
];
Print["Start breaking of Lagrangian"];
LSMEFTcanonicalbroken=GreensSimplify[BrokenPhase[LSMEFTcanonical]];
\[ScriptCapitalL]SMEFTbrokenSimpl=Matchete`SSB`PackagePrivate`PrepareLagrangianForMatching[LSMEFTcanonicalbroken];
\[ScriptCapitalL]SMEFTbrokenSimplEff =\[ScriptCapitalL]SMEFTbrokenSimpl //ReplaceEffectiveCouplings;
Global`cHkin = (Coupling[Global`cHBox,{},0]-1/4*Coupling[Global`cHD,{},0])*Coupling[Global`v,{},0]^2;
Global`vT = Sqrt[2*Coupling[Global`\[Mu]2,{},0]/Coupling[Global`\[Lambda],{},0]] + 3*Coupling[Global`\[Mu]2,{},0]^(3/2)/(Sqrt[2]*Coupling[Global`\[Lambda],{},0]^(5/2))*Coupling[Global`cH,{},0];
Print["Lagrangian is broken. Starting with expansion up to order 1."];
\[ScriptCapitalL]SMEFTbrokenSimplEfforder = \[ScriptCapitalL]SMEFTbrokenSimplEff /. Coupling[sym_Symbol,args___]/;StringTake[SymbolName[sym],1]=="c":>Coupling[sym,args]*order;
LSMEFTfinal =Normal[Series[\[ScriptCapitalL]SMEFTbrokenSimplEfforder,{order,0,1}]]/.order-> 1;
LSMEFTfinal /.Global`c\[Theta]^2 ->1- Global`s\[Theta]^2 //Expand

]
End[]



Begin["Global`"]
DiagonaliseBrokenLagrangian[Lag_]:= Module[{},
(*DefineCoupling[Global`CKM, Indices -> {Global`Flavor,Global`Flavor}];*)

DefineCoupling[Global`me, Indices -> Global`Flavor, DiagonalCoupling -> True];
DefineCoupling[Global`mu, Indices -> Global`Flavor, DiagonalCoupling -> True];
DefineCoupling[Global`md, Indices -> Global`Flavor, DiagonalCoupling -> True];

diagonal = {
Global`Yd[r_,s_] -> 1/Coupling[Global`v,{},0] Global`md[r] Delta[Index[r,Global`Flavor],Index[s,Global`Flavor]],
Bar@Global`Yd[r_,s_] -> 1/Coupling[Global`v,{},0] Global`md[r] Delta[Index[r,Global`Flavor],Index[s,Global`Flavor]],
Global`Ye[r_,s_] -> 1/Coupling[Global`v,{},0] Global`me[r] Delta[Index[r,Global`Flavor],Index[s,Global`Flavor]],
Bar@Global`Ye[r_,s_] -> 1/Coupling[Global`v,{},0] Global`me[r] Delta[Index[r,Global`Flavor],Index[s,Global`Flavor]],
Bar@Global`CKM[p_,r_]Global`Yu[p_,s_] -> 1/Coupling[Global`v,{},0]  Global`mu[r] Delta[Index[r,Global`Flavor],Index[s,Global`Flavor]],
Global`CKM[p_,r_]Bar@Global`Yu[p_,s_] -> 1/Coupling[Global`v,{},0] Global`mu[r] Delta[Index[r,Global`Flavor],Index[s,Global`Flavor]],
Bar@Global`CKM[r_, p_]Global`CKM[r_,s_]:>Delta[Index[p,Global`Flavor],Index[s,Global`Flavor]]
};

Lagdiag = Lag /.diagonal;

Lagdiagcont = Contract@Lagdiag;

Lagdiagcont = Lagdiagcont //. c___*(a___**DiracProduct[ops1___,Proj[-1],ops2___]**b__)+c___*(a___**DiracProduct[ops1___,Proj[1],ops2___]**b___):> c a**DiracProduct[ops1,ops2]**b //Contract;

ExplicitFlavor[Lagdiagcont]

]
End[]

