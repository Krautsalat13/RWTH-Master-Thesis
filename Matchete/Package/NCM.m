(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`NCM`*)


(* ::Subtitle:: *)
(*Contains the definitions of the NonCommutativeMultiply*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection:: *)
(*Exported*)


Unprotect@ NonCommutativeMultiply; (* overwrite NonCommutativeMultiply *)


PackageExport["Bar"]


(* ::Subsubsection:: *)
(*Internal*)


PackageScope["NonBarableQ"]
PackageScope["SetBarable"]
PackageScope["UnsetBarable"]


PackageScope["CommutativeQ"]
PackageScope["Commutative"]


PackageScope["Commutator"]
PackageScope["AntiCommutator"]


PackageScope["FermionTrace"]
PackageScope["CanonizeSpinorLines"]


PackageScope["OpenSpinChainQ"]


PackageScope["NCDot"]


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


NonCommutativeMultiply::usage = "NonCommutativeMultiply[a, b, ...] (alternatively a**b**...) denotes the Dirac product of spinors, Dirac matrices, etc.";


Bar::usage   = "Bar[field] returns the bar of a fermion field or the conjugate otherwise. Bar can also be applied to conjugate indices and charges.";


(* ::Subsubsection::Closed:: *)
(*Internal*)


NonBarableQ::usage  = "Specifies whether a term is trivial under bar or not.";
SetBarable::usage   = "Function for setting complex objects.";
UnsetBarable::usage = "Function for unsetting complex objects.";


CommutativeQ::usage      = "CommutativeQ[obj] returns True if obj is a comutative object and false otherwise.";
Commutative::usage       = "Commutative[obj] is an auxiliary head that can be used to treat the object obj temporarily as commutative.";


Commutator::usage     = "Commutator[X,Y] gives the commutator of X and Y, i.e. [X,Y]=X**Y-Y**X.";
AntiCommutator::usage = "AntiCommutator[X,Y] gives the anti-commutator of X and Y, i.e. {X,Y}=X**Y+Y**X.";


FermionTrace::usage        = "FermionTrace[expr] closes fermion lines in expr and performs Dirac gamma traces.";
CanonizeSpinorLines::usage = "CanonizeSpinorLines[expr] separates out nested spinor lines, breaking up the NCM accordingly.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Bar*)


(* ::Text:: *)
(*Bar denotes complex conjugation except for on fermion lines, where they are the Bar of the Dirac algebra*)


(* ::Subsubsection::Closed:: *)
(*Set rules for what is trivial under a Bar *)


(* ::Text:: *)
(*Default is no action*)


NonBarableQ[f_?NonBarableQ[x___]]:= And@@ NonBarableQ/@ {x};
NonBarableQ[f_@ x___]:= False;
NonBarableQ[_] := True; 
NonBarableQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold]= False;


(* ::Text:: *)
(*Objects with non-trivial behavior under Bar *)


NonBarableQ@ Alternatives[Complex, CG, Index]= False;
NonBarableQ@ _SparseArray= False;


NonBarableQ@ Alternatives[Coupling, DiracProduct, Field, FieldStrength,
	NonCommutativeMultiply, Transp, EvaOp]= False;


(* ::Subsubsection::Closed:: *)
(*General properties of Bar*)


Bar@ Bar@ x_:= x;


Bar@ x_?NonBarableQ:= x;


Bar[expr:Alternatives[_Plus, _Times, _Power, _Log, _List]]:= Bar/@expr;


Bar@ x_SparseArray:= Conjugate@ x;


Bar@ x_Complex:= Conjugate@ x;


(* ::Subsubsection::Closed:: *)
(*Action on Matchete symbols*)


Bar[expr_NonCommutativeMultiply]:= Bar/@ Reverse@ expr;
Bar@ x_DiracProduct:= Reverse@ x/. {Gamma5-> -Gamma5, Proj@ s_:> Proj[-s], g_GammaM:> Reverse@ g};


Bar@ f:Field[label_, Scalar|Vector[_]|Ghost, __]/; $FieldAssociation[label, SelfConjugate]:=f;
Bar@ f:Field[label_, Fermion, __]/; $FieldAssociation[label, SelfConjugate]:= 
	Transp[f]**DiracProduct[GammaCC]; (* Bar for Majorana Fermions *)
Bar@ FieldStrength[label_, lind_, ginds_, cdinds_]/; $FieldAssociation[label, SelfConjugate]:= 
	FieldStrength[label, lind, Bar/@ ginds, cdinds];
Bar@ c:Coupling[label_, __]/; $CouplingAssociation[label, SelfConjugate]:= c;
Bar@ Coupling[label_, ind_List,x__]/; (Head[$CouplingAssociation[label, SelfConjugate]]===List):= Coupling[label, Permute[ind, $CouplingAssociation[label, SelfConjugate]],x];


(*Special case of charge conjugation before Bar*)
Transp/: Bar@ Transp@ f_:= Transp@ Bar@ f;


LF/: Bar@ lf_LF:= lf;
LCTensor/: Bar@ lc_LCTensor:= lc;


Bar@ Index[x_, Lorentz]:= Index[x, Lorentz];


(* ::Subsubsection::Closed:: *)
(*For setting objects complex*)


SetBarable[x_, y__]  := (SetBarable@ x; SetBarable@ y;);
SetBarable[x_] := (NonBarableQ@ x ^= False;);


(* ::Text:: *)
(*And unseting them for removal*)


UnsetBarable[x_Symbol] := (x/: NonBarableQ@ x =.);
UnsetBarable[x:_Symbol@_] := (Evaluate@ Head@ x/: NonBarableQ@ x =.);


(* ::Section:: *)
(*NonCommutativeMultiply functionality*)


(* ::Subsection:: *)
(*Commutative check*)


(* ::Text:: *)
(*Hard-code commutative properties of objects to prevent iterative calls to all subparts: greatly improves performance *)


CommutativeQ@ _DiracProduct= False;
CommutativeQ@ Field[_, Fermion| Ghost, __]= False;


CommutativeQ@ Field[_, Scalar| _Vector, __]= True;
CommutativeQ@ _FieldStrength= True;
CommutativeQ@ _Coupling= True;
CommutativeQ@ _CG= True;
CommutativeQ@ _LF= True;


(* ::Text:: *)
(*Default assumption is that everything is commutative*)


(* Memoization speeds up DeriveSubstitutions by a factor of ~2 or so *)
(* Memoization seems to have very poor performance in EOM and GreensSimplify creating 100,000s of Downvalues *)
(*CommutativeQ[x_NonCommutativeMultiply]:= (CommutativeQ[x]= ClosedSpinChainQ[x]);
CommutativeQ[f_?CommutativeQ[x___]]:= (CommutativeQ[f[x]]= And@@ CommutativeQ/@ {x});*)

CommutativeQ@ x_NonCommutativeMultiply:= ClosedSpinChainQ@ x;
CommutativeQ[f_?CommutativeQ[x___]] := And@@ CommutativeQ/@ {x};

CommutativeQ[f_[x___]] := False;
CommutativeQ[_] := True; 


(* ::Text:: *)
(*Define the non-commutative objects*)


CommutativeQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold] = False;


(* ::Text:: *)
(*Commutative is a head that can be used to temporarily treat an object as commutative*)


CommutativeQ@ Commutative@ _ ^= True; 


(* ::Subsection:: *)
(*Properties of NonCommutativeMultiply*)


(* ::Text:: *)
(*Remove attributes problematic for pattern matching*)


ClearAttributes[NonCommutativeMultiply, {Flat, OneIdentity}];


(* ::Text:: *)
(*Set definitions for the non-commutative product*)


(* ::Subsubsection::Closed:: *)
(*General properties*)


(* ::Text:: *)
(*NCM contractions*)


NonCommutativeMultiply[a___, NonCommutativeMultiply[b__], c___]:= NonCommutativeMultiply[a, b, c] /;FreeQ[List@b,Fermion]
NonCommutativeMultiply[] = 1;
NonCommutativeMultiply@NonCommutativeMultiply@x___:=NonCommutativeMultiply@x;


(* Extracting commuting objects *)
NonCommutativeMultiply[a___, b_?CommutativeQ, c___]:= b * NonCommutativeMultiply[a, c];
NonCommutativeMultiply[a___, b_?CommutativeQ * x_, c___]:= b * NonCommutativeMultiply[a, x, c];
(* Distributivity *)
NonCommutativeMultiply[a___, b_Plus, c___]:= NonCommutativeMultiply[a, #, c] & /@ b 


(*On matrix*)
NCDot[a_, b_]/;(Length@Dimensions@a===2 && Length@Dimensions@b===2 && (Dimensions[a][[2]])===(Dimensions[b][[1]])):=
	Inner[NonCommutativeMultiply,a,b,Plus];
NCDot[a_,b_,c__]:=NCDot[NCDot[a,b],c];
NCDot[0,0]:=0;


(* ::Subsubsection::Closed:: *)
(*Contraction of spin chains*)


NonCommutativeMultiply[a___, x:NonCommutativeMultiply@ b___, c___]:=
	If[ClosedSpinChainQ@ x,
		x NonCommutativeMultiply[a, c],
		NonCommutativeMultiply[a, b, c]
	];


(* ::Subsubsection::Closed:: *)
(*Canonical ordering of operators with transpose*)


(*No terms with two fields transposed*)
NonCommutativeMultiply[A: Transp@Field[_,Fermion,___] , B___/;FreeQ[List@B,Fermion] ,C: Transp@Bar@Field[_,Fermion,___]]:= - Transp[C]** Transp[B] ** Transp[A] 
(*No complex fermion transpose on the left*)
NonCommutativeMultiply[A: Transp@Field[labelA_,Fermion,___]/;!$FieldAssociation[labelA][SelfConjugate], B___/;FreeQ[List@B,Fermion], C: Field[labelC_,Fermion,___]/;$FieldAssociation[labelC][SelfConjugate]]:= 
	- Transp[C]** Transp[B] ** Transp[A] 
(*If there must be a transpose, and no derivatices, choose canonical ordering of labels*)
NonCommutativeMultiply[A:Transp@(Field[labelA_,Fermion,_,{}])/; $FieldAssociation[labelA][SelfConjugate], B___/;FreeQ[List@B,Fermion], C: Field[labelC_,Fermion,_,{}]/;$FieldAssociation[labelC][SelfConjugate]]:= 
	- Transp[C]** Transp[B] ** Transp[A]  /;!OrderedQ[{labelA,labelC}] 
NonCommutativeMultiply[A:Transp@(Field[labelA_,Fermion,_,{}])/;!$FieldAssociation[labelA][SelfConjugate], B___/;FreeQ[List@B,Fermion], C: Field[labelC_,Fermion,_,{}]/;!$FieldAssociation[labelC][SelfConjugate]]:= 
	- Transp[C]** Transp[B] ** Transp[A]  /;!OrderedQ[{labelA,labelC}]
NonCommutativeMultiply[A:Bar@(Field[labelA_,Fermion,_,{}]), B___/;FreeQ[List@B,Fermion], C: Transp@Bar@(Field[labelC_,Fermion,_,{}]) ]:=
	- Transp[C]** Transp[B] ** Transp[A]  /;!OrderedQ[{labelA,labelC}]


(*If same field and same number of derivatives, place fields according to indices*)
(*NonCommutativeMultiply[A:Bar@(Field[label_,Fermion,indA:Except[_Pattern],derA_]), B___/;FreeQ[List@B,Fermion], C: Transp@Bar@(Field[label_,Fermion,indC:Except[_Pattern],derC_]) ]/;(Length[derA]===Length[derC]):=
	- Transp[C]** Transp[B] ** Transp[A]  /;!OrderedQ[{indA,indC}/.{Index[_,Flavor]->Nothing,Index[_Pattern,_]->Nothing}]*)


(*Place field with most derivative on the right*)
NonCommutativeMultiply[A:Transp@Field[_,Fermion,_,Lind1_], B___/;FreeQ[List@B,Fermion], C: Field[_,Fermion,_,Lind2_]]:=
	- Transp[C]** Transp[B] ** Transp[A]  /;(Length[Lind1]>Length[Lind2]);
NonCommutativeMultiply[A:Bar@Field[_,Fermion,_,Lind1_], B___/;FreeQ[List@B,Fermion], C: Transp@Bar@Field[_,Fermion,_,Lind2_]]:=
	- Transp[C]** Transp[B] ** Transp[A]  /;(Length[Lind1]>Length[Lind2]);


(*Place field with EoM head on the right*)
NonCommutativeMultiply[EoM@ Transp[f1:Field[_, Fermion, _, {}]], d_DiracProduct, f2:Field[_, Fermion, _, {}]]:=
	- Transp@ f2 ** Transp@ d ** EoM@ f1;
NonCommutativeMultiply[EoM@ Bar[f1:Field[_, Fermion, _, {}]], d_DiracProduct, Transp[f2:Bar@ Field[_, Fermion, _, {}]]]:=
	- f2 ** Transp@ d ** EoM@ Transp@ Bar@ f1;


(* ::Subsection:: *)
(*Utility functions*)


Commutator[x_, y_]     := x ** y - y ** x;
AntiCommutator[x_, y_] := x ** y + y ** x;


(* ::Subsection:: *)
(*Non-commutative symbols used in Matchete should be declared here*)


(* ::Text:: *)
(*From DiracAlgebra.m*)


CommutativeQ@ Alternatives[GammaM, Gamma5, GammaCC]= False;


(* ::Text:: *)
(*From FunctionalTools.m*)


(* SetNonCommutative[TransposeThisSpinChain] *)
(* For some reason the lie above does not work, but the one below does... *)
CommutativeQ[TransposeThisSpinChain]=False


(* ::Section:: *)
(*Active manipulation*)


(* ::Subsection:: *)
(*Organize fermion lines*)


(* ::Subsubsection::Closed:: *)
(*Fermion trace*)


(* ::Text:: *)
(*Function for closing fermion lines in traces as appear in some supertraces. *)


FermionTrace::cantclose = "The fermion line cannot be closed";
FermionTrace::multopens = "Multiple open fermion lines encountered in one term.";


OpenSpinChainQ= Not@* ClosedSpinChainQ;


FermionTrace@ expr_:= Module[{out= BetterExpand@ expr},	
	(*Trace terms in the sums separately*)
	If[Head@ out === Plus, Return[FermionTrace/@ out];]; 
	
	Switch[Count[{out}, _NonCommutativeMultiply? OpenSpinChainQ, All]
		,0 ,
			4 out
		,1 ,
			out/. x_NonCommutativeMultiply? OpenSpinChainQ:> TraceSpin@ x// ContractMetric
		,_ ,
			Message[FermionTrace::multopens];
			Abort[];
	]
];


TraceSpin@ spinChain_NonCommutativeMultiply:= Module[{spinCount= CumulativeSpin@ spinChain},
	If[Last@ spinCount =!= 0,
		Message[FermionTrace::cantclose];
		Abort[];	
	];
	(*Check if the fermion line can be closed through the trace *)
	If[MemberQ[spinCount, -1],
		(*Sign is from Grassmanian variables*)
		Return[-RotateLeft[spinChain, FirstPosition[spinCount, -1]]];
	];
	(*Otherwise do a Dirac trace*)
	DiracTrace@ spinChain
];


(* ::Subsubsection::Closed:: *)
(*Spinor line canonization*)


(* ::Text:: *)
(*Function for separating out nested spinor lines, breaking up the NCM accordingly*)


CanonizeSpinorLines::lopens = "The Spinor line has >1 L-open spinors.";
CanonizeSpinorLines::ropens = "The Spinor line has >1 R-open spinors.";


CanonizeSpinorLines@ fline_NonCommutativeMultiply:= Module[
		{cumSpinCount, max, min, l, r, expr},
	cumSpinCount= CumulativeSpin@ fline;
	min= Min@ cumSpinCount;
	If[min < -1,
		Message[CanonizeSpinorLines::lopens];
		Abort[];
	];
	If[Last@ cumSpinCount> 1,
		Message[CanonizeSpinorLines::ropens];
		Abort[];
	];
	
	If[min ===-1,
		cumSpinCount[[;; First@ FirstPosition[cumSpinCount, -1] -1]]--;
	]; 
	max= Max@ cumSpinCount;
		
	(*If a closed spinor line appear it is separated out*)
	{l, r}= First@ SequencePosition[cumSpinCount, {max..}, Overlaps-> False]+ {0, 1};
	If[r <= Length@ fline,
		Return[fline[[l;;r]] CanonizeSpinorLines[fline[[;; l-1]]** fline[[r+ 1;;]] ] ];
	];

	fline
];


CanonizeSpinorLines@ expr_:= expr//. x_NonCommutativeMultiply:> CanonizeSpinorLines@ x;
