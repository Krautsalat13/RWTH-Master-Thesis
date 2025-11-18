(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Definitions`*)


(* ::Subtitle:: *)
(*Field and coupling definitions for Matchete paclet.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection:: *)
(*Exported*)


PackageExport["Field"]
PackageExport["DefineField"]
PackageExport["RemoveField"]
PackageExport["ResetFields"]
PackageExport["GetFields"]
PackageExport["GetFieldsByProperty"]


PackageExport["Scalar"]
PackageExport["Fermion"]
PackageExport["Vector"]
PackageExport["Ghost"]


PackageExport["Indices"]
PackageExport["Charges"]
PackageExport["SelfConjugate"]
PackageExport["Mass"]
PackageExport["Heavy"]
PackageExport["Light"]
PackageExport["Type"]
PackageExport["Chiral"]
PackageExport["LeftHanded"]
PackageExport["RightHanded"]
PackageExport["GoldstoneBoson"]
PackageExport["UFO$Options"]


PackageExport["Symmetries"]
PackageExport["SymmetricPermutation"]
PackageExport["AntisymmetricPermutation"]


PackageExport["Coupling"]
PackageExport["DefineCoupling"]
PackageExport["RemoveCoupling"]
PackageExport["GetCouplings"]
PackageExport["ResetCouplings"]
PackageExport["DiagonalCoupling"]


PackageExport["GetFlavorIndices"]
PackageExport["IndexDimension"]
PackageExport["DefineFlavorIndex"]
PackageExport["RemoveFlavorIndex"]
PackageExport["ResetFlavorIndices"]


PackageExport["IndexAlphabet"]


PackageExport["Group"]
PackageExport["Abelian"]
PackageExport["Representations"]
PackageExport["ClebschGordanCoefficients"]
PackageExport["RepresentationProperties"]
PackageExport["Dimension"]
PackageExport["DynkinCoefficients"]
PackageExport["Reality"]
PackageExport["ClebschGordanIndices"]


PackageExport["U1"]
PackageExport["fund"]
PackageExport["tFundf"]
PackageExport["adj"]
PackageExport["gen"]
PackageExport["fStruct"]
PackageExport["eps"]
PackageExport["dSym"]


PackageExport["FieldStrength"]
PackageExport["FS"]


PackageExport["GetGaugeGroups"]
PackageExport["FundAlphabet"]
PackageExport["AdjAlphabet"]
PackageExport["DefineGaugeGroup"]
PackageExport["RemoveGaugeGroup"]
PackageExport["ResetGaugeGroups"]


PackageExport["GetGlobalGroups"]
PackageExport["DefineGlobalGroup"]
PackageExport["RemoveGlobalGroup"]
PackageExport["ResetGlobalGroups"]


PackageExport["ResetAll"]


PackageExport["PlusHc"]


PackageExport["FreeLag"]
PackageExport["LoadModel"]
PackageExport["ParentModel"]
PackageExport["ParameterDefault"]
PackageExport["SetModelOptions"]
PackageExport["GetModels"]
PackageExport["ModelParameters"]


(* ::Subsubsection:: *)
(*Internal*)


PackageExport["FieldGenerators"]


PackageScope["SymmetryOverride"]
PackageScope["CouplingAssociationEntry"]


PackageScope["GroupFromRep"]
PackageScope["GroupFromInd"]
PackageScope["FieldsFromGroup"]
PackageScope["CouplingsFromGroup"]


PackageScope["$FieldAssociation"]
PackageScope["$CouplingAssociation"]
PackageScope["$FlavorIndices"]
PackageScope["$GaugeGroups"]
PackageScope["GetGaugeGroupByProperty"]
PackageScope["$GlobalGroups"]


PackageScope["KinOpLagrangian"]


PackageScope["$DropDiagonalCouplings"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


Field::usage               = "Field[label,type,{indices},{CDerivs}] denotes a field with given label, type and list of {indices}. The argument {CDerivs} is the list of Lorentz indices of the covariant derivatives acting on the field.";
DefineField::usage         = "DefineField[label, type] defines a new field of given type and label. The type can be Scalar,Fermion or Vector. The default options for this function are: Indices->{}, Charges->{}, SelfConjugate->False, Chiral -> False, Mass->Heavy, NiceForm-> Default, UFO$Options -> <||>.";
RemoveField::usage         = "RemoveField[alias] undefines the field with label alias together with its corresponding mass coupling.";
ResetFields::usage         = "Undefines all fields (excluding gauge fields) and field masses."
GetFields::usage           = "GetFields[] returns an association of all fields that have been defined by DefineField and their properties. GetFields[field] returns an association with all the properties of an already defined field.";
GetFieldsByProperty::usage = "GetFieldsByProperty[properties] returns all fields matching the property, list of properties or sequence of properties determined by the argument.";


Scalar::usage  =    "Field[\[Phi],Scalar,{indices},{CDerivs}] denotes a scalar field \[Phi].";
Vector::usage  =    "Field[A,Vector[Index[\[Mu],Lorentz]],{indices},{CDerivs}] denotes a vector field \!\(\*SubscriptBox[\(A\), \(\[Mu]\)]\).";
Fermion::usage =    "Field[\[Psi],Fermion,{indices},{CDerivs}] denotes a fermion field \[Psi].";
Ghost::usage   =    "Field[c,Ghost,{indices},{CDerivs}] denotes a ghost field c.";


Indices::usage       = "Indices is an option for the routine DefineField that specifies a list of representations under which the field transforms. By default, Indices is an empty list {}.";
Charges::usage       = "Charges is an option for the routine DefineField that specifies a list of U(1) groups and corresponding charges of the field. By default, Charges is an empty list {}.";
SelfConjugate::usage = "SelfConjugate is an option for the routine DefineField that specifies whether the field is self-conjugated. By default, SelfConjugate is False.";
Mass::usage          = "Mass is an option for the routine DefineField that specifies the mass and EFT scaling of the field. The allowed values are either Heavy|Light|0, {Light,0}, {Heavy|Light,mlabel}, or {Heavy|Light,mlabel,{Indices}} where mlabel is a custom label for the mass and Indices must be a subset of the field indices. The default value is Mass->Heavy. If no mlabel is given, the default value Mflabel, with flabel being the label of the corresponding field, is used.";
Heavy::usage         = "Heavy specifies that the mass of a field is at the UV scale.";
Light::usage         = "Light specifies that the mass of a field is at the IR scale.";
Type::usage          = "Type specifies the type of a field. It can be Scalar, Vector or Fermion.";
Chiral::usage        = "Chiral is an option for the routine DefineField that specifies if a fermion field is LeftHanded or RightHanded.";
LeftHanded::usage    = "LeftHanded specifies that a field is of left-handed chirality."
RightHanded::usage   = "RightHanded specifies that a field is of right-handed chirality."
GoldstoneBoson::usage= "GoldstoneBoson indicates if the field is a GB associated with a heavy vector field."


Symmetries::usage               = "Symmetries is an option that takes a list of symmetries. The possible option values are: SymmetricIndices[n1, n2, ...], AntisymmetricIndices[n1, n2, ...], SymmetricPermutation[n1, n2, ...] or AntisymmetricPermutation[n1, n2, ...]." ;
SymmetricPermutation::usage     = "SymmetricPermutation is a possible value for the option Symmetries. SymmetricPermutation[n1, n2, ...] denotes symmetry under a particular permutation of indices." ;
AntisymmetricPermutation::usage = "AntisymmetricPermutation is a possible value for the option Symmetries. AntisymmetricPermutation[n1, n2, ...] denotes antisymmetry under a particular permutation of indices." ;


EFTOrder::usage = "EFTOrder is an option for various routines specifying the targeted order in power-counting, as well as defining the power-counting of couplings.";


Coupling::usage         = "Coupling[label,{indices},EFTOrder] denotes the coupling label with a set of indices and an EFTOrder."
DefineCoupling::usage   = "DefineCoupling[label] defines a new coupling with a given label. The default options for this function are: EFTOrder -> 0, Indices->{}, SelfConjugate->False, Symmetries-> {}, DiagonalCoupling -> False, , NiceForm-> Default.";
RemoveCoupling::usage   = "RemoveCoupling[alias] undefines the coupling with label alias.";
GetCouplings::usage     = "GetCouplings[] returns an association of all couplings that have been defined by DefineCoupling and their properties.GetCouplings[coupling] returns an association with all the properties of an already defined coupling.";
ResetCouplings::usage   = "Undefines all couplings (excluding the field masses and the gauge couplings)."
DiagonalCoupling::usage = "DiagonalCoupling is an option for the routine DefineCoupling. Allowed values are a list of booleans or an empty list. The n-th boolean specifies whether the n-th flavor index is a diagonal in flavor space. By default all flavor indices are assumed to be non-diagonal.

DiagonalCoupling is also the corresponding key in the association returned by GetCouplings[]."


GetFlavorIndices::usage   = "GetFlavorIndices[] returns an association of all flavor indices that have been defined by DefineFlavorIndex and their dimension. GetFlavorIndices[GaugeGroup] returns an association with the dimension of an already defined flavor index.";
IndexDimension::usage     = "IndexDimension is a key for the association returned by GetFlavorIndices.";
DefineFlavorIndex::usage  = "DefineFlavorIndex[IndexName, Dimension] defines a flavor index with given IndexName and Dimension. This routine has one optional argument, IndexAlphabet, to define the printing alphabet of the representation.";
RemoveFlavorIndex::usage  = "RemoveFlavorIndex[alias] undefines the flavor index with given alias.";
ResetFlavorIndices::usage = "Undefines all flavor indices.";


IndexAlphabet::usage             = "IndexAlphabet is an option for the routine DefineRepresentation and DefineFlavorIndex that specifies a printing alphabet for the representation or flavor indices.";


Group::usage                     = "Group is a key for the associations returned by GetGaugeGroups and RepresentationProperties.";
Abelian::usage                   = "Abelian is a key for the association returned by GetGaugeGroups.";
Representations::usage           = "Representations is a key for the association returned by GetGaugeGroups.";
ClebschGordanCoefficients::usage = "ClebschGordanCoefficients is a key for the association returned by GetGaugeGroups.";
RepresentationProperties::usage  = "RepresentationProperties[RepName] return the main properties of the representation RepName."
Dimension::usage                 = "Dimension is a key for the association returned by RepresentationProperties."
DynkinCoefficients::usage        = "Representations is a key for the association returned by RepresentationProperties."
Reality::usage                   = "Reality is a key for the association returned by RepresentationProperties."
ClebschGordanIndices::usage      = "ClebschGordanIndices[CGname] return the representation of the indices (ordered) for the Clebsch-Gordan coefficient CGname."


U1::usage      = "Standard name for Abelian U(1) group.";
fund::usage    = "Fundamental representations of groups are refered to by group@ fund."
tFundf::usage  = "Standard head for the Clebsch-Gordan coefficient, CG[tFundf@ group, {i, j, A, B}] associated with the commonly occuring combination \!\(\*SubscriptBox[SuperscriptBox[\(T\), \(Ci\)], \(j\)]\)\!\(\*SuperscriptBox[\(f\), \(CAB\)]\), where T is the generator of the fundamental representation of the group.";
adj::usage     = "Adjoint representations of groups are refered to by group@ adj."
gen::usage     = "Standard head for generator Clebsch-Gordan coefficients. E.g., CG[gen@ rep, {A, i, j}] for \!\(\*SubscriptBox[SuperscriptBox[\(T\), \(Ai\)], \(j\)]\), where A is an adjoint index and i an index of the representation rep.";
fStruct::usage = "Standard head for structure constant Clebsch-Gordan coefficients. E.g., CG[fStruct@ group, {A, B, C}] for \!\(\*SuperscriptBox[\(f\), \(ABC\)]\)";
eps::usage     = "Standard head for Levi-Civita type Clebsch-Gordan coefficients for SU(N) groups. E.g., CG[eps@ group, {a, b, ...}] for \!\(\*SuperscriptBox[\(\[Epsilon]\), \(ab ... \)]\), with indices in the fundamental representation.";
dSym::usage    = "Standard head for the totally symmetric d Clebsch-Gordan coefficients for SU(N\[GreaterEqual]3) groups. E.g., CG[dSym@ group, {A, B, C}] for \!\(\*SuperscriptBox[\(d\), \(ABC\)]\).";


FieldStrength::usage = "FieldStrength[V, {Index[\[Mu],Lorentz], Index[\[Nu],Lorentz]}, {indices}, {CDerivs}] denotes the field-strength tensor \!\(\*SubscriptBox[\(V\), \(\[Mu]\[Nu]\)]\) of the vector field V."
FS::usage            = "FS[V,\[Mu],\[Nu],indices] denotes the field-strength tensor of an already defined vector field V with Lorentz indices \[Mu] and \[Nu], and gauge indices as the last argument."


GetGaugeGroups::usage   = "GetGaugeGroups[] returns an association of all gauge groups that have been defined by DefineGaugeGroup and their properties. GetGaugeGroups[GaugeGroup] returns an association with all the properties of an already defined gauge group.";
FundAlphabet::usage     = "FundAlphabet is an option for the routine DefineGaugeGroup that specifies a printing alphabet for the indices of the fundamental representation.";
AdjAlphabet::usage      = "AdjAlphabet is an option for the routine DefineGaugeGroup that specifies a printing alphabet for the indices of the adjoint representation.";
DefineGaugeGroup::usage = "DefineGaugeGroup[groupName, LieGroup, coupling, fieldName] defines a gauge group labelled groupName from a given simple LieGroup, together with its corresponding coupling and gauge bosons labelled fieldName. DefineGaugeGroup also defines the fundamental and adjoint representations of the corresponding group, with some basic Clebsch-Gordan coefficients. This routine has two optional arguments, FundAlphabet and AdjAlphabet, to define respectively the printing alphabet of the fundamental and adjoint representations.";
RemoveGaugeGroup::usage = "RemoveGaugeGroup[alias] undefines the gauge group with given alias, together with its associated representations, Clebsch-Gordan coefficients, gauge coupling and bosons.";
ResetGaugeGroups::usage = "Undefines all gauge groups and associated bosons, couplings, representations, and Clebsch-Gordan coefficients."


GetGlobalGroups::usage   = "GetGlobalGroups[] returns an association of all global groups that have been defined by DefineGlobalGroup and their properties. GetGlobalGroups[GlobalGroup] returns an association with all the properties of an already defined global group.";
DefineGlobalGroup::usage = "DefineGlobalGroup[groupName, LieGroup] defines a global group labelled groupName from a given simple LieGroup. DefineGlobalGroup also defines the fundamental and adjoint representations of the corresponding group, together with some basic Clebsch-Gordan coefficients. This routine has two optional arguments, FundAlphabet and AdjAlphabet, to define respectively the printing alphabet of the fundamental and adjoint representations.";
RemoveGlobalGroup::usage = "RemoveGlobalGroup[alias] undefines the global group with given alias together with its associated Clebsch-Gordan coefficients.";
ResetGlobalGroups::usage = "Undefines all global groups and associated representations, and Clebsch-Gordan coefficients."


ResetAll::usage = "Undefines all fields, couplings, flavor indices, global and gauge groups."


PlusHc::usage = "PlusHc[expr] returns expr + Bar[expr] and can be used for convenience when writing the Lagrangian."


FreeLag::usage         = "FreeLag[f] returns the free Lagrangian for the fields specified by the argument f, which can be either a single label or a sequence of labels. If the function is called without arguments, it returns the free Lagrangian of all fields which have been defined.";
LoadModel::usage       = "LoadModel[ModelName] loads the model 'ModelName' previously defined in the Models folder in Matchete. The optional arguments lets the user change the name/value of any of the model parameters or add/modify options of fields, couplings, or groups.";
GetModels::usage       = "Returns an association with all available model files that can be used with the LoadModel[...] command.";
ModelParameters::usage = "ModelParameters is an option of the LoadModel routine."
SetModelOptions::usage = "SetModelOptions is an option of the LoadModel routine."
ParentModel::usage     = "ParentModel[<model file>] is used to indicate that a model file builds directly on another."
ParameterDefault::usage= "ParameterDefault[<parameter> -> <default value>] is used to indicate a default value for a parameter in a model file if no other value is given during loading."


(* ::Subsubsection::Closed:: *)
(*Internal*)


$FieldAssociation::usage    = "$FieldAssociation is an association containing the complete information on all defined fields."
$CouplingAssociation::usage = "$CouplingAssociation is an association containing the complete information on all defined couplings."
$FlavorIndices::usage       = "$FlavorIndices is an association containing the complete information on all defined flavor indices."
$GaugeGroups::usage         = "$GaugeGroups is an association containing information on all defined gauge groups. It does not contain neither representations nor Clebsch-Gordan coefficients associated to the group."
$GlobalGroups::usage        = "$GlobalGroups is an association containing information on all defined global groups. It does not contain neither representations nor Clebsch-Gordan coefficients associated to the group."


KinOpLagrangian::usage      = "KinOpLagrangian[field] returns the  the part of the free Lagrangian that reproduces the kinetic piece of the fluctuation operator (see [2012.08506] eq. (2.3))."


GroupFromRep::usage         = "Returns the name of the group associtated to a given representation."
FieldsFromGroup::usage      = "Returns a list of fields that are charged under a given group."
CouplingsFromGroup::usage   = "Returns a list of couplings that are charged under a given group."
FieldGenerators::usage      = "FieldGenerators[Field,GaugeGroup,Indices] returns the generator of the given GaugeGroup in the same representation as the given Field, with given Indices."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Field definitions*)


(* ::Subsubsection::Closed:: *)
(*Association of all fields*)


$FieldAssociation=<||>;
GetFields[FieldName___]:=Return@$FieldAssociation[FieldName];


GetFieldsByProperty[]:=GetFields[]

GetFieldsByProperty[propsSeq__]:=GetFieldsByProperty[List[propsSeq]]

GetFieldsByProperty[props_Association]:= GetFieldsByProperty@ Normal@ props;

GetFieldsByProperty[props:_List|_Rule]:=
	Keys@ Select[$FieldAssociation, MatchQ[#, KeyValuePattern[props]]&];


(* ::Subsection:: *)
(*Defining new fields*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineField::FieldLabel   = "The label '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another label.";
DefineField::FieldType    = "The field type '`1`' is invalid. Please use either Scalar, Fermion or Vector.";
DefineField::Heavy0       = "The mass of the '`1`' cannot be simultaneously `Heavy` and `0`.";
DefineField::LorentzIndex = "Lorentz cannot be part of the Indices list.";
DefineField::MassIndices  = "The indices in the Mass option must be a subset of the field flavor Indices";
DefineField::MassLabel    = "A mass with label '`1`' has already been defined. Please remove the previous definition or choose another label for the '`2`' field.";
DefineField::Chiral       = "The option Chiral can only be different from False for fields of the Fermion type with SelfConjugate->False";
DefineField::SelfGhost    = "The option SelfConjugate->True is not valid for Ghost fields";
DefineField::ChiralMass   = "A chiral fermion with U(1) charges or complex indices can only be created with a vanishing mass, i.e. by adding the option Mass->0."
DefineField::Complex      = "A real field cannot transform as a complex representation or be charged under an U(1) group."
DefineField::LightIndices = "DefineField does not support light scalar/vector masses with flavor indices. Light mass should be set to 0 in DefineField and the squared mass should be defined separately with DefineCoupling with EFTOrder->2. ";
DefineField::DoubleGaugeCharges = "The field contains multiple charges or indices assigned to the same group, only a single charge and a single index is allowed. For non-Abelian groups, define a new representation to have a field transforming under representations other than adjoint and fundamental.";
DefineField::GBscalar     = "Only a scalar field can be a Goldstone boson."
DefineField::pseudorealmassterm = "The chiral fermion is in a pseudoreal representation: a mass term is possible only if it has flavor indices. If no mass term is desired, please provide the option \"Mass -> 0\".";
DefineField::pseudorealheavymass = "The chiral fermion is in a pseudoreal representation and a flavor-diagonal mass term is impossible. Anti-symmetric heavy masses are not presently supported in Matchete.";


(* ::Subsubsection:: *)
(*Options*)


Options[DefineField]={
	Charges->{},
	Chiral-> False,
	GoldstoneBoson-> False,
	Indices->{},
	Mass-> Heavy,
	NiceForm-> Default,
	SelfConjugate->False,
	UFO$Options -> <||>
};


OptionTest[DefineField, Charges]= If[ListQ@ #, And@@ ChargeQ/@ #, ChargeQ@ #]&;
OptionTest[DefineField, Indices]:= If[ListQ@ #,
		And@@ fieldIndQ/@# && DuplicateFreeQ[GroupFromRep/@ Cases[#, Not@* MemberQ[Keys@ $FlavorIndices]]]
	,
		fieldIndQ@ #
	]&;
OptionTest[DefineField, NiceForm]= MatchQ[_String | Default | {_String | Default , _String | Default}];
OptionTest[DefineField, UFO$Options] := AssociationQ;


OptionTest[DefineField, Mass]= MatchQ[Heavy| Light| 0| {Light, 0}| {Heavy| Light, _}|
	{Heavy|Light, _, _? (SubsetQ[Keys@ GetFlavorIndices[], #]&)}];


fieldIndQ= MemberQ[Keys@ $FlavorIndices, #] || MemberQ[
	Join[Keys@ $GlobalGroups, Keys@ $GaugeGroups], GroupFromRep@ #]&


ChargeQ@ symb_Symbol[_Integer | _Rational | _Symbol]:= MemberQ[GetGaugeGroupByProperty[Group-> U1], symb];
ChargeQ@ _= False;


OptionMessage[Charges, DefineField, val_]:= Message[General::optexpectsval, Charges, DefineField, val, "(list of) valid group charge(s)"];
OptionMessage[Indices, DefineField, val_]:= Message[General::optexpectsval, Indices, DefineField, val, "(list of) already defined group representation(s) (only one per group) and/or flavor indice(s)"];
OptionMessage[Mass, DefineField, val_]:= Message[General::optexpectsval, Mass, DefineField, val, "value Heavy, Light, 0, {Light,0}, {Heavy,MassLabel}, {Light,MassLabel}, {Heavy,MassLabel,{FlavorIndex, \[Ellipsis]}} or {Light,MassLabel,{FlavorIndex, \[Ellipsis]}}, with FlavorIndex among the flavor indices of the field, "];
OptionMessage[UFO$Options, DefineField, val_]:=Message[General::optexpectsval, UFO$Options, DefineField, val, "association of necessary UFO parameters, e.g. <|\"pdgcode\" -> ..., ...|>"];


(* ::Subsubsection:: *)
(*DefineField*)


DefineField[fieldLabel_,type_,opts:OptionsPattern[]]? OptionsCheck:=
DefineField[fieldLabel,type,opts]=Module[
	{
		massInfo = OptionValue@Mass,
		scale,
		flavorInds,
		groupInds,
		isPseudoreal,
		massInds,
		massIndCount,
		massLabel,
		fieldInds,
		fieldCharges,
		tmpAssociation
	},

	fieldInds= If[Head@ # === List, #, {#}]&@ OptionValue@ Indices;
	fieldCharges= If[Head@ # === List, #, {#}]&@ OptionValue@ Charges;
	fieldCharges= fieldCharges/. ( _[0|0.] -> Nothing ); (* remove charges explicitly specified to vanish *)

	(*Check that the symbol 'label' is not yet in use.*)
	If[Defined[fieldLabel],
			Message[DefineField::FieldLabel,fieldLabel];
			Abort[]
	];

	(*Check that 'type' is valid.*)
	If[!MatchQ[type,Scalar|Vector|Fermion|Ghost],
		Message[DefineField::FieldType,type];
		Abort[]
	];

	(*Check Lorentz indices are not part on the internal indices *)
	If[MemberQ[fieldInds,Lorentz],
		Message[DefineField::LorentzIndex];
		Abort[]
	];

	(*Check that only complex fermion fields are chiral *)
	If[OptionValue@Chiral=!=False&&(type=!=Fermion||OptionValue@SelfConjugate===True),
		Message[DefineField::Chiral];
		Abort[]
	];

	(*Check that ghosts are not self-conjugate *)
	(*If[OptionValue@SelfConjugate===True && type===Ghost,
		Message[DefineField::SelfGhost];
		Abort[]
	];*)

	(* Check that charges appear only once per group *)
	If[!DuplicateFreeQ[Head/@(fieldCharges)]||
			!DuplicateFreeQ[DeleteCases[GroupFromRep/@ fieldInds, None]],
		Message[DefineField::DoubleGaugeCharges];
		Abort[];
	];

	(* Check that chiral fermion mass is zero if field is charged under any U(1) or has complex indices *)
	If[(type===Fermion && OptionValue@Chiral=!=False
			&& (Length[fieldCharges]=!=0 || Or@@(MemberQ[Select[Keys@$Representations,$Representations[#][Reality]===0&],#]&/@fieldInds)) )
			&&  !MatchQ[massInfo,0|{Light,0}],
		Message[DefineField::ChiralMass];
		Abort[];
	];

	(* Check that indices are not complex when SelfConjugate -> True*)
	If[(OptionValue@SelfConjugate===True
			&& (Length[fieldCharges]=!=0 || Or@@(MemberQ[Select[Keys@$Representations,$Representations[#][Reality]===0&],#]&/@fieldInds))),
		Message[DefineField::Complex];
		Abort[];
	];

	(* Check only scalar fields are given Goldstone boson option*)
	If[OptionValue@ GoldstoneBoson && type =!= Scalar,
		Message[DefineField::GBscalar];
		Abort[];
	];

	(*Check if mass term should have be symmetric or anti-symmetric in flavor indices depending on
		whether the chiral fermions are in a pseudoreal representation or not*)
	isPseudoreal= If[type === Fermion && OptionValue@ Chiral =!= False,
			flavorInds= Intersection[Keys@ $FlavorIndices, fieldInds];
			groupInds= Complement[fieldInds, flavorInds];
			OddQ[Plus@@ (If[$Representations[#, Reality] === -1, 1, 0]&)/@ groupInds]
		,
			False
		];

	(*Interpret mass option*)
	{scale, massLabel, massInds}= Switch[massInfo
		,Heavy|Light,
			If[isPseudoreal,
				Message[DefineField::pseudorealmassterm];
				Abort[];
			];
			{massInfo, Symbol[If[massInfo === Heavy, "M", "m"] <> ToString@ fieldLabel], {}}
		,0| {Light, 0},
			{Light, 0, {}}
		,{Heavy, 0},
			Message[DefineField::Heavy0, fieldLabel];
			Abort[];
		,{_, Except@ 0},
			Append[massInfo, {}]
		,{_, _, _},
			(*Do not allow light scalar and vector masses with flavor indices*)
			If[First@ massInfo === Light && MemberQ[{Scalar, Vector}, type] && Last@ massInfo =!= {},
				Message[DefineField::LightIndices];
				Abort[];
			];

			(* Check that mass indices are consistent with field indices *)
			If[!MatchQ[fieldInds,{OrderlessPatternSequence[Last@ massInfo/.List->Sequence,___]}],
				Message[DefineField::MassIndices];
				Abort[];
			];

			(* Heavy chiral fermions in a pseudoreal are not presently supported *)
			If[isPseudoreal && First@ massInfo === Heavy,
				Message[DefineField::pseudorealheavymass];
				Abort[];
			];
			massInfo
		];

	(*Check if mass label is already in use*)
	If[massLabel =!= 0 && Defined[massLabel],
		Message[DefineField::MassLabel, massLabel, fieldLabel];
		Abort[]
	];

	(*Define mass coupling*)
	Which[
	scale === Heavy && massLabel =!= 0,
		DefineCoupling[massLabel, EFTOrder-> 0,
			SelfConjugate-> True,
			Indices-> massInds,
			DiagonalCoupling-> ConstantArray[True, Length@ massInds],
			NiceForm-> Last[OptionValue[NiceForm], Default]]
	,scale === Light && massLabel =!= 0 && type === Fermion && OptionValue@ Chiral =!= False,
		massIndCount= Length@ massInds;
		DefineCoupling[massLabel, EFTOrder-> 1,
			SelfConjugate-> False,
			Indices-> Join[massInds, massInds],
			Symmetries-> {If[isPseudoreal, AntisymmetricPermutation, SymmetricPermutation]@@
				Join[massIndCount+ Range@ massIndCount, Range@ massIndCount]},
			NiceForm-> Last[OptionValue[NiceForm], Default]]
	,scale === Light && massLabel =!= 0,
		massIndCount= Length@ massInds;
		DefineCoupling[massLabel, EFTOrder-> 1,
			SelfConjugate-> If[massIndCount === 0,
					True
				,
					Join[massIndCount+ Range@ massIndCount, Range@ massIndCount]
				],
			Indices-> Join[massInds, massInds],
			NiceForm-> Last[OptionValue[NiceForm], Default]]
	];

	scale = (scale === Heavy);

	(*Add field to the list of fields*)
	AppendTo[$FieldAssociation,fieldLabel->
		<|
		Type-> type,
		Indices-> fieldInds,
		Charges-> fieldCharges,
		SelfConjugate-> OptionValue@SelfConjugate,
		Chiral-> OptionValue@Chiral,
		Mass-> massLabel,
		Heavy-> scale,
		GoldstoneBoson-> OptionValue@ GoldstoneBoson,
		UFO$Options -> OptionValue @ UFO$Options
		|>
	];

	(* safe NiceForm printing info *)
	If[OptionValue[NiceForm]=!=Default,
		AppendTo[LabelsNiceForm[Field],fieldLabel->First[OptionValue[NiceForm],OptionValue[NiceForm]]];
	];

	(*In case of vectors, add Lorentz to indices list *)
	If[type===Vector, PrependTo[fieldInds,Lorentz]];

	(* Create the usage message for the new field *)
	If[Length[fieldInds]=== 0,
	fieldLabel::usage=ToString[fieldLabel]<>"[]: Gives a " <> ToString[type] <> " field with label "<>ToString[fieldLabel]<>".";
	,
	fieldLabel::usage=ToString[fieldLabel]<>"[indices]: Gives a " <> ToString[type] <> " field with label "<>ToString[fieldLabel]<>". The argument must be a sequence of indices following the order: "<> StringJoin[ToString/@Riffle[fieldInds,","]] <>".";
	];
	If[Length[fieldCharges]=!= 0,fieldLabel::usage=fieldLabel::usage<>" The charges of the field are "<> ToString[fieldCharges,InputForm] <> "."];
	fieldLabel::IndexNumber="Incorrect number of indices specified. Expected `2` indices, but given `1`.";

	(* Define the new field *)
	fieldLabel[indexlist___]:=Module[
		{
			l=fieldLabel,
			t=type,
			inputInd=List[indexlist],
			i=fieldInds
		},

		(*check the length of indices*)
		If[Length[i]!=Length[inputInd],
			Message[fieldLabel::IndexNumber,Length[inputInd],Length[i]];
			Abort[]
		];

		(*Extract Lorentz index for vectors*)
		If[t===Vector,
			t   = Vector@Index[Quiet[First@inputInd,First::normal],Lorentz];
			inputInd = Drop[inputInd,1];
			i = Drop[i,1];
		];

		If[OptionValue@Chiral===False,
			Field[l, t, Thread@Index[inputInd,i], {}],
			If[OptionValue@Chiral===LeftHanded,PL,PR]**Field[l, t, Thread@Index[inputInd,i], {}]
		]
	];

	(*In case of vectors, define field-strength tensor *)
	If[type===Vector,
		FS[fieldLabel, \[Mu]_, \[Nu]_, indexlist___]:=Module[
			{
				inputInd=List[indexlist],
				i=Drop[fieldInds,1]
			},

			(*check the length of indices*)
			If[Length[i]!=Length[inputInd],
				Message[fieldLabel::IndexNumber,Length@inputInd,Length@i];
				Abort[];
			];

			FieldStrength[fieldLabel, {Index[\[Mu],Lorentz], Index[\[Nu],Lorentz]}, Thread@Index[inputInd, i], {}]
		]
	];

	(*Setup of Conj (used to clasify operator types)*)
	If[OptionValue@ SelfConjugate,
		Conj@ fieldLabel= fieldLabel
	];
]


(* ::Subsection:: *)
(*Remove fields*)


(* ::Subsubsection::Closed:: *)
(*Function*)


RemoveField::rmgaugefield="A gauge fields must be removed by removing the associated gauge group with RemoveGaugeGroup[group] or by resetting all gauge groups with ResetGaugeGroups[]. Continuing without removing the gauge field.";


RemoveField[alias_,OptionsPattern[{"gauge"->False}]]:= Module[{m=$FieldAssociation[alias][Mass]},
						If[!OptionValue["gauge"] && GaugeFieldLabelQ[alias],
							Message[RemoveField::rmgaugefield];
							Return[]
						];
						If[KeyExistsQ[LabelsNiceForm[Field],alias],KeyDropFrom[LabelsNiceForm[Field], alias]];
					    If[KeyExistsQ[$FieldAssociation,alias],
							KeyDropFrom[$FieldAssociation, alias];
							(* delete some auxiliary information *)
							If[Head[$currentFieldAssociation]===Association && KeyExistsQ[$currentFieldAssociation,alias],
								KeyDropFrom[$currentFieldAssociation, alias];
							];
							If[KeyExistsQ[$currentHeavyDims,alias],
								KeyDropFrom[$currentHeavyDims, alias];
							];
							$currentEOMs=DeleteCases[$currentEOMs,KeyValuePattern[Field[alias|{alias,_},___]:>_],All];
							Quiet[Conj@ alias=.];
							If[m=!=0 && !MemberQ[Lookup[Mass]@Values@GetFields[],m],RemoveCoupling[m]];
							(* ClearAllValues[alias] is rather slow here, instead clear specific places *)
							ClearAllValues[alias, {
								"DefineField",
								"Matchete`PackageScope`ExpandVectorFluctuations",
								"Matchete`FunctionalTools`PackagePrivate`VarDraw"
							}];
							(* unset $currentLagrangian if necessary *)
							If[!FreeQ[Matchete`Matching`PackagePrivate`$currentLagrangian,alias,All],
								Quiet[Matchete`Matching`PackagePrivate`$currentLagrangian=.]
							];
							ClearAll@alias;
					    ];
					  ];


ResetFields[]:=(RemoveField/@Complement[Keys@GetFields[],If[GetGaugeGroups[]===<||>,{},Values@Transpose[GetGaugeGroups[]][Field]]];);


GaugeFieldLabelQ[alias_]:=If[MemberQ[$GaugeGroups, KeyValuePattern[Field->alias]],
	True,
	False
]


(* ::Section:: *)
(*Coupling definitions*)


(* ::Subsection:: *)
(*Coupling setup*)


(* ::Subsubsection::Closed:: *)
(*Properties*)


Coupling[0, __]= 0;


(* ::Subsubsection::Closed:: *)
(*Association of all couplings*)


$CouplingAssociation=<||>;
GetCouplings[CouplingName___]:=$CouplingAssociation[CouplingName]//Map[ReplaceAll[<|x___,Symmetries-> val_,y___|>:><|x,Symmetries->If[Length[List@@val]===1,{}, val],y|>]]


(* ::Subsubsection::Closed:: *)
(*Remove an already defined coupling*)


RemoveCoupling::rmgaugecoupling="A gauge coupling must be removed by removing the associated gauge group with RemoveGaugeGroup[group] or by resetting all gauge groups with ResetGaugeGroups[]. Continuing without removing the gauge coupling.";


RemoveCoupling[alias_,OptionsPattern[{"gauge"->False}]]:= (
	If[!OptionValue["gauge"] && GaugeCouplingLabelQ[alias],
		Message[RemoveCoupling::rmgaugecoupling];
		Return[]
	];
	If[KeyExistsQ[LabelsNiceForm[Coupling],alias],KeyDropFrom[LabelsNiceForm[Coupling], alias]];
	If[KeyExistsQ[$CouplingAssociation,alias],
		KeyDropFrom[$CouplingAssociation, alias];
		(* ClearAllValues[alias] is rather slow here, hence specify where to clean up *)
		ClearAllValues[alias, {
			"DefineCoupling",
			"DefineGaugeGroup",
			"EvaluateLoopFunctions"
		}];
		(* unset $currentLagrangian if necessary *)
		If[!FreeQ[Matchete`Matching`PackagePrivate`$currentLagrangian,alias,All],
			Quiet[Matchete`Matching`PackagePrivate`$currentLagrangian=.]
		];
		ClearAll@alias;
		(*regenerates the replacement rule to drop diagonal flavor indices on couplings*)
		UpdateDropDiagonalCouplings[];
	];
)


ResetCouplings[]:=Module[{mList1,mList2,cList},
					mList1=GetFields[#][Mass]&/@(Keys@GetFields[]);
					mList2=If[GetGaugeGroups[]===<||>,{},Values@Transpose[GetGaugeGroups[]][Coupling]];
					cList=Keys@GetCouplings[];
					RemoveCoupling/@Complement[cList,Join[mList1,mList2]];
				  ];


GaugeCouplingLabelQ[alias_]:=If[MemberQ[$GaugeGroups, KeyValuePattern[Coupling->alias]],
	True,
	False
]


(* ::Subsection:: *)
(*Defining new couplings*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineCoupling::CouplingLabel=    "The label '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another label.";
DefineCoupling::LorentzIndex=     "Lorentz cannot be part of the Indices list.";
DefineCoupling::DiagonalCoupling= "The Option DiagonalCoupling must be (list of) boolean(s) of the same length as the number of indices carried by the coupling.";


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[DefineCoupling]= {
		DiagonalCoupling-> Default,
		EFTOrder-> 0,
		Indices-> {},
		NiceForm-> Default,
		SelfConjugate-> False,
		Symmetries-> {}
	};


(* ::Text:: *)
(*Option checks and messages  *)


OptionTest[DefineCoupling, DiagonalCoupling]= MatchQ[{_?BooleanQ...}| _?BooleanQ | Default];
OptionTest[DefineCoupling, EFTOrder]= MatchQ[_Integer? NonNegative];
OptionTest[DefineCoupling, Indices]= flavorIndQ@ # || (ListQ[#]&& And@@ flavorIndQ/@#)&;
OptionTest[DefineCoupling, NiceForm]= MatchQ[_String | Default];
OptionTest[DefineCoupling, Symmetries]= (ListQ[#] && And@@((MatchQ[Head[#1],SymmetricIndices|AntisymmetricIndices|SymmetricPermutation|AntisymmetricPermutation] && VectorQ[List@@#1,Positive])&)/@#) || (MatchQ[#,_SymmetryOverride]) &;


flavorIndQ= MemberQ[Keys@ $FlavorIndices, #]|| MemberQ[Keys@$GlobalGroups, GroupFromRep@ #]&


OptionMessage[DiagonalCoupling, func:DefineCoupling, val_]:= Message[General::optexpectsval, DiagonalCoupling, func, val, "(list of) boolean(s) or Default"];
OptionMessage[EFTOrder, DefineCoupling, val_]:= Message[General::optexpectsval, EFTOrder, DefineCoupling, val, "non-negative integer"];
OptionMessage[Indices, func:DefineCoupling, val_]:= Message[General::optexpectsval, Indices, func, val, "(list of) already defined flavor or global indice(s)"];
OptionMessage[Symmetries, func:DefineCoupling, val_]:= Message[General::optexpectsval, Symmetries, func, val, "a list of terms of the form SymmetricIndices[n1,n2,..], AntisymmetricIndices[n1,n2,..], SymmetricPermutation[n1,n2,..], or AntisymmetricPermutation[n1,n2,..], with n1,n2,... being positive integers indicating the index positions, or an empty list (in case of no symmetries)"];


(* ::Subsubsection::Closed:: *)
(*DefineCoupling*)


(* this creates the coupling association entry both used when defining a new coupling, but also for the duplicate check in IntroduceEffectiveCouplings *)
CouplingAssociationEntry[order_,indInternal_, sc_,symmetries_,dc_]:=Module[{},
	<|
			EFTOrder-> order,
			Indices-> indInternal,
			SelfConjugate-> sc,
			Symmetries-> symmetries,
			DiagonalCoupling-> dc
		|>
]


(* ::Text:: *)
(*Main function *)


DefineCoupling[label_,opts:OptionsPattern[]]? OptionsCheck:=
DefineCoupling[label,opts]= Module[
	{
		order = OptionValue@ EFTOrder,
		indInternal,
		symmetries,
		tmpAssociation,
		diagCoup= OptionValue@ DiagonalCoupling
	},

	indInternal= If[Head@ # === List, #, {#}]&@ OptionValue@ Indices;

	(*Check that the symbol 'label' is not yet in use.*)
	If[Defined[label],
		Message[DefineCoupling::CouplingLabel,label];
		Abort[]
	];

	(*Check that DiagonalCoupling and Indices have the same length*)
	diagCoup= Switch[diagCoup
		,Default, ConstantArray[False, Length@ indInternal]
		,_?BooleanQ, {diagCoup}
		,_, diagCoup
		];
	If[Length[diagCoup] =!= Length[indInternal],
		Message[DefineCoupling::DiagonalCoupling]; Abort[]
	];

	(*Check Lorentz indices are not part on the internal indices *)
	If[MemberQ[indInternal,Lorentz],
		Message[DefineCoupling::LorentzIndex];
		Abort[]
	];

	(*Perform check of coupling input*)
	symmetries= CouplingSymmetries[Length@ indInternal, OptionValue@ Symmetries];

	(*Add coupling to the list of couplings*)
	AppendTo[$CouplingAssociation, label-> CouplingAssociationEntry[order,
		indInternal, OptionValue@ SelfConjugate, symmetries, diagCoup]];

	(* safe NiceForm printing info *)
	If[OptionValue[NiceForm]=!=Default,
		AppendTo[LabelsNiceForm[Coupling],label->OptionValue[NiceForm]];
	];

	(* Create the usage message for the new coupling *)
	If[Length[indInternal]== 0,
	label::usage=ToString[label]<>"[]: Gives a coupling with label "<>ToString[label]<>" and EFTOrder "<>ToString[order]<>".";
	,
	label::usage=ToString[label]<>"[indices]: Gives a coupling with label "<>ToString[label]<>" and EFTOrder "<>ToString[order]<>". The argument must be a sequence of indices following the order: "<> StringJoin[ToString/@Riffle[indInternal,","]] <>".";
	];
	label::IndexNumber="Incorrect number of indices specified. Expected `2` indices, but given `1`.";

	(* Define the new coupling *)
	label[indexlist___]:=Module[
		{
			l=label,
			ind=List[indexlist],
			i=indInternal
		},
		(*check the length of indices*)
		If[Length[i]!=Length[ind],
			Message[label::IndexNumber,Length[ind],Length[i]];
			Abort[]
		];

		Coupling[l, Thread@ Index[ind,i],order]
	];

	(*regenerates the replacement rule to drop diagonal flavor indices on couplings*)
	UpdateDropDiagonalCouplings[];
]


(* ::Subsubsection::Closed:: *)
(*Determine all symmetric permutation of coupling indices *)


(* ::Text:: *)
(*Returns an association with viable coupling permutations with their sign*)


CouplingSymmetries[noInds_, symmetries_]:=Module[{couplingPermutations, evenGens, perm, oddGens, newPerm},
	couplingPermutations=<|Range@ noInds-> 1|>;

	evenGens={}; oddGens={};
	Do[
		Switch[Head@ perm
		,SymmetricIndices,
			evenGens= evenGens~ Join~ {Cycles@ List[List@@ perm], Cycles@ List[List@@ perm[[{1, 2}]]]};
		,AntisymmetricIndices,
			AppendTo[oddGens, Cycles@ List[List@@ perm[[{1, 2}]]]];
			(*oddGens= oddGens~ Join~ {Cycles@ List[List@@ perm], Cycles@ List[List@@ perm[[{1, 2}]]]};*)
			If[OddQ@ Length@ perm,
				AppendTo[evenGens, Cycles@ List[List@@ perm] ];,
				AppendTo[oddGens, Cycles@ List[List@@ perm] ];
			]
		,SymmetricPermutation,
			AppendTo[evenGens, FindPermutation@ perm];
		,AntisymmetricPermutation,
			AppendTo[oddGens, FindPermutation@ perm];
		]
	,{perm, symmetries}];

	(*Compose the symmetries *)
	newPerm= couplingPermutations;
	While[Length@ newPerm > 0,
		newPerm= Association@@ Join[
				KeyValueMap[Function[{key, val}, Permute[key, #]-> val], newPerm]&/@ evenGens,
				KeyValueMap[Function[{key, val}, Permute[key, #]-> -val], newPerm]&/@ oddGens
			];
		KeyDropFrom[newPerm, Keys@ couplingPermutations];
		couplingPermutations= Join[couplingPermutations, newPerm];
	];
	couplingPermutations
]


(* ::Text:: *)
(*Override to directly give the the symmetries in the internal form*)


CouplingSymmetries[ _ , SymmetryOverride[x_]]:=Module[{out},
	Switch[Head@x,
		List, out = Association @@ x,
		Association, out = x,
		Rule, out = Association @ x,
		_, Association[ {1,2,3,4}->1 ]
	];
	Return[out]
]


(* ::Subsubsection::Closed:: *)
(*Rule for dropping diagonal indices on Couplings*)


(* ::Text:: *)
(*Rule to drop diagonal flavor indices from couplings when relabeling indices. Defined here once globally for better performance.*)


UpdateDropDiagonalCouplings[]:= Block[{},
	$DropDiagonalCouplings= Join[{_FlavorSum->1}, Table[
			With[{tmp= Position[GetCouplings[x][DiagonalCoupling],True], lab= x},
				If[tmp==={},
					Nothing,
					Coupling[x,ind_List,ord_]:> Coupling[lab, ReplacePart[ind,tmp:>Unique[]],ord] (* Unique ensures that 1/(m[p]-m[r]) yields 1/(m[]-m[])=ComplexInfinity *)
				]
			]
			,
			{x,Keys@GetCouplings[]}
		]
	];
];


UpdateDropDiagonalCouplings[]


(* ::Section:: *)
(*Flavor indices*)


(* ::Subsection:: *)
(*Flavor indices setup*)


(* ::Subsubsection::Closed:: *)
(*Association of all flavor indices*)


$FlavorIndices= <||>;
GetFlavorIndices[FlavorIndex___]:=Return@$FlavorIndices[FlavorIndex];


(* ::Subsubsection::Closed:: *)
(*Removing  flavor indices*)


RemoveFlavorIndex[alias_]:= Module[{},
					       If[KeyExistsQ[$FlavorIndices,alias],
					          RemoveCoupling/@CouplingsFromGroup[alias];
					          RemoveField/@FieldsFromGroup[alias];
					          KeyDropFrom[$IndexAlphabets, alias];
							  KeyDropFrom[$FlavorIndices, alias];
							  Bar[Index[ind_,alias]]=.;
							  ClearAllValues[alias];
							  ClearAll@alias;
					       ]
					      ]


ResetFlavorIndices[]:=(RemoveFlavorIndex/@Keys[$FlavorIndices[]];);


(* ::Subsection:: *)
(*Defining flavor indices*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineFlavorIndex::FlavorName     = "The flavor index name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another flavor index name.";
DefineFlavorIndex::IndexDimension = "The flavor index dimension '`1`' is not a positive integer larger than 1 or an undefined Symbol."


(* ::Subsubsection::Closed:: *)
(*Options*)


(* Define the option patterns for DefineFlavorIndices *)
Options[DefineFlavorIndex]={IndexAlphabet -> None};


(* ::Subsubsection::Closed:: *)
(*DefineFlavorIndices*)


(* ::Text:: *)
(*Function to define a flavor indices*)


DefineFlavorIndex[flavorName_,indexDim_,opts:OptionsPattern[]]? OptionsCheck:=
DefineFlavorIndex[flavorName,indexDim,opts]= Module[{},
	(*Check if 'Name' is already defined as a flavor index.*)
	If[Defined[flavorName]||KeyExistsQ[$FlavorIndices,flavorName]||KeyExistsQ[$GaugeGroups,flavorName]||KeyExistsQ[$GlobalGroups,flavorName],
		Message[DefineFlavorIndex::FlavorName,flavorName];
		Abort[]
	];

	(*Check if 'IndexDimension' is an undefined symbol or a positive integer larger than 1.*)
	If[Defined[indexDim]&&!(IntegerQ[indexDim] && indexDim > 1),
		Message[DefineFlavorIndex::IndexDimension,indexDim];
		Abort[]
	];

	(*Update list of flavor indices*)
	AppendTo[$FlavorIndices, flavorName-> <|
			IndexDimension-> indexDim
	|>];

	(* If provided, add index alphabet for the representation *)
	If[OptionValue@IndexAlphabet=!=None,AppendTo[$IndexAlphabets,flavorName -> BuildIndexAssoc[OptionValue@IndexAlphabet]]];

	(* Define a usage message for this object *)
	flavorName::usage=ToString[flavorName]<>" is a label for a flavor index.";

	(* Define flavor indices to be neutral under Bar[]*)
	Bar@Index[ind_,flavorName]:=Index[ind,flavorName];
]


(* ::Section:: *)
(*Gauge and global groups*)


(* ::Subsubsection::Closed:: *)
(*Representation and  Clebsch Gordan info*)


(* ::Text:: *)
(*Shows the main representation properties*)


RepresentationProperties[RepName_]:=Module[{real= $Representations[RepName][Reality]},
									<|
										Group -> $Representations[RepName][GroupName],
										DynkinCoefficients ->  $Representations[RepName][DynkinCoefficients],
										Dimension -> $Representations[RepName][RepDimension],
										Reality -> Which[real==1,"Real",real==-1,"PseudoReal",real==0,"Complex"]
									|>

]


(* ::Text:: *)
(*Shows the representation indices for a given Clebsch Gordan*)


ClebschGordanIndices::ClebschGordanName   = "The Clebsch-Gordan coefficient has not been defined. Please use another name or use DefineCG to define it."


ClebschGordanIndices[CGname_]:= Module[{},
								(*Check that the Clebsch-Gordan is defined.*)
								If[!KeyExistsQ[$CGproperties,CGname],
									Message[ClebschGordanIndices::ClebschGordanName,CGname];
									Abort[]
								];
							    $CGproperties[CGname][Indices]
							  ]


(* ::Subsubsection::Closed:: *)
(*GroupFromRep*)


(* ::Text:: *)
(*For finding the group associated to a representation*)


GroupFromRep@ Bar@ rep_:= GroupFromRep@ rep;
GroupFromRep@ rep_:= Module[{gr},
	gr= Lookup[$Representations, rep, None]@ GroupName;
	If[KeyExistsQ[$GaugeGroups, gr]||KeyExistsQ[$GlobalGroups, gr], gr, None]
];


(* ::Text:: *)
(*Returns the group associated to an index*)


GroupFromInd@ Bar@ ind_:= GroupFromInd@ ind;
GroupFromInd@ Index[_, rep_]:= GroupFromRep@ rep;


(* ::Subsubsection::Closed:: *)
(*FieldsFromGroup*)


FieldsFromGroup[group_]:= Module[{tmp},
	If[MemberQ[Keys@$FlavorIndices,group],
		Return[Select[Keys@$FieldAssociation,MemberQ[Transpose[$FieldAssociation][Indices][#],group]&]];
		,
		If[MatchQ[Join[$GaugeGroups,$GlobalGroups][group][Group], U1],
			tmp=Select[Keys@$FieldAssociation,MemberQ[Head/@Transpose[$FieldAssociation][Charges][#],group]&];
			Return[If[MemberQ[Keys@$GaugeGroups,group],AppendTo[tmp,$GaugeGroups[group][Field]],tmp]];
			,
			Return[Select[Keys@$FieldAssociation,MemberQ[GroupFromRep/@Transpose[$FieldAssociation][Indices][#],group]&]];
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*CouplingsFromGroup*)


CouplingsFromGroup[group_]:= Module[{},
	If[MemberQ[Keys@$FlavorIndices,group],
		Return[Select[Keys@$CouplingAssociation,MemberQ[Transpose[$CouplingAssociation][Indices][#],group]&]];
		,
		If[MatchQ[Join[$GaugeGroups,$GlobalGroups][group][Group], U1],
			Return[Select[Keys@$CouplingAssociation,MemberQ[Head/@Transpose[$CouplingAssociation][Charges][#],group]&]];
			,
			Return[Select[Keys@$CouplingAssociation,MemberQ[GroupFromRep/@Transpose[$CouplingAssociation][Indices][#],group]&]];
		]
	]
]


(* ::Subsection:: *)
(*Get field generators*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


FieldGenerators::Field       = "The field '`1`' has not been defined. Please define the field using DefineField or use another field name.";
FieldGenerators::GroupName   = "The group '`1`' has not been defined. Please define the group using DefineGaugeGroup or DefineGlobalGroup, or use another group name.";
FieldGenerators::Abelian     = "The group '`1`' is Abelian so no indices should be provided. Please use FieldGenerators[`2`, `1`] instead or use a different group name."
FieldGenerators::NotAbelian  = "The group '`1`' is not Abelian so generator indices should be provided. Please use FieldGenerators[`2`, `1`,Indices] instead or use a different group name."
FieldGenerators::Charged     = "The field '`1`' is not charged under the group '`2`'. Please provide another group name."
FieldGenerators::Indices     = "The indices argument '`1`' is incorrect. Please provide a list of three Symbols";


(* ::Subsubsection::Closed:: *)
(*FieldGenerators*)


(* ::Text:: *)
(*For finding the generators associated to a given field and gauge group*)


FieldGenerators[Bar@field_,GroupName_,indices_]:=Bar@FieldGenerators[field,GroupName,indices]
FieldGenerators[field_,GroupName_,indices_]:= Module[{rep,generators},

	(*Check that the field is already defined*)
	If[!KeyExistsQ[$FieldAssociation,field],
		Message[FieldGenerators::Field,field];
		Abort[]
	];

	(*Check that the group is already defined*)
	If[!KeyExistsQ[$GaugeGroups,GroupName] && !KeyExistsQ[$GlobalGroups,GroupName],
		Message[FieldGenerators::GroupName,GroupName];
		Abort[]
	];

	(* Check that the group is not Abelian *)
	If[$GaugeGroups[GroupName][Abelian] || $GlobalGroups[GroupName][Abelian],
		Message[FieldGenerators::Abelian,GroupName,field];
		Abort[]
	];

	(*Check that the field is charged under the given group*)
	If[!MemberQ[DeleteCases[GroupFromRep/@$FieldAssociation[field][Indices],None],GroupName],
		Message[FieldGenerators::Charged,field,GroupName];
		Abort[]
	];

	(*Check that indices is a list of three indices*)
	If[!MatchQ[indices,{_Symbol,_Symbol,_Symbol}],
		Message[FieldGenerators::Indices,indices];
		Abort[]
	];

	rep=First@Select[GetFields[field][Indices],GroupFromRep[#]==GroupName &];
	If[Head@rep===Bar,
		Bar@CG[gen[rep[[1]]],indices]
	,
		CG[gen[rep],indices]
	]
];


(* ::Text:: *)
(*In the case of Abelian fields it returns the charge *)


FieldGenerators[Bar@field_,GroupName_]:=Bar@FieldGenerators[field,GroupName]
FieldGenerators[field_,GroupName_]:= Module[{generators},

	(*Check that the field is already defined*)
	If[!KeyExistsQ[$FieldAssociation,field],
		Message[FieldGenerators::Field,field];
		Abort[]
	];

	(*Check that the group is already defined*)
	If[!KeyExistsQ[$GaugeGroups,GroupName] && !KeyExistsQ[$GlobalGroups,GroupName],
		Message[FieldGenerators::GroupName,GroupName];
		Abort[]
	];

	(* Check that the group is Abelian *)
	If[!$GaugeGroups[GroupName][Abelian] || !$GlobalGroups[GroupName][Abelian],
		Message[FieldGenerators::NotAbelian,GroupName,field];
		Abort[]
	];

	(*Check that the field is charged under the given group*)
	If[!MemberQ[Head/@$FieldAssociation[field][Charges],GroupName],
		Message[FieldGenerators::Charged,field,GroupName];
		Abort[]
	];

	Select[$FieldAssociation[field][Charges],Head[#]==GroupName&][[1]][[1]]
];


(* ::Section:: *)
(*Gauge groups*)


(* ::Subsection:: *)
(*Gauge setup*)


(* ::Subsubsection::Closed:: *)
(*Field-strength properties*)


FieldStrength[_, {\[Mu]_, \[Mu]_}, __]:= 0;
FieldStrength[_, _, {a_, a_}, __]:= 0;
FieldStrength[_, _, {a_, Bar@ a_}, __]:= 0;
FieldStrength[_, _, {Bar@ a_, a_}, __]:= 0;


(* Order the Lorentz indices canonically *)
FieldStrength[label_, {\[Mu]_,\[Nu]_}, rest___] := -FieldStrength[label, {\[Nu],\[Mu]}, rest] /; !OrderedQ[{\[Mu],\[Nu]}]


(* ::Subsubsection::Closed:: *)
(*Association of all gauge groups*)


(* ::Text:: *)
(*Basic gauge group information*)


$GaugeGroups= <||>;
GetGaugeGroups[GaugeGroupName___]:=Module[{$tmpGauge},
	$tmpGauge=$GaugeGroups;
	Table[
		AppendTo[$tmpGauge[i], ClebschGordanCoefficients -> Keys@Select[$CGproperties,#[GroupName]==i&]];
		AppendTo[$tmpGauge[i], Representations -> Keys@Select[$Representations,#[GroupName]==i&]];
	,{i,Keys@$GaugeGroups}];
	$tmpGauge[GaugeGroupName]
];


GetGaugeGroupByProperty[]:= Keys@ $GaugeGroups;

GetGaugeGroupByProperty[propsSeq__]:= GetGaugeGroupByProperty@ List@ propsSeq;

GetGaugeGroupByProperty[props_Association]:= GetGaugeGroupByProperty@ Normal@ props;

GetGaugeGroupByProperty[props: _List| _Rule]:=
	Keys@ Select[$GaugeGroups, MatchQ[#, KeyValuePattern[props]]&];
(* GetGaugeGroupByProperty[props: (_List| _Rule)]:= Module[{gg= $GaugeGroups},
	Table[If[MatchQ[gg[k], KeyValuePattern[props]], k, Nothing], {k, Keys[gg]}]
] *)


(* ::Subsubsection:: *)
(*Removing gauge group*)


RemoveGaugeGroup[groupName_]:= 
	If[KeyExistsQ[$GaugeGroups, groupName],
		RemoveCoupling[$GaugeGroups[groupName, Coupling], "gauge"-> True];
		RemoveField[#, "gauge"-> True]&/@ FieldsFromGroup[groupName];
		RemoveCG/@ GetGaugeGroups[groupName][ClebschGordanCoefficients];
		RemoveRepresentation/@ GetGaugeGroups[groupName][Representations];
		KeyDropFrom[$IndexAlphabets, Keys@ Select[$Representations, #[GroupName] === groupName &]];
		If[!MatchQ[$GaugeGroups[groupName, Group], U1], KeyDropFrom[$Groups, groupName];];
		KeyDropFrom[$GaugeGroups, groupName];
		ClearAllValues@ groupName;
		ClearAll@ groupName;
	];


ResetGaugeGroups[]:=(RemoveGaugeGroup/@ Keys[$GaugeGroups[]];);


(* ::Subsection:: *)
(*Defining gauge group *)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineGaugeGroup::CouplingLabel   = "The coupling '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another symbol for the gauge coupling.";
DefineGaugeGroup::GroupName       = "The group name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another group name.";
DefineGaugeGroup::GaugeFieldLabel = "The gauge field name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another symbol for the gauge field.";
DefineGaugeGroup::LieGroup        = "'`1`' is not a known simple group.";


(* ::Subsubsection:: *)
(*DefineGaugeGroup*)


(* ::Text:: *)
(*Function to define a gauge group together with some basic representations and Clebsch-Gordan coefficients*)


Options[DefineGaugeGroup]= {AdjAlphabet-> None, FundAlphabet -> None, NiceForm-> Default, UFO$Options -> <||>}


OptionTest[DefineGaugeGroup, NiceForm]= MatchQ[Default | {_String | Default, _String | Default}];
OptionTest[DefineGaugeGroup, UFO$Options] := AssociationQ;


OptionMessage[NiceForm, DefineGaugeGroup, val_]:= Message[General::optexpectsval, NiceForm, DefineGaugeGroup, val, "value Default or a list with two elements each of which is either Default or a string"];
OptionMessage[UFO$Options, DefineGaugeGroup, val_]:=Message[General::optexpectsval, UFO$Options, DefineGaugeGroup, val, "association of necessary UFO parameters, e.g. <|\"pdgcode\" -> ..., ...|>"];


DefineGaugeGroup[grName_,lieGroup_, coupling_, gaugeField_,opts:OptionsPattern[]]? OptionsCheck:=
DefineGaugeGroup[grName,lieGroup, coupling, gaugeField,opts]= Module[{fundRep, adjRep, rank, nfCoupling, nfField},
	(*Check that the symbol 'coupling' is not yet in use.*)
	If[Defined[coupling] && !KeyExistsQ[$GaugeGroups,coupling],
		Message[DefineGaugeGroup::CouplingLabel,coupling];
		Abort[]
	];

	(*Check that the symbol 'grName' is not yet in use.*)
	If[Defined[grName]||KeyExistsQ[$FlavorIndices,grName]||KeyExistsQ[$GaugeGroups,grName]||KeyExistsQ[$GlobalGroups,grName],
		Message[DefineGaugeGroup::GroupName,grName];
		Abort[]
	];

	(*Check that the symbol 'gaugeField' is not yet in use.*)
	If[Defined[gaugeField],
		Message[DefineGaugeGroup::GaugeFieldLabel,gaugeField];
		Abort[]
	];

	(*Check that the symbol 'lieGroup' is a valid simple Lie group.*)
	If[!MatchQ[lieGroup, U1|Alg["A",n_?Positive]|Alg["B",n_?(#>1&)]|Alg["C",n_?(#>1&)]|Alg["D",n_?(#>2&)]|Alg["E",n_?(MatchQ[#,6|7|8]&)]|Alg["F",4]|Alg["G",2]],
		Message[DefineGaugeGroup::LieGroup,lieGroup];
		Abort[]
	];

	(*Update list of gauge groups*)
	AppendTo[$GaugeGroups, grName-> <|
			Group-> lieGroup,
			Coupling-> coupling,
			Field-> gaugeField,
			Abelian-> lieGroup === U1
	|>];

	(*Define usage message*)
	grName::usage=ToString[grName]<>" is the label of a gauge group with gauge field "<>ToString[gaugeField]<> " and gauge coupling "<>ToString[coupling]<>".";

	(*Set up Lie group for non-Ablelian gauge groups*)
	If[lieGroup=!= U1,
		(*Add group*)
		DefineCGGroup[grName, lieGroup];
		adjRep= AdjointRepresentation@ lieGroup;
		fundRep= FundamentalRepresentation@ lieGroup;
		rank= Last@ lieGroup;

		(*Add common representations*)
		DefineRepresentation[grName@ adj, grName, adjRep, IndexAlphabet-> OptionValue@AdjAlphabet];
		If[lieGroup=!=Alg["E",8], (*For E_8, the fundamental conincide with the Adjoint*)
			DefineRepresentation[grName@ fund, grName, fundRep,
				IndexAlphabet -> OptionValue@ FundAlphabet];
		];

		(*Structure constatns and identification with adj generators*)
		DefineCG[fStruct@grName,{grName@adj,grName@adj,grName@adj},StructureConstants@lieGroup];
		CG[gen[grName[adj]], indices_]:=-I CG[fStruct@ grName, indices];

		(*Define the commonly occuring t^ai_j f^{abc}*)
		DefineCompositeCG[tFundf@ grName, {gen@ grName@ fund, fStruct@ grName}, {{$a, $i, $j}, {$a, $b, $c}}];

		(*Add common CGs for SU(N) algebras (limited to N\[LessEqual]5 for performance reasons)*)
		If[MatchQ[lieGroup, Alg["A", _?(#>0 && #<5 &)]],
			DefineCG[eps@grName, ConstantArray[grName@fund, rank+ 1],
				First@InvariantTensors[lieGroup, Table[fundRep,rank+ 1],
					AntisymmetricIndices-> Range[rank+ 1],
					Normalization-> (rank+ 1)!]];
		];

		If[MatchQ[lieGroup, Alg["A", _?(#>1 &)]],
			DefineCG[dSym@grName,{grName@adj, grName@adj, grName@adj},
				First@ InvariantTensors[lieGroup, {adjRep, adjRep, adjRep},
					SymmetricIndices-> {1, 2, 3},
					Normalization-> ((rank+1)^2-1)((rank+1)^2-4)/(rank+1)]];
		];

		(*Add two-index invariant CGs for Sp(N) algebras*)
		If[MatchQ[lieGroup, Alg["C", _]],
			DefineCG[eps@ grName, Table[grName@ fund, 2],
				First@ InvariantTensors[lieGroup, Table[FundamentalRepresentation@lieGroup, 2],
					AntisymmetricIndices-> {1, 2}, Normalization-> 2* lieGroup[[2]]]];
		];
	];

	nfCoupling=First[OptionValue[NiceForm],Default];
	nfField=Last[OptionValue[NiceForm],Default];

	(*Add gauge coupling*)
	DefineCoupling[coupling, SelfConjugate -> True, NiceForm->nfCoupling];

	(*Add gauge field*)
	DefineField[gaugeField, Vector, Indices-> If[lieGroup=!= U1,{grName@ adj},{}], Mass-> 0, SelfConjugate-> True, NiceForm->nfField, UFO$Options -> OptionValue @ UFO$Options];
];


(* ::Text:: *)
(*Order indices of structure constants*)


CG[fStruct[group_], indices_List]/;(!OrderedQ[indices]):=Signature[indices]CG[fStruct[group],Sort@indices]


(* ::Section:: *)
(*Global groups*)


(* ::Subsection:: *)
(*Global groups setup*)


(* ::Subsubsection::Closed:: *)
(*Association of all global groups*)


(* ::Text:: *)
(*Basic global group information*)


$GlobalGroups= <||>;
GetGlobalGroups[GlobalGroupName___]:=Module[{$tmpGlobal},
	$tmpGlobal=$GlobalGroups;
	Table[
		AppendTo[$tmpGlobal[i], ClebschGordanCoefficients -> Keys@Select[$CGproperties,#[GroupName]==i&]];
		AppendTo[$tmpGlobal[i], Representations -> Keys@Select[$Representations,#[GroupName]==i&]];
	,{i,Keys@$GlobalGroups}];
	$tmpGlobal[GlobalGroupName]
];


(* ::Subsubsection::Closed:: *)
(*Removing global group*)


RemoveGlobalGroup[groupName_]:= 
	If[KeyExistsQ[$GlobalGroups, groupName],
		RemoveCG/@ GetGlobalGroups[groupName][ClebschGordanCoefficients];
		RemoveRepresentation/@ GetGlobalGroups[groupName][Representations];
		KeyDropFrom[$IndexAlphabets, Keys@ Select[$Representations, #[GroupName] === groupName &]];
		If[!MatchQ[$GlobalGroups[groupName, Group], U1], KeyDropFrom[$Groups, groupName];];
		KeyDropFrom[$GlobalGroups, groupName];
		ClearAllValues@ groupName;
		ClearAll@ groupName;
	];


ResetGlobalGroups[]:=(RemoveGlobalGroup/@Keys[$GlobalGroups[]];);


(* ::Subsection:: *)
(*Defining global group *)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineGlobalGroup::GroupName = "The group name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another group name.";
DefineGlobalGroup::LieGroup  = "'`1`' is not a known simple group.";


(* ::Subsubsection::Closed:: *)
(*Options*)


(* Define the option patterns for DefineGlobalGroup *)
Options[DefineGlobalGroup] = {FundAlphabet -> None, AdjAlphabet -> None}


(* ::Subsubsection::Closed:: *)
(*DefineGlobalGroup*)


(* ::Text:: *)
(*Function to define a global group together with some basic representations and Clebsch-Gordan coefficients*)


DefineGlobalGroup[grName_,lieGroup_,opts:OptionsPattern[]]? OptionsCheck:=
DefineGlobalGroup[grName,lieGroup,opts]=Module[{},
	(*Check that the symbol 'grName' is not yet in use.*)
	If[Defined[grName]||KeyExistsQ[$FlavorIndices,grName]||KeyExistsQ[$GaugeGroups,grName]||KeyExistsQ[$GlobalGroups,grName],
		Message[DefineGlobalGroup::GroupName,grName];
		Abort[]
	];

	(*Check that the symbol 'lieGroup' is a valid simple Lie group.*)
	If[!MatchQ[lieGroup, U1|Alg["A",n_?Positive]|Alg["B",n_?(#>1&)]|Alg["C",n_?(#>1&)]|Alg["D",n_?(#>2&)]|Alg["E",n_?(MatchQ[#,6|7|8]&)]|Alg["F",4]|Alg["G",2]],
		Message[DefineGlobalGroup::LieGroup,lieGroup];
		Abort[]
	];

	(*Update list of global groups*)
	AppendTo[$GlobalGroups, grName-> <|
			Group-> lieGroup,
			Abelian-> lieGroup === U1
	|>];

	(*Define a usage message*)
	grName::usage=ToString[grName]<> " is the label of a global symmetry.";

	(*Set up Lie group for non-Ablelian gauge groups*)
	If[lieGroup=!= U1,
		(*Add group*)
		DefineCGGroup[grName, lieGroup];

		(*Add common representations*)
		If[lieGroup===Alg["E",8],
			DefineRepresentation[grName@adj, grName, AdjointRepresentation@lieGroup,IndexAlphabet -> OptionValue@AdjAlphabet];
		,
			(*It's important to add adj before fund*)
			DefineRepresentation[grName@ adj, grName, AdjointRepresentation@lieGroup,IndexAlphabet -> OptionValue@AdjAlphabet];
			DefineRepresentation[grName@ fund, grName, FundamentalRepresentation@lieGroup,IndexAlphabet -> OptionValue@FundAlphabet];
		];

		(*Add structure constants and identification with adj generators*)
		DefineCG[fStruct@grName, {grName@ adj, grName@ adj, grName@ adj}, StructureConstants@ lieGroup];
		CG[gen@ grName@ adj, indices_]:= -I CG[fStruct@ grName, indices];

		(*Add common CGs for SU(N) algebras (limited to N\[LessEqual]5 for performance reasons)*)
		If[MatchQ[lieGroup, Alg["A",n_?(#>0 && #<5 &)]],
			DefineCG[eps@grName,Table[grName@fund,lieGroup[[2]]+1],First@InvariantTensors[lieGroup, Table[FundamentalRepresentation@lieGroup,lieGroup[[2]]+1],AntisymmetricIndices->Range[lieGroup[[2]]+1],Normalization->(lieGroup[[2]]+1)!]];
		];
		If[MatchQ[lieGroup, Alg["A",n_?(#>1 && #<5 &)]],
			DefineCG[dSym@grName,{grName@adj,grName@adj,grName@adj},First@InvariantTensors[lieGroup, {AdjointRepresentation@lieGroup, AdjointRepresentation@lieGroup, AdjointRepresentation@lieGroup},SymmetricIndices->{1,2,3},Normalization->((lieGroup[[2]]+1)^2-1)((lieGroup[[2]]+1)^2-4)/(lieGroup[[2]]+1)]];
		];

		(*Add two-index invariant CGs for Sp(N) algebras*)
		If[MatchQ[lieGroup, Alg["C", _]],
			DefineCG[eps@ grName, Table[grName@ fund, 2],
				First@ InvariantTensors[lieGroup, Table[FundamentalRepresentation@lieGroup, 2],
					AntisymmetricIndices-> {1, 2}, Normalization-> 2* lieGroup[[2]]]];
		];
	];
];


(* ::Section:: *)
(*Reset all*)


(* ::Subsubsection::Closed:: *)
(*Remove all global data*)


ResetAll[]:= Block[{},
		Quiet[Matchete`Matching`PackagePrivate`$currentLagrangian=.];
		RemoveAssociatedDownValues@ LoadModel[_String, ___]; (* Clear LoadModel cache*)
		ResetAuxiliaryMatchingInformation[];
		ResetGaugeGroups[];
		ResetCouplings[];
		ResetFields[];
		ResetFlavorIndices[];
		ResetGlobalGroups[];
		ResetIndexAlphabets[];
		ResetInternalCouplings[];
		ResetOperatorAssociations[];
		ResetSubstitutionInformation[];
		ResetTempCouplings[];
		ResetEvanescentOperators[];

		ClearGroups[];(*needed to reset the CG numbering*)
	];


(* ::Section:: *)
(*Removing values associated with all Matchete` symbols*)


(* ::Subsubsection::Closed:: *)
(*Function for removing all Up-, Down-, and SubValues off all functions that contain a given symbol*)


PackageScope["ClearAllValues"]


ClearAllValues[var_]:=Do[
	If[NameQ[function],
		With[{func=function},
			If[!FreeQ[DownValues[func],var,All],
				DownValues[func]=DeleteCases[DownValues[func],x_/;!FreeQ[First@x,var,All]];
			];
			If[!FreeQ[UpValues[func],var,All],
				UpValues[func]=DeleteCases[UpValues[func],x_/;!FreeQ[First@x,var,All]];
			];
			If[!FreeQ[SubValues[func],var,All],
				SubValues[func]=DeleteCases[SubValues[func],x_/;!FreeQ[First@x,var,All]];
			];
			(* OwnValues cannot contain var on the LHS *)
			(*OwnValues[func]=DeleteCases[OwnValues[func],x_/;!FreeQ[First@x,var,All]];*)
		]
	];
	,
	{function,Names["Matchete`"~~___]}
]


(* ::Subsubsection::Closed:: *)
(*Function clearing all Up-, Down-, and SubValues off a specific functions that contain a given symbol*)


(* ::Text:: *)
(*Note: func must be given as a String including the full Context*)


ClearAllValues[var_, func_]:=If[NameQ[func],
	DownValues[func] = DeleteCases[DownValues[func], x_/;!FreeQ[First@x,var,All]];
	UpValues[func]   = DeleteCases[UpValues[func],   x_/;!FreeQ[First@x,var,All]];
	SubValues[func]  = DeleteCases[SubValues[func],  x_/;!FreeQ[First@x,var,All]];
	(* OwnValues cannot contain var on the LHS *)
	(*OwnValues[func]=DeleteCases[OwnValues[func],x_/;!FreeQ[First@x,var,All]];*)
]


(* ::Text:: *)
(*Note: functions must be given as a List of Strings including the full Context*)


ClearAllValues[var_, functions_List]:=Do[
	If[NameQ[func],
		With[{f=func},
			DownValues[f] = DeleteCases[DownValues[f], x_/;!FreeQ[First@x,var,All]];
			UpValues[f]   = DeleteCases[UpValues[f],   x_/;!FreeQ[First@x,var,All]];
			SubValues[f]  = DeleteCases[SubValues[f],  x_/;!FreeQ[First@x,var,All]];
			(* OwnValues cannot contain var on the LHS *)
			(*OwnValues[f]=DeleteCases[OwnValues[f],x_/;!FreeQ[First@x,var,All]];*)
		]
	]
	,
	{func,functions}
]


(* ::Subsubsection::Closed:: *)
(*Function for checking in which definitions a symbol appears*)


PackageScope["GetAllValues"]

GetAllValues::usage="GetAllValues[var] lists all DownValues, UpValues, SubValues, and OwnValues of any symbol inside the Matchete` context or its subcontexts that contains the symbol var.";

GetAllValues[var_]:=Do[
	If[NameQ[function],
		With[{func=function},
			Cases[DownValues[func], x_/;!FreeQ[x,var,All] :> Echo[x, "DownValue"] ];
			Cases[UpValues[func],   x_/;!FreeQ[x,var,All] :> Echo[x, "UpValue"]   ];
			Cases[SubValues[func],  x_/;!FreeQ[x,var,All] :> Echo[x, "SubValue"]  ];
			Cases[OwnValues[func],  x_/;!FreeQ[x,var,All] :> Echo[x, "OwnValue"]  ];
		]
	]
	,
	{function,Names["Matchete`"~~___]}
]


(* ::Section:: *)
(*Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*Hermitian conjugate*)


(* ::Text:: *)
(*Ease of life function for adding '+H.c.' terms to the Lagrangian*)


PlusHc[x_]:=x+Bar@x


(* ::Subsubsection::Closed:: *)
(*Free Lagrangian*)


FreeLag::UndefinedField = "The field '`1`' has not been defined.";


FreeLag[field_Symbol]:= Module[
	{fieldType, isChiralFermion, massInd, normalization, fieldInd, indK, ind1, ind2, indM, m, i, j, a, \[Mu], \[Nu]},

	(* Check if the field has been defined *)
	If[!KeyExistsQ[$FieldAssociation, field], Message[FreeLag::UndefinedField, field]; Abort[]];

	isChiralFermion= GetFields[field, Chiral] =!= False;

	massInd= GetCouplings[$FieldAssociation[field, Mass]][Indices]// DeleteDuplicates;
	fieldInd= $FieldAssociation[field, Indices];
	fieldType= $FieldAssociation[field, Type];

	indK = ConstantArray[a, Length@ fieldInd]/. List-> Sequence;

	If[GetFields[field, Mass] =!= 0,
		ind1= Sequence@@ ConstantArray[i, Length@ fieldInd];
		ind2= Sequence@@ ConstantArray[j, Length@ fieldInd];
		m= $FieldAssociation[field, Mass]@@ Join[ConstantArray[i, Length@ massInd],
			If[GetFields[field, Heavy], {}, ConstantArray[j, Length@ massInd]] ];
		(*Currently light boson masses with indices are not supported, so we may square the mass*)
		If[GetFields[field, Type] =!= Fermion, m= m^2];
		(*Insert the group structure to make the mass term invariant*)
		m*= Times@@ Table[
				If[!isChiralFermion || MemberQ[Keys@ $FlavorIndices, indType] ||
					$Representations[indType/. Bar-> Identity, Reality] =!= -1
				,
					Delta[Index[i, indType], Index[j, indType]]
				,
					If[Head@ indType === Bar, Identity, Bar]@ CG[eps@ GroupFromRep@ indType, {i, j}]
				]
			, {indType, Complement[fieldInd, If[GetFields[field, Heavy], {}, massInd]]}]
	,
		ind1= indK;
		ind2= indK;
		m= 0;
	];

	If[$FieldAssociation[field, SelfConjugate] === True, 1/2, 1]* Switch[fieldType,
		Scalar|Ghost,
			(Bar[CD[\[Mu],field[indK]]]**CD[\[Mu],field[indK]] - m Bar[field[ind1]]**field[ind2]),
		Fermion,
			If[$FieldAssociation[field][Chiral] === False,
				(I Bar[field[indK]]** \[Gamma][\[Mu]]** CD[\[Mu], field[indK]] - m Bar[field[ind1]]** field[ind2])
			,
				(I Bar[field[indK]]** \[Gamma][\[Mu]]** CD[\[Mu], field[indK]] - (m/2 Bar[CConj[field[ind1]]]** field[ind2]+ Bar@ m/2 Bar[field[ind1]]** CConj[field[ind2]]))
			],
		Vector,
			(*Gauge fields are normalized with their couplings*)
			normalization= First[Query[Select[#@ Field === field&], Key@ Coupling]@ $GaugeGroups, 1&][]^-2;
			(-1/2 normalization* Bar[FS[field,\[Mu],\[Nu],indK]]**FS[field,\[Mu],\[Nu],indK] + m Bar[field[\[Mu],ind1]]**field[\[Mu],ind2])
	]//Contract//RelabelIndices
]


FreeLag[field_, fields__]:=Plus@@FreeLag/@List[field, fields]
FreeLag[]:=Sum[FreeLag[f],{f,Keys[GetFields[]]}]


(* ::Subsubsection::Closed:: *)
(*Kinetic Lagrangian*)


(* ::Text:: *)
(*This is the part of the free Lagrangian that reproduces the Kinetic piece, \[CapitalDelta], of the fluctuation operator (see [2012.08506] eq. (2.3)).*)


KinOpLagrangian[field_, fields__]:= Plus@@ KinOpLagrangian/@ {field, fields};


KinOpLagrangian[]= 0;
KinOpLagrangian[field_Symbol]:=Module[
		{FieldType, MassInd, FieldInd, indK, ind1, ind2, indM, m, i, j, a, \[Mu], \[Nu],fluctuation},
	FieldInd= GetFields[field, Indices];
	FieldType= GetFields[field, Type];

	indK= Sequence@@ Table[a, {n, Length@ FieldInd}];

	If[GetFieldsUpdated[field, Heavy],
		MassInd= DeleteDuplicates@ GetCouplings[GetFieldsUpdated[field, Mass], Indices];
		ind1 = Sequence@@ Table[If[MemberQ[MassInd, FieldInd[[n]]], i, a], {n, Length@ FieldInd}];
		indM = Sequence@@ Table[i, Length@ MassInd];
		ind2 = ind1;
		m= GetFieldsUpdated[field, Mass]@ indM;
	,
		ind1= ind2= indK;
		m= 0;
	];

	If[GetFields[field, SelfConjugate], 1/2, 1]Switch[FieldType
		,Scalar| Ghost,
			Bar[BackgroundCD[\[Mu],field[indK]]]**BackgroundCD[\[Mu],field[indK]] - m^2 Bar[field[ind1]]**field[ind2]
		,Fermion,
			If[GetFields[field][Chiral]===False,
				I Bar[field[indK]]**\[Gamma][\[Mu]]**BackgroundCD[\[Mu],field[indK]] - m Bar[field[ind1]]**field[ind2]
			,
				I Bar[field[indK]]**\[Gamma][\[Mu]]**BackgroundCD[\[Mu],field[indK]] - m/2 (Bar[CConj[field[ind1]]]**field[ind2]+ Bar[field[ind1]]**CConj[field[ind2]])
			]
		,Vector,
			If[GaugeFieldQ@ field,
				-GetFieldsUpdated[field, Coupling]* BackgroundCD[\[Mu], field[\[Nu], indK]]^2
			,
				-Bar@ BackgroundCD[\[Mu], field[\[Nu], indK]]** BackgroundCD[\[Mu], field[\[Nu], indK]] + m^2 Bar[field[\[Mu],ind1]]**field[\[Mu],ind2]
			]
	]//Expand//RelabelIndices
];


(* ::Subsection:: *)
(*Routine to load Lagrangians*)


(* ::Subsubsection::Closed:: *)
(*Options*)


OptionTest[LoadModel, SetModelOptions]=
	MatchQ[{{_Symbol, _Symbol, _Symbol-> _}...}| {_Symbol, _Symbol, _Symbol-> _}];


OptionMessage[SetModelOptions, LoadModel, val_]:=
	Message[General::optexpectsval, SetModelOptions, LoadModel, val, "list of sets of the kind {<DefinitionFunction>, <Model symbol>, <Option> -> <New value>}"];


(* ::Subsubsection::Closed:: *)
(*Load model*)


LoadModel::str= "The model name must be a string."


Options[LoadModel]={ModelParameters -> {}, SetModelOptions-> {}, Verbose-> False}


LoadModel[fileName_,opts:OptionsPattern[]]? OptionsCheck:= LoadModel[fileName, opts]= Module[{modelPath},

	(* Check that the file name is a string and that the model file exits *)
	If[Head@ fileName =!= String, Message[LoadModel::str]; Abort[]; ];

	modelPath= DetermineModelPath[fileName, OptionValue@ Verbose];

	ExecuteModelDefinition[modelPath, Sequence@@ FilterRules[{opts}, Options@ ExecuteModelDefinition]]
]


(* ::Subsubsection::Closed:: *)
(*Determine model path*)


(* ::Text:: *)
(*Determine which model file is indicated by a string and return full path *)


LoadModel::fileext= "The extension \"`1`\" of model name \"`2`\" is not valid. Use .m or .wl (or leave balnk).";
LoadModel::nofiles= "No files matches the model name \"`1`\".";
LoadModel::multfiles= "Multiple files `1` matches the model name \"`2`\". Please specify an extension and/or a valid model path.";


DetermineModelPath[fileName_String, verbose_]:= Block[
		{dir, extension, matches, models, modName, modPossibilities},
	{dir, modName}= Through@ {DirectoryName, FileNameTake}@ fileName;

	modName= StringReplace[modName, RegularExpression["\\.$"]-> ""];
	modPossibilities= Switch[extension= FileExtension@ modName
		,"",
			{modName <> ".m", modName <> ".wl"}
		,"m"|"wl",
			{modName}
		,_,
			Message[LoadModel::fileext, extension, fileName];
			Abort[];
	];

	(*If directory is specified*)
	If[dir =!= "",
		matches= FileNames[modPossibilities, dir];
		Switch[Length@ matches
		,0,
			Message[LoadModel::nofiles, fileName];
			Abort[];
		,1,
			If[verbose, Print["Loading model ", First@ matches]];
			Return@ First@ matches;
		,_,
			Message[LoadModel::multfiles, matches, fileName];
			Abort[];
		];
	];

	(*Full file match*)
	models= Normal@ GetModels[Path-> True];
	matches= Position[models, Alternatives@@ modPossibilities];
	If[Length@ matches === 0,
		(*Partial match*)
		matches= Position[models, mod_String/; StringContainsQ[mod, fileName]];
	];

	matches=  FileNameJoin/@ Transpose@
		{Extract[models, {First@ #, 1}&/@ matches], Extract[models, matches]};
	Switch[Length@ matches
	,0,
		Message[LoadModel::nofiles, fileName];
		Abort[];
	,1,
		If[verbose, Print["Loading model ", First@ matches]];
		First@ matches
	,_,
		Message[LoadModel::multfiles, matches, fileName];
		Abort[];
	]

]


(* ::Text:: *)
(*Get models*)


GetDirectories[]:= Block[{notebookDir},
	(*Returns $Failed if the notebook is unsaved*)
	notebookDir= Quiet[NotebookDirectory[], NotebookDirectory::nosv];
	DeleteDuplicatesBy[Association@@ {
		FileNameJoin@ {"Matchete", "Models"}-> FileNameJoin@{$MatchetePath, "Models"},
		If[notebookDir =!= $Failed, "NotebookDirectory"-> NotebookDirectory[], Nothing],
		"Directory"-> Directory[]
	}, FileNameDrop[#, 0]&]
];


Options@ GetModels= {Path-> False};


GetModels[OptionsPattern[]]? OptionsCheck:= Block[{models, dirs},
	dirs= GetDirectories[];
	models= DeleteCases[(FileNameTake/@ FileNames[{"*.m", "*.wl"}, #]&)/@ dirs, {}];
	If[OptionValue@ Path,
		models= KeyMap[#/. dirs&, models];
	];
	models
];


Quiet[NotebookDirectory[],{NotebookDirectory::nosv}]


(* ::Subsubsection:: *)
(*Read and execute model file*)


(* ::Text:: *)
(*Read model file, change parameters and options accordingly, and execute definition.*)


LoadModel::defsymb= "The symbol `1` is already in use and cannot be used as a new model parameter.";
LoadModel::defaultpar= "The LHS of the default value `1` from the (parent) model file must be a Symbol.";
LoadModel::modelpar= "The parameter `1` does not appear in the model file and cannot be renamed.";
LoadModel::localvar= "The redefined symbol(s) \"`1`\" appear as local variables in Module or Block in the (parent) model. This redefinition is likely to cause an error and is not allowed.";
LoadModel::nulllag= "The Lagrangian returned by the model file `1` is Null (Parent file level `2`). Ensure that the Lagrangian expression is the last statement in the model file and that it DOES NOT terminate with a semicolon ' ; '.";
LoadModel::optscall= "The error happens with the Option call SpecifyOptions\[Rule] {\[Ellipsis], `1` ,\[Ellipsis]}.";
LoadModel::parcall= "The error happens with the Option call ModelParameters\[Rule] {\[Ellipsis], `1` ,\[Ellipsis]}.";
LoadModel::unkwnfunc= "The model file does not define `1`[`2`, \[Ellipsis]], which you are attempting to change options for. Perhaps you changed the model parameter `2`.";
LoadModel::unkwnopt= "Attempt to modify the option `1` for a definition with `2`. This option does not exist. ";
LoadModel::unkwnoptfunc= "Cannot change options associated with function `1`. SetModelOptions should only be used to change options associated with functions `2`.";


Options@ ExecuteModelDefinition= {
	ModelParameters-> {},
	SetModelOptions-> {}
};


ExecuteModelDefinition[file_String, OptionsPattern[]]:= Module[
		{change, changedPars, childLag, defaultParams, lag, modelDef, modelOpts, modelLocalVars, n, overlap,
		parameterSubs, symb},

	(*Load held Model definition*)
	modelDef= ReadList[file, Hold[Expression]];
	modelDef= DeleteCases[modelDef, Hold@ Null];

	modelDef= ReadParentModels@ modelDef;

	(*Determine which are the local variables in the model file(s)*)
	modelLocalVars= Join@@ Cases[modelDef, HoldPattern@ (Module|Block)[pars_List, _]:> pars, All]// DeleteDuplicates;

	(*Check parameter changes*)
	parameterSubs= OptionValue[ModelParameters];
	Do[
		If[FreeQ[modelDef, symb[[1]]],
			Message[LoadModel::modelpar, symb[[1]]];
			Message[LoadModel::parcall, symb];
			Abort[];
		];
		(* the checks below are problematic in light of parameter name overloading done consistently between UV and EFT *)
		(*
		If[MatchQ[symb[[2]], _Symbol] && Defined@ symb[[2]],
			Message[LoadModel::defsymb, symb[[2]]];
			Message[LoadModel::parcall, symb];
			Abort[];
		];
		*)
	
		overlap= Intersection[modelLocalVars, List@@ symb];
		If[Length@ overlap =!= 0,
			Message[LoadModel::localvar, StringJoin@@ Riffle[ToString/@ overlap, ", "] ];
			Message[LoadModel::parcall, symb];
		];
	, {symb, parameterSubs}];

	(*Remove Niceforms where parameter changes are performed*)
	changedPars= parameterSubs[[;;, 1]];
	modelDef= modelDef/. {
			HoldPattern@ DefineCoupling[p_, l___, Rule[NiceForm, _], r___]/; MemberQ[changedPars, p]:>
				DefineCoupling[p, l, r],
			HoldPattern@ DefineField[p_, l___, Rule[NiceForm, _], r___]/; MemberQ[changedPars, p]:>
				DefineField[p, l, r],
			HoldPattern@ DefineGaugeGroup[s:PatternSequence[_, _, coup_, field_], l___, Rule[NiceForm, _], r___]/; SubsetQ[changedPars, {coup, field}]:>
				DefineGaugeGroup[s, l, r],
			HoldPattern@ DefineGaugeGroup[s:PatternSequence[_, _, _, field_], l___, Rule[NiceForm, {nfCoup_, _}], r___]/; MemberQ[changedPars, field]:>
				DefineGaugeGroup[s, l, r, NiceForm-> {nfCoup, Default}],
			HoldPattern@ DefineGaugeGroup[s:PatternSequence[_, _, coup_, _], l___, Rule[NiceForm, {_, nfField_}], r___]/; MemberQ[changedPars, coup]:>
				DefineGaugeGroup[s, l, r, NiceForm-> {Default, nfField}]
		};

	(*Incorporate defualt paramters*)
	defaultParams= Cases[modelDef, ParameterDefault[par_-> val_]:>
		If[MatchQ[par, _Symbol],
			par-> val
		,
			Message[LoadModel::defaultpar, FullForm@ par-> val];
			Abort[];
		]
		, All];
	modelDef= DeleteCases[modelDef, _ParameterDefault, All];

	(*Use value from ModelParameters or the value from the most junior model that provides a value*)
	parameterSubs= DeleteDuplicatesBy[Join[parameterSubs, Reverse@ defaultParams], First];

	(*Change parameter names*)
	modelDef= modelDef/. parameterSubs;

	(*Set/change specific options*)
	modelOpts= OptionValue@ SetModelOptions;
	If[MatchQ[modelOpts, {_Symbol, __}], modelOpts= {modelOpts}];

	Do[
		With[{func= change[[1]], var= change[[2]], opt= change[[3, 1]], newVal= change[[3, 2]]},
			If[!MemberQ[$setModelOptionsFuncs, func],
				Message[LoadModel::unkwnoptfunc, func, $setModelOptionsFuncs];
				Message[LoadModel::optscall, change];
				Abort[];
			];

			If[FreeQ[modelDef, HoldPattern@ func[var, ___]],
				Message[LoadModel::unkwnfunc, func, var];
				Message[LoadModel::optscall, change];
				Abort[];
			];

			If[!MemberQ[Options[func][[;;, 1]], opt],
				Message[LoadModel::unkwnopt, opt, func];
				Message[LoadModel::optscall, change];
				Abort[];
			];

			modelDef= If[FreeQ[modelDef, HoldPattern@ func[var, ___, Rule[opt, _], ___]],
					modelDef/. HoldPattern@ func[var, rest___]:> func[var, rest, opt-> newVal]
				,
					modelDef/. HoldPattern@ func[var, rest1___, Rule[opt, _], rest2___]:>
						func[var, rest1, opt-> newVal, rest2]
				];
		]
	, {change, modelOpts}];

	(*Run model definition over consecutive child models*)
	lag= 0;
	Do[
		childLag= Last@ ReleaseHold@ modelDef[[n]];
		If[childLag === Null,
		Message[LoadModel::nulllag, FileBaseName@ file, Length@ modelDef- n];
			Abort[];
		];
		lag+= childLag;
	, {n, Length@ modelDef}];

	lag
]


(* ::Text:: *)
(*Valid functions to change with SetModelOptions*)


$setModelOptionsFuncs= {DefineCoupling, DefineField, DefineFlavorIndex, DefineGaugeGroup, DefineGlobalGroup};


(* ::Subsubsection::Closed:: *)
(*Determine parent model*)


(* ::Text:: *)
(*Recursively inserts the model definition from ParentModel files into the Model definition *)


LoadModel::findparent= "The parentmodel '`1`' could not be found.";
LoadModel::parentarg= "The argument/option `1` in `2` is invalid.";
LoadModel::parentopt= "The option \"Use Lagrangian\" in `1` takes a boolean value.";
LoadModel::parfile= "The model in ParentModel[`1`] must be a String.";
LoadModel::parmoduse= "ParentModel[<model>] must appear as the first expression in the model file.";


Options@ ParentModel= {
	"Use Lagrangian"-> True
};


ReadParentModels@ modelDef_List:= Module[
		{childDef, parentDef, parentFile, parentOps, parentPath, useParentLag= True},
	If[FreeQ[modelDef, _ParentModel],
		Return@ {modelDef};
	];

	parentFile= First@ modelDef;
	childDef= Rest@ modelDef;

	{parentFile, parentOps}= Switch[parentFile
	, Hold[ParentModel[_, ___] ],
		{parentFile[[1, 1]], List@@ parentFile[[1, 2;;]]}
	, Hold[CompoundExpression[ParentModel[_, ___], Null]], (* If file reads ParentModel[...]; *)
		{parentFile[[1, 1, 1]], List@@ parentFile[[1, 1, 2;;]]}
	, _,
		Message[LoadModel::parmoduse];
		Abort[];
	];

	If[!MatchQ[parentFile, _String],
		Message[LoadModel::parfile, parentFile];
		Abort[];
	];

	(*Check parent options*)
	Do[
		Switch[opt
		, "Use Lagrangian"->  _?BooleanQ,
			useParentLag= opt[[2]];
		, "Use Lagrangian"->  _,
			Message[LoadModel::parentopt, ParentModel[parentFile, Sequence@@ parentOps]];
			Abort[];
		, _,
			Message[LoadModel::parentarg, opt, ParentModel[parentFile, Sequence@@ parentOps]];
			Abort[];
		];
	, {opt, parentOps}];

	(*Read parent model definitions*)
	parentPath= CheckAbort[DetermineModelPath[parentFile, False],
			Message[LoadModel::findparent, parentFile];
			Abort[];
		];

	parentDef= ReadList[parentPath, Hold[Expression]];
	parentDef= DeleteCases[parentDef, Hold@ Null];

	(*Return 0 if Lagrangian should not be used. This way maintains all definitions that may
		be performed in the last file expression, e.g., in a module.*)
	If[!useParentLag,
		AppendTo[parentDef, 0];
	];

	(*Recursively check for parents*)
	parentDef= ReadParentModels@ parentDef;

	Append[parentDef, childDef]
]
