(* Created with the Wolfram Language : www.wolfram.com *)
<|"Model" -> "VLF_toy_model", "Version" -> "0.3.0", 
 "Date" -> DateObject[{2025, 5, 30}, "Day"], "Time (Match)" -> 3.498521, 
 "Time (GreensSimplify)" -> 0.040421, "Time (EOMSimplify)" -> 0.515094, 
 "Time (MapEffectiveCouplings)" -> None, "UV Lagrangian" -> 
  -1/2*(Coupling[m, {}, 1]^2*Field[\[Phi], Scalar, {}, {}]^2) + 
   Field[\[Phi], Scalar, {}, {Index[d$$1, Lorentz]}]^2/2 - 
   FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2/
    (4*Coupling[e, {}, 0]^2) - Coupling[M, {}, 0]*
    Bar[Field[\[CapitalPsi], Fermion, {}, {}]]**Field[\[CapitalPsi], Fermion, 
      {}, {}] + I*Bar[Field[\[Psi], Fermion, {}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]]]**Field[\[Psi], Fermion, {}, 
      {Index[d$$1, Lorentz]}] - Coupling[y, {}, 0]*
    Field[\[Phi], Scalar, {}, {}]*Bar[Field[\[Psi], Fermion, {}, {}]]**
     DiracProduct[Proj[1]]**Field[\[CapitalPsi], Fermion, {}, {}] + 
   I*Bar[Field[\[CapitalPsi], Fermion, {}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]]]**Field[\[CapitalPsi], 
      Fermion, {}, {Index[d$$1, Lorentz]}] - Bar[Coupling[y, {}, 0]]*
    Field[\[Phi], Scalar, {}, {}]*Bar[Field[\[CapitalPsi], Fermion, {}, {}]]**
     DiracProduct[Proj[-1]]**Field[\[Psi], Fermion, {}, {}], 
 "Off-shell EFT Lagrangian" -> 
  (hbar*Bar[Coupling[y, {}, 0]]^3*Coupling[y, {}, 0]^3*
     Field[\[Phi], Scalar, {}, {}]^6)/(3*Coupling[M, {}, 0]^2) + 
   (13*hbar*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
     Field[\[Phi], Scalar, {}, {}]^3*Field[\[Phi], Scalar, {}, 
      {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}])/
    (18*Coupling[M, {}, 0]^2) + (hbar*Bar[Coupling[y, {}, 0]]*
     Coupling[y, {}, 0]*Field[\[Phi], Scalar, {}, {Index[d$$1, Lorentz], 
       Index[d$$1, Lorentz]}]*Field[\[Phi], Scalar, {}, 
      {Index[d$$2, Lorentz], Index[d$$2, Lorentz]}])/
    (3*Coupling[M, {}, 0]^2) + (hbar*Bar[Coupling[y, {}, 0]]*
     Coupling[y, {}, 0]*Field[\[Phi], Scalar, {}, {}]^2*
     FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
      2)/(3*Coupling[M, {}, 0]^2) - 
   (2*hbar*FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
      {Index[d$$2, Lorentz]}]*FieldStrength[A, {Index[d$$1, Lorentz], 
       Index[d$$3, Lorentz]}, {}, {Index[d$$3, Lorentz]}])/
    (15*Coupling[M, {}, 0]^2) + 
   FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2*
    (-1/4*1/Coupling[e, {}, 0]^2 + hbar*(-1/3*1/\[Epsilon] - 
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2]/3)) + 
   hbar*Field[\[Phi], Scalar, {}, {}]^4*
    (-((Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2)/\[Epsilon]) - 
     Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
      Log[\[Mu]bar2/Coupling[M, {}, 0]^2]) - 
   Field[\[Phi], Scalar, {}, {Index[d$$1, Lorentz]}]^2*
    (-1/2 + hbar*(-1/2*(Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]) - 
       (Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0])/\[Epsilon] - 
       Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2])) + 
   Field[\[Phi], Scalar, {}, {}]^2*(-1/2*Coupling[m, {}, 1]^2 + 
     hbar*((-2*Bar[Coupling[y, {}, 0]]*Coupling[M, {}, 0]^2*
         Coupling[y, {}, 0])/\[Epsilon] + Coupling[M, {}, 0]^2*
        (-2*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0] - 
         2*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2]))) + 
   I*Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]]]**Field[\[Psi], Fermion, {}, 
      {Index[d$$1, Lorentz]}] + (7*hbar*Bar[Coupling[y, {}, 0]]*
     Coupling[y, {}, 0]*FieldStrength[A, {Index[d$$1, Lorentz], 
       Index[d$$2, Lorentz]}, {}, {Index[d$$2, Lorentz]}]*
     Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[\[Psi], Fermion, {}, 
       {}])/(36*Coupling[M, {}, 0]^2) + 
   hbar*(((3*I)/4)*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0] + 
     ((I/2)*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0] + 
       (I*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*Coupling[y, {}, 0])/
        Coupling[M, {}, 0]^2)/\[Epsilon] + (I/2)*Bar[Coupling[y, {}, 0]]*
      Coupling[y, {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2] + 
     (((3*I)/2)*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*
        Coupling[y, {}, 0] + I*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*
        Coupling[y, {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^2)*Bar[Field[\[Psi], Fermion, {}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] + 
   (((I/2)*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0])/Coupling[M, {}, 0]^2 + 
     hbar*(((-I)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2)/
        (\[Epsilon]*Coupling[M, {}, 0]^2) + 
       (((-5*I)/4)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2 - 
         I*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2))*
    (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
     Field[\[Phi], Scalar, {}, {}]^2*
      Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[\[Psi], Fermion, {}, {}]) + 
   (hbar*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
     (-(FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
         {}]*Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz], Index[d$$2, Lorentz]], 
          GammaM[Index[d$$3, Lorentz]], Proj[-1]]**Field[\[Psi], Fermion, {}, 
          {Index[d$$3, Lorentz]}]) + 
      FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
       Bar[Field[\[Psi], Fermion, {}, {Index[d$$3, Lorentz]}]]**
        DiracProduct[GammaM[Index[d$$3, Lorentz]], 
         GammaM[Index[d$$1, Lorentz], Index[d$$2, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}]))/(8*Coupling[M, {}, 0]^2) + 
   ((I/6)*hbar*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
     (Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[\[Psi], Fermion, {}, {Index[d$$2, Lorentz], 
         Index[d$$2, Lorentz]}] - Bar[Field[\[Psi], Fermion, {}, 
         {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]**
       DiracProduct[GammaM[Index[d$$2, Lorentz]], Proj[-1]]**
       Field[\[Psi], Fermion, {}, {Index[d$$2, Lorentz]}]))/
    Coupling[M, {}, 0]^2, "On-shell EFT Lagrangian" -> 
  -1/2*(Coupling[m, {}, 1]^2*Field[\[Phi], Scalar, {}, {}]^2) + 
   (hbar*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*Coupling[y, {}, 0]*
     Field[\[Phi], Scalar, {}, {}]^2)/2 + 
   (hbar*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*Coupling[y, {}, 0]*
     Field[\[Phi], Scalar, {}, {}]^2)/\[Epsilon] + 
   (hbar*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^4*Coupling[y, {}, 0]*
     Field[\[Phi], Scalar, {}, {}]^2)/(3*Coupling[M, {}, 0]^2) - 
   (hbar*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
     Field[\[Phi], Scalar, {}, {}]^4)/\[Epsilon] - 
   (13*hbar*Bar[Coupling[y, {}, 0]]^2*Coupling[m, {}, 1]^2*
     Coupling[y, {}, 0]^2*Field[\[Phi], Scalar, {}, {}]^4)/
    (18*Coupling[M, {}, 0]^2) + (hbar*Bar[Coupling[y, {}, 0]]^3*
     Coupling[y, {}, 0]^3*Field[\[Phi], Scalar, {}, {}]^6)/
    (3*Coupling[M, {}, 0]^2) + 
   Field[\[Phi], Scalar, {}, {Index[d$$1, Lorentz]}]^2/2 - 
   (hbar*FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
       {}]^2)/(3*\[Epsilon]) - 
   FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2/
    (4*Coupling[e, {}, 0]^2) + (hbar*Bar[Coupling[y, {}, 0]]*
     Coupling[y, {}, 0]*Field[\[Phi], Scalar, {}, {}]^2*
     FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
      2)/(3*Coupling[M, {}, 0]^2) + hbar*Bar[Coupling[y, {}, 0]]*
    Coupling[m, {}, 1]^2*Coupling[y, {}, 0]*Field[\[Phi], Scalar, {}, {}]^2*
    Log[\[Mu]bar2/Coupling[M, {}, 0]^2] - hbar*Bar[Coupling[y, {}, 0]]^2*
    Coupling[y, {}, 0]^2*Field[\[Phi], Scalar, {}, {}]^4*
    Log[\[Mu]bar2/Coupling[M, {}, 0]^2] - 
   (hbar*FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
       {}]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/3 + 
   hbar*(-2*Bar[Coupling[y, {}, 0]]*Coupling[M, {}, 0]^2*Coupling[y, {}, 0]*
      Field[\[Phi], Scalar, {}, {}]^2 - (2*Bar[Coupling[y, {}, 0]]*
       Coupling[M, {}, 0]^2*Coupling[y, {}, 0]*Field[\[Phi], Scalar, {}, {}]^
        2)/\[Epsilon] - 2*Bar[Coupling[y, {}, 0]]*Coupling[M, {}, 0]^2*
      Coupling[y, {}, 0]*Field[\[Phi], Scalar, {}, {}]^2*
      Log[\[Mu]bar2/Coupling[M, {}, 0]^2]) + 
   I*Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]]]**Field[\[Psi], Fermion, {}, 
      {Index[d$$1, Lorentz]}] - (2*hbar*Coupling[e, {}, 0]^4*
     (Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
        GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[\[Psi], Fermion, {}, 
        {}])^2)/(15*Coupling[M, {}, 0]^2) + 
   (7*hbar*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*Coupling[y, {}, 0]*
     (Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
        GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[\[Psi], Fermion, {}, 
        {}])^2)/(36*Coupling[M, {}, 0]^2) - 
   (4*hbar*Coupling[e, {}, 0]^4*Bar[Field[\[Psi], Fermion, {}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[\[Psi], Fermion, {}, {}]*Bar[Field[\[Psi], Fermion, {}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[\[Psi], Fermion, {}, {}])/(15*Coupling[M, {}, 0]^2) + 
   (7*hbar*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*Coupling[y, {}, 0]*
     Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[\[Psi], Fermion, {}, 
       {}]*Bar[Field[\[Psi], Fermion, {}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[\[Psi], Fermion, {}, {}])/(36*Coupling[M, {}, 0]^2) - 
   (2*hbar*Coupling[e, {}, 0]^4*(Bar[Field[\[Psi], Fermion, {}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[\[Psi], Fermion, {}, {}])^2)/(15*Coupling[M, {}, 0]^2), 
 "SuperTraces" -> <|"hFermion-lScalar" -> 
    (7*hbar*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
       FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
        {Index[d$$2, Lorentz]}]*Bar[Field[\[Psi], Fermion, {}, {}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}])/(36*Coupling[M, {}, 0]^2) + 
     hbar*(((3*I)/4)*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0] + 
       ((I/2)*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0])/\[Epsilon] + 
       (I/2)*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2])*
      Bar[Field[\[Psi], Fermion, {}, {}]]**DiracProduct[
        GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[\[Psi], Fermion, {}, 
        {Index[d$$1, Lorentz]}] + (hbar*Bar[Coupling[y, {}, 0]]*
       Coupling[y, {}, 0]*(-(FieldStrength[A, {Index[d$$1, Lorentz], 
            Index[d$$2, Lorentz]}, {}, {}]*Bar[Field[\[Psi], Fermion, {}, 
             {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz], 
             Index[d$$2, Lorentz]], GammaM[Index[d$$3, Lorentz]], Proj[-1]]**
           Field[\[Psi], Fermion, {}, {Index[d$$3, Lorentz]}]) + 
        FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
          {}]*Bar[Field[\[Psi], Fermion, {}, {Index[d$$3, Lorentz]}]]**
          DiracProduct[GammaM[Index[d$$3, Lorentz]], 
           GammaM[Index[d$$1, Lorentz], Index[d$$2, Lorentz]], Proj[-1]]**
          Field[\[Psi], Fermion, {}, {}]))/(8*Coupling[M, {}, 0]^2) + 
     ((I/6)*hbar*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
       (Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[\[Psi], Fermion, {}, {Index[d$$2, Lorentz], 
           Index[d$$2, Lorentz]}] - Bar[Field[\[Psi], Fermion, {}, 
           {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]**
         DiracProduct[GammaM[Index[d$$2, Lorentz]], Proj[-1]]**
         Field[\[Psi], Fermion, {}, {Index[d$$2, Lorentz]}]))/
      Coupling[M, {}, 0]^2, "hFermion-lFermion" -> 
    (hbar*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*Field[\[Phi], Scalar, 
        {}, {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]*
       Field[\[Phi], Scalar, {}, {Index[d$$2, Lorentz], 
         Index[d$$2, Lorentz]}])/(3*Coupling[M, {}, 0]^2) + 
     (hbar*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
       Field[\[Phi], Scalar, {}, {}]^2*
       FieldStrength[A, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
        2)/(3*Coupling[M, {}, 0]^2) - 
     hbar*Field[\[Phi], Scalar, {}, {Index[d$$1, Lorentz]}]^2*
      (-1/2*(Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]) - 
       (Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0])/\[Epsilon] - 
       Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2]) + 
     hbar*Field[\[Phi], Scalar, {}, {}]^2*
      ((-2*Bar[Coupling[y, {}, 0]]*Coupling[M, {}, 0]^2*Coupling[y, {}, 0])/
        \[Epsilon] + Coupling[M, {}, 0]^2*(-2*Bar[Coupling[y, {}, 0]]*
          Coupling[y, {}, 0] - 2*Bar[Coupling[y, {}, 0]]*Coupling[y, {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])), 
   "hFermion-lVector" -> hbar*((((-7*I)/2)*Bar[Coupling[y, {}, 0]]*
        Coupling[e, {}, 0]^2*Coupling[y, {}, 0])/
       (\[Epsilon]*Coupling[M, {}, 0]^2) + 
      (((-7*I)/4)*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*
         Coupling[y, {}, 0] - ((7*I)/2)*Bar[Coupling[y, {}, 0]]*
         Coupling[e, {}, 0]^2*Coupling[y, {}, 0]*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2)*
     (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
      Field[\[Phi], Scalar, {}, {}]^2*
       Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}]), "hFermion-lScalar-lScalar" -> 
    hbar*((I*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*Coupling[y, {}, 0])/
       (\[Epsilon]*Coupling[M, {}, 0]^2) + 
      (((3*I)/2)*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*
         Coupling[y, {}, 0] + I*Bar[Coupling[y, {}, 0]]*Coupling[m, {}, 1]^2*
         Coupling[y, {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
       Coupling[M, {}, 0]^2)*Bar[Field[\[Psi], Fermion, {}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}], 
   "hFermion-lScalar-lFermion" -> 
    hbar*(((-1/2*I)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2)/
       (\[Epsilon]*Coupling[M, {}, 0]^2) + 
      ((-1/2*I)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2 - 
        (I/2)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2)*
     (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
      Field[\[Phi], Scalar, {}, {}]^2*
       Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}]), "hFermion-lFermion-lScalar" -> 
    hbar*(((-1/2*I)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2)/
       (\[Epsilon]*Coupling[M, {}, 0]^2) + 
      ((-1/2*I)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2 - 
        (I/2)*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2)*
     (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
      Field[\[Phi], Scalar, {}, {}]^2*
       Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}]), "hFermion-lFermion-lVector" -> 
    hbar*((((3*I)/2)*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*
        Coupling[y, {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^2) + 
      (((3*I)/4)*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*
         Coupling[y, {}, 0] + ((3*I)/2)*Bar[Coupling[y, {}, 0]]*
         Coupling[e, {}, 0]^2*Coupling[y, {}, 0]*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2)*
     (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
      Field[\[Phi], Scalar, {}, {}]^2*
       Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}]), "hFermion-lVector-lFermion" -> 
    hbar*((((3*I)/2)*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*
        Coupling[y, {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^2) + 
      (((3*I)/4)*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*
         Coupling[y, {}, 0] + ((3*I)/2)*Bar[Coupling[y, {}, 0]]*
         Coupling[e, {}, 0]^2*Coupling[y, {}, 0]*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2)*
     (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
      Field[\[Phi], Scalar, {}, {}]^2*
       Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}]), 
   "hFermion-lScalar-hFermion-lScalar" -> 0, 
   "hFermion-lScalar-hFermion-lFermion" -> 
    ((-1/4*I)*hbar*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
      (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
       Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[\[Psi], Fermion, {}, 
          {}]))/Coupling[M, {}, 0]^2, 
   "hFermion-lFermion-hFermion-lFermion" -> 
    (13*hbar*Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
       Field[\[Phi], Scalar, {}, {}]^3*Field[\[Phi], Scalar, {}, 
        {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}])/
      (18*Coupling[M, {}, 0]^2) + hbar*Field[\[Phi], Scalar, {}, {}]^4*
      (-((Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2)/\[Epsilon]) - 
       Bar[Coupling[y, {}, 0]]^2*Coupling[y, {}, 0]^2*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2]), 
   "hFermion-lFermion-lVector-lFermion" -> 
    hbar*(((I/2)*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*
        Coupling[y, {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^2) + 
      ((I/4)*Bar[Coupling[y, {}, 0]]*Coupling[e, {}, 0]^2*
         Coupling[y, {}, 0] + (I/2)*Bar[Coupling[y, {}, 0]]*
         Coupling[e, {}, 0]^2*Coupling[y, {}, 0]*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2)*
     (Field[\[Phi], Scalar, {}, {}]^2*Bar[Field[\[Psi], Fermion, {}, {}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}] - 
      Field[\[Phi], Scalar, {}, {}]^2*
       Bar[Field[\[Psi], Fermion, {}, {Index[d$$1, Lorentz]}]]**
        DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[\[Psi], Fermion, {}, {}]), 
   "hFermion-lFermion-hFermion-lFermion-hFermion-lFermion" -> 
    (hbar*Bar[Coupling[y, {}, 0]]^3*Coupling[y, {}, 0]^3*
      Field[\[Phi], Scalar, {}, {}]^6)/(3*Coupling[M, {}, 0]^2)|>, 
 "Matching Conditions" -> None|>
