(* Created with the Wolfram Language : www.wolfram.com *)
<|"Model" -> "Singlet_Scalar_Extension", "Version" -> "0.3.0", 
 "Date" -> DateObject[{2025, 5, 30}, "Day"], "Time (Match)" -> 4.765902, 
 "Time (GreensSimplify)" -> 0.276364, "Time (EOMSimplify)" -> 3.927217, 
 "Time (MapEffectiveCouplings)" -> 12.86903, 
 "UV Lagrangian" -> Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}] + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}] - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/2 - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[\[Phi], Scalar, {}, {}] - 
   (Coupling[M, {}, 0]^2*Field[\[Phi], Scalar, {}, {}]^2)/2 - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[\[Phi], Scalar, {}, {}]^2)/2 - 
   (Coupling[\[Mu]\[Phi], {}, 0]*Field[\[Phi], Scalar, {}, {}]^3)/6 - 
   (Coupling[\[Lambda]\[Phi], {}, 0]*Field[\[Phi], Scalar, {}, {}]^4)/24 + 
   Field[\[Phi], Scalar, {}, {Index[d$$1, Lorentz]}]^2/2 - 
   FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2/
    (4*Coupling[gY, {}, 0]^2) - 
   FieldStrength[G, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU3c[adj]]}, {}]^2/(4*Coupling[gs, {}, 0]^2) - 
   FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]^2/(4*Coupling[gL, {}, 0]^2) + 
   I*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
      {Index[d$$1, Lorentz]}] - 
   Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   I*Bar[Field[e, Fermion, {Index[d$$1, Flavor]}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[e, Fermion, {Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
   Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] - 
   Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$2, Flavor]}, 
      {}] + I*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
      {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
      {Index[d$$1, Lorentz]}] - 
   Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}] + 
   I*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
    Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}] + 
   I*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
      {Index[d$$1, Lorentz]}] - 
   Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
      Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
     {Index[d$$2, SU2L[fund]]}, {}]*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}], "Off-shell EFT Lagrangian" -> 
  (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$2, Lorentz], 
       Index[d$$2, Lorentz]}])/(6*Coupling[M, {}, 0]^4) - 
   FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2/
    (4*Coupling[gY, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
      2)/(12*Coupling[M, {}, 0]^4) - 
   FieldStrength[G, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU3c[adj]]}, {}]^2/(4*Coupling[gs, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}])/(3*Coupling[M, {}, 0]^4) - 
   FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]^2/(4*Coupling[gL, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
       {Index[d$$1, SU2L[adj]]}, {}]^2)/(12*Coupling[M, {}, 0]^4) + 
   hbar*(-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {}]*FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
        {}, {Index[d$$2, Lorentz]}]) + 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
       {Index[d$$2, Lorentz]}])*(((-1/12*I)*Coupling[A, {}, 0]^2)/
      (\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((-17*I)/72)*Coupling[A, {}, 0]^2 - (I/12)*Coupling[A, {}, 0]^2*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
   hbar*(-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {Index[d$$1, Lorentz]}]]*CG[gen[SU2L[fund]], 
        {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
         Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]*FieldStrength[W, 
        {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
        {Index[d$$1, SU2L[adj]]}, {Index[d$$2, Lorentz]}]) + 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
        Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
       {Index[d$$1, SU2L[adj]]}, {Index[d$$2, Lorentz]}])*
    (((-1/6*I)*Coupling[A, {}, 0]^2)/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((-17*I)/36)*Coupling[A, {}, 0]^2 - (I/6)*Coupling[A, {}, 0]^2*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((-3*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/
      (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     ((-5*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/4 - 
       (3*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^4) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    (-(Coupling[A, {}, 0]^2/Coupling[M, {}, 0]^4) + 
     hbar*(-1/12*Coupling[\[Kappa], {}, 0]^2/Coupling[M, {}, 0]^2 + 
       ((2*Coupling[A, {}, 0]^4)/Coupling[M, {}, 0]^6 - 
         (Coupling[A, {}, 0]^2*(9*Coupling[gL, {}, 0]^2 + 
            3*Coupling[gY, {}, 0]^2 + 4*(-4*Coupling[\[Kappa], {}, 0] + 
              Coupling[\[Lambda]\[Phi], {}, 0])))/(4*Coupling[M, {}, 0]^4))/
        \[Epsilon] + ((Coupling[A, {}, 0]^2*(26*Coupling[A, {}, 0]^2 - 
            26*Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0] + 
            11*Coupling[\[Mu]\[Phi], {}, 0]^2))/12 + 2*Coupling[A, {}, 0]^4*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
       (-1/24*(Coupling[A, {}, 0]*(Coupling[A, {}, 0]*
             (45*Coupling[gL, {}, 0]^2 + 15*Coupling[gY, {}, 0]^2 - 
              100*Coupling[\[Kappa], {}, 0] + 24*Coupling[\[Lambda]\[Phi], 
                {}, 0]) + 20*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], 
              {}, 0])) - (Coupling[A, {}, 0]^2*(9*Coupling[gL, {}, 0]^2 + 
            3*Coupling[gY, {}, 0]^2 + 4*(-4*Coupling[\[Kappa], {}, 0] + 
              Coupling[\[Lambda]\[Phi], {}, 0]))*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4)/Coupling[M, {}, 0]^4)) + 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
         Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
        Index[d$$1, Lorentz]}])*
    (-1/2*Coupling[A, {}, 0]^2/Coupling[M, {}, 0]^4 + 
     hbar*(-1/24*Coupling[\[Kappa], {}, 0]^2/Coupling[M, {}, 0]^2 + 
       ((3*Coupling[A, {}, 0]^4)/(2*Coupling[M, {}, 0]^6) - 
         (Coupling[A, {}, 0]^2*(3*Coupling[gL, {}, 0]^2 + 
            3*Coupling[gY, {}, 0]^2 + 4*(-4*Coupling[\[Kappa], {}, 0] + 
              Coupling[\[Lambda], {}, 0] + Coupling[\[Lambda]\[Phi], {}, 
               0])))/(8*Coupling[M, {}, 0]^4))/\[Epsilon] + 
       ((Coupling[A, {}, 0]^2*(80*Coupling[A, {}, 0]^2 - 
            32*Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0] + 
            11*Coupling[\[Mu]\[Phi], {}, 0]^2))/24 + 
         (3*Coupling[A, {}, 0]^4*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/
        Coupling[M, {}, 0]^6 + 
       (-1/48*(Coupling[A, {}, 0]*(Coupling[A, {}, 0]*
             (15*Coupling[gL, {}, 0]^2 + 15*Coupling[gY, {}, 0]^2 - 
              136*Coupling[\[Kappa], {}, 0] + 108*Coupling[\[Lambda], {}, 
                0] + 24*Coupling[\[Lambda]\[Phi], {}, 0]) + 
            20*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])) - 
         (Coupling[A, {}, 0]^2*(3*Coupling[gL, {}, 0]^2 + 
            3*Coupling[gY, {}, 0]^2 + 4*(-4*Coupling[\[Kappa], {}, 0] + 
              Coupling[\[Lambda], {}, 0] + Coupling[\[Lambda]\[Phi], {}, 0]))*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/8)/Coupling[M, {}, 0]^4)) - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
    (-1 + hbar*(-1/2*Coupling[A, {}, 0]^2/Coupling[M, {}, 0]^2 + 
       (Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((5*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2])/2 + 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4)) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    (Coupling[\[Mu]2, {}, 2] + hbar*(Coupling[A, {}, 0]^2 - 
       (Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/2 + 
       (Coupling[A, {}, 0]^2 + (Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 
            0])/2 - (Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2])/
          Coupling[M, {}, 0]^2 + (Coupling[A, {}, 0]^2*
           Coupling[\[Mu]2, {}, 2]^2)/Coupling[M, {}, 0]^4 - 
         (Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/2)/\[Epsilon] + 
       (Coupling[A, {}, 0]^2 - (Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 
            0])/2)*Log[\[Mu]bar2/Coupling[M, {}, 0]^2] + 
       Coupling[M, {}, 0]^2*(Coupling[\[Kappa], {}, 0]/2 + 
         (Coupling[\[Kappa], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2) + 
       (-(Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]) - 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2 + 
       (Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2 + 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4)) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    (Coupling[A, {}, 0]^2/(2*Coupling[M, {}, 0]^2) - 
     Coupling[\[Lambda], {}, 0]/2 + 
     hbar*((Coupling[\[Kappa], {}, 0]^2/4 + (Coupling[A, {}, 0]^2*
           (-8*Coupling[\[Kappa], {}, 0] + 12*Coupling[\[Lambda], {}, 0] + 
            Coupling[\[Lambda]\[Phi], {}, 0]))/(4*Coupling[M, {}, 0]^2) + 
         (-2*Coupling[A, {}, 0]^4 + 3*Coupling[A, {}, 0]^2*
            (Coupling[\[Kappa], {}, 0] - 2*Coupling[\[Lambda], {}, 0])*
            Coupling[\[Mu]2, {}, 2])/Coupling[M, {}, 0]^4 + 
         (Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*
           (6*Coupling[A, {}, 0] - Coupling[\[Mu]\[Phi], {}, 0]))/
          Coupling[M, {}, 0]^6)/\[Epsilon] + 
       (Coupling[\[Kappa], {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4 + 
       (-1/4*(Coupling[A, {}, 0]*(Coupling[A, {}, 0]*(12*Coupling[\[Kappa], 
                {}, 0] - 12*Coupling[\[Lambda], {}, 0] - Coupling[
               \[Lambda]\[Phi], {}, 0]) - 2*Coupling[\[Kappa], {}, 0]*
             Coupling[\[Mu]\[Phi], {}, 0])) + (Coupling[A, {}, 0]^2*
           (-8*Coupling[\[Kappa], {}, 0] + 12*Coupling[\[Lambda], {}, 0] + 
            Coupling[\[Lambda]\[Phi], {}, 0])*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4)/Coupling[M, {}, 0]^2 + 
       (-1/4*(Coupling[A, {}, 0]^2*(12*Coupling[A, {}, 0]^2 - 
            16*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 2] + 
            24*Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]2, {}, 2] - 
            4*Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0] + 
            Coupling[\[Mu]\[Phi], {}, 0]^2)) + (-2*Coupling[A, {}, 0]^4 + 
           3*Coupling[A, {}, 0]^2*(Coupling[\[Kappa], {}, 0] - 
             2*Coupling[\[Lambda], {}, 0])*Coupling[\[Mu]2, {}, 2])*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
       (2*Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*
          (4*Coupling[A, {}, 0] - Coupling[\[Mu]\[Phi], {}, 0]) + 
         Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*(6*Coupling[A, {}, 0] - 
           Coupling[\[Mu]\[Phi], {}, 0])*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^6)) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
    (-1/2*(Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
       Coupling[M, {}, 0]^4 + (Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 
        0])/(6*Coupling[M, {}, 0]^6) + 
     hbar*(-1/12*Coupling[\[Kappa], {}, 0]^3/Coupling[M, {}, 0]^2 + 
       ((Coupling[A, {}, 0]^2*(12*Coupling[\[Kappa], {}, 0]^2 + 
            36*Coupling[\[Lambda], {}, 0]^2 - Coupling[\[Kappa], {}, 0]*
             (36*Coupling[\[Lambda], {}, 0] + Coupling[\[Lambda]\[Phi], {
                }, 0])))/(4*Coupling[M, {}, 0]^4) + 
         (9*Coupling[A, {}, 0]^6 - (7*Coupling[A, {}, 0]^5*
             Coupling[\[Mu]\[Phi], {}, 0])/2)/Coupling[M, {}, 0]^8 + 
         (Coupling[A, {}, 0]^3*(18*Coupling[A, {}, 0]*
             (7*Coupling[\[Kappa], {}, 0] - 12*Coupling[\[Lambda], {}, 0]) + 
            (-12*Coupling[\[Kappa], {}, 0] + 36*Coupling[\[Lambda], {}, 0] + 
              Coupling[\[Lambda]\[Phi], {}, 0])*Coupling[\[Mu]\[Phi], {}, 
              0]))/(12*Coupling[M, {}, 0]^6))/\[Epsilon] + 
       ((Coupling[A, {}, 0]*(Coupling[A, {}, 0]*(22*Coupling[\[Kappa], {}, 0]^
                2 + 36*Coupling[\[Lambda], {}, 0]^2 - 2*Coupling[\[Kappa], 
                {}, 0]*(24*Coupling[\[Lambda], {}, 0] + Coupling[
                 \[Lambda]\[Phi], {}, 0])) - Coupling[\[Kappa], {}, 0]^2*
             Coupling[\[Mu]\[Phi], {}, 0]))/4 + 
         (Coupling[A, {}, 0]^2*(12*Coupling[\[Kappa], {}, 0]^2 + 
            36*Coupling[\[Lambda], {}, 0]^2 - Coupling[\[Kappa], {}, 0]*
             (36*Coupling[\[Lambda], {}, 0] + Coupling[\[Lambda]\[Phi], {
                }, 0]))*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4)/
        Coupling[M, {}, 0]^4 + 
       ((Coupling[A, {}, 0]^3*(86*Coupling[A, {}, 0]^3 - 
            45*Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0] + 
            6*Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]^2 - 
            Coupling[\[Mu]\[Phi], {}, 0]^3))/6 + 
         (9*Coupling[A, {}, 0]^6 - (7*Coupling[A, {}, 0]^5*
             Coupling[\[Mu]\[Phi], {}, 0])/2)*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8 + 
       ((Coupling[A, {}, 0]^2*(3*Coupling[A, {}, 0]^2*
             (37*Coupling[\[Kappa], {}, 0] - 48*Coupling[\[Lambda], {}, 0] - 
              Coupling[\[Lambda]\[Phi], {}, 0]) + 2*Coupling[A, {}, 0]*
             (-15*Coupling[\[Kappa], {}, 0] + 18*Coupling[\[Lambda], {}, 0] + 
              Coupling[\[Lambda]\[Phi], {}, 0])*Coupling[\[Mu]\[Phi], {}, 
              0] + 3*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]^
              2))/6 + (Coupling[A, {}, 0]^3*(18*Coupling[A, {}, 0]*
             (7*Coupling[\[Kappa], {}, 0] - 12*Coupling[\[Lambda], {}, 0]) + 
            (-12*Coupling[\[Kappa], {}, 0] + 36*Coupling[\[Lambda], {}, 0] + 
              Coupling[\[Lambda]\[Phi], {}, 0])*Coupling[\[Mu]\[Phi], {}, 0])*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/12)/Coupling[M, {}, 0]^6)) - 
   Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
         Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
          Index[d$$3, Flavor]}, 0] + 
       Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
          Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   hbar*(((I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((5*I)/8)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
          0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
        {}] - Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
        {}]) + I*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[d, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]]*CG[eps[SU2L], 
     {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
    (((-1/2*I)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((-5*I)/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
          0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] - 
       (I/2)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   hbar*(((I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     ((I/8)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$1, Lorentz]}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {Index[d$$1, Lorentz]}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
        Proj[1]]**Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$1, Flavor]}, {}]) - 
   Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
         Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
          Index[d$$3, Flavor]}, 0] + 
       Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
          Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[l, Fermion, 
      {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}] + 
   hbar*(((I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((5*I)/8)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
          0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[e, Fermion, {Index[d$$1, Flavor]}, {}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[e, Fermion, {Index[d$$1, Flavor]}, {}]) + 
   I*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[e, Fermion, {Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}] + 
   hbar*(((I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     ((I/8)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[e, Fermion, {Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[e, Fermion, {Index[d$$1, Flavor]}, {}]) - 
   Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
      {}] + I*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
      {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
         Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
          Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
          Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
      {}] + 
   hbar*(((-1/4*I)*Bar[Coupling[Ye, {Index[d$$1, Flavor], 
          Index[d$$3, Flavor]}, 0]]*Coupling[A, {}, 0]^2*
       Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((-5*I)/8)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
          0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0] - 
       (I/4)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
        {}] - Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
        {}]) + 
   hbar*(((I/4)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
         Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     ((I/8)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$1, Lorentz]}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {Index[d$$1, Lorentz]}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
        Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}]) - 
   Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   I*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
         Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
          Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
          Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
    Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   hbar*(((-1/4*I)*Coupling[A, {}, 0]^2*
       (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] + 
        Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]))/
      (\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((-5*I)/8)*Coupling[A, {}, 0]^2*
        (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] + 
         Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]) - 
       (I/4)*Coupling[A, {}, 0]^2*
        (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] + 
         Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}]) + 
   hbar*(((I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
         Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((5*I)/8)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
          0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}]) + 
   hbar*(((I/4)*Coupling[A, {}, 0]^2*
       (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] - 
        Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]))/
      (\[Epsilon]*Coupling[M, {}, 0]^4) + 
     ((I/8)*Coupling[A, {}, 0]^2*
        (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] - 
         Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]) + 
       (I/4)*Coupling[A, {}, 0]^2*
        (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] - 
         Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])*
        Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}]) + 
   hbar*(((I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
         Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     ((I/8)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
          Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}]) + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
         Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
          Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
          Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] - 
   Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
      Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
     {Index[d$$2, SU2L[fund]]}, {}]*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
      Bar[Index[d$$3, SU2L[fund]]]}]*Field[H, Scalar, 
     {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
     {Index[d$$3, SU2L[fund]]}, {}]*
    ((Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$4, Flavor], 
         Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$4, Flavor], 
          Index[d$$3, Flavor]}, 0] + 
       Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$4, Flavor], 
          Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   hbar*CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
      Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
     {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
     {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
    (((I/2)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((5*I)/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
          0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       (I/2)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   hbar*(((-1/4*I)*Bar[Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$2, Flavor]}, 0]]*Coupling[A, {}, 0]^2*
       Coupling[Yu, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
      (\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (((-5*I)/8)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
          0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] - 
       (I/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
        {}] - Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
        {}]) + I*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[u, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] + 
   hbar*(((I/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
         0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
         Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     ((I/8)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
          Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$1, Lorentz]}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
         {Index[d$$1, Lorentz]}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
        Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$1, Flavor]}, {}]), "On-shell EFT Lagrangian" -> 
  Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}] - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/
    (2*Coupling[M, {}, 0]^2) + 
   (8*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}] - 
   (7*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (2*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (2*Coupling[M, {}, 0]^2) - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^2 - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/(4*\[Epsilon]) - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/2 + 
   (7*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (2*Coupling[M, {}, 0]^2) + 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (4*Coupling[M, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (4*\[Epsilon]*Coupling[M, {}, 0]^2) + 
   (19*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^6 + 
   (10*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^6) - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^4 - 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (36*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (36*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (29*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^4) + 
   (7*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^2) - 
   (40*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^4) - 
   (8*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^4 - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (2*Coupling[M, {}, 0]^2) - 
   (14*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*Coupling[\[Mu]\[Phi], {}, 
      0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*Coupling[\[Mu]\[Phi], {}, 
      0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^6) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 2]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (4*Coupling[M, {}, 0]^4) + 
   (11*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Coupling[\[Mu]\[Phi], {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}])/(12*Coupling[M, {}, 0]^6) + 
   (19*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^4*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^6) + 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^4*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^6) - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/Coupling[M, {}, 0]^4 - 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(12*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(36*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(6*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (18*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (3*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (25*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(6*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[\[Kappa], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^2) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^4 - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (13*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^3*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(6*Coupling[M, {}, 0]^6) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^4) + 
   (11*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[\[Mu]\[Phi], {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^6) + 
   (175*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^6*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^8) + 
   (16*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^6*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^8) - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^6 - 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(36*Coupling[M, {}, 0]^6) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(6*\[Epsilon]*Coupling[M, {}, 0]^6) - 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(36*Coupling[M, {}, 0]^6) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(6*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (371*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^6) + 
   (37*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^6) - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   (65*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]^3*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^2) - 
   (87*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (2*Coupling[M, {}, 0]^6) - 
   (28*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^4 + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (36*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (36*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (53*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^4) - 
   (13*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^2) + 
   (41*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^4) + 
   (10*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^6 - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(2*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(4*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (149*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^8) - 
   (7*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^8) + 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^6) - 
   (41*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (4*Coupling[M, {}, 0]^4) + 
   (26*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^6) + 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Lambda]\[Phi], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (3*Coupling[M, {}, 0]^6) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Lambda]\[Phi], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (12*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (29*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Mu]\[Phi], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (12*Coupling[M, {}, 0]^8) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(2*Coupling[M, {}, 0]^6) - 
   (11*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}])/(12*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]^3*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*Coupling[M, {}, 0]^8) - 
   FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2/
    (4*Coupling[gY, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
      2)/(12*Coupling[M, {}, 0]^4) - 
   FieldStrength[G, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU3c[adj]]}, {}]^2/(4*Coupling[gs, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}])/(3*Coupling[M, {}, 0]^4) - 
   FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]^2/(4*Coupling[gL, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
       {Index[d$$1, SU2L[adj]]}, {}]^2)/(12*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4 + 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^2) + 
   (10*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) + 
   (7*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (8*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*Coupling[\[Mu]\[Phi], {}, 
      0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^4*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(3*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
   (16*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^6*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8 - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (6*Coupling[M, {}, 0]^6) - 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (6*Coupling[M, {}, 0]^6) + 
   (37*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4) - 
   (28*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) - 
   (13*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
   (10*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (4*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^4 - 
   (7*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^8) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Lambda]\[Phi], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(12*Coupling[M, {}, 0]^6) + 
   hbar*(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}] + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/
      \[Epsilon] + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/2 + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/(2*\[Epsilon]) - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^2 - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^2) + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/
      Coupling[M, {}, 0]^4 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {}]]*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/
      (\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}])/2 - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}])/(2*\[Epsilon]) - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      Coupling[M, {}, 0]^4 - (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      (\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      (2*Coupling[M, {}, 0]^2) - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      (2*\[Epsilon]*Coupling[M, {}, 0]^2) + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}])/Coupling[M, {}, 0]^6 + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      Coupling[M, {}, 0]^8 - (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      (\[Epsilon]*Coupling[M, {}, 0]^8) + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      (2*Coupling[M, {}, 0]^4) + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
      (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Log[\[Mu]bar2/Coupling[M, {}, 0]^2] + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2 - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {}]]*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      2 - (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^2) + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^6 - (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8 + 
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4)) - 
   Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[M, {}, 0]^2) - 
   (17*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (33*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(8*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (17*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]^2*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[M, {}, 0]^2) + 
   (29*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) - 
   (11*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) - 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (((5*I)/8)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[d, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/216)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/36)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((5*I)/8)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/216)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/36)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/36)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/36)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   I*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] - 
   (((5*I)/4)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/2)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*CG[eps[SU2L], 
      {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   ((I/2)*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]]*CG[eps[SU2L], 
      {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[M, {}, 0]^2) - 
   (17*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (33*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(8*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (17*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Kappa], {}, 0]^2*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[M, {}, 0]^2) + 
   (29*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) - 
   (11*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[M, {}, 0]^6) - 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) - 
   (3*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (((5*I)/8)*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((5*I)/8)*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   I*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[e, Fermion, {Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}] - 
   Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
      {}] + (hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (4*Coupling[M, {}, 0]^2) - (17*hbar*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*Coupling[M, {}, 0]^4) - 
   (hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[M, {}, 0]^4) + 
   I*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] + 
   (hbar*Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$1, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(6*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$1, Flavor]}, 0]*Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(6*Coupling[M, {}, 0]^4) - 
   (33*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(8*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(4*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (17*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(6*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Kappa], {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(24*Coupling[M, {}, 0]^2) + 
   (29*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(3*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*Coupling[M, {}, 0]^4) - 
   (11*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(24*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[M, {}, 0]^4) - 
   (((5*I)/8)*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((5*I)/8)*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   (((17*I)/144)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/144)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/144)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/144)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - Coupling[Yd, {Index[d$$2, Flavor], 
      Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
     {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
     DiracProduct[Proj[1]]**Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
       Index[d$$1, Flavor]}, {}] + (hbar*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (4*Coupling[M, {}, 0]^2) - (17*hbar*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   I*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}] + 
   (hbar*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
       Index[d$$1, Flavor]}, 0]*Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (6*Coupling[M, {}, 0]^4) + 
   (hbar*CG[eps[SU2L], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[Yu, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$2, Flavor]}, {}])/(6*Coupling[M, {}, 0]^4) - 
   (33*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (8*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (4*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (17*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (6*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Kappa], {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (24*Coupling[M, {}, 0]^2) + 
   (29*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^3*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (3*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) - 
   (11*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (24*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^4*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
    Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (4*Coupling[M, {}, 0]^2) - 
   (17*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) - 
   (((5*I)/8)*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((5*I)/8)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((5*I)/8)*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((5*I)/8)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/12)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   (((5*I)/8)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/144)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/432)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((5*I)/8)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], 
        Index[d$$3, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/144)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/432)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   ((I/24)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   ((I/72)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 - 
   (33*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^4*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(8*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^4*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (17*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (6*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Kappa], {}, 0]^2*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[M, {}, 0]^2) + 
   (29*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^3*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (3*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) - 
   (11*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0]^2*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (24*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^4*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Kappa], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
       Index[d$$1, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (6*Coupling[M, {}, 0]^4) + 
   (hbar*CG[eps[SU2L], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[Yu, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (6*Coupling[M, {}, 0]^4) - 
   Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
      Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
     {Index[d$$2, SU2L[fund]]}, {}]*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[M, {}, 0]^2) - 
   (17*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) - 
   (33*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^4*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(8*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^4*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) + 
   (Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (31*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (17*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[\[Kappa], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[M, {}, 0]^2) + 
   (29*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (4*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^3*
     Coupling[\[Mu]\[Phi], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) - 
   (11*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Mu]\[Phi], {}, 0]^2*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[M, {}, 0]^6) - 
   (3*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^4*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^6) + 
   (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[M, {}, 0]^4) + 
   (3*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) - 
   (2*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[M, {}, 0]^4 + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda]\[Phi], {}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[M, {}, 0]^4) + 
   (((5*I)/4)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*CG[Bar[eps[SU2L]], 
      {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]]}]*
     Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/2)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   ((I/2)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Coupling[Yd, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   (((5*I)/8)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[u, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((17*I)/108)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/18)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) + 
   (((5*I)/8)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   (((17*I)/108)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/18)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[M, {}, 0]^4) - 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
     Coupling[Yu, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/18)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   ((I/4)*hbar*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 - 
   ((I/18)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*Delta[Index[d$$1, Flavor], 
      Index[d$$2, Flavor]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[M, {}, 0]^4 + 
   I*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(6*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
       Index[d$$2, Flavor]}, 0]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(6*Coupling[M, {}, 0]^4), 
 "SuperTraces" -> 
  <|"hScalar" -> hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      (-1/2*(Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]) + 
       (Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 0] - 
         Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/(2*\[Epsilon]) - 
       (Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2 + Coupling[M, {}, 0]^2*
        (Coupling[\[Kappa], {}, 0]/2 + (Coupling[\[Kappa], {}, 0]*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (((Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2)/
          Coupling[M, {}, 0]^6 - (Coupling[A, {}, 0]*
           (Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/
          Coupling[M, {}, 0]^4)/\[Epsilon] + 
       (Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2 + 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
       (-(Coupling[A, {}, 0]*(Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], 
              {}, 0] + Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 
              0])) - Coupling[A, {}, 0]*(Coupling[A, {}, 0]*
            Coupling[\[Lambda]\[Phi], {}, 0] + Coupling[\[Kappa], {}, 0]*
            Coupling[\[Mu]\[Phi], {}, 0])*Log[\[Mu]bar2/Coupling[M, {}, 0]^
             2])/Coupling[M, {}, 0]^4) + 
     hbar*(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}])*
      (((Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2)/
          (2*Coupling[M, {}, 0]^6) - (Coupling[A, {}, 0]*
           (Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/
          (2*Coupling[M, {}, 0]^4))/\[Epsilon] + 
       ((Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2)/2 + 
         (Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^6 + 
       (-1/2*(Coupling[A, {}, 0]*(Coupling[A, {}, 0]*Coupling[
              \[Lambda]\[Phi], {}, 0] + Coupling[\[Kappa], {}, 0]*
             Coupling[\[Mu]\[Phi], {}, 0])) - 
         (Coupling[A, {}, 0]*(Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], 
              {}, 0] + Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 
              0])*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/
        Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((-1/4*(Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2)/
           Coupling[M, {}, 0]^4 + (Coupling[A, {}, 0]*
           (Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/
          (4*Coupling[M, {}, 0]^2))/\[Epsilon] + 
       (-1/4*(Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2) - 
         (Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4)/Coupling[M, {}, 0]^4 + 
       ((Coupling[A, {}, 0]*(Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 
              0] + 2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/
          4 + (Coupling[A, {}, 0]*(Coupling[A, {}, 0]*Coupling[
              \[Lambda]\[Phi], {}, 0] + 2*Coupling[\[Kappa], {}, 0]*
             Coupling[\[Mu]\[Phi], {}, 0])*Log[\[Mu]bar2/Coupling[M, {}, 0]^
              2])/4)/Coupling[M, {}, 0]^2) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      ((-1/4*(Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]^3)/
           Coupling[M, {}, 0]^8 - (Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 
            0]*(Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/
          (2*Coupling[M, {}, 0]^4) + (Coupling[A, {}, 0]^2*
           Coupling[\[Mu]\[Phi], {}, 0]*(4*Coupling[A, {}, 0]*
             Coupling[\[Lambda]\[Phi], {}, 0] + 9*Coupling[\[Kappa], {}, 0]*
             Coupling[\[Mu]\[Phi], {}, 0]))/(12*Coupling[M, {}, 0]^6))/
        \[Epsilon] + (-1/4*(Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]^
            3) - (Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]^3*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4)/Coupling[M, {}, 0]^8 + 
       (-1/2*(Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*
           (Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])) - 
         (Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*
           (Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^4 + 
       ((Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]*
           (4*Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            9*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/12 + 
         (Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]*
           (4*Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            9*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/12)/Coupling[M, {}, 0]^6), 
   "hScalar-hScalar" -> hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Coupling[\[Kappa], {}, 0]^2/4 - (Coupling[A, {}, 0]*
           Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/
          (2*Coupling[M, {}, 0]^2) + (Coupling[A, {}, 0]^2*
           Coupling[\[Mu]\[Phi], {}, 0]^2)/(4*Coupling[M, {}, 0]^4))/
        \[Epsilon] + (Coupling[\[Kappa], {}, 0]^2*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4 - 
       (Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], 
          {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        (2*Coupling[M, {}, 0]^2) + (Coupling[A, {}, 0]^2*
         Coupling[\[Mu]\[Phi], {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        (4*Coupling[M, {}, 0]^4)) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      (((Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]^3)/
          (4*Coupling[M, {}, 0]^8) + (Coupling[A, {}, 0]*Coupling[\[Kappa], 
            {}, 0]*(Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
            2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/
          (4*Coupling[M, {}, 0]^4) - (Coupling[A, {}, 0]^2*
           Coupling[\[Mu]\[Phi], {}, 0]*(Coupling[A, {}, 0]*
             Coupling[\[Lambda]\[Phi], {}, 0] + 3*Coupling[\[Kappa], {}, 0]*
             Coupling[\[Mu]\[Phi], {}, 0]))/(4*Coupling[M, {}, 0]^6))/
        \[Epsilon] + (Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]^3*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^8) + 
       (Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*
         (Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
          2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^4) - 
       (Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]*
         (Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0] + 
          3*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^6)) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/12*Coupling[\[Kappa], {}, 0]^2/Coupling[M, {}, 0]^2 + 
       ((Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], 
            {}, 0])/Coupling[M, {}, 0]^4 - (Coupling[A, {}, 0]^2*
           Coupling[\[Mu]\[Phi], {}, 0]^2)/Coupling[M, {}, 0]^6)/\[Epsilon] + 
       ((Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], 
            {}, 0])/6 + Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*
          Coupling[\[Mu]\[Phi], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4 + 
       (-1/12*(Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2) - 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6) + 
     hbar*(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}])*(-1/24*Coupling[\[Kappa], {}, 0]^2/
         Coupling[M, {}, 0]^2 + 
       ((Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], 
            {}, 0])/(2*Coupling[M, {}, 0]^4) - 
         (Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2)/
          (2*Coupling[M, {}, 0]^6))/\[Epsilon] + 
       ((Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], 
            {}, 0])/12 + (Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*
           Coupling[\[Mu]\[Phi], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
          2)/Coupling[M, {}, 0]^4 + 
       (-1/24*(Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2) - 
         (Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^6), 
   "hScalar-lScalar" -> 
    (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
      (2*Coupling[M, {}, 0]^2) + 
     (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*Coupling[A, {}, 0]^2*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$2, Lorentz], 
         Index[d$$2, Lorentz]}])/(6*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
        2)/(12*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
         Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[A, {}, 0]^2*
       Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
       FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
       FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
        {Index[d$$1, SU2L[adj]]}, {}])/(3*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[A, {}, 0]^2*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
       FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
         {Index[d$$1, SU2L[adj]]}, {}]^2)/(12*Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      (Coupling[A, {}, 0]^2 + Coupling[A, {}, 0]^2/\[Epsilon] + 
       Coupling[A, {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2]) + 
     hbar*(-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*FieldStrength[B, 
          {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
          {Index[d$$2, Lorentz]}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
         {Index[d$$2, Lorentz]}])*(((-1/12*I)*Coupling[A, {}, 0]^2)/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((-17*I)/72)*Coupling[A, {}, 0]^2 - (I/12)*Coupling[A, {}, 0]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*(-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*CG[gen[SU2L[fund]], 
          {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
           Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
          {Index[d$$2, SU2L[fund]]}, {}]*FieldStrength[W, 
          {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
          {Index[d$$1, SU2L[adj]]}, {Index[d$$2, Lorentz]}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
          Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
         {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
         {Index[d$$1, SU2L[adj]]}, {Index[d$$2, Lorentz]}])*
      (((-1/6*I)*Coupling[A, {}, 0]^2)/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((-17*I)/36)*Coupling[A, {}, 0]^2 - (I/6)*Coupling[A, {}, 0]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((-2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^2) + 
       (-2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0] - 
         2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2) + 
     hbar*(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}])*
      ((2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((5*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/2 + 
         2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((4*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (4*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0] + 
         4*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      (((3*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2)/
          Coupling[M, {}, 0]^4 - (Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 
            0]*Coupling[\[Mu]\[Phi], {}, 0])/Coupling[M, {}, 0]^6)/
        \[Epsilon] + (3*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2 + 
         3*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
       (-(Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 0]*
           Coupling[\[Mu]\[Phi], {}, 0]) - Coupling[A, {}, 0]^3*
          Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6), 
   "hScalar-hScalar-hScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     (-1/12*Coupling[\[Kappa], {}, 0]^3/Coupling[M, {}, 0]^2 + 
      (Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]^2*Coupling[\[Mu]\[Phi], 
         {}, 0])/(4*Coupling[M, {}, 0]^4) - 
      (Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
        Coupling[\[Mu]\[Phi], {}, 0]^2)/(4*Coupling[M, {}, 0]^6) + 
      (Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0]^3)/
       (12*Coupling[M, {}, 0]^8))*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}], 
   "hScalar-hScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      (-((Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
         Coupling[M, {}, 0]^2) + (Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], 
          {}, 0])/Coupling[M, {}, 0]^4)*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*
      ((Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
        (6*Coupling[M, {}, 0]^4) - (13*Coupling[A, {}, 0]^3*
         Coupling[\[Mu]\[Phi], {}, 0])/(6*Coupling[M, {}, 0]^6))*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
     hbar*((Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
        (3*Coupling[M, {}, 0]^4) - (4*Coupling[A, {}, 0]^3*
         Coupling[\[Mu]\[Phi], {}, 0])/(3*Coupling[M, {}, 0]^6))*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}]) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      ((2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2)/
        Coupling[M, {}, 0]^4 + (Coupling[A, {}, 0]^4*
         Coupling[\[Mu]\[Phi], {}, 0]^2)/(2*Coupling[M, {}, 0]^8) - 
       (Coupling[A, {}, 0]^3*(Coupling[A, {}, 0]*Coupling[\[Lambda]\[Phi], 
            {}, 0] + 6*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 
            0]))/(2*Coupling[M, {}, 0]^6))*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {}]*Field[H, Scalar, 
       {Index[d$$3, SU2L[fund]]}, {}], "hScalar-lScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (Coupling[A, {}, 0]^4/(\[Epsilon]*Coupling[M, {}, 0]^6) + 
       (Coupling[A, {}, 0]^4/2 + Coupling[A, {}, 0]^4*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6) + 
     hbar*(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}])*
      ((Coupling[A, {}, 0]^4/Coupling[M, {}, 0]^6 - 
         (Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0])/
          (2*Coupling[M, {}, 0]^4))/\[Epsilon] + 
       ((3*Coupling[A, {}, 0]^4)/2 + Coupling[A, {}, 0]^4*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
       ((-9*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0])/4 - 
         (Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      (-((Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2])/
         (\[Epsilon]*Coupling[M, {}, 0]^2)) + 
       (-(Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]) - 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2) - 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]*((Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 
          2])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((5*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2])/2 + 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (((3*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0])/
          Coupling[M, {}, 0]^2 + (-Coupling[A, {}, 0]^4 + 
           2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
            Coupling[\[Mu]2, {}, 2])/Coupling[M, {}, 0]^4)/\[Epsilon] + 
       (3*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0] + 
         3*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2 + 
       (-Coupling[A, {}, 0]^4 + 2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 
           0]*Coupling[\[Mu]2, {}, 2] + (-Coupling[A, {}, 0]^4 + 
           2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
            Coupling[\[Mu]2, {}, 2])*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      (((7*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0])/
          (2*Coupling[M, {}, 0]^6) - (6*Coupling[A, {}, 0]^2*
           Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda], {}, 0])/
          Coupling[M, {}, 0]^4 - (Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], 
            {}, 0])/(2*Coupling[M, {}, 0]^8))/\[Epsilon] + 
       ((7*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0])/2 + 
         (7*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0]*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^6 + 
       (-6*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
          Coupling[\[Lambda], {}, 0] - 6*Coupling[A, {}, 0]^2*
          Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
       (-1/2*(Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]) - 
         (Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^8), 
   "hScalar-hScalar-hScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     ((Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2)/
       (2*Coupling[M, {}, 0]^4) - (Coupling[A, {}, 0]^3*
        Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/
       Coupling[M, {}, 0]^6 + (Coupling[A, {}, 0]^4*
        Coupling[\[Mu]\[Phi], {}, 0]^2)/(2*Coupling[M, {}, 0]^8))*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}], 
   "hScalar-hScalar-lScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (((Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 
            2])/Coupling[M, {}, 0]^4 - (Coupling[A, {}, 0]^3*
           Coupling[\[Mu]2, {}, 2]*Coupling[\[Mu]\[Phi], {}, 0])/
          Coupling[M, {}, 0]^6)/\[Epsilon] + 
       (2*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 
           2] + Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
          Coupling[\[Mu]2, {}, 2]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4 + (-2*Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 
           2]*Coupling[\[Mu]\[Phi], {}, 0] - Coupling[A, {}, 0]^3*
          Coupling[\[Mu]2, {}, 2]*Coupling[\[Mu]\[Phi], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      (((-3*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
           Coupling[\[Lambda], {}, 0])/Coupling[M, {}, 0]^4 - 
         (Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0])/
          Coupling[M, {}, 0]^8 + (Coupling[A, {}, 0]^3*
           (Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0] + 
            3*Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]))/
          Coupling[M, {}, 0]^6)/\[Epsilon] + 
       (-6*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
          Coupling[\[Lambda], {}, 0] - 3*Coupling[A, {}, 0]^2*
          Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
       (-2*Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0] - 
         Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8 + 
       (2*Coupling[A, {}, 0]^3*(Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 
             0] + 3*Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]\[Phi], {}, 
             0]) + Coupling[A, {}, 0]^3*(Coupling[A, {}, 0]*
            Coupling[\[Kappa], {}, 0] + 3*Coupling[\[Lambda], {}, 0]*
            Coupling[\[Mu]\[Phi], {}, 0])*Log[\[Mu]bar2/Coupling[M, {}, 0]^
             2])/Coupling[M, {}, 0]^6), "hScalar-lScalar-hScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-(Coupling[A, {}, 0]^4/(\[Epsilon]*Coupling[M, {}, 0]^4)) + 
       (-2*Coupling[A, {}, 0]^4 - Coupling[A, {}, 0]^4*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}])*(Coupling[A, {}, 0]^4/
        (2*\[Epsilon]*Coupling[M, {}, 0]^6) + 
       ((11*Coupling[A, {}, 0]^4)/6 + (Coupling[A, {}, 0]^4*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^6) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (Coupling[A, {}, 0]^4/(\[Epsilon]*Coupling[M, {}, 0]^6) + 
       ((5*Coupling[A, {}, 0]^4)/3 + Coupling[A, {}, 0]^4*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      ((4*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^6) + 
       (8*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0] + 
         4*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6), 
   "hScalar-lScalar-lScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      ((Coupling[A, {}, 0]^6/Coupling[M, {}, 0]^8 - 
         (6*Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0])/
          Coupling[M, {}, 0]^6 + (9*Coupling[A, {}, 0]^2*
           Coupling[\[Lambda], {}, 0]^2)/Coupling[M, {}, 0]^4)/\[Epsilon] + 
       (Coupling[A, {}, 0]^6 + Coupling[A, {}, 0]^6*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8 + 
       (-6*Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0] - 
         6*Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
       (9*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]^2 + 
         9*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (((2*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2])/Coupling[M, {}, 0]^
           6 - (6*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
           Coupling[\[Mu]2, {}, 2])/Coupling[M, {}, 0]^4)/\[Epsilon] + 
       (2*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2] + 
         2*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
       (-6*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
          Coupling[\[Mu]2, {}, 2] - 6*Coupling[A, {}, 0]^2*
          Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      ((Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2)/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2 + 
         Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4), 
   "hScalar-lScalar-lFermion-lScalar" -> 
    hbar*(((I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
           0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((5*I)/8)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, 
              Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yd, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
          {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
         {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
          {}] - Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
          {}]) + hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*CG[eps[SU2L], 
       {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
      (((-1/2*I)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
           0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((-5*I)/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, 
              Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yu, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] - 
         (I/2)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
          Index[d$$2, Flavor]}, {}]]**DiracProduct[
        GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[u, Fermion, 
        {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
     hbar*(((I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, 
             Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yd, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((I/8)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[d, Fermion, 
          {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]) + 
     hbar*(((I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, 
             Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((5*I)/8)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, 
              Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
          {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
         {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[e, Fermion, {Index[d$$1, Flavor]}, {}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[e, Fermion, {Index[d$$1, Flavor]}, {}]) + 
     hbar*(((I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, 
             Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((I/8)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[e, Fermion, {Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[e, Fermion, 
          {Index[d$$1, Flavor]}, {}]) + 
     hbar*(((-1/4*I)*Bar[Coupling[Ye, {Index[d$$1, Flavor], 
            Index[d$$3, Flavor]}, 0]]*Coupling[A, {}, 0]^2*
         Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((-5*I)/8)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, 
              Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Ye, 
           {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] - 
         (I/4)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
            Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
          {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
         {Index[d$$2, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}] - Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}]) + 
     hbar*(((I/4)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, 
             Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Ye, 
          {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((I/8)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
            Index[d$$3, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
            Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]) + 
     hbar*(((-1/4*I)*Coupling[A, {}, 0]^2*
         (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] + 
          Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((-5*I)/8)*Coupling[A, {}, 0]^2*
          (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
            Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] + 
           Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
            Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]) - 
         (I/4)*Coupling[A, {}, 0]^2*
          (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
            Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] + 
           Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
            Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
          {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
         {Index[d$$2, SU2L[fund]]}, {}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
           Index[d$$1, Flavor]}, {}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
           Index[d$$1, Flavor]}, {}]) + 
     hbar*(((I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, 
             Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yu, 
          {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((5*I)/8)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, 
              Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yu, 
           {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
            Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
          {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
         {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
           Index[d$$1, Flavor]}, {}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
           Index[d$$1, Flavor]}, {}]) + 
     hbar*(((I/4)*Coupling[A, {}, 0]^2*
         (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] - 
          Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((I/8)*Coupling[A, {}, 0]^2*(Bar[Coupling[Yd, {Index[d$$1, Flavor], 
               Index[d$$3, Flavor]}, 0]]*Coupling[Yd, {Index[d$$2, Flavor], 
              Index[d$$3, Flavor]}, 0] - 
           Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
            Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]) + 
         (I/4)*Coupling[A, {}, 0]^2*
          (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
            Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0] - 
           Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
            Coupling[Yu, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
           Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
           Index[d$$1, Flavor]}, {}]) + 
     hbar*(((I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, 
             Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yu, 
          {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((I/8)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
            Index[d$$3, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
            Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
           Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
             SU2L[fund]], Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}]]**
         DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
           Index[d$$1, Flavor]}, {}]) + 
     hbar*CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
        Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
       {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      (((I/2)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
           0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((5*I)/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, 
              Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yd, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] + 
         (I/2)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
          Index[d$$2, Flavor]}, {}]]**DiracProduct[
        GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[d, Fermion, 
        {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
     hbar*(((-1/4*I)*Bar[Coupling[Yu, {Index[d$$3, Flavor], 
            Index[d$$2, Flavor]}, 0]]*Coupling[A, {}, 0]^2*
         Coupling[Yu, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (((-5*I)/8)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, 
              Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yu, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] - 
         (I/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
          {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
         {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
          {}] - Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
          {}]) + 
     hbar*(((I/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, 
             Flavor]}, 0]]*Coupling[A, {}, 0]^2*Coupling[Yu, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((I/8)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 
            0]]*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
         Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[u, Fermion, 
          {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]), 
   "hScalar-lScalar-lVector-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
       {}]*((-3*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/
        (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((-5*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/4 - 
         (3*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2)/Coupling[M, {}, 0]^4) + 
     hbar*(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}])*
      ((-3*Coupling[A, {}, 0]^2*(Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^
           2))/(8*\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((-5*Coupling[A, {}, 0]^2*(Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^
             2))/16 - (3*Coupling[A, {}, 0]^2*(Coupling[gL, {}, 0]^2 + 
            Coupling[gY, {}, 0]^2)*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/8)/
        Coupling[M, {}, 0]^4) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Bar[Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {}]]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((-3*Coupling[A, {}, 0]^2*(3*Coupling[gL, {}, 0]^2 + 
          Coupling[gY, {}, 0]^2))/(4*\[Epsilon]*Coupling[M, {}, 0]^4) + 
       ((-5*Coupling[A, {}, 0]^2*(3*Coupling[gL, {}, 0]^2 + 
            Coupling[gY, {}, 0]^2))/8 - (3*Coupling[A, {}, 0]^2*
           (3*Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/4)/Coupling[M, {}, 0]^4), 
   "hScalar-hScalar-lScalar-hScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     (((2*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0])/
         Coupling[M, {}, 0]^6 - (2*Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], 
           {}, 0])/Coupling[M, {}, 0]^8)/\[Epsilon] + 
      (5*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0] + 
        2*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0]*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
      (-5*Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0] - 
        2*Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8), 
   "hScalar-hScalar-lScalar-lVector-lScalar" -> 0, 
   "hScalar-lScalar-hScalar-lScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
      (((4*Coupling[A, {}, 0]^6)/Coupling[M, {}, 0]^8 - 
         (12*Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0])/
          Coupling[M, {}, 0]^6)/\[Epsilon] + 
       (6*Coupling[A, {}, 0]^6 + 4*Coupling[A, {}, 0]^6*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8 + 
       (-18*Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0] - 
         12*Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((4*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2])/
        (\[Epsilon]*Coupling[M, {}, 0]^6) + 
       (6*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2] + 
         4*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6), 
   "hScalar-lScalar-lScalar-lVector-lScalar" -> 0, 
   "hScalar-lScalar-lFermion-lFermion-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
           Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
            Index[d$$3, Flavor]}, 0] + 
         Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Yd, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$4, Flavor], 
            Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
          Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
           Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
            Index[d$$3, Flavor]}, 0] + 
         Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$4, Flavor], 
            Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[-1]]**Field[l, Fermion, 
        {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
           Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
            Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0] + 
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$2, Flavor], 
            Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
          Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
       Field[e, Fermion, {Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
           Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
            Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0] + 
         Bar[Coupling[Yd, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
            Index[d$$4, Flavor]}, 0]*Coupling[Yd, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
          Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[1]]**Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
           Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
            Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0] + 
         Bar[Coupling[Yu, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$2, Flavor], 
            Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$3, Flavor], 
            Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
          Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
        Bar[Index[d$$3, SU2L[fund]]]}]*Field[H, Scalar, 
       {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
       {Index[d$$3, SU2L[fund]]}, {}]*
      ((Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$4, Flavor], 
           Index[d$$3, Flavor]}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
       (Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$4, Flavor], 
            Index[d$$3, Flavor]}, 0] + 
         Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$4, Flavor], 
            Index[d$$3, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
        Coupling[M, {}, 0]^4)*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
          Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}], 
   "hScalar-lScalar-lFermion-lVector-lScalar" -> 0, 
   "hScalar-lScalar-lVector-lScalar-lScalar" -> 0, 
   "hScalar-lScalar-lVector-lFermion-lScalar" -> 0, 
   "hScalar-lScalar-lVector-lVector-lScalar" -> 0, 
   "hScalar-lScalar-hScalar-lScalar-hScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     ((4*Coupling[A, {}, 0]^6)/(\[Epsilon]*Coupling[M, {}, 0]^8) + 
      ((22*Coupling[A, {}, 0]^6)/3 + 4*Coupling[A, {}, 0]^6*
         Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^8), 
   "hScalar-lScalar-hScalar-lScalar-lVector-lScalar" -> 0, 
   "hScalar-lScalar-lVector-lScalar-lVector-lScalar" -> 0|>, 
 "Matching Conditions" -> {Coupling[gL, {}, 0] -> Coupling[gL, {}, 0], 
   Coupling[gs, {}, 0] -> Coupling[gs, {}, 0], Coupling[gY, {}, 0] -> 
    Coupling[gY, {}, 0], Coupling[Yd, {Index[i1_, Flavor], 
      Index[i2_, Flavor]}, 0] -> 
    Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0] - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^2) + 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (12*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
   Coupling[Ye, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0] - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^2) + 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (12*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
   Coupling[Yu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0] - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^2) + 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (12*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
   Coupling[\[Lambda], {}, 0] -> (7*hbar*Coupling[A, {}, 0]^4)/
      Coupling[M, {}, 0]^4 + (4*hbar*Coupling[A, {}, 0]^4)/
      (\[Epsilon]*Coupling[M, {}, 0]^4) - Coupling[A, {}, 0]^2/
      Coupling[M, {}, 0]^2 + (6*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0])/Coupling[M, {}, 0]^2 + 
     (4*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
      (\[Epsilon]*Coupling[M, {}, 0]^2) - (hbar*Coupling[\[Kappa], {}, 0]^2)/
      (2*\[Epsilon]) + Coupling[\[Lambda], {}, 0] - 
     (7*hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0])/
      Coupling[M, {}, 0]^2 - (6*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Lambda], {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^2) - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0])/
      (2*Coupling[M, {}, 0]^2) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Lambda]\[Phi], {}, 0])/(2*\[Epsilon]*
       Coupling[M, {}, 0]^2) - (92*hbar*Coupling[A, {}, 0]^4*
       Coupling[\[Mu]2, {}, 2])/(3*Coupling[M, {}, 0]^6) - 
     (16*hbar*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2])/
      (\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (31*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[\[Mu]2, {}, 2])/(9*Coupling[M, {}, 0]^4) - 
     (10*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[\[Mu]2, {}, 2])/(3*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (11*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Coupling[\[Mu]2, {}, 2])/Coupling[M, {}, 0]^4 - 
     (6*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Coupling[\[Mu]2, {}, 2])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (80*hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
       Coupling[\[Mu]2, {}, 2])/(3*Coupling[M, {}, 0]^4) + 
     (16*hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
       Coupling[\[Mu]2, {}, 2])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (2*hbar*Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0])/
      Coupling[M, {}, 0]^4 - (hbar*Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 
        0]*Coupling[\[Mu]\[Phi], {}, 0])/Coupling[M, {}, 0]^2 + 
     (5*hbar*Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*
       Coupling[\[Mu]\[Phi], {}, 0])/Coupling[M, {}, 0]^6 + 
     (2*hbar*Coupling[A, {}, 0]^3*Coupling[\[Mu]2, {}, 2]*
       Coupling[\[Mu]\[Phi], {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^6) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2)/
      (2*Coupling[M, {}, 0]^4) + (4*hbar*Coupling[A, {}, 0]^4*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
     (4*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2 - 
     (hbar*Coupling[\[Kappa], {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      2 - (6*hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2 - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^2) - 
     (16*hbar*Coupling[A, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 - 
     (10*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[\[Mu]2, {}, 2]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (3*Coupling[M, {}, 0]^4) - (6*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]2, {}, 2]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
     (16*hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
       Coupling[\[Mu]2, {}, 2]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4 + (2*hbar*Coupling[A, {}, 0]^3*
       Coupling[\[Mu]2, {}, 2]*Coupling[\[Mu]\[Phi], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6, 
   Coupling[\[Mu]2, {}, 2] -> hbar*Coupling[A, {}, 0]^2 + 
     (hbar*Coupling[A, {}, 0]^2)/\[Epsilon] + 
     (hbar*Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 0])/2 + 
     (hbar*Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 0])/(2*\[Epsilon]) + 
     Coupling[\[Mu]2, {}, 2] - (3*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Mu]2, {}, 2])/(2*Coupling[M, {}, 0]^2) - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2])/
      (\[Epsilon]*Coupling[M, {}, 0]^2) + 
     (11*hbar*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2)/
      (3*Coupling[M, {}, 0]^4) + (2*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Mu]2, {}, 2]^2)/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (hbar*Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/2 - 
     (hbar*Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/(2*\[Epsilon]) + 
     hbar*Coupling[A, {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2] + 
     (hbar*Coupling[M, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2 - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^2 + 
     (2*hbar*Coupling[A, {}, 0]^2*Coupling[\[Mu]2, {}, 2]^2*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
     (hbar*Coupling[A, {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/2, 
   Coupling[cllHH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cdB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cdd, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cdG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cdH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) + 
     (3*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (29*hbar*Coupling[A, {}, 0]^4*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(12*Coupling[M, {}, 0]^6) - 
     (hbar*Coupling[A, {}, 0]^4*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (31*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (36*Coupling[M, {}, 0]^4) - (5*hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(6*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (3*hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Kappa], {}, 0])/
      (4*Coupling[M, {}, 0]^4) + (29*hbar*Coupling[A, {}, 0]^2*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[\[Lambda], {}, 0])/(12*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
      (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^3*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/
      (4*Coupling[M, {}, 0]^6) + 
     (3*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^4*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^6) - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
   Coupling[cduq, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cduu, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cdW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[ceB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[ced, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cee, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[ceH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) + 
     (3*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (29*hbar*Coupling[A, {}, 0]^4*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(12*Coupling[M, {}, 0]^6) - 
     (hbar*Coupling[A, {}, 0]^4*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (31*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (36*Coupling[M, {}, 0]^4) - (5*hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(6*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (3*hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Kappa], {}, 0])/
      (4*Coupling[M, {}, 0]^4) + (29*hbar*Coupling[A, {}, 0]^2*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[\[Lambda], {}, 0])/(12*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
      (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^3*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/
      (4*Coupling[M, {}, 0]^6) + 
     (3*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^4*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^6) - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
   Coupling[ceu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[ceW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cG, {}, 0] -> 0, Coupling[cGt, {}, 0] -> 0, 
   Coupling[cH, {}, 0] -> (19*hbar*Coupling[A, {}, 0]^6)/
      Coupling[M, {}, 0]^8 + (10*hbar*Coupling[A, {}, 0]^6)/
      (\[Epsilon]*Coupling[M, {}, 0]^8) + 
     (31*hbar*Coupling[A, {}, 0]^4*Coupling[gL, {}, 0]^2)/
      (18*Coupling[M, {}, 0]^6) + (5*hbar*Coupling[A, {}, 0]^4*
       Coupling[gL, {}, 0]^2)/(3*\[Epsilon]*Coupling[M, {}, 0]^6) + 
     (83*hbar*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0])/
      (4*Coupling[M, {}, 0]^6) + (21*hbar*Coupling[A, {}, 0]^4*
       Coupling[\[Kappa], {}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
      (2*Coupling[M, {}, 0]^4) + (11*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0]^2)/(2*Coupling[M, {}, 0]^4) + 
     (3*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2)/
      (\[Epsilon]*Coupling[M, {}, 0]^4) - (hbar*Coupling[\[Kappa], {}, 0]^3)/
      (12*Coupling[M, {}, 0]^2) - (100*hbar*Coupling[A, {}, 0]^4*
       Coupling[\[Lambda], {}, 0])/(3*Coupling[M, {}, 0]^6) - 
     (20*hbar*Coupling[A, {}, 0]^4*Coupling[\[Lambda], {}, 0])/
      (\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (31*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[\[Lambda], {}, 0])/(18*Coupling[M, {}, 0]^4) - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[\[Lambda], {}, 0])/(3*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (27*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Coupling[\[Lambda], {}, 0])/(2*Coupling[M, {}, 0]^4) - 
     (9*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Coupling[\[Lambda], {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (41*hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda], {}, 0]^2)/
      (3*Coupling[M, {}, 0]^4) + (10*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Lambda], {}, 0]^2)/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (hbar*Coupling[A, {}, 0]^4*Coupling[\[Lambda]\[Phi], {}, 0])/
      (2*Coupling[M, {}, 0]^6) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0])/
      (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0])/
      (4*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (33*hbar*Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0])/
      (4*Coupling[M, {}, 0]^8) - (7*hbar*Coupling[A, {}, 0]^5*
       Coupling[\[Mu]\[Phi], {}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^8) + 
     (Coupling[A, {}, 0]^3*Coupling[\[Mu]\[Phi], {}, 0])/
      (6*Coupling[M, {}, 0]^6) - (5*hbar*Coupling[A, {}, 0]^3*
       Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/
      Coupling[M, {}, 0]^6 - (hbar*Coupling[A, {}, 0]^3*
       Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/
      (\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (hbar*Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]^2*
       Coupling[\[Mu]\[Phi], {}, 0])/(4*Coupling[M, {}, 0]^4) + 
     (13*hbar*Coupling[A, {}, 0]^3*Coupling[\[Lambda], {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0])/(2*Coupling[M, {}, 0]^6) + 
     (3*hbar*Coupling[A, {}, 0]^3*Coupling[\[Lambda], {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^6) + 
     (hbar*Coupling[A, {}, 0]^3*Coupling[\[Lambda]\[Phi], {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0])/(3*Coupling[M, {}, 0]^6) + 
     (hbar*Coupling[A, {}, 0]^3*Coupling[\[Lambda]\[Phi], {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0])/(12*\[Epsilon]*Coupling[M, {}, 0]^6) + 
     (hbar*Coupling[A, {}, 0]^4*Coupling[\[Mu]\[Phi], {}, 0]^2)/
      Coupling[M, {}, 0]^8 + (hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]^2)/
      (2*Coupling[M, {}, 0]^6) - (hbar*Coupling[A, {}, 0]^3*
       Coupling[\[Mu]\[Phi], {}, 0]^3)/(6*Coupling[M, {}, 0]^8) + 
     (10*hbar*Coupling[A, {}, 0]^6*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^8 + (5*hbar*Coupling[A, {}, 0]^4*
       Coupling[gL, {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (3*Coupling[M, {}, 0]^6) + (21*hbar*Coupling[A, {}, 0]^4*
       Coupling[\[Kappa], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (2*Coupling[M, {}, 0]^6) + (3*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4 - (20*hbar*Coupling[A, {}, 0]^4*
       Coupling[\[Lambda], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^6 - (5*hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Coupling[\[Lambda], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(3*Coupling[M, {}, 0]^4) - 
     (9*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Coupling[\[Lambda], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4 + (10*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Lambda], {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^4 - (hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^4) - 
     (7*hbar*Coupling[A, {}, 0]^5*Coupling[\[Mu]\[Phi], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^8) - 
     (hbar*Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      Coupling[M, {}, 0]^6 + (3*hbar*Coupling[A, {}, 0]^3*
       Coupling[\[Lambda], {}, 0]*Coupling[\[Mu]\[Phi], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 + 
     (hbar*Coupling[A, {}, 0]^3*Coupling[\[Lambda]\[Phi], {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (12*Coupling[M, {}, 0]^6), Coupling[cHB, {}, 0] -> 
    (hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/
     (12*Coupling[M, {}, 0]^4), Coupling[cHBox, {}, 0] -> 
    (19*hbar*Coupling[A, {}, 0]^4)/(12*Coupling[M, {}, 0]^6) + 
     (hbar*Coupling[A, {}, 0]^4)/(\[Epsilon]*Coupling[M, {}, 0]^6) - 
     Coupling[A, {}, 0]^2/(2*Coupling[M, {}, 0]^4) - 
     (31*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2)/
      (24*Coupling[M, {}, 0]^4) - (5*hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2)/(4*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (31*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/
      (72*Coupling[M, {}, 0]^4) - (5*hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2)/(12*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (25*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
      (12*Coupling[M, {}, 0]^4) + (2*hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Kappa], {}, 0])/(\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (hbar*Coupling[\[Kappa], {}, 0]^2)/(24*Coupling[M, {}, 0]^2) - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0])/
      (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[\[Lambda]\[Phi], {}, 0])/(2*\[Epsilon]*
       Coupling[M, {}, 0]^4) - (13*hbar*Coupling[A, {}, 0]^3*
       Coupling[\[Mu]\[Phi], {}, 0])/(12*Coupling[M, {}, 0]^6) - 
     (5*hbar*Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]*
       Coupling[\[Mu]\[Phi], {}, 0])/(12*Coupling[M, {}, 0]^4) + 
     (11*hbar*Coupling[A, {}, 0]^2*Coupling[\[Mu]\[Phi], {}, 0]^2)/
      (24*Coupling[M, {}, 0]^6) + (hbar*Coupling[A, {}, 0]^4*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^4) - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(12*Coupling[M, {}, 0]^4) + 
     (2*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
     (hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
   Coupling[cHBt, {}, 0] -> 0, 
   Coupling[cHd, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(8*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(4*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (216*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (36*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (4*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(36*Coupling[M, {}, 0]^4), 
   Coupling[cHD, {}, 0] -> (-31*hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2)/(18*Coupling[M, {}, 0]^4) - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/
      (3*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(3*Coupling[M, {}, 0]^4), 
   Coupling[cHe, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-5*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(8*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(4*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (72*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (12*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (4*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(12*Coupling[M, {}, 0]^4), 
   Coupling[cHG, {}, 0] -> 0, Coupling[cHGt, {}, 0] -> 0, 
   Coupling[cHl1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(8*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (144*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (24*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (8*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(24*Coupling[M, {}, 0]^4), 
   Coupling[cHl3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(8*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (144*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (24*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (8*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(24*Coupling[M, {}, 0]^4), 
   Coupling[cHq1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(8*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (5*hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(8*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (432*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (72*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (8*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (8*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(72*Coupling[M, {}, 0]^4), 
   Coupling[cHq3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(8*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (5*hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0])/(8*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (144*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (24*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (8*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (8*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(24*Coupling[M, {}, 0]^4), 
   Coupling[cHu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(8*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(4*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (108*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (18*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (4*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
       Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(18*Coupling[M, {}, 0]^4), 
   Coupling[cHud, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (2*Coupling[M, {}, 0]^4), Coupling[cHW, {}, 0] -> 
    (hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2)/
     (12*Coupling[M, {}, 0]^4), Coupling[cHWB, {}, 0] -> 
    (hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]*Coupling[gY, {}, 0])/
     (6*Coupling[M, {}, 0]^4), Coupling[cHWt, {}, 0] -> 0, 
   Coupling[cHWtB, {}, 0] -> 0, 
   Coupling[cld, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cle, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 
    -1/12*(hbar*Bar[Coupling[Ye, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
   Coupling[cledq, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    (hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[i3, Flavor]}, 0]]*
      Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
        Index[i2, Flavor]}, 0])/(6*Coupling[M, {}, 0]^4), 
   Coupling[clequ1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/6*(hbar*Coupling[A, {}, 0]^2*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Yu, {Index[i3, Flavor], 
         Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
   Coupling[clequ3, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cll, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[clq1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[clq3, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[clu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqd1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/36*(hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
   Coupling[cqd8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/6*(hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
   Coupling[cqe, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqq1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqq3, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqqq, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqqu, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqu1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/36*(hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
   Coupling[cqu8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/6*(hbar*Bar[Coupling[Yu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
   Coupling[cquqd1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    (hbar*Coupling[A, {}, 0]^2*Coupling[Yd, {Index[i3, Flavor], 
        Index[i4, Flavor]}, 0]*Coupling[Yu, {Index[i1, Flavor], 
        Index[i2, Flavor]}, 0])/(6*Coupling[M, {}, 0]^4), 
   Coupling[cquqd8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cuB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cud1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cud8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cuG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cuH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (5*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) + 
     (3*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (29*hbar*Coupling[A, {}, 0]^4*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(12*Coupling[M, {}, 0]^6) - 
     (hbar*Coupling[A, {}, 0]^4*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(2*\[Epsilon]*Coupling[M, {}, 0]^6) - 
     (31*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (36*Coupling[M, {}, 0]^4) - (5*hbar*Coupling[A, {}, 0]^2*
       Coupling[gL, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(6*\[Epsilon]*Coupling[M, {}, 0]^4) - 
     (3*hbar*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Kappa], {}, 0])/
      (4*Coupling[M, {}, 0]^4) + (29*hbar*Coupling[A, {}, 0]^2*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[\[Lambda], {}, 0])/(12*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
      (2*\[Epsilon]*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^3*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Mu]\[Phi], {}, 0])/
      (4*Coupling[M, {}, 0]^6) + 
     (3*hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[A, {}, 0]^2*Coupling[Yu, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[Yu, {Index[i1, Flavor], 
         Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
      (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^4*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^6) - 
     (5*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(6*Coupling[M, {}, 0]^4) + 
     (hbar*Coupling[A, {}, 0]^2*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
       Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
   Coupling[cuu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cuW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cW, {}, 0] -> 0, Coupling[cWt, {}, 0] -> 0}|>
