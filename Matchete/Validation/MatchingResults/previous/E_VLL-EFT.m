(* Created with the Wolfram Language : www.wolfram.com *)
<|"Model" -> "E_VLL", "Version" -> "0.3.0", 
 "Date" -> DateObject[{2025, 5, 30}, "Day"], "Time (Match)" -> 11.435209, 
 "Time (GreensSimplify)" -> 0.204514, "Time (EOMSimplify)" -> 3.985093, 
 "Time (MapEffectiveCouplings)" -> 15.525172, 
 "UV Lagrangian" -> Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}] + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}] - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/2 - 
   FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2/
    (4*Coupling[gY, {}, 0]^2) - 
   FieldStrength[G, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU3c[adj]]}, {}]^2/(4*Coupling[gs, {}, 0]^2) - 
   FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]^2/(4*Coupling[gL, {}, 0]^2) - 
   Coupling[ME, {}, 0]*Bar[Field[EE, Fermion, {}, {}]]**
     Field[EE, Fermion, {}, {}] + 
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
       Index[d$$1, Flavor]}, {}] + I*Bar[Field[EE, Fermion, {}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]]]**Field[EE, Fermion, {}, 
      {Index[d$$1, Lorentz]}] - Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[EE, Fermion, {}, {}]]**DiracProduct[Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}] - 
   Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$2, Flavor]}, 
      {}] - Coupling[yE, {Index[d$$1, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[EE, Fermion, {}, {}] + 
   I*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
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
  (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz], 
        Index[d$$1, Lorentz]}]]*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$2, Lorentz], 
       Index[d$$2, Lorentz]}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     (Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0] + 
      6*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])*
     (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
      Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz], 
         Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
        {}]))/(6*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
      2)/(8*Coupling[ME, {}, 0]^2) - 
   (2*hbar*FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
      {Index[d$$2, Lorentz]}]*FieldStrength[B, {Index[d$$1, Lorentz], 
       Index[d$$3, Lorentz]}, {}, {Index[d$$3, Lorentz]}])/
    (15*Coupling[ME, {}, 0]^2) - 
   FieldStrength[G, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU3c[adj]]}, {}]^2/(4*Coupling[gs, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}])/(3*Coupling[ME, {}, 0]^2) - 
   FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]^2/(4*Coupling[gL, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
       {Index[d$$1, SU2L[adj]]}, {}]^2)/(24*Coupling[ME, {}, 0]^2) + 
   FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2*
    (-1/4*1/Coupling[gY, {}, 0]^2 + hbar*(-1/3*1/\[Epsilon] - 
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]/3)) - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
    (-1 + hbar*(-1/2*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$1, Flavor]}, 0]) - 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$1, Flavor]}, 0])/\[Epsilon] - 
       Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
         {Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])) + 
   hbar*(-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {}]*FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
        {}, {Index[d$$2, Lorentz]}]) + 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
       {Index[d$$2, Lorentz]}])*
    (((I/6)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (((13*I)/36)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$1, Flavor]}, 0] + 
       (I/6)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$1, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2) + 
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
    (((-1/3*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (((-5*I)/18)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$1, Flavor]}, 0] - 
       (I/3)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$1, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2) + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*
         (-(Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[yE, 
             {Index[d$$1, Flavor]}, 0]) + 
          Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]))/2 + 
       Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2) + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
    ((2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*
         (-2*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$1, Flavor]}, 0] + 
          9*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]))/3 + 
       2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    (-1/2*Coupling[\[Lambda], {}, 0] + 
     hbar*(-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0] - 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*(Bar[Coupling[yE, {Index[d$$2, Flavor]}, 
             0]]*Coupling[yE, {Index[d$$1, Flavor]}, 0] + 
          2*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]))/
        \[Epsilon] - Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*
        (Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[yE, 
           {Index[d$$1, Flavor]}, 0] + 
         2*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])) + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
    ((-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
       Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
              {Index[d$$2, Flavor]}, 0] + 6*Bar[Coupling[Ye, {Index[d$$3, 
                 Flavor], Index[d$$4, Flavor]}, 0]]*Coupling[Ye, 
              {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]) - 
          6*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
           Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]))/3 - 
       2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
        Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    (Coupling[\[Mu]2, {}, 2] + 
     hbar*((-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[ME, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0])/
        \[Epsilon] + Coupling[ME, {}, 0]^2*
        (-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$1, Flavor]}, 0] - 
         2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
           {Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^
             2]))) - Bar[Coupling[Yd, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]]*Bar[Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   I*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] - 
   Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$3, Flavor]}, 0]*(3*Coupling[gY, {}, 0]^2 + 
        2*Coupling[\[Lambda], {}, 0]))/(2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*(Coupling[gY, {}, 0]^2 + 
          2*Coupling[\[Lambda], {}, 0]))/2 + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*(3*Coupling[gY, {}, 0]^2 + 
          2*Coupling[\[Lambda], {}, 0])*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
        2)/Coupling[ME, {}, 0]^2)*
    Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   hbar*(((-1/4*I)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((I/8)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] - 
       (I/4)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
   hbar*(((I/4)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (((3*I)/8)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] + 
       (I/4)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
      {}] - hbar*FieldStrength[B, {Index[d$$1, Lorentz], 
      Index[d$$2, Lorentz]}, {}, {Index[d$$2, Lorentz]}]*
    ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0])/(12*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     ((25*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/72 + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/12)/
      Coupling[ME, {}, 0]^2)*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}] + 
   (I*Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]] + 
     hbar*(((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
       ((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         (I*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
          Coupling[ME, {}, 0]^2)/\[Epsilon] + 
       (I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] + 
       (((-3*I)/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2] - 
         I*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
           {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2))*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
      {Index[d$$1, Lorentz]}] - hbar*CG[gen[SU2L[fund]], 
     {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
      Bar[Index[d$$2, SU2L[fund]]]}]*FieldStrength[W, 
     {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {Index[d$$1, SU2L[adj]]}, 
     {Index[d$$2, Lorentz]}]*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0])/(6*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     ((11*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/36 + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/6)/
      Coupling[ME, {}, 0]^2)*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}] - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     (-(FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
         {}]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
            Index[d$$2, Flavor]}, {}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz], Index[d$$2, Lorentz]], 
          GammaM[Index[d$$3, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$3, Lorentz]}]) + 
      FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
       Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
          {Index[d$$3, Lorentz]}]]**DiracProduct[
         GammaM[Index[d$$3, Lorentz]], GammaM[Index[d$$1, Lorentz], 
          Index[d$$2, Lorentz]], Proj[-1]]**Field[l, Fermion, 
         {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]))/
    (8*Coupling[ME, {}, 0]^2) + 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     (Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
         {Index[d$$1, Lorentz]}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
        Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {Index[d$$2, Lorentz], 
         Index[d$$2, Lorentz]}] - 
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
         {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]**
       DiracProduct[GammaM[Index[d$$2, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$2, Lorentz]}]))/Coupling[ME, {}, 0]^2 + 
   hbar*(-1/4*(Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 
          0]]*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
        Coupling[Yd, {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
        Coupling[yE, {Index[d$$4, Flavor]}, 0])/(\[Epsilon]*
        Coupling[ME, {}, 0]^2) + 
     ((-3*Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[Yd, 
          {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0])/8 - 
       (Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[Yd, 
          {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*
         Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}] + 
   hbar*(-1/4*(Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0])/
       (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((-3*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0])/8 - 
       (Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
         Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
    Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}] + 
   hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
       (3*Coupling[gY, {}, 0]^2 + 2*Coupling[\[Lambda], {}, 0]))/
      (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0]*(Coupling[gY, {}, 0]^2 + 
          2*Coupling[\[Lambda], {}, 0]))/2 + 
       (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0]*(3*Coupling[gY, {}, 0]^2 + 
          2*Coupling[\[Lambda], {}, 0])*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
        2)/Coupling[ME, {}, 0]^2)*
    Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
      {}] + (((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0])/Coupling[ME, {}, 0]^2 + 
     hbar*((((-3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 + 3*Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       ((-1/16*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (11*Coupling[gL, {}, 0]^2 + 9*Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         ((3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 + 3*Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2))*
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
   hbar*(((-1/4*I)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (3*Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
          2*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]) + 
        2*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((-1/8*I)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (11*Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
           4*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
            (-(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
                {Index[d$$2, Flavor]}, 0]) + Bar[Coupling[Ye, 
                {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
              Coupling[Ye, {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0])) + 
         4*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) - 
       (I/4)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (3*Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
           2*Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]) + 
         2*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
    (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
         {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {}]*Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
           Index[d$$2, Flavor]}, {}]]**DiracProduct[
         GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
         {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]) + 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
        {}]) + (((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0])/Coupling[ME, {}, 0]^2 + 
     hbar*(((-1/2*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*(3*Coupling[gY, {}, 0]^2 - 
          Coupling[\[Lambda], {}, 0]))/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((-1/4*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*(2*Coupling[gY, {}, 0]^2 - 
           3*Coupling[\[Lambda], {}, 0]) - 
         (I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*(3*Coupling[gY, {}, 0]^2 - 
           Coupling[\[Lambda], {}, 0])*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
        Coupling[ME, {}, 0]^2))*
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
         Index[d$$1, Flavor]}, {}]) + 
   hbar*(((I/4)*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0] + 
        Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
          2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 
            0])))/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((I/8)*(3*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0] + 
         Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (-2*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
             {Index[d$$3, Flavor]}, 0] + 3*Bar[Coupling[Ye, {Index[d$$3, 
                Flavor], Index[d$$4, Flavor]}, 0]]*Coupling[yE, 
             {Index[d$$3, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
              Index[d$$4, Flavor]}, 0] + 6*Coupling[yE, 
             {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])) + 
       (I/4)*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0] + 
         Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 
             0]))*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      Coupling[ME, {}, 0]^2)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$1, Lorentz]}] - 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {Index[d$$1, Lorentz]}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
        Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}]) + 
   hbar*(((I/4)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
        Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (((5*I)/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
         Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) + 
       (I/4)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
         Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
    (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$1, Lorentz]}] + 
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {Index[d$$1, Lorentz]}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
        Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}]) + 
   hbar*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0])/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (-1/8*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*
         (Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$3, Flavor]}, 0] - 
          3*Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[Ye, {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0])) + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
         Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}] - 
   Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   I*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}] - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
    Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   hbar*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$3, Flavor]}, 0]*
       (Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[Yd, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0] + 
        Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[Yu, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]))/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$3, Flavor]}, 0]*
         (Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[Yd, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0] + 
          Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[Yu, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]))/8 + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$3, Flavor]}, 0]*
         (Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[Yd, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0] + 
          Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[Yu, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0])*
         Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$4, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$2, Flavor]}, {}] + 
   hbar*(-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
         {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0])/
       (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((-3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
          {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0])/8 - 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
          {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
         Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$4, Flavor]}, {}]]**DiracProduct[
      GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
       Index[d$$2, Flavor]}, {}] - 
   Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
      Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
     {Index[d$$2, SU2L[fund]]}, {}]*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   I*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] + 
   hbar*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
        {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0])/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     ((3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
          {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0])/8 + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
          {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]*
         Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}], 
 "On-shell EFT Lagrangian" -> 
  Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}] - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/2 - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/\[Epsilon] + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}] - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/\[Epsilon] - 
   2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
    Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Coupling[yE, {Index[d$$2, Flavor]}, 0]*
    Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] - 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/\[Epsilon] - 
   (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/2 + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/2 + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/\[Epsilon] - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (15*Coupling[ME, {}, 0]^2) - 
   (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (36*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (36*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/Coupling[ME, {}, 0]^2 - 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (3*Coupling[ME, {}, 0]^2) - 
   (2*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(15*Coupling[ME, {}, 0]^2) + 
   (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (18*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (3*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}])/(2*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (2*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^4*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
    (15*Coupling[ME, {}, 0]^2) - 
   (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
    (12*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
    (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
    (36*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
    (6*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}])/(3*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
    Coupling[ME, {}, 0]^2 + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$2, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/Coupling[ME, {}, 0]^2 - 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/Coupling[ME, {}, 0]^2 - 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (15*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (36*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (36*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (6*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (3*Coupling[ME, {}, 0]^2) - 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/Coupling[ME, {}, 0]^2 + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]^2*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
    (3*Coupling[ME, {}, 0]^2) - 
   (hbar*FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
       {}]^2)/(3*\[Epsilon]) - 
   FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^2/
    (4*Coupling[gY, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
      2)/(8*Coupling[ME, {}, 0]^2) - 
   FieldStrength[G, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU3c[adj]]}, {}]^2/(4*Coupling[gs, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}])/(3*Coupling[ME, {}, 0]^2) - 
   FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]^2/(4*Coupling[gL, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$1, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
       {Index[d$$1, SU2L[adj]]}, {}]^2)/(24*Coupling[ME, {}, 0]^2) - 
   hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] - 
   hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
    Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Coupling[yE, {Index[d$$1, Flavor]}, 0]*
    Coupling[yE, {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, 
     {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
     {Index[d$$2, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] - 
   2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
    Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Coupling[yE, {Index[d$$2, Flavor]}, 0]*
    Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] + 
   hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
    Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
    Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2 - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(2*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2 - 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2 + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) - 
   (hbar*FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
       {}]^2*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/3 + 
   hbar*(-2*Bar[Coupling[yE, {Index[u49025, Flavor]}, 0]]*
      Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Coupling[ME, {}, 0]^2*Coupling[yE, {Index[u49025, Flavor]}, 0]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}] - 
     (2*Bar[Coupling[yE, {Index[u49026, Flavor]}, 0]]*
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[ME, {}, 0]^2*Coupling[yE, {Index[u49026, Flavor]}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}])/\[Epsilon] - 
     2*Bar[Coupling[yE, {Index[u49027, Flavor]}, 0]]*
      Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Coupling[ME, {}, 0]^2*Coupling[yE, {Index[u49027, Flavor]}, 0]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]) - 
   Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/4 + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(30*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (13*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[-1]]**
      Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/2 + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[ME, {}, 0]^2) - 
   (((2*I)/45)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[gY, {}, 0]^4*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((13*I)/108)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/18)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((2*I)/45)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((13*I)/108)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/18)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   ((I/18)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/18)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   I*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] - (2*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, {}]*
     Bar[Field[d, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, {}])/
    (135*Coupling[ME, {}, 0]^2) - 
   Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/4 + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/8 + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*\[Epsilon]) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) - 
   (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(30*Coupling[ME, {}, 0]^2) + 
   (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (13*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[ME, {}, 0]^2) - 
   (19*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(16*Coupling[ME, {}, 0]^2) - 
   (9*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(8*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (37*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(24*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/2 + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/4 - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[ME, {}, 0]^2) - 
   (9*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(8*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(4*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/24)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/4)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((2*I)/15)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[gY, {}, 0]^4*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   (((13*I)/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   ((I/24)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/4)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((2*I)/15)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   (((13*I)/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/4)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/4)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   I*Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[e, Fermion, {Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}] - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*FieldStrength[B, 
      {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz], Index[d$$2, Lorentz]], 
       Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz], Index[d$$2, Lorentz]], 
       Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   (4*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$4, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$4, Flavor]}, {}])/
    (45*Coupling[ME, {}, 0]^2) - (2*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$3, Flavor]}, {}]*
     Bar[Field[e, Fermion, {Index[d$$4, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$4, Flavor]}, {}])/
    (15*Coupling[ME, {}, 0]^2) - 
   Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
      {}] + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/4 + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*\[Epsilon]) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/8 + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(4*\[Epsilon]) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(3*Coupling[ME, {}, 0]^2) - 
   (3*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(4*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/2 + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/4 - (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[ME, {}, 0]^2) + 
   I*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz], 
        Index[d$$2, Lorentz]], Proj[1]]**Field[e, Fermion, 
       {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
      {Index[d$$1, SU2L[adj]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz], 
        Index[d$$2, Lorentz]], Proj[1]]**Field[e, Fermion, 
       {Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   (hbar*Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$4, Flavor], Index[d$$1, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$4, Flavor], Index[d$$1, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(3*Coupling[ME, {}, 0]^2) - 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (8*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (25*hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[d, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (216*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[1]]**Field[d, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (36*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (2*hbar*Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Delta[Index[d$$2, Flavor], Index[d$$4, Flavor]]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (45*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (4*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (36*Coupling[ME, {}, 0]^2) - 
   (3*hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (8*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (25*hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (2*hbar*Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Delta[Index[d$$2, Flavor], Index[d$$4, Flavor]]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (15*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (4*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (12*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(30*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (13*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(6*Coupling[ME, {}, 0]^2) + 
   (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[ME, {}, 0]^2) + 
   (37*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(24*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/Coupling[ME, {}, 0]^2 - 
   (19*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(16*Coupling[ME, {}, 0]^2) - 
   (9*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(8*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$4, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(4*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(3*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(2*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(12*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(4*Coupling[ME, {}, 0]^2) - 
   (9*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/(8*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}])/Coupling[ME, {}, 0]^2 + 
   ((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((77*I)/144)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((7*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((9*I)/16)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((9*I)/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((5*I)/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((3*I)/4)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((5*I)/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((77*I)/144)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((7*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((9*I)/16)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((9*I)/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((5*I)/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((3*I)/4)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((5*I)/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((7*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((9*I)/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((3*I)/4)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((7*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((9*I)/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((3*I)/4)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + (((187*I)/144)*hbar*
     Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((17*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((25*I)/144)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/15)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[gY, {}, 0]^4*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((5*I)/72)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((13*I)/72)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((187*I)/144)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((17*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((25*I)/144)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   ((I/15)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((5*I)/72)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((13*I)/72)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((17*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((17*I)/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/24)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (8*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (8*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (11*hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (4*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[ME, {}, 0]^2) + 
   (11*hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (144*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (24*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (25*hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (144*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (24*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Delta[Index[d$$2, Flavor], Index[d$$4, Flavor]]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (30*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (24*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}])/
    (24*Coupling[ME, {}, 0]^2) - 
   Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    4 + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]) - (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    2 + I*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
     DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
     Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$2, Flavor]}, {Index[d$$1, Lorentz]}] + 
   (hbar*Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$4, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     CG[eps[SU2L], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, {Index[d$$1, Flavor]}, 
       {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$2, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (2*hbar*Coupling[gY, {}, 0]^4*Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$3, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$4, Flavor]}, {}])/(45*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (30*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (13*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (6*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[Yd, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[Yd, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[ME, {}, 0]^2) - 
   Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
    CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
    Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
    Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    4 + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (2*\[Epsilon]) - (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    2 + (((5*I)/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((5*I)/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/45)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[gY, {}, 0]^4*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   (((5*I)/72)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((13*I)/216)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/45)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   (((5*I)/72)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (((13*I)/216)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 - 
   ((I/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/12)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   ((I/36)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   (3*hbar*Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(8*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(8*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (11*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(4*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(4*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(12*Coupling[ME, {}, 0]^2) - 
   (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(8*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (11*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(144*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(24*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (25*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(432*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(72*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (hbar*Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Delta[Index[d$$2, Flavor], Index[d$$4, Flavor]]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(45*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(4*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(24*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[gY, {}, 0]^4*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(30*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (13*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yu, {Index[d$$2, Flavor], 
       Index[d$$1, Flavor]}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
      {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
      {}]*Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], 
         Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]]**
      DiracProduct[Proj[1]]**Field[u, Fermion, {Index[d$$1, SU3c[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[Yu, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (12*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     Coupling[Yd, {Index[d$$4, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     CG[eps[SU2L], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]]}]*
     Coupling[Yd, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (3*Coupling[ME, {}, 0]^2) + (2*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$4, Flavor]}, {}])/(135*Coupling[ME, {}, 0]^2) - 
   (hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$3, Flavor]}, {}]*
     Bar[Field[q, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$4, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$2, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$4, Flavor]}, {}])/(270*Coupling[ME, {}, 0]^2) - 
   Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
      Bar[Index[d$$2, SU2L[fund]]]}]*Field[H, Scalar, 
     {Index[d$$2, SU2L[fund]]}, {}]*
    Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
      {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
       Index[d$$1, Flavor]}, {}] + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/4 + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(2*\[Epsilon]) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[\[Mu]2, {}, 2]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/2 + 
   (hbar*Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[gY, {}, 0]^4*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(30*Coupling[ME, {}, 0]^2) + 
   (5*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (13*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(72*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$4, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[yE, {Index[d$$4, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(6*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[yE, {Index[d$$4, Flavor]}, 0]*
     Coupling[Ye, {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/Coupling[ME, {}, 0]^2 + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Coupling[\[Lambda], {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
      {}]*Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[gL, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$2, SU2L[fund]]], 
       Bar[Index[d$$3, SU2L[fund]]]}]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Field[H, Scalar, 
      {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
      {Index[d$$3, SU2L[fund]]}, {}]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(12*Coupling[ME, {}, 0]^2) + 
   (((4*I)/45)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]]*Coupling[gY, {}, 0]^4*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   (((13*I)/54)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 - 
   ((I/9)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (((4*I)/45)*hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   (((13*I)/54)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/9)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   ((I/9)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   ((I/9)*hbar*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
     Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
     Delta[Index[d$$1, Flavor], Index[d$$2, Flavor]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}])/
    Coupling[ME, {}, 0]^2 + 
   I*Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
       {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
     Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, 
      {Index[d$$1, Lorentz]}] + 
   (hbar*Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$5, Flavor]}, 0]*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (8*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
   (25*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (108*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
       {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (18*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
   (4*hbar*Coupling[gY, {}, 0]^4*Delta[Index[d$$1, Flavor], 
      Index[d$$3, Flavor]]*Delta[Index[d$$2, Flavor], Index[d$$4, Flavor]]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (45*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]*
     Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (4*Coupling[ME, {}, 0]^2) - 
   (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[yE, {Index[d$$3, Flavor]}, 0]*Delta[Index[d$$2, Flavor], 
      Index[d$$4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]*
     Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
      Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}])/
    (18*Coupling[ME, {}, 0]^2) + (8*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
      DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[e, Fermion, {Index[d$$3, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, {}])/
    (45*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[yE, {Index[d$$5, Flavor]}, 0]*Coupling[Yu, 
      {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[yE, {Index[d$$5, Flavor]}, 0]]*
     Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0]]*
     CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], 
       Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$5, Flavor]}, 0]*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
        Index[d$$1, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[Proj[-1]]**Field[q, Fermion, 
       {Index[d$$2, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$2, Flavor]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
   (8*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, {}])/
    (135*Coupling[ME, {}, 0]^2) - (4*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$3, Flavor]}, {}]]**DiracProduct[
       GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[q, Fermion, 
       {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
        Index[d$$3, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, {}])/
    (135*Coupling[ME, {}, 0]^2) - (8*hbar*Coupling[gY, {}, 0]^4*
     Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, {}]*
     Bar[Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, 
        {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
      Field[u, Fermion, {Index[d$$2, SU3c[fund]], Index[d$$4, Flavor]}, {}])/
    (135*Coupling[ME, {}, 0]^2), "SuperTraces" -> 
  <|"hFermion-lScalar" -> 
    -(hbar*FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
        {Index[d$$2, Lorentz]}]*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0])/(12*\[Epsilon]*
          Coupling[ME, {}, 0]^2) + 
        ((25*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
             {Index[d$$2, Flavor]}, 0])/72 + 
          (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
             {Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
           12)/Coupling[ME, {}, 0]^2)*
       Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
          {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
         {}]) + hbar*(((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
       ((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/\[Epsilon] + 
       (I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$1, Lorentz]}] - hbar*CG[gen[SU2L[fund]], 
       {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
        Bar[Index[d$$2, SU2L[fund]]]}]*FieldStrength[W, 
       {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
       {Index[d$$1, SU2L[adj]]}, {Index[d$$2, Lorentz]}]*
      ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0])/(6*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((11*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$2, Flavor]}, 0])/36 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
          6)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
        {}] - (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*
       (-(FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
           {}]*Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
              Index[d$$2, Flavor]}, {}]]**DiracProduct[
            GammaM[Index[d$$1, Lorentz], Index[d$$2, Lorentz]], 
            GammaM[Index[d$$3, Lorentz]], Proj[-1]]**Field[l, Fermion, 
            {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
            {Index[d$$3, Lorentz]}]) + FieldStrength[B, 
          {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
         Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, 
              Flavor]}, {Index[d$$3, Lorentz]}]]**DiracProduct[
           GammaM[Index[d$$3, Lorentz]], GammaM[Index[d$$1, Lorentz], 
            Index[d$$2, Lorentz]], Proj[-1]]**Field[l, Fermion, 
           {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]))/
      (8*Coupling[ME, {}, 0]^2) + 
     ((I/6)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*
       (Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$2, Lorentz], Index[d$$2, Lorentz]}] - 
        Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]**
         DiracProduct[GammaM[Index[d$$2, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$2, Lorentz]}]))/Coupling[ME, {}, 0]^2, 
   "hFermion-lFermion" -> (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz], 
          Index[d$$1, Lorentz]}]]*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$2, Lorentz], 
         Index[d$$2, Lorentz]}])/(3*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*
       FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]^
        2)/(8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       CG[gen[SU2L[fund]], {Index[d$$1, SU2L[adj]], Index[d$$1, SU2L[fund]], 
         Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[yE, {Index[d$$1, Flavor]}, 
        0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
       FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, {}]*
       FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
        {Index[d$$1, SU2L[adj]]}, {}])/(3*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*
       FieldStrength[W, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, 
         {Index[d$$1, SU2L[adj]]}, {}]^2)/(24*Coupling[ME, {}, 0]^2) - 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]*
      (-1/2*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$1, Flavor]}, 0]) - 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$1, Flavor]}, 0])/\[Epsilon] - 
       Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
         {Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      ((-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[ME, {}, 0]^2*
         Coupling[yE, {Index[d$$1, Flavor]}, 0])/\[Epsilon] + 
       Coupling[ME, {}, 0]^2*(-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$1, Flavor]}, 0] - 
         2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
           {Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])) + 
     hbar*(-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*FieldStrength[B, 
          {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
          {Index[d$$2, Lorentz]}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        FieldStrength[B, {Index[d$$1, Lorentz], Index[d$$2, Lorentz]}, {}, 
         {Index[d$$2, Lorentz]}])*
      (((I/6)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$1, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       (((13*I)/36)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$1, Flavor]}, 0] + 
         (I/6)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2) + 
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
      (((-1/3*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$1, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       (((-5*I)/18)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$1, Flavor]}, 0] - 
         (I/3)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2), 
   "hFermion-lVector" -> 
    hbar*((((-7*I)/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((-7*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         ((7*I)/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*((((-7*I)/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 
           0]]*Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((-7*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         ((7*I)/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lScalar-lScalar" -> 
    hbar*(((-I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((-3*I)/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2] - 
         I*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
           {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
        {Index[d$$1, Lorentz]}] + 
     hbar*(((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0] + 
         (I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
     hbar*(((I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0] + 
         (I/2)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lScalar-lFermion" -> 0, "hFermion-lScalar-lVector" -> 
    hbar*((((5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((11*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         ((5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*(((I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lFermion-lScalar" -> 0, "hFermion-lFermion-lFermion" -> 0, 
   "hFermion-lFermion-lVector" -> 
    hbar*((I*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         I*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*((((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         ((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lVector-lScalar" -> 
    hbar*((((5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((11*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         ((5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*(((I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lVector-lFermion" -> 
    hbar*((I*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         I*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*((((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         ((3*I)/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lScalar-lScalar-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
        (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
          2 + (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/2)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
        (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
             Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/2 + 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
             Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/2)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], "hFermion-lScalar-hFermion-lScalar" -> 
    -1/8*(hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$4, Flavor]}, 0]*
       Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
          {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
       Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
          {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
        Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
         {}])/Coupling[ME, {}, 0]^2, "hFermion-lScalar-hFermion-lFermion" -> 
    ((I/2)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$3, Flavor]}, 0]*
       (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
            {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
           {Index[d$$1, SU2L[fund]]}, {}]*
          Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
               Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
            Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
             Index[d$$1, Flavor]}, {}]) + 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {}]))/Coupling[ME, {}, 0]^2 - 
     ((I/4)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$3, Flavor]}, 0]*
       (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {Index[d$$1, Lorentz]}]]**DiracProduct[
           GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
           {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]))/
      Coupling[ME, {}, 0]^2, "hFermion-lScalar-lFermion-lScalar" -> 
    hbar*(-1/4*(Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 
            0]]*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
          Coupling[Yd, {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
          Coupling[yE, {Index[d$$4, Flavor]}, 0])/(\[Epsilon]*
          Coupling[ME, {}, 0]^2) + 
       ((-3*Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[Yd, 
            {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
           Coupling[yE, {Index[d$$4, Flavor]}, 0])/8 - 
         (Bar[Coupling[Yd, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[Yd, 
            {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
           Coupling[yE, {Index[d$$4, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$3, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[d, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, Flavor]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}] + hbar*(-1/4*(Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((-3*Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0])/8 - 
         (Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$5, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$5, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$3, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[e, Fermion, {Index[d$$1, Flavor]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$4, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$2, Flavor]}, 
        {}] + hbar*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0])/
        (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}]*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$4, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, 
        {}] + hbar*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*
         (Bar[Coupling[Yd, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[Yd, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0] + 
          Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[Yu, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]))/
        (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$3, Flavor]}, 0]*(Bar[Coupling[Yd, {Index[d$$2, Flavor], 
                Index[d$$5, Flavor]}, 0]]*Coupling[Yd, {Index[d$$4, Flavor], 
               Index[d$$5, Flavor]}, 0] + Bar[Coupling[Yu, {Index[d$$2, 
                 Flavor], Index[d$$5, Flavor]}, 0]]*Coupling[Yu, 
              {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]))/8 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$3, Flavor]}, 0]*(Bar[Coupling[Yd, {Index[d$$2, Flavor], 
                Index[d$$5, Flavor]}, 0]]*Coupling[Yd, {Index[d$$4, Flavor], 
               Index[d$$5, Flavor]}, 0] + Bar[Coupling[Yu, {Index[d$$2, 
                 Flavor], Index[d$$5, Flavor]}, 0]]*Coupling[Yu, 
              {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0])*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$4, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$1, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}] + 
     hbar*(-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
           {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((-3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
            {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0])/8 - 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Yu, {Index[d$$2, Flavor], Index[d$$5, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
            {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
      Bar[Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, 
           SU2L[fund]], Index[d$$4, Flavor]}, {}]]**
       DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[q, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, SU2L[fund]], 
         Index[d$$2, Flavor]}, {}] + 
     hbar*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
          {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0])/
        (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
            {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Yu, {Index[d$$5, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Yu, 
            {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$3, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
       Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*
      Bar[Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$4, Flavor]}, 
         {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[1]]**
       Field[u, Fermion, {Index[d$$1, SU3c[fund]], Index[d$$2, Flavor]}, {}], 
   "hFermion-lScalar-lFermion-lFermion" -> 
    hbar*(((-1/4*I)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
          Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((-1/4*I)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) - 
         (I/4)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}]) + hbar*(((I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
          Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/16)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) + 
         (I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]) + 
     hbar*(((I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
          Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((5*I)/16)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) + 
         (I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lScalar-lVector-lScalar" -> 
    hbar*(((I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       (((3*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*(((I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}]) + hbar*(((I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       (((3*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
     hbar*(((-1/4*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((-3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         (I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lScalar-lVector-lFermion" -> 
    hbar*((((-5*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       (((-11*I)/32)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] - ((5*I)/16)*
          Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*((((-5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 
           0]]*Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((-11*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         ((5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}]) - ((I/16)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*
       (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
            Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {Index[d$$1, Lorentz]}]]**DiracProduct[
           GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
           {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]))/
      Coupling[ME, {}, 0]^2 + 
     ((I/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
       (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {Index[d$$1, Lorentz]}]]**DiracProduct[
           GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
           {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]))/
      Coupling[ME, {}, 0]^2, "hFermion-lFermion-lScalar-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
        (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
          2 + (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/2)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
          {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
           Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
        (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
             Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/2 + 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[yE, 
            {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[d$$3, Flavor], 
             Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/2)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], "hFermion-lFermion-lScalar-lFermion" -> 
    hbar*(((-1/4*I)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((I/8)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] - 
         (I/4)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
     hbar*(((I/4)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/8)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$4, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {Index[d$$1, Flavor]}, {}]), "hFermion-lFermion-lScalar-lVector" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
           2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
            2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], "hFermion-lFermion-hFermion-lFermion" -> 
    -1/2*(hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {Index[d$$1, Lorentz]}]]*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, 
         {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}])/
       Coupling[ME, {}, 0]^2 + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0]*
       (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
           {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz], 
           Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {}]))/(6*Coupling[ME, {}, 0]^2) - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
       Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
         {Index[d$$1, Lorentz]}]]*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Field[H, Scalar, 
        {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, 
        {Index[d$$2, SU2L[fund]]}, {Index[d$$1, Lorentz]}])/
      (3*Coupling[ME, {}, 0]^2) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*Coupling[yE, 
           {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
         \[Epsilon]) - Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
         {Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]), 
   "hFermion-lFermion-lFermion-lScalar" -> 
    hbar*(((-1/4*I)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
          Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((-1/4*I)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) - 
         (I/4)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}]) + hbar*(((I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
          Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((3*I)/16)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) + 
         (I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] + 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]) + 
     hbar*(((I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
          Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]))/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((5*I)/16)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]) + 
         (I/8)*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$3, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$2, Flavor], Index[d$$4, Flavor]}, 0] - 
           Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
            Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]]*
            Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
             {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0])*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lFermion-lFermion-lFermion" -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
           {Index[d$$1, Lorentz], Index[d$$1, Lorentz]}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz], 
           Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
          {}]))/Coupling[ME, {}, 0]^2 + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0] - 
       (2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/\[Epsilon] - 
       2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
        Log[\[Mu]bar2/Coupling[ME, {}, 0]^2]) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
       {}]*((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/2 + 
         Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2) + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
        {Index[d$$1, Lorentz]}]]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, 
       {Index[d$$1, Lorentz]}]*
      ((2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (3*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0] + 
         2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2), 
   "hFermion-lFermion-lFermion-lVector" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0] + 
         2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((2*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] + 
         2*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], "hFermion-lFermion-lVector-lScalar" -> 
    hbar*((((-5*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       (((-11*I)/32)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] - ((5*I)/16)*
          Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*((((-5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 
           0]]*Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((-11*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         ((5*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}]) - ((I/16)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*
       (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
            Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {Index[d$$1, Lorentz]}]]**DiracProduct[
           GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
           {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]))/
      Coupling[ME, {}, 0]^2 + 
     ((I/8)*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
       (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {Index[d$$1, Lorentz]}] - 
        Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
         Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {Index[d$$1, Lorentz]}]]**DiracProduct[
           GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
           {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]))/
      Coupling[ME, {}, 0]^2, "hFermion-lFermion-lVector-lFermion" -> 
    hbar*(((I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       (((-3*I)/16)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
          {}]) + hbar*(((I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (((-3*I)/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (-(Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
           {Index[d$$1, Lorentz]}]]*Field[H, Scalar, 
          {Index[d$$1, SU2L[fund]]}, {}]*
         Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, 
              Flavor]}, {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], 
           Proj[-1]]**Field[l, Fermion, {Index[d$$2, SU2L[fund]], 
            Index[d$$1, Flavor]}, {}]) + 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {Index[d$$1, Lorentz]}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {}]) + hbar*(((-1/8*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$2, Flavor]}, 0])/(\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       ((-1/16*I)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0] - 
         (I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          (Coupling[gL, {}, 0]^2 - Coupling[gY, {}, 0]^2)*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
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
     hbar*(((I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((I/8)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0] + 
         (I/4)*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {}]]**DiracProduct[GammaM[Index[d$$1, Lorentz]], Proj[-1]]**
         Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, 
          {Index[d$$1, Lorentz]}] - 
       Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
        Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
        Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
           {Index[d$$1, Lorentz]}]]**DiracProduct[
          GammaM[Index[d$$1, Lorentz]], Proj[-1]]**Field[l, Fermion, 
          {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]), 
   "hFermion-lVector-lScalar-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
           2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
            2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], "hFermion-lVector-lFermion-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0] + 
         2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((2*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0] + 
         2*Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lScalar-hFermion-lFermion-lFermion" -> 0, 
   "hFermion-lScalar-lFermion-hFermion-lFermion" -> 0, 
   "hFermion-lScalar-lFermion-lFermion-lFermion" -> 0, 
   "hFermion-lScalar-lVector-lScalar-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$3, Flavor]}, 0])/(8*\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*(Coupling[gL, {}, 0]^2 + 
          Coupling[gY, {}, 0]^2)*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
         Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (8*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*(Coupling[gL, {}, 0]^2 + 
            Coupling[gY, {}, 0]^2)*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
           Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lScalar-lVector-lFermion-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
           2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
            2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lFermion-lScalar-hFermion-lFermion" -> 0, 
   "hFermion-lFermion-lScalar-lFermion-lFermion" -> 0, 
   "hFermion-lFermion-lScalar-lVector-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$3, Flavor]}, 0])/(8*\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*(Coupling[gL, {}, 0]^2 + 
          Coupling[gY, {}, 0]^2)*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
         Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (8*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*(Coupling[gL, {}, 0]^2 + 
            Coupling[gY, {}, 0]^2)*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
           Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           (Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lFermion-lScalar-lVector-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$3, Flavor]}, 0])/(8*\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*(-Coupling[gL, {}, 0]^2 + 
          Coupling[gY, {}, 0]^2)*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
         Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (8*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lFermion-hFermion-lFermion-lFermion" -> 0, 
   "hFermion-lFermion-lFermion-lScalar-lFermion" -> 0, 
   "hFermion-lFermion-lFermion-lFermion-lScalar" -> 0, 
   "hFermion-lFermion-lFermion-lFermion-lFermion" -> 0, 
   "hFermion-lFermion-lFermion-lVector-lScalar" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
           2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       (-1/4*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
            2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]) - 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lFermion-lFermion-lVector-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2)) + 
       (-1/2*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]) - 
         Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2)) + 
       (-1/2*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
            2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]) - 
         Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lFermion-lVector-lScalar-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
         (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
         Coupling[yE, {Index[d$$3, Flavor]}, 0])/(8*\[Epsilon]*
         Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$3, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*(-Coupling[gL, {}, 0]^2 + 
          Coupling[gY, {}, 0]^2)*Coupling[yE, {Index[d$$2, Flavor]}, 0]*
         Coupling[Ye, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
        (8*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
       ((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/8 + 
         (Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
           (-Coupling[gL, {}, 0]^2 + Coupling[gY, {}, 0]^2)*
           Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
           Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/8)/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lFermion-lVector-lFermion-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-((Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2)) + 
       (-1/2*(Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
           Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
           Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]) - 
         Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
          Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]]*
          Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[e, Fermion, {Index[d$$2, Flavor]}, {}]]**
       DiracProduct[Proj[-1]]**Field[l, Fermion, {Index[d$$1, SU2L[fund]], 
         Index[d$$1, Flavor]}, {}] + 
     hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      (-((Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
         (\[Epsilon]*Coupling[ME, {}, 0]^2)) + 
       (-1/2*(Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^
            2*Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
            {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]) - 
         Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
          Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
           {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
          Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)*
      Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, 
         {}]]**DiracProduct[Proj[1]]**Field[e, Fermion, 
        {Index[d$$1, Flavor]}, {}], 
   "hFermion-lFermion-hFermion-lFermion-hFermion-lFermion" -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
      Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
      Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
       {Index[d$$2, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0]*
      Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/
     (3*Coupling[ME, {}, 0]^2), 
   "hFermion-lFermion-hFermion-lFermion-lFermion-lFermion" -> 
    (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
      Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
      Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
      Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
       {Index[d$$3, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
        Index[d$$4, Flavor]}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, 
       {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
      Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}])/Coupling[ME, {}, 0]^2, 
   "hFermion-lFermion-lFermion-hFermion-lFermion-lFermion" -> 0, 
   "hFermion-lFermion-lFermion-lFermion-lFermion-lFermion" -> 
    hbar*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*
     Bar[Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]]*
     Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*
     Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]*
     ((-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
        Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
         {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
        Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0])/
       (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
      (-2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
         Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0] - 
        2*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
         Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
         Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
          {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
         Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
         Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2)|>, 
 "Matching Conditions" -> {Coupling[gL, {}, 0] -> -Coupling[gL, {}, 0], 
   Coupling[gs, {}, 0] -> Coupling[gs, {}, 0], Coupling[gY, {}, 0] -> 
    -Coupling[gY, {}, 0] + (2*hbar*Coupling[gY, {}, 0]^3)/(3*\[Epsilon]) + 
     (2*hbar*Coupling[gY, {}, 0]^3*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/3, 
   Coupling[Yd, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0] - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/4 - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/(2*\[Epsilon]) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (3*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/2, 
   Coupling[Ye, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0])/8 - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0])/(4*\[Epsilon]) + 
     Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0] - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0])/4 - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0])/(2*\[Epsilon]) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (4*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (3*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/4 - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/2 + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(2*Coupling[ME, {}, 0]^2), 
   Coupling[Yu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0] - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0])/4 - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0])/(2*\[Epsilon]) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      (3*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/2, Coupling[\[Lambda], {}, 0] -> 
    (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0])/\[Epsilon] + 
     4*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
      Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
       {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0] + 
     (4*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/\[Epsilon] + 
     Coupling[\[Lambda], {}, 0] - 
     hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0] - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
      \[Epsilon] - (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2])/(9*Coupling[ME, {}, 0]^2) - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2])/(3*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      Coupling[ME, {}, 0]^2 + 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2])/Coupling[ME, {}, 0]^2 + 
     (4*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2])/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (4*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
       Coupling[\[Mu]2, {}, 2])/(3*Coupling[ME, {}, 0]^2) + 
     2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
       {Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] + 
     4*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
      Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
       {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
      Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] - 
     2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
      Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (3*Coupling[ME, {}, 0]^2) + 
     (4*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[\[Mu]2, {}, 2]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      Coupling[ME, {}, 0]^2, Coupling[\[Mu]2, {}, 2] -> 
    -2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[ME, {}, 0]^2*
      Coupling[yE, {Index[d$$1, Flavor]}, 0] - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[ME, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0])/
      \[Epsilon] + Coupling[\[Mu]2, {}, 2] - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/2 - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2])/
      \[Epsilon] + (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]^2)/
      (3*Coupling[ME, {}, 0]^2) - 
     2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[ME, {}, 0]^2*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*
      Log[\[Mu]bar2/Coupling[ME, {}, 0]^2] - 
     hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Mu]2, {}, 2]*
      Log[\[Mu]bar2/Coupling[ME, {}, 0]^2], 
   Coupling[cllHH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cdB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cdd, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> (-2*hbar*Coupling[gY, {}, 0]^4*
      Delta[Index[i1, Flavor], Index[i2, Flavor]]*Delta[Index[i3, Flavor], 
       Index[i4, Flavor]])/(135*Coupling[ME, {}, 0]^2), 
   Coupling[cdG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cdH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[Yd, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0]*Coupling[yE, {Index[d$$1, Flavor]}, 0])/
      (36*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/(6*\[Epsilon]*
       Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0])/(2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
      (3*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2, 
   Coupling[cduq, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cduu, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cdW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[ceB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]*
      Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
       {Index[d$$1, Flavor], Index[i2, Flavor]}, 0])/
     (12*Coupling[ME, {}, 0]^2), 
   Coupling[ced, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> (-4*hbar*Coupling[gY, {}, 0]^4*
      Delta[Index[i1, Flavor], Index[i2, Flavor]]*Delta[Index[i3, Flavor], 
       Index[i4, Flavor]])/(45*Coupling[ME, {}, 0]^2), 
   Coupling[cee, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 
    -1/15*(hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
         Index[i4, Flavor]]*Delta[Index[i2, Flavor], Index[i3, Flavor]])/
       Coupling[ME, {}, 0]^2 - (hbar*Coupling[gY, {}, 0]^4*
       Delta[Index[i1, Flavor], Index[i2, Flavor]]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(15*Coupling[ME, {}, 0]^2), 
   Coupling[ceH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[yE, 
        {Index[i1, Flavor]}, 0]*Coupling[Ye, {Index[d$$1, Flavor], 
         Index[i2, Flavor]}, 0])/(2*Coupling[ME, {}, 0]^2) - 
     (19*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[i1, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
         Index[i2, Flavor]}, 0])/(16*Coupling[ME, {}, 0]^2) - 
     (9*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[i1, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
         Index[i2, Flavor]}, 0])/(8*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[Ye, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0])/
      (4*Coupling[ME, {}, 0]^2) + 
     (37*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[Ye, {Index[i1, Flavor], Index[d$$3, Flavor]}, 0])/
      (24*Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[Ye, {Index[i1, Flavor], Index[d$$3, Flavor]}, 0])/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (36*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (6*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0]*Coupling[Ye, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 
        0])/(2*Coupling[ME, {}, 0]^2) + 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 
        0])/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 
        0])/(3*Coupling[ME, {}, 0]^2) - 
     (9*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[i1, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
         Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (8*Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[Ye, {Index[i1, Flavor], Index[d$$3, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(4*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2 + 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 
        0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2, 
   Coupling[ceu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> (8*hbar*Coupling[gY, {}, 0]^4*
      Delta[Index[i1, Flavor], Index[i2, Flavor]]*Delta[Index[i3, Flavor], 
       Index[i4, Flavor]])/(45*Coupling[ME, {}, 0]^2), 
   Coupling[ceW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]*
      Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
       {Index[d$$1, Flavor], Index[i2, Flavor]}, 0])/
     (24*Coupling[ME, {}, 0]^2), Coupling[cG, {}, 0] -> 0, 
   Coupling[cGt, {}, 0] -> 0, Coupling[cH, {}, 0] -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0]*Coupling[yE, {Index[d$$3, Flavor]}, 0])/
      (3*Coupling[ME, {}, 0]^2) + 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$3, Flavor]}, 0]*Coupling[Ye, {Index[d$$2, Flavor], 
         Index[d$$4, Flavor]}, 0])/Coupling[ME, {}, 0]^2 - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
       Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0])/
      Coupling[ME, {}, 0]^2 - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
       Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Coupling[\[Lambda], {}, 0])/(18*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
      (3*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
      Coupling[ME, {}, 0]^2 + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[\[Lambda], {}, 0])/Coupling[ME, {}, 0]^2 + 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[\[Lambda], {}, 0])/(\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]^2)/
      (3*Coupling[ME, {}, 0]^2) - 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$4, Flavor], Index[d$$5, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$5, Flavor]}, 0]*
       Coupling[Ye, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2 - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[\[Lambda], {}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(3*Coupling[ME, {}, 0]^2) + 
     (2*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[\[Lambda], {}, 0]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      Coupling[ME, {}, 0]^2, Coupling[cHB, {}, 0] -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
      Coupling[yE, {Index[d$$1, Flavor]}, 0])/(8*Coupling[ME, {}, 0]^2), 
   Coupling[cHBox, {}, 0] -> -1/30*(hbar*Coupling[gY, {}, 0]^4)/
       Coupling[ME, {}, 0]^2 - 
     (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0])/
      (24*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/(4*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0])/
      (72*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/(12*\[Epsilon]*
       Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0])/(3*Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(4*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(12*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2, 
   Coupling[cHBt, {}, 0] -> 0, 
   Coupling[cHd, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (2*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(45*Coupling[ME, {}, 0]^2) - 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (108*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(18*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (18*Coupling[ME, {}, 0]^2), Coupling[cHD, {}, 0] -> 
    (-2*hbar*Coupling[gY, {}, 0]^4)/(15*Coupling[ME, {}, 0]^2) + 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0])/
      (18*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/(3*\[Epsilon]*
       Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0])/(2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(3*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2, 
   Coupling[cHe, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0])/
      (24*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0])/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (2*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(15*Coupling[ME, {}, 0]^2) - 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (36*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(6*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(4*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (6*Coupling[ME, {}, 0]^2), Coupling[cHG, {}, 0] -> 0, 
   Coupling[cHGt, {}, 0] -> 0, 
   Coupling[cHl1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    -1/4*(Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
        Coupling[yE, {Index[i1, Flavor]}, 0])/Coupling[ME, {}, 0]^2 - 
     (33*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(32*Coupling[ME, {}, 0]^2) - 
     (9*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(16*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     (31*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(288*Coupling[ME, {}, 0]^2) + 
     (25*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(48*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[yE, 
        {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[i1, Flavor]}, 0])/
      (16*Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[yE, 
        {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[i1, Flavor]}, 0])/
      (8*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
      (2*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
      (2*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[d$$2, Flavor]}, 0])/
      (2*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[d$$2, Flavor]}, 0])/
      (2*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (15*Coupling[ME, {}, 0]^2) - 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (72*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (9*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(16*Coupling[ME, {}, 0]^2) + 
     (25*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(48*Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[yE, 
        {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[i2, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(2*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[d$$2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(2*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (12*Coupling[ME, {}, 0]^2), 
   Coupling[cHl3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    -1/4*(Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
        Coupling[yE, {Index[i1, Flavor]}, 0])/Coupling[ME, {}, 0]^2 + 
     (77*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(288*Coupling[ME, {}, 0]^2) + 
     (7*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(48*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     (9*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(32*Coupling[ME, {}, 0]^2) + 
     (9*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(16*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[yE, 
        {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[i1, Flavor]}, 0])/
      (16*Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[yE, 
        {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[i1, Flavor]}, 0])/
      (8*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (72*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (7*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(48*Coupling[ME, {}, 0]^2) + 
     (9*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(16*Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[yE, 
        {Index[d$$1, Flavor]}, 0]*Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (12*Coupling[ME, {}, 0]^2), 
   Coupling[cHq1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    -1/45*(hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
         Index[i2, Flavor]])/Coupling[ME, {}, 0]^2 + 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (216*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(36*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (36*Coupling[ME, {}, 0]^2), 
   Coupling[cHq3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (72*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (12*Coupling[ME, {}, 0]^2), 
   Coupling[cHu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-4*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(45*Coupling[ME, {}, 0]^2) + 
     (13*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gY, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Delta[Index[i1, Flavor], Index[i2, Flavor]])/
      (54*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(9*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (9*Coupling[ME, {}, 0]^2), 
   Coupling[cHud, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cHW, {}, 0] -> (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0])/
     (24*Coupling[ME, {}, 0]^2), Coupling[cHWB, {}, 0] -> 
    -1/6*(hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]*Coupling[gY, {}, 0]*Coupling[yE, 
        {Index[d$$1, Flavor]}, 0])/Coupling[ME, {}, 0]^2, 
   Coupling[cHWt, {}, 0] -> 0, Coupling[cHWtB, {}, 0] -> 0, 
   Coupling[cld, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 
    (-3*hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(4*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     (25*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(216*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(36*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (2*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      (45*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[Yd, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(4*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (36*Coupling[ME, {}, 0]^2), 
   Coupling[cle, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 
    (-3*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0])/
      (8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0])/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i4, Flavor]}, 0])/
      (6*Coupling[ME, {}, 0]^2) + 
     (25*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(72*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (2*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      (15*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(4*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (12*Coupling[ME, {}, 0]^2), 
   Coupling[cledq, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    (hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[i3, Flavor]}, 0]]*
      Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
       {Index[i1, Flavor], Index[i2, Flavor]}, 0])/(3*Coupling[ME, {}, 0]^2), 
   Coupling[clequ1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/3*(hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Coupling[Yu, {Index[i3, Flavor], Index[i4, Flavor]}, 0])/
      Coupling[ME, {}, 0]^2, Coupling[clequ3, {Index[i1_, Flavor], 
      Index[i2_, Flavor], Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cll, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 
    -1/8*(hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
        Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[yE, 
         {Index[i1, Flavor]}, 0]*Coupling[yE, {Index[i3, Flavor]}, 0])/
       Coupling[ME, {}, 0]^2 + 
     (3*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[d$$1, Flavor]}, 0])/
      (8*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[d$$1, Flavor]}, 0])/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (11*hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(144*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(24*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (25*hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(144*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]])/(24*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (11*hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i2, Flavor], 
        Index[i3, Flavor]])/(72*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i2, Flavor], 
        Index[i3, Flavor]])/(12*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      (30*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Coupling[Ye, 
        {Index[i1, Flavor], Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(4*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (24*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (24*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i4, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i2, Flavor], 
        Index[i3, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (12*Coupling[ME, {}, 0]^2), 
   Coupling[clq1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    (3*hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(16*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(8*\[Epsilon]*
       Coupling[ME, {}, 0]^2) - 
     (3*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0])/
      (16*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0])/
      (8*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (25*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(432*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(72*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
       Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      (45*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (72*Coupling[ME, {}, 0]^2), 
   Coupling[clq3, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    (3*hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(16*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0])/(8*\[Epsilon]*
       Coupling[ME, {}, 0]^2) + 
     (3*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0])/
      (16*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0])/
      (8*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (11*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(144*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(24*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[Yd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[Yd, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(8*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(8*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (24*Coupling[ME, {}, 0]^2), 
   Coupling[clu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> 
    (3*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0])/
      (8*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0])/
      (4*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (25*hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(108*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]])/(18*\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (4*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      (45*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Coupling[Yu, 
        {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(4*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
       Coupling[yE, {Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], 
        Index[i4, Flavor]]*Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/
      (18*Coupling[ME, {}, 0]^2), 
   Coupling[cqd1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/18*(hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
        Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Coupling[Yd, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
        Coupling[yE, {Index[d$$1, Flavor]}, 0])/Coupling[ME, {}, 0]^2 + 
     (2*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      (135*Coupling[ME, {}, 0]^2), 
   Coupling[cqd8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/3*(hbar*Bar[Coupling[Yd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[Yd, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0])/Coupling[ME, {}, 0]^2, 
   Coupling[cqe, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> (2*hbar*Coupling[gY, {}, 0]^4*
      Delta[Index[i1, Flavor], Index[i2, Flavor]]*Delta[Index[i3, Flavor], 
       Index[i4, Flavor]])/(45*Coupling[ME, {}, 0]^2), 
   Coupling[cqq1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/270*(hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      Coupling[ME, {}, 0]^2, 
   Coupling[cqq3, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqqq, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqqu, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cqu1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/18*(hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
        Bar[Coupling[Yu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
        Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
         {Index[i1, Flavor], Index[i4, Flavor]}, 0])/Coupling[ME, {}, 0]^2 - 
     (4*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], 
        Index[i2, Flavor]]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
      (135*Coupling[ME, {}, 0]^2), 
   Coupling[cqu8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    -1/3*(hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Yu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i4, Flavor]}, 0])/Coupling[ME, {}, 0]^2, 
   Coupling[cquqd1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
      Coupling[Yd, {Index[i3, Flavor], Index[i4, Flavor]}, 0]*
      Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
       {Index[i1, Flavor], Index[i2, Flavor]}, 0])/(3*Coupling[ME, {}, 0]^2), 
   Coupling[cquqd8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cuB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cud1, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
    (8*hbar*Coupling[gY, {}, 0]^4*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
      Delta[Index[i3, Flavor], Index[i4, Flavor]])/
     (135*Coupling[ME, {}, 0]^2), 
   Coupling[cud8, {Index[i1_, Flavor], Index[i2_, Flavor], 
      Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
   Coupling[cuG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cuH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
    (-5*hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[gL, {}, 0]^2*Coupling[yE, {Index[d$$1, Flavor]}, 0]*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (36*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (6*\[Epsilon]*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[yE, {Index[d$$2, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[yE, 
        {Index[d$$2, Flavor]}, 0]*Coupling[Yu, {Index[i1, Flavor], 
         Index[i2, Flavor]}, 0])/(2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (2*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
      (\[Epsilon]*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*Coupling[\[Lambda], {}, 
        0])/(3*Coupling[ME, {}, 0]^2) - 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*Coupling[gL, {}, 0]^2*
       Coupling[yE, {Index[d$$1, Flavor]}, 0]*Coupling[Yu, 
        {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/(6*Coupling[ME, {}, 0]^2) + 
     (hbar*Bar[Coupling[yE, {Index[d$$1, Flavor]}, 0]]*
       Bar[Coupling[Ye, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
       Coupling[yE, {Index[d$$2, Flavor]}, 0]*Coupling[Ye, 
        {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
       Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
       Log[\[Mu]bar2/Coupling[ME, {}, 0]^2])/Coupling[ME, {}, 0]^2, 
   Coupling[cuu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
      Index[i4_, Flavor]}, 0] -> (-8*hbar*Coupling[gY, {}, 0]^4*
      Delta[Index[i1, Flavor], Index[i2, Flavor]]*Delta[Index[i3, Flavor], 
       Index[i4, Flavor]])/(135*Coupling[ME, {}, 0]^2), 
   Coupling[cuW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
   Coupling[cW, {}, 0] -> 0, Coupling[cWt, {}, 0] -> 0}|>
