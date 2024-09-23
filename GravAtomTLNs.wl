(* ::Package:: *)

BeginPackage["GravAtomTLNs`"];


Lovenumber::usage = "Lovenumber[n,li,mi,l] computes the tidal Love number of multipolar index l of a gravitational atom of indices {n,li,mi}"


Mc;
\[Alpha];
MBH;


Begin["`Private`"];


C1[l_,m_,li_,mi_,k_]:=(-1)^(mi+m) Sqrt[((2l+1)(2li+1)(2Abs[l-li]+4k+1))/(4\[Pi])]ThreeJSymbol[{l,0},{li,0},{Abs[l-li]+2k,0}]ThreeJSymbol[{l,-m},{li,-mi},{Abs[l-li]+2k,m+mi}]
C2[l_,m_,li_,mi_,k_]:=(-1)^m Sqrt[((2l+1)(2li+1)(2Abs[l-li]+4k+1))/(4\[Pi])]ThreeJSymbol[{l,0},{li,0},{Abs[l-li]+2k,0}]ThreeJSymbol[{l,-m},{li,-mi},{Abs[l-li]+2k,m-mi}]


IntG1[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2l+2k+2) Hypergeometric1F1[l+2k-n-2li,2l-2li+4k+2,x]HypergeometricU[-n,2li+2,x],{x,0,y},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntG2[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2l+2k+2) HypergeometricU[l+2k-n-2li,2l-2li+4k+2,x]HypergeometricU[-n,2li+2,x],{x,y,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntegrandG[y_?NumericQ,k_,n_,li_,l_]:=Exp[-y]y^(2l+2k+2) HypergeometricU[-n,2li+2,y](HypergeometricU[l+2k-n-2li,2l-2li+4k+2,y]IntG1[y,k,n,li,l]+Hypergeometric1F1[l+2k-n-2li,2l-2li+4k+2,y]IntG2[y,k,n,li,l])
Gf[k_,n_,li_,l_]:=(4\[Pi])/(2l+1) 1/(n!(n+li+1)(n+2li+1)!) Gamma[l+2k-n-2li]/Gamma[2+2l-2li+4k] NIntegrate[IntegrandG[y,k,n,li,l],{y,0,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]


IntF1[y_?NumericQ,n_,li_,l_]:=NIntegrate[E^-x x^(n+2li+l+3) HypergeometricPFQ[{1,1},{2,2n+2li+3},x]HypergeometricU[-n,2li+2,x],{x,0,y},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntF2[y_?NumericQ,n_,li_,l_]:=NIntegrate[E^-x x^(n+2li+l+2) (Sum[Binomial[2n+2li+1,-s+2n+2li+1]Gamma[s]x^-s,{s,1,2n+2li+1}]-Log[x])HypergeometricU[-n,2li+2,x],{x,y,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntegrandF[y_?NumericQ,n_,li_,l_]:=E^-y y^(n+2li+l+2) HypergeometricU[-n,2li+2,y]((-5-4n-4li+(2+2n+2li)EulerGamma+y) (Gamma[n+2li+l+3]Gamma[n+l+2])/((2+2n+2li)Gamma[l+2])+(Gamma[n+2li+l+4]Gamma[n+l+3])/((2+2n+2li)Gamma[l+3])-1/(2+2n+2li) IntF1[y,n,li,l]-IntF2[y,n,li,l]+(-1)^(n+1) (Sum[Binomial[2n+2li+1,-s+2n+2li+1]Gamma[s]y^-s,{s,1,2n+2li+1}]-Log[y])Sum[Binomial[n,p]Pochhammer[2li+2+p,n-p](-1)^p ((n+2li+l+p+2)!-Sum[(n+2li+l+p+2)!/q! y^q E^-y,{q,0,n+2li+l+p+2}])
,{p,0,n}]+(-1)^(n+1) y/(2+2n+2li) HypergeometricPFQ[{1,1},{2,2n+2li+3},y]Sum[Binomial[n,p]Pochhammer[2li+2+p,n-p](-1)^p Sum[(n+2li+l+p+2)!/q! y^q E^-y,{q,0,n+2li+l+p+2}],{p,0,n}])
Ff[n_,li_,l_]:=(4\[Pi])/(2l+1) 1/(n!(n+li+1)(n+2li+1)!) 1/Gamma[2+2n+2li] NIntegrate[IntegrandF[y,n,li,l],{y,0,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]


Af[k_,n_,li_,l_]:=(8\[Pi])/(2l+1) 1/(Factorial[n](n+li+1)Factorial[n+2li+1]) (1/(Gamma[2+2li-2l+4k](n+li+1)Factorial[n+l-2k]))(Sum[Sum[(-1)^s Binomial[n+l-2k,p]Binomial[n,s-p](Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2li-2l+4k+2,p] Factorial[2li+2k+s+2]Sum[Sum[(-1)^j Binomial[n+l-2k,v]Binomial[n,j-v]Pochhammer[2li-2l+4k+2+v,n+l-2k-v]Pochhammer[2li+2+j-v,n-j+v](Factorial[2li+2k+2+j]-Sum[Factorial[2li+2k+2+q+j]/(2^(2li+2k+3+q+j)Factorial[q]),{q,0,2li+2k+s+2}]),{v,0,j}],{j,0,2n+l-2k}]+Pochhammer[2li-2l+4k+2+p,n+l-2k-p]Pochhammer[2li+2+s-p,n-s+p]Factorial[2li+2k+s+2]Sum[Sum[(-1)^j Binomial[n+l-2k,v]Binomial[n,j-v] Pochhammer[2li+2+j-v,n-j+v]/Pochhammer[2li-2l+4k+2,v] Sum[Factorial[2li+2k+2+q+j]/(2^(2li+2k+3+q+j)Factorial[q]),{q,0,2li+2k+s+2}],{v,0,j}],{j,0,2n+l-2k}]),{p,0,s}],{s,0,2n+l-2k}])


IntC1[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2li+2k+2)Hypergeometric1F1[2k-n-l,2li-2l+4k+2,x]HypergeometricU[-n,2li+2,x],{x,0,y},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntC2[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2li+2k+2)HypergeometricU[2k-n-l,2li-2l+4k+2,x]HypergeometricU[-n,2li+2,x],{x,y,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntegrandC[y_?NumericQ,k_,n_,li_,l_]:=Exp[-y]y^(2li+2k+2)HypergeometricU[-n,2li+2,y](HypergeometricU[2k-n-l,2li-2l+4k+2,y]IntC1[y,k,n,li,l]+Hypergeometric1F1[2k-n-l,2li-2l+4k+2,y]IntC2[y,k,n,li,l])
Cf[k_,n_,li_,l_]:=(4\[Pi])/(2l+1) 1/(Factorial[n](n+li+1)Factorial[n+2li+1]) Gamma[2k-n-l]/Gamma[2+2li-2l+4k] NIntegrate[IntegrandC[y,k,n,li,l],{y,0,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]


Hf[k_,n_,li_,l_]:=(8\[Pi])/(2l+1) 1/(Factorial[n](n+li+1)Factorial[n+2li+1]) (1/(Gamma[2+2l-2li+4k](n+li+1)Factorial[n+2li-l-2k]))(Sum[Sum[(-1)^s Binomial[n+2li-l-2k,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2l-2li+4k+2,p] Factorial[2l+2k+s+2](Sum[(-1)^j Binomial[n+2li-l-2k,j]Pochhammer[2l-2li+4k+2+j,n+2li-l-2k-j]Sum[(-1)^v Binomial[n,v]Pochhammer[2li+2+v,n-v]Factorial[2l+2k+2+j+v],{v,0,n}],{j,0,n+2li-l-2k}]-Sum[Sum[Binomial[n+2li-l-2k,q]Pochhammer[2l-2li+4k+2+q,n+2li-l-2k-q] (-1)^q/Factorial[j-q] Sum[(-1)^v Binomial[n,v]Pochhammer[2li+2+v,n-v] Factorial[2l+2k+2+j+v]/2^(2l+2k+3+j+v),{v,0,n}],{q,0,j}],{j,0,n+2li+l+s+2}]),{p,0,s}],{s,0,2n+2li-l-2k}]+Sum[Sum[(-1)^s Binomial[n+2li-l-2k,p]Binomial[n,s-p]Pochhammer[2l-2li+4k+2+p,n+2li-2k-p]Pochhammer[2li+2+s-p,n-s+p]Factorial[2l+2k+s+2]Sum[Sum[Binomial[n+2li-l-2k,q]((-1)^q/(Pochhammer[2l-2li+4k+2,q]Factorial[j-q]))Sum[(-1)^v Binomial[n,v]Pochhammer[2li+2+v,n-v] Factorial[2l+2k+2+j+v]/2^(2l+2k+3+j+v),{v,0,n}],{q,0,j}],{j,0,n+2li+l+s+2}],{p,0,s}],{s,0,2n+2li-l-2k}])


IntJ1[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2l+2k+2) (D[Hypergeometric1F1[a,2l-2li+4k+2,x],a]/.a->-n-2li+l+2k)HypergeometricU[-n,2li+2,x],{x,0,y},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntJ2[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2l+2k+2) (D[HypergeometricU[a,2l-2li+4k+2,x],a]/.a->-n-2li+l+2k)HypergeometricU[-n,2li+2,x],{x,y,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntegrandJ[y_?NumericQ,k_,n_,li_,l_]:=Exp[-y]y^(2l+2k+2) HypergeometricU[-n,2li+2,y](1/(n+li+1) HypergeometricU[-n-2li+l+2k,2l-2li+4k+2,y](2l-2li+4k+5/2-y/2+PolyGamma[n+2li-l-2k+1])Sum[Sum[(-1)^s Binomial[n+2li-l-2k,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2l-2li+4k+2,p] (Factorial[2l+2k+s+2]-Sum[Factorial[2l+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+2}]),{p,0,s}],{s,0,2n+2li-l-2k}]-1/(2(n+li+1)) HypergeometricU[-n-2li+l+2k,2l-2li+4k+2,y]Sum[Sum[(-1)^s Binomial[n+2li-l-2k,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2l-2li+4k+2,p] (Factorial[2l+2k+s+3]-Sum[Factorial[2l+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+3}]),{p,0,s}],{s,0,2n+2li-l-2k}]+(-n-2li+l+2k)/((2l-2li+4k+2)(n+li+1)) HypergeometricU[-n-2li+l+2k,2l-2li+4k+2,y]Sum[Sum[(-1)^s Binomial[n+2li-l-2k-1,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2l-2li+4k+3,p] (Factorial[2l+2k+s+3]-Sum[Factorial[2l+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+3}]),{p,0,s}],{s,0,2n+2li-l-2k-1}]+(-1)^n HypergeometricU[-n-2li+l+2k,2l-2li+4k+2,y]IntJ1[y,k,n,li,l]+((n+2li-l-2k)/(n+li+1) y HypergeometricU[-n-2li+l+2k+1,2l-2li+4k+3,y]+(D[HypergeometricU[a,2l-2li+4k+2,y],a]/.a->-n-2li+l+2k))Sum[Sum[(-1)^s Binomial[n+2li-l-2k,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2l-2li+4k+2,p] (Factorial[2l+2k+s+2]-Sum[Factorial[2l+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+2}]),{p,0,s}],{s,0,2n+2li-l-2k}]+1/(n+li+1) Hypergeometric1F1[-n-2li+l+2k,2l-2li+4k+2,y](2l-2li+4k+5/2-y/2+PolyGamma[n+2li-l-2k+1])Sum[Sum[(-1)^(n+l+s) Binomial[n+2li-l-2k,p]Binomial[n,s-p]Pochhammer[2l-2li+4k+2+p,n+2li-l-2k-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2l+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+2}],{p,0,s}],{s,0,2n+2li-l-2k}]-1/(2(n+li+1)) Hypergeometric1F1[-n-2li+l+2k,2l-2li+4k+2,y]Sum[Sum[(-1)^(n+l+s) Binomial[n+2li-l-2k,p]Binomial[n,s-p]Pochhammer[2l-2li+4k+2+p,n+2li-l-2k-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2l+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+3}],{p,0,s}],{s,0,2n+2li-l-2k}]+((-n-2li+l+2k)/(n+li+1) y Hypergeometric1F1[-n-2li+l+2k+1,2l-2li+4k+3,y]+(D[Hypergeometric1F1[a,2l-2li+4k+2,y],a]/.a->-n-2li+l+2k))Sum[Sum[(-1)^(n+l+s) Binomial[n+2li-l-2k,p]Binomial[n,s-p]Pochhammer[2l-2li+4k+2+p,n+2li-l-2k-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2l+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+2}],{p,0,s}],{s,0,2n+2li-l-2k}]+(n+2li-l-2k)/(n+li+1) Hypergeometric1F1[-n-2li+l+2k,2l-2li+4k+2,y]Sum[Sum[(-1)^(n+l+1+s) Binomial[n+2li-l-2k-1,p]Binomial[n,s-p]Pochhammer[2l-2li+4k+3+p,n+2li-l-2k-1-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2l+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2l+2k+s+3}],{p,0,s}],{s,0,2n+2li-l-2k-1}]+(-1)^n Hypergeometric1F1[-n-2li+l+2k,2l-2li+4k+2,y]IntJ2[y,k,n,li,l])
Jf[k_,n_,li_,l_]:=(4\[Pi])/(2l+1) (-1)^l/(Factorial[n](n+li+1)Factorial[n+2li+1]) (1/(Gamma[2+2l-2li+4k]Factorial[n+2li-l-2k]))NIntegrate[IntegrandJ[y,k,n,li,l],{y,0,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]


IntB2[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2li+2k+2)(D[HypergeometricU[a,2li-2l+4k+2,x],a]/.a->-n-l+2k)HypergeometricU[-n,2li+2,x],{x,y,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntB1[y_?NumericQ,k_,n_,li_,l_]:=NIntegrate[Exp[-x]x^(2li+2k+2)(D[Hypergeometric1F1[a,2li-2l+4k+2,x],a]/.a->-n-l+2k)HypergeometricU[-n,2li+2,x],{x,0,y},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]
IntegrandB[y_?NumericQ,k_,n_,li_,l_]:=Exp[-y]y^(2li+2k+2)HypergeometricU[-n,2li+2,y](1/(n+li+1) HypergeometricU[-n-l+2k,2li-2l+4k+2,y](2li-2l+4k+5/2-y/2+PolyGamma[n+l-2k+1])Sum[Sum[(-1)^(n+s)Binomial[n+l-2k,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2li-2l+4k+2,p] (Factorial[2li+2k+s+2]-Sum[Factorial[2li+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+2}]),{p,0,s}],{s,0,2n+l-2k}]-1/(2(n+li+1)) HypergeometricU[-n-l+2k,2li-2l+4k+2,y]Sum[Sum[(-1)^(n+s)Binomial[n+l-2k,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2li-2l+4k+2,p] (Factorial[2li+2k+s+3]-Sum[Factorial[2li+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+3}]),{p,0,s}],{s,0,2n+l-2k}]+(-n-l+2k)/((n+li+1)(2li-2l+4k+2)) HypergeometricU[-n-l+2k,2li-2l+4k+2,y]Sum[Sum[(-1)^(n+s)Binomial[n+l-2k-1,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2li-2l+4k+3,p] (Factorial[2li+2k+s+3]-Sum[Factorial[2li+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+3}]),{p,0,s}],{s,0,2n+l-2k-1}]+HypergeometricU[-n-l+2k,2li-2l+4k+2,y]IntB1[y,k,n,li,l]+((n+l-2k)/(n+li+1) y HypergeometricU[-n-l+2k+1,2li-2l+4k+3,y]+(D[HypergeometricU[a,2li-2l+4k+2,y],a]/.a->-n-l+2k))Sum[Sum[(-1)^(n+s)Binomial[n+l-2k,p]Binomial[n,s-p] Pochhammer[2li+2+s-p,n-s+p]/Pochhammer[2li-2l+4k+2,p] (Factorial[2li+2k+s+2]-Sum[Factorial[2li+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+2}]),{p,0,s}],{s,0,2n+l-2k}]+1/(n+li+1) Hypergeometric1F1[-n-l+2k,2li-2l+4k+2,y](2li-2l+4k+5/2-y/2+PolyGamma[n+l-2k+1])Sum[Sum[(-1)^(l+s)Binomial[n+l-2k,p]Binomial[n,s-p]Pochhammer[2li-2l+4k+2+p,n+l-2k-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2li+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+2}],{p,0,s}],{s,0,2n+l-2k}]-1/(2(n+li+1)) Hypergeometric1F1[-n-l+2k,2li-2l+4k+2,y]Sum[Sum[(-1)^(l+s)Binomial[n+l-2k,p]Binomial[n,s-p]Pochhammer[2li-2l+4k+2+p,n+l-2k-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2li+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+3}],{p,0,s}],{s,0,2n+l-2k}]+((-n-l+2k)/(n+li+1) y Hypergeometric1F1[-n-l+2k+1,2li-2l+4k+3,y]+(D[Hypergeometric1F1[a,2li-2l+4k+2,y],a]/.a->-n-l+2k))Sum[Sum[(-1)^(l+s)Binomial[n+l-2k,p]Binomial[n,s-p]Pochhammer[2li-2l+4k+2+p,n+l-2k-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2li+2k+s+2]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+2}],{p,0,s}],{s,0,2n+l-2k}]+(n+l-2k)/(n+li+1) Hypergeometric1F1[-n-l+2k,2li-2l+4k+2,y]Sum[Sum[(-1)^(l+s+1)Binomial[n+l-2k-1,p]Binomial[n,s-p]Pochhammer[2li-2l+4k+3+p,n+l-2k-1-p]Pochhammer[2li+2+s-p,n-s+p]Sum[Factorial[2li+2k+s+3]/Factorial[q] y^q Exp[-y],{q,0,2li+2k+s+3}],{p,0,s}],{s,0,2n+l-2k-1}]+Hypergeometric1F1[-n-l+2k,2li-2l+4k+2,y]IntB2[y,k,n,li,l])
Bf[k_,n_,li_,l_]:=(4\[Pi])/(2l+1) (-1)^(n+l)/(Factorial[n](n+li+1)Factorial[n+2li+1]) (1/(Gamma[2+2li-2l+4k](n+li+1)Factorial[n+l-2k]))NIntegrate[IntegrandB[y,k,n,li,l],{y,0,\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0},PrecisionGoal->4]


fnd[]:=Print["Undefined"]
f1[n_,li_,mi_,l_]:=((n+li+1)/2)^(2l+2) Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Bf[k,n,li,l],{k,0,l}] Mc/(MBH \[Alpha]^(4l+2))
f2[n_,li_,mi_,l_]:=((n+li+1)/2)^(2l+2) Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Jf[k,n,li,l],{k,0,li}] Mc/(MBH \[Alpha]^(4l+2))
f3[n_,li_,mi_,l_]:=Which[EvenQ[n],((n+li+1)/2)^(2n+2) (Sum[C1[n,0,li,mi,k]C1[n,0,li,mi,k]Jf[k,n,li,n],{k,0,li-1}]-C1[n,0,li,mi,li]C1[n,0,li,mi,li]Ff[n,li,n]) Mc/(MBH \[Alpha]^(4n+2)),OddQ[n],fnd[]]
f4[n_,li_,mi_,l_]:=Which[EvenQ[n],((n+li+1)/2)^(2l+2) (Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Jf[k,n,li,l],{k,0,li-(l-n)/2-1}]-C1[l,0,li,mi,li+(n-l)/2]C1[l,0,li,mi,li+(n-l)/2]Ff[n,li,l]+Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Gf[k,n,li,l],{k,li-(l-n)/2+1,li}]) Mc/(MBH \[Alpha]^(4l+2)),OddQ[n],((n+li+1)/2)^(2l+2) (Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Jf[k,n,li,l],{k,0,(2li+n-l-1)/2}]+Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Gf[k,n,li,l],{k,(2li+n-l+1)/2,li}]) Mc/(MBH \[Alpha]^(4l+2))]
f5[n_,li_,mi_,l_]:=Which[EvenQ[n],fnd[],OddQ[n],f4[n,li,mi,l]]
f6[n_,li_,mi_,l_]:=Which[EvenQ[n],((n+li+1)/2)^(2n+4li+2) (Sum[C1[n+2li,0,li,mi,k]C1[n+2li,0,li,mi,k]Gf[k,n,li,n+2li],{k,1,li}]-C1[n+2li,0,li,mi,0]C1[n+2li,0,li,mi,0]Ff[n,li,n+2li]) Mc/(MBH \[Alpha]^(4n+8li+2)),OddQ[n],fnd[]]
f7[n_,li_,mi_,l_]:=((n+li+1)/2)^(2l+2) Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Gf[k,n,li,l],{k,0,li}] Mc/(MBH \[Alpha]^(4l+2))
f9[n_,li_,mi_,l_]:=Which[EvenQ[n],((n+li+1)/2)^(2n+2) (Sum[C1[n,0,li,mi,k]C1[n,0,li,mi,k]Bf[k,n,li,n],{k,0,n-1}]-C1[n,0,li,mi,n]C1[n,0,li,mi,n]Ff[n,li,n]) Mc/(MBH \[Alpha]^(4n+2)),OddQ[n],fnd[]]
f21[n_,li_,mi_,l_]:=Which[EvenQ[n],((n+li+1)/2)^(2l+2) (Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Bf[k,n,li,l],{k,0,(n+l)/2-1}]-C1[l,0,li,mi,(n+l)/2]C1[l,0,li,mi,(n+l)/2]Ff[n,li,l]+Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Cf[k,n,li,l],{k,(n+l)/2+1,l}]) Mc/(MBH \[Alpha]^(4l+2)),OddQ[n],((n+li+1)/2)^(2l+2) (Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Bf[k,n,li,l],{k,0,(n+l-1)/2}]+Sum[C1[l,0,li,mi,k]C1[l,0,li,mi,k]Cf[k,n,li,l],{k,(n+l+1)/2,l}]) Mc/(MBH \[Alpha]^(4l+2))]
f23[n_,li_,mi_,l_]:=f21[n,li,mi,l]


f27[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],True,f7[n,li,mi,l]]
f28[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f3[n,li,mi,l],
l>2,f7[n,li,mi,l]]
f29[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],(l>=2&&l<n),f2[n,li,mi,l],
l==n,f3[n,li,mi,l],
l>n,f7[n,li,mi,l]]
f30[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f3[n,li,mi,l],
l==4,f6[n,li,mi,l],
l>4,f7[n,li,mi,l]]
f31[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],(l>=2&&l<n),f2[n,li,mi,l],
l==n,f3[n,li,mi,l],
l==n+1,f5[n,li,mi,l],
l==n+2,f6[n,li,mi,l],
l>n+2,f7[n,li,mi,l]]
f32[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f1[n,li,mi,l],
l==4||l==6,f4[n,li,mi,l],
l>7,f7[n,li,mi,l]]
f33[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l>=2&&l<=li,f1[n,li,mi,l],
l>li&&l<n,f2[n,li,mi,l],
l==n,f3[n,li,mi,l],
l>n&&l<n+2li,f4[n,li,mi,l],
l==n+2li,f6[n,li,mi,l],
l>n+2li,f7[n,li,mi,l]]
f34[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],True,f7[n,li,mi,l]]
f35[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f4[n,li,mi,l],
l>2,f7[n,li,mi,l]]
f36[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f9[n,li,mi,l],
l==4,f4[n,li,mi,l],
l==6,f6[n,li,mi,l],
l>6,f7[n,li,mi,l]]
f37[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l>=2&&l<n,f1[n,li,mi,l],
l==n,f9[n,li,mi,l],
l>n&&l<3n,f4[n,li,mi,l],
l==3n,f6[n,li,mi,l],
l>3n,f7[n,li,mi,l]]
f38[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f21[n,li,mi,l],
l==4,f6[n,li,mi,l],
l>4,f7[n,li,mi,l]]
f39[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l>=2&&l<=li,f21[n,li,mi,l],
l>li&&l<n+2li,f4[n,li,mi,l],
l==n+2li,f6[n,li,mi,l],
l>n+2li,f7[n,li,mi,l]]
f40[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f9[n,li,mi,l],
l>2&&l<=li,f21[n,li,mi,l],
l>li&&l<n+2li,f4[n,li,mi,l],
l==n+2li,f6[n,li,mi,l],
l>n+2li,f7[n,li,mi,l]]
f41[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],2<=l&&l<n,f1[n,li,mi,l],
l==n,f9[n,li,mi,l],
l>n&&l<=li,f21[n,li,mi,l],
l>li&&l<n+2li,f4[n,li,mi,l],
l==n+2li,f6[n,li,mi,l],
l>n+2li,f7[n,li,mi,l]]
f42[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f6[n,li,mi,l],
l>2,f7[n,li,mi,l]]
f43[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f23[n,li,mi,l],
l==4,f4[n,li,mi,l],
l>4,f7[n,li,mi,l]]
f44[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l>=2&&l<=li,f23[n,li,mi,l],
l>li&&l<n+2li,f4[n,li,mi,l],
l==n+2li,f6[n,li,mi,l],
l>n+2li,f7[n,li,mi,l]]
f45[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f9[n,li,mi,l],
l==4||l==6,f4[n,li,mi,l],
l==8,f6[n,li,mi,l],
l>8,f7[n,li,mi,l]]
f46[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l==2,f9[n,li,mi,l],
l>2&&l<=li,f23[n,li,mi,l],
l>li&&l<2+2li,f4[n,li,mi,l],
l==2+2li,f6[n,li,mi,l],
l>2+2li,f7[n,li,mi,l]]
f47[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l>=2&&l<n,f1[n,li,mi,l],
l==n,f9[n,li,mi,l],
l==li,f23[n,li,mi,l],
l>li&&l<n+2li,f4[n,li,mi,l],
l==n+2li,f6[n,li,mi,l],
l>n+2li,f7[n,li,mi,l]]
f48[n_,li_,mi_,l_]:=Which[OddQ[l],fnd[],l>=2&&l<n,f1[n,li,mi,l],
l==n,f9[n,li,mi,l],
l>n&&l<=li,f23[n,li,mi,l],
l>li&&l<n+2li,f4[n,li,mi,l],
l==n+2li,f6[n,li,mi,l],
l>n+2li,f7[n,li,mi,l]]


Lovenumber[n_,li_,mi_,l_]:=Which[
(n==0&&li==0&&mi==0),f34[n,li,mi,l],
(n==0&&li==1),f42[n,li,mi,l],
(n==0&&li==2),f38[n,li,mi,l],
(n==0&&li>1&&OddQ[li]),f44[n,li,mi,l],
(n==0&&li>2&&EvenQ[li]),f39[n,li,mi,l],

(n==1&&li==0&&mi==0),f27[n,li,mi,l],
(n==1&&li==1),f35[n,li,mi,l],
(n==1&&li==2),f43[n,li,mi,l],
(n==1&&li>1&&OddQ[li]),f39[n,li,mi,l],
(n==1&&li>2&&EvenQ[li]),f44[n,li,mi,l],

(n==2&&li==0&&mi==0),f28[n,li,mi,l],
(n==2&&li==1),f30[n,li,mi,l],
(n==2&&li==2),f36[n,li,mi,l],
(n==2&&li==3),f45[n,li,mi,l],
(n==2&&li>3&&OddQ[li]),f46[n,li,mi,l],
(n==2&&li>2&&EvenQ[li]),f40[n,li,mi,l],

(n==3&&li==2),f32[n,li,mi,l],

(n>2&&li==0&&mi==0),f29[n,li,mi,l],
(n>2&&li==1),f31[n,li,mi,l],
(n>3&&li>=2&&n>li),f33[n,li,mi,l],
(n==li&&li>2),f37[n,li,mi,l],
(n>2&&li>n&&EvenQ[li]==EvenQ[n]),f41[n,li,mi,l],
(n>2&&li>n+1&&EvenQ[li]==OddQ[n]),f48[n,li,mi,l],
(n>2&&li==n+1),f47[n,li,mi,l]
]


End[];


EndPackage[];
