(* ::Package:: *)

BeginPackage["LyapunovCalc`"]
 
CalculateLCE::usage="CalculateLCE[F, xinit, tmax, tstep]
Finds LCEs of arbitrary dimensioned system"
CalculateLCE::BadDimensions="One or all of arguments have different number of dimensions"

Begin["Private`"]


CalculateLCE[F_, xinit_List, tmax_, tstep_,\[CurlyEpsilon]_:1*^-3 ,tWindow_Integer:10]:=Module[{x0, \[CapitalPhi]0, f0,
             dim, \[CapitalPhi], eqs, \[Phi], \[Xi], \[CapitalXi], sol, tT=0, ws={}, lces={}, cnorm},
 x0=xinit;
 dim=Length[xinit];
 \[CapitalXi][t_]=Table[\[Xi][i][t],{i,dim}];
 (*Jacobi matrice*)
 \[CapitalPhi][t_]=Table[\[Phi][i,j][t],{i,dim},{j,dim}];
 eqs:={D[\[CapitalXi][t],t]==F@@\[CapitalXi][t], D[\[CapitalPhi][t],t]==\[CapitalPhi][t].Transpose@Outer[D, F@@\[CapitalXi][t], \[CapitalXi][t]],\[CapitalXi][0]==x0,\[CapitalPhi][0]==\[CapitalPhi]0};
 \[CapitalPhi]0=IdentityMatrix[dim];
 While[tT<tmax&&(Length[ws]<tWindow||Apply[And,StandardDeviation@#>=\[CurlyEpsilon] Abs@Mean@#&/@Transpose[lces[[-tWindow;;]]]]),
 sol=First[NDSolve[eqs,Join[\[CapitalXi][t],Flatten[\[CapitalPhi][t]]],{t,0,tstep}]];tT+=tstep;
 x0=\[CapitalXi][t]/.sol/.t->tstep;
 (*Stupidly obfuscated Gram Schmidt orthogonalization*)
 \[CapitalPhi]0=Last@FoldList[Append[#1,#2-If[Length@#1!=0,Plus@@(#1.#2 #1/(Plus@@@(#1 #1))),0]]&,{},\[CapitalPhi][t]/.sol/.t->tstep];
 cnorm=Norm/@\[CapitalPhi]0;
 ws=Append[ws,cnorm];
 \[CapitalPhi]0/=cnorm;
 lces=Append[lces,(Plus@@Log[ws])/tT];
 ];
 lces
];

End[]

EndPackage[]






