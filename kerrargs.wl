(* ::Package:: *)

(*The maths describing the path of an individual photon*)
KerrMap[\[Phi]c_,\[Theta]c_,\[Delta]_,M_,r0_,\[Theta]0_,\[Phi]0_,a_]:=
Module[{\[CapitalSigma]0,\[CapitalDelta]0,td0,R,\[CapitalTheta],\[CapitalPhi],r,\[Theta],\[Phi],t,\[Sigma],rd0,\[Theta]d0,\[Phi]d0,tconst,\[Phi]const,test,domain,k,end},
\[CapitalSigma]0=r0^2+a^2 Cos[\[Theta]0]^2;
\[CapitalDelta]0=r0^2-2r0+a^2;

(*Initial Velocities*)
td0=1;

rd0=-Cos[\[Phi]c] Sin[\[Theta]c] \[Sqrt]((\[CapitalDelta]0 (r0^2 (-2 r0+\[CapitalSigma]0)+a^2 (-r0+\[CapitalSigma]0)-a^2 r0 Cos[2 \[Theta]0]) (16 a r0 Abs[-2 r0+\[CapitalSigma]0] Sin[\[Phi]c] \[Sqrt](\[CapitalSigma]0 (r0^2 (-2 r0+\[CapitalSigma]0)+a^2 (-r0+\[CapitalSigma]0)-a^2 r0 Cos[2 \[Theta]0])) Sin[\[Theta]0] Sin[\[Theta]c]+(2 r0-\[CapitalSigma]0) (-2 a^2 r0^2+4 a^2 r0 \[CapitalSigma]0+8 r0^3 \[CapitalSigma]0-4 a^2 \[CapitalSigma]0^2-4 r0^2 \[CapitalSigma]0^2+2 a^2 r0^2 Cos[2 \[Phi]c]+a^2 r0 Cos[2 \[Theta]0] (2 r0+4 \[CapitalSigma]0-2 r0 Cos[2 \[Theta]c]+r0 Cos[2 \[Theta]c-2 \[Phi]c]-2 r0 Cos[2 \[Phi]c]+r0 Cos[2 (\[Theta]c+\[Phi]c)])+4 a^2 r0^2 Cos[2 \[Theta]c] Sin[\[Phi]c]^2)))/(\[CapitalSigma]0 (-2 (a^2+r0^2) (2 r0-\[CapitalSigma]0) \[CapitalSigma]0+a^2 r0 (-2 r0+4 \[CapitalSigma]0+2 r0 Cos[2 \[Theta]c]-r0 Cos[2 \[Theta]c-2 \[Phi]c]+2 r0 Cos[2 \[Phi]c]-r0 Cos[2 (\[Theta]c+\[Phi]c)]) Sin[\[Theta]0]^2)^2));
\[Theta]d0=-Cos[\[Theta]c] \[Sqrt](((r0^2 (-2 r0+\[CapitalSigma]0)+a^2 (-r0+\[CapitalSigma]0)-a^2 r0 Cos[2 \[Theta]0]) (16 a r0 Abs[-2 r0+\[CapitalSigma]0]Sin[\[Phi]c] \[Sqrt](\[CapitalSigma]0 (r0^2 (-2 r0+\[CapitalSigma]0)+a^2 (-r0+\[CapitalSigma]0)-a^2 r0 Cos[2 \[Theta]0])) Sin[\[Theta]0] Sin[\[Theta]c]+(2 r0-\[CapitalSigma]0) (-2 a^2 r0^2+4 a^2 r0 \[CapitalSigma]0+8 r0^3 \[CapitalSigma]0-4 a^2 \[CapitalSigma]0^2-4 r0^2 \[CapitalSigma]0^2+2 a^2 r0^2 Cos[2 \[Phi]c]+a^2 r0 Cos[2 \[Theta]0] (2 r0+4 \[CapitalSigma]0-2 r0 Cos[2 \[Theta]c]+r0 Cos[2 \[Theta]c-2 \[Phi]c]-2 r0 Cos[2 \[Phi]c]+r0 Cos[2 (\[Theta]c+\[Phi]c)])+4 a^2 r0^2 Cos[2 \[Theta]c] Sin[\[Phi]c]^2)))/(\[CapitalSigma]0 (-2 (a^2+r0^2) (2 r0-\[CapitalSigma]0) \[CapitalSigma]0+a^2 r0 (-2 r0+4 \[CapitalSigma]0+2 r0 Cos[2 \[Theta]c]-r0 Cos[2 \[Theta]c-2 \[Phi]c]+2 r0 Cos[2 \[Phi]c]-r0 Cos[2 (\[Theta]c+\[Phi]c)]) Sin[\[Theta]0]^2)^2));
\[Phi]d0=(Csc[\[Theta]0]^2 Sin[\[Theta]c] (Abs[-2 r0+\[CapitalSigma]0] Sin[\[Phi]c] Sqrt[\[CapitalSigma]0 (r0^2 (-2 r0+\[CapitalSigma]0)+a^2 (-r0+\[CapitalSigma]0)-a^2 r0 Cos[2 \[Theta]0])] Csc[\[Theta]0]+2 a r0 (-2 r0+\[CapitalSigma]0) Sin[\[Theta]c] Sin[\[Phi]c]^2))/(4 a^2 r0^2 Cos[\[Phi]c]^2-(a^2+r0^2) (2 r0-\[CapitalSigma]0) \[CapitalSigma]0 Csc[\[Theta]0]^2+2 a^2 r0 (-2 r0+\[CapitalSigma]0+2 r0 Cos[\[Theta]c]^2 Sin[\[Phi]c]^2));
tconst=-2(1-(2 r0)/\[CapitalSigma]0)td0-(4 r0 a Sin[\[Theta]0]^2)/\[CapitalSigma]0 \[Phi]d0;
\[Phi]const=2(r0^2+a^2+(2 r0 a^2)/\[CapitalSigma]0 Sin[\[Theta]0]^2)Sin[\[Theta]0]^2 \[Phi]d0-(4 r0 a Sin[\[Theta]0]^2)/\[CapitalSigma]0 td0;

(*Test is a Boolean: 1 if geodesic doesn't cross event horizon,
0 becomes 0 if geodesic crosses event horizon*)
test=1;

{R,\[CapitalTheta],\[CapitalPhi]}={r,\[Theta],\[Phi]}/.NDSolve[{
t'[\[Sigma]]-(2 a r[\[Sigma]] (\[Phi]const+a tconst Sin[\[Theta][\[Sigma]]]^2)+tconst (a^2+r[\[Sigma]]^2) \[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]])/(4 r[\[Sigma]] (a^2 Cos[\[Theta][\[Sigma]]]^2+r[\[Sigma]]^2)-2 (a^2+r[\[Sigma]]^2) \[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]])==0,
-((\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]] r'[\[Sigma]]^2 \[CapitalDelta]'[r[\[Sigma]]])/\[CapitalDelta][r[\[Sigma]]]^2)+(4 a Sin[\[Theta][\[Sigma]]]^2 t'[\[Sigma]] \[Phi]'[\[Sigma]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]+(2 \[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]] r''[\[Sigma]])/\[CapitalDelta][r[\[Sigma]]]-(r'[\[Sigma]]^2 \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalDelta][r[\[Sigma]]]-\[Theta]'[\[Sigma]]^2 \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]]-(4 a r[\[Sigma]] Sin[\[Theta][\[Sigma]]]^2 t'[\[Sigma]] \[Phi]'[\[Sigma]] \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]^2-t'[\[Sigma]]^2 (2/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]-(2 r[\[Sigma]] \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]^2)-Sin[\[Theta][\[Sigma]]]^2 \[Phi]'[\[Sigma]]^2 (2 r[\[Sigma]]+(2 a^2 Sin[\[Theta][\[Sigma]]]^2)/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]-(2 a^2 r[\[Sigma]] Sin[\[Theta][\[Sigma]]]^2 \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]^2)+(2 r'[\[Sigma]] (\[Theta]'[\[Sigma]] \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]]+r'[\[Sigma]] \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]]))/\[CapitalDelta][r[\[Sigma]]]==0,
(8 a Cos[\[Theta][\[Sigma]]] r[\[Sigma]] Sin[\[Theta][\[Sigma]]] t'[\[Sigma]] \[Phi]'[\[Sigma]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]-2 Cos[\[Theta][\[Sigma]]] Sin[\[Theta][\[Sigma]]] (a^2+r[\[Sigma]]^2+(2 a^2 r[\[Sigma]] Sin[\[Theta][\[Sigma]]]^2)/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]) \[Phi]'[\[Sigma]]^2+2 \[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]] \[Theta]''[\[Sigma]]-(r'[\[Sigma]]^2 \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalDelta][r[\[Sigma]]]+(2 r[\[Sigma]] t'[\[Sigma]]^2 \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]^2-\[Theta]'[\[Sigma]]^2 \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]]-(4 a r[\[Sigma]] Sin[\[Theta][\[Sigma]]]^2 t'[\[Sigma]] \[Phi]'[\[Sigma]] \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]^2-Sin[\[Theta][\[Sigma]]]^2 \[Phi]'[\[Sigma]]^2 ((4 a^2 Cos[\[Theta][\[Sigma]]] r[\[Sigma]] Sin[\[Theta][\[Sigma]]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]-(2 a^2 r[\[Sigma]] Sin[\[Theta][\[Sigma]]]^2 \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])/\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]^2)+2 \[Theta]'[\[Sigma]] (\[Theta]'[\[Sigma]] \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]]+r'[\[Sigma]] \!\(\*SuperscriptBox[\(\[CapitalSigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[r[\[Sigma]],\[Theta][\[Sigma]]])==0,
\[Phi]'[\[Sigma]]-(2 a tconst r[\[Sigma]]+\[Phi]const Csc[\[Theta][\[Sigma]]]^2 (2 r[\[Sigma]]-\[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]]))/(4 r[\[Sigma]] (a^2 Cos[\[Theta][\[Sigma]]]^2+r[\[Sigma]]^2)-2 (a^2+r[\[Sigma]]^2) \[CapitalSigma][r[\[Sigma]],\[Theta][\[Sigma]]])==0,
r'[0]==rd0,
\[Theta]'[0]==\[Theta]d0,
r[0]==r0,
t[0]==0,
\[Theta][0]==\[Theta]0,
\[Phi][0]==\[Phi]0,
WhenEvent[r[\[Sigma]]-(1+Sqrt[1-a^2]+\[Delta])==0,(Clear[test];test:=0;Throw[end=\[Sigma],"StopIntegration"])]},
{r,\[Theta],\[Phi]},{\[Sigma],-200,0},

Method->{"EquationSimplification"->"Residual"}
][[1]];

If[test==1,
{If[Mod[\[CapitalTheta][-200],2\[Pi]]<=\[Pi],Mod[-\[CapitalPhi][-200],2\[Pi]],Mod[\[CapitalPhi][-200],2\[Pi]]+\[Pi]],\[Pi]-If[Mod[\[CapitalTheta][R["Domain"][[1,1]]],2\[Pi]]<=\[Pi],Mod[\[CapitalTheta][R["Domain"][[1,1]]],2\[Pi]],2\[Pi]-Mod[\[CapitalTheta][R["Domain"][[1,1]]],2\[Pi]]]},
{10\[Pi],10\[Pi]}
]

];


(*pic is the picture you want to warp*)
pic=Import[$CommandLine[[4]]];
{width,height}=ImageDimensions[pic];
asprat=height/width;


(*M: Black hole mass. Can be left as 1. This physics is independent of scale*)
M=1;
(*Observer Position in spherical polars*)
r0=ToExpression[$CommandLine[[5]]];
\[Theta]0=ToExpression[$CommandLine[[6]]]\[Pi];
\[Phi]0=ToExpression[$CommandLine[[7]]]\[Pi];
(*How close we get to event horizon before ending integration*)
\[Delta]=.1;
(*Spin parameter: 0=> Schwarzschild black hole
.999 ... max spin*)
a=ToExpression[$CommandLine[[8]]];


\[CapitalSigma][r_,\[Theta]_]:=r^2+a^2 Cos[\[Theta]]^2;
\[CapitalDelta][r_]:=r^2-2r+a^2;

(*SplitTransform: warps image*)
(*rows, columns: split image into sections to process in parallel*)
(*hmargin, vmargin: cut of fraction of sides 
(i.e. zooming in on image, cuts down on processing time)*)
SplitTransform[image_,rows_,columns_,hmargin_,vmargin_]:=
ParallelTable[ImageTransformation[image,KerrMap[2\[Pi] #[[1]],\[Pi]  #[[2]]/asprat,\[Delta],M,r0,\[Theta]0,\[Phi]0,a]{1/(2\[Pi]),asprat/\[Pi]}&,
PlotRange->{{hmargin+n (1-2hmargin)/columns,hmargin+(n+1) (1-2hmargin)/columns},asprat{1-vmargin-(m+1)((1-2vmargin)/rows),1-vmargin-m((1-2vmargin)/rows)}}],{m,0,rows-1,1},{n,0,columns-1,1},
Method->"CoarsestGrained"]//ImageAssemble;


warp=SplitTransform[pic,8,5,1/3,1/4]


Export[$CommandLine[[9]],warp]
