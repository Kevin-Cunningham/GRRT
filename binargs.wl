(* ::Package:: *)

(*The maths describing the path of an individual photon*)
(*See arXiv:1603.04469v3 [gr-qc] 9 Jun 2016*)
U[x_,y_,z_,a_,\[Alpha]_]:=1+1/Sqrt[(x+a/2 Cos[\[Alpha]])^2+(y+a/2 Sin[\[Alpha]])^2+(z)^2]+1/Sqrt[(x-a/2 Cos[\[Alpha]])^2+(y-a/2 Sin[\[Alpha]])^2+(z)^2];
binmap[\[Phi]_,\[Theta]_,x0_,y0_,z0_,a_,\[Alpha]_]:=
Module[{t,x,y,z,px,py,pz,T,X,Y,Z,eventtest,dom},
eventtest=1;
{T,X,Y,Z}={t,x,y,z}/.NDSolve[{t'[\[Lambda]]==U[x[\[Lambda]],y[\[Lambda]],z[\[Lambda]],a,\[Alpha]]^2,
x'[\[Lambda]]==1/U[x[\[Lambda]],y[\[Lambda]],z[\[Lambda]],a,\[Alpha]]^2 px[\[Lambda]],px'[\[Lambda]]==D[U[x[\[Lambda]],y[\[Lambda]],z[\[Lambda]],a,\[Alpha]]^2,x[\[Lambda]]],
y'[\[Lambda]]==1/U[x[\[Lambda]],y[\[Lambda]],z[\[Lambda]],a,\[Alpha]]^2 py[\[Lambda]],py'[\[Lambda]]==D[U[x[\[Lambda]],y[\[Lambda]],z[\[Lambda]],a,\[Alpha]]^2,y[\[Lambda]]],
z'[\[Lambda]]==1/U[x[\[Lambda]],y[\[Lambda]],z[\[Lambda]],a,\[Alpha]]^2 pz[\[Lambda]],pz'[\[Lambda]]==D[U[x[\[Lambda]],y[\[Lambda]],z[\[Lambda]],a,\[Alpha]]^2,z[\[Lambda]]],
x[0]==x0,y[0]==y0,z[0]==z0,t[0]==0,
px[0]==Sin[\[Theta]]Sin[\[Phi]],py[0]==Sin[\[Theta]]Cos[\[Phi]],pz[0]==-Cos[\[Theta]]
},{t,x,y,z},{\[Lambda],0,100},
Method->{"EventLocator",
"Event"->{Abs[x[\[Lambda]]-a/2 Cos[\[Alpha]]]-.1<0&&Abs[y[\[Lambda]]-a/2 Sin[\[Alpha]]]-.1<0&&Abs[z[\[Lambda]]]-.1<0,Abs[x[\[Lambda]]+a/2 Cos[\[Alpha]]]-.1<0&&Abs[y[\[Lambda]]+a/2 Sin[\[Alpha]]]-.1<0&&Abs[z[\[Lambda]]]-.1<0},
"EventAction":>{
{Clear[eventtest];eventtest:=0;Throw[end=\[Sigma],"StopIntegration"]},
{Clear[eventtest];eventtest:=0;Throw[end=\[Sigma],"StopIntegration"]}}
}
][[1]];
dom=X["Domain"][[1,2]];
If[eventtest==1,{If[X[dom]<0,2\[Pi]-ArcCos[Y[dom]/Sqrt[X[dom]^2+Y[dom]^2]],ArcCos[Y[dom]/Sqrt[X[dom]^2+Y[dom]^2]]],ArcCos[-Z[dom]/Sqrt[X[dom]^2+Y[dom]^2+Z[dom]^2]]},{10,10}]
];


(*pic is the picture you want to warp*)
pic=Import[$CommandLine[[4]]];
{width,height}=ImageDimensions[pic];
asprat=height/width;


(*Observer position*)
x0=ToExpression[$CommandLine[[5]]];
y0=ToExpression[$CommandLine[[6]]];
z0=ToExpression[$CommandLine[[7]]];
(*Distance between black holes*)
a=ToExpression[$CommandLine[[8]]];
(*Angle black holes make with x-axis (I think, it might have the y-axis)*)
\[Alpha]=ToExpression[$CommandLine[[9]]];


(*SplitTransform: warps image*)
(*rows, columns: split image into sections to process in parallel*)
(*hmargin, vmargin: cut of fraction of sides 
(i.e. zooming in on image, cuts down on processing time)*)
SplitTransform[image_,rows_,columns_,hmargin_,vmargin_,\[Alpha]_]:=
ParallelTable[ImageTransformation[image,binmap[2\[Pi] #[[1]],\[Pi]  #[[2]]/asprat,x0,y0,z0,a,\[Alpha]]{1/(2\[Pi]),asprat/\[Pi]}&,PlotRange->{{hmargin+n (1-2hmargin)/columns,hmargin+(n+1) (1-2hmargin)/columns},asprat{1-vmargin-(m+1)((1-2vmargin)/rows),1-vmargin-m((1-2vmargin)/rows)}}],{m,0,rows-1,1},{n,0,columns-1,1},Method->"CoarsestGrained"]//ImageAssemble;


warp = SplitTransform[pic,8,5,1/3,1/4,\[Alpha]];


Export[$CommandLine[[10]],warp]
