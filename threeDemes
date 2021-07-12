(*
 Program that models the likelihood of a gene drive spillover from a \
targeted population to a non target population
In this case, we examine the case of three populations arranged in a \
line, with symmetrical migration between populations
Produces a graph that shows the relationship between fitness cost (s) \
and maximum tolerated migration rate {m^*}
*)

$HistoryLength=0;

(*Given values of fitness cost, CRISPR conversion rate, gene dominance, and migration rate, solves the system of equations of allele frequency in each population and returns a table showing all equilibrium solutions*)
eqpoints1[s_,c_,h_,m_]:=Module[{eeq1,eeq2,eeq3,sn,sc},(*Gives the equilibrium solutions for s,c,h,m*)
sn=1/2 (1-c)(1-h s);
sc=c(1-s);
eeq1=q1==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q1+m*q2)};
eeq2=q2==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2m)q2+m*q1+m*q3)};
eeq3=q3==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q3+m*q2)};
sols=NSolve[{eeq1,eeq2,eeq3},{q1,q2,q3},Reals,Method->"Newton", WorkingPrecision->24];
z1=Table[
If[
Length[ss]!=3,"N",
If[ss[[1,2]]<0-10^-5||ss[[2,2]]<0-10^-5||ss[[3,2]]<0-10^-5||ss[[1,2]]>1.+10^-5||ss[[2,2]]>1.+10^-5||ss[[3,2]]>1.+10^-5,"N",{ss[[1,2]],ss[[2,2]],ss[[3,2]]}]
](*Here we remove all solution that are not in the interval [0,1] with some small error margin due to numeric accuracy*)
,{ss,sols}];
z2=DeleteCases[z1,"N"];
z3=DeleteDuplicates[Round[z2,10^-5]];(*Due to numeric accuracy, some solutions are counted twice sometimes. Here we remove duplicate solutions*)
N[Sort[z3]]
]

(*System of equations representing allele frequency in each population*)
eq1[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q1+m*q2)};
eq2[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2*m)q2+m*q1+m*q3)};
eq3[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q3+m*q2)};

(*Jacobian matrix*)
J[m_,s_,sc_,sn_] := {{D[eq1[m,s,sc,sn], q1],D[eq1[m,s,sc,sn], q2],D[eq1[m,s,sc,sn], q3]}, {D[eq2[m,s,sc,sn], q1], D[eq2[m,s,sc,sn], q2], D[eq2[m,s,sc,sn], q3]},{D[eq3[m,s,sc,sn], q1], D[eq3[m,s,sc,sn], q2], D[eq3[m,s,sc,sn], q3]}};

(*Given values of fitness cost, CRISPR conversion rate, gene dominance, and migration rate, selects the equilibria which represent Differential Targeting (gene drive allele frequency above 50% in the target population and below 50% in all non-target populations)
Then, calculates stability of equilibria by finding eigenvalues of the Jacobian matrix
Returns list of stable DTEs, for both central targeting case and peripheral targeting case*)
DTE[c_,h_,s_,m_]:=Module[{},
sc=c*(1-s); sn=(1-c)*(1-h*s);
nn=eqpoints1[s,c,h,m];
nna=Select[nn,#[[1]]>0.5&&#[[2]]<0.5&&#[[3]]<0.5&];(*potential DTE for target pop peripheral*)
nnb=Select[nn,#[[1]]<0.5&&#[[2]]>0.5&&#[[3]]<0.5&];(*potential DTE for target pop central*)
staa=If[nna!={},Table[J[m,s,sc,sn]/.{q1->nna[[i,1]],q2->nna[[i,2]],q3->nna[[i,3]]},{i,Length[nna]}],"None1"];
sta=If[!StringQ[staa],Table[Eigenvalues[SetPrecision[staa[[j]],MachinePrecision],1,Method->"Arnoldi"],{j,Length[staa]}],"None1"];
stbb=If[nnb!={},Table[J[m,s,sc,sn]/.{q1->nna[[i,1]],q2->nna[[i,2]],q3->nna[[i,3]]},{i,Length[nnb]}],"None1"];
stb=If[!StringQ[stbb],Table[Eigenvalues[SetPrecision[stbb[[j]],MachinePrecision],1,Method->"Arnoldi"],{j,Length[stbb]}],"None1"];
DTEsp=If[!StringQ[sta],DeleteCases[Table[If[Abs[sta[[i,1]]]<1(*&&Abs[sta[[i,2]]]<1&&Abs[sta[[i,3]]]<1*),nna[[i]],"N"],{i,Length[sta]}],"N"],sta];(*DTEs for target pop peripheral*)
DTEsc=If[!StringQ[stb],DeleteCases[Table[If[Abs[stb[[i,1]]]<1(*&&Abs[stb[[i,2]]]<1&&Abs[stb[[i,3]]]<1*),nnb[[i]],"N"],{i,Length[stb]}],"N"],stb];(*DTEs for target pop central*)
If[StringQ[DTEsp],DTEsp="None2"];
If[StringQ[DTEsc],DTEsc="None2"];
{DTEsp,DTEsc}
]

(*Given a value for fitness cost, find the maximum tolerated migration rate for both central targeting case and peripheral targeting case*)
mqstar[ss_]:=Module[{maxm,jumps,t0,t1,t2,t3,outm}, 
cc = 1;
hh = 0;
maxm=0.15;(*Maximal m for searching m^* *)
jumps=10^-3;(*Resolution in m for searching m^*  *)


t0 = ParallelTable[{N[mm],DTE[cc,hh,ss,mm]},{mm,0,maxm,jumps}];
t1=Table[{t[[1]],t[[2]]},{t,t0}];
t1p=Table[{t[[1]],Length[t[[2,1]]]},{t,t1}];
t1c=Table[{t[[1]],Length[t[[2,2]]]},{t,t1}];
mstarp=If[Last[t1p][[2]]==0,Max[{Min[Flatten[Take[Select[t1p,#[[2]]==0&],All,{1}]]]-jumps,0}],"Non"];
mstarc=If[Last[t1c][[2]]==0,Max[{Min[Flatten[Take[Select[t1c,#[[2]]==0&],All,{1}]]]-jumps,0}],"Non"];
{mstarp,mstarc}
]

c=1;h=0; time=100; th=0.01;
dss = 0.01;

SetDirectory["/home/lab-heavy/Yonatan"];
Directory[]
DateString[]
(*Export results*)
res06 = Table[{ss,mqstar[ss]}, {ss, 0.5, 0.6, dss}];
Export["Dropbox/res06.mx", res06];
res07 = Table[{ss,mqstar[ss]}, {ss, 0.6+dss, 0.7, dss}];
Export["Dropbox/res07.mx", res07];
res08 = Table[{ss,mqstar[ss]}, {ss, 0.7+dss, 0.8, dss}];
Export["Dropbox/res08.mx", res08];
res09 = Table[{ss,mqstar[ss]}, {ss, 0.8+dss, 0.9, dss}];
Export["Dropbox/res09.mx", res09];
res10 = Table[{ss,mqstar[ss]}, {ss, 0.9+dss, 1, dss}];
Export["Dropbox/res10.mx", res10];

(*Import results*)
res6 = Import["res06.mx"];
res7 = Import["res07.mx"];
res8 = Import["res08.mx"];
res9 = Import["res09.mx"];
res10 = Import["res10.mx"];

twoResImp = Import["twoDemesRes.mx"];
res = Join[res6, res7, res8, res9, res10];

twoRes = Table[{r[[1]],r[[2]]}, {r,twoResImp}];
resp = Table[{r[[1]],N[Round[r[[2,1]],10^-5]]}, {r,res}];
resp2 = resp/.{{0.64,0.033}->Nothing,{0.8,0.04}->Nothing,{0.87,0.011}->Nothing,{0.88,0.01}->Nothing,{0.93,0.}->Nothing,{0.9400000000000001`,0.006`}->Nothing};
resc = Table[{r[[1]],r[[2,2]]}, {r,res}];

(*Graphics for plot*)
pMax = {0.69,0.13};
cMax = {0.67,0.06};
tMax = {0.72,0.107};
fiveCMax = {0.65,0.067};
pLine=Line[{{0.69,0},pMax}];
cLine=Line[{{0.67,0},cMax}];
tLine=Line[{{0.72,0},tMax}];
fiveCLine=Line[{{0.65,0},fiveCMax}];
g1 = Graphics[{Black, Circle[{3,0}], Orange, Circle[{6.5,0}], Black, Circle[{10,0}], Arrowheads[{-0.03,0.03}], Thickness[0.005], Arrow[{{4,0},{5.5,0}}], Arrow[{{7.5,0},{9,0}}], (*Text[Style[\[Alpha],FontSize->Medium],{3,0}],*) Text[Style[\[Alpha],FontSize->Medium],{6.5,0}]}];
g2 = Graphics[{Black, Circle[{3,0}], Circle[{6.5,0}], Brown, Circle[{10,0}], Black, Circle[{13.5,0}], Circle[{17,0}], Arrowheads[{-0.03,0.03}], Thickness[0.005], Arrow[{{4,0},{5.5,0}}], Arrow[{{7.5,0},{9,0}}],Arrow[{{11,0},{12.5,0}}], Arrow[{{14.5,0},{16,0}}], (*Text[Style[\[Alpha],FontSize->Medium],{3,0}],*) Text[Style[\[Beta],FontSize->Medium],{10,0}]}];
g3 = Graphics[{Thickness[0.01],Blue, Circle[{3,0}], Red, Circle[{6.5,0}], Green, Circle[{10,0}], Black, Arrowheads[{-0.03,0.03}], Thickness[0.005], Arrow[{{4,0},{5.5,0}}], Arrow[{{7.5,0},{9,0}}]}];

(*Produce and export plot*)
p1 = ListPlot[{resc, resp2)}, Joined->True, PlotRange->{{0.5,1},{0,0.075}}, Frame->{True,True,False,False}, FrameLabel->{"Fitness (s)","Maximum Tolerated Migration (\!\(\*SuperscriptBox[\(m\), \(*\)]\))"},PlotLabel->"Linear 3 Demes",PlotStyle->{Orange, Red}, PlotLegends->{"Center of 3", "Peripheral Targeting")},Epilog->{(*{Directive[Brown, Dashed], fiveCLine},*){Directive[Orange, Dashed], cLine},Inset[g1, {0.712,0.067},{0,0},{.18,.36}]},AspectRatio->1];
UseFrontEnd[Export["threeDemes.png", p1]];
p1
Export["Dropbox/threeDemes.png", p1];
