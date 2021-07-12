(*
 Program that models the likelihood of a gene drive spillover from a targeted population to a non target population
In this case, we examine the case of five populations arranged in a line, with symmetrical migration between populations
Produces a graph that shows the relationship between fitness cost (s) and maximum tolerated migration rate {m^*}
*)

$HistoryLength = 0;

(*Given values of fitness cost, CRISPR conversion rate, gene dominance, and migration rate, solves the system of equations of allele frequency in each population and returns a table showing all equilibrium solutions*)
eqpoints1[s_,c_,h_,m_]:=Module[{eeq1,eeq2,eeq3,eeq4,eeq5,sn,sc},(*Gives the equilibrium solutions for s,c,h,m*)
sn=1/2 (1-c)(1-h s);
sc=c(1-s);
eeq1=q1==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q1+m*q2)};
eeq2=q2==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2m)q2+m*q1+m*q3)};
eeq3=q3==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2m)q3+m*q2+m*q4)};
eeq4=q4==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2m)q4+m*q3+m*q5)};
eeq5=q5==(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q5+m*q4)};
sols=NSolve[{eeq1,eeq2,eeq3,eeq4,eeq5},{q1,q2,q3,q4,q5},Reals,Method->"Newton", WorkingPrecision->30];

z1=Table[
If[
Length[sols]==0||Length[ss]!=5,"N",
If[ss[[1,2]]<0-10^-5||ss[[2,2]]<0-10^-5||ss[[3,2]]<0-10^-5||ss[[4,2]]<0-10^-5||ss[[5,2]]<0-10^-5 ||ss[[1,2]]>1.+10^-5||ss[[2,2]]>1.+10^-5||ss[[3,2]]>1.+10^-5 ||ss[[4,2]]>1.+10^-5||ss[[5,2]]>1.+10^-5,"N",{ss[[1,2]],ss[[2,2]],ss[[3,2]],ss[[4,2]],ss[[5,2]]}]
](*Here we remove all solution that are not in the interval [0,1] with some small error margin due to numeric accuracy*)
,{ss,sols}];
z2=DeleteCases[z1,"N"];
z3=DeleteDuplicates[Round[z2,10^-5]];(*Due to numeric accuracy, some solutions are counted twice sometimes. Here we remove duplicate solutions*)
N[Sort[z3]]
]

(*System of equations representing allele frequency in each population*)
eq1[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q1+m*q2)};
eq2[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2*m)q2+m*q1+m*q3)};
eq3[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2*m)q3+m*q2+m*q4)};
eq4[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-2*m)q4+m*q3+m*q5)};
eq5[m_,s_,sc_,sn_]:=(q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q5+m*q5)};

(*Jacobian matrix*)
J[m_,s_,sc_,sn_] := {{D[eq1[m,s,sc,sn], q1],D[eq1[m,s,sc,sn], q2],D[eq1[m,s,sc,sn], q3],D[eq1[m,s,sc,sn], q4],D[eq1[m,s,sc,sn], q5]}, {D[eq2[m,s,sc,sn], q1],D[eq2[m,s,sc,sn], q2],D[eq2[m,s,sc,sn], q3],D[eq2[m,s,sc,sn], q4],D[eq2[m,s,sc,sn], q5]},{D[eq3[m,s,sc,sn], q1],D[eq3[m,s,sc,sn], q2],D[eq3[m,s,sc,sn], q3],D[eq3[m,s,sc,sn], q4],D[eq3[m,s,sc,sn], q5]},{D[eq4[m,s,sc,sn], q1],D[eq4[m,s,sc,sn], q2],D[eq4[m,s,sc,sn], q3],D[eq4[m,s,sc,sn], q4],D[eq4[m,s,sc,sn], q5]},{D[eq5[m,s,sc,sn], q1],D[eq5[m,s,sc,sn], q2],D[eq5[m,s,sc,sn], q3],D[eq5[m,s,sc,sn], q4],D[eq5[m,s,sc,sn], q5]}};

(*Given values of fitness cost, CRISPR conversion rate, gene dominance, and migration rate, selects the equilibria which represent Differential Targeting (gene drive allele frequency above 50% in the target population and below 50% in all non-target populations)
Then, calculates stability of equilibria by finding eigenvalues of the Jacobian matrix
Returns list of stable DTEs, for central targeting case, peripheral targeting, and 2-peripheral targeting cases*)
DTE[c_,h_,s_,m_]:=Module[{},
sc=c*(1-s); sn=(1-c)*(1-h*s);
nn=eqpoints1[s,c,h,m];
nna=Select[nn,#[[1]]>0.5&&#[[2]]<0.5&&#[[3]]<0.5&&#[[4]]<0.5&&#[[5]]<0.5&];(*potential DTE for target pop peripheral*)
nnb=Select[nn,#[[1]]<0.5&&#[[2]]>0.5&&#[[3]]<0.5&&#[[4]]<0.5&&#[[5]]<0.5&];(*potential DTE for target pop 2-peripheral*)
nnc=Select[nn,#[[1]]<0.5&&#[[2]]<0.5&&#[[3]]>0.5&&#[[4]]<0.5&&#[[5]]<0.5&];(*potential DTE for target pop central*)

staa=If[nna!={},Table[J[m,s,sc,sn]/.{q1->nna[[i,1]],q2->nna[[i,2]],q3->nna[[i,3]],q4->nna[[i,4]],q5->nna[[i,5]]},{i,Length[nna]}],"None1"];
sta=If[!StringQ[staa],Table[Eigenvalues[SetPrecision[staa[[j]],MachinePrecision],1,Method->"Arnoldi"],{j,Length[staa]}],"None1"];
stbb=If[nnb!={},Table[J[m,s,sc,sn]/.{q1->nnb[[i,1]],q2->nnb[[i,2]],q3->nnb[[i,3]],q4->nnb[[i,4]],q5->nnb[[i,5]]},{i,Length[nnb]}],"None1"];
stb=If[!StringQ[stbb],Table[Eigenvalues[SetPrecision[stbb[[j]],MachinePrecision],1,Method->"Arnoldi"],{j,Length[stbb]}],"None1"];
stcc=If[nnc!={},Table[J[m,s,sc,sn]/.{q1->nnc[[i,1]],q2->nnc[[i,2]],q3->nnc[[i,3]],q4->nnc[[i,4]],q5->nnc[[i,5]]},{i,Length[nnc]}],"None1"];
stc=If[!StringQ[stcc],Table[Eigenvalues[SetPrecision[stcc[[j]],MachinePrecision],1,Method->"Arnoldi"],{j,Length[stcc]}],"None1"];

DTEsp=If[!StringQ[sta],DeleteCases[Table[If[Abs[sta[[i,1]]]<1(*&&Abs[sta[[i,2]]]<1&&Abs[sta[[i,3]]]<1&&Abs[sta[[i,4]]]<1&&Abs[sta[[i,5]]]<1*),nna[[i]],"N"],{i,Length[sta]}],"N"],sta];(*DTEs for target pop peripheral*)
DTEst=If[!StringQ[stb],DeleteCases[Table[If[Abs[stb[[i,1]]]<1(*&&Abs[stb[[i,2]]]<1&&Abs[stb[[i,3]]]<1&&Abs[stb[[i,4]]]<1&&Abs[stb[[i,5]]]<1*),nnb[[i]],"N"],{i,Length[stb]}],"N"],stb];(*DTEs for target pop 2-peripheral*)
DTEsc=If[!StringQ[stc],DeleteCases[Table[If[Abs[stc[[i,1]]]<1(*&&Abs[stc[[i,2]]]<1&&Abs[stc[[i,3]]]<1&&Abs[stc[[i,4]]]<1&&Abs[stc[[i,5]]]<1*),nnc[[i]],"N"],{i,Length[stc]}],"N"],stc];(*DTEs for target pop central*)

If[StringQ[DTEsp],DTEsp="None2"];
If[StringQ[DTEst],DTEst="None2"];
If[StringQ[DTEsc],DTEsc="None2"];

{DTEsp,DTEst,DTEsc}
]

(*Given a value for fitness cost, find the maximum tolerated migration rate for both central targeting case and peripheral targeting case*)
mqstar[ss_]:=Module[{maxm,jumps,t0,t1,t2,t3,outm}, (*For a given c, h, and s, returns m^* and Subscript[q, 2]^* *)
Print[ss];
cc = 1;
hh = 0;
maxm=0.14;(*Maximal m for searching m^* *)
jumps=10^-3;(*Resolution in m for searching m^*  *)


t0 = ParallelTable[{N[mm],DTE[cc,hh,ss,mm]},{mm,0,maxm,jumps}];

t1=Table[{t[[1]],t[[2]]},{t,t0}];
t1p=Table[{t[[1]],Length[t[[2,1]]]},{t,t1}];
t1t=Table[{t[[1]],Length[t[[2,2]]]},{t,t1}];
t1c=Table[{t[[1]],Length[t[[2,3]]]},{t,t1}];
mstarp=If[Last[t1p][[2]]==0,Max[{Min[Flatten[Take[Select[t1p,#[[2]]==0&],All,{1}]]]-jumps,0}],"Non"];
mstart=If[Last[t1t][[2]]==0,Max[{Min[Flatten[Take[Select[t1t,#[[2]]==0&],All,{1}]]]-jumps,0}],"Non"];
mstarc=If[Last[t1c][[2]]==0,Max[{Min[Flatten[Take[Select[t1c,#[[2]]==0&],All,{1}]]]-jumps,0}],"Non"];
{mstarp,mstart,mstarc}
]

c=1;h=0; time=100; th=0.01;
dss = 0.005;

SetDirectory["Yonatan"];
Directory[]
DateString[]

(*Export results*)
res06 = Table[{ss,mqstar[ss]}, {ss, 0.5, 0.6, dss}];
Export["fiveRes06HR2.mx", res06];
Export["Dropbox/fiveRes06HR2.mx", res06];
res07 = Table[{ss,mqstar[ss]}, {ss, 0.6+dss, 0.7, dss}];
Export["fiveRes07HR2.mx", res07];
Export["Dropbox/fiveRes07HR2.mx", res07];
res08 = Table[{ss,mqstar[ss]}, {ss, 0.7+dss, 0.8, dss}];
Export["fiveRes08HR2.mx", res08];
Export["Dropbox/fiveRes08HR2.mx", res08];
res09 = Table[{ss,mqstar[ss]}, {ss, 0.8+dss, 0.9, dss}];
Export["fiveRes09HR2.mx", res09];
Export["Dropbox/fiveRes09HR2.mx", res09];
res10 = Table[{ss,mqstar[ss]}, {ss, 0.9+dss, 1, dss}];
Export["fiveRes10HR2.mx", res10];
Export["Dropbox/fiveRes10HR2.mx", res10];


(*Import results and produce graph*)
res62 = Import["fiveRes06HR2.mx"];
res72 = Import["fiveRes07HR2.mx"];
res82 = Import["fiveRes08HR2.mx"];
res92 = Import["fiveRes09HR2.mx"];
res102 = Import["fiveRes10HR2.mx"];

res2 = Join[res62, res72, res82, res92, res102];

resp2 = Table[{r[[1]],r[[2,1]]}, {r,res2}];
rest2 = Table[{r[[1]],r[[2,2]]}, {r,res2}];
resc2 = Table[{r[[1]],r[[2,3]]}, {r,res2}];

resp2
rest2
resc2
UseFrontEnd[Export["fiveResPFinal.mx", resp2]];
Export["Dropbox/fiveResPFinal.mx", resp2];
UseFrontEnd[Export["fiveResTFinal.mx", rest2]];
Export["Dropbox/fiveResTFinal.mx", rest2];
UseFrontEnd[Export["fiveResCFinal.mx", resc2]];
Export["Dropbox/fiveResCFinal.mx", resc2];

pMax = {0.69,0.13};
tMax = {0.665,0.061};
cMax = {0.65,0.067};

pLine=Line[{{0.69,0},pMax}];
tLine=Line[{{0.665,0},tMax}];
cLine=Line[{{0.65,0},cMax}];

g1 = Graphics[{Black, Circle[{3,0}],Blue, Circle[{6.5,0}], Black, Circle[{10,0}], Black, Circle[{13.5,0}], Circle[{17,0}], Arrowheads[{-0.03,0.03}], Thickness[0.005], Arrow[{{4,0},{5.5,0}}], Arrow[{{7.5,0},{9,0}}], Arrow[{{11,0},{12.5,0}}], Arrow[{{14.5,0},{16,0}}], (*Text[Style[\[Alpha],FontSize->Medium],{3,0}], Text[Style[\[Beta],FontSize->Medium],{6.5,0}],*) Text[Style[\[Gamma], FontSize->Medium],{6.5,0}], Arrowheads[0.05]}];
g2 = Graphics[{Black, Circle[{3,0}], Blue, Circle[{6.5,0}], Orange, Circle[{10,0}], Black, Circle[{13.5,0}], Circle[{17,0}], Arrowheads[{-0.03,0.03}], Thickness[0.005], Arrow[{{4,0},{5.5,0}}], Arrow[{{7.5,0},{9,0}}], Arrow[{{11,0},{12.5,0}}], Arrow[{{14.5,0},{16,0}}], Text[Style[\[Alpha],FontSize->Medium],{6.5,0}], Text[Style[\[Beta], FontSize->Medium],{10,0}], Arrowheads[0.05]}];
g3 = Graphics[{Black, Circle[{3,0}], Blue, Circle[{6.5,0}], Black, Circle[{10,0}],(* Black, Circle[{13.5,0}], Circle[{17,0}], Circle[{20.5,0}],*) Arrowheads[{-0.03,0.03}], Thickness[0.005], Arrow[{{4,0},{5.5,0}}], Arrow[{{7.5,0},{9,0}}], (*Arrow[{{11,0},{12.5,0}}], Arrow[{{14.5,0},{16,0}}], Arrow[{{18,0},{19.5,0}}],*) Text[Style[\[Alpha],FontSize->Medium],{6.5,0}], Arrowheads[0.05]}]
UseFrontEnd[Export["threeDemeMap2.png", g3]];
Export["Dropbox/threeDemeMap2.png", g3];

p1 = ListPlot[{resp2,rest2,resc2}, Joined->True, Frame->{True,True,False,False}, FrameLabel->{"Fitness (s)","Maximum Tolerated Migration (m^*)"},PlotLabel->"Maximum Tolerated Migration In Linear 5-Deme Network", PlotRange->{{0.5,1},{0,0.155}},PlotStyle->{Blue,Orange,Brown}, PlotLegends->{"Targeting \[Alpha]","Targeting \[Beta]","Targeting \[Gamma]"},Epilog -> {{Directive[Blue, Dashed], pLine},{Directive[Orange, Dashed], tLine},{Directive[Brown, Dashed], cLine},Text[pMax, Offset[{0, 10}, pMax]],Text[cMax, Offset[{-27, 10}, cMax]],Text[tMax, Offset[{30, 10}, tMax], Background->White], Inset[g2, {0.64,0.145},{0,0},{.33,.66}]}, AspectRatio->1];
p1
