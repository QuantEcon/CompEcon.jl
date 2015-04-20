function demopt02

%% DEMOPT02 Changes in Nelder-Mead simplex

% Preliminary tasks
demosetup(mfilename)

color0  = [.7 .7 .7];
color1  = [.2 .2 .2];
rvector = ones(3,1);
rfactor = 2;

S  = [[0.5;0.5] [.65;1.75] [1.5;0.35]];
x1 = S*rvector-rfactor*S(:,1);        % reflection
S1 = S; S1(:,1) = x1;
x2 = 1.5*x1-0.5*S(:,1);              	% expansion
S2 = S; S2(:,1) = x2;
x3 = 0.75*S(:,1)+0.25*x1;            	% contraction
S3 = S; S3(:,1) = x3;
S4 = S/2+S(:,2)*(ones(1,3)/2);        % shrinkage

figure
subplot(2,2,1)
hold on
patch(S(1,:),S(2,:),color0)
patch(S1(1,:),S1(2,:),color1)
layout('Reflection')

subplot(2,2,2)
hold on
patch(S(1,:),S(2,:),color0)
patch(S2(1,:),S2(2,:),color1)
layout('Expansion')

subplot(2,2,3)
hold on
patch(S(1,:),S(2,:),color0)
patch(S3(1,:),S3(2,:),color1)
layout('Contraction')

subplot(2,2,4)
hold on
patch(S(1,:),S(2,:),color0)
patch(S4(1,:),S4(2,:),color1)
layout('Shrinkage')

% Save Plots as EPS Files
printfigures(mfilename,1,1,0)


function layout(name)
axis([0.5 2.225 0.35 2.15])
axis off
text(0.35,0.4,'A')
text(0.55,1.9,'B')
text(1.50,0.2,'C')
title(name)