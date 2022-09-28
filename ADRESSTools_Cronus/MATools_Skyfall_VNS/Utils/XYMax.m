function [x,y]=XYMax(PTS)
% [x,y]=xymax(PTS) quadraticly interpolates a 3-point vector to find the maximum.
% The points are assumed located at [0 1 2]. The maximum may be outside. 
% Based on QUAD2 from the Optimization toolbox.

% input check
if size(PTS,1)>size(PTS,2); pts=PTS; else pts=PTS'; end
    
c=pts(1);
ab=[-1 0.5; 2 -0.5]*(pts(2:3)-c*ones(2,1));
x=-ab(2)/(2*ab(1));
y=ab(1)*x^2+ab(2)*x+c;


