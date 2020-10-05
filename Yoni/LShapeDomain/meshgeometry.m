%Make domain for mesh in this funciton

function [gd sf ns]=meshgeometry(w)

%pdepoly([-0.5,0.5,0.5,0,-0.5],[-0.5,-0.5,0,0.5,0.5]);

%Current: Pentagon. First use pdepoly in Matlab to get geometry. Then
%copy-paste in this file the three variables that you get when you do "export geometry
%data" in the drawing tool box.



%Lshape with bottom left corner (0 0) uniform mesh.

gd = [2.0000    2.0000    2.0000;
    4.0000    4.0000    4.0000;
         0    0.5000         0;
    0.5000    1.0000    0.5000;
    0.5000    1.0000    0.5000;
         0    0.5000         0;
         0         0    0.5000;
         0         0    0.5000;
    0.5000    0.5000    1.0000;
    0.5000    0.5000    1.0000];


sf = 'P1+P2+P3';

ns = [80    80    80;
      49    50    51];


end
