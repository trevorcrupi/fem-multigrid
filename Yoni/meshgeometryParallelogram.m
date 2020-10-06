%Make domain for mesh in this funciton

function [gd sf ns]=meshgeometryParallelogram(w)

%pdepoly([-0.5,0.5,0.5,0,-0.5],[-0.5,-0.5,0,0.5,0.5]);

%Use pdepoly in Matlab to get geometry. Then
%copy-paste into this file the three variables that you get when you do "export geometry
%data" in the drawing tool box.


%Domain: [0,1] times [0,1]

gd =[2;
     4;
     0;
     1;
     1;
     0;
     0;
     0;
     1;
     1];

sf = 'P1';

ns = [80; 49];


%[0,4] times [-4 4] for axisymmetric example
%{
gd = [2     2;
      4     4;
     0     0;
     4     4;
     4     4;
     0     0;
     0     0;
     0     0;
     4    -4;
     4    -4];


ns = [80    80;
      49    50];


sf = 'P1+P2';

%}

%Pentagon on my first draft of paper.
%{
gd=[2.0000; 
    5.0000;
   -0.5000;
    0.5000;
    0.5000;
         0;
   -0.5000;
   -0.5000;
   -0.5000;
         0;
    0.5000;
    0.5000];

sf='P1';
ns=[80;49];
%}

%Uniform mesh of pentagon.
%{
gd = [2.0000    2.0000    2.0000    2.0000;
    4.0000    4.0000    3.0000    4.0000;
   -0.5000   -0.5000         0         0;
         0         0    0.5000         0;
         0         0         0    0.5000;
   -0.5000   -0.5000         0    0.5000;
   -0.5000         0         0   -0.5000;
   -0.5000         0    0.5000         0;
         0    0.5000         0         0;
         0    0.5000         0   -0.5000];



ns = [ 80    80    80    80;
       49    50    51    52];


sf = 'P1+P2+P3+P4';
%}



%%unit square containing the origin.
%{
gd =[2     2     2     2;
     4     4     4     4;
     1     0     0     1;
     1     0     0     1;
     0    -1    -1     0;
     0    -1    -1     0;
     0     0    -1    -1;
     1     1     0     0;
     1     1     0     0;
     0     0    -1    -1];

sf='P1+P2+P3+P4';

ns=[80    80    80    80;
    49    50    51    52];
%}

%unit Square whose left bottom vertex is the origin.
%{
gd = [2.0000    2.0000    2.0000    2.0000;
    4.0000    4.0000    4.0000    4.0000;
         0    0.5000    0.5000         0;
    0.5000    1.0000    1.0000    0.5000;
    0.5000    1.0000    1.0000    0.5000;
         0    0.5000    0.5000         0;
         0         0    0.5000    0.5000;
         0         0    0.5000    0.5000;
    0.5000    0.5000    1.0000    1.0000;
    0.5000    0.5000    1.0000    1.0000];

sf = 'P1+P2+P3+P4';

ns = [80    80    80    80;
      49    50    51    52];

%}

%Lshape with bottom left corner (0 0) uniform mesh.
%{
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
%}

%Lshape mesh quasi uniform, not uniform
%{
gd = [2.0000;
    6.0000;
         0;
    1.0000;
    1.0000;
    0.5000;
    0.5000;
         0;
         0;
         0;
    0.5000;
    0.5000;
    1.0000;
    1.0000];


ns =[ 80;
    49];



sf = 'P1';
%}

%Standard Unit Square Standard Uniform Mesh
%{
gd = [2     2;
     3     3;
    -1    -1;
     1     1;
     1    -1;
    -1    -1;
    -1     1;
     1     1];


sf = 'P1+P2';

ns = [80    80;
      49    50];
%}






%[-4 4]\times [-4 4]
%{
gd =[2;
     4;
    -4;
     4;
     4;
    -4;
    -4;
    -4;
     4;
     4];

sf = 'P1';


ns=[80; 49];
%}

% pdepoly([0,1,1],[0,0,1]) + pdepoly([0,1,0],[0,0,-1])
gd = [2     2;
     3     3;
     0     0;
     1     1;
     1     0;
     0     0;
     0     0;
     1    -1];




ns = [80    80;
      49    50];
    
sf = 'P1+P2';



end
