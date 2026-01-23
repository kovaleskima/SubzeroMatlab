function c2_boundary =initialize_boundaries()

%Boundaries are defined as a set of contours e.g. [x1 NaN x2 ; y1 NaN y2]

%Adding walls around the domain
Lx=5e4; Ly=5e4;
x=[-1 -1 1 1 -1]*Lx; 
y=[-1 1 1 -1 -1]*Ly;

c2_boundary = [x; y];

end