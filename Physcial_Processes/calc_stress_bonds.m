function [floe] = calc_stress_bonds(floe,Nfrac,dt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a = floe.bonds.interactions;
global Modulus
if ~isempty(a)
    
    %% Initialize stress calc by creating floes out of subfloes
    RIDGING=false;
    
    PERIODIC=false;
    
    COLLISION = true;
        
    RAFTING = false;
    
    doInt.flag = false;
    doInt.step = 10;
    rmax = floe.rmax;
    
    ocean.Xo = [-10 10]*rmax; ocean.Yo= [-10 10]*rmax; ocean.Uocn= 0; ocean.Vocn = 0;
    ocean.fCoriolis=1.4e-4; ocean.r = 1/(2.5*24*3600); 
    winds.u = 0; winds.v = 0;
    HFo = 0; min_floe_size = 0; Nb = 0; Nx = 0; Ny = 0;dissolvedNEW = 0;
    
    forces = a(:,3:4)/Nfrac;
    c2_boundary = floe.c0;
    c2_boundary_poly = polyshape(c2_boundary');
    height.mean = floe.h; height.delta = 0;
    c2_border = polyshape(2*rmax*[-1 -1 1 1; -1 1 1 -1]'); c2_border = subtract(c2_border, c2_boundary_poly);
    floebound = initialize_floe_values(c2_border, height,1);
    Floe = []; Xk = zeros(length(floe.bonds.Vert),1); Yk = zeros(length(floe.bonds.Vert),1);
    for kk = 1:length(floe.bonds.Vert)
        poly = polyshape(floe.bonds.Vert{kk});
        floenew = initialize_floe_values(poly, height,1);
        Xk(kk) = floe.bonds.Xs(kk);
        Yk(kk) = floe.bonds.Ys(kk);
        in = inpolygon(a(:,1),a(:,2),poly.Vertices(:,1),poly.Vertices(:,2));
        floenew.FxOA =  sum(forces(in,1))/floenew.area; %force by interactions with other floes
        floenew.FyOA =  sum(forces(in,2))/floenew.area; %force by interactions with other floes      
        Floe = [Floe floenew];
    end
    
    Xi = cat(1,Floe.Xi); Yi = cat(1,Floe.Yi);
    
    %% Run for desired number of steps to move subfloes about
    for kk = 1:50
        [Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, 0, 0, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
    end

    %%Treat as elastic by finding displacement from Xk, Yk
    Xi_new = cat(1,Floe.Xi); Yi_new = cat(1,Floe.Yi);
    for kk = 1:length(Floe)        
        A_rot=[cos(Floe(kk).alpha_i) -sin(Floe(kk).alpha_i); sin(Floe(kk).alpha_i) cos(Floe(kk).alpha_i)]; %rotation matrix
        Voronoi_new = A_rot*[Xi_new(kk);Yi_new(kk)];
        dx = Voronoi_new(1)-Xi(kk); dy = Voronoi_new(2)-Yi(kk);
        Xk_new(kk) = Xk(kk)+dx;
        Yk_new(kk) = Yk(kk)+dy;
    end

    q = 5.2; SigC = 250e3;
    Sig1 = (1/q+1)*SigC/(1/q-q);
    Sig2 = q*Sig1+SigC;
    Sig11 = 4e4;
    Sig22 = q*Sig11+SigC;
    MohrX = [Sig1; Sig11; Sig22];
    MohrY = [Sig2; Sig22; Sig11];
    Mohr = polyshape(-MohrX,-MohrY);
    thet = 45*pi/180;
    A_rot2 = [cos(thet) -sin(thet); sin(thet) cos(thet)]; %rotation matrix
    for kk = 1:length(Floe)
        A_rot=[cos(Floe(kk).alpha_i) -sin(Floe(kk).alpha_i); sin(Floe(kk).alpha_i) cos(Floe(kk).alpha_i)]; %rotation matrix
        Voronoi_new = A_rot*[Xi_new(kk);Yi_new(kk)];
        dx(kk) = Voronoi_new(1)-Xi(kk); dy(kk) = Voronoi_new(2)-Yi(kk);
        bonds =floe.bonds.Num{kk};
        L = floe.bonds.d{kk};
        clear Princ1; clear Princ2;
        for jj = 1:length(bonds)
            theta = atan((Yk(bonds(jj))-Yk(kk))/(Xk(bonds(jj))-Xk(kk)));
            theta_new = atan((Yk_new(bonds(jj))-Yk_new(kk))/(Xk_new(bonds(jj))-Xk_new(kk)));
            A_rot=[cos(theta_new-theta) -sin(theta_new-theta); sin(theta_new-theta) cos(theta_new-theta)]; %rotation matrix
            d0 = sqrt((Xk(kk)-Xk(bonds(jj))).^2+(Yk(kk)-Yk(bonds(jj))).^2);
            d1 = sqrt((Xk_new(kk)-Xk_new(bonds(jj))).^2+(Yk_new(kk)-Yk_new(bonds(jj))).^2);
            h1 = Floe(kk).h; h2 = Floe(bonds(jj)).h;
            r1 = sqrt(Floe(kk).area); r2 = sqrt(Floe(bonds(jj)).area);
            Force_factor=Modulus*(h1*h2)/(h1*r2+h2*r1);
            Force = Force_factor*(d0-d1);
            Stress = A_rot*[Force; 0]/(floe.h*L(jj));%use elastic model to find force between subfloes that are bonded together...
            Princ = eig(A_rot2*[Stress(1) Stress(2); Stress(2) 0]);
            Princ1(jj) = max(Princ); Princ2(jj) = min(Princ);
        end
        [inpoly,~] = inpolygon(Princ1,Princ2,Mohr.Vertices(:,1), Mohr.Vertices(:,2));
        if min(inpoly) == 0
            xx = 1; xx(1) =[1 2];
        end
    end
end
floe.bonds.interactions = [];
end