function [ocean, heat_flux, h0]=initialize_ocean(dt,nDTOut)

% defining ocean currents
ocean.fCoriolis=1.45e-4; % Coriolis parameter.
ocean.r = 1/(2.5*24*3600); % Damping parameter
ocean.H = 15; %ocean mix layer depth

ocean.U = 0.5; ocean.V = ocean.U;

ocean.turn_angle=-15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% ocean grid;
dXo=500; % in meters

Xo=-3e5:dXo:3e5; Yo=-3e5:dXo:3e5; 
[Xocn, Yocn]=meshgrid(Xo,Yo);
ocean.Xocn = Xocn; ocean.Yocn = Yocn;

%defining ocean streamfunction with some eddies
transport=5e3; % horizontal transport, in m^2/s (controls ocean currents) 
psi_ocean=transport*sin(2*pi*Xocn/40e4).*cos(2*pi*Yocn/50e4); 

%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=zeros(size(Xocn));

ocean.Xo=Xo;
ocean.Yo=Yo;
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;
ocean.Uocn_p=Uocn;
ocean.Vocn_p=Vocn;

Tice = -20; Tocean = 2*ones(size(Xocn));
heat_flux = 7.4*10^(-4)*(Tice-Tocean)/(72); %cm^2/s
heat_flux = heat_flux/100^2; %m^2/s
h0 = real(sqrt(-2*dt*heat_flux*nDTOut));
h0 = mean(h0(:));

ocean.TauAtmX_p=zeros(size(Xocn));
ocean.TauIceX_p=zeros(size(Xocn));
ocean.TauAtmY_p=zeros(size(Xocn));
ocean.TauIceY_p=zeros(size(Xocn));

end
