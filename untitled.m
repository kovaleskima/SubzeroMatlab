% computes M
load(‘buoy-ERA5-2011.mat’);
time = serialTime; % units of days, 1 hourly increments
longitude(longitude < 0) = 360 + longitude(longitude < 0);
[time, indexUnique] = unique(time);
latitude = latitude(indexUnique);
longitude = longitude(indexUnique);
ice = icebuoy(indexUnique);
uwind = uwind(indexUnique);
vwind = vwind(indexUnique);
dt = 1/24;
timeGrid = time(1):dt:time(end);
latitudeGrid = interp1(time,latitude,timeGrid);
longitudeGrid = interp1(time,longitude,timeGrid);
iceGrid = interp1(time,ice,timeGrid);
uwindGrid = interp1(time,uwind,timeGrid);
vwindGrid = interp1(time,vwind,timeGrid);
Uwind = uwindGrid +1i*vwindGrid; % store wind velocity as complex number
%% ------------------------
% Now compute velocities:
DEG2RAD = pi/180;
OMEGA = 7.292e-5*((60*60*24)/(2*pi)); % cycles/day
f0 = 2*OMEGA*sin(latitudeGrid*DEG2RAD);
% stereographic projection
x = 110.949*(90-latitudeGrid).*cos(longitudeGrid.*DEG2RAD); y = 110.949*(90-latitudeGrid).*sin(longitudeGrid.*DEG2RAD); % km
% velocity
u = diff(x)./dt; v = diff(y)./dt; % km/d
x = x(1:end-1); y = y(1:end-1); iceGrid = iceGrid(1:end-1);
timeGrid = timeGrid(1:end-1);
U = u +1i*v; % store velocity as complex number
dx = gc_dist(latitudeGrid(1:end-1),longitudeGrid(1:end-1),latitudeGrid(2:end),longitudeGrid(2:end));
speed = dx/dt; %km/day
speedms = speed*1000/(24*60*60); %m/s
% Now compute M:
% Cf. Torrence and Compo (1998)
N = length(timeGrid);
dj = 1/12; % default from ‘wt’, gives good resolution
w0 = 6; % non-dimensional frequency, value 6 to satisfy ‘admissibility condition’, Terrence and Compo ’98
s0 = 2*dt; % base scale
MaxScale = N*0.17*2*dt; % after wt
J = round((log2(MaxScale/s0))/dj);
j = 0:1:J;
s = s0* 2.^(j*dj); % scales, in units of ‘hours’, this is important
fT = (4*pi*s)/(w0 + sqrt(2+w0^2)); % associated fourier periods
fF = 1./fT; % and fourier frequencies
%% WTIF
sig.val = U; % signal structure, input to cwtft wavelet transform
sig.period = dt;
cwtstruct = cwtft(sig,‘scales’,s); % continuous wavelet transform
coef = cwtstruct.cfs; % coefficients from transform
% find closest fourier frequency to f_0, take that coefficient
inCoef = NaN(1,N);
for j = 1:N
    fInd = find(abs(f0(j)-fF)==min(abs(f0(j)-fF)));
    inCoef(j)=coef(fInd,j);% coefficient from transform
end
C = abs(inCoef); % Magnitude of coefficient
Ubar = smooth(abs(U),3/dt)’; % Magnitude of velocity, smoothed at 3 days
M = C./Ubar; % M ratio