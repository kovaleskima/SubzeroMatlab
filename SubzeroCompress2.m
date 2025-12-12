close all; clear all;

%% Set Flags

RIDGING=false; 

FRACTURES=true;

PERIODIC=false;

PACKING = false;

WELDING = false;

CORNERS = false;

COLLISION = true;

AVERAGE = false;

RAFTING = false;

KEEP_MIN = false;

ifPlot = false; %Plot floe figures or not?

ifPlotStress = false;

ifPlotStressStrain = false;

%% Initialize model vars

dt=10; %Time step in sec
height.mean = 0.25;
height.delta = 0;
% h0 = 0.1; %thickness of ice that gets packed in

%Define ocean currents
nDTpack = 4500;
rho_ice=920;
[ocean, HFo, h0]=initialize_ocean(dt,nDTpack);

%Define 10m winds
winds=[0 0];

%Define boundaries
%c2_boundary=initialize_boundaries();
% Lx=3.5e4; Ly=3.5e4; f=2*pi/(24*3600);
Lx=1e5; Ly=1e5; f=2*pi/(24*3600);
% Lx=5e4; Ly=5e4; f=2*pi/(24*3600);
x=[-1 -1 1 1 -1]*Lx; 
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
uright = 0; uleft = 0;
min_floe_size = 4*Lx*Ly/10000;

%Initialize Floe state
target_concentration = 1;%[1;1;0];
[Floe, bonds, Nbond, Nbound] = initial_concentration(c2_boundary,target_concentration,height,100,1,min_floe_size);
x=[-1 -1 1 1 -1]*Lx; 
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
c2_boundary_poly = polyshape(c2_boundary');
c2_border = polyshape(2*[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly]'); c2_border = subtract(c2_border, c2_boundary_poly);
floebound = initialize_floe_values(c2_border, height, 1);
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
% load floetest
min_floe_size = (4*Lx*Ly-sum(cat(1,Floe(1:Nbound).area)))/20000;
global Modulus
%load Modulus
collide = 0;
Modulus = 1.5e3*(mean(sqrt(cat(1,Floe.area)))+min(sqrt(cat(1,Floe.area))));
save('Modulus.mat','Modulus');
fdir.dir =[0;0];fdir.loc =[0 0];
save('F_dir.mat','fdir')
frac = 0;

%%

dhdt = 1;

nDTOut=40; %Output frequency (in number of time steps)

nSnapshots=20; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nSimp = 20;

target_concentration=1;
tStart = tic; 
doInt.flag = true;
doInt.step = 10;

% specify coarse grid size
LxO= 2*max(ocean.Xo);LyO= 2*max(ocean.Yo);
Nx=10; Ny=10;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);

%Initiailize Eulearian Data
[eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nbound,c2_boundary,dt,PERIODIC);
Vd = zeros(Ny,Nx,2);
Vdnew=zeros(Ny, Nx);
SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
DSigX = 0; DSigY= 0; DSig1= 0; DSig2= 0;
DivSigX = zeros(Ny, Nx); DivSig1 = zeros(Ny, Nx);
DivSigY = zeros(Ny, Nx); DivSig2 = zeros(Ny, Nx);
Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
U = zeros(Ny, Nx); V = zeros(Ny, Nx);
dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
Fx = zeros(Ny, Nx); Fy = zeros(Ny, Nx);
Sig = zeros(Ny, Nx); mass = zeros(Ny,Nx);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
FloeOld = Floe;
[Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, uright, 0, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nbound, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);

%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end
if isempty(dir('FloesCom2')); disp('Creating folder: FloesCom2'); mkdir('FloesCom2'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=figure('Position',[100 100 1000 500],'visible','on');
    fig3 = figure;
end
%i_step = (im_num-1)*nDTOut;

for ii =1:length(Floe)
    poly(ii) = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end
%% Solving for floe trajectories
tic; 
%while Time<120000
while im_num<500
% while im_num<nSnapshots || area(intersect(poly(1),poly(2)))>0

    if mod(i_step,10)==0
        dXo=20000; transport=5e4; % horizontal transport, in m^2/s (controls ocean currents) 
        disp(' ');
        toc
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        disp(['number of collisions: ' num2str(numCollisions)]);
        disp(' ');
        tic
        doInt.flag=true;
    else
        doInt.flag=false;
    end


    if mod(i_step,nSimp)==0
    FloeOld = Floe;
    parfor j=1:length(FloeOld)
        FloeOld(j).poly = polyshape(FloeOld(j).c_alpha'+[FloeOld(j).Xi FloeOld(j).Yi]);
    end
        floenew = [];
        for ii = 1:length(Floe)
            ParFloes(ii).floenew = [];
            ParFloes(ii).kill = [];
            ParFloes(ii).verts = 0;
        end
        if isfield(Floe,'poly')
            Floe=rmfield(Floe,{'poly'});
        end
        parfor ii = 1:length(Floe)
            floe = Floe(ii);
            if length(Floe(ii).c0) > 30
                [floe2,kill] = FloeSimplify(Floe(ii),0,FloeOld);
                if isfield(floe2,'poly')
                    floe2=rmfield(floe2,{'poly'});
                end
                if isempty(kill)
                    ParFloes(ii).kill = [ParFloes(ii).kill kill];
                end
                for jj = 1:length(floe2)
                    if jj == 1
                        Floe(ii) = floe2(jj);
                        ParFloes(ii).verts=length(floe2(jj).c0(1,:));
                    else
                        ParFloes(ii).floenew = [ParFloes(ii).floenew floe2(jj)];
                    end
                end
            end
        end
        if isfield(Floe,'poly')
            Floe=rmfield(Floe,{'poly'});
        end
        floenew =[]; kill = [];
        for ii = 1:length(ParFloes)
            floenew = [floenew ParFloes(ii).floenew];
            kill = [kill ParFloes(ii).kill];
        end
        kill = unique(kill(kill>0));
        live = cat(1,Floe.alive);
        Floe(live==0)=[]; 
        Floe =[Floe floenew];
    end

    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        

        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nbound,Nbond,c2_boundary,dt,PERIODIC);
        if ifPlot
            [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nbound,PERIODIC);
        end
        

        if AVERAGE
            SigXXa = SigXX/fix(nDTOut); SigYXa = SigYX/fix(nDTOut);
            SigXYa = SigXY/fix(nDTOut); SigYYa = SigYY/fix(nDTOut);
            DivSigXa = DivSigX/fix(nDTOut); DivSig1a = DivSig1/fix(nDTOut);
            DivSigYa = DivSigY/fix(nDTOut); DivSig2a = DivSig2/fix(nDTOut);
            Eux = Eux/fix(nDTOut); Evx = Evx/fix(nDTOut);
            Euy = Euy/fix(nDTOut); Evy = Evy/fix(nDTOut);
            U = U/fix(nDTOut); V = V/fix(nDTOut);
            dU = dU/fix(nDTOut); dV = dV/fix(nDTOut);
            Fx = Fx/fix(nDTOut); Fy = Fy/fix(nDTOut);
            Sig = Sig/fix(nDTOut); 
            mass = mass/fix(nDTOut);
        else
            SigXXa = squeeze(eularian_data.stressxx); SigYXa = squeeze(eularian_data.stressyx);
            SigXYa = squeeze(eularian_data.stressxy); SigYYa = squeeze(eularian_data.stressyy);
            DivSigXa = DSigX; DivSig1a = DSig1;
            DivSigYa = DSigY; DivSig2a = DSig2;
            Eux = squeeze(eularian_data.strainux); Evx = squeeze(eularian_data.strainvx);
            Euy = squeeze(eularian_data.strainuy); Evy = squeeze(eularian_data.strainvy);
            U = U+squeeze(eularian_data.u);V = V+squeeze(eularian_data.v);
            dU = dU+squeeze(eularian_data.du);dV = dV+squeeze(eularian_data.dv);
            Fx = Fx+squeeze(eularian_data.force_x);Fy = Fy+squeeze(eularian_data.force_x);
            Sig = Sig+squeeze(eularian_data.stress);
        end
        
    end
    
    FloeOld2 = FloeOld; FloeOld = Floe;
    if PACKING && h0 > 0
        if mod(i_step,nDTpack)==0
            height.mean = 0.2;%h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice_new(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean, height, min_floe_size, PERIODIC,3,3);
        end
    end
    
    live = cat(1,Floe.alive);
    Floe(live==0)=[];     
    if mod(i_step,nDTOut)==0
      save(['./FloesCom2/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','eularian_data','SigXXa','SigXYa', 'SigYYa','DivSigXa','DivSigYa','DivSig1a','DivSig2a','U','dU','Fx','mass','c2_boundary_poly');
        SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
        SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
        DivSigX = zeros(Ny, Nx); DivSig1 = zeros(Ny, Nx);
        DivSigY = zeros(Ny, Nx); DivSig2 = zeros(Ny, Nx);
        Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
        Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
        U = zeros(Ny, Nx); V = zeros(Ny, Nx);
        dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
        Fx = zeros(Ny, Nx); Fy = zeros(Ny, Nx);
        Sig = zeros(Ny, Nx); mass = zeros(Ny, Nx);
        
        M = cat(1,Floe.mass);
        Mtot(im_num) = sum(M)+sum(Vdnew(:));
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    

    FloeOld2 = FloeOld; FloeOld = Floe;
    %Calculate forces and torques and intergrate forward
%     tic

    
    [Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, uright, 0, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nbound, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING);
    
    if AVERAGE
        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nbound,c2_boundary,dt,PERIODIC);
        SigXX = SigXX+squeeze(eularian_data.stressxx); SigYX = SigYX+squeeze(eularian_data.stressyx);
        SigXY = SigXY+squeeze(eularian_data.stressxy); SigYY = SigYY+squeeze(eularian_data.stressyy);
        Eux = Eux+squeeze(eularian_data.strainux); Evx = Evx+squeeze(eularian_data.strainvx);
        Euy = Euy+squeeze(eularian_data.strainuy); Evy = Evy+squeeze(eularian_data.strainvy);
        U = U+squeeze(eularian_data.u);V = V+squeeze(eularian_data.v);
        dU = dU+squeeze(eularian_data.du);dV = dV+squeeze(eularian_data.dv); 
        Fx = Fx+squeeze(eularian_data.force_x);Fy = Fy+squeeze(eularian_data.force_y);
        Sig = Sig+squeeze(eularian_data.stress);
        mass = mass+squeeze(eularian_data.Mtot);
        [DSig1, DSig2, DSigX, DSigY] = Calc_Stress(eularian_data,dt, c2_boundary);
        DivSigX = DivSigX+DSigX; DivSig1 = DivSig1+DSig1;
        DivSigY = DivSigY+DSigY; DivSig2 = DivSig2+DSig2;
    end
    
    FloeOld = Floe;
    if WELDING && mod(i_step,25)==0
        weldrate = 150;%Set rate at which floes will meld
        A=cat(1,Floe.area);
        if max(A) > Amax
           Amax = max(A);
        end
        FloeOld = Floe;
        if mod(i_step,5000)==0
            weldrate = 150;
            Floe = Weld_Floes_par2(Floe,Nbound,weldrate,c2_boundary,Amax/2,1,1);
        elseif mod(i_step,500)==0
            Floe = Weld_Floes_par2(Floe,Nbound,weldrate,c2_boundary,Amax/3,2,2);
        else
            Floe = Weld_Floes_par2(Floe,Nbound,weldrate,c2_boundary,Amax/3,3,3);
        end
    end
    live = cat(1,Floe.alive);
    Floe(live==0)=[]; 

    FloeOld = Floe;
    if FRACTURES && mod(i_step,50)==0 %&& im_num>150 
        compactness = sum(cat(1,Floe.area))/area(c2_boundary_poly);
        [Floe,Princ] = FracMohr2(Floe,Nbound,min_floe_size,compactness);

        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        if KEEP_MIN
            keep(cat(1,Floe.area)<min_floe_size)=1;
        end
    end
    live = cat(1,Floe.alive);
    Floe(live==0)=[]; 
    
    if CORNERS && mod(i_step,10)==0

        keep = rand(length(Floe),1)>0.7;
        keep(1:Nbound) = ones(Nbound,1);
        if KEEP_MIN
            keep(cat(1,Floe.area)<min_floe_size)=1;
        end
         overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
         keep(overlapArea>0.15) = 0;
        fracturedFloes=corners(Floe(~keep),Nbound);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
    end   
    
    
    %Advect the dissolved mass
    Area=cat(1,Floe.area);
    if ~KEEP_MIN
        dissolvedNEW = calc_dissolved_mass(Floe(Area<min_floe_size),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
    end
    Vdnew = Vd(:,:,1)+dissolvedNEW;
    dissolvedNEW=zeros(Ny,Nx);
    Vd(:,:,2) = Vd(:,:,1);
    Vd(:,:,1) = Vdnew;
    
    Area=cat(1,Floe.area);
    if ~KEEP_MIN
        if sum(Area<min_floe_size)>0, display(['num of small floes killed:' num2str(sum(Area<min_floe_size))]); end
        Floe=Floe(Area> min_floe_size);
    end
    live = cat(1,Floe.alive);
    Floe(live == 0) = [];
    
    if sum(cat(1,Floe.area))/area(c2_boundary_poly)>1.05
       % xx = 1; xx(1) = [1 2];
    end
    
    Time=Time+dt; i_step=i_step+1; %update time index
    if i_step < 750 %&& mod(i_step,10)==0 
        xb = c2_boundary(1,:);
        yb = c2_boundary(2,:);
        yb = yb - 4*[-1 1 1 -1 -1];
        c2_boundary = [xb; yb];
        Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
        c2_boundary_poly = polyshape(c2_boundary');
        c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
        floebound = initialize_floe_values(c2_border, height);
    end
    
end

tEnd = toc(tStart)
%%
