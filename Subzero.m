close all; clearvars;

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

KEEP_MIN = true;

SIMPLIFY = false;

PLOT = true; %Plot floe figures or not?

%% Initialize model vars

%add paths
paths

dt=10; %Time step in sec
height.mean = 0.25;
height.delta = 0;

%Define ocean currents
nDTpack = 5500;
rho_ice=920;
[ocean, HFo, h0]=initialize_ocean(dt,nDTpack);

%Define 10m winds
winds.x = ocean.Xocn; winds.y = ocean.Yocn;
U0 = 0; V0 = 0;
winds.u=U0*ones(size(ocean.Xocn));
winds.v=V0*ones(size(ocean.Xocn));

%Define boundaries an d floes that act as boundaries
side = 1;
c2_boundary=initialize_boundaries();
Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
c2_boundary_poly = polyshape(c2_boundary');
c2_border = polyshape(2*[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly]'); c2_border = subtract(c2_border, c2_boundary_poly);

uright = -7; uleft = 0; %Define speeds that boundaries might be moving with
min_floe_size = 2*Lx*Ly/10000;% Define the minimum floe size you want in initial configuration

%Initialize Floe state
target_concentration = 1;
[Floe,bonds, Nb,Nbond] = initial_concentration(c2_boundary,target_concentration,height,100,1,min_floe_size);

Floe0 = Floe;

if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end

min_floe_size = (4*Lx*Ly-sum(cat(1,Floe(1:Nb).area)))/20000; %define minimum floe size that can exist during model run

Ly = max(c2_boundary(2,:));Lx = 1.25*max(c2_boundary(1,:));
c2_boundary =[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly];
c2_boundary_poly = polyshape(c2_boundary');
c2_border = polyshape(2*[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly]'); c2_border = subtract(c2_border, c2_boundary_poly);
floebound = initialize_floe_values(c2_border, height,1);

%Define Modulus for floe interactions
global Modulus r_mean L_mean
Modulus = 2.5e3*(mean(sqrt(cat(1,Floe.area)))+min(sqrt(cat(1,Floe.area))));
r_mean = mean(sqrt(cat(1,Floe.area)));
L = [];
for ii = 1:length(Floe)
    L_tmp=cat(1,Floe(ii).bonds.L);
    L = [L; L_tmp];
end
% bond_area = 0;%max(cat(1,Floe(1:Nbond).area));
L_mean = median(L);
save('Modulus.mat','Modulus','r_mean','L_mean');


dhdt = 1; %Set to 1 for ice to grow in thickness over time

nDTOut=10; %Output frequency (in number of time steps)

nSnapshots=5;  %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nSimp = 20;

% use pc=(...) to set local write directory for parpool on HPC cluster
pc = parcluster('Processes')
pc.JobStorageLocation = 'matlab_jobs';
nPar = 6; %Number of workers for parfor
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(pc, nPar);
else
    delete(poolobj);
    parpool(pc, nPar);
end

target_concentration=1;%Setting new target concentration to one if ice will freeze during winter for packing

tStart = tic; 

%set flags to be used for updating ice ocean interactions
doInt.flag = false; 
doInt.step = 10;

% specify coarse grid size
LxO= 2*max(ocean.Xo);LyO= 2*max(ocean.Yo);
Nx=10; Ny=10;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);

%Initialize Eulearian Data
[eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,Nbond,c2_boundary,dt,PERIODIC);
Vd = zeros(Ny,Nx,2);
Vdnew=zeros(Ny, Nx);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
[Floe,dissolvedNEW] = floe_interactions_all(Floe,floebound, uright, 0, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);


%% Initialize time and other stuff to zero
if isempty(dir('FloesOut')); disp('Creating folder: FloesOut'); mkdir('FloesOut'); end
if isempty(dir('./FloesOut/figs')); disp('Creating folder: figs'); mkdir('./FloesOut/figs'); end
if isempty(dir('./FloesOut/Floes')); disp('Creating folder: Floes'); mkdir('./FloesOut/Floes'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=figure('Position',[100 100 1000 500],'visible','on');
    fig3 = figure;
end

%% Solving for floe trajectories
tic; 
while side < 2.5
% while im_num<nSnapshots

    if mod(i_step,10)==0        
        disp(' ');
        toc
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        disp(['number of collisions: ' num2str(numCollisions)]);
        disp(' ');
        tic
        doInt.flag=false;
    else
        doInt.flag=false;
    end

    %Simplify the floe shapes if there are to many vertices
    if mod(i_step,nSimp)==0 && SIMPLIFY
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
        parfor ii = 1+Nb:length(Floe)
            floe = Floe(ii);
            if length(Floe(ii).c0) > 30
                [floe2,kill] = FloeSimplify(Floe(ii),0,FloeOld,polyAU);
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
        clear FloeOld
        for ii = 1:length(ParFloes)
            floenew = [floenew ParFloes(ii).floenew];
            kill = [kill ParFloes(ii).kill];
        end
        kill = unique(kill(kill>Nb));
        live = cat(1,Floe.alive);
        live(kill) = 0;
        Floe(live==0)=[];
        Floe =[Floe floenew];
    end

    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps

        % MOVE BOUNDARIES %
        xb = c2_boundary(1,:);
        yb = c2_boundary(2,:);
        yb = yb - 10*[-1 1 1 -1];
        c2_boundary = [xb; yb];
        Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
        c2_boundary_poly = polyshape(c2_boundary');
        new_boundary = subtract(c2_border, c2_boundary_poly);
        floebound = initialize_floe_values(new_boundary, height, 1);

        % if plotting eulerian data use this but I'm not for now
        %[eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,Nbond,c2_boundary,dt,PERIODIC);
        if PLOT
            fig = figure;
            axis equal
            axis manual
            [fig] =plot_basic_bonds(fig,Floe,ocean, Lx, Ly, c2_boundary_poly,Nb,Nbond,PERIODIC);
            exportgraphics(fig,['./FloesOut/figs/fig' num2str((i_step/10),'%03.f') '.jpg']);
            close(fig);
            save(['./FloesOut/Floes/Floe' num2str(i_step/10, '%03.f') '.mat'], 'Floe', 'bonds', 'Nbond', 'Nb');
        end
        
    end
    
    if PACKING && h0 > 0
        if mod(i_step,nDTpack)==0
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice_new(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean, height, min_floe_size, PERIODIC,3,3,Nb);
        end
    end
       
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW] = floe_interactions_all(Floe,floebound, uright, 0, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING);
    
    %Perform welding if selected to at designated rate
    if WELDING && mod(i_step,25)==0
        weldrate = 150;%Set rate at which floes will meld
        A=cat(1,Floe.area);
        if max(A) > Amax
           Amax = max(A);
        end
        %Do welding over different sized grid boxes
        if mod(i_step,5000)==0
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/2,1,1);
        elseif mod(i_step,500)==0
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/3,2,2);
        else
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/3,3,3);
        end
    end


    if FRACTURES && mod(i_step,10)==0 %&& im_num > 40
        compactness = sum(cat(1,Floe.area))/area(c2_boundary_poly);
        [Floe,Princ] = FracMohr(Floe,Nb,min_floe_size,compactness);
    end
 
    if CORNERS && mod(i_step,10)==0
        keep = rand(length(Floe),1)>0.7;
        keep(1:Nb) = ones(Nb,1);
        if KEEP_MIN
            keep(cat(1,Floe.area)<min_floe_size)=1;
        end
         overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
         keep(overlapArea>0.15) = 0;
        fracturedFloes=corners(Floe(~keep),Nb);
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
    
        
    Time=Time+dt; i_step=i_step+1; %update time index
    
end

tEnd = toc(tStart)
