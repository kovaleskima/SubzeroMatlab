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

KEEP_MIN = true;

SIMPLIFY = false;

ifPlot = false; %Plot floe figures or not?

ifPlotStress = false;

ifPlotStressStrain = false;

justfrac = false;

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
floebound = initialize_floe_values(c2_border, height,1);
uright = 0; uleft = 0; %Define speeds that boundaries might be moving with
min_floe_size = 2*Lx*Ly/10000;% Define the minimum floe size you want in initial configuration

%Initialize Floe state
target_concentration = 1;
%[Floe,bonds, Nb,Nbond] = initial_concentration(c2_boundary,target_concentration,height,1500,1,min_floe_size);
%save('FloeInit.mat','Floe','bonds','Nbond','Nb');
load FloeInit2
%load('./Floes_bnds2/Floe0000001.mat','Floe','Nbond','Nb');
Floe0 = Floe;
Nums = cat(1,Floe.num);
for ii = 1:length(Floe)
    floe = Floe(ii);
    bnds = unique(cat(1,Floe(ii).bonds.Num));
    bonds1 = cat(1,Floe(ii).bonds.Num);
    for jj = 1:length(bnds)
        Lia = ismember(Nums,bnds(jj)); floe2 = Floe(Lia); floeNum = floe2.num;
        Num1 = sum(ismember(bonds1,floeNum)); BondNum2 = cat(1,floe2.bonds.Num); Num2 = sum(ismember(BondNum2,Floe(ii).num));
        if ~Num1
            xx = 1; xx(1) =[1 2];
        elseif ~Num2
            xx = 1; xx(1) =[1 2];
        elseif sum(Num1)~=sum(Num2)
            xx = 1; xx(1) =[1 2];
        end
    end
end
%Floe(1).Ui = 0.5; %Floe(2).Ui = -0.5;
%Floe(1).Vi = 0.1; Floe(2).Vi = -0.5;
% Floe(1).ksi_ice = 0.01; Floe(2).ksi_ice = 0;
% Floe(2) = [];
% save('Floe0.mat','Floe','Nb','Nbond');
% xx = 1; xx(1) =[1 2];
% Floe0 = Floe;
% save('Floe0.mat','Floe0');
%load('./Floes_bnds/Floe0000001.mat','Floe','Nbond','Nb');
% for jj = 1+Nbond:length(Floe)
%     BondNum = Floe(jj).bonds.BondNum;
%     if ~isempty(Floe(jj).bonds.poly)
%         R = regions(polyshape(Floe(jj).bonds.poly));
%         if ~(length(R)==length(BondNum))
%             xx = 1; xx(1) =[1 2];
%         end
%     end
% end

if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
% Floe = Floe(Nbond+1:end);
min_floe_size = (4*Lx*Ly-sum(cat(1,Floe(1:Nb).area)))/20000; %define minimum floe size that can exist during model run

% load Floe0; Nbond = 1;
% Floe(2).Ui = 0.15;
% Floe(3).Ui = -0.15;

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

%%

dhdt = 1; %Set to 1 for ice to grow in thickness over time

nDTOut=150; %Output frequency (in number of time steps)

nSnapshots=800; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nSimp = 20;

nPar = 6; %Number of workers for parfor
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(nPar);
else
    delete(poolobj);
    parpool(nPar);
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

%Initiailize Eulearian Data
[eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,Nbond,c2_boundary,dt,PERIODIC);
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
[Floe,dissolvedNEW] = floe_interactions_all(Floe,floebound, uright, 0, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);


%% Initialize time and other stuff to zero
if isempty(dir('Floes_bnds2')); disp('Creating folder: Floes_bnds2'); mkdir('Floes_bnds2'); end
if isempty(dir('./Floes_bnds2/figs')); disp('Creating folder: figs'); mkdir('./Floes_bnds2/figs'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=figure('Position',[100 100 1000 500],'visible','on');
    fig3 = figure;
end
% im_num=51;
% load(['./Floes_bnds2/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','c2_boundary');
% i_step = (im_num-1)*nDTOut;
% xb = c2_boundary(1,:);
% yb = c2_boundary(2,:);
% % yb = yb - 2.5*[-1 1 1 -1];
% % c2_boundary = [xb; yb];
% Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
% c2_boundary_poly = polyshape(c2_boundary');
% c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
% floebound = initialize_floe_values(c2_border, height, 1);
% Time = i_step*dt;
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
        

        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,Nbond,c2_boundary,dt,PERIODIC);
        if ifPlot
            [fig] =plot_basic_bonds(fig,Floe,ocean,c2_boundary_poly,Nb,Nbond,PERIODIC);
%            [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nb,PERIODIC);
            exportgraphics(fig,['./figs/' num2str(im_num,'%03.f') '.jpg']);
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
        
                
        if ifPlotStress
            [fig] =plot_basic_stress(fig, Time,Floe,ocean,c2_boundary_poly,Nb,bonds);
            saveas(fig,['./Floes_bnds2/figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
%             figure(2)
%             plot(Xc,SigXXa,'kx','linewidth',2); title(['Time = ' num2str(Time/3600) ' hours'],'fontsize',24);
            drawnow
        end
        
        if ifPlotStressStrain
            fig2 = figure(fig2);
            SigO = Sig;
            subplot(2,4,1); imagesc(Xc,Yc,SigXXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{xx}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,2); imagesc(Xc,Yc,SigYXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{yx}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,5); imagesc(Xc,Yc,SigXYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{xy}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,6); imagesc(Xc,Yc,SigYYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{yy}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,3); imagesc(Xc,Yc,Eux); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{11}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,4); imagesc(Xc,Yc,Evx); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{21}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,7); imagesc(Xc,Yc,Euy); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{12}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,8); imagesc(Xc,Yc,Evy); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{22}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            saveas(fig2,['./figs/' num2str(im_num,'Stress%03.f') '.jpg'],'jpg');
            fig3 = figure(fig3);
            subplot(2,2,1); imagesc(Xc,Yc,DivSigXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{x}~Homogenized$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            subplot(2,2,2); imagesc(Xc,Yc,DivSigYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{y}~Homogenized$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            subplot(2,2,3); imagesc(Xc,Yc,DivSig1a); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{x}~Momentum$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            subplot(2,2,4); imagesc(Xc,Yc,DivSig2a); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{y}~Momentum$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            saveas(fig3,['./figs/' num2str(im_num,'DivStress%03.f') '.jpg'],'jpg');
        end
        
    end
    
    if PACKING && h0 > 0
        if mod(i_step,nDTpack)==0
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice_new(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean, height, min_floe_size, PERIODIC,3,3,Nb);
        end
    end
    
    
    if mod(i_step,nDTOut)==0
        save(['./FloesOut/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','Nb','Nbond','eularian_data','SigXXa','SigXYa', 'SigYYa','U','dU','Fx','mass','c2_boundary');
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
    
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW] = floe_interactions_all(Floe,floebound, uright, 0, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING);

    
    if AVERAGE
        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,Nbond,c2_boundary,dt,PERIODIC);
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


    if FRACTURES && mod(i_step,50)==0 %&& im_num > 40
%         for ii = 1:length(bonds)
%             for jj = 1:length(bonds(ii).bond)
%                 bonds(ii).bond(jj).Stress = bonds(ii).bond(jj).Stress/250;
%             end
%         end
        compactness = sum(cat(1,Floe.area))/area(c2_boundary_poly);
%        yb = c2_boundary(2,:);
%        if max(yb) <= 49000
%            save('FloeFail.mat','Floe','Nb','min_floe_size','compactness');
%            xx = 1; xx(1) =[1 2];
%        end
        [Floe,Princ] = FracMohr(Floe,Nb,min_floe_size,compactness);
%         Areas = cat(1,Floe.area);
%         Nbond = length(Areas(Areas<bond_area+1));
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
%     Vdnew = Advect_Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt);
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
    
    if max(c2_boundary(2,:)) > 48000 && mod(i_step,20)==0 && side == 1
        xb = c2_boundary(1,:);
        yb = c2_boundary(2,:);
        xb = xb + 2.5*[-1 -1 1 1];
        yb = yb - 2.5*[-1 1 1 -1];%     y=[-1 1 1 (-1-0.1*Time*f/pi*sin(f*Time)) -1]*Ly;
        %     xb = xb + [1 1 -1 -1 1];
        %    Ly = max(y); %Lx = area(c2_boundary_poly)/(4*Ly);
        %   x=[-1 -1 1 1 -1]*Lx;
        c2_boundary = [xb; yb];
        Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
        c2_boundary_poly = polyshape(c2_boundary');
        c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
        floebound = initialize_floe_values(c2_border, height, 1);
        if max(yb) <= 48000
            Adomain = area(c2_boundary_poly);
            for jj = 1:length(Floe)
                xmax(jj) = max(abs(Floe(jj).Xi+Floe(jj).c_alpha(1,:)));
            end
            Lx = max(xmax); Ly = Adomain/(Lx*4);
            xb = Lx*[-1 -1 1 1]; 
            yb = Ly*[-1 1 1 -1];
            c2_boundary_poly = polyshape(c2_boundary');
            c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
            floebound = initialize_floe_values(c2_border, height, 1);
            side = 2;
        end
    elseif max(c2_boundary(1,:)) > 48000 && mod(i_step,20)==0 && side == 2
        xb = c2_boundary(1,:);
        yb = c2_boundary(2,:);
        xb = xb - 2.5*[-1 -1 1 1];
        yb = yb + 2.5*[-1 1 1 -1];%     y=[-1 1 1 (-1-0.1*Time*f/pi*sin(f*Time)) -1]*Ly;
        %     xb = xb + [1 1 -1 -1 1];
        %    Ly = max(y); %Lx = area(c2_boundary_poly)/(4*Ly);
        %   x=[-1 -1 1 1 -1]*Lx;
        c2_boundary = [xb; yb];
        Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
        c2_boundary_poly = polyshape(c2_boundary');
        c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
        floebound = initialize_floe_values(c2_border, height, 1);
        if max(xb) <= 48000
            side = 3;
        end
    end
    
    Time=Time+dt; i_step=i_step+1; %update time index
    
end

tEnd = toc(tStart)