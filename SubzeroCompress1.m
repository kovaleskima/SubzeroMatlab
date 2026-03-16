close all; clear all;
paths
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
nDTpack = 4500; Nfrac = 10;
rho_ice=920;
% [ocean, HFo, h0]=initialize_ocean(dt,nDTpack);
[ocean, HFo, h0]=initialize_ocean(dt,nDTpack);
% [ocean, HFo, h0]=initialize_ocean_Nares(dt,nDTpack);

Lxo = max(ocean.Xo);

%Define 10m winds
clear winds
clear Time
[nyo,nxo] = size(ocean.Xocn);
U0 = 0; V0 = 0;
winds.x = ocean.Xocn; winds.y = ocean.Yocn;
winds.u=U0*ones(nyo,nxo);%.*(exp(-ocean.Yo.^2/(300*max(ocean.Yo)))'*exp(-ocean.Xo.^2/(300*max(ocean.Xo))));
winds.v=V0*ones(nyo,nxo);%.*(exp(-ocean.Yo.^2/(300*max(ocean.Yo)))'*exp(-ocean.Xo.^2/(300*max(ocean.Xo))));

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
min_floe_size = 4*Lx*Ly/10000;%4*Lx*Ly/250000;

%Initialize Floe state
%Floe=initialize_Floe('FloeShapes.mat');
% target_concentration = ones(1,1);
target_concentration = 1;%[1;1;0];
%[Floe, Nb] = initial_concentration_Nares(c2_boundary,target_concentration,height,50,min_floe_size);
[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,100,50,min_floe_size);
xx = 1; xx(1) =[1 2];
%[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,100,min_floe_size);
%load(['./FloesS3/Floe0000001.mat']); Nb = 0;
x=[-1 -1 1 1 -1]*Lx; 
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
%Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
c2_boundary_poly = polyshape(c2_boundary');
c2_border = polyshape(2*[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly]'); c2_border = subtract(c2_border, c2_boundary_poly);
floebound = initialize_floe_values(c2_border, height,1);
floebound = initialize_floe_values(c2_border, height,1);
% load Floe1; Nb = 0;
% Floe = Floe(1:200);
%im_num=94;
%load(['./Floes/Floe' num2str(im_num,'%07.f') '.mat']);
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
% load floetest
min_floe_size = (4*Lx*Ly-sum(cat(1,Floe(1:Nb).area)))/20000;
global Modulus
%load Modulus
collide = 0;
Modulus = 1.5e3*(mean(sqrt(cat(1,Floe.area)))+min(sqrt(cat(1,Floe.area))));
save('Modulus.mat','Modulus');
%load(['./FloesS3/Floe0000404.mat']);
%load('PackedFloesFullDomain.mat','Floe');
%Floe= create_packed_domain();
% load F_dir
% clear fdir
fdir.dir =[0;0];fdir.loc =[0 0];
% f_hist_t =[0;0];
save('F_dir.mat','fdir')
frac = 0;

% load Fcompress1.mat
% x = c2_boundary_poly.Vertices(:,1)';
% y = c2_boundary_poly.Vertices(:,2)';
% c2_boundary = [x; y];
% Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
% c2_border = polyshape(2*[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly]'); c2_border = subtract(c2_border, c2_boundary_poly);
% floebound = initialize_floe_values(c2_border, height);

%%

dhdt = 1;

nDTOut=250; %Output frequency (in number of time steps)

nSnapshots=20; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nSimp = 20;

%nPar = 36; %Number of workers for parfor
%poolobj = gcp('nocreate'); % If no pool, do not create new one.
%if isempty(poolobj)
%    parpool(nPar);
%else
%    delete(poolobj);
%    parpool(nPar);
%end

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
[eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
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
[Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, uright, 0, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);

%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end
if isempty(dir('FloesCom1')); disp('Creating folder: FloesCom1'); mkdir('FloesCom1'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=figure('Position',[100 100 1000 500],'visible','on');
    fig3 = figure;
end
% im_num = 116;
% i_step = (im_num-1)*nDTOut;

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
%         psi_ocean=transport/1*(sin(2*ocean.kx*ocean.Xocn).*cos(2*ocean.ky*ocean.Yocn)+cos(ocean.fCoriolis*Time)*sin(ocean.kx*ocean.Xocn).*cos(ocean.ky*ocean.Yocn)); 
% 
%         Uocn=zeros(size(ocean.Xocn)); Vocn=zeros(size(ocean.Xocn));
%         ocean.Uocn(2:end,:)=-(psi_ocean(2:end,:)-psi_ocean(1:end-1,:))/dXo;
%         ocean.Vocn(:,2:end)=(psi_ocean(:,2:end)-psi_ocean(:,1:end-1))/dXo;
        
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


%     if mod(i_step,nSimp)==0
%     FloeOld = Floe;
%     parfor j=1:length(FloeOld)
%         FloeOld(j).poly = polyshape(FloeOld(j).c_alpha'+[FloeOld(j).Xi FloeOld(j).Yi]);
%     end
%         floenew = [];
%         for ii = 1:length(Floe)
%             ParFloes(ii).floenew = [];
%             ParFloes(ii).kill = [];
%             ParFloes(ii).verts = 0;
%         end
%         if isfield(Floe,'poly')
%             Floe=rmfield(Floe,{'poly'});
%         end
%         parfor ii = 1:length(Floe)
%             floe = Floe(ii);
%             if length(Floe(ii).c0) > 30
%                 [floe2,kill] = FloeSimplify(Floe(ii),0,FloeOld);
%                 if isfield(floe2,'poly')
%                     floe2=rmfield(floe2,{'poly'});
%                 end
%                 if isempty(kill)
%                     ParFloes(ii).kill = [ParFloes(ii).kill kill];
%                 end
%                 for jj = 1:length(floe2)
%                     if jj == 1
%                         Floe(ii) = floe2(jj);
%                         ParFloes(ii).verts=length(floe2(jj).c0(1,:));
%                     else
%                         ParFloes(ii).floenew = [ParFloes(ii).floenew floe2(jj)];
%                     end
%                 end
%             end
%         end
%         if isfield(Floe,'poly')
%             Floe=rmfield(Floe,{'poly'});
%         end
%         floenew =[]; kill = [];
%         for ii = 1:length(ParFloes)
%             floenew = [floenew ParFloes(ii).floenew];
%             kill = [kill ParFloes(ii).kill];
%         end
%         kill = unique(kill(kill>0));
%         if ~isempty(kill)
%             save('FloeF.mat','Floe','kill','FloeOld')
%             xx = 1; xx(1) =[1 2];
%         end
%         live = cat(1,Floe.alive);
%         Floe(live==0)=[]; 
%         Floe =[Floe floenew];
%         h = cat(1,Floe.h); h2 = cat(1,FloeOld.h);
%         uice = cat(1,Floe.Ui); vice = cat(1,Floe.Vi);
%         if max(abs(uice)) >100 || max(abs(vice)) >100
% %        if sum(cat(1,Floe(h>0.1).area))/sum(cat(1,FloeOld(h2>0.1).area)) < 0.98
%             save('FloeF.mat','Floe','FloeOld')
%             xx = 1; xx(1) =[1 2];
%         end
%     end

    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        

        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
%         [eularian_data] = calc_eulerian_data(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
        if ifPlot
            [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nb,PERIODIC);
%            img =  getframe(gcf);
%            savepng(img.cdata,['./figs/' num2str(im_num,'%03.f') '.jpg']);
%             [fig] =plot_Nares(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
%              saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
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
            [fig] =plot_basic_stress(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
            figure(2)
            plot(Xc,SigXXa,'kx','linewidth',2); title(['Time = ' num2str(Time/3600) ' hours'],'fontsize',24);
            drawnow
        end
        
        if ifPlotStressStrain
            fig2 = figure(fig2);
%             if im_num==3
%                 xx = 1; xx(1) =[1 2];
%             end
            SigO = Sig;
%             if Time > 3900
%                 xx = 1; xx(1) =[1 2];
%             end
%             imagesc(Xc,Yc,abs(Sig)); colorbar;%title('$\sigma_{xx}$','interpreter','latex','fontsize',16); colorbar; caxis([0 1e3])
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
        
%         if im_num > 60
%             xx = 1;
%             xx(1) =[1 2];
%         end
    end
    
    FloeOld2 = FloeOld; FloeOld = Floe;
    if PACKING && h0 > 0
        if mod(i_step,nDTpack)==0
            height.mean = 0.2;%h0;
            height.delta = 0;
%             xx = 1; xx(1) =[1 2];
            [Floe,Vd] = pack_ice_new(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean, height, min_floe_size, PERIODIC,3,3);
        end
%         for ii = 1+Nb:length(Floe)
%             h1 = Floe(ii).mass/(Floe(ii).area*rho_ice);
%             if abs(Floe(ii).h -h1) > 0.05
%                 xx = 1;
%                 xx(1) = [1 2];
%             end
%         end
    end
    
    live = cat(1,Floe.alive);
    Floe(live==0)=[]; 
    h = cat(1,Floe.h); h2 = cat(1,FloeOld.h);
    uice = cat(1,Floe.Ui); vice = cat(1,Floe.Vi);
    if max(abs(uice)) >100 || max(abs(vice)) >100
%    if sum(cat(1,Floe(h>0.1).area))/sum(cat(1,FloeOld(h2>0.1).area)) < 0.98
        save('FloeF.mat','Floe','FloeOld')
        %xx = 1; xx(1) =[1 2];
    end
    
    if mod(i_step,nDTOut)==0
      save(['./FloesCom1/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','eularian_data','SigXXa','SigXYa', 'SigYYa','DivSigXa','DivSigYa','DivSig1a','DivSig2a','U','dU','Fx','mass','c2_boundary_poly');
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

    if sum(cat(1,Floe.area))/area(c2_boundary_poly)>1.05
        %xx = 1; xx(1) = [1 2];
    end
    [Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, uright, 0, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING);
    h = cat(1,Floe.h); h2 = cat(1,FloeOld.h);
    uice = cat(1,Floe.Ui); vice = cat(1,Floe.Vi);
    if max(abs(uice)) >100 || max(abs(vice)) >100
%    if sum(cat(1,Floe(h>0.1).area))/sum(cat(1,FloeOld(h2>0.1).area)) < 0.93
        save('FloeF.mat','Floe','FloeOld')
        xx = 1; xx(1) =[1 2];
    end

%     toc
%     for ii =1:length(Floe)
%         poly(ii) = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
%     end
%     if abs(Floe(1).Fx)>0
%         collide = 1;
%     end


%     if Time == 24620
%         xx = 1; xx(1) =[1 2];
%     end
    
    if AVERAGE
        %         [eularian_data] = calc_eulerian_stress(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
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
    
    if sum(cat(1,Floe.area))/area(c2_boundary_poly)>1.05
        %xx = 1; xx(1) = [1 2];
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
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/2,1,1);
        elseif mod(i_step,500)==0
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/3,2,2);
        else
            Floe = Weld_Floes_par2(Floe,Nb,weldrate,c2_boundary,Amax/3,3,3);
        end
    end
%    if WELDING && mod(i_step,7)==0
%        weldrate = 75;%Set rate at which floes will meld
%        A=cat(1,Floe.area);
%        if max(A) > Amax
%           Amax = max(A);
%        end
%        FloeOld = Floe;
%        keep = zeros(1,length(Floe));
%        keep(A>mean(A)/5)=1;
%        keep(1:Nb) = ones(Nb,1);
%        keep = logical(keep);
%        Floe = Weld_Floes(Floe,Nb,weldrate,Amax,mean(A)/5);
%        if ~isempty(FloesW)
%            Floe=[Floe(keep) FloesW];
%        end
%    end
    live = cat(1,Floe.alive);
    Floe(live==0)=[]; 
    h = cat(1,Floe.h); h2 = cat(1,FloeOld.h);
    uice = cat(1,Floe.Ui); vice = cat(1,Floe.Vi);
    if max(abs(uice)) >100 || max(abs(vice)) >100
%    if sum(cat(1,Floe(h>0.1).area))/sum(cat(1,FloeOld(h2>0.1).area)) < 0.98
        save('FloeF.mat','Floe','FloeOld')
        xx = 1; xx(1) =[1 2];
    end

    if sum(cat(1,Floe.area))/area(c2_boundary_poly)>1.05
        %xx = 1; xx(1) = [1 2];
    end

    FloeOld = Floe;
    if FRACTURES && mod(i_step,Nfrac)==0 && i_step>50 
        for ii = 1:length(Floe)
            [fracFloe] = calc_stress_bonds(Floe,Nfrac,dt);
            if length(fracFloe)>1
                Floe(ii) = [];
                Floe = [Floe fracFloe];
            else
                Floe(ii) = fracFloe;
            end
        end
        compactness = sum(cat(1,Floe.area))/area(c2_boundary_poly);
%         [Floe] = FracMohrSubfloes(Floe,Nb,min_floe_size);
        if sum(cat(1,Floe.area))/(4*Lx^2)-1 > 0.005
            xx = 1; xx(1) =[1 2];
        end
%        save(['./FloesCom/PrincO' num2str(im_num,'%07.f') '.mat'],'Princ');
%         clear Stress, clear A;
%         for ii = 1:length(Floe)
%             Stress(ii) = max(abs(eig(Floe(ii).Stress)));
%             A(ii) = Floe(ii).area;
%         end
%         n = 6:0.5:ceil(max(log10(A)));
%         BinEdges = 10.^n;
%         for ii = 1:length(n)-1
%             Stresses = Stress;
%             Stresses(A<BinEdges(ii)) = nan; Stresses(A>BinEdges(ii+1)) = nan;
%             Mstress = max(Stress(~isnan(Stresses)));
%             Stress(~isnan(Stresses)) = Stress(~isnan(Stresses))/Mstress;
%         end
%         [B,TF] = rmoutliers(Stress);
%         keep=rand(length(Floe),1)>Stress'/max(Stress(TF));
%         keep=rand(length(Floe),1)<1.05-Stress'/max(Stress);
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
%        keep=rand(length(Floe),1)>0.001;%2*overlapArea;
%        keep=rand(length(Floe),1)>overlapArea;
        if KEEP_MIN
            keep(cat(1,Floe.area)<min_floe_size)=1;
        end
%        keep(1:Nb) = ones(Nb,1);
%        fracturedFloes=fracture_floe(Floe(~keep),3);
% %         fracturedFloes=fracture_leads(Floe(~keep),Nx,Ny,c2_boundary,eularian_data);
%        if ~isempty(fracturedFloes)
%            Floe=[Floe(keep) fracturedFloes];
%        end
% %         for ii = 1+Nb:length(Floe)
% %             h1 = Floe(ii).mass/(Floe(ii).area*rho_ice);
% %             if abs(Floe(ii).h -h1) > 0.05
% %                 xx = 1;
% %                 xx(1) = [1 2];
% %             end
% %         end
%         %         [Floe] = FracLeads(Floe,Ny,Nx,Nb,c2_boundary,eularian_data);
%         %         [Floe] = FracIso(Floe,Ny,Nx,Nb,c2_boundary,SigO);
%        Princ = zeros(length(Floe),2);
%        Pstar = 7.25e5; C = 20;
%        Amohr=cat(1,Floe.area);
%        concentration = sum(cat(1,Floe.area))/area(c2_boundary_poly);
%        h = mean(cat(1,Floe.h));
%        P = Pstar*h*exp(-C*(1-concentration));
%        t = linspace(0,2*pi) ;
%        a = P*sqrt(2)/2 ; b = a/2 ;
%        x = a*cos(t) ;
%        y = b*sin(t) ;
%        Mohr = polyshape(x,y);
%        Mohr = rotate(Mohr,45);
%        Mohr = translate(Mohr,[-P/2, -P/2]);
%        for ii = 1:length(Floe)
%            Stress = eig(Floe(ii).Stress);
%            Princ(ii,1) = max(Stress);
%            Princ(ii,2) = min(Stress);
%            Princ(ii,3) = Floe(ii).area;
%        end
%        Princ1 = Princ(:,1); Princ2 = Princ(:,2); 
%        [in,out] = inpolygon(Princ1,Princ2,Mohr.Vertices(:,1), Mohr.Vertices(:,2));
%        keep = zeros(length(Floe),1);
%        keep(in) = 1;
%        keep(Amohr<min_floe_size)=1;
%        if length(keep)-sum(keep)>0
%            save('FloeFailMohr.mat','FloeOld','Nb','min_floe_size','c2_boundary_poly')
%        end
%        save(['./FloesCom/PrincN' num2str(im_num,'%07.f') '.mat'],'Princ');
    end
    live = cat(1,Floe.alive);
    Floe(live==0)=[]; 
    h = cat(1,Floe.h); h2 = cat(1,FloeOld.h);
    uice = cat(1,Floe.Ui); vice = cat(1,Floe.Vi);
    if max(abs(uice)) >100 || max(abs(vice)) >100
    %if sum(cat(1,Floe(h>0.1).area))/sum(cat(1,FloeOld(h2>0.1).area)) < 0.98
        save('FloeF.mat','Floe','FloeOld')
        xx = 1; xx(1) =[1 2];
    end
    
    if CORNERS && mod(i_step,10)==0
%         clear Stress, clear A;
%         for ii = 1:length(Floe)
%             Stress(ii) = max(abs(eig(Floe(ii).Stress)));
% %             A(ii) = Floe(ii).area;
%         end
%         [B,TF] = rmoutliers(Stress);
%         if max(Stress) == 0
%             keep = ones(length(Floe),1);
%         elseif max(Stress(~TF))==0
%             grind = rand(length(Floe),1)<2*Stress'/max(Stress);
%             keep = ~grind;
%         elseif max(Stress(~TF))<max(Stress(TF))
%             keep=rand(length(Floe),1)<2*Stress'/max(Stress(~TF));
%         else
%             grind = rand(length(Floe),1)<2*Stress'/max(Stress);
%             keep = ~grind;
%         end
%         n = 6:0.5:ceil(max(log10(A)));
%         BinEdges = 10.^n;
%         for ii = 1:length(n)-1
%             Stresses = Stress;
%             Stresses(A<BinEdges(ii)) = nan; Stresses(A>BinEdges(ii+1)) = nan;
%             Mstress = max(Stress(~isnan(Stresses)));
%             Stress(~isnan(Stresses)) = Stress(~isnan(Stresses))/Mstress;
%         end
%         for ii = 1+Nb:length(Floe)
%             h1 = Floe(ii).mass/(Floe(ii).area*rho_ice);
%             if abs(Floe(ii).h -h1) > 0.05
%                 xx = 1;
%                 xx(1) = [1 2];
%             end
%         end
%         Floe2 = Floe;
%         stress = zeros(length(Floe),1);
%         for ii = 1:length(Floe)
%             stress(ii) = trace(abs(Floe(ii).Stress));
%         end
%         if max(stress)>0
%             stresses=stress/max(stress);
%         else
%             stresses = stress;
%         end
%         keep=stresses<4*rand(length(Floe),1);
%         keep=rand(length(Floe),1)/2<1.05-Stress';
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
%         for ii = 1+Nb:length(Floe)
%             h1 = Floe(ii).mass/(Floe(ii).area*rho_ice);
%             if abs(Floe(ii).h -h1) > 0.05
%                 xx = 1;
%                 xx(1) = [1 2];
%             end
%         end
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
    
    if sum(cat(1,Floe.area))/area(c2_boundary_poly)>1.05
       % xx = 1; xx(1) = [1 2];
    end
    
    Time=Time+dt; i_step=i_step+1; %update time index
%     x=[-1 -1 1 (1-0.1*Time*f/pi*cos(f*Time)) -1]*Lx;
    if i_step < 100000 && mod(i_step,50)==0 
        xb = c2_boundary(1,:);
        yb = c2_boundary(2,:);
        yb = yb - 2.5*[-1 1 1 -1 -1];%     y=[-1 1 1 (-1-0.1*Time*f/pi*sin(f*Time)) -1]*Ly;
        %     xb = xb + [1 1 -1 -1 1];
        %    Ly = max(y); %Lx = area(c2_boundary_poly)/(4*Ly);
        %   x=[-1 -1 1 1 -1]*Lx;
        c2_boundary = [xb; yb];
        Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
        c2_boundary_poly = polyshape(c2_boundary');
        c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
        floebound = initialize_floe_values(c2_border, height, 1);
    end
    
end

% xx = 1; xx(1) =[1 2];

% clear Modulus
tEnd = toc(tStart)
%%
