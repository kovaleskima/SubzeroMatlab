function [fig] =plot_stress(fig, Time,Floe,ocean,c2_boundary_poly,Nb,Nbond)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 
px = [-Lx -Lx Lx Lx];py = [Ly Ly+10000 Ly+10000 Ly];
p1 = polyshape([px;py]'); p2 = polyshape([px;-py]');
live = cat(1,Floe.alive);
Floe(live == 0) = [];

%% Set Up The Plots

ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig)
clf(fig);

dn=1; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

% title(['Time = ' num2str(Time/3600) ' hours'],'fontsize',24);
Floe_bonds = Floe(1:Nbond);
Floes = Floe(Nbond+1:end);
for ii =1:length(Floe_bonds)
    Floe_bonds(ii).poly = polyshape(Floe_bonds(ii).c_alpha'+[Floe_bonds(ii).Xi Floe_bonds(ii).Yi]);
    Stress_bnds = eig(Floe_bonds(ii).Stress);
    Princ1_bnds(ii) = max(abs(Stress_bnds));
end
for ii =1:length(Floes)
    Floes(ii).poly = polyshape(Floes(ii).c_alpha'+[Floes(ii).Xi Floes(ii).Yi]);
    Tmax(ii) = sqrt((Floes(ii).Stress(1,1)-Floes(ii).Stress(2,2)).^2/4+Floes(ii).Stress(1,2).^2);
%     Stress = eig(Floes(ii).Stress);
%     Princ1(ii) = max(abs(Stress));
end
% Stress = cat(1,Floe.MaxShear);

%% Plot the Floes and Ghost Floes
Shear_bnds = abs(Princ1_bnds);
S1_bnds = mean(Shear_bnds);
C2_bnds = Shear_bnds/S1_bnds;
C2_bnds(C2_bnds<2)=0;
C2_bnds(C2_bnds>1)=1;

% Shear = abs(Princ1);
S1 = 20000;%mean(Shear);
C2 = Tmax/S1;
% C2(C2<2)=0;
C2(C2>1)=1;
%else
%    C2 = (Shear-min(Shear))/(max(Shear)-min(Shear));
%end
% if max(Shear)>0

plot([Floe_bonds.poly],'FaceColor',[0 0 0],'EdgeColor',[1 1 1]*0.2,'linestyle','none');
for i = 1:length(Floes)
    plot(Floes(i).poly,'FaceColor',[1 0 0]*C2(i),'EdgeColor',[1 1 1]*0.2);
end

% plot([Floes(C2==0).poly],'FaceColor',[0 0 0],'EdgeColor',[1 1 1]*0.2,'linestyle','none');
% plot([Floes(C2==1).poly],'FaceColor',[0 1 1]*0.7,'EdgeColor',[1 1 1]*0.2,'linestyle','none');
% plot([Floe_bonds(C2_bnds==0).poly],'FaceColor',[0 0 0],'EdgeColor',[1 1 1]*0.2,'linestyle','none');
%plot([Floe_bonds(C2_bnds==1).poly],'FaceColor',[1 0 0],'EdgeColor',[1 1 1]*0.2,'linestyle','none','FaceAlpha',0.6);%     for i = 1:length(Floe)

%         plot(Floe(i).poly,'FaceColor',[1 0 0]*C2(i),'EdgeColor',[1 1 1]*0.2);
% %         plot(Floe(i).poly,'FaceColor',[1 0 0]*Stress(i)/max(Stress),'EdgeColor',[1 1 1]*0.2);
%     end
% else
%     plot([Floe.poly],'FaceColor',[1 0 0]*0,'EdgeColor',[1 1 1]*0.2);
% %     plot([Floe.poly],'FaceColor','r','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
% end
if Nb > 0
    plot([Floe(1:Nb).poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
end
plot(p1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]*0.2); plot(p1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]*0.2);
plot(p2,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]*0.2); plot(p2,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]*0.2);

X = cat(1,Floes.Xi); Y = cat(1,Floes.Yi);
U = cat(1,Floes.Ui); V = cat(1,Floes.Vi);
quiver(X,Y,U,V,'k','linewidth',2)


set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
% plot(xb,yb, 'k-','linewidth',2);


colormap('gray'); caxis([0 1]);
axis([-Lx+Lx/10 Lx-Lx/10 -5e4 5e4])
%axis([-2e5 2e5 -Ly-Ly/10 Ly+Ly/10])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
set(gca,'xtick',[])
set(gca,'ytick',[])

drawnow

warning('on',id)
warning('on',id3)
end