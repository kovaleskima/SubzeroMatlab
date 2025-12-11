function [fig] =plot_Nares(fig, Time,Floe,ocean,c2_boundary_poly,Nb,eularian_data)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 

%% Set Up The Plots

% ratio=max(ocean.Yo)/max(ocean.Xo);
ratio = abs(-1.75e5/2)/(Lx); %Ly/Lx;%
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[10 10 200 200*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig)
clf(fig);

dn=1; % plot every dn'th velocity vector
% quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

% title(['Time = ' num2str(Time) ' s'],'fontsize',24);
for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
    Stress(ii) = max(abs(eig(Floe(ii).Stress)));
end

A = cat(1,Floe.area);

plot([Floe(1+Nb:end).poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
if Nb > 0
    plot([Floe(1:Nb).poly],'FaceColor',[0 0.2 0],'FaceAlpha',0.75,'EdgeColor',[1 1 1]*0.2);
end

[Ny, Nx]=size(eularian_data.u);%fix(Nx*LyO/LxO);
c2_boundary =c2_boundary_poly.Vertices';
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;
inP= zeros(Ny,Nx); 
%for jj = 1:Nb
[Xcg,Ycg] = meshgrid(Xc,Yc);
for jj = 1:Nb
    in = inpolygon(Xcg(:),Ycg(:),Floe(jj).poly.Vertices(:,1),Floe(jj).poly.Vertices(:,2));
    inP(:) = inP(:) + in;
end
quiver(Xcg(inP==0),Ycg(inP==0),eularian_data.u(inP==0),eularian_data.v(inP==0),'b','autoscalefactor',0.5)


set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
plot(xb,yb, 'k-','linewidth',2);
%title(['Time = ' num2str(Time) ' s']);
dim = [.71 .125 .175 .03];
Time = round(10*Time/(24*3600))/10;
%str = 'Time = ' num2str(Time) ' days';
annotation('textbox',dim,'String',['Time = ' num2str(Time) ' days'],'fontsize',16,'Color',[1, 1 ,1]);

colormap('gray'); caxis([0 1]);
axis([-Lx Lx -Ly Ly])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
% set(gca,'xtick',[])
% set(gca,'ytick',[])
ylim([-1e5 0.75e5])

drawnow
end