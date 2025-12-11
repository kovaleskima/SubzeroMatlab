function [fig] =plot_basic_bonds(fig, Floe,ocean,c2_boundary_poly,Nbound,Nbond,PERIODIC)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 
live = cat(1,Floe.alive);
Floe(live == 0) = [];

N0=length(Floe);
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    translation = [];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            parent=[parent  i];
            translation = [translation; -2*Lx*sign(x(i)) 0];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            parent=[parent  i];
            translation = [translation; 0 -2*Ly*sign(y(i))];

        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%Find length of new Floe variable including the ghost floes
N=length(Floe);

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

grid
ax = gca;
ax.GridAlpha = 0.5;

%% Plot the Floes and Ghost Floes
FloeNumbers = 1:length(Floe);
Floes = Floe(Nbound+1:end);
Nums = cat(1,Floe.num);
clear f1; clear f2
Xc = []; Xb = [];
Yc = []; Yb = [];
Xi = cat(1,Floe.Xi); Yi = cat(1,Floe.Yi);
count = 1;
for ii = 1:length(Floes)
    BondNum = unique(cat(1,Floes(ii).bonds.Num));
    Xb = [Xb; cat(1,Floe(ii).bonds.Xb)+Xi(ii)]; 
    Xc = [Xc; cat(1,Floe(ii).bonds.Xc)+Xi(ii)]; 
    Yb = [Yb; cat(1,Floe(ii).bonds.Yb)+Yi(ii)]; 
    Yc = [Yc; cat(1,Floe(ii).bonds.Yc)+Yi(ii)]; 
    for jj = 1:length(BondNum)
        f1(count) = Floes(ii).num;
        f2(count) = BondNum(jj);
        count = count+1;
    end
end
f2(f1 ==0) = [];f1(f1 ==0) = [];
G = digraph(f1,f2);
[bins,binsizes] = conncomp(G,'Type','weak');
Lia = ismember(FloeNumbers,Nums);Lia(1:Nbond) = 0;
clear poly
for iii = 1:length(Floe)
    poly(iii) = polyshape(Floe(iii).c_alpha'+[Floe(iii).Xi Floe(iii).Yi]);
end
count = 1;
clear p
for ii = 1+Nbond:length(binsizes)
    if Lia(ii)
        Lia1 = ismember(bins,ii);
        Lia2 = ismember(Nums,FloeNumbers(Lia1));
        ptmp = union([poly(Lia2)]);
        p(count) = ptmp;%rmholes(ptmp);
        count = count+1;  
    end
end
plot(poly(1+Nbound:length(Floe)),'FaceColor','none','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
plot(Xc,Yc,'r.','linewidth',3)

hold on
colors = distinguishable_colors(length(p));
for ii = 1:length(p)
    plot(p(ii),'FaceColor',colors(ii,:) )
end


set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
plot(xb,yb, 'k-','linewidth',2);

colormap('gray'); caxis([0 1]);
axis([-65000 65000 -65000 65000])
% axis([000 20000 00 20000])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
% set(gca,'xtick',[])
% set(gca,'ytick',[])

% drawnow
end