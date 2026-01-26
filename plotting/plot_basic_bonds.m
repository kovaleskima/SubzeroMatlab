function [fig] =plot_basic_bonds(fig, Floe, ocean, Lx, Ly, c2_boundary_poly,Nbound,Nbond,PERIODIC)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
Lxmax = Lx;
Lymax = Ly;
Lx = max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly = max(c2_boundary_poly.Vertices(:,2)); 
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
for iii = 1:Nbound
    poly_bound(iii) = polyshape(Floe(iii).c_alpha'+[Floe(iii).Xi Floe(iii).Yi]);
end
%% Set Up The Plots

ratio=1.616;%max(ocean.Yo)/max(ocean.Xo);
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
Floes = Floe(Nbound+1:end);
FloeNumbers = 1:length(Floes);
Nums = cat(1,Floes.num);
clear f1; clear f2
Xbond = []; Xbroken = [];
Ybond = []; Ybroken = [];
Xi = cat(1,Floes.Xi); Yi = cat(1,Floes.Yi);
count = 1;

for ii = 1:length(Floes)

    % Extract bonds
    bonds = Floes(ii).bonds;

    if isempty(bonds)
        continue
    end

    % differentiate intact / broken bonds
    isBroken = arrayfun(@(b) ~isempty(b.broken) && logical(b.broken), bonds);
    intactBnds  = bonds(~isBroken);
    brokenBnds  = bonds(isBroken);

    % rotation matrix
    A_rot = [ cos(Floes(ii).alpha_i) -sin(Floes(ii).alpha_i)
              sin(Floes(ii).alpha_i)  cos(Floes(ii).alpha_i) ];

    % store intact bond positions
    if ~isempty(intactBnds) && ~isempty(cat(1,intactBnds.Xc))
        pos = A_rot * [cat(1,intactBnds.Xc) cat(1,intactBnds.Yc)]';
        Xbond = [Xbond; pos(1,:)' + Xi(ii)];
        Ybond = [Ybond; pos(2,:)' + Yi(ii)];
    end

    % store broken bond positions
    if ~isempty(brokenBnds) && ~isempty(cat(1,brokenBnds.Xc))
        pos = A_rot * [cat(1,brokenBnds.Xc) cat(1,brokenBnds.Yc)]';
        Xbroken = [Xbroken; pos(1,:)' + Xi(ii)];
        Ybroken = [Ybroken; pos(2,:)' + Yi(ii)];
    end

    % store bond number for intact bonds to plot
    BondNum = unique([intactBnds.Num]);

    for jj = 1:length(BondNum)
        f1(count) = Floes(ii).num;
        f2(count) = BondNum(jj);
        count = count + 1;
    end
end


clear poly
for iii = 1:length(Floes)
    poly(iii) = polyshape(Floes(iii).c_alpha'+[Floes(iii).Xi Floes(iii).Yi]);
end

count = 1;
clear p

axis manual

plot(Xbond,Ybond,'r.','linewidth',3)
plot(Xbroken,Ybroken,'bx','linewidth',2)
hold on

plot(poly,'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2,'linestyle','none');
if Nbound>0
    plot(poly_bound,'FaceColor',[0 0.2 0],'FaceAlpha',1,'EdgeAlpha',0.5);
end

set(0,'CurrentFigure',fig);

colormap('gray'); caxis([0 1]);
axis([-Lxmax Lxmax  -Lymax Lymax])

set(gca,'Ydir','normal');

end
