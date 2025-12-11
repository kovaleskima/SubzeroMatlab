kk = 1;
fig = 0; num_chunks = zeros(1,427);Nbound = 0;edges = 0:2:102;
load(['/dat1/bpm5026/Kennedy/Floes_bnds5/Floe' num2str(kk,'%07.f') '.mat'])
c2_boundary_poly = polyshape(c2_boundary'); fig = 0; %Nbound = Nb;
%[fig] =plot_basic_bonds(fig,Floe,ocean,c2_boundary_poly,Nb,Nbond,PERIODIC);

close all
FloeNumbers = 1:length(Floe);
Floes = Floe;
Nums = cat(1,Floe.num);
clear f1; clear f2
Xc = []; Xb = [];
Yc = []; Yb = [];
Xi = cat(1,Floe.Xi); Yi = cat(1,Floe.Yi);
count = 1;
for ii = 1:length(Floes)
    BondNum = [FloeNumbers(ii); unique(cat(1,Floes(ii).bonds.Num))];
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
plot(G,'Layout','layered')
[bins,binsizes] = conncomp(G,'Type','weak');
Lia = ismember(FloeNumbers,Nums);Lia(1:Nbond) = 0;
clear poly
for iii = 1:length(Floe)
    poly(iii) = polyshape(Floe(iii).c_alpha'+[Floe(iii).Xi Floe(iii).Yi]);
end
count = 1;
p = polyshape(); 
for ii = 1:length(binsizes)
    Lia1 = ismember(bins,ii);
    Lia2 = ismember(Nums,FloeNumbers(Lia1));
    if sum(Lia2)> 0
        ptmp = union([poly(Lia2)]);
        p(count) = ptmp;%rmholes(ptmp);
        count = count+1;
    end
end
%num_chunks(kk) = length(p);

close all
ratio=2.1;
fig=figure('Position',[100 100 1000*ratio 1000],'visible','on');
set(fig,'PaperSize',12*[ratio 1],'PaperPosition',12*[0 0 ratio 1]);
figure(fig)
clf(fig);
subplot(1,2,1)
plot(poly(1:Nb),'FaceColor',[0 0.4 0],'FaceAlpha',1,'EdgeAlpha',0.5);
hold on
plot(poly(1+Nb:length(Floe)),'FaceColor','none','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
plot(Xc,Yc,'r.','linewidth',3);
colors = distinguishable_colors(length(p));
for ii = 1:length(p)
    plot(p(ii),'FaceColor',colors(ii,:) )
end

set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
plot(xb,yb, 'k-','linewidth',2);
colormap('gray'); caxis([0 1]);
axis equal
axis([-65000 65000 -65000 65000])
set(gca,'Ydir','normal');
box on

subplot(1,2,2)
plot(poly(1:Nb),'FaceColor',[0 0.2 0],'FaceAlpha',1,'EdgeAlpha',0.5);
hold on
plot(p,'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
plot(xb,yb, 'k-','linewidth',2);
colormap('gray'); caxis([0 1]);
axis equal
axis([-65000 65000 -65000 65000])
set(gca,'Ydir','normal');
box on
fig = figure(1);
exportgraphics(fig,['./plotting/bonds_colors.jpg'] ,'resolution',300);  
