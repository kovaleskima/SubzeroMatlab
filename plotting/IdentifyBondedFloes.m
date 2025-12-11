close all
FloeNumbers = 1:length(Floe);
Floes = Floe(Nbound+1:end);
Nums = cat(1,Floe.num);
clear f1; clear f2
count = 1;
for ii = 1:length(Floes)
    BondNum = unique(cat(1,Floes(ii).bonds.Num));
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

close all
figure
hold on
colors = distinguishable_colors(length(p));
for ii = 1:length(p)
    plot(p(ii),'FaceColor',colors(ii,:) )
end

Abox = area(c2_boundary_poly);
a = sqrt(area(p));
hFSD = histogram(a/1e3);
edges = hFSD.BinEdges;
val = hFSD.Values;
clear FSD
count3 = 1;
%Loop through to find floe sizes per km^2
for ii = 1:length(val)
    FSD(count3) = sum(val(ii:end))/Abox*1e6;
    count3 = count3+1;
end
%Plot FSD
bins = (edges(2:end)+edges(1:end-1))/2;
bins(FSD<1e-4) = []; FSD(FSD<1e-4) = [];
figure;
fig(1) = loglog(bins,FSD,'linewidth',2);
min_size = 1;
binsUpper = bins(bins>min_size); slopes1 = -2;
hold on
fig(2) = loglog([binsUpper(1) 50], 5*1e-1*[binsUpper(1) 50].^(slopes1),'k--','linewidth',2);
set(gca, 'YScale', 'log')
set(gca, 'xScale', 'log')
xlabel('Floe Size (km)','fontsize',20,'interpreter','latex')
ylabel({'FSD (floes per km$^2$)'},'fontsize',20,'interpreter','latex')
fig = figure(2);

%% 
fig = 0;
for ii = 1:427
    load(['./Floes_bnds1/Floe' num2str(ii,'%07.f') '.mat'])
    c2_boundary_poly = polyshape(c2_boundary');
    [fig] =plot_basic_bonds(fig,Floe,ocean,c2_boundary_poly,Nb,Nbond,PERIODIC);
    exportgraphics(fig,['./figs/' num2str(ii,'%03.f') '.jpg'] ,'resolution',300);
end
