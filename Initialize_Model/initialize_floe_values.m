function FloeNEW = initialize_floe_values(poly1, height,Nsubfloes)
%%This function populates all the fields of the floes upon their creation

rho_ice=920;
dX = 480;

polyout = sortregions(poly1,'area','descend');
R = regions(polyout);
poly1 = R(1);
polya = rmholes(poly1);
h=height.mean+(-1)^randi([0 1])*rand*height.delta;

FloeNEW.poly = poly1;
FloeNEW.num = 0;
[Xi,Yi] = centroid(FloeNEW.poly);
FloeNEW.area = area(FloeNEW.poly);
FloeNEW.h = h;
FloeNEW.mass = FloeNEW.area*h*rho_ice;
FloeNEW.c_alpha = [(polya.Vertices-[Xi Yi])' [polya.Vertices(1,1)-Xi; polya.Vertices(1,2)-Yi]];
FloeNEW.c0 = FloeNEW.c_alpha;
FloeNEW.inertia_moment = PolygonMoments(FloeNEW.c0',h);
FloeNEW.angles = polyangles(polya.Vertices(:,1),polya.Vertices(:,2));
FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
% n=(fix(FloeNEW.rmax/dX)+1); n=dX*(-n:n);
% FloeNEW.Xg = n;
% FloeNEW.Yg = n;
% [X, Y]= meshgrid(n, n);
% FloeNEW.X = X;
% FloeNEW.Y = Y;
FloeNEW.Stress = [0 0; 0 0];
FloeNEW.strain = [0 0; 0 0];
FloeNEW.StressH = zeros(2,2,10);
FloeNEW.StressCount = 1;
%FloeNEW.strainCount = 0;
%FloeNEW.MaxShear = 0;
FloeNEW.Fx = 0; FloeNEW.Fy = 0;
FloeNEW.FxOA = 0; FloeNEW.FyOA = 0; FloeNEW.torqueOA = 0;
FloeNEW.alpha_i = 0; FloeNEW.Ui = 0; FloeNEW.Vi = 0;

err = 1;
count = 1;
while err > 0.1
    FloeNEW.X = FloeNEW.rmax*(2*rand(50,1) - 1);
    FloeNEW.Y = FloeNEW.rmax*(2*rand(50,1) - 1);
    FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
    err = abs((sum(FloeNEW.A)/50*4*FloeNEW.rmax^2-area(polya)))/area(polya);
    count = count+1; if count>10; err = 0; FloeNEW.alive = 0; end
end
% [in] = inpolygon(FloeNEW.X(:)+Xi, FloeNEW.Y(:)+Yi,FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
% FloeNEW.A=reshape(in,length(FloeNEW.X),length(FloeNEW.X));
if Nsubfloes > 1
    polyOrig = polyshape(FloeNEW.c_alpha');
    x = [min(FloeNEW.c_alpha(1,:)) max(FloeNEW.c_alpha(1,:))]; dx = x(2)-x(1);
    y = [min(FloeNEW.c_alpha(2,:)) max(FloeNEW.c_alpha(2,:))]; dy = y(2) -y(1);
    
    inn = 1;
    while inn < 3
        XX = 0.975*dx/2*(2*rand(ceil(Nsubfloes*dx*dy/FloeNEW.area),1)-1)+(x(2)+x(1))/2;
        YY = 0.975*dy/2*(2*rand(ceil(Nsubfloes*dx*dy/FloeNEW.area),1)-1)+(y(2)+y(1))/2;
        % X = 0.975*FloeNEW.rmax*(2*rand(ceil(Nsubfloes*FloeNEW.rmax^2/FloeNEW.area),1)-1);
        % Y = 0.975*FloeNEW.rmax*(2*rand(ceil(Nsubfloes*FloeNEW.rmax^2/FloeNEW.area),1)-1);
        in = inpolygon(XX,YY,FloeNEW.c_alpha(1,:)',FloeNEW.c_alpha(2,:)');
        Ys = YY(in); Xs = XX(in);
        [d_min1] = p_poly_dist(Xs, Ys,FloeNEW.c_alpha(1,:)', FloeNEW.c_alpha(2,:)');
        Ys = Ys(abs(d_min1)>FloeNEW.rmax/10); Xs = Xs(abs(d_min1)>FloeNEW.rmax/10);
        inn = length(Ys);
    end
    FloeNEW.bonds.Xs = Xs; FloeNEW.bonds.Ys = Ys;
    FloeNEW.bonds.interactions = [];
    in = inpolygon(Xs,Ys,FloeNEW.c_alpha(1,:)',FloeNEW.c_alpha(2,:)');
    if sum(in)/length(in)<1
        xx = 1; xx(1) =[1 2];
    end
    [~, b,~,~,~] = polybnd_voronoi([Xs Ys],FloeNEW.c_alpha');
    clear nanb
    for kk = 1:length(b)
        nanb(kk) = isnan(max(max(b{kk})));
    end
    b(nanb == 1) = []; Xs(nanb == 1) = []; Ys(nanb == 1) = [];
    for jj = 1:length(b); poly(jj) = polyshape(b{jj}); poly(jj) = intersect(poly(jj),polyOrig); end
    areas = area(poly); Xs=Xs(areas>0); Ys=Ys(areas>0); b = b(areas>0); 
    areas = areas(areas>0);
    for kk = 1:length(areas)
        FloeNEW.bonds.Vert{kk} = poly(kk).Vertices;
    end
    
    for ii = 1:length(b)
        FloeNEW.bonds.Num{ii,1} = [];
        FloeNEW.bonds.d{ii,1} = [];
%         FloeNEW.bonds.Theta{ii,1} = [];
%         poly = polyshape(b{ii}); polyNEW = subtract(polyOrig,poly);
%         FloeNEW.bonds.In(ii,1) = logical(polyNEW.NumHoles);
    end
%     FloeNEW.bonds.Vert = b';
    A_rot=[cos(FloeNEW.alpha_i) -sin(FloeNEW.alpha_i); sin(FloeNEW.alpha_i) cos(FloeNEW.alpha_i)]; %rotation matrix

    Subfloes = A_rot*[FloeNEW.bonds.Xs'; FloeNEW.bonds.Ys'];
    in = inpolygon(Subfloes(1,:),Subfloes(2,:),FloeNEW.c_alpha(1,:)',FloeNEW.c_alpha(2,:)');
    if sum(in)/length(in)<1
        xx = 1; xx(1) =[1 2];
    end

    for ii = 1:length(b)-1
        for jj = ii+1:length(b)
            [k,dist] = dsearchn(b{ii},b{jj});
            if length(k(dist<1)) > 1
                verts = b{ii}; verts = verts(k(dist<1),:);
                d = sqrt((verts(1,1)-verts(2,1))^2+(verts(1,2)-verts(2,2))^2);
                FloeNEW.bonds.Num{ii} = [FloeNEW.bonds.Num{ii} jj];
                FloeNEW.bonds.d{ii} = [FloeNEW.bonds.d{ii} d];
                %FloeNEW.bonds.Theta{ii} = [FloeNEW.bonds.Theta{ii} atan((Ys(jj)-Ys(ii))/(Xs(jj)-Xs(ii)))];
                FloeNEW.bonds.Num{jj} = [FloeNEW.bonds.Num{jj} ii];
                FloeNEW.bonds.d{jj} = [FloeNEW.bonds.d{jj} d];
                %FloeNEW.bonds.Theta{jj} = [FloeNEW.bonds.Theta{jj} atan((Ys(jj)-Ys(ii))/(Xs(jj)-Xs(ii)))];
            end
        end
    end
else
    FloeNEW.bonds.bowtie = 0;
    FloeNEW.bonds.L = 0;
    FloeNEW.bonds.FloeNum = 0;
    FloeNEW.bonds.BondNum = 0;
    FloeNEW.bonds.interactions = 1;
%     FloeNEW.bonds.Xs = 0; FloeNEW.bonds.Ys = 0;
%     FloeNEW.bonds.Num = 0; FloeNEW.bonds.Theta = 0;
%     FloeNEW.bonds.d = 0; FloeNEW.bonds.interactions = [];
end

if Nsubfloes > 1
    for ii = 1:length(FloeNEW.bonds.Num)
        if abs(length(FloeNEW.bonds.Num{ii})-length(FloeNEW.bonds.d{ii})) > 0
            xx = 1; xx(1) =[1 2];
        end
    end
end

FloeNEW.Xi = Xi; FloeNEW.Yi = Yi; FloeNEW.alive = 1;
FloeNEW.dXi_p = 0; FloeNEW.dYi_p = 0;
FloeNEW.dUi_p = 0; FloeNEW.dVi_p = 0;
FloeNEW.dalpha_i_p = 0; FloeNEW.ksi_ice = 0;
FloeNEW.dksi_ice_p = 0;
FloeNEW.interactions = [];
FloeNEW.potentialInteractions = [];
FloeNEW.collision_force = 0;
%         FloeNEW.fracture_force = 0;
FloeNEW.collision_torque = 0;
FloeNEW.OverlapArea = 0;



if isnan(FloeNEW.inertia_moment)
    xx = 1;
    xx(1) = [1 2];
end

end
