function Floes=fracture_floe_subfloes(floe,Nsubfloes,Floe0,polynew)
%%This function takes an input of floes and fractures each floe into a
%%specified number of smaller ones using Voronoi Tesselation
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Floes=[]; rho_ice = 920; 
% save('floefail.mat','Floe','N','Floe0')
[~,~,Nz] = size(Floe0(1).StressH);

a = floe.interactions;

% %Permanently deform shapes during fracture
% XInew = [];
% if ~isempty(a)
%     k1 = isinf(a(:,1));
%     a(k1,:) = [];
% end
% if ~isempty(a)
%     [~,k] = max(a(:,7));
%     X1new = floe.c_alpha(1,:)+floe.Xi;
%     Y1new = floe.c_alpha(2,:)+floe.Yi;
%     clip = a(k,1);
%     if clip < length(Floe0)+1
%         X2new = Floe0(clip).c_alpha(1,:)+Floe0(clip).Xi;
%         Y2new = Floe0(clip).c_alpha(2,:)+Floe0(clip).Yi;
%         [XInew,YInew] = polyclip([X1new' Y1new'],[X2new' Y2new'],'int');
%     end
% end
% if ~isempty(XInew)
%     Xt = XInew{1}; Yt = YInew{1};
%     [xm,ym] = centroid(polyshape(Xt,Yt));
%     [d_min1] = p_poly_dist(xm, ym,Xt', Yt');
%     F = vecnorm(a(k,2:3));
%     xs = a(k,2)*abs(d_min1)/2/F;
%     ys = a(k,3)*abs(d_min1)/2/F;
%     X2new = X2new + xs; Y2new = Y2new+ys;
%     [XInew,YInew] = polyclip([X1new' Y1new'],[X2new' Y2new'],'dif');
%     if ~isempty(XInew)
%         Xt = XInew{1}; Yt = YInew{1};
%         Anew = polyarea(Xt,Yt);
%         if Anew/floe.area > 0.9
%             [xm,ym] = centroid(polyshape(Xt,Yt));
%             floe.Xi = xm; floe.Yi = ym;
%             floe.area = polyarea(Xt,Yt);
%             floe.c_alpha =[ Xt'-xm;Yt'-ym];
%         end
%     end
% end

for i =1:length(polynew)
    a=regions(polynew(i));
    
    %%Loop through all the new shapes to calculate the new properties of
    %%each
    for p=1:length(a)
        FloeNEW.poly = rmholes(a(p));
        [Xi,Yi] = centroid(FloeNEW.poly);
        FloeNEW.area = area(FloeNEW.poly);
        FloeNEW.mass = floe.mass*area(a(p))/floe.area;
        FloeNEW.h = floe.mass*area(a(p))/(rho_ice*FloeNEW.area*floe.area);
        FloeNEW.c_alpha = [(FloeNEW.poly.Vertices-[Xi Yi])' [FloeNEW.poly.Vertices(1,1)-Xi; FloeNEW.poly.Vertices(1,2)-Yi]];
        FloeNEW.c0 = FloeNEW.c_alpha;
        FloeNEW.inertia_moment = PolygonMoments(FloeNEW.c0',FloeNEW.h);

        FloeNEW.angles = polyangles(FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
        FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
        
        FloeNEW.strain = floe.strain;
        FloeNEW.Stress = zeros(2);
        FloeNEW.StressH = zeros(2,2,Nz);
        FloeNEW.StressCount = 1;
        FloeNEW.Fx = floe.Fx; FloeNEW.Fy = floe.Fy;
        FloeNEW.FxOA = 0; FloeNEW.FyOA = 0; FloeNEW.torqueOA = 0;
        
        err = 1;
        while err > 0.1
            FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1);
            FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1);
            FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
            err = (sum(FloeNEW.A)/1000*4*FloeNEW.rmax^2-FloeNEW.area)/FloeNEW.area;
        end
        
        FloeNEW.Xi = floe.Xi+Xi; FloeNEW.Yi = floe.Yi+Yi; FloeNEW.alive = 1;
        FloeNEW.alpha_i = 0; FloeNEW.Ui = floe.Ui; FloeNEW.Vi = floe.Vi;
        FloeNEW.dXi_p = floe.dXi_p; FloeNEW.dYi_p = floe.dYi_p;
        FloeNEW.dUi_p = floe.dUi_p; FloeNEW.dVi_p = floe.dVi_p;
        FloeNEW.dalpha_i_p = 0; FloeNEW.ksi_ice = floe.ksi_ice;
        FloeNEW.dksi_ice_p = floe.dksi_ice_p;
        FloeNEW.interactions = [];
        FloeNEW.potentialInteractions = [];
        FloeNEW.collision_force = 0;
        FloeNEW.collision_torque = 0;
        FloeNEW.OverlapArea = 0;
        FloeNEW.Fx = floe.Fx*area(a(p))/floe.area;
        FloeNEW.Fy = floe.Fy*area(a(p))/floe.area;
        
        polyOrig = polyshape(FloeNEW.c_alpha');
        x = [min(FloeNEW.c_alpha(1,:)) max(FloeNEW.c_alpha(1,:))]; dx = x(2)-x(1);
        y = [min(FloeNEW.c_alpha(2,:)) max(FloeNEW.c_alpha(2,:))]; dy = y(2) -y(1);
        XX = 0.975*dx/2*(2*rand(ceil(Nsubfloes*dx*dy/FloeNEW.area),1)-1)+(x(2)+x(1))/2;
        YY = 0.975*dy/2*(2*rand(ceil(Nsubfloes*dx*dy/FloeNEW.area),1)-1)+(y(2)+y(1))/2;
        in = inpolygon(XX,YY,FloeNEW.c_alpha(1,:)',FloeNEW.c_alpha(2,:)');
        Ys = YY(in); Xs = XX(in);
        [d_min1] = p_poly_dist(Xs, Ys,FloeNEW.c_alpha(1,:)', FloeNEW.c_alpha(2,:)');
        Ys = Ys(abs(d_min1)>FloeNEW.rmax/10); Xs = Xs(abs(d_min1)>FloeNEW.rmax/10);
        Afloes = 1; c2 = 1;
        while Afloes > 0.01
            nmin = 3;count = 1; NYs = 1;
            while NYs<nmin
                XX = 0.975*dx/2*(2*rand(ceil(Nsubfloes*dx*dy/FloeNEW.area),1)-1)+(x(2)+x(1))/2;
                YY = 0.975*dy/2*(2*rand(ceil(Nsubfloes*dx*dy/FloeNEW.area),1)-1)+(y(2)+y(1))/2;
                in = inpolygon(XX,YY,FloeNEW.c_alpha(1,:)',FloeNEW.c_alpha(2,:)');
                Ys = YY(in); Xs = XX(in);
                [d_min1] = p_poly_dist(Xs, Ys,FloeNEW.c_alpha(1,:)', FloeNEW.c_alpha(2,:)');
                Ys = Ys(abs(d_min1)>FloeNEW.rmax/(10^count)); Xs = Xs(abs(d_min1)>FloeNEW.rmax/(10^count));
                NYs = length(Ys);
                count = count+1;
            end
            [~, b,~,~,~] = polybnd_voronoi([Xs Ys],FloeNEW.c_alpha');
            clear nanb
            for kk = 1:length(b)
                nanb(kk) = isnan(max(max(b{kk})));
            end
            b(nanb == 1) = []; Xs(nanb == 1) = []; Ys(nanb == 1) = [];
            clear poly
            for jj = 1:length(b); poly(jj) = polyshape(b{jj}); poly(jj) = intersect(poly(jj),polyOrig); end
            areas = area(poly); Xs=Xs(areas>0); Ys=Ys(areas>0); poly = poly(areas>0); b = b(areas>0);
            areas = areas(areas>0);
            Afloes = sum(areas)/FloeNEW.area-1;
            c2 = c2+1;
            if c2 > 10
                xx = 1; xx(1) =[1 2];
            end
        end
        if sum(areas)/FloeNEW.area-1>0.01
            xx = 1; xx(1) =[1 2];
        end
        FloeNEW.bonds.Xs = Xs; FloeNEW.bonds.Ys = Ys;
        for kk = 1:length(areas)
            FloeNEW.bonds.Vert{kk} = poly(kk).Vertices;
        end
%         FloeNEW.bonds.Vert = b';
        for kk = 1:length(b)
            FloeNEW.bonds.Num{kk,1} = [];
            FloeNEW.bonds.d{kk,1} = [];
%            FloeNEW.bonds.Theta{kk,1} = [];
%             poly = polyshape(b{kk}); polyNEW = subtract(polyOrig,poly);
%             FloeNEW.bonds.In(kk,1) = logical(polyNEW.NumHoles);
        end
        for kk = 1:length(b)-1
            for jj = kk+1:length(b)
                [k,dist] = dsearchn(b{kk},b{jj});
                if length(k(dist<1)) > 1
                    verts = b{ii}; verts(k(dist<1),:);
                    d = sqrt((verts(1,1)-verts(2,1))^2+(verts(1,2)-verts(2,2))^2);
                    FloeNEW.bonds.Num{ii} = [FloeNEW.bonds.Num{ii} jj];
                    FloeNEW.bonds.d{ii} = [FloeNEW.bonds.d{ii} d];
                    %FloeNEW.bonds.Theta{ii} = [FloeNEW.bonds.Theta{ii} atan((Ys(jj)-Ys(ii))/(Xs(jj)-Xs(ii)))];
                    FloeNEW.bonds.Num{jj} = [FloeNEW.bonds.Num{jj} ii];
                    FloeNEW.bonds.d{ii} = [FloeNEW.bonds.d{ii} d];
                    %FloeNEW.bonds.Theta{jj} = [FloeNEW.bonds.Theta{jj} atan((Ys(jj)-Ys(ii))/(Xs(jj)-Xs(ii)))];
                end
            end
        end

        
        Floes = [Floes FloeNEW];
        clear FloeNEW
    end
    
end

if sum(cat(1,Floes.area))/floe.area-1 > 0.1
    xx = 1;
    xx(1) =[1 2];
end

if ~isempty(Floes)
    Floes=rmfield(Floes,{'poly'});
end

warning('on',id)
warning('on',id3)
end
