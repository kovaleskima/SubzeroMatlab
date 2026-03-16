function [Floe,bonds,Nbound,Nbond] = initial_concentration(c2_boundary,target_concentration,height, NumFloes, NumSubfloes, min_floe_size)
%% This function is used to generate the initial floe field

%Identify the grids to align with the concentrations specified by the input
[Ny, Nx] = size(target_concentration);
c = flipud(target_concentration);
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
dx = x(2)-x(1);
dy = y(2)-y(1);

%Create floes that act as boundaries and wont move
% x1 = [Lx/2 Lx/2 1.5e4 1.5e4]; y1 = [-Ly/2 -1e4 -Ly/2+3e5 -Ly/2];
% B1 = polyshape(x1, y1);
% B2 = polyshape(-x1, y1);
% x2 = [Lx/2 Lx/2 -Lx/2 -Lx/2]; y2 = [-1e5 -Ly/2 -Ly/2 -1e5];
% B3 = polyshape(x2,y2);
% Floe1 = initialize_floe_values(B1,height,0,SUBFLOES);
% Floe2 = initialize_floe_values(B2,height,0,SUBFLOES);
% bound = subtract(c2_boundary_poly, B1);
% bound = subtract(bound, B2);
%Floe = [Floe1 Floe2];

nx=40; ny=4;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;
[xx,yy] = meshgrid(Xc,Yc); X = xx(:); Y = yy(:);%[yc(1) yc(2) yc(2) yc(1)];

Nbound = 0;
Floe = [];
FloeBnds = [];
FloeNum = 1;

%Loop through all the regions of the domain to create new floes
for jj = 1:Ny
    for ii = 1:Nx
        if c(jj,ii)>0
            boundary = polyshape([x(ii) x(ii) x(ii+1) x(ii+1)], [y(jj) y(jj+1) y(jj+1) y(jj)]);
            boundary = intersect(boundary,c2_boundary_poly);
            N = ceil(NumFloes*area(boundary)/area(c2_boundary_poly)/c(jj,ii));
%             poly = intersect(bound,boundary); %Use these when having
%             boundaries
%             N = 4*ceil(NumFloes*area(poly)/area(bound)/c(jj,ii)); %Use these when having
%             boundaries
            if N == 1
                Floe = initialize_floe_values(c2_boundary_poly,height,NumSubfloes);
            else
                X = 0.975*dx/2*(2*rand(N,1)-1)+(x(ii)+x(ii+1))/2;
                Y = 0.975*dy/2*(2*rand(N,1)-1)+(y(jj)+y(jj+1))/2;
                in = inpolygon(X,Y,boundary.Vertices(:,1),boundary.Vertices(:,2));
                X = X(in); Y = Y(in);
                %             for i = 1:2%nx
                %                 X = [(-1)^i*dx/2 (-1)^i*dx/2 0 0];%xx(:,i:i+1); X = X(:);
                %                 Y = [-dy/2 dy/2 dy/2 -dy/2];%Y = yy(:,i:i+1); Y = Y(:);
                %                 b{i} = [X',Y'];
                %             end
                [~, b,~,~,~] = polybnd_voronoi([X Y],boundary.Vertices);
                %             [~, b,~,~,~] = polybnd_voronoi([X Y],boundary.Vertices); %%Use these when having
                %             boundaries
                
                for iii = 1:length(b)
                    bonds(iii).bond = [];
                end
                r = sqrt(min_floe_size);
                for m = 1:length(b)
                    poly(m) = polyshape(b{m});
                end
%                 xx = 1; xx(1) =[1 2];
%                 load('./Initialize_Model/polystart.mat','poly');
%                load('./Initialize_Model/pstart.mat','poly');
%                 poly = polys;
                Nf = 1:length(poly);%randperm(length(b));
                for iii = 1:length(poly)-1
                    for jjj = iii+1:length(poly)
                        [k,dist] = dsearchn(poly(iii).Vertices,poly(jjj).Vertices);
                        if length(k(dist<0.1)) > 1
                            verts = poly(iii).Vertices; verts = verts(k(dist<0.1),:);
                            Xi = (verts(1,1)+verts(2,1))/2; Yi = (verts(1,2)+verts(2,2))/2;
                            Xi1 = (verts(1,1)+Xi)/2; Yi1 = (verts(1,2)+Yi)/2;
                            Xi2 = (verts(2,1)+Xi)/2; Yi2 = (verts(2,2)+Yi)/2;
                            d = sqrt((verts(1,1)-verts(2,1))^2+(verts(1,2)-verts(2,2))^2);
                            dx=(verts(1,1)-verts(2,1))/d; dy = (verts(1,2)-verts(2,2))/d;
                            bonds_i.Num = jjj; bonds_j.Num = iii;
                            vertsy = [Yi1, Yi1+r/2*(dy*sqrt(3)-dx),Yi1+r/2*(-dy*sqrt(3)-dx)];
                            vertsx = [Xi1, Xi1+r/2*(dx*sqrt(3)+dy),Xi1+r/2*(-dx*sqrt(3)+dy)];
                            vertsy2 = [Yi1, Yi1+r/2*(dy*sqrt(3)+dx),Yi1+r/2*(-dy*sqrt(3)+dx)];
                            vertsx2 = [Xi1, Xi1+r/2*(dx*sqrt(3)-dy),Xi1+r/2*(-dx*sqrt(3)-dy)];
                            poly_bnd1 = polyshape(vertsx,vertsy);
                            poly_bnd1 = translate(poly_bnd1,-r/10*[dy -dx]);
                            poly_bnd2 = polyshape(vertsx2,vertsy2);
                            poly_bnd2 = translate(poly_bnd2,r/10*[dy -dx]);
                            vertsy = [Yi2, Yi2+r/2*(dy*sqrt(3)-dx),Yi2+r/2*(-dy*sqrt(3)-dx)];
                            vertsx = [Xi2, Xi2+r/2*(dx*sqrt(3)+dy),Xi2+r/2*(-dx*sqrt(3)+dy)];
                            vertsy2 = [Yi2, Yi2+r/2*(dy*sqrt(3)+dx),Yi2+r/2*(-dy*sqrt(3)+dx)];
                            vertsx2 = [Xi2, Xi2+r/2*(dx*sqrt(3)-dy),Xi2+r/2*(-dx*sqrt(3)-dy)];
                            poly_bnd3 = polyshape(vertsx,vertsy);
                            poly_bnd3 = translate(poly_bnd3,-r/10*[dy -dx]);
                            poly_bnd4 = polyshape(vertsx2,vertsy2);
                            poly_bnd4 = translate(poly_bnd4,r/10*[dy -dx]);
                            poly_bnd = union(poly_bnd1,poly_bnd2);
                            if area(intersect(poly_bnd,poly(iii)))>area(intersect(poly_bnd,poly(jjj)))
                                bonds_i.c(1) = poly_bnd1; bonds_i.c(2) = poly_bnd3;
                                bonds_j.c(1) = poly_bnd2; bonds_j.c(2) = poly_bnd4;
                                bonds_i.r_bnd = sqrt(area(poly_bnd1)); bonds_j.r_bnd = sqrt(area(poly_bnd2));
                            else
                                bonds_i.c(1) = poly_bnd2; bonds_i.c(2) = poly_bnd4;
                                bonds_j.c(1) = poly_bnd1; bonds_j.c(2) = poly_bnd3;
                                bonds_i.r_bnd = sqrt(area(poly_bnd2)); bonds_j.r_bnd = sqrt(area(poly_bnd1));
                            end
%                                 bonds_i.c = poly_bnd; %bonds_i.piece = 0;
%                                 bonds_j.c = poly_bnd; %bonds_j.piece =  1;
                            bonds_i.L = d; bonds_j.L = d;
                            bonds_i.Xb = []; bonds_j.Xb = [];
                            bonds_i.Yb = []; bonds_j.Yb = [];
                            bonds_i.Xc(1) = Xi1; bonds_j.Xc(1) = Xi1;
                            bonds_i.Xc(2) = Xi2; bonds_j.Xc(2) = Xi2;
                            bonds_i.Yc(1) = Yi1; bonds_j.Yc(1) = Yi1;
                            bonds_i.Yc(2) = Yi2; bonds_j.Yc(2) = Yi2;
                            bonds_i.Fx_p = 0; bonds_j.Fx_p = 0;
                            bonds_i.Fy_p = 0; bonds_j.Fy_p = 0;
                            bonds_i.Stress = [0;0]; bonds_j.Stress = [0;0];
                            bonds(iii).bond = [bonds(iii).bond bonds_i];
                            bonds(jjj).bond = [bonds(jjj).bond bonds_j];
%                             xx = 1; xx(1)=[1 2];
%                             poly_bond = union(poly_bnd,poly_bnd2);
%                             [xb,yb] = centroid(poly_bond); pbnd2 = translate(poly_bond,-[xb,yb]);
%                             pbnd2 = scale(pbnd2,1); poly_bond= translate(pbnd2,[xb,yb]);
%                             p1 = subtract(poly_bond,poly(iii)); a1 = area(intersect(poly_bond,poly(iii)));
%                             p2 = subtract(poly_bond,poly(jjj)); a2 = area(intersect(poly_bond,poly(iii)));
%                             a = area(poly_bond);
%                             if (a1 + a2)/a > 0.999 && p1.NumRegions == 1 && p2.NumRegions == 1
%                                 floe_bnd = initialize_floe_values(poly_bond, height,1);
%                                 floe_bnd.num = FloeNum; 
%                                 floe_bnd.bonds.bowtie = 1;
%                                 floe_bnd.bonds.L = d;
%                                 floe_bnd.bonds.FloeNum = [iii jjj];
%                                 FloeBnds = [FloeBnds floe_bnd];
%                                 bonds.poly{iii} = [bonds.poly{iii} poly_bond];
%                                 bonds.poly{jjj} = [bonds.poly{jjj} poly_bond];
%                                 bonds.BndNum{iii} = [bonds.BndNum{iii} FloeNum];
%                                 bonds.BndNum{jjj} = [bonds.BndNum{jjj} FloeNum];
%                                 FloeNum = FloeNum+1;
%                             end
                        end
                    end
                end
%                 bnd_size = FloeBnds(1).area; 
%                 for iii = 1:length(FloeBnds)
%                     p = intersect(FloeBnds(iii).poly,[FloeBnds.poly]);
%                     a = area(p)/bnd_size; a(iii) = 0;
%                     if sum(a)>0.25
%                         xx = 1; xx(1) =[1 2];
%                     end
%                 end
                Atot = 0;
                count = 1;
                polyFloes = [];
%                 while Atot/area(boundary)<=c(jj,ii)
%                     if ~isnan(b{Nf(count)})
%                         pnew = poly(Nf(count));
%                         polyFloes = [polyFloes pnew];
%                         Atot = Atot+area(polyFloes(count));
%                     end
%                     count = count+1;
%                     if count > length(Nf)
%                         Atot = area(boundary)+1;
%                     end
%                 end
                polyFloes = poly;
                for kk= 1:length(polyFloes)
                    polyFloesNew = polyFloes(kk);
                    bond_tmp = bonds(kk).bond;
                    new_bnds = [];
%                     for iii=1:length(bond_tmp)
%                             hole = intersect(bond_tmp(iii).c,polyFloesNew);
%                             bond_tmp(iii).hole = hole.Vertices;
% %                         if bond_tmp(iii).piece
% %                             bond_tmp(iii).c = subtract(bond_tmp(iii).c,polyFloesNew);
% %                         else
%                             polyFloesNew = subtract(polyFloesNew,bond_tmp(iii).c);
% %                             bond_tmp(iii).c = intersect(bond_tmp(iii).c,polyFloes(kk));
% %                             bond_tmp(iii).c = intersect(bond_tmp(iii).c,polyFloes(kk));
%                             if polyFloesNew.NumHoles > 0
%                                 [k,dist] = dsearchn(bond_tmp(iii).c.Vertices,polyFloesNew.Vertices);
%                                 xx = 1; xx(1) =[1 2];
%                                 for jjj = 1:length(dist)
%                                     if dist(jjj) < 1
%                                         polyFloesNew.Vertices(jjj,:) = bond_tmp(iii).c.Vertices(k(jjj),:);
%                                     end
%                                 end
%                             end
% %                         end
%                     end
                    if polyFloesNew.NumRegions > 1
                        polyout = sortregions(polyFloesNew,'area','descend');
                        R = regions(polyout);
                        polyFloesNew = R(1);
                    end
                    floenew = initialize_floe_values(polyFloesNew,height,NumSubfloes);
                    floenew.num = Nf(kk);
                    if ~isempty(bond_tmp)
                        for iii = 1:length(bond_tmp)
                            for jjj = 1:length(bond_tmp(iii).c)
                                bond_tmp2 = bond_tmp(iii);
                                [Xb,Yb] = centroid(bond_tmp(iii).c(jjj));  bond_tmp2.Xb = Xb-floenew.Xi; bond_tmp2.Yb = Yb-floenew.Yi;
                                bond_tmp2.Xc = bond_tmp(iii).Xc(jjj)-floenew.Xi; bond_tmp2.Yc = bond_tmp(iii).Yc(jjj)-floenew.Yi;
                                new_bnds = [new_bnds; bond_tmp2];
                            end
%                             Verts = [bond_tmp(iii).c.Vertices(:,1)-floenew.Xi bond_tmp(iii).c.Vertices(:,2)-floenew.Yi];
%                             bond_tmp(iii).c = Verts;
                        end
                        new_bnds=rmfield(new_bnds,{'c'});
%                         floenew.bonds = bonds(kk);
%                         bonds(kk).poly = Verts;%[Xvert Yvert];%pbnd.Vertices;
%                         clear R
%                         for iii = 1:length(bond_tmp)
%                             R(iii) = polyshape(bond_tmp(iii).c);%regions(polyshape(floenew.bonds.poly));
%                         end
%                         polynew = polyshape(floenew.c0');
%                         poly = union(polynew,union(R)); %Fix that this is not ii but needs to be FloeNumber field
%                         if poly.NumRegions > 1
%                             for jjj = 1:length(R)
%                                 [k,dist] = dsearchn(R(jjj).Vertices,polynew.Vertices);
%                                 for iii = 1:length(dist)
%                                     if dist(iii) < 1
%                                         polynew.Vertices(iii,:) = R(jjj).Vertices(k(iii),:);
%                                     end
%                                 end
%                             end
%                             poly = union(polynew,union(R)); %Fix that this is not ii but needs to be FloeNumber field
%                             if poly.NumRegions > 1
%                                 xx = 1; xx(1) =[1 2];
%                             end
%                         end
                        floenew.bonds = new_bnds;
                    else
                        xx = 1; xx(1) =[1 2];
                    end
                    if abs(Nf(kk)-kk)>0
                        xx = 1; xx(1) =[1 2];
                    end
                    Floe = [Floe floenew];
                end
            end
        end
    end
end
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

% for ii = 1:length(bonds)
%     for jj = 1:length(bonds(ii).bond)
% %         if bonds(ii).bond(jj).piece
%             bonds(ii).bond(jj).h = Floe(ii).h;
% %         else
% %             bonds(ii).bond(jj).h = Floe(bonds(ii).bond(jj).Num).h;
% %         end
%     end
% end


% areas = cat(1,Floe.area);
% small = Nums(areas<min_floe_size);
% for ii = 1:length(small)
%     bnds = unique(cat(1,Floe(small(ii)).bonds.Num));
%     for jj= length(bnds):-1:1
%         Lia = ismember(Nums,bnds(jj));
%         BndNums = cat(1,Floe(Lia).bonds.Num);
%         Lia2 = ismember(BndNums,Floe(small(ii)).num);
%         Floe(Lia).bonds(Lia2) = [];
%     end
% end
% Floe(areas<min_floe_size)=[];
Nbond = 0;



% Nbond = length(FloeBnds); 
% for ii = 1:Nbond
%     FloeBnds(ii).bonds.FloeNum = FloeBnds(ii).bonds.FloeNum + Nbond;
% end
% Floe = [FloeBnds Floe];
% 
% for ii = 1:Nbond
%     Fnum = Floe(ii).bonds.FloeNum;
%     Nums = cat(1,Floe.num);
%     Numbers = [1:length(Floe)]';
%     for jj = 1:length(Fnum)
%         kk = Numbers(Nums==Fnum(jj));
%         BondNum = Floe(kk).bonds.BondNum;
%         clear R
%         for iii = 1:length(Floe(kk).bonds.poly)
%             R(iii) = polyshape(Floe(kk).bonds.poly{iii});%regions(polyshape(floenew.bonds.poly));
%         end   
%         polynew = polyshape(Floe(kk).c0');
%         poly = union(polynew,R(BondNum==Floe(ii).num)); %Fix that this is not ii but needs to be FloeNumber field
%         if poly.NumRegions > 1
%             [k,dist] = dsearchn(R(BondNum==Floe(ii).num).Vertices,polynew.Vertices);
%             for iii = 1:length(dist)
%                 if dist(iii) < 1
%                     polynew.Vertices(iii,:) = R(BondNum==Floe(ii).num).Vertices(k(iii),:);
%                 end
%             end
%             poly = union(polynew,R(BondNum==Floe(ii).num)); %Fix that this is not ii but needs to be FloeNumber field
%             if poly.NumRegions > 1
%                 xx = 1; xx(1) =[1 2];
%             end
%         end
%     end
% end
end

