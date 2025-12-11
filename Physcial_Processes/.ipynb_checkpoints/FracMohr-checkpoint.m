function [Floe,Princ] = FracMohr(Floe,Nb,min_floe_size,concentration)
%Use Mohr's cone to determine which floes are fractured
%   If Principal stresses are outside the cone then the floes are fractured
rho_ice=920;

        Pstar = 2.25e5; C = 20;
        h = mean(cat(1,Floe.h));
        P = Pstar*h*exp(-C*(1-concentration));
        t = linspace(0,2*pi) ;
        a = P*sqrt(2)/2 ; b = a/2 ;
        x = [0 -1.5e5 -5.5e5];%a*cos(t) ;
        y = [0 -5.5e5 -1.5e5];%b*sin(t) ;
        Mohr = polyshape(x,y);
        Mohr = translate(Mohr, [100, 100]);
%        Mohr = rotate(Mohr,45);
%        Mohr = translate(Mohr,[-P/2, -P/2]);
        A = cat(1,Floe.area);
%        q = 5.2; SigC = 250e3;
%        Sig1 = (1/q+1)*SigC/(1/q-q);
%        Sig2 = q*Sig1+SigC;
%        Sig11 = 1e8;
%        Sig22 = q*Sig11+SigC;
%        MohrX = [Sig1; Sig11; Sig22];
%        MohrY = [Sig2; Sig22; Sig11];
        %MohrX = [Sig1; Sig11; Sig22]*0.45;*1.1
        %MohrY = [Sig2; Sig22; Sig11]*0.45;
%        Mohr = polyshape(MohrX,MohrY);
q = 5.2; SigC = 250e3;
Sig1 = (1/q+1)*SigC/(1/q-q);
Sig2 = q*Sig1+SigC;
Sig11 = -3e4;
Sig22 = q*Sig11+SigC;
MohrX = [Sig1; Sig11; Sig22];
MohrY = [Sig2; Sig22; Sig11];
Mohr = polyshape(-MohrX,-MohrY);
        p1 = -59520; p2 = 4e8;
        for ii = 1:length(Floe)
            Stress = eig(Floe(ii).Stress);
            Princ1(ii) = max(Stress); Princ2(ii) = min(Stress); 
            Princ(ii,1) = max(Stress);
            Princ(ii,2) = min(Stress);
            Princ(ii,3) = Floe(ii).area;
            Stresses(ii) = vecnorm([Princ1 Princ2]);
%            Stresses(ii) = vecnorm([Princ(ii,1:2)]);
        end
        Princ1 = Princ(:,1); Princ2 = Princ(:,2); 
        %[d_min, ~, ~] = p_poly_dist(Princ1, Princ2, Mohr.Vertices(:,1), Mohr.Vertices(:,2),true);
        [in,out] = inpolygon(Princ1,Princ2,Mohr.Vertices(:,1), Mohr.Vertices(:,2));
%         [~,x_d_min, y_d_min] = p_poly_dist(Princ1, Princ2,Mohr.Vertices(:,1), Mohr.Vertices(:,2),true);
%         [d_min2] = p_poly_dist(Princ1, Princ2,[p1; p2], [p1; p2],false);
%         [d_min3] = p_poly_dist(x_d_min, y_d_min,[p1; p2], [p1; p2],false);
count = 1; Sig1 = []; Sig2 = [];
FloeNums = cat(1,Floe.num);
for ii = 1:length(Floe)
%     A_rot=[cos(Floe(ii).alpha_i) -sin(Floe(ii).alpha_i); sin(Floe(ii).alpha_i) cos(Floe(ii).alpha_i)]; %rotation matrix
%     a = Floe(ii).interactions;
    if ~isempty(Floe(ii).bonds)
        bnds = unique(cat(1,Floe(ii).bonds.Num));
%         b = sqrt(a(:,2).^2+a(:,3).^2);
%         a(b==0,:) = []; [Ny,~] = size(a);
%         bonds = zeros(Ny,1); bonds(a(:,7)==0) = 1;
%         a(~logical(bonds),:) = [];
%         nums = a(:,1);
        for jj = length(bnds):-1:1
            Lia = ismember(cat(1,Floe(ii).bonds.Num),bnds(jj));
            bnd_tmp = Floe(ii).bonds(Lia);
            Sig1(count) = 0; Sig2(count) = 0;
            for kk = 1:length(bnd_tmp)
                Sig1(count) = Sig1(count)+bnd_tmp(kk).Stress(1)/50; Sig2(count) = Sig2(count)+bnd_tmp(kk).Stress(2)/50;
                bnd_tmp(kk).Stress = [0; 0];
            end
            Sig1(count) = Sig1(count)/kk; Sig2(count) = Sig2(count)/kk;
%                 xx= 1; xx(1) =[1 2];
                %Floe(ii).bonds(jj).Stress=[0;0];
%                 F_comp0 = [bnd_tmp.Xc-bnd_tmp.Xb; bnd_tmp.Yc-bnd_tmp.Yb];F_comp0 = F_comp0/vecnorm(F_comp0);
%                 F_comp =A_rot*F_comp0;
%                 F_bnd_dir = [a(jj,2); a(jj,3)]/sqrt(a(jj,2).^2+a(jj,3).^2);
%                 F_bnd_c = sqrt(a(jj,2).^2+a(jj,3).^2)*dot(F_bnd_dir,F_comp);
%                 F_bnd_s = sqrt(a(jj,2).^2+a(jj,3).^2-F_bnd_c.^2);
%                 Sig1(count) = F_bnd_c/(bnd_tmp.L*Floe(ii).h); Sig2(count) = F_bnd_s/(bnd_tmp.L*Floe(ii).h);
            if abs(Sig1(count))>1.65e4
                Floe(ii).bonds(Lia) = [];
                Lia = ismember(FloeNums,bnds(jj));
                BndNums = cat(1,Floe(Lia).bonds.Num);
                Lia2 = ismember(BndNums,Floe(ii).num);
                Floe(Lia).bonds(Lia2) = [];
            elseif abs(Sig2(count))>1.65e4
                Floe(ii).bonds(Lia) = [];
                Lia = ismember(FloeNums,bnds(jj));
                BndNums = cat(1,Floe(Lia).bonds.Num);
                Lia2 = ismember(BndNums,Floe(ii).num);
                Floe(Lia).bonds(Lia2) = [];
            else
                Floe(ii).bonds(Lia) = bnd_tmp;
            end
            count = count+1;
%             end
        end
    end
end
%         xx = 1; xx(1) =[1 2];
        %keep=rand(length(Floe),1)>d_min;
%         keep=rand(length(Floe),1)>0.1*(d_min2./d_min3).^2;
        p = rand(length(Floe),1);
        StressNorm = Stresses/max(Stresses);%rmoutliers(Stresses));
%        in(p>StressNorm'/2) = 1;
        keep = zeros(length(Floe),1);
        keep(in) = 1;
        keep(A<min_floe_size)=1;
        keep(1:Nb) = ones(Nb,1);
        
        
        %%%%%just to only fracture bonds
        keep = ones(length(Floe),1);%(Nbond+1:length(Floe)) = ones(length(Floe)-Nbond,1);
        
        
        keep = logical(keep);
        for ii = 1:length(Floe)
            FracFloes(ii).floenew = [];
        end
%         if sum(keep(1:Nbond))<length(keep(1:Nbond))
%             live = logical(cat(1,Floe.alive));
%             live(Nb+1:Nb+Nbond)=logical(keep(Nb+1:Nb+Nbond));
%             keep(Nb+1:Nb+Nbond)=1;
%             for ii = 1:length(live)
%                 if ~live(ii)
%                     Floe(ii).alive=0;
%                     Fnum = Floe(ii).bonds.FloeNum;
%                     Nums = cat(1,Floe.num);
%                     Numbers = 1:length(Floe);
%                     for jj = 1:length(Fnum)
%                         kk = Numbers(Nums==Fnum(jj));
%                         BondNum = Floe(kk).bonds.BondNum;
%                         clear R
%                         for iii = 1:length(Floe(kk).bonds.poly)
%                             R(iii) = polyshape(Floe(kk).bonds.poly{iii});%regions(polyshape(floenew.bonds.poly));
%                         end
%                         polynew = polyshape(Floe(kk).c0');
%                         poly = union(polynew,R(BondNum==Floe(ii).num)); %Fix that this is not ii but needs to be FloeNumber field
%                         if poly.NumRegions > 1
%                             [k,dist] = dsearchn(R(BondNum==Floe(ii).num).Vertices,polynew.Vertices);
%                             for iii = 1:length(dist)
%                                 if dist(iii) < 1
%                                     polynew.Vertices(iii,:) = R(BondNum==Floe(ii).num).Vertices(k(iii),:);
%                                 end
%                             end
%                             poly = union(polynew,R(BondNum==Floe(ii).num)); %Fix that this is not ii but needs to be FloeNumber field
%                             if poly.NumRegions > 1
%                                 xx = 1; xx(1) =[1 2];
%                             end
%                         end
%                         if poly.NumHoles > 0
%                             poly = rmholes(poly);
%                         end
%                         mass = area(R(BondNum==Floe(ii).num))*rho_ice*Floe(ii).h;
%                         if ~isempty(R)
%                             if ~(length(R)==length(BondNum))
%                                 xx = 1; xx(1) =[1 2];
%                             end
%                         end
%                         R(BondNum==Floe(ii).num) = []; BondNum(BondNum==Floe(ii).num) = [];
%                         Floe(kk).bonds.BondNum = BondNum;
%                         if ~isempty(R)
%                             Verts = {};
%                             for iii = 1:length(R)
%                                 Verts{iii} = R(iii).Vertices;
%                             end
%                             Floe(kk).bonds.poly = Verts;
% %                             pbnd = union(R); Floe(kk).bonds.poly = pbnd.Vertices;
%                             if ~(length(R)==length(BondNum))
%                                 xx = 1; xx(1) =[1 2];
%                             end
%                         else
%                             Floe(kk).bonds.poly = [];
%                         end
%                         [Xi_new,Yi_new] = centroid(poly);
%                         pold = polyshape(Floe(kk).c0'); [Xi_old,Yi_old] = centroid(pold);
%                         dx = Xi_new-Xi_old; dy = Yi_new-Yi_old;
%                         Floe(kk).Xi = Floe(kk).Xi+dx; Floe(kk).Yi = Floe(kk).Yi+dy;
%                         Floe(kk).area = area(poly);
%                         Floe(kk).mass = Floe(kk).mass+mass;
%                         Floe(kk).h = Floe(kk).mass/(rho_ice*Floe(kk).area);
%                         Floe(kk).c0 = poly.Vertices';
%                         A_rot=[cos(Floe(kk).alpha_i) -sin(Floe(kk).alpha_i); sin(Floe(kk).alpha_i) cos(Floe(kk).alpha_i)]; %rotation matrix
%                         Floe(kk).c_alpha=A_rot*Floe(kk).c0; %rotate floe contour
%                         Floe(kk).inertia_moment = PolygonMoments(Floe(kk).c0',Floe(kk).h);
%                         Floe(kk).angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
%                     end
%                 end
%             end
%             Nbond = sum(live(Nb+1:Nb+Nbond));
%         end
        parfor ii = 1:length(keep)
            if ~keep(ii)
                FracFloes(ii).floenew=fracture_floe(Floe(ii),3,Floe,1);
            end
        end
        fracturedFloes =[];
        for ii = 1:length(FracFloes)
            fracturedFloes = [fracturedFloes FracFloes(ii).floenew];
        end
        if isfield(fracturedFloes,'potentialInteractions')
            fracturedFloes=rmfield(fracturedFloes,{'potentialInteractions'});
        end
        if ~isempty(fracturedFloes)
%            save('FloeOld.mat','Floe')
 %           xx = 1; xx(1) =[1 2];
            Floe=[Floe(keep) fracturedFloes];
        end
    live = cat(1,Floe.alive);
    Floe(live == 0) = [];
end

