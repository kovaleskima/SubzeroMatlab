function [Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, ubound, ucell, ocean, winds,c2_boundary, dt, HFo, min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING)

id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
rho_ice=920;
global Modulus r_mean L_mean 
% d = sqrt(bond_area);

Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
live = cat(1,Floe.alive);
Floe(live==0)=[];

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
x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
u=cat(1,Floe.Ui);
v=cat(1,Floe.Vi);
ksi=cat(1,Floe.ksi_ice);
rmax=cat(1,Floe.rmax);
alive=cat(1,Floe.alive);

for i=1+Nb:N  %do interactions with boundary in a separate parfor loop
    
    Floe(i).interactions=[];
    
    Floe(i).OverlapArea = 0;
    
    Floe(i).potentialInteractions=[];
    
    Floe(i).collision_force=[0 0];
    
    Floe(i).Stress=zeros(2);
    
    Floe(i).collision_torque=0;
    
    k=1;
    
    if ( alive(i) && ~isnan(x(i)) ) && COLLISION
        for j=1:N
            if j>i && alive(j) && sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j)) % if floes are potentially overlapping
                Floe(i).potentialInteractions(k).floeNum=j;
                Floe(i).potentialInteractions(k).Num=Floe(j).num;
                Floe(i).potentialInteractions(k).c=[Floe(j).c_alpha(1,:)+x(j); Floe(j).c_alpha(2,:)+y(j)];
                Floe(i).potentialInteractions(k).Ui=u(j);
                Floe(i).potentialInteractions(k).Vi=v(j);
                Floe(i).potentialInteractions(k).h=Floe(j).h;
                Floe(i).potentialInteractions(k).area=Floe(j).area;
                Floe(i).potentialInteractions(k).Xi=x(j);
                Floe(i).potentialInteractions(k).Yi=y(j);
                Floe(i).potentialInteractions(k).ksi_ice = ksi(j);
                Floe(i).potentialInteractions(k).alpha=Floe(j).alpha_i;
%                 Floe(i).potentialInteractions(k).bonds = cat(1,bonds(i).bond.Num);
                Floe(i).potentialInteractions(k).bonds = Floe(j).bonds;
                k=k+1;
            end
            
        end
        
%         A_rot=[cos(Floe(i).alpha_i) -sin(Floe(i).alpha_i); sin(Floe(i).alpha_i) cos(Floe(i).alpha_i)]; %rotation matrix
%         for ii = 1:length(bonds(i).bond)
% %             if bonds(i).bond(ii).piece
%                 bonds(i).bond(ii).h = Floe(i).h;
% %             else
% %                 bonds(i).bond(ii).h = Floe(bonds(i).bond(ii).Num).h;
% %             end
%             bonds(i).bond(ii).Ui = u(i);
%             bonds(i).bond(ii).Xi = x(i);
%             bonds(i).bond(ii).Vi = v(i);
%             bonds(i).bond(ii).Yi = y(i);
% %             vel_bond = ([u(i) v(i)]+ ksi(i)*([bonds(i).bond(ii).Xb bonds(i).bond(ii).Yb]-[x(i) y(i)]));
% %             Stress = [bonds(i).bond(ii).Ub-vel_bond(1) bonds(i).bond(ii).Ub-vel_bond(1)]/dt * rho_ice *bonds(i).bond(ii).area/d;
% %             bonds(i).bond(ii).Stress = bonds(i).bond(ii).Stress+Stress;
% %             bonds(i).bond(ii).Ub = vel_bond(1); bonds(i).bond(ii).Vb = vel_bond(2);
%             c0 = [bonds(i).bond(ii).c(:,1)'- bonds(i).bond(ii).Xb; bonds(i).bond(ii).c(:,2)'- bonds(i).bond(ii).Yb];
%             c_alpha=A_rot*c0;
%             bonds(i).bond(ii).c_alpha = [c_alpha(1,:)+bonds(i).bond(ii).Xb;c_alpha(2,:)+bonds(i).bond(ii).Yb];
%             bonds(i).bond(ii).ksi_ice = ksi(i);
%             
%         end
        
    end
end

kill = zeros(1,N0); transfer = kill;
% xx = 1; xx(1) =[1 2];
%for i=1+Nb:N  %now the interactions could be calculated in a parfor loop!
parfor i=1+Nb:N  %now the interactions could be calculated in a parfor loop!


    c1=[Floe(i).c_alpha(1,:)+x(i); Floe(i).c_alpha(2,:)+y(i)];
%     bondNums = cat(1,bonds(i).bond.Num);
    if ~isempty(Floe(i).potentialInteractions)
        
        for k=1:length(Floe(i).potentialInteractions)
            
            floeNum=Floe(i).potentialInteractions(k).floeNum;
            
            [force_j,P_j, overlap] = floe_interactions_poly_con2(Floe(i),Floe(i).potentialInteractions(k),c2_boundary,PERIODIC,Modulus,dt,r_mean, L_mean);
            
            if sum(abs(force_j(:)))~=0
                Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1) overlap'];
                Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
            elseif isinf(overlap)
                if i <= N0 && sign(overlap)>0
                    kill(i) = i;
                    transfer(i)=floeNum;
                elseif floeNum <= N0
                    kill(i) = floeNum;
                end
            end
            
%             Lia = ismember(bondNums,Floe(i).potentialInteractions(k).Num);
%             if sum(Lia)>0 && sum(abs(force_j(:)))~=0
%                 [force_j,P_j, overlap] = floe_interactions_poly_con2(bonds(i).bond(Lia),Floe(i).potentialInteractions(k),c2_boundary,PERIODIC,Modulus,dt,r_mean, r_bond);
% 
%                 if sum(abs(force_j(:)))~=0
%                     Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1) overlap'];
%                     Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
%                     P_j = [P_j(:,1)-bonds(i).bond(Lia).Xb-x(i) P_j(:,2)-bonds(i).bond(Lia).Yb-y(i)];
%                     bonds(i).bond(Lia).interactions=[floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1) overlap'];
%                 end
%             end
            
        end
        
    end
    if ~PERIODIC
        [force_b, P_j, overlap] = floe_interactions_poly_con2(Floe(i), floebound,c2_boundary,PERIODIC,Modulus,dt,r_mean,L_mean);
        in = inpolygon(x(i),y(i),c2_boundary(1,:)',c2_boundary(2,:)');
        if ~in
            Floe(i).alive = 0;
        end
        
%     if ~worked, display(['contact points issue for (' num2str(i) ', boundary)' ]); end
        if sum(abs(force_b(:)))~=0
            [mm,~] = size(P_j);
            for ii =1:mm
                if abs(P_j(ii,2)) == Ly
                    force_b(ii,1) = 0;
                end
            end
            % boundary will be recorded as floe number Inf;
            Floe(i).interactions=[Floe(i).interactions ; Inf*ones(size(force_b,1),1) force_b P_j zeros(size(force_b,1),1) overlap'];
            Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
            Floe(i).potentialInteractions(end+1).floeNum = Inf;
            Floe(i).potentialInteractions(end).c = c2_boundary;
        end
    end
    
end
for i = 1:length(kill);
    if abs(kill(i)-i)>0 && kill(i)>0;
        transfer(kill(i))=i;
    end
end

if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end

if isfield(Floe,'potentialInteractions')
    Floe=rmfield(Floe,{'potentialInteractions'});
end

%Fill the lower part of the interacton matrix (floe_i,floe_j) for floes with j<i
for i=1:N %this has to be done sequentially
      
    if ~isempty(Floe(i).interactions)
        
        a=Floe(i).interactions;
        
        indx=a(:,1);
        
        for j=1:length(indx)
            
            if indx(j)<=N && indx(j)>i
                Floe(indx(j)).interactions=[Floe(indx(j)).interactions; i -a(j,2:3) a(j,4:5) 0 a(j,7)];   % 0 is torque here that is to be calculated below
                Floe(indx(j)).OverlapArea = Floe(indx(j)).OverlapArea + a(j,7);
%                 m = size(Floe(indx(j)).potentialInteractions,2);
%                 Floe(indx(j)).potentialInteractions(m+1).floeNum=i;
%                 Floe(indx(j)).potentialInteractions(m+1).c=[Floe(i).c_alpha(1,:)+x(i); Floe(i).c_alpha(2,:)+y(i)];
%                 Floe(indx(j)).potentialInteractions(m+1).Ui=Floe(i).Ui;
%                 Floe(indx(j)).potentialInteractions(m+1).Vi=Floe(i).Vi;
%                 Floe(indx(j)).potentialInteractions(m+1).Xi=x(i);
%                 Floe(indx(j)).potentialInteractions(m+1).Yi=y(i);
%                 Floe(indx(j)).potentialInteractions(m+1).alpha=Floe(i).alpha_i;
%                 Floe(indx(j)).potentialInteractions(m+1).ksi_ice = Floe(i).ksi_ice;
            end
            
        end
    end

end

% calculate all torques from forces
if PERIODIC
    
   parfor i=N0+1:N %do this in parfor
        
        if ~isempty(Floe(i).interactions)
            
            a=Floe(i).interactions;
            r=[x(i) y(i)];
            for k=1:size(a,1)
                floe_Rforce=a(k,4:5);
                floe_force=a(k,2:3);
                floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
                Floe(i).interactions(k,6)=floe_torque(3);
            end
            
            Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1);
            Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1);
            
        end
        
    end
     %add forces and torques from ghost floes to their parents; ghost floes
    %begin with the index N0+1
    for i=1:length(parent)
        Floe(parent(i)).collision_force =Floe(parent(i)).collision_force +Floe(N0+i).collision_force;
        Floe(parent(i)).collision_torque=Floe(parent(i)).collision_torque+Floe(N0+i).collision_torque;
    end
end

keep = ones(1,N0);
%parfor i=1+Nb:N0
for i=1+Nb:N0
    
    if ~isempty(Floe(i).interactions)
        
        a=Floe(i).interactions;
        r=[x(i) y(i)];
        for k=1:size(a,1)
            floe_Rforce=a(k,4:5);
            floe_force=a(k,2:3);
            floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
            Floe(i).interactions(k,6)=floe_torque(3);
        end
        
        Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1)+Floe(i).collision_force;
        Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1)+Floe(i).collision_torque;
        
    end
    
    if PERIODIC
        
        if abs(Floe(i).Xi)>Lx %if floe got out of periodic bounds, put it on the other end
            Floe(i).Xi=Floe(i).Xi-2*Lx*sign(Floe(i).Xi);
        end
        
        if abs(Floe(i).Yi)>Ly %if floe got out of periodic bounds, put it on the other end
            Floe(i).Yi=Floe(i).Yi-2*Ly*sign(Floe(i).Yi);
        end
        
    end
    
    if Floe(i).alive
        [tmp, frac,Fx,Fy] =calc_trajectory(dt,ocean,winds,Floe(i),HFo,doInt,c2_boundary, ubound,ucell,i);
        if (isempty(tmp) || isnan(x(i)) ), kill(i)=i;xx = 1; xx(1) =[1 2]; elseif frac == 1, keep(i) = 0; else; Floe(i)=tmp; Floe(i).Fx = Fx; Floe(i).Fy = Fy; end
    end
end


% nums = cat(1,Floe.num);
% bond0 = bonds;
% parfor ii = 1:length(bonds)
%     for jj = 1:length(bonds(ii).bond)
%         if ~isempty(bonds(ii).bond(jj).interactions)
%             Lia = ismember(nums,bonds(ii).bond(jj).Num);
%             bnd_tmp = bond0(Lia).bond;
%             bond_num = cat(1,bnd_tmp.Num);
%             Lia2 = ismember(bond_num,nums(ii));
%             a=bonds(ii).bond(jj).interactions;
%             a2 = bnd_tmp(Lia2).interactions;
%             a = [a; a2];
%             r=[bonds(ii).bond(jj).Xi bonds(ii).bond(jj).Yi];
%             Stress =1/(2*bonds(ii).bond(jj).area*bonds(ii).bond(jj).h)*([sum((a(:,4)-r(1)).*a(:,2)) sum((a(:,5)-r(2)).*a(:,2)); sum((a(:,4)-r(1)).*a(:,3)) sum((a(:,5)-r(2)).*a(:,3))]...
%                 +[sum(a(:,2).*(a(:,4)-r(1))) sum(a(:,3).*(a(:,4)-r(1))); sum(a(:,2).*(a(:,5)-r(2))) sum(a(:,3).*(a(:,5)-r(2)))]);
%             bonds(ii).bond(jj).Stress = bonds(ii).bond(jj).Stress+Stress;
%         end
%     end
% end

% Areas = cat(1,Floe.area);
% Nbond = length(Areas(Areas<bond_area+1));
% for ii = 1:Nbond
%     if Floe(ii).bonds.interactions >10
%         Floe(ii).alive = 0;
%         Fnum = Floe(ii).bonds.FloeNum;
%         Nums = cat(1,Floe.num);
%         Numbers = [1:length(Floe)]';
%         for jj = 1:length(Fnum)
%             kk = Numbers(Nums==Fnum(jj));
%             BondNum = Floe(kk).bonds.BondNum;
%             clear R
%             for iii = 1:length(Floe(kk).bonds.poly)
%                 R(iii) = polyshape(Floe(kk).bonds.poly{iii});%regions(polyshape(floenew.bonds.poly));
%             end
%             polynew = polyshape(Floe(kk).c0');
%             poly = union(polynew,R(BondNum==Floe(ii).num)); %Fix that this is not ii but needs to be FloeNumber field
%             if poly.NumRegions > 1
%                 [k,dist] = dsearchn(R(BondNum==Floe(ii).num).Vertices,polynew.Vertices);
%                 for iii = 1:length(dist)
%                     if dist(iii) < 1
%                         polynew.Vertices(iii,:) = R(BondNum==Floe(ii).num).Vertices(k(iii),:);
%                     end
%                 end
%                 poly = union(polynew,R(BondNum==Floe(ii).num)); %Fix that this is not ii but needs to be FloeNumber field
%                 if poly.NumRegions > 1
%                     xx = 1; xx(1) =[1 2];
%                 end
%             end
%             if poly.NumHoles > 0
%                 poly = rmholes(poly);
%             end
%             mass = area(R(BondNum==Floe(ii).num))*rho_ice*Floe(ii).h;
%             if ~isempty(R)
%                 if ~(length(R)==length(BondNum))
%                     xx = 1; xx(1) =[1 2];
%                 end
%             end
%             R(BondNum==Floe(ii).num) = []; BondNum(BondNum==Floe(ii).num) = [];
%             Floe(kk).bonds.BondNum = BondNum;
%             if ~isempty(R)
%                 Verts = {};
%                 for iii = 1:length(R)
%                     Verts{iii} = R(iii).Vertices;
%                 end
%                 Floe(kk).bonds.poly = Verts;
% %                 pbnd = union(R); Floe(kk).bonds.poly = pbnd.Vertices;
%                 if ~(length(R)==length(BondNum))
%                     xx = 1; xx(1) =[1 2];
%                 end
%             else
%                 Floe(kk).bonds.poly = [];
%             end
%             [Xi_new,Yi_new] = centroid(poly);
%             pold = polyshape(Floe(kk).c0'); [Xi_old,Yi_old] = centroid(pold);
%             dx = Xi_new-Xi_old; dy = Yi_new-Yi_old;
%             Floe(kk).Xi = Floe(kk).Xi+dx; Floe(kk).Yi = Floe(kk).Yi+dy;
%             Floe(kk).area = area(poly);
%             Floe(kk).mass = Floe(kk).mass+mass;
%             Floe(kk).h = Floe(kk).mass/(rho_ice*Floe(kk).area);
%             Floe(kk).c0 = poly.Vertices';
%             A_rot=[cos(Floe(kk).alpha_i) -sin(Floe(kk).alpha_i); sin(Floe(kk).alpha_i) cos(Floe(kk).alpha_i)]; %rotation matrix
%             Floe(kk).c_alpha=A_rot*Floe(kk).c0; %rotate floe contour
%             Floe(kk).inertia_moment = PolygonMoments(Floe(kk).c0',Floe(kk).h);
%             Floe(kk).angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
%         end
%     else
%         Floe(ii).alive = 1;
%     end
% end

%% 

floenew = [];
Ridge = zeros(1,length(Floe));
if RIDGING && doInt.flag
    %Create a function to control probability that ridging will occur
    h = cat(1,Floe.h);
    keepR=rand(length(Floe),1)<0.05;
    for ii=1+Nb:N0
        
        if Floe(ii).alive && ~isempty(Floe(ii).interactions)
            a = Floe(ii).interactions;
            c1 = Floe(ii).c_alpha+[Floe(ii).Xi; Floe(ii).Yi];
            abound = zeros(1+Nb,1);
            if ~isempty(a)
                if ~isempty(InterX(c1,c2_boundary))
                    abound(1+Nb) = 1;
                end
                a(isinf(a(:,1)),:)=[];
            end

            if  ~keepR(ii) && h(ii)<5  && ~isempty(a)
                clear overlap;
                for jj = 1:size(a,1)
                    if a(jj,1) < length(Floe)+1
                        overlap(jj) = a(jj,7)/min([Floe(ii).area Floe(a(jj,1)).area]);
                    else
                        overlap(jj) = a(jj,7)/Floe(ii).area;
                    end
                end
                overlap(overlap<1e-6) = nan; overlap(overlap>0.95) = nan;
                overlappingFloes = a(~isnan(overlap),1);
                overlappingFloes = unique(overlappingFloes);
                abound(overlappingFloes<Nb+1) = 1;
                for jj = length(overlappingFloes):-1:1
                    if Ridge(overlappingFloes(jj))
                        overlappingFloes(jj)=[];
                    end
                end
                for jj = 1:length(overlappingFloes)
                    if Floe(overlappingFloes(jj)).h < 5
                        [Floe1, Floe2] = ridging(Floe(ii),Floe(overlappingFloes(jj)),c2_boundary_poly,PERIODIC,min_floe_size);
                        if length(Floe1) > 1
                            Floe(ii) = Floe1(1);
                            Ridge(ii) = 1;
                            floenew = [floenew Floe1(2:end)];
                        else
                            Floe(ii) = Floe1;
                            if Floe1.alive == 0
                                kill(ii) = ii;
                            end
                        end
                        if length(Floe2) > 1
                            Floe(overlappingFloes(jj)) = Floe2(1);
                            Ridge(overlappingFloes(jj)) = 1;
                            floenew = [floenew Floe2(2:end)];
                        else
                            Floe(overlappingFloes(jj)) = Floe2;
                            if Floe2.alive == 0 && overlappingFloes(jj) <= N0
                                kill(overlappingFloes(jj)) = overlappingFloes(jj);
                            end
                        end
                    end
                end

            end
            if sum(abound)>0 && h(ii)<1.25 && Floe(ii).area > min_floe_size
                for jj = 1:length(abound)
                    if abound(jj) == 1 && jj == Nb+1
                        [Floe1, ~] = ridging(Floe(ii),floebound,c2_boundary_poly,PERIODIC,min_floe_size);
                    elseif abound(jj)  == 1
                        poly = polyshape([Floe(abound(jj)).c_alpha(1,:)+Floe(abound(jj)).Xi; Floeabound(jj).c_alpha(2,:)+Floe(abound(jj)).Yi]');
                        [Floe1, ~] = ridging(Floe(ii),Floe(abound(jj)),poly,PERIODIC,min_floe_size);
                    end
                    if length(Floe1) > 1
                        Floe(ii) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(ii) = Floe1;
                        if Floe1.alive == 0
                            kill(ii) = ii;
                        end
                    end
                end
            end
        end
    end
end


Raft = zeros(1,length(Floe));
if RAFTING && doInt.flag
    %Create a function to control probability that ridging will occur
    h = cat(1,Floe.h);
    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
    keepR=rand(length(Floe),1)>0.5*overlapArea;
    for ii=1+Nb:N0
        
        if Floe(ii).alive && ~isempty(Floe(ii).interactions)
            a = Floe(ii).interactions;
            c1 = Floe(ii).c_alpha+[Floe(ii).Xi; Floe(ii).Yi];
            abound = zeros(1+Nb,1);
            if ~isempty(a)
                if ~isempty(InterX(c1,c2_boundary))
                    abound(1+Nb) = 1;
                end
                a(isinf(a(:,1)),:)=[];
            end

            if  ~keepR(ii) && h(ii)<0.25  && ~isempty(a)
                clear overlap;
                for jj = 1:size(a,1)
                    if a(jj,1) < length(Floe)+1
                        overlap(jj) = a(jj,7)/min([Floe(ii).area Floe(a(jj,1)).area]);
                    else
                        overlap(jj) = a(jj,7)/Floe(ii).area;
                    end
                end
                overlap(overlap<1e-6) = nan; overlap(overlap>0.95) = nan;
                overlappingFloes = a(~isnan(overlap),1);
                overlappingFloes = unique(overlappingFloes);
                abound(overlappingFloes<Nb+1) = 1;
                for jj = length(overlappingFloes):-1:1
                    if Raft(overlappingFloes(jj))
                        overlappingFloes(jj)=[];
                    end
                end
                for jj = 1:length(overlappingFloes)
                    if Floe(overlappingFloes(jj)).h < 0.25
                        [Floe1, Floe2] = rafting(Floe(ii),Floe(overlappingFloes(jj)),c2_boundary_poly,PERIODIC,min_floe_size);
                        
                        if length(Floe1) > 1
                            Floe(ii) = Floe1(1);
                            Raft(ii) = 1;
                            floenew = [floenew Floe1(2:end)];
                        else
                            Floe(ii) = Floe1;
                            if Floe1.alive == 0
                                kill(ii) = ii;
                            end
                        end
                        if length(Floe2) > 1
                            Floe(overlappingFloes(jj)) = Floe2(1);
                            Raft(overlappingFloes(jj)) = 1;
                            floenew = [floenew Floe2(2:end)];
                        else
                            Floe(overlappingFloes(jj)) = Floe2;
                            if Floe2.alive == 0 && overlappingFloes(jj) <= N0
                                kill(overlappingFloes(jj)) = overlappingFloes(jj);
                            end
                        end
                    end
                end

            end
            if sum(abound)>0 && h(ii)<0.25 && Floe(ii).area > min_floe_size
                for jj = 1:length(abound)
                    if abound(jj) == 1 && jj == Nb+1
                        [Floe1, ~] = rafting(Floe(ii),floebound,c2_boundary_poly,PERIODIC,min_floe_size);
                    elseif abound(jj)  == 1
                        poly = polyshape([Floe(abound(jj)).c_alpha(1,:)+Floe(abound(jj)).Xi; Floeabound(jj).c_alpha(2,:)+Floe(abound(jj)).Yi]');
                        [Floe1, ~] = rafting(Floe(ii),Floe(abound(jj)),poly,PERIODIC,min_floe_size);
                    end
                    if length(Floe1) > 1
                        Floe(ii) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(ii) = Floe1;
                        if Floe1.alive == 0
                            kill(ii) = ii;
                        end
                    end
                end
            end
        end
    end
end

Floe=Floe(1:N0); % ditch the ghost floes.

killO = kill; transferO = transfer; FloeO = Floe;

if ~isempty(kill(kill>0)) 
    kill(kill>N0)=0;
    transfer(transfer>N0) = 0;
    transfer = transfer(kill>0);
    kill = kill(kill>0);
    for ii = 1:length(kill)
        if Floe(kill(ii)).alive>0 
            if transfer(ii)>0 && (Floe(kill(ii)).area>2e4 || Floe(kill(ii)).area>2e4)
                Floe1 = Floe(kill(ii));
                Floe2 = Floe(transfer(ii));
                Floe1.poly = polyshape(Floe1.c_alpha'+[Floe1.Xi Floe1.Yi]);
                Floe2.poly = polyshape(Floe2.c_alpha'+[Floe2.Xi Floe2.Yi]);
                floes  = Floe(transfer(ii));%FuseFloes(Floe1,Floe2);
                if isfield(floes,'poly')
                    floes=rmfield(floes,{'poly'});
                end
                if length(floes)>1
                    Floe(transfer(ii)) = floes(1);
                    for jj = 2:length(floes)
                        floenew = [floenew floes(2:end)];
                    end
                elseif isempty(floes)
                    save('FloeTransfer.mat','Floe','kill','transfer','Floe1','Floe2')
                else
                    Floe(transfer(ii)) = floes(1);
                end
            end
            dissolvedNEW = dissolvedNEW+calc_dissolved_mass(Floe(kill(ii)),Nx,Ny,c2_boundary_poly);
            Floe(kill(ii)).alive = 0;
        end
    end
end
Floe = [Floe floenew];
live = cat(1,Floe.alive);

Floe(live==0)=[]; %remove any floes that got dissolved so they do not take up space

if isfield(Floe,'potentialInteractions')
    Floe=rmfield(Floe,{'potentialInteractions'});
end

warning('on',id)
warning('on',id3)
end
