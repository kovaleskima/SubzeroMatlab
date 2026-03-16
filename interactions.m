% u=cat(1,Floe.Ui);
% v=cat(1,Floe.Vi);
% h=cat(1,Floe.h);
% alpha=cat(1,Floe.alpha_i);
% areas=cat(1,Floe.area);
% ksi_ice=cat(1,Floe.ksi_ice);
% alive = cat(1,Floe.alive);
% floenums = 1:N;
% nums = cat(1,Floe.num);
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
        
%         live = logical(ones(length(Floe),1));
%         live(1:i)=0; 
%         live(~alive) = 0; 
%         live(sqrt((x(i)-x).^2 + (y(i)-y).^2)>(rmax(i)+rmax))=0;
%         floenew = Floe(live);
%         Floe(i).potentialInteractions.floeNum=num2cell(floenums(live))';
%         Floe(i).potentialInteractions.Num=num2cell(nums(live))';
%         Floe(i).potentialInteractions.Ui=num2cell(u(live))';
%         Floe(i).potentialInteractions.Vi=num2cell(v(live))';
%         Floe(i).potentialInteractions.h=num2cell(h(live))';
%         Floe(i).potentialInteractions.area=num2cell(areas(live))';
%         Floe(i).potentialInteractions.Xi=num2cell(x(live))';
%         Floe(i).potentialInteractions.Yi=num2cell(y(live))';
%         Floe(i).potentialInteractions.ksi_ice = num2cell(ksi(live))';
%         Floe(i).potentialInteractions.alpha=num2cell(alpha(live))';
%         for j = 1:length(floenew)
%             Floe(i).potentialInteractions(j).c=[floenew(j).c_alpha(1,:)+floenew(j).Xi; floenew(j).c_alpha(2,:)+floenew(j).Yi];
%             Floe(i).potentialInteractions(j).bonds = floenew(j).bonds;
%         end
    end
end
%% 

kill = zeros(1,N0); transfer = kill;

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