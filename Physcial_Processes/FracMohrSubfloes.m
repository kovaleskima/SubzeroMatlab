function [Floe] = FracMohrSubfloes(Floe,Nb,min_floe_size)
%Use Mohr's cone to determine which floes are fractured
%   If Principal stresses are outside the cone then the floes are fractured
Floe0 = Floe;
A = cat(1,Floe.area);
q = 5.2; SigC = 250e3;
Sig1 = (1/q+1)*SigC/(1/q-q);
Sig2 = q*Sig1+SigC;
Sig11 = 4e4;
Sig22 = q*Sig11+SigC;
MohrX = [Sig1; Sig11; Sig22];
MohrY = [Sig2; Sig22; Sig11];
Mohr = polyshape(-MohrX,-MohrY);
% xx = 1; xx(1) =[1 2];
keep = logical(ones(length(Floe),1));
for ii = 1:length(Floe)
    FracFloes(ii).floenew = [];
end
for ii = 1:length(Floe)
%     Stress = Floe(ii).Stress;
    Stress = eig(Floe(ii).Stress);
    Princ1 = max(Stress); Princ2 = min(Stress);
    for jj = 1:length(Floe(ii).bonds.Num)
        bonds = Floe(ii).bonds.Num{jj};
        theta = Floe(ii).bonds.Theta{jj};
%         P1 = zeros(1,length(theta)); P2 = zeros(1,length(theta)); 
        in = logical(zeros(1,length(theta)));
        for kk = 1:length(theta)
%             Q = [cos(Floe(ii).alpha_i + theta(kk)) sin(Floe(ii).alpha_i + theta(kk)); -sin(Floe(ii).alpha_i + theta(kk)) cos(Floe(ii).alpha_i + theta(kk))];
%             PrincSig = Q * Stress * Q.';
%             P1 = [PrincSig(1,1); P2(kk) = PrincSig(2,2);
            Mohr2 = rotate(Mohr,abs(theta(kk)*180/pi),[0 0]);
            [inpoly,~] = inpolygon(Princ1,Princ2,Mohr2.Vertices(:,1), Mohr2.Vertices(:,2));
            in(kk) = inpoly;
        end
%         [in,~] = inpolygon(P1,P2,Mohr.Vertices(:,1), Mohr.Vertices(:,2));
        Floe(ii).bonds.Num{jj} = bonds(in);
        Floe(ii).bonds.Theta{jj} = theta(in);
    end
    t = []; s = []; count = 1;
    for jj = 1:length(Floe(ii).bonds.Num)
        bonds = Floe(ii).bonds.Num{jj};
        if isempty(bonds)
%             xx = 1; xx(1) =[1 2];
            s = [s jj];
            t = [t length(Floe(ii).bonds.Xs)+count];
            count = count+1;
        end
        s = [s jj*ones(1,length(bonds))];
        t = [t bonds];
    end
    G = digraph(s,t);
    floe = Floe(ii);
    polyOrig = polyshape(floe.c_alpha');
    clear poly
    [bins,binsize] = conncomp(G,'Type','weak');
    b=floe.bonds.Vert;
    A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)];
    for jj = 1:length(b); poly(jj) = polyshape((A_rot*b{jj}')'); poly(jj) = intersect(poly(jj),polyOrig); end
    if length(binsize)>1
        for jj = length(bins):-1:length(Floe(ii).bonds.Num)+1
            binsize(bins(jj)) = binsize(bins(jj))-1;
            if binsize(bins(jj)) == 0
                bins(bins > jj) = bins(bins > jj)-1;

            end
            bins(jj) = [];
        end
        H = histogram(bins); binsize = H.Values;
        Nbins = length(binsize);
        clear polyFrac
        for jj = Nbins:-1:1
            polyNEW = union(poly(bins==jj));
            polyFrac(jj) = polyNEW;
        end
        NholeMax = 1;
        while NholeMax > 0
            clear Nholes
            clear polyFrac; clear polyNEW
            Nbins = length(binsize);
            for jj = Nbins:-1:1
                polyNEW = union(poly(bins==jj));
                polyFrac(jj) = polyNEW;
                Nholes(jj) = polyNEW.NumHoles;
                if polyNEW.NumHoles > 0
                    P = [Floe(ii).bonds.Xs Floe(ii).bonds.Ys];
                    bink = 1:length(bins);
                    bink = bink(bins==jj);
                    P(bins == jj,:) = 1e7*ones(size(P(bins == jj,:)));
                    PQ = [Floe(ii).bonds.Xs(bins == jj) Floe(ii).bonds.Ys(bins == jj)];
                    [k,~] = dsearchn(P,PQ);
                    for kk = 1:length(k)
                        bins(bink(kk)) = bins(k(kk));
                    end
                    bins(bins > jj) = bins(bins > jj)-1;
                    H = histogram(bins); binsize = H.Values;
                end
            end
            NholeMax = max(Nholes);
            if length(binsize) ==1
                NholeMax = 0;
            end
        end
    end
    Nbins = length(binsize);
    polyNEW = [];
    for jj = Nbins:-1:1
        polyFrac = union(poly(bins==jj));
        R = regions(polyFrac);
        polyNEW = [polyNEW; R];
    end
    if abs(length(bins)-sum(binsize))>0
        xx = 1; xx(1) =[1 2];
    elseif sum(area(polyNEW))/floe.area-1 > 0.1
        xx = 1;
        xx(1) =[1 2];
    end
    if length(polyNEW)>1 && ii > Nb && Floe(ii).area > min_floe_size
        keep(ii) = 0;
            FracNEW=fracture_floe_subfloes(Floe(ii),30,Floe0,polyNEW);
        FracFloes(ii).floenew=FracNEW;
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
    Floe=[Floe(keep) fracturedFloes];
end
if sum(cat(1,Floe.area))/sum(cat(1,Floe0.area))-1 > 0.005
    xx = 1;
    xx(1) =[1 2];
end
Lx = 1e5;
if sum(cat(1,Floe.area))/(4*Lx^2)-1 > 0.005
    xx = 1; xx(1) =[1 2];
end
end

