% -------------------------------------------------------------------------
% Whitened AntiPSIICOS projection
% -------------------------------------------------------------------------
% _________________________________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

function [UcorrRnk, Wps] = ProjectionAwayFromCorrelationWhitened(Gain_virt_sens, ...
    fixed, bLoad, PwrRnk, CTF)

    % Span vectors for the correlation subspace
    if fixed == 0
        Nsrc = size(Gain_virt_sens,2)/2;
        Swp = zeros(4);
        Swp(2,3) = 1;
        Swp(3,2) = 1;
        Swp(1,1) = 1; 
        Swp(4,4) = 1;
        RankG = size(Gain_virt_sens,1);
        C_re = zeros(RankG^2);
        

        if(bLoad)
            if CTF == 0
                load C_re.mat;
            else
                load C_re_CTF.mat;
            end
        else
            parfor i=1:Nsrc
                 range_i = i*2-1:i*2;
                 ai = Gain_virt_sens(:,range_i);
                 rng = 1:4;
                 X = zeros(RankG^2,4*(Nsrc-i), 'single');
                 for j = i+1:Nsrc
                    range_j = j*2-1:j*2;
                    aj = Gain_virt_sens(:,range_j);
                    X(:,rng) = kron(ai,aj)+kron(aj,ai)*Swp;
                    rng = rng+4;
                 end
                C_re = C_re + X*X';
                i
            end
            save C_re C_re
        end
    else
        [RankG, Nsrc] = size(Gain_virt_sens);
        C_re = zeros(RankG^2);

        if(bLoad)
            load('C_re_fixed.mat');
        else
            parfor i = 1:Nsrc
                 gi = Gain_virt_sens(:,i);
                 X = zeros(RankG^2,(Nsrc-i), 'single');
                 k = 1;
                 for j = (i+1):Nsrc
                    gj = Gain_virt_sens(:,j);
                    X(:,k) = kron(gi,gj)+kron(gj,gi);
                    k = k+1;
                 end
                C_re = C_re + X*X';
                i
            end
            save C_re_fixed C_re
        end
    end
 
    % create correlation matrix for power subspace
    if fixed == 0
        A = zeros(RankG^2,3*Nsrc);
        for i = 1:Nsrc
             range_i = i*2-1:i*2;
             gi = Gain_virt_sens(:,2*i-1);
             v = gi*gi';
             A(:,3*i-2) = v(:)/norm(v(:));
             gj = Gain_virt_sens(:,2*i);
             v = gj*gj';
             A(:,3*i-1) = v(:)/norm(v(:));
             v = gi*gj' + gj*gi';
             A(:,3*i) = v(:)/norm(v(:));
        end
        P = A*A';
    else
        A = zeros(RankG^2,Nsrc);
        for i = 1:Nsrc
             gi = Gain_virt_sens(:,i);
             v = gi*gi';
             A(:,i) = v(:)/norm(v(:));
        end
        P = A*A';
    end

    % whitening transformation matrix
    Wps = sqrtm(inv(P + 0.01*trace(P)/(RankG^2)*eye(size(P))));

    WCreW = Wps*C_re*Wps';
    [UcorrRnk s] = eigs(WCreW, PwrRnk);
   
end
