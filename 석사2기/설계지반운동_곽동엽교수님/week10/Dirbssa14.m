function p = Dirbssa14(M,Vs30,Rjb,ftype,region,dz1,ab)
    Co = readmatrix('BSSA14_Coefficients_071314_Revisedf4_071514.csv');
    
% Example Data
%     % M
%     M = 6.5;
%     
%     % region 
%     % 1 = globalCATWNZ 2 = ChinaTurkey 3 = Dc3ItalyJapan
%     region = 1;
%     
%     % dz1
%     dz1 = 0;
%     
%     % Rjb (km)
%     Rjb = 10;
%     
%     %Vs30 (m/s)
%     Vs30 = 760;
% 
%     %ftype
%     ftype = 4;
    
  
% Coefficients
    period = Co(:,1);
    
    % Magnitude Scaling
    e0  = Co(:,2);
    e1	= Co(:,3);
    e2  = Co(:,4);
    e3  = Co(:,5);
    e4  = Co(:,6);
    e5  = Co(:,7);
    e6  = Co(:,8);
    Mh  = Co(:,9);
    
    % Distance Scaling
    c1  = Co(:,10);
    c2  = Co(:,11);
    c3  = Co(:,12);
    Mref= Co(:,13);
    Rref= Co(:,14);
    h	= Co(:,15);
    
    % Anelastic attenuation
    Dc3globalCATWNZ = Co(:,16);
    Dc3ChinaTurkey	= Co(:,17);
    Dc3ItalyJapan	= Co(:,18);
    
    % Linear site term
    clin= Co(:,19);
    Vc	= Co(:,20);
    Vref= Co(:,21);
    
    % Ninlinear site term
    f1	= Co(:,22);
    f3	= Co(:,23);
    f4	= Co(:,24);
    f5	= Co(:,25);
    
    % Basin depth
    f6	= Co(:,26);
    f7	= Co(:,27);
    
    % Aleatory Uncertainty
    R1	= Co(:,28);
    R2	= Co(:,29);
    DfR	= Co(:,30);
    DfV	= Co(:,31);
    V1	= Co(:,32);
    V2	= Co(:,33);
    phi1= Co(:,34);
    phi2= Co(:,35);
    tau1= Co(:,36);
    tau2= Co(:,37);
    
% ftype
% 1 = unspecified / 2 = strike slip / 3 = normal slip / 4 = reverse slip
    if ftype == 1
        US = 1;
        SS = 0;
        NS = 0;
        RS = 0;
    elseif ftype == 2
        US = 0;
        SS = 1;
        NS = 0;
        RS = 0;
    elseif ftype == 3
        US = 0;
        SS = 0;
        NS = 1;
        RS = 0;
    else 
        US = 0;
        SS = 0;
        NS = 0;
        RS = 1;
    end
  
% Calculate the PGAr
    if   region == 1
        Dc3 = Dc3globalCATWNZ;
    elseif  region == 2
        Dc3 = Dc3ChinaTurkey;
    else 
        Dc3 = Dc3ItalyJapan;
    end
  
    i=2;
    
    % source
    if M <= Mh(i)
        Fm1 = e0(i) .* US + e1(i) .* SS + e2(i) .* NS + e3(i) .* RS + e4(i) .* (M - Mh(i)) + e5(i) .* ((M - Mh(i)).^2);
    else 
        Fm1 = e0(i) .* US + e1(i) .* SS + e2(i) .* NS + e3(i) .* RS + e6(i) .* (M - Mh(i));
    end

    %path
    R = sqrt(Rjb^2 + (h(i).^2));
    
    Fd1 = (c1(i) + c2(i) .* (M - Mref(i))) .* log(R ./ Rref(i)) + (c3(i) + Dc3(i)) .* (R - Rref(i));
    
    pgar = exp(Fm1 + Fd1);

% Calcuation with PGAr
% Magnitude
    FM = zeros(107,1);
    for i = 1 : 1 : 107
        if M <= Mh(i)
            FM(i) = e0(i) .* US + e1(i) .* SS + e2(i) .* NS + e3(i) .* RS + e4(i) .* (M - Mh(i)) + e5(i) .* ((M - Mh(i)).^2);
        else 
            FM(i) = e0(i) .* US + e1(i) .* SS + e2(i) .* NS + e3(i) .* RS + e6(i) .* (M - Mh(i));
        end
    end

    
% Path
    R = sqrt(Rjb .^ 2 + (h .^ 2)) ;
    Fp = (c1 + c2 .* (M - Mref)).* log(R ./ Rref) + (c3 + Dc3) .* (R - Rref);

% Site
% basin depth on fround motion amplitude
% We realize that in many applications
% z1 may be unknown; in such cases we recommend using the default value of ??z1 = 0.0, which
% turns off this adjustment factor (i.e., F??z1 = 0).
    Fbd1 = zeros(107,1);
    f76 = f7./f6;
    for i = 1 : 1 : 107
        if(dz1 == 0 && period(i) < 0.65)
            Fbd1(i) = 0;
        elseif(period(i) >= 0.65 && dz1 <= f76(i))
            Fbd1(i) = f6(i) * dz1;
        elseif(period(i) >= 0.65 && dz1 > f76(i))
            Fbd1(i) = f7(i);
        else 
            Fbd1(i) =  0;
        end
    end


% linear component of site amplification 
    Flin = zeros(107,1);
    for i = 1 : 1 : 107
        if Vs30 <= Vc(i)
            Flin(i) = clin(i) * log(Vs30 / Vref(i));  
        else 
            Flin(i) = clin(i) * log(Vc(i) / Vref(i));
        end
    end


% nonlinear component of site amplification 
    Fnl = zeros(107,1);
    f2 = zeros(107,1);
    for i = 1 : 1 : 107
        if(Vs30 == Vref(i))
            Fnl(i) = 0;
        else 
            f2(i) = f4(i) * (exp(f5(i) * (min(Vs30, Vref(i)) - 360)) - exp(f5(i) * (Vref(i) - 360)));
            Fnl(i) = f1(i) + f2(i) * log((pgar + f3(i)) / f3(i));
        end
    end

    
% ALEATORY-UNCERTAINTY FUNCTION
    % between event variability Tau (??)
    if M <= 4.5
        tauM = tau1;
    elseif(M > 4.5 && M<5.5)
        tauM = tau1+(tau2-tau1)*(M-4.5);
    else 
        tauM = tau2;
    end

    
    % within-event variability (??)
    if M <= 4.5
        phiM = phi1;
    elseif(M > 4.5 && M<5.5)
        phiM= phi1+(phi2-phi1).*(M-4.5);
    else 
        phiM = phi2;
    end
    
    phiMR = zeros(107,1);
    for i = 1 : 1 : 107
        if Rjb <= R1(i)
            phiMR(i) =  phiM(i);
        elseif(Rjb > R1(i) && Rjb <= R2(i))
            phiMR(i) = phiM(i) + DfR(i)*(log(Rjb/R1(i))/log(R2(i)/R1(i)));
        else 
            phiMR(i) = phiM + DfR(i);
        end
    end

    
    phiMRV = zeros(107,1);
    for i = 1 : 1 : 107
        if Vs30 >= V2(i)
            phiMRV(i) = phiMR(i);
        elseif(Vs30 >= V1(i) && Vs30 < V2(i))
            phiMRV(i) = phiMR(i) - DfV(i)*(log(V2(i)/Vs30)/log(V2(i)/V1(i)));
        else 
            phiMRV(i) = phiMR(i) - DfV(i);
        end
    end

    
    % total standard deviation ??
    sigMRV = sqrt((phiMRV).^2 + tauM.^2);
    

%% Evaluate the Isochrone Directivity Parameter(IDP)
h = 10;

a1x = -24;
a1y = 30;
a2x = -20;
a2y = 40;
as = sqrt((a1x)^2 + (a1y)^2);
aRrup = sqrt((a1x-a2x)^2 + (a1y-a2y)^2);
aRhyp = sqrt((a2x)^2 + (a2y)^2 + 10^2);
aRt = 0.9;
aRu = 0.3;
aD = sqrt((a1x)^2 + (a1y)^2 + 10^2);

b1x = -80/41;
b1y = -5/4 * b1x;
b2x = 20;
b2y = 20;
bs = sqrt((b1x)^2 + (b1y)^2);
bRrup = sqrt((b1x-b2x)^2 + (b1y-b2y)^2);
bRhyp = sqrt((b2x)^2 + (b2y)^2 + 10^2);
bRt = 0.1;
bRu = 1;
bD = 10;
%% Calculate AIDP
cbar = (1/0.8 - (aRhyp - aRrup)/aD )^(-1);
C = (min(cbar, 2.45)-0.8) / (2.45 - 0.8);
S = log(min(75,max(as,h)));
ARri = max(sqrt(aRu^2 + aRt^2),0.2);
AIDP = C * S * ARri;
%% Calculate BIDP
cbar = (1/0.8 - (bRhyp - bRrup)/bD )^(-1);
C = (min(cbar, 2.45)-0.8) / (2.45 - 0.8);
S = log(min(75,max(bs,h)));
BRri = max(sqrt(bRu^2 + bRt^2),0.2);
BIDP = C * S * BRri;
%% IDPbar
a = 1.192885;
b = 0.909248;
c = 0.155513;
z = 0;
aIDPbar = a + b/(cosh(c*max(aRrup - z,0)));
bIDPbar = a + b/(cosh(c*max(bRrup - z,0)));
%% qMT
C4 = -1.1736;
C5 = 0.2971;
g = 0.6132;
qMT = -(log10(period) - (C4 + C5*M)).^2/(2*g^2) ;
%% bMT
C2 = 0.0823;
C3 = 0.1665;
C1 = 5.7;
bMT = (C2 + C3*max(M-C1,0))*exp(qMT);
%% frRRR Distance Taper
afrRRR = max(0,1-max(0,aRrup - 40)/30);
bfrRRR = max(0,1-max(0,bRrup - 40)/30);
%% fD
afD = afrRRR*bMT*(AIDP-aIDPbar);
bfD = bfrRRR*bMT*(BIDP-bIDPbar);

if ab == 1
    Fd = afD;
else
    Fd = bfD;
end
%% Results of the model
    
    lnY = Fp + FM + Fbd1 + Fnl + Flin + Fd;
    %+ sigMRV
    expY = exp(lnY);

    period(1,:) = [];
    lnY(1,:) = [];
    expY(1,:) = [];
    
    period(1,:) = [];
    lnY(1,:) = [];
    expY(1,:) = [];
    
p = [period lnY expY];
end

