%%%%%%%%%%%   Compare how theory and data react to minimization
%%%% We will try to make the theoretical prob. distribution just as dense
%%%% and with the same cuttoff to make a fair comparison

clear all


%%%% LOAD DATA %%%
%erg=load('C:\Users\customer\Desktop\AL-FOLDER\TOMOGRAPHY\scramble_case\Feb08_trial2.csv');
erg=load('C:\Users\customer\Desktop\AL-FOLDER\TOMOGRAPHY\Data2007\May_07\loopy_48.csv');
erg=erg(1:3:300,1:10)';
%erg=erg(1:1:144,1:10)';
si= size(erg);

%%%%  RENORMALIZE  to a probability distribution
 ergs = sum(erg) - erg(1,:);
 for i1 = 2:si(1)
     erg(i1,:) = erg(i1,:)./ergs;
 end




%%%% CUTOFF PARAMETERS %%%
alpha_cut=floor(erg(1,si(2)))+1;  %% alpha_cut is the cuttoff value in phase space, 

NrPh=50;                 %% nr. of photons up to which we expand the POVMs %%%

DiNr=20;        %% cutoff Nr or photons for display                 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  THE THEORETICAL POVMs %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LoopyLoss=loss(theory(9,NrPh-1),0.5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      THE THEORETICAL Probability Distributions    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density_f=floor(si(2)/alpha_cut)+1;



for nalpha = 1:Density_f*alpha_cut
    PD_LoopyLoss(1,nalpha) = nalpha./Density_f;
    v(1:NrPh) = exp(-PD_LoopyLoss(1,nalpha))*PD_LoopyLoss(1,nalpha).^(0:(NrPh-1))./factorial(0:(NrPh-1));
    for bi = 0:7
        PD_LoopyLoss(2+bi,nalpha) = LoopyLoss(:,bi+1)'*v';
    end
    PD_LoopyLoss(2+8,nalpha) = 1 - sum(PD_LoopyLoss(2:9,nalpha));
end
%
Terg=PD_LoopyLoss(:,1:1:si(2));  %saturates around alpha^2=120 
%note that for alpha =120, P(alpha^2=120, n=7, loss=50) = 0.005)

Tsi = size(Terg); %Determines size of the data set



figure(1)
for i=2:10
 plot(PD_LoopyLoss(1,:),PD_LoopyLoss(i,:))
 hold on
end
hold on
for i=2:10
 plot(erg(1,:),erg(i,:),':r')
 hold on
end
title 'Loopy and loss (theory & Data)'


 
% erg(1,:)=erg(1,:)*0.95; % Apply scaling factor due to unknown calibration
                         %of power meter 
mmax = 20+floor(erg(1,si(2))*exp(1)); % Determine cut-off in Fock space
imax = si(2);                         % Number of data points

% Create binomial coefficients, factorials and gamma function
% I do this to speed up code as these operations are computationally
% very intensive in the determination of F matrix
%
for i1 = 1:2*(mmax+1)
    bino(i1,1) = 1;
    for i2 = 2:i1
        m1=i1-1;
        m2=i2-1;
        bino(i1,i2) = bino(i1,i2-1)*(m1-m2+1)/m2;
    end
end
fac(1)=1;
for i1=2:2*mmax
    fac(i1)=fac(i1-1)*(i1-1);
end
for i1 = 0:mmax+1
    ga(i1+1)=gamma(i1+1/2);
end

% In contrast to Martin.m we now carry out an anti-Kenny
% average. This affects the definition of F which changes
% from 
% F_{mi} = exp(-a_i^2)*a_i^(2m)/m!

%%%% Pure incoming coherent state case %%%
Fpu=[];
    for i1=1:mmax;
        Fpu(i1, 1:imax)=exp(-Terg(1,1:imax)).*Terg(1,1:imax).^(i1-1)./factorial(i1-1);
    end

    
%
% to and average of this expression with a Gaussian around
% a_i with probability density proportional to exp(-(a-a0)^2/2/s/s))
%
% Index m is photon number in Fock state expansion (should be >amax^2*e)
% Index i numbers measured alpha values
tic
F = zeros(mmax,imax);
for i1 = 1:mmax
%    i1/mmax % Show where I am 
    for i2 = 1:imax
        s = max(sqrt(0.02*erg(1,i2)),0.0012); % Standard deviation, treat 
                                              % int=0 point separately
        F(i1,i2) = 0;
        m = i1-1;
        for r = 0:m
%            c1 = factorial(2*m)/factorial(2*r)/factorial(2*(m-r));
            c1 = bino(2*m+1,2*r+1);
            c2 = (2*s*s/(1+2*s*s))^(r)*sqrt(2/(1+2*s*s));
            c3 = (sqrt(erg(1,i2))/(1+2*s*s))^(2*(m-r));
%            F(i1,i2) = F(i1,i2) + c1*c2*c3*gamma(r+1/2);
            F(i1,i2) = F(i1,i2) + c1*c2*c3*ga(r+1);
        end
%        F(i1,i2) = F(i1,i2)*(sqrt(1/2/pi)/factorial(m))*exp(-erg(1,i2)/(1+2*s*s));
        F(i1,i2) = F(i1,i2)*(sqrt(1/2/pi)/fac(i1))*exp(-erg(1,i2)/(1+2*s*s));
    end
end
toc


%%%% to get rid of NaN  %%%
Fpu(isnan(Fpu))=0;
F(isnan(F))=0;














% Now we determine the POVM, m-th column gives the 
% POVM elements for given number of photons coming in
% Columns need to add to 1


povm = sdpvar(si(1)-1,mmax);
Tpovm = sdpvar(si(1)-1,mmax);

tic

m=zeros(mmax,mmax);
for j=1:mmax
    m(j,j)=1;
end



%%% to improve ill conditionning %%

totinvT = (Terg(2:si(1),:) - Tpovm*Fpu);
%totinv = (erg(2:si(1),:)*pinv(Fpu) - povm)*m;
totinv = (erg(2:si(1),:) - povm*Fpu);

CON = set(Tpovm(:) > 0);
CON = CON + set(sum(Tpovm)==1);
solvesdp(CON,norm(totinvT,'fro'))
Tpovmend = double(Tpovm);

CON = set(povm(:) > 0);
CON = CON + set(sum(povm)==1);
solvesdp(CON,norm(totinv,'fro'))
povmend = double(povm);


figure(51)
bar3(Tpovmend(1:9,1:12)')
title('reconstructed Loopy Loss')

figure(52)
bar3(povmend(1:9,1:12)')
title('reconstructed Measured Detector Feb_set2')

