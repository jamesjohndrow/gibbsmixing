clear; clc; close all;
if ismac
    boxname = 'Dropbox (Personal)';
else
    boxname = 'Dropbox';
end

[~, name] = system('hostname');
isshark = ~isempty(strfind(name,'sailfish'));

if isshark
    addpath(strcat('/ytmp/JEJ/',boxname,'/Utilities/Matlab Functions'));
    addpath(strcat('/ytmp/JEJ/',boxname,'/Utilities/Matlab Functions/Tensor'));
else
    addpath(strcat('~/',boxname,'/utilities/Matlab Functions'));
    addpath(strcat('~/',boxname,'/utilities/Matlab Functions/Tensor'));
end
usecanned = false; makefigures = false; saving = true;
nmc = 10000; burn = 0; thin = 1;
center = true; 


load('Data/convdat.mat');
Ytr = full(Ytr);
Ytr = Ytr(Ntr>0,:);
Ntr = Ntr(Ntr>0);
%Ytr = Ytr(:,[46 62]);
Ytr = Ytr(:,46);


disp_int = 1000;
n = length(Ntr); 
N = Ntr; Y = Ytr;

mu = log(((1+Y)./N)./(1-(1+Y)./N));
mu(N==1) = mean(mu(N>1));
mu0 = mean(mu);

ACC = zeros(n,1);

Tau2 = 5; Tau = sqrt(5);
NY = N-Y;

MU = zeros((nmc-burn)./thin,n);
TAU = zeros((nmc-burn)./thin,1);
MU0 = zeros((nmc-burn)./thin,1);
v = 9.*ones(n,1);
mu00 = -12; tau20 = 36; 
a0 = 0; b0 = 0;
tic;
for t=1:nmc+burn
    if mod(t,disp_int) == 0 %|| (mod(t,10) == 0 && t<100)
        disp(num2str(t)); toc; tic;
        if t > burn && makefigures
            if t <= 500; tstart = 1; else tstart = 500; end
            figure(3);plot(MU(tstart:(t-burn)/thin-1,1));drawnow;
            figure(4);plot(MU(tstart:(t-burn)/thin-1,2));drawnow;
            figure(6);plot(TAU(tstart:(t-burn)/thin-1));drawnow;
            figure(5);hist(ACC./t,30);drawnow;
            figure(7);plot(MU0(tstart:(t-burn)/thin-1));drawnow;
            figure(8);hist(mu,30);drawnow;
        end
    end


%     if t < 100 || mod(t,10) == 0
%         opts = optimoptions('fminunc','display','off','gradobj','on','hessian','on','TolFun',1e-6,'TolX',1e-6);
%         [muhats,Hs] = bayeslogitfit(Y,N,mu,sqrt(Tau2),mu0.*ones(n,1),opts,usecanned);
%     end
%     if usecanned
%         muprop = muhats + sqrt(v./diag(Hs)).*trnd(5);
%     else
%         muprop = muhats + sqrt(v./Hs).*trnd(5);
%     end
    muprop = normrnd(mu,1);

    if center
        pprop = exp(muprop)./(1+exp(muprop));
        pcurr = exp(mu)./(1+exp(mu));
        logar = Y.*log(pprop)+(N-Y).*log(1-pprop)-Y.*log(pcurr)-(N-Y).*log(1-pcurr);
        logpr = -.5./(Tau2).*(mu0-muprop).^2+.5./(Tau2).*(mu0-mu).^2;
    else
        pprop = exp(muprop+mu0)./(1+exp(muprop+mu0));
        pcurr = exp(mu+mu0)./(1+exp(mu+mu0));
        logar = Y.*log(pprop)+(N-Y).*log(1-pprop)-Y.*log(pcurr)-(N-Y).*log(1-pcurr);
        logpr = -.5./(Tau2).*(muprop).^2+.5./(Tau2).*(mu).^2;
    end

    accrat = rand(n,1)<exp(logar+logpr);
    mu(accrat) = muprop(accrat);
    ACC = accrat+ACC;
    
    if center
        res = mu - mu0;
        ssr = sum(res.^2);
        %Tau2 = 1./gamrnd((n-1)./2,2./(ssr));
        Tau2 = 1./gamrnd((n+a0)./2,2./(ssr+b0));
        Tau = sqrt(Tau2);
        
        s = (n./Tau2+1./tau20)^(-1);
        m = s.*(sum(mu)./Tau2+mu00/tau20);
        mu0 = normrnd(m,sqrt(s));
    else 
        ssr = sum(mu.^2);
        %Tau2 = 1./gamrnd((n-1)./2,2./(ssr));
        Tau2 = 1./gamrnd((n+a0)./2,2./(ssr+b0));
        Tau = sqrt(Tau2);
        
        mu0prop = normrnd(mu0,.5);
        pprop = exp(mu+mu0prop)./(1+exp(mu+mu0prop));
        pcurr = exp(mu+mu0)./(1+exp(mu+mu0));
        logar = sum(Y.*log(pprop)+(N-Y).*log(1-pprop)-Y.*log(pcurr)-(N-Y).*log(1-pcurr));
        logpr = -.5./(tau20).*(mu0-mu00).^2+.5./(tau20).*(mu0-mu00).^2;
        accrat0 = rand <exp(logar+logpr);
        if accrat0
           mu0 = mu0prop; 
        end
        
        
    end
    
     if t>burn && mod(t,thin)==0
          TAU((t-burn)./thin) = Tau2;
          MU((t-burn)./thin,:) = mu';
          MU0((t-burn)./thin) = mu0;
     end
end

AR = zeros(size(MU,2),1);
for j=1:size(MU,2)
    if mod(j,1000)==0
        disp(num2str(j));
    end
    tmp = ar(MU(5000:end,j),1);
    AR(j) = tmp;
end

if saving
    save('Outputs/binomhmmh_simple.mat');
    ctr = 0;maxidx = 0;
    for j=5000:5000:size(MU,2)
        ctr = ctr + 1;
        maxidx = j;
        MUsv = MU(:,j-5000+1:j);
        save(strcat('Outputs/binomhmmh_mu',num2str(ctr),'.mat'),'MUsv');
    end
    MUsv = MU(:,maxidx+1:end);
    save(strcat('Outputs/binomhmmh_mu',num2str(ctr+1),'.mat'),'MUsv');


    MUsml = MU(:,1:100:n);
    save('Outputs/binomhmmh_simple_small.mat','MUsml','MU0','AR','TAU');
end


