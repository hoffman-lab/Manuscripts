function [ Index ] = TypicalSparsityMeasures( x, opt )
% Input: x...signal to be evaluated (a vector)
%        opt...sparsity indicator name and its parameters
%        opt.crit..name of sparsity indicator
%        opt.a..parameter a of Generalized Gini Index (GGI), opt.a>0
%        opt.p..parameter p of Fully Generalized Gini Index (FGGI) and
%               Power function-based Gini Index 1, 2 and 3 (PFGI1, PFGI2 and PFGI3), opt.p>0
% Output: Index...sparsity indicator value
%--------------------------------------------------------------------------
% Note: (1) opt.a and opt.p are set to 0 when opt.crt={'Kurt','NE','DN','L2/L1','HI','GI','GI2','GI3'}
%       (2) opt.a is set to 0 when opt.crt={'PFGI1','PFGI2','PFGI3'}
%       (3) opt.p is set to 0 when opt.crt={'GGI'}
% Such as:
% opt.crt={'Kurt','NE','DN','L2/L1','HI','GI','GI2','GI3','GGI','FGGI','PFGI1','PFGI2','PFGI3'};
% opt.p=[0,0,0,0,0,0,0,0,0,3,3,3,3];
% opt.a=[0,0,0,0,0,0,0,0,2,3,0,0,0];
%**************************************************************************

% Example of a setting:
% opt.crit='Kurt'; opt.a = 0; opt.p=0;
% [ Index ] = SparsityMeasures( x, opt );

% Example 1 -------------
% Index=[];nn=[];
% IndStr={'NE','DN','L2/L1','GI','GI2'}; % 'Kurt','HI','MSI','GI3'
% figure(1),clf,
% for j=1:length(IndStr)
%     opt.crit=IndStr{j}; opt.a=0; opt.p=0;
%     for i=1:100
%         x=zeros(1000,1);
%         x(1:1*i:1000)=1;
%         Index(i,j)=TypicalSparsityMeasures( x, opt );
%         nn(i)=length(find(x==1));
%     end
% end
% 
% plot(nn,Index,'linewidth',1.5)
% grid on
% legend(IndStr)
% xlabel('Number of Pulses')
% ylabel('Sparse Measurers')

% Example 2 -------------
% Index=[];nn=[];
% IndStr={'GGI','FGGI','PFGI1','PFGI2','PFGI3'};
% paraA=[0.5,0.8,0,0,0];
% paraP=[0,2,0.5,1.5,2];
% figure(2),clf,
% for j=1:length(IndStr)
%     opt.crit=IndStr{j}; opt.a=paraA(j); opt.p=paraP(j);
%     for i=1:100
%         x=zeros(1000,1);
%         x(1:1*i:1000)=1;
%         Index(i,j)=TypicalSparsityMeasures( x, opt );
%         nn(i)=length(find(x==1));
%     end
% end
% 
% plot(nn,Index,'linewidth',1.5)
% grid on
% legend(IndStr)
% xlabel('Number of Pulses')
% ylabel('Sparse Measurers')
%
%--------------------------------------------------------------------------
% Authors: Bingyan Chen, Fengshou GU 2023
% Date of revision: 29/03/2023
%--------------------------------------------------------------------------
%
% Reference 1: FGGI
% MSSP A full generalization of the Gini index for bearing condition monitoring
% https://doi.org/10.1016/j.ymssp.2022.109998
% Reference 2: PFGI1, PFGI2, PFGI3
% SHM Power function-based Gini indices: New sparsity measures using power function-based 
% quasi-arithmetic means for bearing condition monitoring
% https://doi.org/10.1177/14759217221149745
% Reference 3: PFGI1
% JDMD IGIgram: An Improved Gini Index-Based Envelope Analysis for Rolling Bearing Fault Diagnosis
% https://doi.org/10.37965/jdmd.2022.65
% Reference 4: MSI
% MSSP Investigations on improved Gini indices for bearing fault feature characterization 
% and condition monitoring https://doi.org/10.1016/j.ymssp.2022.109165
%==========================================================================

N = length(x);
x = x(:);
%x = x(:) - mean(x);

%x = abs(hilbert(x)); % envelope signal
%x = abs(hilbert(x)).^2; % squared envelope signal

ax = sort(abs(x),'ascend');
%ax = sort(x,'ascend');
if strcmp(opt.crit,'Kurt') == 1
    %x = x(:) - mean(x);
    Index = mean(x.^4)/(mean(x.^2)^2);    % kurtosis
elseif strcmp(opt.crit,'NE') == 1
    nx = (x.^2)/mean(x.^2)+10^(-5);
    Index = mean(nx.*log(nx));  % negentropy
elseif strcmp(opt.crit,'DN') == 1
    Index = max(abs(x))/sqrt(sum(x.^2));  % D-norm
elseif strcmp(opt.crit,'L2/L1') == 1
    Index = sqrt(N)*sqrt(sum(x.^2))/sum(abs(x));  % L2/L1 norm
    %Index = (sqrt(N)*sqrt(sum(x.^2))/sum(abs(x))-1)/(sqrt(N)-1);    % normalized L2/L1 norm
elseif strcmp(opt.crit,'HI') == 1
    Index = (sqrt(N)-sum(abs(x))/sqrt(sum(x.^2)))/(sqrt(N)-1);  % Hoyer Index
elseif strcmp(opt.crit,'MSI') == 1
    %Index = mean(abs(x))/exp(mean(log(abs(x))));  % reciprocal of smoothness index
    Index = 1-exp(mean(log(abs(x)+10^(-5))))/mean(abs(x));  % modified smoothness index
    %Index = 1-exp(mean(log(abs(x))))/mean(abs(x));  % modified smoothness index
elseif strcmp(opt.crit,'GI') == 1
    weight = (N-1+0.5 : -1 :0.5)/N;
    Index = 1-2*weight*ax/sum(abs(x));  % Gini index
elseif strcmp(opt.crit,'GI2') == 1
    weight1 = (2*N-1 : -2 : 1)/(N^2);
    weight2 = (1 : 2 : 2*N-1)/(N^2);
    Index = 1-(weight1*ax)/(weight2*ax);  % Gini index 2
elseif strcmp(opt.crit,'GI3') == 1
    weight = (1 : 2 : 2*N-1)/(N^2);
    Index = 1-mean(abs(x))/(weight*ax);  % Gini index 3
elseif strcmp(opt.crit,'GGI') == 1
    weight = (2*N-1 : -2 : 1).^opt.a;
    nw = weight/sum(weight);
    Index = 1-nw*ax/mean(abs(x));  % generalized Gini index
elseif strcmp(opt.crit,'FGGI') == 1
    weight = (2*N-1 : -2 : 1).^opt.a;
    nw = weight/sum(weight);
    axp = ax.^opt.p;
    Index = 1-(nw*axp/mean(axp))^(1/opt.p);  % Fully generalized Gini index
elseif strcmp(opt.crit,'PFGI1') == 1
    weight = (2*N-1 : -2 : 1)/(N^2);
    axp = ax.^opt.p;
    Index = 1-(weight*axp/mean(axp))^(1/opt.p);  % Power function-based Gini index 1
elseif strcmp(opt.crit,'PFGI2') == 1
    weight1 = (2*N-1 : -2 : 1)/(N^2);
    weight2 = (1 : 2 : 2*N-1)/(N^2);
    axp = ax.^opt.p;
    Index = 1-((weight1*axp)/(weight2*axp))^(1/opt.p);  % Power function-based Gini index 2
elseif strcmp(opt.crit,'PFGI3') == 1
    weight = (1 : 2 : 2*N-1)/(N^2);
    axp = ax.^opt.p;
    Index = 1-(mean(axp)/(weight*axp))^(1/opt.p);  % Power function-based Gini index 3
else
    disp('parmeter is not set right!')
end
    
end




   