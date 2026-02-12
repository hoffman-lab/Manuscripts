%% Introduction
%
% This demo introduces usage of 'TypicalSparsityMeasures' on a data
% sequence for evaluating the sparsity
% -------------------------------------------------------------------------
% This case shows the evaluation of data sequence sparsity using a set of
% numerical simulation signals.
%
% -------------------------------------------------------------------------
% Author: Bingyan Chen, Fengshou Gu
% Date of revision: 29/03/2023
% -------------------------------------------------------------------------


clc;clear
close all


%% 1. generate a set signals with different sparsity

N = 1000;
X = zeros(N,100);
nn = zeros(1,100);
for i=1:100
    x=zeros(N,1);
    X(1:1*i:N,i)=1;
    nn(i)=length(find(X(:,i)==1)); % number of pulses
end

% Display three signals with different sparsity

figure(1);
subplot(3,1,1);plot(1:N,X(:,90));
xlabel('Samples');ylabel('Amplitude');
title(['Number of Pulses: ',num2str(nn(90))]);
subplot(3,1,2);plot(1:N,X(:,50));
xlabel('Samples');ylabel('Amplitude')
title(['Number of Pulses: ',num2str(nn(50))]);
subplot(3,1,3);plot(1:N,X(:,20));
xlabel('Samples');ylabel('Amplitude')
title(['Number of Pulses: ',num2str(nn(20))]);


%% 2. traditional sparsity measures (the value range is greater than 1)
% such as Kurtosis, Negentropy, L2/L1 norm

Index=[];%nn=[];opt=[];
IndStr={'Kurt','NE','L2/L1'};
figure(2),clf,
for j=1:length(IndStr)
    opt.crit=IndStr{j}; opt.a=0; opt.p=0;
    for i=1:100
        x = X(:,i);
        Index(i,j) = TypicalSparsityMeasures( x, opt );
    end
end

plot(nn,Index,'linewidth',1.5);grid on
legend(IndStr)
xlabel('Number of Pulses')
ylabel('Sparsity Measures')


%% 3. sparsity measures (the value range from 0 to 1)

Index=[];%nn=[];opt=[];
IndStr={'DN','HI','MSI','GI','GI2','GI3'}; %,'GI2','GI3','GGI','PFGI2','PFGI3'}
figure(3),clf,
for j=1:length(IndStr)
    opt.crit=IndStr{j}; opt.a=0; opt.p=0;
    for i=1:100
        x = X(:,i);
        Index(i,j) = TypicalSparsityMeasures( x, opt );
    end
end

plot(nn,Index,'linewidth',1.5);grid on
legend(IndStr)
xlabel('Number of Pulses')
ylabel('Sparsity Measures')


%% 4. novel sparsity measures with adjustable parameter (the value range from 0 to 1)

Index=[];%nn=[];opt=[];
IndStr={'GGI','FGGI','PFGI1','PFGI2','PFGI3'};
paraA=[0.5,0.8,0,0,0];
paraP=[0,2,0.5,1.5,2];
figure(4),clf,
for j=1:length(IndStr)
    opt.crit=IndStr{j}; opt.a=paraA(j); opt.p=paraP(j);
    for i=1:100
        x = X(:,i);
        Index(i,j) = TypicalSparsityMeasures( x, opt );
    end
end

plot(nn,Index,'linewidth',1.5);grid on
legend(IndStr)
xlabel('Number of Pulses')
ylabel('Sparsity Measures')

