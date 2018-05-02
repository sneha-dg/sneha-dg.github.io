%This code refers to the graphs generated under Expanding the Model in the 
%Research Report
param.cities = 5;             %The number of cities in the network of populations
%[low crime, high crime]
param.numvar = 3;


param.alpha = 100*[2 1 0 0.05 0.75];     % crime contagion rates
param.beta =  [0.25  2 1 0.35 0];       % rate of C -> X
param.gamma = [1 0.5 2 0.05 0.05];           % incarceration rate
param.epsilon = [2 0.2 0 0.1 0.22];   % rate of I -> X
param.delta = [2 1 0.5 0.05 2] ;    % rate of I -> C
param.alphap = [2 0.02 0.03 0.15 1];   % natural propensity rate to do crime

    %for varying params
% param.alpha = 100*[1 1 1 1 1];     % crime contagion rates
% param.beta =  [1 1 1 1 1];       % rate of C -> X
% param.gamma = [1 1 1 1 1];           % incarceration rate
% param.epsilon = [1 1 1 1 1];   % rate of I -> X
% param.delta = [1 1 1 1 1] ;    % rate of I -> C
% param.alphap = [1 1 1 1 1];   % natural propensity rate to do crime

theta = 50;
%param.T = theta*zeros(param.cities,param.cities);
param.T = theta*(ones(param.cities, param.cities) - eye(param.cities));
%param.T = theta*[0 ones(1,param.cities-1); ones(param.cities-1,1) zeros(param.cities-1,param.cities-1)];
%param.T = theta*(ones(param.cities, param.cities) - eye(param.cities));
P0 = [10 25 45 45 50;
      40 30 25 35 40
      50 45 30 20 10]/500;

i = 1;

numparam=100;
max_param =  100;
final_C = NaN(numparam, param.cities);

ep0 = 50; %initial conditions for epsilon1, and delta1
d0 = 1;

k1 = -1;
k2 = 2;
for gamma = linspace(0, 100,numparam)

    %final time
    tf = 10;
    eps1 = ep0 + k1 * gamma;
    d1 = d0 + k2 * gamma;
    if d1 < 0
        d1 = 0;
    end
    if eps1 < 0
        eps1 = 0;
    end
    
    param.gamma = [gamma 0.5 2 0.05 0.05];           % incarceration rate
    param.epsilon = [eps1 0.2 0 0.1 0.22];   % rate of I -> X
    param.delta = [d1 1 0.5 0.05 2] ;    % rate of I -> C
    
    [t,XCI] = ode23s( @(t,x) equations(t,x,param), [0,tf], P0(:) );
    C1 = XCI(:,2); C2 = XCI(:,5); C3 = XCI(:, 8); C4 = XCI(:, 11); C5 = XCI(:, 14);
    final_C(i,1) = C1(end);
    final_C(i,2) = C2(end);
    final_C(i,3) = C3(end);
    final_C(i,4) = C4(end);
    final_C(i,5) = C5(end);
    i = i + 1;
    disp([num2str(i-1) '/' num2str(numparam)])
end

figure()
hold on
plot(linspace(0,max_param, numparam), final_C(:,1))
plot(linspace(0,max_param, numparam), final_C(:,2))
plot(linspace(0,max_param, numparam), final_C(:,3))
plot(linspace(0,max_param, numparam), final_C(:,4))
plot(linspace(0,max_param, numparam), final_C(:,5))
xlabel('Gamma')
ylabel('Criminally active Population')
title('Criminal activity when k2 > 0 and k1 < 0')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
legend('1', '2', '3', '4', '5')  
hold off

function dxdt = equations(t,x,param)
% function that computes the right-hand side of the differential equation
    Y = reshape(x,3,param.cities);
    X = Y(1,:);
    C = Y(2,:);
    I = Y(3,:); % extract the variables
    
    dxdt =  [param.beta .* C - param.alpha .* X .* C - param.alphap .* X + param.epsilon .* I+ (X*param.T) - sum(param.T,1) .* X; ... dX/dt
           -param.beta  .* C + param.alpha .* X .* C + param.alphap .* X + param.delta .* I - param.gamma .* C + (C*param.T) - sum(param.T,1) .* C; ... dC/dt
            param.gamma .* C - (param.delta + param.epsilon) .*I        ... dI/dt
           ];
	dxdt = dxdt (:);
end