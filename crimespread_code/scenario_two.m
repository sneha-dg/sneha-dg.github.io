%This code refers to the graphs generated under Scenario 2 in the
%Research Paper
param.cities = 5;             %The number of cities in the network of populations
%[low crime, high crime]
param.numvar = 3;

% param.alpha = [1 0.95 0.77 0.65 0.6];     % crime contagion rates [Toronto Halton Durham Peel York] 
% param.beta =  [0.85  0.9 0.95 1 1];       % rate of C -> X
% param.gamma = [0.95 0.8 0.75 0.67 0.55];           % incarceration rate
% param.epsilon = [0.5 0.5 0.5 0.5 0.5];   % rate of I -> X
% param.delta = [0.45 0.85 0.39 0.40 0.6] ;    % rate of I -> C
% param.alphap = [0.6 0.8 0.5 0.3 0.3];         % natural propensity rate to do crime
%Make sure all values in the vectors add to 1
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
for a = linspace(0, 100,numparam)
    b = linspace(0, 100, numparam);
    c = linspace(100, 0, numparam);
    %final time
    tf = 10;
    
        %param.T controls the network structure. Comment out the network structure
        %you don't want to use in your simulation
        %Direct Transit to and from, assuming transit has no delays and is constant
        %between all regions
    %param.T = theta(i)*(ones(param.cities, param.cities) - eye(param.cities));
        %Movement between neighbouring regions
    %param.T = theta*diag(ones(param.cities-1,1),1)+ theta*diag(ones(param.cities-1,1), -1); 
        %People normally travel from Toronto back to other regions, assume Toronto is
        %Region 1 
    %param.T = theta(i)*[0 ones(1,param.cities-1); ones(param.cities-1,1) zeros(param.cities-1,param.cities-1)];
        %No travel
    %param.T = theta*zeros(param.cities,param.cities);
        %Travel in one direction in a circle
    %param.T  = theta*diag(ones(param.cities-1,1),1);
    
    % to vary different parameters
    param.alpha = 100*[a 1 0 0.05 0.75];     % crime contagion rates
    %param.beta =  [a  2 1 0.35 0];       % rate of C -> X
    %param.beta =  [a  0.25 0.25 0.35 0];
%     param.gamma = [a 0.5 2 0.05 0.05]; 
%     param.epsilon = [c(i) 0.2 0 0.1 0.22];   % rate of I -> X
%     param.delta = [b(i) 1 0.5 0.05 2] ;    % rate of I -> C
%     
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
xlabel('Alpha')
ylabel('Criminally active Population')
title('Criminal Population vs Alpha')
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