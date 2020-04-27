function out_data = EVM(Data,Mw,theta,e,s,m,lb,ub)

% Written by Bruno V. Krieger 2019/2020
% 
% Data must be a MATLAB structure containing:
%   f10 = [0.1; 0.1; 0.3; 0.3]; Column vector containing the initial mole fraction of each
%                               data set. The value repeats for every data point in each
%                               set
%   Xw = [0.01; 0.02; 0.03; 0.04]; Column vector containing the degree of mass conversion
%                                  for every data point for each data set
%   F1 = [0.23; 0.22; 0.41; 0.42]; Column vector containing the cumulative copolymer 
%                                  compositions of all data sets by order of initial mole 
%                                  fraction followed by order of conversion
%   Data example:
%    f10 = [0.0254;	   Xw = [0.0111;	F1 = [0.2711;
%           0.0254;          0.0296;          0.4504;
%           0.0254;          0.0669;          0.3009;
%           0.0254;          0.0995;          0.2586;
%           0.0254;          0.1480;          0.1762;
%           0.0596;          0.0371;          0.4832;
%           0.0596;          0.0843;          0.3907;
%           0.0596;          0.1102;    	  0.4234;
%           0.0596;          0.1511;          0.3882;
%           0.0596;          0.2061;          0.3009;
%           0.0882;          0.1012;          0.4410;
%           0.0882;          0.2164;          0.3922;
%           0.0882;          0.2961;          0.3091;
%           0.0882;          0.3943;          0.2291;
%           0.0882];         0.4664];         0.1925];
% In this example there are 3 f10 levels 0.0254, 0.0596, and 0.0882, these correspond to 
% three experiments with different monomer 1 concentrations. For each experiment the 
% different lines represent a diferent mass conversion level with an associated cumulative
% copolymer composition, these must be placed in order inside each inital mole fraction 
% group.
% 
% Other input arguments:
% Mw = [Monomer_1_Mw Monomer_2_Mw]; Line vector containing monomer molecular weights
% e = [e_F; e_f; e_X]; Column vector of associated measurement errors.
% theta = [.01 .3]; Line vector containing the initial reactivity ratio estimates.
% s = [15 15 15]; Line vector containing the number of data sets and how many measurements 
%                 there are in each one. In the example given we have 3 data sets, each  
%                 level has 15 data points.
% m = 1; Method code, can either be 1 for fmincon or 2 for SCE
% lb = 0*ones(size(theta,2),1)'; Lower bound line vector containing the lower bound 
%                                values for each theta used
% ub = 5*ones(size(theta,2),1)'; Upper bound line vector containing the upper bound 
%                                values for each theta used
% 
% Output data: Is a MATLAB structure contaning the relevant output result
% data. It contains: the final parameter estimate, f_theta; minimized
% objective function value, phi; the true value of the inputed variables,
% xi; values of the error function g, and its derivative B; and finally the
% corrected conversion values, Xn.
% 
% This function uses global variables internally.
% 

clear global
% Setting global variable spaces for result estimates and internal values:
global ests ints 
ests.phi = [];
ests.theta = [];
ints.g = [];
ints.xi = [];
ints.B = [];

% Variable reallocation:
F1 = Data.F1;
f10 = Data.f10;
Xw = Data.Xw;
% Variance co-variance matrix:
V = (diag(e)^2)/3;
% Number of f10 levels in the data:
c = size(s,2); 

%  Inner loop call:
switch m
    case 1
        % Using fmincon:
        options = optimoptions('fmincon');
        options = optimoptions(options,'Display', 'iter');
        options = optimoptions(options,'FunValCheck', 'off');
        options = optimoptions(options,'Algorithm', 'sqp');
        options = optimoptions(options,'FiniteDifferenceType', 'central');
        options = optimoptions(options,'UseParallel', false);
        options = optimoptions(options,'StepTolerance', 1e-10);
        options = optimoptions(options,'OptimalityTolerance', 1e-8);
        options = optimoptions(options,'MaxFunctionEvaluations', 1000);
        A = [];
        b = [];
        [theta(2,:),phi,~,output,~,grad,hessian] = fmincon(@(theta) ...
            EVM_Inner(theta,F1,f10,Xw,s,Mw,V,c),theta,A,b,[],[],lb,ub,[],options);
    case 2
        % Using SCE:
        maxn=10000;
        kstop=10;
        pcento=0.1; 
        peps=0.001;
        iseed=-1;
        iniflg=0;
        ngs = 20;
        [theta(2,:),phi] = sceua_n(theta,lb,ub,maxn,kstop,pcento,peps,ngs,iseed, ...
                                 iniflg,F1,f10,Xw,s,Mw,V,c);
        hessian = [];
    otherwise
        error('Method code must be between 1 for fmincon or 2 for SCE.')
end

min_phi = phi;
min_index = find(min_phi==ests.phi); % Finds the index of min phi
g = ints.g(:,:,:,min_index); % Uses min phi index to find the correct g
xi = ints.xi(:,:,:,min_index); % Uses min phi index to find the correct xi
B = ints.B(:,:,:,min_index); % Uses min phi index to find the correct B
Xn = Xw.*(Mw(1)*f10+(1-f10)*Mw(2))./(Mw(1)*(xi(:,:,end))+(1-(xi(:,:,end)))*Mw(2));

ODEoptions = odeset('RelTol',1e-10,'AbsTol',1e-12);
for i = 1:c
    [iXn,f1(:,i)] = ode45(@(Xn_,f1_) odefun(Xn_,f1_,theta(2,:)), ...
                                            linspace(0,.9999,1000), ...
                                            f10(1+sum(s(1,1:i-1),2),1),ODEoptions);
    F1_ctot(:,i) = (f1(1,i)-f1(:,i).*(1-iXn))./iXn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JCR plotting sub-fuction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = input('Would you like to plot the Joint Confidence Regions? (y/n) \n','s');
isOK = false;
while ~isOK
    switch p
        case 'y'
            isOK = true;
            if isempty(hessian)
                hessian = hessian_n_alt(@(theta) EVM_Inner(theta,F1,f10,Xw,s,Mw, ...
                                                           V,c),theta(2,:));
            end      
            JCR = error_ellipse(inv(hessian),theta(2,:),'conf',0.95);
            hold on
            grid on
            set(JCR,'LineStyle','-.','Color','k')
            sct = scatter(theta(2,1),theta(2,2),'s','MarkerEdgeColor','k','Marker', ...
                          'square');
            out_data.JCR = JCR;
        case 'n'
            isOK =true;
        otherwise
            p =  input('Anwser must be either y or n. \n','s');
    end
end

% Outputting:
out_data.f_theta = theta(2,:);
out_data.phi = phi;
out_data.g = g;
out_data.xi = xi;
out_data.B = B;
out_data.Xn = Xn;
out_data.iXn = iXn;
out_data.F1_ctot = F1_ctot;
if ~isempty(hessian)
    out_data.hessian = hessian;
else
    hessian = hessian_n_alt(@(theta) EVM_Inner(theta,F1,f10,Xw,s,Mw,V,c),theta(2,:));
    out_data.hessian = hessian;
end

clear global
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mayo-Lewis equation sub-function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df1dXn = odefun(x,y,theta)
r1 = theta(1);
r2 = theta(2);
f1 = y;
f2 = 1-f1;
Xn = x;
F11 = ( r1*f1^2 + f1*f2 )/( r1*f1^2 + 2*f1*f2 + r2*f2^2 );
% Vector of equations:
df1dXn = (f1-F11)/(1-Xn);
end