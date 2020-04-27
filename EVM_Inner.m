function phi = EVM_Inner(theta,F1,f10,Xw,s,Mw,V,c)

% Written by Bruno V. Krieger 2019/2020
% Should be used alongside EVM.m for the calculation of the objetive
% function values for the estimation of the reactivity ratios in a
% copolymerization problem. Uses global variables.

global ests ints
xi = (F1);
tol = 1e-6;

ODEoptions = odeset('RelTol',1e-10,'AbsTol',1e-12);
k = 1;
intp = 1000;

Xn = Xw.*(Mw(1)*f10+(1-f10)*Mw(2))./(Mw(1)*(F1(:,1,end))+(1-(F1(:,1,end)))*Mw(2));
    
for i = 1:c
    sol(i) = ode45(@(Xn_,f1_) odefun(Xn_,f1_,theta),linspace(0,.9999,intp), ...
                                     f10(1+sum(s(1,1:i-1),2),1),ODEoptions);
end

while 1
    clearvars F1_c f1
    F1_c = zeros(size(F1,1),1);
    f1 = zeros(size(f10,1),1);
    for i = 1:c
        idx = (1+sum(s(1,1:i-1),2):sum(s(1,1:i),2));
        f1(idx,1) = deval(sol(i),Xn(idx,1))';
        Xn_ = Xn(idx,1);
        F1_c(idx,1) = (f10(1+sum(s(1,1:i-1),2),1)-f1(idx,1).*(1-Xn_))./Xn_;
    end
    
    clearvars g B Xn_
    for i = 1:size(F1,1)
        g(:,:,i) = ((xi(i,1,end)) - (F1_c(i,1)));
%         B(:,:,i) = (xi(i,1,end)); % B vector for multiplicative error
        B(:,:,i) = 1; % B vector for additive error
    end
    
    for i = 1:size(F1,1)
        C(:,:,i) = B(:,:,i)*V*B(:,:,i)';
        b(:,:,i) = (g(:,1,i)+B(:,:,i)*((F1(i,:))-(xi(i,:,k)))');
        S(:,:,i) = chol(C(:,:,i));
        h(:,:,i) = S(:,:,i)'\b(:,:,i);
        t(:,:,i) = S(:,:,i)\h(:,:,i);
        xi(i,:,k+1) = ((F1(i,:)) - (V*B(:,:,i)'*t(:,:,i))');
        phi(1,i) = h(:,:,i)'*h(:,:,i)/2;
    end    
    
    if abs(xi(:,:,k+1)-xi(:,:,k)) < tol
        Xn = Xw.*(Mw(1)*f10+(1-f10)*Mw(2))./(Mw(1)*(xi(:,1,end))+(1-(xi(:,1,end)))*Mw(2));
        break
    end
    
    Xn = Xw.*(Mw(1)*f10+(1-f10)*Mw(2))./(Mw(1)*(xi(:,1,end))+(1-(xi(:,1,end)))*Mw(2));
    k = k+1;
end

phi = sum(phi);
ests.phi = cat(1,ests.phi,phi);
ests.theta = cat(1,ests.theta,theta);
ints.g = cat(4,ints.g,g);
ints.xi = cat(4,ints.xi,(xi(:,:,end)));
ints.B = cat(4,ints.B,B);
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