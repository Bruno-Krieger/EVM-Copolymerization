function phi = MaLe_Inner(theta,F1,f10,Xw,s,Mw,V,c)

% Written by Bruno V. Krieger 2020
% Mayo-Lewis instantaneous copolymer composition model
% Should be used alongside EVM.m for the calculation of the objetive
% function values for the estimation of the reactivity ratios in a
% copolymerization problem. Uses global variables. For it to be used
% replace the name EVM_Inner by MaLe_Inner in the fmincon function 
% call in the EVM.m file. If you want to use SCE go into the 
% sceua_n.m and cceua_n.m files and replace the name EVM_Inner by 
% MaLe_Inner where it appears, there are multiple locations.

global ests ints

xi = [f10 F1];
tol = 1e-6;
k = 1;
r1 = theta(1);
r2 = theta(2);

while 1
    clearvars F1_
    F1_ = zeros(size(F1,1),1);
    for i = 1:c
        idx = (1+sum(s(1,1:i-1),2):sum(s(1,1:i),2));
        f = xi(1+sum(s(1,1:i-1),2),1,end);
        F1_(idx,1) = (r1*f^2 + f*(1-f))/(r1*f^2 + 2*f*(1-f) + r2*(1-f)^2);
    end
    
    clearvars g B
    for i = 1:size(F1,1)
        g(:,:,i) = ((xi(i,2,end)) - (F1_(i,1)));
        
% %       B calculation for multiplicative error:
%         f = xi(i,1,end);
%         B(:,1,i) = f*(-r2*(1-f)*((2*r1-1)*f+1)-r1*f^2)/(f*((r1-2)*f+2)+r2*(1-f)^2)^2;
%         B(:,2,i) = (xi(i,2,end));

%       B calculation for additive error:
        f = xi(i,1,end);
        B(:,1,i) = (-r2*(1-f)*((2*r1-1)*f+1)-r1*f^2)/(f*((r1-2)*f+2)+r2*(1-f)^2)^2;
        B(:,2,i) = 1;

    end
    
    for i = 1:size(F1,1)
        C(:,:,i) = B(:,:,i)*V*B(:,:,i)';
        b(:,:,i) = (g(:,1,i)+B(:,:,i)*([(f10(i,:)) (F1(i,:))]-(xi(i,:,k)))');
        S(:,:,i) = chol(C(:,:,i));
        h(:,:,i) = S(:,:,i)'\b(:,:,i);
        t(:,:,i) = S(:,:,i)\h(:,:,i);
        xi(i,:,k+1) = ([(f10(i,:)) (F1(i,:))] - (V*B(:,:,i)'*t(:,:,i))');
        phi(1,i) = h(:,:,i)'*h(:,:,i)/2;
    end
    
    if abs(xi(:,:,k+1)-xi(:,:,k)) < tol
        break
    end
    
    k = k+1;
end

phi = sum(phi);
ests.phi = cat(1,ests.phi,phi);
ests.theta = cat(1,ests.theta,theta);
ints.g = cat(4,ints.g,g);
ints.xi = cat(4,ints.xi,(xi(:,:,end)));
ints.B = cat(4,ints.B,B);
end