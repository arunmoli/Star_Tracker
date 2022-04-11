function  [q_est, P, dist_index, IFAIL] = q_method(W, V, nstar, dist_index) 

ERR_POS = 0.3;

[nstar, ~] = size(W);
dist_index = 0;
wgt = zeros(nstar,1);
for i=1:1:nstar-1
    dist_index = dist_index + W(i,:)*W(i+1,:)';
    wgt(i) = 1/nstar;
end
dist_index = dist_index + W(nstar-1,:)*W(1,:)';

B = W'.*diag(wgt)* V;
rho = trace(B); 
S = B + B';
Z = zeros(3,1);

for i=1:1:nstar
Z = Z + wgt(i).*getSuperCross(W(i,:))*V(i,:)';
end
lambda = 0;   
[q_est, lambda] = fastq(S, Z, rho);
A = getAfromQ(q_est);  

P = eye(3)-(A*V'*diag(wgt)*V*A');

dev_tot = 0.0 ;
for i=1:1:nstar 
   if (mag(i) >= 2.0) 
       std_dev = (4.5+0.7*(mag(i)-2.)+ERR_POS)/3600.*pi/180.;
   end
   if (mag(i) <  2.0) 
       std_dev = (4.5+ERR_POS)/3600.*pi/180. ; 
   end
   dev_tot = dev_tot + 1/(std_dev*std_dev) ;
end

% P(2,1) = (P(1,2)+P(2,1))/2.0 ;
% P(1,2) = P(2,1);
% P(3,1) = (P(1,3)+P(3,1))/2.0 ;
% P(1,3) = P(3,1);
% P(3,2) = (P(2,3)+P(3,2))/2.0 ;
% P(2,3) = P(3,2);
% P 
P = (P + P')/2;
[P, IFAIL] = SYMINV3(P);
if (IFAIL == 1)   
      error('Matrix inversion error at q_method()') ;
      return;
else
    
end
% P(2,1) = (P(1,2)+P(2,1))/2.0 ;
% P(3,1) = (P(1,3)+P(3,1))/2.0 ;
% P(3,2) = (P(2,3)+P(3,2))/2.0 ;
% P(1,3) = P(3,1);
% P(1,2) = P(2,1);
% P(2,3) = P(3,2);

%     for i=1:1:3
%         for j=1:1:3
%             P(i,j) = P(i,j)/dev_tot ; % P_dthe_dthe = 4*P_dQ_dQ */
%         end
%     end
%     
    P = P/dev_tot;

P_range = zeros(3,1);    
P_range(1) = sqrt(4*P(1,1))*180/pi*3600;
P_range(2) = sqrt(4*P(2,2))*180/pi*3600;
P_range(3) = sqrt(4*P(3,3))*180/pi*3600;    
        
end 