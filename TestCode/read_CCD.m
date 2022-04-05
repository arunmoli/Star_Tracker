function [t, x, y, xmag,m_star, nstar] = read_CCD(t, FOVs_meas, Ns)

% sm = 0 ; % 0:no problem, 1:bad data, 2:no data

% READ THROUGH FOVs_meas FOR MATCHING t
% RETURN ROW NUMBER r OR return
for r = 1:size(FOVs_meas,1)
   %End of FOVs_meas
    if r == size(FOVs_meas,1)
       printf("No more CCD measurements (read_CCD:0).\n");
       nstar = -999;
       return;
   %Equivalent t found within FOVs_meas
   elseif t == FOVs_meas(r,1) && r < size(FOVs_meas,1)
       r = r + 1;
       break;
   end
end

t = FOVs_meas(r,1);

ii = 2; %Star's starting index within FOVs_meas
for i = 1:Ns
    num1 = FOVs_meas(r,ii);
    num2 = FOVs_meas(r,ii+1);
    num3 = FOVs_meas(r,ii+2);
    cx(i) = num1;
    cy(i) = num2;
    cm(i) = num3;
    mag(i) = num3;
    ii = ii+3;
end

nstar = 1;

for i = 1:Ns
    for j = 1:3
        m_star(i).L(j) = 0.0;
        m_star(i).mag = mag(i);
    end
    
    x(i) = 0;
    y(i) = 0;
    xmag(i) = mag(i);
    phi(nstar) = cx(i);
    lambda(nstar) = cy(i);
    mag(nstar) = mag(i);
    
    if mag(nstar) >= -2.0 && mag(nstar) < 7.5
        nstar = nstar +1;
    end
end



for i = 1:nstar
   x(i) = phi(i);
   y(i) = lambda(i);
   xmag(i) = mag(i);
   phi(i) = phi(i)*pi/180.0; % degree to radian
   lambda(i) = lambda(i)*pi/180.0; % degree to radian
   s(1) = phi(i);
   s(2) = lambda(i);
   m_star(i).L(3) = 1./sqrt(1+s(1)*s(1)+s(2)*s(2));
   m_star(i).L(1) = s(1)*m_star(i).L(3) ;
   m_star(i).L(2) = s(2)*m_star(i).L(3) ;
   m_star(i).mag = mag(i);  
end

% sort(nstar, x, y, xmag);
for i = 1:nstar
    for j = i+1:nstar
        if m_star(j).mag < m_star(i).mag
            temp = m_star(i);
            m_star(i) = m_star(j);
            m_star(j) = temp;
            
            tmp_x = x(i);
            tmp_y = y(i);
            tmp_mag = xmag(i);
            x(i) = x(j);
            y(i) = y(j);
            xmag(i) = xmag(j);
            x(j) = tmp_x;
            y(j) = tmp_y;
            xmag(j) = tmp_mag;            
        end
    end
end

return;
end

