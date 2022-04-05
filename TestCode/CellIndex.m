function Index = CellIndex(az, el)
% %int       n, j, No, cellsu, jmax ;
% %    double    el_dif, az_dif, de ;
% 
No = 14 ; % 15 ellinamtion zones 
% 
if (el <  90.00 && el >=  83.79)
    cellsu =  1 ;
elseif (el <  83.79 && el >=  71.38)
    cellsu =  5 ;
elseif (el <  71.38 && el >=  58.97)
    cellsu =  9 ;
elseif (el <  58.97 && el >=  46.55)
    cellsu = 13 ;
elseif (el <  46.55 && el >=  34.14)  
    cellsu = 17 ;
elseif (el <  34.14 && el >=  21.72)  
    cellsu = 21 ;
elseif (el <  21.72 && el >=   9.31)  
    cellsu = 25 ;
elseif (el <   9.31 && el >=  -3.10)  
    cellsu = 29 ;
elseif (el <  -3.10 && el >= -15.52)  
    cellsu = 27 ;
elseif (el < -15.52 && el >= -27.93)  
    cellsu = 23 ;
elseif (el < -27.93 && el >= -40.34)  
    cellsu = 19 ;
elseif (el < -40.34 && el >= -52.76)  
    cellsu = 15 ;
elseif (el < -52.76 && el >= -65.17)  
    cellsu = 11 ;
elseif (el < -65.17 && el >= -77.59)  
    cellsu =  7 ;
elseif (el < -77.59 && el >= -90.00)  
    cellsu =  3 ;
else
end
 
el_dif = (83.79+90.00)/No ;
az_dif = 360./cellsu ;

j    = ceil(az/az_dif - 0.5) ;
jmax = ceil(360/az_dif - 0.5) ;

de = 90 - el ;
 
if (de < 93.1001295051)
    n = 2*ceil(de/el_dif - 0.5) ;
else
    n = 2*No + 1 - 2*ceil(de/el_dif - 0.5) ;
end

if ( (n*n+j+1) == (n*n+jmax+1) )
    Index = n*n+1;
else
    Index = n*n+j+1 ; 
end
end
    