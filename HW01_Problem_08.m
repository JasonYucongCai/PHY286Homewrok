function [  ] = HW01_Problem_08(T,vmin,vmax)
%the function accepts three inputs: the temp in Fahrenheit (T), a minimum 
%wind speed (Vmin), and maximum wind speed (Vmax).
%You may enter the temperature in Fahreneit, Wind Speed in miles per hour.
%For example: HW01_Problem_08(15,2,160)


%---------------------------------------------------------------------
if nargin ~=3 %nargin -- new concept! see help nargin
  disp('Need exactly three inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%input confirm, structure comfirm.
%-------------------------------------------------------------------------





A=linspace(vmin,vmax,100);
   
for i=1:1:100
   A(1,i);
   WCF(1,i)=35.7+0.6*T-35.7*(A(1,i)^0.16)+0.43*T*(A(1,i)^0.16);


end

figure
plot(A,WCF,'b--')

legend('Wind Chill Factor');

xlabel('Wind Speed(miles per hour)');
ylabel('WCF');

a=vmin;
b=vmax;

xlim([a b]);
title('Wind Chill Factor in given temperature');

end

















