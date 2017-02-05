function [  ] = HW01_Problem_08(T,vmin,vmax)
%A wind chill factor (WCF) describes how cold it “feels” for a given 
%temperature T, in Fahrenheit, and a given wind speed V(in miles per hour). 
%
%The function accepts three inputs: the temp in Fahrenheit (T), a minimum 
%wind speed (Vmin), and maximum wind speed (Vmax). Use a for loop to 
%compute and appropriately display the wind chill factor (WCF) for 100 
%speeds, evenly spaced across the the wind speed range specified and plot 
%your result as WCF vs. wind speeds.
%
%The function accepts three inputs: the temp in Fahrenheit (T), a minimum 
%wind speed (Vmin), and maximum wind speed (Vmax).
%You may enter the temperature in Fahreneit, Wind Speed in miles per hour.
%For example: HW01_Problem_08(15,2,160)
%Yucong Cai

%---------------------------------------------------------------------
if nargin ~=3 %nargin -- new concept! see help nargin
  disp('Need exactly three inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%input confirm, structure comfirm.
%-------------------------------------------------------------------------

if vmin<0
    disp('The minium wind speed must be equal or bigger than 0');
return;
else if vmin>=vmax
    disp('The minium wind speed must smaller than the maxmiun wind speed');
    return;
    end
end


if T< (-459.67)
   disp('The minium temeprature must be bigger than -459.67 Fahrenheit.'); 
return;
end



%input comfirm, legal comfirm.
%-------------------------------------------------------------------------
    




A=linspace(vmin,vmax,100);
   
for i=1:1:100
        A(1,i);
        WCF(1,i)=35.7+0.6*T-35.7*(A(1,i)^0.16)+0.43*T*(A(1,i)^0.16);


end





%Calculation
%-------------------------------------------------------------------------
figure
plot(A,WCF,'b--')

legend('Wind Chill Factor');

xlabel('Wind Speed(miles per hour)');
ylabel('WCF');

a=vmin;
b=vmax;
xlim([a b]);%range of x axis

title('Wind Chill Factor in given temperature');


%plots
%-------------------------------------------------------------------------
end


%-------------------------------------------------------------------------










