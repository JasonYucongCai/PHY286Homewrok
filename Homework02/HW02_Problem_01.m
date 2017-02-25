
function[]=HW02_Problem_01(N0,a,b,T,n)

%check
%-------------------------------------------------------------------------
if nargin ~=5 %nargin -- new concept! see help nargin
  disp('Need exactly five inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------
if n<0 || T<0 || N0<0
    error('check your input, the population, time, or step could not be negative.');
end
%Input comfirm, legal adjusting.
%-------------------------------------------------------------------------

xaxis=linspace(1,n,n);
delta=T/n;
N(1)=N0;


%setup
%-------------------------------------------------------------------------


for i=1:n;

%RK
Nrk(i)=
t=



%loop RK
%-------------------------------------------------------------------------


%calculate


end
