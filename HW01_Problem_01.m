function[]=HW01_Problem_01(text,n)
% The function takes two input arguments, the first input as a string 
% and the second an integer. The function should print the string on 
% separate lines, as many times as indicated by the second number input.
%
% You must enter the text with quotation mark, 
% for example HW01_Problem_01('hellow world',5)
% the text will be printed in separated lines.
% Yucong Cai

%-------------------------------------------------------------------------
if nargin ~=2 %nargin -- new concept! see help nargin
  disp('Need exactly two inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%Input comfirm, Structure comfirm
%-------------------------------------------------------------------------


for i=1:1:n
    fprintf('%s\n',text);
end

%-------------------------------------------------------------------------
