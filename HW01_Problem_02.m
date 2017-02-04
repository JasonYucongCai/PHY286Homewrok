function[]=HW01_Problem_02(text,n)
% you must enter the text with quotation mark, 
%for example HW01_Problem_01('hellow world',5)
%the text will be printed in one line
%Yuocng Cai

%-------------------------------------------------------------------------
if nargin ~=2 %nargin -- new concept! see help nargin
  disp('Need exactly two inputs. See help');
    %error('Need exactly two inputs. See help');
  return;
end
%-------------------------------------------------------------------------

for i=1:1:n
    fprintf('%s',text);
end