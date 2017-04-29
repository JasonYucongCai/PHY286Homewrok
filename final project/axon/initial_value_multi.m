function[v_initial,n_initial,m_initial,h_initial,I_ki,I_ni,I_Li]=initial_value_multi(n_in,m_in,h_in,dt,v_input,tolerance)
 
n_initial_temp_out=n_in;
m_initial_temp_out=m_in;
h_initial_temp_out=h_in;
v_initial_temp_out=v_input;
 
tolerance_value=tolerance*dt;
%tolerance_value=tolerance;
tolerance_comput=tolerance_value+1;

while tolerance_comput>=(tolerance_value)
 
v_initial_temp_in=v_initial_temp_out;
n_initial_temp_in=n_initial_temp_out;
m_initial_temp_in=m_initial_temp_out;
h_initial_temp_in=h_initial_temp_out;
 
[v_initial_temp_out,n_initial_temp_out,m_initial_temp_out,h_initial_temp_out,I_ki,I_ni,I_Li]=main_calculation_multi(n_initial_temp_in,m_initial_temp_in,h_initial_temp_in,dt,v_initial_temp_in,0);
 
 
v_initial_temp_compare=abs(v_initial_temp_in-v_initial_temp_out);
n_initial_temp_compare=abs(n_initial_temp_in-n_initial_temp_out);
m_initial_temp_compare=abs(m_initial_temp_in-m_initial_temp_out);
h_initial_temp_compare=abs(h_initial_temp_in-h_initial_temp_out);
 
tolerance_comput=max([v_initial_temp_compare,n_initial_temp_compare,m_initial_temp_compare,h_initial_temp_compare]);

fprintf('ddbug tolorance %d\n',tolerance_comput );

end
 
v_initial=v_initial_temp_out;
n_initial=n_initial_temp_out;
m_initial=m_initial_temp_out;
h_initial=h_initial_temp_out;
 
 
end%initial_value
