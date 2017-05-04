function[v_main_out,n_out, m_out,h_out,Ik,In,Il]=main_calculation_multi(n_main_in,m_main_in,h_main_in,dt,v_main_in,I_main_in)



C_m = 1;        %uF / cm^2
G_k = 36;     %S / cm^2
V_k = -12;      %mV
G_na = 120;     %S / cm^2
V_na = 115;     %mV
G_leak = .3;    %mS / cm^2
V_leak = 10.6;  %mV
 

    
n_main=n_main_in;
m_main=m_main_in;
h_main=h_main_in;
v_in=v_main_in;
I_in=I_main_in;
 
dn_rk=alpha_n(v_in)*(1-n_main)-beta_n(v_in)*n_main;
n_rk=n_main+dn_rk*dt;
 
 
dm_rk=alpha_m(v_in)*(1-m_main)-beta_m(v_in)*m_main;
m_rk=m_main+dm_rk*dt;
 
 
dh_rk=alpha_h(v_in)*(1-h_main)-beta_h(v_in)*h_main;
h_rk=h_main+dh_rk*dt;
 
 
Ik_rk=G_k*(n_rk^4)*(v_in-V_k);
 
In_rk=G_na*(m_rk^3)*h_rk*(v_in-V_na);
 
Il_rk=G_leak*(v_in-V_leak);
 
 
 
dV_rk=(1/C_m)*(I_in-Ik_rk-In_rk-Il_rk);
 
V_rk=v_in+0.5*dV_rk*dt;
 
 
 
%------------------------------Rk mediate----------------------------------------
 
 
 
dn=alpha_n(V_rk)*(1-n_main)-beta_n(V_rk)*n_main;
n_main=n_main+dn*dt;
 
 
dm=alpha_m(V_rk)*(1-m_main)-beta_m(V_rk)*m_main;
m_main=m_main+dm*dt;
 
 
dh=alpha_h(V_rk)*(1-h_main)-beta_h(V_rk)*h_main;
h_main=h_main+dh*dt;
 
Ik=G_k*((n_main)^4)*(V_rk-V_k);
 
In=G_na*((m_main)^3)*h_main*(V_rk-V_na);
 
Il=G_leak*(V_rk-V_leak);
 
dV=(1/C_m)*(I_in-Ik-In-Il);
v_main_out=v_in+dV*dt;

n_out=n_main;
m_out=m_main;
h_out=h_main;




function[an_out]=alpha_n(v_suv_in)
an_out=(0.01*(10-v_suv_in))/(exp((10-v_suv_in)/10)-1);
end
 
 
function[bn_out]=beta_n(v_suv_in)
%bn_out=0.125*exp(-v_suv_in/20);
bn_out=0.125*exp(-v_suv_in/80);
end
 
function[am_out]=alpha_m(v_suv_in)
am_out=(0.1*(25-v_suv_in))/(exp((25-v_suv_in)/10)-1);
end
 
function[bm_out]=beta_m(v_suv_in)
bm_out=4*exp(-v_suv_in/18);
end
 
 
function[ah_out]=alpha_h(v_suv_in)
ah_out=0.07*exp(-v_suv_in/20);
end
 
 
function[bh_out]=beta_h(v_suv_in)
bh_out=1/(exp((30-v_suv_in)/10)+1);
end



end%mainfunciton
