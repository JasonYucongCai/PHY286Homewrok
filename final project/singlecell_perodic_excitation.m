
% fp2(0.5,1,30,30,0.01,210,1E-15);

function [] = singlecell_perodic_excitation(I_start_time,I_duration,I_separation,I_induct_value,dt,tmax,tolerance)

%I_start_time, I_end_time,I_induct_value

N_IC=(tmax-I_start_time)/(I_duration+I_separation);

n_in=0;
m_in =0;
h_in=1;
v_input=0;


time_axis = linspace(0, tmax, round(tmax/dt));
%{
I_start_time;
I_end_time;
I_induct_value;

fp(0,0,1,0,2,4,0.01,0.01,10,0.01);


%[v_initial,n_initial,m_initial,h_initial]=initial_value(n_in,m_in,h_in,dt,v_input,tolerance)
%[v_main_out,n_out, m_out,h_out,Ik,In,Il]=
main_calculation(n_main_in,m_main_in,h_main_in,dt,v_main_in,I_main_in)
%}
N=length(time_axis);


%{
[set,n(1),m(1),h(1)]=initial_value(n_in,m_in,h_in,dt,v_input,tolerance);
v=ones(1,N)*set;
%}
[v(1),n(1),m(1),h(1),I_k,I_n,I_L]=initial_value(n_in,m_in,h_in,dt,v_input,tolerance);
time_axis=0;
Induct=[0];
output_Ik=I_k;
output_In=I_n;
output_Il=I_L;
h1(1)=0;
h2(1)=0;

TF_change_I_induct=0;

I_end_time=I_start_time+I_duration;

for i=1:N-1
    
    if (time_axis(i)>=I_start_time)&& (time_axis(i)<=I_end_time);
        I_induct=I_induct_value;
        %N_IC=(tmax-I_start_time)/(I_duration+I_separation);
        TF_change_I_induct=1;
        
    else
        I_induct=0;
        if TF_change_I_induct==1
        I_start_time=I_start_time+I_duration+I_separation;
        I_end_time=I_end_time+I_duration+I_separation;
        
        TF_change_I_induct=0;
        end
    end
    

    [v(i+1),n(i+1), m(i+1),h(i+1),I_k,I_n,I_L]=main_calculation(n(i),m(i),h(i),dt,v(i),I_induct);
    
  %  fprintf('ddbug main %d\n', v(i+1));
%{
       plot(time_axis, v, '.b');
       legend('V_plot');
       title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
       drawnow;
%}
    
    Induct=[Induct,I_induct];
    output_Ik=[output_Ik,I_k];
    output_In=[output_In,I_n];
    time_axis=[time_axis, time_axis(i)+dt];
    output_Il=[output_Il,I_L];
    
    
    Graph_Induct=Induct./1000;
     Graph_output_Ik=output_Ik./1000;
     Graph_output_In=output_In./1000;
     
     Graph_output_Il=output_Il./1000;
     Graph_output_compare=Graph_output_In+Graph_output_Ik;
    
     
      Graph_output_add=Graph_output_compare+Graph_output_Il;
      
      
      h1=[h1,sum(Graph_Induct)];
      h2=[h2,sum(Graph_output_add)];
%     Graph_Induct=Induct;
%     Graph_output_Ik=output_Ik;
%     Graph_output_In=output_In;
if mod(i,1000)==1
figure (24)
plot(time_axis,h1,'.b',time_axis,h2,'.k');
       legend('Induct current Integrate','Added current Integrate');
       title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
       drawnow;


    figure (21)
       plot(time_axis,v, '.b');
       legend('V_plot');
       title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
       drawnow;
    figure(22)
        plot(time_axis,Graph_Induct,'.k',time_axis,Graph_output_In,'.y',time_axis,Graph_output_Ik,'.g');
       legend('Induct current','Na current','Ik current');
       %ylim([-0.01 0.02]);
       title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
       drawnow;
       
    figure(23)
    plot(time_axis,Graph_output_Il, '.b',time_axis,Graph_output_compare, '.y',time_axis,Graph_output_add,'.k');
    legend('Leaking current','Sum of Na, K current','Added curent');
    title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
    drawnow;
end



end

figure (24)
plot(time_axis,h1,'.b',time_axis,h2,'.k');
       legend('Induct current Integrate','Added current Integrate');
       title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
       drawnow;


    figure (21)
       plot(time_axis,v, '.b');
       legend('V_plot');
       title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
       drawnow;
    figure(22)
        plot(time_axis,Graph_Induct,'.k',time_axis,Graph_output_In,'.y',time_axis,Graph_output_Ik,'.g');
       legend('Induct current','Na current','Ik current');
       %ylim([-0.01 0.02]);
       title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
       drawnow;
       
    figure(23)
    plot(time_axis,Graph_output_Il, '.b',time_axis,Graph_output_compare, '.y',time_axis,Graph_output_add,'.k');
    legend('Leaking current','Sum of Na, K current','Added curent');
    title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
    drawnow;
    
    
    
    
    
    
    
    
    


function[v_initial,n_initial,m_initial,h_initial,I_ki,I_ni,I_Li]=initial_value(n_in,m_in,h_in,dt,v_input,tolerance)
 
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
 
[v_initial_temp_out,n_initial_temp_out,m_initial_temp_out,h_initial_temp_out,I_ki,I_ni,I_Li]=main_calculation(n_initial_temp_in,m_initial_temp_in,h_initial_temp_in,dt,v_initial_temp_in,0);
 
 
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








function[v_main_out,n_out, m_out,h_out,Ik,In,Il]=main_calculation(n_main_in,m_main_in,h_main_in,dt,v_main_in,I_main_in)


%original
C_m = 1;        %uF / cm^2
G_k = 36;     %S / cm^2
V_k = -12;      %mV
G_na = 102;     %S / cm^2
V_na = 115;     %mV
G_leak = .3;    %mS / cm^2
V_leak = 10.6;  %mV

% a=0.95;
% b=0.38;
% 
% a=1;
% b=1;
% C_m = 1;        %uF / cm^2
% G_k = 36*a;     %S / cm^2
% V_k = -12*b;      %mV
% G_na = 120*a;     %S / cm^2
% V_na = 115*b;     %mV
% G_leak = 0.3*a;    %mS / cm^2
% V_leak = 10.6*b;  %mV
    
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



end%total funcion
