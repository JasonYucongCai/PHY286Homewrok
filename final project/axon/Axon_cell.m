

%Axon_cell(0,0,1,0,1,2,2,0.01,50,1E-15,12);
function [] = Axon_cell(n_in, m_in, h_in,v_input,I_start_time, I_end_time,I_induct_value,dt,tmax,tolerance,N_axon)


time_axis = linspace(0, tmax, round(tmax/dt));

N=length(time_axis);


[v_axon,n_axon,m_axon,h_axon,I_k_axon,I_n_axon,I_L_axon]=initial_value_multi(n_in,m_in,h_in,dt,v_input,tolerance);
% function[v_initial,n_initial,m_initial,h_initial,I_ki,I_ni,I_Li]
%     =initial_value_multi(n_in,m_in,h_in,dt,v_input,tolerance)
 

for i=1:N_axon
v(i,1)=v_axon;
n(i,1)=n_axon;
m(i,1)=m_axon;
h(i,1)=h_axon;
I_k(i,1)=I_k_axon;
I_n(i,1)=I_n_axon;
I_L(i,1)=I_L_axon;
%for j=1:N_axon
I_induct(i,1)=0;
%end
end

time_axis=0;



for i=1:N-1
    
    if (time_axis(i)>=I_start_time)&& (time_axis(i)<=I_end_time);
        I_induct_initial_compartment=I_induct_value;
    else
        I_induct_initial_compartment=0;
    end
    
    I_induct(1,i+1)= I_induct_initial_compartment;
    [v(1,i+1),n(1,i+1), m(1,i+1),h(1,i+1),I_k(1,i+1),I_n(1,i+1),I_L(1,i+1)]=main_calculation_multi(n(1,i),m(1,i),h(1,i),dt,v(1,i),I_induct_initial_compartment);
    
    I_induct(2,i)=I_k(1,i+1)+I_n(1,i+1)+I_L(1,i+1);
    
    for j=1:N_axon
       
        [v(j,i+1),n(j,i+1), m(j,i+1),h(j,i+1),I_k(j,i+1),I_n(j,i+1),I_L(j,i+1)]=main_calculation_multi(n(j,i),m(j,i),h(j,i),dt,v(j,i),I_induct(j,i));
        
        I_induct(j+1,i)=I_k(j,i+1)+I_n(j,i+1)+I_L(j,i+1);
        
        
    end
    
   
    time_axis=[time_axis, time_axis(i)+dt];    
    
    figure (1)

for j=1:N_axon
    
    
    
    subplot(4,ceil(N_axon/4),j);
       plot(time_axis,v(j,:), '.b');
       legend('V_plot');
        if j==1
    title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
    end   
     
end
  drawnow;
figure (2)

for j=1:N_axon
    subplot(4,ceil(N_axon/4),j);
       plot(time_axis,I_induct(j,:),'.k',time_axis,I_n(j,:),'.b',time_axis,I_k(j,:),'.g');
       legend('Induct current','Na current','Ik current');
       %ylim([-0.01 0.02]);
       
           if j==1
    title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
    end
end
      drawnow;
    
%     
%     
% for j=1:N_axon
%     figure (j)
%        plot(time_axis,v(j,:), '.b');
%        legend('V_plot');
%        title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
%        drawnow;
% 
%     figure(j+N_axon)
%         plot(time_axis,I_induct(j,:),'.k',time_axis,I_n(j,:),'.b',time_axis,I_k(j,:),'.g');
%        legend('Induct current','Na current','Ik current');
%        %ylim([-0.01 0.02]);
%        title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
%        drawnow;
% 
%        
% end



end

 
end%axon_cell
