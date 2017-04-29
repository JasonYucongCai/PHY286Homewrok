

%Axon_cell_par(0,0,1,0,1,2,2,0.01,50,1E-15,12);
function [] = Axon_cell_par(n_in, m_in, h_in,v_input,I_start_time, I_end_time,I_induct_value,dt,tmax,tolerance,N_axon)
tic;

%parpool('local',2);%--------------------------------------------------------------------------------------------------------------
parpool('local');%shutdonw the pool by delete(gcp);
time_axis = linspace(0, tmax, round(tmax/dt));

N=length(time_axis);


[v_axon,n_axon,m_axon,h_axon,I_k_axon,I_n_axon,I_L_axon]=initial_value_multi(n_in,m_in,h_in,dt,v_input,tolerance);
% function[v_initial,n_initial,m_initial,h_initial,I_ki,I_ni,I_Li]
%     =initial_value_multi(n_in,m_in,h_in,dt,v_input,tolerance)
 
ParaPool_size = parcluster('local');
P_N=ParaPool_size.NumWorkers;

parfor i=1:N_axon
for j=1:2%P_N    
v(i,j)=v_axon;
n(i,j)=n_axon;
m(i,j)=m_axon;
h(i,j)=h_axon;
I_k(i,j)=I_k_axon;
I_n(i,j)=I_n_axon;
I_L(i,j)=I_L_axon;
%for j=1:N_axon
I_induct(i,j)=0;
%end
end
end

time_axis=0;

if P_N~=0
   %for i=1:2%P_N
    time_axis=[time_axis, time_axis(1)+dt];
   %end
end

for j=1:N_axon
[v(j,2),n(j,2), m(j,2),h(j,2),I_k(j,2),I_n(j,2),I_L(j,2)]=main_calculation_multi(n(j,1),m(j,1),h(j,1),dt,v(j,1),I_induct(j,1));
        
I_induct(j+1,1)=I_k(j,2)+I_n(j,2)+I_L(j,2);
end



for i=2:N-1
    

    
    if (time_axis(i)>=I_start_time)&& (time_axis(i)<=I_end_time);
        I_induct_initial_compartment=I_induct_value;
    else
        I_induct_initial_compartment=0;
    end
    
    I_induct(1,i+1)= I_induct_initial_compartment;
    [v(1,i+1),n(1,i+1), m(1,i+1),h(1,i+1),I_k(1,i+1),I_n(1,i+1),I_L(1,i+1)]=main_calculation_multi(n(1,i),m(1,i),h(1,i),dt,v(1,i),I_induct_initial_compartment);
    
    I_induct(2,i)=I_k(1,i+1)+I_n(1,i+1)+I_L(1,i+1);
    
%        fixed_i=i;
%        ti1(j)=n(j,fixed_i-1);
%        ti2(j)=m(j,fixed_i-1);
%        ti3(j)=h(j,fixed_i-1);
%        ti4(j)=dt;
%        ti5(j)=v(j,fixed_i-1);
%        ti6(j)=I_induct(j,fixed_i-1);
%     
     parfor j=1:N_axon
         
       ti1(j)=n(j,i-1);
       ti2(j)=m(j,i-1);
       ti3(j)=h(j,i-1);
       ti4(j)=dt;
       ti5(j)=v(j,i-1);
       ti6(j)=I_induct(j,i-1);
    
     end
     
     parfor j=1:N_axon
   
        [to1(j),to2(j),to3(j),to4(j),to5(j),to6(j),to7(j)]=main_calculation_multi(ti1(j),ti2(j),ti3(j),ti4(j),ti5(j),ti6(j));

     end
     parfor j=1:N_axon
        v(j,i)=to1(j);
        n(j,i)=to2(j); 
        m(j,i)=to3(j);
        h(j,i)=to4(j);
        I_k(j,i)=to5(j);
        I_n(j,i)=to6(j);
        I_L(j,i)=to7(j);
        I_induct(j+1,i)=to5(j)+to6(j)+to7(j);
     end
       % I_induct(j+1,fixed_i)=I_k(j,fixed_i)+I_n(j,fixed_i)+I_L(j,fixed_i);
       
    
%     parfor j=1:N_axon %#ok<*PFUNK>
%        ti1=n(j,i-1);
%        ti2=m(j,i-1);
%        ti3=h(j,i-1);
%        ti4=dt;
%        ti5=v(j,i-1);
%        ti6=I_induct(j,i-1);
%         
%         
%         [to1,to2,to3,to4,to5,to6,to7]=main_calculation_multi(ti1,ti2,ti3,ti4,ti5,ti6);
%         v(j,i)=to1;
%         n(j,i)=to2; 
%         m(j,i)=to3;
%         h(j,i)=to4;
%         I_k(j,i)=to5;
%         I_n(j,i)=to6;
%         I_L(j,i)=to7;
%         
%         I_induct(j+1,i)=I_k(j,i)+I_n(j,i)+I_L(j,i);
%         
%         
%     end
    
    

    time_axis=[time_axis, time_axis(i)+dt];    

    
if mod(i,100)==2
    figure (3)

for j=1:N_axon

    subplot(4,ceil(N_axon/4),j);
       plot(time_axis,v(j,:), '.b');
       legend('V_plot');
        if j==1
    title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
    end   
     
end
  drawnow;
figure (5)

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
end
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

toc;
figure (1)
for j=1:N_axon

    subplot(4,ceil(N_axon/4),j);
       plot(time_axis,v(j,:), '.b');
       legend('V_plot');
        if j==1
    title(sprintf('Time step %d of %d at time %d s', i, N,(i*dt)));
    end   
     
end
  
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

delete(gcp);%------------------------------------------------------------------------------------------------------------------
end%axon_cell
