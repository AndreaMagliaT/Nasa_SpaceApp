%%Simulation of an epidemic outbreak model in a closed space with the 
%%deployment of DAIANA protocol, developed for SpaceApps Challenge COVID-19 edition

clear all
close all 
clc

%% parameters of the model
N_pop = 1000;   %population
N_days = 300;   %days considered
r_meet0 = 10;   %base unilateral encounter rate
prob0 = 0.02;   %base infection probability

%tresholds considered for the protocol
tresh1 = N_pop * 0.01;	
tresh2 = N_pop * 0.10;
tresh3 = N_pop * 0.4;

%Diana activation
active = false;

%% recovery time cdf
    recovery_mean = 30; %average of death/recovery time
    recovery_std = 2;   %death/recovery standard deviation
    recovery_array= 0:1:N_days;
    recovery_cdf = cdf('norm',recovery_array,recovery_mean,recovery_std);

%% detection time cdf
    if(active == true)
        detection_mean = 15; %average detection time, it is lower with DAIANA
    else
        detection_mean = 20;
    end
    detection_std = 2;  %detection standard deviation
    detection_array = 0:1:N_days;
    detection_distribution = cdf('norm',detection_array,detection_mean,detection_std);

%% repeated simulation of the model in order to achieve mean behavioural indexes
N_sim = 100;
sim_cases = zeros(N_sim,N_days);    %total cases matrix
sim_infected = zeros(N_sim,N_days); %infected matrix
sim_detected = zeros(N_sim,N_days); %detected matrix

for k=1:N_sim
    
    %lockdown flags for simulating less flexible protocols
    lock2 = false;
    lock3 = false;
    
    %initialisation of the variables
    lock_count = 0;
    N_det = 0 ;
    
    health_state = zeros(1,N_pop);  %i-th element is 1 if i-th person is infected, it is 0 otherwise
    infected = zeros(1,N_pop);      %i-th element is 1 if i-th person is infected, it is 0 otherwise
    detected = zeros (1,N_pop);     %i-th element is 1 if i-th person is detected as infected, it is 0 otherwise
    recovered = zeros(1,N_pop);     %i-th element is 1 if i-th person is recovered (or dead), it is 0 otherwise
    
    total_cases = zeros(1,N_days);      %i-th element contains the number of total cases on day i 
    infected_counter = zeros(1,N_days); %i-th element contains the number of infected on day i
    detected_counter = zeros(1,N_days);  %i-th element contains the number of detected on day i
    
    detection_time = zeros (1,N_days);  %dynamical detection time
    r = zeros(1,N_days);                %dynamical unilateral encounter rate
    prob = zeros(1,N_days);             %dynamical infection probability
    phase = zeros(1,N_days);            %dynamical threat-level

    %first case of infection
    health_state(ceil(rand*1000)) = 1 ;
    
    %simulation of a model
    for d=1:1:N_days
      day_before = N_det ; 
      N_det = sum(detected) - sum(and(recovered,detected));
      difference = N_det-day_before ;
      
      %phase 0
      if(lock2 == false && lock3 == false)
          phase(d) = 0;    
          detection_time(d) = 14 ;
          r(d) = r_meet0;
          prob(d) = prob0;
      end

      if(active)    %DAIANA deployment if active is true
          %phase 1
          if((N_det > tresh1 && N_det <= tresh2 || (30>difference && difference>20)) && (lock3==false && lock2==false))
              phase(d) = 1;
              r(d) = r_meet0 - 0.3*r_meet0;
              prob(d) = prob0 - 0.2*prob0;
          end
          %phase 2
          if((N_det > tresh2 && N_det <= tresh3|| (30>difference && difference>40) || lock2 == true) && lock3==false)
              if(lock2 == true)
                  lock_count = lock_count + 1;
              end
              phase(d) = 2;
              detection_time(d) = 10; 
              r(d) = r_meet0 - 0.8*r_meet0;
              prob(d) = prob0 - 0.5*prob0;
              lock2 = true;
              if(lock_count >= 30)
                  lock2 = false;
                  lock_count = 0;
              end
          end
          %phase 3
          if(N_det >= tresh3 || difference>40 || lock3 == true)
            phase(d) = 3;
            r(d) = 0.1*r_meet0;
            prob(d) = 0.2*prob0;
            detection_time(d) = 10 ;
            lock3 = true;
            if(N_det == 0)
                lock3 = false;
            end
          end
      end

    clear numeri_casuali
        for i=1:1:N_pop

            if (health_state(i) == 1 )
               infected (i) = infected(i)+1;          
            end
            %%isolamenti
            if (infected(i) > 0)
               if(rand <= detection_distribution(infected(i)))
                   detected (i) = 1;
               end
            end
            if(infected(i) > 0)
                if(rand <= recovery_cdf(infected(i)))
                   recovered(i) = 1;
                end
            end
            encount_numb= ceil(r(d)*rand) ;
            rand_vect=ceil(N_pop*rand(1,encount_numb)) ;  

            for x=1:1:encount_numb
                if (health_state(i)==1 || health_state(rand_vect(x))==1) 
                    if (recovered(i) == 0)
                        if( detected(rand_vect(x))==1  && rand < (0.5*prob(d)+~active*0.5*prob(d)))  %put 0.5*prob(d) if diana is active
                            health_state(rand_vect(x))= 1 ;  
                            health_state(i)= 1 ;
                        end
                        if(detected(rand_vect(x))==0  && rand < (0.5*prob(d)+~active*0.5*prob(d)))  %put 0.5*prob(d) if diana is active
                             health_state(rand_vect(x))= 1 ;  
                            health_state(i)= 1 ;
                        end
                    end 
                end
             end  
        end
        infected_counter(d)=sum(health_state)-sum(recovered) ; 
        total_cases(d)= sum(health_state) ; 
        detected_counter(d)=N_det ;
    end 

    sim_cases(k,:) = total_cases;
    sim_infected(k,:) = infected_counter;
    sim_detected(k,:) = detected_counter;
end

% figure (1)
% hold on 
% plot (infected_counter,'black')
% xlabel ('Days')
% ylabel('cases')
% % ylim([0 1000]);
% grid on
% plot(total_cases,'r+')
% plot(detected_counter,'b.')
% legend ('infected cases','total cases','known infected cases') ;
% 
% figure(2)
% plot(phase);
% xlabel ('Days')
% ylabel('DAIANA threat level')
% grid on

%mean behaviour
mean_total_cases = mean(sim_cases);
max_sim = max(sim_infected');
mean_infected = mean(sim_infected);
std_infected = std(sim_infected);
mean_detected = mean(sim_detected);
std_detected = std(sim_detected);

figure(3)
plot(mean_total_cases);
xlabel('days');
ylabel('mean cases');

figure(4)
plot(max_sim');
xlabel('days');
ylabel('max infected');

if(active == true)
    figure(5)
    plot(mean_infected,'b');
    hold on 
    xax=1:1:length(mean_infected) ;
    xlabel('days');
    ylabel('mean infected');
    title('Mean infected curve with DAIANA protocol');
    hold on
    plot(mean_detected,'r');
    errorbar(xax(1:10:end),mean_infected(1:10:end),2*std_infected(1:10:end),'b.')
    legend('mean infected','mean detected','confidence interval at 97.5%');
    grid
else
    figure(5)
    plot(mean_infected,'b');
    hold on 
    xax=1:1:length(mean_infected) ;
    xlabel('days');
    ylabel('mean infected');
    title('Mean infected curve with no protocol');
    hold on
    plot(mean_detected,'r');
    errorbar(xax(1:10:end),mean_infected(1:10:end),2*std_infected(1:10:end),'b.')
    legend('mean infected','mean detected','confidence interval at 97.5%');
    grid
end


