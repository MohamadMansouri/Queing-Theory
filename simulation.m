%% COPYRIGHTS && PLEASE READ ME
% This code has copyrights to Mohamad MANSOURI
% Please dont distribute it
% THIS CODE WILL SIMULATE 3 queing policies...
% The first is Shortest Job First (SJF)
% The second is First Come First Serve (FCFS)
% The last is a Custom policy that i would like to call it Donation


%% EXPLANATION OF THE DONATION POLICY %%
% it is a combination between the shortest job first policy and the
% processpr sharing policy... it works as follows
% we assume we have 2 types of jobs BIG and SMALL ones
% each time we get a big job we do as SJF do, we through it to the back of
% all the small jobs but this time we will let each small job donate this
% big job an amount X of time, This mean the server wills serve the big job
% for an amount X then serve the small job...
% this happens repeatidly so as the small jobs are served the big job is
% getting smaller and smaller...
% the amount of donation X is an hyperparameter that i tuned to find an
% optimal value, it appear to be 2.



%% CLEAR THE DIRT %%

clc 
clear
tic
%% CONSTANTS %%
r=1; 
x_axis=0.02:0.02:0.18; % lamda values from 0.02 to 0.18 taking 0.02 steps this mean raw goes from 0.1 --> 0.9
NSamples=10000; % for each simulation take "NSample" arrivals

seed1 = RandStream.create('mcg16807','Seed',3); % make 2 different
seed2 = RandStream.create('mcg16807','Seed',4); % seeds

U1 = rand(seed1, 1, NSamples); % uniform distributed numbers between 0 and 1 (will be used to form the interarrival time vector)
U2 = rand(seed2, 1, NSamples); % uniform distributed numbers between 0 and 1 (will be used to form the service time vector )

fprintf('THIS PROGRAM WILL DISPLAY 3 FIGURES DESCRIBING THE PERFORMANCE OF FCFS, SJF AND A CUSTOM QUEING SYSTEMS... SO BE READY\n')

%% ITERATE %% 
for j = x_axis % loop on all the 9 values of lambda 
            
    fprintf('calculation are going for rho = %f...\n',j)
    lambda=j; % lambda current value 
    donate=0; % these is only an initialization u will see how this variable used later
    donated=0;
    %%----CACULATE S----%%
    %%%%%%%%%%%%%%%%%%%%%%
    for k=1:NSamples;if U2(k)<=0.02 S(k)=201;else S(k)=1;end;end  % S is the service time, let it equal 1 sec with probability 0.98 and equal 201 sec with probability 0.02
                                                                  % E(S) = 5 , E(S^2)= 809
    SSJF=S;
    SC=S;
    %%%%%%%%%%%%%%%%%%%%%%	
    %%------------------%%


    %%---CACULATE tau---%%
    %%%%%%%%%%%%%%%%%%%%%%

    tau = -1/lambda*log(1-U1); % tau is inter-arrival time, we transformed the uniform random numbers to exponential ( 1 - exp(-lambda*x)=U1) 

    %%%%%%%%%%%%%%%%%%%%%%
    %%------------------%%    

	%%---FCFS INITIALIZATIONS---%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	T = cumsum(tau); % T(i) is the time since the system has started till the arrival of i event
    D = zeros(NSamples, 1); % D(i) is the departure time of the ith event
    D(1) = S(1) + tau(1); % initialize D; D(1) is the departural time of the first event ( no queuing time)
    W = zeros(NSamples, 1); % W(i) is the waiting time of the ith event
    W(1) = 0; % W(1) = 0 since the first event will not wait since the queue is empty    

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%--------------------------%%



	%%----SJF INITIALIZATIONS---%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	TSJF = cumsum(tau); % TJSF(i) is the time since the system has started till the arrival of i event ( SJF VECTOR )
	DSJF = zeros(NSamples, 1); % DSJF(i) is the departure time of the ith event ( SJF VECTOR )
	DSJF(1) = S(1) + tau(1); % initialize DSJF; DSJF(1) is the departural time of the first event ( no queuing time) ( SJF VECTOR )	
	WSJF = zeros(NSamples, 1); % W(i) is the waiting time of the ith event ( SJF VECTOR )
	WSJF(1) = 0; % WSJF(1) = 0 since the first event will not wait since the queue is empty ( SJF VECTOR )     

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%--------------------------%%
    

	%%---CUSTOM POLICY INITIALIZATIONS---%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	TC = cumsum(tau); % T(i) is the time since the system has started till the arrival of i event
    DC = zeros(NSamples, 1); % D(i) is the departure time of the ith event
    DC(1) = SC(1) + tau(1); % initialize D; D(1) is the departural time of the first event ( no queuing time)
    WC = zeros(NSamples, 1); % W(i) is the waiting time of the ith event
    WC(1) = 0; % W(1) = 0 since the first event will not wait since the queue is empty    

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%--------------------------%%
    
    
    jobs_in_SJF=0;
    jobs_in_FCFS=0;
    
	for i = 2:NSamples % iterate on the samples  
       %% SJF CALCUALTION PART %% 
        k=1; 
        while k<=NSamples-i && TSJF(i+k)< DSJF(i-1)%check if thier are jobs in the SJF que 
            if SSJF(i+k)<SSJF(i-1+k)% check if thier is jobs in queue smaller than the job that will be server now
                % replace the current job with the smaller job
                TSJF([i-1+k i+k])= TSJF([i+k i-1+k]); % swap the indices of time of arrival vector
                SSJF([i-1+k i+k])=SSJF([i+k i-1+k]); % and swap the indices of size
            end
            k=k+1;
        end
        
        if TSJF(i)> DSJF(i-1) k=0; end % check if the system is idle

        jobs_in_SJF=[jobs_in_SJF k]; % save the number of jobs in the system

        
        
        DSJF(i) = max(TSJF(i), DSJF(i-1)) + SSJF(i); % departure time will be equal to the
									 	 % Service time 
                                         %  + 
                                         % Deparure time of the previous event (if the event is arrived BEFORE previous one is finished)
                                         % Arrival time of this event (if the event is arrived AFTER previous one is finished (i.e. empty queue))
        WSJF(i) = max(DSJF(i-1)-TSJF(i), 0);  % wait is either 0 (i.e. empty queue) or Departure time of the prev event - arrival of the new one if this value is +ve
    	
    
    
        %% FCFS CALCULATION PART %%
        k=0;
        while k<=NSamples-i && T(i+k)< D(i-1) % count the number of jobs in the FCFS system
            
            k=k+1;
        end
        jobs_in_FCFS=[jobs_in_FCFS k];
            
        D(i) = max(T(i), D(i-1)) + S(i); % departure time will be equal to the
									 	 % Service time 
                                         %  + 
                                         % Deparure time of the previous event (if the event is arrived BEFORE previous one is finished)
                                         % Arrival time of this event (if the event is arrived AFTER previous one is finished (i.e. empty queue))
        W(i) = max(D(i-1)-T(i), 0);  % wait is either 0 (i.e. empty queue) or Departure time of the prev event - arrival of the new one if this value is +ve
    	
	
       %% DONATION POLICY CALCULATION PART %%
        
        k=1; 
        
        while k<=NSamples-i && TC(i+k)< DC(i-1)%check if thier are jobs in the SJF que 
            if SC(i+k)<SC(i-1+k) && SC(i-1+k)>donate% check if thier is jobs in queue smaller than the job that will be server now and that our Big value is grater than the donation that he will take
                donate=2; % set donation value to 2 ( this is found as optimal after tuning)
                SC(i-1+k)=SC(i-1+k)-donate; % give the big value 2 sec of service as donation 
                
                donated=k; % count how many jobs are donating the big job
                % replace the current job with the small one
                TC([i-1+k i+k])= TC([i+k i-1+k]); % swap the indices of time of arrival vector
                SC([i-1+k i+k])=SC([i+k i-1+k]); % and swap the indices of size
                
            end
            k=k+1;
        end
        if donated ==0 donate=0;end % check if thier all jobs that donated this big job have finished
        
        DC(i) = max(TC(i), DC(i-1)) + SC(i)+ donate; % departure time will be equal to the
									 	 % Service time 
                                         %  + 
                                         % Deparure time of the previous event (if the event is arrived BEFORE previous one is finished)
                                         % Arrival time of this event (if the event is arrived AFTER previous one is finished (i.e. empty queue))
                                         %  +
                                         % The donated time that it gave to
                                         % the big job (if i didnt donate the donate value will be 0) 
        donated = donated -1; % decremant the number of jobs that donated the big job till 0 
        if donated<0 donated=0;end
        WC(i) = max(DC(i-1)-TC(i), 0);  % wait is either 0 (i.e. empty queue) or Departure time of the prev event - arrival of the new one if this value is +ve
     	
          
        
    end
  
    Delay_Time = zeros(NSamples, 1); % the delay time vector – in–queue wait time and service time
    Delay_Time_SJF = zeros(NSamples, 1); % the delay time vector – in–queue wait time and service time
    Delay_Time_C = zeros(NSamples, 1); % the delay time vector – in–queue wait time and service time
    for i = 1:NSamples,
        Delay_Time(i) = W(i) + S(i); % Delay time is adding the wait time and the service time of each event
        Delay_Time_SJF(i) = WSJF(i) + SSJF(i); % Delay time is adding the wait time and the service time of each event
        Delay_Time_C(i) = WC(i) + SC(i); % Delay time is adding the wait time and the service time of each event
						    
	end
    
    delay(r)=mean(Delay_Time); % calcule the mean delay time and append it to a list for each value of lambda
    delay_SJF(r)=mean(Delay_Time_SJF);% calcule the mean delay time and append it to a list for each value of lambda
    delay_C(r)=mean(Delay_Time_C);% calcule the mean delay time and append it to a list for each value of lambda
    r=r+1;
    
    
% draw graph of nb of jobs over iteration
figure(1)
subplot(3,3,j*5*10)
hold on;
title(j*5)
plot(1:size(jobs_in_SJF,2),jobs_in_SJF,'b')
plot(1:size(jobs_in_FCFS,2),jobs_in_FCFS,'g--')
xlabel('Iteration')
ylabel('# of jobs')
legend('SJF','FCFS')
hold off;
a = axes;
t1 = title('# of jobs in system after each iteration');
a.Visible = 'off';
t1.Visible = 'on';


%draw histograms of nb of jobs
figure(2)
subplot(3,3,j*5*10)
hold on;

nbins = 0:j*500;
title(j*5)
histogram(jobs_in_SJF,nbins ,'FaceColor','b','facealpha',1,'EdgeColor','k')
histogram(jobs_in_FCFS,nbins,'FaceColor','g','facealpha',0.7)
legend('SJF','FCFS')
hold off
a = axes;
t1 = title('jobs in system frequency');
a.Visible = 'off';
t1.Visible = 'on';

end

% draw the performone while increasing utilization
figure(3);
title('Delay Time')


mg1=80.9.*x_axis./(1/5-x_axis)+1/5; % this is the theoretical value of response time using E(S)=5 and E(S^2)=809 NB: C^2 = 31
mm1=1./(1/5-x_axis);
hold on;
plot(x_axis.*5,delay,'g') % plot the delay time for each lambda
plot(x_axis.*5,mg1,'r') % plot the theoretical curve 
plot(x_axis.*5,delay_SJF,'b') % plot the delay time for each lambda
plot(x_axis.*5,delay_C) % plot the delay time for each lambda
plot(x_axis.*5,mm1)
xlabel('\rho')
ylabel({'Delay Time','(sec)'})

hold off;
legend('FCFS','theoretic FCFS( mg1 )','SJF','CUSTOM( Donation )','theoretic PS( mm1 )')

fprintf('\n\n')
fprintf('we can see from the graph of the performance\n')
fprintf('that the FCFS system and the SJF system have roughly same performance at law utilization\n')
fprintf('but it appears that the SJF system does better as the utilization approaches 1\n')
fprintf('also the other graphs shows the same result in the 2 systems\n')
fprintf('moreover the donation policy is doing best performance at high load since it\nis close to sharing the processor which make it close to M/M/1 \n')
fprintf('\n\n')


toc
    