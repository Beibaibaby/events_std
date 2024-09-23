%{
This code calls a mex .cpp file to simulate a birth-death process
rate model with plasticity on the connections emerging from the exictatory
population.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basic system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syssize=325;
Ne = syssize;
Ni = syssize;
ttoss=0;
tlength=3e4+ttoss;
dt = 0.1;

%connection strengths
a = sqrt(1000/syssize);
jee  = 2*a;
jei  = 1*a;
jie  = 5*a;
jii  = 2*a;
%input drives
Ie   = -0.12*a;
Ii   = -0.2*a;
%inhibitory timescale
taui = 1;

%synapse parameters
trie = 40;
tdie = 10;
threshie=0.2;
magie=2;
slie=50;

tree = 40;
tdee = 10;
threshee=0.5;
magee = 2;
slee = 50;

%for ran2 inside the simulation
seed = -round(10000*rand);

%define parameter vectors for passing to our function
times = [tlength, ttoss, dt];
params = [jee, jei, jie, jii, Ie, Ii, taui];
ie_params = [trie,tdie,threshie,magie,slie];
ee_params = [tree,tdee,threshee,magee,slee];

%choose initial condition
N = [Ne, Ni];
re0 = 0.6;
ri0 = 0.9;
n0 = [round(re0*Ne), round(ri0*Ni)];

%run the simulation
[t,ve,vi,pE,pI] = markov_2D(times,N,n0,params,ee_params,ie_params,seed);

%scale time
t=t/100;

%plot example traces
clf();
subplot(2,1,1);
hold on;
plot(t,pE,'Color','red');
plot(t,pI,'Color','blue');
xlim([ttoss/100 tlength/100]);
ylim([-0.01 1.01]);
ylabel('Plasticity');
    
subplot(2,1,2);
hold on;
plot(t,ve,'Color','red');
plot(t,-vi,'Color','blue');
xlabel('Time (s)');
ylabel('Activity');





