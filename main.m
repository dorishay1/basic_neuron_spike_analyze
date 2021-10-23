clear; close all; clc;
 % Ex - 1


 %% Get the data

load ('S1.mat');
load ('S2.mat');

%% Properties
si = S1;

fs = 10000;      %sampling rate
dt = 1/fs;       %time step
N = length(si);  %number of components in si
t = dt*(1:N);    %time vector
TH = -30;        %threshold
t_cur_step = 0.2; %duration of current step in sec
t_seg = 0.3*fs;   %number of samples in one segment

%% (1)Plot

plot(t,si);
title('voltage of the neuron as function of time')
xlabel('time (sec)');
ylabel('voltage (mV)');
ylim([-90,60]);
xlim([0,7]);
hold

%% (2)Finding the spike times
SaTH = (si > TH);               %logical vector (1 for above the TH)
SaTHdiff = diff (SaTH);         %marks where there are changes in the previous vector

L2H = find(SaTHdiff == 1);      %marks low to high as 1
H2L = find(SaTHdiff == -1);     %marks high to low as -1

LM = zeros(2,length(H2L));      %setting a matrix for the data - 1st row - time in ms, 2nd row - the voltage

for i = 1:length(LM)
    range = L2H(i):H2L(i);      %set the range for the local maxima
    A = si(range);              %take si only in the range 
    [value, index]=max(A);      %extract value and index
    LM(1,i) = range(index);
    LM(2,i) = value;
end

%% (3)Finding the spike rate per segment

seg = (t_seg:t_seg:length(si));     %setting times for segments
SC = zeros(1,length(seg));          %preparing segment count vector



for j = 1:length(SC)
    if j==1
        Dlimit = 0;                             %first segment needs a down limit of 0 (from start)
    else
    Dlimit = seg(j-1);                          %all others using the previous upper limit
    end
    Ulimit = seg(j);                            %setting upper limit for segment
    B = (LM(1,:) > Dlimit & LM(1,:) < Ulimit);  %mark with ones the local maximas only within limits
    SC(j) = sum(B);                             %count number of local max
end 

R = SC/t_cur_step;                              %spikes per second per segment

%% (4) Adjust plot

%scatter L2H, H2L, LM & firring rate
LM_t = LM(1,:)*dt;                      
LL = si(L2H);
MM = si(H2L);

scatter(LM_t,LM(2,:));
scatter(L2H*dt,LL);
scatter(H2L*dt,MM);

a= (t_seg/2:t_seg:N);
R_text = num2str(R);
R_t = split(R_text);
text(a*dt,ones(1,length(SC))*50,R_t);
legend({'voltage', 'local maxima', 'low to high', 'hight to low'},'FontSize',7, 'Location', 'northeastoutside')
text(6.4,50,'\leftarrow firing rate (Hz)', 'FontSize', 7)

%setting figure size
x0=10;
y0=10;
width=1850;
height=400;
set(gcf,'position',[x0,y0,width,height])
