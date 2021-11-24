


% say 100Hz
srate = 100;
tt = [1:5000]/srate;

% say we have a 1 Hz response
a = sin(tt*2*pi*.8);

figure('Position',[0 0 800 300]),hold on
plot(tt,a)
% say our TR  = 1.5 s
sample_points = find(mod(tt,1.5)==0);
plot(tt(sample_points),a(sample_points),'r*')
plot(tt(sample_points),a(sample_points),'r-')
xlabel('time (s)')

