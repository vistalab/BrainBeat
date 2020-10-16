

t = [1:1000]/1000;

% s1 = log((10:100).^2);%10+sin(t*2*pi);
% s2 = log(10:100);%5+cos(t*2*pi);
s1 = 10+sin(t*2*pi);
s2 = 5+cos(t*2*pi);

te1 = 8;
te2 = 28;

r = (s2 - s1)/(te2 - te1);
s0 = s1 - (s2 - s1)*(te1/(te2 - te1));

figure

subplot(2,1,1),hold on
plot(s1,'ro')
s1_test = s0 + r*te1;
plot(s1_test,'k-','LineWidth',2)


subplot(2,1,2),hold on
plot(s2,'bo')
s2_test = s0 + r*te2;
plot(s2_test,'g-','Linewidth',2)


%% adding the log

t = [1:1000]/1000;

s1 = log(10+sin(t*2*pi));
s2 = log(5+cos(t*2*pi));

te1 = 8;
te2 = 28;

r = (s2 - s1)/(te2 - te1);
lns0 = s1 - (s2 - s1)*(te1/(te2 - te1));

t2s = -1./r;
s0 = exp(lns0);

figure

subplot(2,1,1),hold on
plot(exp(s1),'ro')
% s1_test = exp(lns0).*exp(r*te1);
s1_test = s0.*exp(-te1./t2s);
plot(s1_test,'k-','LineWidth',2)

subplot(2,1,2),hold on
plot(exp(s2),'bo')
s2_test = s0.*exp(-te2./t2s);
plot(s2_test,'g-','Linewidth',2)



