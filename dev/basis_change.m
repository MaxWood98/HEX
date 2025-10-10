%basis change
clearvars
clc

%basis a
ax = [1,0];
ay = [0,1];

%basis b (defined by vectors in a)
bx = [1,0.5];
by = [-0.5,1];

bx = bx/norm(bx);
by = by/norm(by);


vtest = [1,-2];
vtest = vtest/norm(vtest);



Mb2a = [bx(1),by(1);bx(2),by(2)]
Ma2b = inv(Mb2a)

m = Mb2a;
det = m(1,1)*m(2,2) - m(1,2)*m(2,1);
Ma2bt = (1.0/det)*[m(2,2),-m(1,2);-m(2,1),m(1,1)];


%plot
len = 1.0;
pref = [0,0];
cla reset
hold on
plot_vector(pref,ax,len,'r')
plot_vector(pref,ay,len,'r')

plot_vector(pref,bx,len,'b')
plot_vector(pref,by,len,'b')

plot_vector(pref,vtest,len,'g')


vtestb = Ma2b*vtest';
vtesta = Mb2a*vtestb;

vtestt = vtestb(1)*bx + vtestb(2)*by;

plot([pref(1),pref(1)+len*vtestt(1)],[pref(2),pref(2)+len*vtestt(2)],'m--')



hold off
axis equal



%% functions 

function plot_vector(pref,vector,len,style)
    plot([pref(1),pref(1)+len*vector(1)],[pref(2),pref(2)+len*vector(2)],style)
end