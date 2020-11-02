%READ IN THE DATA
fid=fopen('trace_output.dat','r');
lt=fread(fid,1,'integer*8');
lrec=lt*7;   %length of the data records
data=fread(fid,lrec,'real*8');
data=reshape(data, [7,lt])';    %transpose!


%ORGANIZE
t=data(:,1);
x=data(:,2);
y=data(:,3);
z=data(:,4);
vx=data(:,5);
vy=data(:,6);
vz=data(:,7);


%MAKE SOME PLOTS
close all;

%COMET PLOTS WORK IN MATLAB

figure;
comet3(x,y,z);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');



figure;
MS=20;
FS=24;

plot3(x,y,z);
hold on;
plot3(x(1),y(1),z(1),'*','MarkerSize',MS,'MarkerFaceColor','white');
plot3(x(lt),y(lt),z(lt),'o','MarkerSize',MS,'MarkerFaceColor','white');
set(gca,'FontSize',FS);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
print -depsc2 path.eps    %postscript
%print -dpng path.png    %png

