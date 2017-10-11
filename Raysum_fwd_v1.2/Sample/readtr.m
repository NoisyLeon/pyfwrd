clear all
format compact
format short g

filename='sample.tr';
fid=fopen(filename,'r');

% Read header:
line=fgetl(fid);
while line(1) == '#'
  line=fgetl(fid);
end
dum=str2num(line);
ntr=dum(1); nsamp=dum(2); dt=dum(3); align=dum(4); shift=dum(5);
[ntr,nsamp,dt,align,shift]

% Read each trace
line=fgetl(fid);
for itr=1:ntr
  itr
  while line(1) == '#'
    line=fgetl(fid);
  end
  for isamp=1:nsamp
    traces(:,isamp,itr)=str2num(line);
    line=fgetl(fid);
  end
end

tvec=-shift:dt:(dt*(nsamp-1)-shift);
size(traces)

fclose(fid)

subplot(3,1,1)
contplot(squeeze(traces(2,:,:))',[1:ntr]',dt,1,nsamp,1,ntr,1);
title('SV')
subplot(3,1,2)
contplot(squeeze(traces(3,:,:))',[1:ntr]',dt,1,nsamp,1,ntr,1);
title('SH');
subplot(3,1,3)
contplot(squeeze(traces(1,:,:))',[1:ntr]',dt,1,nsamp,1,ntr,1);
title('P');
xlabel('Trace number');
ylabel('Time (s)');
