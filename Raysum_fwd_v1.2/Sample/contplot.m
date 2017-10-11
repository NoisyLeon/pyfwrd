% contplot.m: Colour plot of traces
% Usage: grd=contplot(T,baz,dt,minn,maxn,t1,t2,normflag)

function grd=contplot(T,baz,dt,minn,maxn,t1,t2,normflag)

tshift = 0;

maxn=min(maxn,size(T,2));
ntraces=t2-t1+1;

% Extract data subset
for tr=t1:t2
  avec=zeros(1,maxn);
  avec=T(tr,minn:maxn);

% Remove average
  avec=avec-mean(avec);

  if normflag == 0

    % Trace-normalize
    maxtr=max(abs(avec));
    if maxtr < 1E-10  % zero trace
      maxtr=1;
    end
    avec=avec/maxtr;

  end

  grd(tr,:)=avec;

end

if normflag == 1

  % Section-normalize
  grd=grd/(max(max(abs(grd))));

end


npts=size(grd,2);
ntr=size(grd,1);

% Because of the way shading works, need to add a dummy trace:
grd(ntr+1,:)=zeros(1,npts);
grd(:,npts+1)=zeros(ntr+1,1);

tvec=[(minn-1)*dt-dt/2:dt:(maxn-1)*dt+dt/2] + tshift;
dbaz=baz(2)-baz(1);
bazvec=[(baz(t1)-dbaz/2):dbaz:(baz(t2)+dbaz/2)];


% -- Colour stuff:
%        r1=[(0:31)/31,ones(1,32)];
%        g1=[(0:31)/31,(31:-1:0)/31];
%        b1=[ones(1,32),(31:-1:0)/31];
%        rwb=[r1',g1',b1'];
%colmin=0; colmax=1; colmid=0.65; colsteps=100;
%cold1=(colmid-colmin)/(colsteps-1); cold2=(colmax-colmid)/(colsteps-1);
%colvec=[colmin:cold1:(colmid-cold1),colmid:cold2:colmax]';
%rwb=[ colvec, colvec, colvec ];
%colormap(rwb);

colormap(gray);


pcolor(bazvec,tvec,grd');
grid off;
shading flat;
axis ij;
axis([baz(t1) baz(t2) (minn-1)*dt (maxn-1)*dt]);
caxis([-1 1]);





