function plot_mouse

% Copyright 2007 John L Baker. All rights reserved.
% This software is provided AS IS under the terms of the Open Source
% MIT License. See http://www.opensource.org/licenses/mit-license.php.

% Plot results from a simulated mouse in a maze experiments

global DATADIR SRCDIR CELLSZ MAZETYPE MINTIME DM NCOMP TSTRT TEND 
global PLOTGRAY PLOTCONTOUR

% Identify which target cell model is in use (see also load_dm)
% For largely historical reasons, the number of compartments is
% used to identify which morphology data should be used.
DM=[]; %Clear dendrite morphology table (loaded on demand)
NCOMP=368; % L56a 50 micron morphology

% Locate the top level directory for data and source files
drv='C:\BNSF_1_0\'; % top level directory for source and data
SRCDIR=[drv 'Src\Baker\']; % used in DM load (.csv morphology files)

% Pick the directory where data results files are found
DATADIR=[drv 'Data\Default\']; % misc data results

CELLSZ=5; % size of spatial grid cell in cm 
MINTIME=0.5; % minimum time in cell in sec to measure firing rate
MAZETYPE=1; % 0-circle maze, 1=square maze
PLOTGRAY=0; % Use gray scale for colormaps (selected plots only)
PLOTCONTOUR=0; % Plot firing maps using filled contours

TSTRT=0e3; % Starting time for reports (in ms)
TEND=9999e3; % Ending time for reports (if any)

% -------------------------------------------------------------
% standard results plots

plot_swim_path
plot_firing_rates(0,0,'')
plot_spike_train(0,0,'')

plot_weight_summary('EC','PP','')
plot_weight_summary('CA3','AC','')

if MAZETYPE==0 % circle maze
   plot_theta_phase(0,[-25 25],[1 11],-1,[0 180;180 360;360 540;540 720;720 900],'')
   plot_linear_rates(0,[-25 25],[1 11],-1,[0 180;180 360;360 540;540 720;720 900],'')
   plot_weight_by_dist('EC','PP','',-1,-12,6,0.01)
   plot_weight_by_dist('CA3','AC','',-1,-12,6,0.1)
end

if MAZETYPE==1 % square maze
   plot_theta_phase(0,[-25 25],[-20 -15],1,[0 180;180 360;360 540;540 720;720 900],'')
   plot_linear_rates(0,[-25 25],[-20 -15],1,[0 180;180 360;360 540;540 720;720 900],'')
   plot_weight_by_dist('EC','PP','',-1,0,-17.5,0.01)
   plot_weight_by_dist('CA3','AC','',-1,0,-17.5,0.05)
end

if 0 % Optionally plot afferent drive (slow)
   plot_weighted_rates('EC','PP',-1)
   plot_weighted_rates('CA3','AC',-1)
end


% -------------------------------------------------------------
% specialized ad hoc plots -- customize as necessary

%plot_dendrite_spikes([-180 180])

%plot_theta_phase(10070,[-25 25],[-20 -15],1,[0 540],'subset')
%plot_theta_phase(20021,[-25 25],[-20 -15],1,[0 540],'subset')
%plot_theta_phase(32546,[-25 25],[-20 -15],1,[0 540],'CA3')
%plot_theta_phase(51001,[-25 25],[-20 -15],1,[0 540],'IN')

%plot_linear_rates(10070,[-25 25],[-20 -15],1,[0 540],'subset')
%plot_linear_rates(20021,[-25 25],[-20 -15],1,[0 540],'subset')
%plot_linear_rates(32546,[-25 25],[-20 -15],1,[0 540],'CA3')
%plot_linear_rates(51001,[-25 25],[-20 -15],1,[0 540],'IN')

%plot_weight_by_comp('EC','PP',-1)
%plot_weight_by_comp('CA3','AC',-1)

%plot_afferent_weight(10070,'PP')
%plot_afferent_weight(32546,'AC')

%plot_firing_rates(10001,10050,'subset') % EC
%plot_firing_rates(20001,20043,'subset') % DG
%plot_firing_rates(30001,30100,'subset') % CA3

%plot_layer_firing_rates('EC')
%plot_layer_firing_rates('DG')
%plot_layer_firing_rates('CA3')

%plot_spike_train(10001,19999,'EC')
%plot_spike_train(20001,29999,'DG')
%plot_spike_train(30001,49999,'CA3')

%plot_spike_train(50001,50999,'IN') % axo-axonic
%plot_spike_train(51001,51999,'IN') % somatic
%plot_spike_train(52001,52999,'IN') % bistatified
%plot_spike_train(53001,53999,'IN') % OL-M

%plot_field_centers('EC','',10001,19999,1,0,-17.5,0)
%plot_field_centers('DG','',20001,29999,1,0,-17.5,0)
%plot_field_centers('CA3','',30001,49999,1,0,-17.5,0)

% -------------------------------------------------------------

function plot_swim_path
global DATADIR CELLSZ MAZETYPE MINTIME TSTRT TEND PLOTGRAY

S=dlmread([DATADIR 'test-baker-mouse-states.txt'],',',1,1);
S=S(find(S(:,1)>=TSTRT & S(:,1)<=TEND),:);

cellSize = CELLSZ;

T=S(:,1);
X=S(:,2);
Y=S(:,3);
hm=(T(end)-T(1))/(length(T)-1)/1000;

last_time_read = T(end)

% Set up for the appropriate maze - circle or rectangle
if MAZETYPE==1
   xsize=25;
   ysize=25;   
   Wx = [-xsize:xsize xsize*ones(1,1+2*ysize) xsize:-1:-xsize -xsize*ones(1,1+2*ysize)];
   Wy = [-ysize*ones(1,1+2*xsize) -ysize:ysize ysize*ones(1,1+2*xsize) ysize:-1:-ysize];
else
   theta = 2*pi*(0:.01:1);
   xsize=25;
   ysize=25;
   Wx=25*cos(theta);
   Wy=25*sin(theta);
end
radius = 25;

Qt=zeros(2,2);
for k=1:length(T)
   q1=(X(k)>0)+1;
   q2=(Y(k)<0)+1;
   Qt(q2,q1)=Qt(q2,q1)+1;
end

times_in_quandrant = Qt

figure
ls1='--';
ls2='-';
if PLOTGRAY
   ls1='k--';
   ls2='k-';
end
plot(Wx,Wy,ls1,X,Y,ls2)
title('Swim path')
xlabel('X (cm)')
ylabel('Y (cm)')
axis([-radius-cellSize/2 radius+cellSize/2 -radius-cellSize/2 radius+cellSize/2]);
axis square

n=1+ceil(2*radius/cellSize);
NC=zeros(n,n);

for k=1:length(T)
   nx = ceil((radius+X(k))/cellSize);
   ny = ceil((radius+Y(k))/cellSize);
   NC(ny,nx)=NC(ny,nx)+1;
end

NCx=-radius-cellSize+cellSize*(1:size(NC,2));
NCy=-radius-cellSize+cellSize*(1:size(NC,1));

if 1
   figure
   plot_color(NCx,NCy,hm*NC,radius,16)
   title(['Cell occupancy time (sec)'])
end

% ----------------------------------------------------------

function plot_field_centers(layerId,setId,fromId,toId,activeOnly,x,y,dist)
global DATADIR MAZETYPE

% Read place field centers
% cols id,x,y,peakrate,isactive
PC = dlmread([DATADIR 'test-baker-' layerId '-pcloc' setId '.txt'],',',1,0);
idx = find(PC(:,1)>=fromId & PC(:,1)<=toId);
ID=PC(idx,1);
X=PC(idx,2);
Y=PC(idx,3);

if activeOnly
   idx=find(PC(idx,5)==1);
   ID=ID(idx);
   X=X(idx);
   Y=Y(idx);
end

if dist>0
   % find points near x and y
   idx=find(dist^2>=(X-x).^2+(Y-y).^2);
   disp(['Place fields near x=' num2str(x) ' y=' num2str(y)]);
   id_pcx_pcy=[ID(idx) X(idx) Y(idx)]
end

% Set up for the appropriate maze - circle or rectangle
if MAZETYPE==1
   xsize=25;
   ysize=25;   
   Wx = [-xsize:xsize xsize*ones(1,1+2*ysize) xsize:-1:-xsize -xsize*ones(1,1+2*ysize)];
   Wy = [-ysize*ones(1,1+2*xsize) -ysize:ysize ysize*ones(1,1+2*xsize) ysize:-1:-ysize];
else
   theta = 2*pi*(0:.01:1);
   xsize=25;
   ysize=25;
   Wx=25*cos(theta);
   Wy=25*sin(theta);
end
radius = 25;

% Do the plot
figure
plot(Wx,Wy,'--',X,Y,'o')
title([layerId ' Place Field Centers'])
axis([-radius-5 radius+5 -radius-5 radius+5]);
axis square

% ----------------------------------------------------------

function plot_firing_rates(fromId,toId,fileId)
global DATADIR CELLSZ MINTIME TSTRT TEND MAZETYPE
MStates=dlmread([DATADIR 'test-baker-mouse-states.txt'],',',1,1);
last_time_read = MStates(end,1)

cellSize = CELLSZ;
plotPath = 0 & fromId==toId                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         & fromId==toId;
plotSpikes = 0 & fromId==toId;
markTime = 0;

idx=find(MStates(:,1)>=TSTRT & MStates(:,1)<=TEND);
MStates=MStates(idx,:);
if isempty(MStates)
   disp('Time range not found')
   return;
end

T=MStates(:,1);
X=MStates(:,2);
Y=MStates(:,3);

hm=(T(end)-T(1))/(length(T)-1)/1000;

radius = 25;

if isempty(fileId)
   fn=[DATADIR 'test-baker-spikes.txt'];
else   
   fn=[DATADIR 'test-baker-' fileId '-spikes.txt'];
end

Spikes=dlmread(fn,',',1,0);
if isempty(Spikes)
   disp('No spikes found')
   return
end

idx = find(Spikes(:,1)>=fromId & Spikes(:,1)<=toId ...
   & Spikes(:,2)>=TSTRT & Spikes(:,2)<TEND);
Spikes=Spikes(idx,:);
if isempty(Spikes)
   disp('No spikes found in time range')
   return
end

Ctr=[];
FRpeak=[];

totalSpikes=0;
n=1+ceil(2*radius/cellSize);
NSall=zeros(n,n);
NCall=1e-8*ones(n,n);

rptNbr=0;
for rptId=fromId:toId
   
   n=1+ceil(2*radius/cellSize);
   NS=zeros(n,n);
   NC=1e-8*ones(n,n);
   
   Xsp=[];
   Ysp=[];
   ns=0;
   swimIdx=1;
   
   for k=1:size(Spikes,1)
      
      id=Spikes(k,1);
      t=Spikes(k,2);
      
      if id==rptId  
         if t>T(end)
            break
         end
         
         while t>T(swimIdx)
            x=X(swimIdx);
            y=Y(swimIdx);         
            swimIdx =swimIdx+1;
            nx = ceil((radius+x)/cellSize);
            ny = ceil((radius+y)/cellSize);
            NC(ny,nx)=NC(ny,nx)+1;            
            NCall(ny,nx)=NCall(ny,nx)+1;            
         end
         
         ns=ns+1;
         Xsp(ns)=x;
         Ysp(ns)=y;
         NS(ny,nx)=NS(ny,nx)+1;
         NSall(ny,nx)=NSall(ny,nx)+1;
      end
   end
   
   [sbits,skinfo,skspar]=infoMetrics(NS,NC);
   
   FR=(NC>=MINTIME/hm).*NS./(NC*hm);
   FRx=-radius-cellSize+cellSize*(1:size(FR,2));
   FRy=-radius-cellSize+cellSize*(1:size(FR,1));
   
   totalSpikes=totalSpikes+ns;
   if ns>0
      avg_spike_rate = sum(sum(NS))/((T(end)-T(1))/1000);
      peak_spike_rate = max(max(FR));
      FRpeak(rptId-fromId+1)=peak_spike_rate;
      FRmean(rptId-fromId+1)=avg_spike_rate;
   else
      avg_spike_rate = 0;
      FRpeak(rptId-fromId+1)=0;
      FRmean(rptId-fromId+1)=0;
   end
   
   if ns>0
      Ctr(rptId-fromId+1,:)=[mean(Xsp) mean(Ysp)];
      Rcov=cov([Xsp' Ysp']);
   else
      Ctr(rptId-fromId+1,:)=[NaN NaN];
      Rcov=zeros(2);
   end
   pcradius=abs(det(Rcov))^.25;
   
   if max(max(FR))>0
      disp(['id=' int2str(rptId)])
      avg_peak_rates_ctr = [avg_spike_rate peak_spike_rate Ctr(rptId-fromId+1,:) ]
      pcradius_sbits_skin_spar = [pcradius sbits skinfo skspar]
   end
   
   maxz=ceil(1e-4+max(max(FR)));
   plotDone=0;
   
   if 1 & ~plotPath & max(max(FR))>=2
      nc=min(2,ceil(sqrt(toId-fromId+1)));
      rptNbr=rptNbr+1;
      if rptNbr>nc*nc
         rptNbr=1;
      end
      if rptNbr==1
         figure
      end
      if fromId~=toId
         subplot(nc,nc,rptNbr)
      end
      plot_color(FRx,FRy,FR,radius,16)
      title(['Cell ' num2str(rptId) ' firing by location (Hz)'])
      plotDone=1;
   end
   
   if plotPath
      figure
      maxz=max(max(FR));
      plot_color(FRx,FRy,FR,radius,16)
      title(['Firing rate for cell ' num2str(toId) ' (Hz)'])
      hold on
      % plot firing locations
      plot3(X,Y,maxz*ones(size(X)),'w') 
      xc=X(1:20:end);
      yc=Y(1:20:end);
      % plot end of the path
      plot3(X(end),Y(end),maxz,'wo')
      if markTime
         % plot time marks
         plot3(xc,yc,maxz*ones(size(xc)),'w.')
      end
      plotDone=1;
   end
   
   if plotSpikes &~isempty(Xsp) & plotDone
      % plot spikes by location
      hold on
      plot3(Xsp,Ysp,maxz*ones(size(Xsp)),'w^')
   end
end

if toId>fromId
   average_firing_rate = totalSpikes/(T(end)-T(1))*1000/(toId-fromId+1)
   mean_peak_rate=mean(FRpeak)
   median_peak_rate=median(FRpeak)
   max_peak_rate = max(FRpeak)
   min_peak_rate = min(FRpeak)
end

if 1 & toId>fromId
   FRall=(NCall>=MINTIME/hm).*NSall./(NCall*hm);
   FRx=-radius-cellSize+cellSize*(1:size(FR,2));
   FRy=-radius-cellSize+cellSize*(1:size(FR,1));
   maxz=ceil(1e-4+max(max(FRall)));
   figure
   plot_color(FRx,FRy,FRall,radius,16)
   title(['Cells ' num2str(fromId) ' to ' num2str(toId) ' firing rate per cell (Hz)'])
end

if 1 & toId>fromId
   actfrp=FRpeak(find(FRpeak>=2));
   actfrm=FRmean(find(FRpeak>=2));
   if ~isempty(actfrp)
      figure
      plot((1:length(actfrp))/(1+length(actfrp)),sort(actfrp),'.-')
      title('Peak firing rates')
      ylabel('firing rate (Hz)')
      xlabel('quantile')
      mean_active_rate = mean(actfrm)
      mean_active_peak = mean(actfrp)
      axis([0 1 0 ceil(1.1*max(actfrp))])
   end
end

% -------------------------------------------------------------------

function plot_theta_phase(cellId,xlimits,ylimits,dir,tlim,fileId)
global DATADIR CELLSZ MINTIME TSTRT TEND MAZETYPE

% Passed params are:
% cellId = id of the cell to plot
% xlimits = array [xlow xhigh] of region to plot
% ylimits = array [ylow yhigh] of region to plot
% dir = +1 for left to right, -1 for right to left
% tlim = matrix where each row is a from and to time interval
% fileId = id of file to read for spikes (see code below)

% Hard coded params are:
xsz = 5; % x bin size in cm 
tpsz = 90; % theta phase bin size in degrees
pfx = 0; % place field center x coordinate
hxlim = dir; % heading limit x coord (-1=right to left, 1=left to right)
hylim = 0; % heading limit y coord
hdlim = 0; % minimum of dot product of heading and [hxlim hylim]
useCircularFit = 1; % Use circular fit (1) or -90-deg wrap point (0)
useNegSlope = 2; % 0= use best fit slope, 1= <0 always, 2= <0 if tie
textY=-260; % Y-coord location of intra-plot labels

% Statistical parameters
tStatP01 = 2.576; % two-sided t-stat value for P<.01 df=inf (normal)
tStatP05 = 1.960; % two-sided t-stat value for P<.05 df=inf (normal)

MStates=dlmread([DATADIR 'test-baker-mouse-states.txt'],',',1,1);
T=MStates(:,1);
X=MStates(:,2);
Y=MStates(:,3);
Hx=MStates(:,4);
Hy=MStates(:,5);

if 0
   % Note theta phase when crossing place field center as a test 
   % for spurious theta phase correlations based on regularity of the path.
   idx=find(Y>=ylimits(1) & Y<=ylimits(2) & X>=xlimits(1) & X<=xlimits(2));
   tpx0=mod(T(idx(find(X(idx(1:end-1))<=pfx & X(idx(2:end))>pfx)))',125)/125*360;
   figure
   xhc=histc(tpx0,0:72:360); % 72 deg = 25 ms
   bar(36:72:360,xhc(1:end-1),'k')
   title('Theta phase at place field center crossing')
   xlabel('Theta phase (deg)')
   ylabel('Count')
   
   % optionally, stop now
   %return
end

hm=(T(end)-T(1))/(length(T)-1)/1000;

if isempty(fileId)
   fn=[DATADIR 'test-baker-spikes.txt'];
else   
   fn=[DATADIR 'test-baker-' fileId '-spikes.txt'];
end

Spikes=dlmread(fn,',',1,0);
spikes_read=size(Spikes)
Spikes=Spikes(find(Spikes(:,1)==cellId),:);
if isempty(Spikes)
   disp('No spikes found')
   return
end

xlnum=ceil((xlimits(2)-xlimits(1))/xsz);
ntlim=size(tlim,1);

if ntlim>1
   figure
end

for n=1:ntlim
   tstrt=tlim(n,1);
   tend=tlim(n,2);
   time_range=[tstrt tend]
   
   Tspk=Spikes(find(Spikes(:,2)>=tstrt*1000 & Spikes(:,2)<=tend*1000),2);
   Xt=zeros(length(Tspk),1);
   Yt=zeros(length(Tspk),1);
   Hxt=zeros(length(Tspk),1);
   Hyt=zeros(length(Tspk),1);
   Theta=zeros(length(Tspk),1);
   swimidx=1;
   for k=1:length(Tspk)
      while swimidx<length(X) & Tspk(k)>T(swimidx)
         swimidx=swimidx+1;
      end
      Xt(k)=X(swimidx);
      Yt(k)=Y(swimidx);
      Hxt(k)=Hx(swimidx);
      Hyt(k)=Hy(swimidx);
      Theta(k)=360*(mod(Tspk(k)+125/2,125))/125-180;
   end
   
   % Pare down data to fit in x and y limits
   idx=find(Yt>=ylimits(1) & Yt<=ylimits(2) & ...
      Xt>=xlimits(1) & Xt<=xlimits(2) & ...
      Hxt*hxlim+Hyt*hylim>=hdlim);
   Xt=Xt(idx);
   Theta=Theta(idx);
   nSpike=length(Xt);
   
   % Set limits for plotting
   if nSpike==0
      disp('No spikes found in range. Stopping.')
      break;
   end
   
   if 0 % adaptively select plotting limits or not
      xlim1=xlimits(1)+xsz*floor((min(Xt)-xlimits(1))/xsz);
      xlim2=xlimits(1)+xsz*ceil((max(Xt)-xlimits(1))/xsz);
   else
      xlim1=xlimits(1);
      xlim2=xlimits(2);
   end
   
   % Find the fit for Theta=m*Xt+b+err where Theta is circular
   % First do a brute force search of precession rates getting
   % the circular variance of err in each case.
   if useCircularFit & length(Xt)>1
      k=0;
      M=[];
      CircVar=[];
      for m=-30.01:.01:30.01
         k=k+1;
         [cm,cr,cv]=cmean(Theta-m*Xt,360);         
         M(k)=m;
         CircVar(k)=cv;
      end
      
      % Look for local minimums to see if there is more than one
      % that is close to the global minimum
      idx=find(CircVar(1:end-2)>CircVar(2:end-1) & CircVar(2:end-1)<=CircVar(3:end));
      mopt=M(1+idx);
      cvopt=CircVar(1+idx);
      multOptFnd=0;
      if length(mopt)>1
         if length(find(min(cvopt)./cvopt>0.90))>1
            disp('Multiple near optimal precession rates found')
            multOptFnd=1;
            mopt
            cvopt
         end
      end
      
      if 1 & ntlim==1 % show fit magnitude as a function of m
         figure
         plot(M,CircVar)
         xlabel('Theta precession rate (deg/cm)')
         ylabel('Circular variance')
      end
      
      % Get the m with least variance. 
      m=mopt(find(cvopt==min(cvopt)));
      m=m(1); % pick an m if there are ties
      if useNegSlope==1 & m>0 & length(find(mopt<0))>0
         disp('Forcing use of negative slope')
         nmi=find(mopt<=0);
         m=mopt(find(cvopt==min(cvopt(nmi))));
         m=m(1)
      end
      if m>0 & length(find(mopt<0))>0
         if useNegSlope==1
            disp('Forcing use of negative slope')
            nmi=find(mopt<=0);
            m=mopt(find(cvopt==min(cvopt(nmi))));
            m=m(1)
         end
         if m>0 & useNegSlope==2 & multOptFnd
            [svar sidx]=sort(cvopt);
            if mopt(sidx(1))>0 & mopt(sidx(2))<0
               disp('For multiple optima, selecting negative slope')
               m=mopt(sidx(2));
            end
         end
      end  
      % And the final answer is (will the winner stand up please)
      m
      
      % Get b from the error mean and move it by full
      % cycles to be nearest 0 at the graph midpoint.
      b=cmean(Theta-m*Xt,360);   
      b=b-360*round((m*(xlim1+xlim2)/2+b)/360)
      
      % Unwrap theta for best alignment with fit
      ThetaUnwrapped=Theta-360*round((Theta-m*Xt-b)/360);
      
      % Optionally, redo fit using unwrapped theta
      if 1
         A=[Xt ones(length(Xt),1)];
         mb=A\(ThetaUnwrapped);
         m=mb(1)
         b=mb(2)
      end
      
   end
   
   if ~useCircularFit & length(Xt)>1
      % Use a fixed unwrap point (near CA1 peak activity)
      ThetaUnwrapped=mod(Theta+90,360)-90;
      A=[Xt ones(length(Xt),1)];
      mb=A\(ThetaUnwrapped);
      m=mb(1)
      b=mb(2)
   end
   
   % Produce some statistics
   nXt=length(Xt);
   cmean_theta_phase=cmean(ThetaUnwrapped,360)
   
   % Compute the sample correlation coefficient (r)
   % This is the fraction of unexplained variance in ThetaUnwrapped
   % assuming the fit is a least squares linear regression, which is
   % the case only if the optional linear fit is done when a circular
   % fit is used for unwrapping. The value is close to a measure of
   % unexplained variance in any case though.
   mean_X=mean(Xt);
   mean_theta_phase=mean(ThetaUnwrapped)
   Sxy = sum((Xt-mean_X).*(ThetaUnwrapped-mean_theta_phase));
   Sx=sum((Xt-mean_X).^2);
   Sy=sum((ThetaUnwrapped-mean_theta_phase).^2);
   r_statistic = Sxy/sqrt(Sx*Sy)
   
   df=length(Xt)-2
   t_statistic = sqrt(df)*r_statistic/sqrt(1-r_statistic^2)
   
   % Plot residuals
   if 0 & ntlim==1
      err=sort(ThetaUnwrapped-m*Xt-b);
      serr=sort(err);
      qerr=((1:length(serr))-0.5)/(length(serr)+1);
      zerr=(serr-mean(serr))./std(serr);
      qnorm=0.5+erf(zerr/sqrt(2))/2;
      
      figure
      plot(Xt,err,'o')
      xlabel('X (cm)')
      ylabel('residual (deg)')
      
      figure
      plot(qnorm,qerr,[0 1],[0 1],'--')
      xlabel('normal quantile')
      ylabel('residual quantile')
   end
   
   % Calculate a mutual information metric.
   % This is sensitive to the bin sizes used and
   % probably could use some clever data smoothing.
   xrng=xlimits(2)-xlimits(1);
   xcnt=zeros(1,1+ceil(xrng/xsz));
   tpcnt=zeros(1+ceil(360/tpsz),1);
   tpxcnt=zeros(length(tpcnt),length(xcnt));
   for k=1:nSpike
      xi=ceil(mod(Xt(k)-mean_X,xrng)/xsz);
      yi=ceil(mod(Theta(k)-cmean_theta_phase,360)/tpsz);
      xcnt(xi)=xcnt(xi)+1;
      tpcnt(yi)=tpcnt(yi)+1;
      tpxcnt(yi,xi)=tpxcnt(yi,xi)+1;
   end
   
   tpxinfo=0;      
   for xi=1:length(xcnt)
      for yi=1:length(tpcnt)
         ptpx=tpxcnt(yi,xi)/nSpike;
         if ptpx>0
            px=xcnt(xi)/nSpike;
            ptp=tpcnt(yi)/nSpike;
            tpxinfo=tpxinfo+ptpx*log(ptpx/(px*ptp));
         end
      end
   end
   theta_phase_x_info=tpxinfo/log(2)
   
   if 1 & ntlim==1 
      figure
      pcolor(tpxcnt)
      colorbar
      title('Spike counts by theta phase & location')
      ylabel('Theta phase bin')
      xlabel('Spatial location bin')
   end
   
   % Plot the fit with theta unwrapped
   if ntlim>1
      subplot(1,ntlim,n)
   else
      figure
   end
   
   % Do the plot
   plot(Xt,ThetaUnwrapped,'o')  
   title([num2str(tstrt) ' to ' num2str(tend) ' sec'])
   xlabel('X (cm)')
   axis([xlim1 xlim2 -360 360])
   
   if n==1
      ylabel('Theta phase (deg)')
   end
   
   box off
   hold on
   lh=plot(xlimits,m*xlimits+b,'r-');
   set(lh,'LineWidth',2);
   
   % label with the r value, slope, and information
   sigFlag=' ';
   if abs(t_statistic)>tStatP05
      sigFlag='*';
   end
   if abs(t_statistic)>tStatP01
      sigFlag='**';
   end
   
   text(xlim1+2,textY,['r = ' num2str(round(100*r_statistic)/100)],'FontSize',10)
   text(xlim1+2,textY-30,['m = ' num2str(round(100*m)/100) sigFlag],'FontSize',10)
   
end

% -------------------------------------------------------------------

function plot_linear_rates(cellId,xlimits,ylimits,dir,tlim,fileId)
global DATADIR CELLSZ MINTIME TSTRT TEND MAZETYPE PLOTGRAY

% Passed params are:
% cellId = id of the cell to plot
% xlimits = array [xlow xhigh] of region to plot
% ylimits = array [ylow yhigh] of region to plot
% tlim = matrix where each row is a from and to time interval
% fileId = id of file to read for spikes (see code below)

% Hard coded params are:
csz = CELLSZ; % x bin in cm for computing firing rates
hxlim = dir; % heading limit x coord (-1=right to left, 1=left to right)
hylim = 0; % heading limit y coord
hdlim = 0; % minimum of dot product of heading and [hxlim hylim]

MStates=dlmread([DATADIR 'test-baker-mouse-states.txt'],',',1,1);
T=MStates(:,1);
X=MStates(:,2);
Y=MStates(:,3);
Hx=MStates(:,4);
Hy=MStates(:,5);

hm=(T(end)-T(1))/(length(T)-1)/1000;

if isempty(fileId)
   fn=[DATADIR 'test-baker-spikes.txt'];
else   
   fn=[DATADIR 'test-baker-' fileId '-spikes.txt'];
end

Spikes=dlmread(fn,',',1,0);
Spikes=Spikes(find(Spikes(:,1)==cellId),:);
if isempty(Spikes)
   disp('No spikes found')
   return
end

figure
if PLOTGRAY
   plotColors='kkkkk';
else
   plotColors='bmgrk';
end
plotSymbols='osd^v';

for tidx=1:min(5,size(tlim,1))
   
   n=1+ceil((xlimits(2)-xlimits(1)-csz)/csz);
   NS=zeros(n,1);
   NC=1e-8*ones(n,1);
   swimIdx=1;
   
   for k=1:size(Spikes,1)
      
      id=Spikes(k,1);
      t=Spikes(k,2);
      
      if t>T(end) 
         break
      end
      
      nx=-1;
      x=nan;
      y=nan;
      while t>T(swimIdx)
         x=X(swimIdx);
         y=Y(swimIdx);
         hx=Hx(swimIdx);
         hy=Hy(swimIdx);
         swimIdx =swimIdx+1; 
         if xlimits(1)<=x & x<=xlimits(2) & ylimits(1)<=y & y<=ylimits(2) ...
               & Hx(swimIdx)*hxlim+hy*hylim>=hdlim ...
               & tlim(tidx,1)<=T(swimIdx)/1000 & T(swimIdx)/1000<=tlim(tidx,2)
            nx = 1+floor((x-xlimits(1))/csz);
            NC(nx)=NC(nx)+1;
         end
      end
      
      if nx>0 & xlimits(1)<=x & x<=xlimits(2) ...
            & ylimits(1)<=y & y<=ylimits(2) ...
            & hx*hxlim+hy*hylim>=hdlim
         if tlim(tidx,1)<=T(swimIdx)/1000 & T(swimIdx)/1000<=tlim(tidx,2)
            NS(nx)=NS(nx)+1;
         end
      end
   end
   
   FR=NS./(NC*hm);
   Xbin=xlimits(1)+csz/2+csz*(0:n-1);
   
   hold on
   plot(Xbin,FR,[plotColors(tidx) plotSymbols(tidx) '-'])
   box off
end
title(['Linear track firing rate for cell ' num2str(cellId)])
xlabel('X (cm)')
ylabel('Firing rate (Hz)')

if size(tlim,1)==2
   legend(['T=' num2str(tlim(1,1)) '-' num2str(tlim(1,2)) ' sec'], ...
      ['T=' num2str(tlim(2,1)) '-' num2str(tlim(2,2)) ' sec']);
end
if size(tlim,1)==3
   legend(['T=' num2str(tlim(1,1)) '-' num2str(tlim(1,2)) ' sec'], ...
      ['T=' num2str(tlim(2,1)) '-' num2str(tlim(2,2)) ' sec'], ...
      ['T=' num2str(tlim(3,1)) '-' num2str(tlim(3,2)) ' sec']);
end
if size(tlim,1)==4
   legend(['T=' num2str(tlim(1,1)) '-' num2str(tlim(1,2)) ' sec'], ...
      ['T=' num2str(tlim(2,1)) '-' num2str(tlim(2,2)) ' sec'], ...
      ['T=' num2str(tlim(3,1)) '-' num2str(tlim(3,2)) ' sec'], ...
      ['T=' num2str(tlim(4,1)) '-' num2str(tlim(4,2)) ' sec']);
end
if size(tlim,1)==5
   legend(['T=' num2str(tlim(1,1)) '-' num2str(tlim(1,2)) ' sec'], ...
      ['T=' num2str(tlim(2,1)) '-' num2str(tlim(2,2)) ' sec'], ...
      ['T=' num2str(tlim(3,1)) '-' num2str(tlim(3,2)) ' sec'], ...
      ['T=' num2str(tlim(4,1)) '-' num2str(tlim(4,2)) ' sec'], ...
      ['T=' num2str(tlim(5,1)) '-' num2str(tlim(5,2)) ' sec']);
end


% -------------------------------------------------------------------

function plot_layer_firing_rates(fileId)
global DATADIR CELLSZ MINTIME
MStates=dlmread([DATADIR 'test-baker-mouse-states.txt'],',',1,1);

cellSize = CELLSZ;

Tm=[0; MStates(:,1)];
X=[0; MStates(:,2)];
Y=[0; MStates(:,3)];
hm=(Tm(end)-Tm(1))/(length(Tm)-1)/1000;

radius = max(sqrt(X.^2+Y.^2));
radius = max(radius,max(abs(X)));
radius = max(radius,max(abs(Y)));
radius = ceil(radius);

if isempty(fileId)
   fn=[DATADIR 'test-baker-spikes.txt'];
else   
   fn=[DATADIR 'test-baker-' fileId '-spikes.txt'];
end

Spikes=dlmread(fn,',',1,0);
minId=min(Spikes(:,1))
maxId=max(Spikes(:,1))
numId=length(unique(Spikes(:,1)))

Ctr=[];
FRpeak=[];
Cir=radius*[cos(-pi:.02:pi)' sin(-pi:.02:pi)'];

totalSpikes=0;

n=1+ceil(2*radius/cellSize);
NS=zeros(n,n);
NC=1e-8*ones(n,n);

swimIdx=1;

for k=1:size(Spikes,1)
   t=Spikes(k,2);
   
   if t>Tm(end)
      break
   end
   
   while t>Tm(swimIdx)
      swimIdx =swimIdx+1;
      x=X(swimIdx-1);
      y=Y(swimIdx-1);         
      nx = ceil((radius+x)/cellSize);
      ny = ceil((radius+y)/cellSize);
      NC(ny,nx)=NC(ny,nx)+1;
   end
   
   NS(ny,nx)=NS(ny,nx)+1;
   totalSpikes = totalSpikes + 1;
end

totalSpikes

FR=(NC>=MINTIME/hm).*NS./(NC*hm)/(numId);
FRx=-radius-cellSize+cellSize*(1:size(FR,2));
FRy=-radius-cellSize+cellSize*(1:size(FR,1));

avg_spike_rate = totalSpikes/((Tm(end)-Tm(1))/1000*(numId))
max_spike_rate = max(max(FR))
min_spike_rate = min(min(FR))

figure
plot_color(FRx,FRy,FR,radius,16)
title([fileId ' firing rate per cell (Hz)'])

% -------------------------------------------------------------------

function plot_weighted_rates(layerId,pathwayId,asOf)
global DATADIR CELLSZ MINTIME TSTRT TEND MAZETYPE
disp(['plotting weighted rates of ' pathwayId])

disp('reading mouse states')
MStates=dlmread([DATADIR 'test-baker-mouse-states.txt'],',',1,1);
idx=find(MStates(:,1)>=TSTRT & MStates(:,1)<=TEND);
MStates=MStates(idx,:);

cellSize = CELLSZ;

Tm=[0; MStates(:,1)];
X=[0; MStates(:,2)];
Y=[0; MStates(:,3)];
hm=(Tm(end)-Tm(1))/(length(Tm)-1)/1000;

radius = max(sqrt(X.^2+Y.^2));
radius = max(radius,max(abs(X)));
radius = max(radius,max(abs(Y)));
radius = ceil(radius);
Cir=radius*[cos(-pi:.02:pi)' sin(-pi:.02:pi)'];

disp('reading weights')
% cols = postId, time, preId, comp, w1, w2, ...
SW = dlmread([DATADIR 'test-baker-' pathwayId '-synapse-weights.txt'],',',1,0);
if asOf>0
   asOfTime=asOf*1000
else
   asOfTime = SW(end,2)
end

disp('reading spikes')
fn=[DATADIR 'test-baker-' layerId '-spikes.txt'];
Spikes=dlmread(fn,',',1,0);

idx = find(Spikes(:,2)>=TSTRT & Spikes(:,2)<TEND);
Spikes=Spikes(idx,:);
if isempty(Spikes)
   disp('No spikes found in range')
   return
end

for tw=asOfTime
   
   disp(['processing as of time =' num2str(tw)])
   Ctr=[];
   FRpeak=[];
   totalSpikes=0;
   
   SWasof = SW(find(SW(:,2)>=tw-1 & SW(:,2)<=tw+1),:);
   if isempty(SWasof)
      disp(['No weights found at time ' num2str(tw)]);
      WstrtId=0;
      WendId=0;
      W=[];
   else
      WstrtId=min(SWasof(:,3))-1;
      WendId=max(SWasof(:,3));
      W=zeros(WendId-WstrtId,1);
      W(SWasof(:,3)-WstrtId)=SWasof(:,5);
   end
   
   n=1+ceil(2*radius/cellSize);
   NSW=zeros(n,n);
   NC=1e-8*ones(n,n);
   wspike=0;
   
   swimIdx=1;
   for k=1:size(Spikes,1)
      t=Spikes(k,2);
      preId=Spikes(k,1);
      
      if t>Tm(end)
         break
      end
      
      while t>Tm(swimIdx)
         swimIdx =swimIdx+1;
         x=X(swimIdx-1);
         y=Y(swimIdx-1);         
         nx = ceil((radius+x)/cellSize);
         ny = ceil((radius+y)/cellSize);
         NC(ny,nx)=NC(ny,nx)+1;
      end
      
      if WstrtId<=preId & preId<=WendId
         NSW(ny,nx)=NSW(ny,nx)+W(preId-WstrtId);
         wspike=wspike+W(preId-WstrtId);
      end
   end
   FR=(NC>=MINTIME/hm).*NSW./(NC*hm);
   FRx=-radius-cellSize+cellSize*(1:size(FR,2));
   FRy=-radius-cellSize+cellSize*(1:size(FR,1));
   
   max_min_spike_rate = [max(max(FR)) min(min(FR))]
   
   [sbits,skinfo,skspar]=infoMetrics(NSW,NC);
   sbits_skinfo_sparsity=[sbits skinfo skspar]
   
   if MAZETYPE==0
      figure
      plot_color(FRx,FRy,FR,radius,16)
      title(['Total ' layerId '-' pathwayId ' weighted afferent firing rate (Hz)'])
      %colormap(cool)
      %caxis([0 ceil(max(max(FR)))])
      colorbar
   end
   
   if MAZETYPE==1
      %figure
      rowsum=sum(FR,2);
      rowidx=find(rowsum==max(rowsum))
      plot(FRx(1:end-1)+cellSize/2,FR(rowidx,1:end-1),'o-')
      xlabel('X (cm)')
      ylabel('Weighted afferent firing rate (Hz)')
   end
   
end

% -------------------------------------------------------------------

function plot_weight_summary(layerId,pathwayId,setId)
global DATADIR CELLSZ MINTIME TSTRT TEND

plotActiveOnly = 1; % only plot weights of active cells

% read synapse weights
% cols = postId, time, preId, comp, w1, w2, ...
SW = dlmread([DATADIR 'test-baker-' pathwayId '-synapse-weights.txt'],',',1,0);
SW=SW(find(SW(:,2)>=TSTRT & SW(:,2)<TEND),:);
if isempty(SW)
   disp('No weights found in time range')
   return
end
T =unique(SW(:,2));

activeTitle='Synaptic';
if plotActiveOnly
   % read place field centers 
   % cols = id,x,y,peakrate,isactive
   PC = dlmread([DATADIR 'test-baker-' layerId '-pcloc' setId '.txt'],',',1,0);
   id0=PC(1,1)-1;
   
   % Figure out which weights go with active cells
   id=SW(:,3)-id0;
   isactive=find(PC(id,5)==1);
   activeTitle='Active synaptic';
end

% Find weights at each time point and summarize them.
% This is not a very efficient algorithm.
for k=1:length(T)
   t=T(k);
   if plotActiveOnly
      idx = find(SW(:,2)==t & PC(id,5)==1);
   else
      idx = find(SW(:,2)==t);
   end
   Wmean(k)=mean(SW(idx,5));
   Wmax(k)=max(SW(idx,5));
   Wstd(k)=std(SW(idx,5));
end

layerId
ending_weight_mean = Wmean(end)
ending_weight_max = Wmax(end)

if ending_weight_mean>0
   ending_max_mean_ratio = Wmax(end)/Wmean(end)
end

figure
[ax,h1,h2]=plotyy(T/1000,Wmean,T/1000,Wmax);
title([activeTitle ' weights for afferent layer ' layerId]);
xlabel('Time (sec)')
set(ax(1),'TickDir','out')
set(ax(2),'TickDir','out')
axes(ax(1))
ylabel('Mean Weight')
axes(ax(2))
ylabel('Max Weight')
% -------------------------------------------------------------------

function plot_afferent_weight(cellId,pathwayId)
global DATADIR TSTRT TEND

% read synapse weights
% cols = postId, time, preId, comp, w1, w2, ...
SW = dlmread([DATADIR 'test-baker-' pathwayId '-synapse-weights.txt'],',',1,0);
SW=SW(find(SW(:,2)>=TSTRT & SW(:,2)<=TEND & SW(:,3)==cellId),:);
if isempty(SW)
   disp('No weights found in time range for cell indicated')
   return
end

figure
plot(SW(:,2)/1000,SW(:,5));
title([pathwayId ' synaptic weight for afferent cell ' num2str(cellId)]);
xlabel('Time (sec)')
ylabel('Synaptic Weight')


% -------------------------------------------------------------------

function plot_weight_by_comp(layerId,pathwayId,asOf)
global DATADIR CELLSZ MINTIME TEND DM

disp(['Weights by dendrite compartment for ' pathwayId]);

load_dm;

% read weights
% cols = postId, time, preId, comp, w1, w2, ...
SW = dlmread([DATADIR 'test-baker-' pathwayId '-synapse-weights.txt'],',',1,0);
if asOf>0
   asOfTime=asOf*1000;
else
   asOfTime =max(SW(find(SW(:,2)<=TEND),2));
end
SW = SW(find(SW(:,2)>=asOfTime-1 & SW(:,2)<=asOfTime+1),:);
if isempty(SW)
   disp('As of time not found in weights dataset')
   return
end

Wtotal=zeros(max(SW(:,4)),1);
Wnsyn=zeros(length(Wtotal),1);
Wmax=zeros(length(Wtotal),1);
for k=1:size(SW,1)
   d=SW(k,4);
   w=SW(k,5);
   Wmax(d)=max(w,Wmax(d));
   Wtotal(d)=Wtotal(d)+SW(k,5);
   Wnsyn(d)=Wnsyn(d)+1;
end
Wdend=find(Wnsyn>0);
Wmean=Wtotal(Wdend)./Wnsyn(Wdend);
Wmax=Wmax(Wdend);

if 1
   figure
   plot(DM(Wdend,9),Wmean,'o');
   title(['Mean weights by compartment for afferent layer ' layerId])
   ylabel('Mean synaptic weight')
   xlabel('Dendrite compartment laminar coordinate (micron)')
end

if 1
   figure
   plot(DM(Wdend,9),Wmax,'o');
   title(['Max weights by compartment for afferent layer ' layerId])
   ylabel('Maximum synaptic weight')
   xlabel('Dendrite compartment laminar coordinate (micron)')
end

% -------------------------------------------------------------------

function plot_weight_by_dist(layerId,pathwayId,setId,asOf,pcx,pcy,wmin)
global DATADIR CELLSZ MINTIME TEND PLOTGRAY

disp(['Weights by distance for pathway ' pathwayId]);

plotActiveOnly = 1; % only plot weights of active cells

% read place field centers
% cols = id,x,y,peakrate,isactive
PC = dlmread([DATADIR 'test-baker-' layerId '-pcloc' setId '.txt'],',',1,0);
id0=PC(1,1)-1;

% read weights
% cols = postId, time, preId, comp, w1, w2, ...
SW = dlmread([DATADIR 'test-baker-' pathwayId '-synapse-weights.txt'],',',1,0);
if asOf>0
   asOfTime=asOf*1000;
else
   asOfTime =max(SW(find(SW(:,2)<=TEND),2));
end
SW = SW(find(SW(:,2)>=asOfTime-1 & SW(:,2)<=asOfTime+1),:);
if isempty(SW)
   disp('As of time not found in weights dataset')
   return
end
id=SW(:,3)-id0;
dist=sqrt((PC(id,2)-pcx).^2+(PC(id,3)-pcy).^2);

isactive=find(PC(id,5)==1);
widx=find(SW(:,5)>wmin & PC(id,5)==1);
active_weight_median = median(SW(isactive,5))
active_weight_mean = mean(SW(isactive,5))
active_weight_std = std(SW(isactive,5))
active_cells = length(isactive)
weights_above_wmin = length(widx)

% Make a linear fit of the form dist=b+m*weight+err
% for weights of active cells with weights above wmin
if isempty(widx)
   disp('No weights found above threshold')
   wfit=[0 -1];
   dfit=[-1 0];
else
   W=ones(length(widx),2);
   W(:,1)=SW(widx,5);
   mb=W\dist(widx);
   m=mb(1); b=mb(2);
   wfit=[-b/m 0];
   dfit=[0 b];
   % Get a measure of fit. Note that this is a
   % fit of distance as predicted by synaptic weight.
   fit_slope_weight_over_dist = 1/m % slope of weight~distance
   r_statistic = m*std(SW(widx,5))/std(dist(widx))
end

figure
ls1='.';
ls2='r-';
if PLOTGRAY
   ls1='k.';
   ls2='k-';
end
if wfit(2)<0 | dfit(1)<0
   if plotActiveOnly
      plot(dist(isactive),SW(isactive,5),ls1);
   else
      plot(dist,SW(:,5),ls1);
   end
   title([pathwayId ' synaptic weight vs center distance'])
   xlabel('Distance between place field centers (cm)')
   ylabel('Synaptic weight')
else
   if plotActiveOnly
      ph=plot(dist(isactive),SW(isactive,5),ls1,dfit,wfit,ls2);
   else
      ph=plot(dist,SW(:,5),ls1,dfit,wfit,ls2);
   end   
   set(ph(2),'LineWidth',2)
   xlim([0 ceil(max(dist))])
   title([pathwayId ' synaptic weight vs center distance'])
   xlabel('Distance between place field centers (cm)')
   ylabel('Synaptic weight')
   legend('actual weight',['linear fit for w>' num2str(wmin)])
   
   text(35,wfit(1)*0.75,['r = ' num2str(round(r_statistic*100)/100)],'FontSize',14)
   
end

% -------------------------------------------------------------------

function plot_spike_train(fromId,toId,fileId)
global DATADIR TSTRT TEND PLOTGRAY

ls0='o';
ls1='-';
if PLOTGRAY
   ls0='ko';
   ls1='k-';
end

if isempty(fileId)
   fn=[DATADIR 'test-baker-spikes.txt'];
else   
   fn=[DATADIR 'test-baker-' fileId '-spikes.txt'];
end

Spikes=dlmread(fn,',',1,0);
if isempty(Spikes)
   disp('No spikes read')
   return
end
idx=find(Spikes(:,1)>=fromId & Spikes(:,1)<=toId);
Spikes=Spikes(idx,:);

idx=find(Spikes(:,2)>=TSTRT & Spikes(:,2)<=TEND);
Spikes=Spikes(idx,:);

if isempty(Spikes)
   fromId
   toId
   disp('No spikes found')
   return
end

if fromId==toId
   idTitle=['cell ' num2str(fromId)];
else
   idTitle=['cells ' num2str(fromId) ' to ' num2str(toId)];
end

minId=min(Spikes(:,1))
maxId=max(Spikes(:,1))
spikes_read = size(Spikes,1)
Theta=mod(Spikes(:,2)+125/2,125)/125*360-180;

if 0 & fromId~=toId
   figure
   plot(Spikes(:,2)/1000,Spikes(:,1),'k.');
   title('Spike train scatter plot')
   xlabel('Time (s)')
   ylabel('Cell numeric identifier')
end

if 1 & fromId==toId
   figure
   plot(Spikes(:,2)/1000,Theta,ls0)
   ylim([-180 180])
   title('Spike scatter plot')
   ylabel('Theta phase (degree)')
   xlabel('Time (sec)')
end

tmax=max(Spikes(:,2));
tmin=min(Spikes(:,2));
avgRate = (spikes_read-1)/(tmax-tmin)*1000
ratePerCell=avgRate/(maxId-minId+1)
mean_theta_phase=cmean(Theta,360)

if 1 
   if toId-fromId<50
      binSize = 500;
   else
      binSize=125;
   end
   
   PR=zeros(ceil(tmax/binSize),1);
   
   for k=1:size(Spikes,1)
      n=floor(Spikes(k,2)/binSize)+1;
      PR(n)=PR(n)+1;
   end
   
   PR=PR/(maxId-minId+1)*1000/binSize;
   PRstrt = max(1,floor(TSTRT/binSize));
   PRend = min(length(PR),ceil(TEND/binSize));
   
   figure
   bar(binSize/2000+binSize/1000*(PRstrt:PRend),PR(PRstrt:PRend),ls1)
   title(['Firing rate for ' idTitle]);
   ylabel(['Firing rate (Hz) using bin size ' num2str(binSize) ' ms'])
   xlabel('Time (s)')
end

if 1 & fromId==toId & size(Spikes,1)>1
   frtau=125;
   ortau=60e3;
   h=5;
   e1=exp(-h/frtau);
   e2=exp(-h/ortau);
   s1=0;
   s2=0;
   nstep=floor((Spikes(end,2)-Spikes(1,2))/h);
   TR=zeros(nstep,1);
   FR=zeros(nstep,1);
   OR=zeros(nstep,1);
   
   k=0;
   nsp=1;
   or=0;
   for t=Spikes(1,2)+h:h:Spikes(end,2)
      s2=e1*(s2+h/frtau*s1);
      s1=s1*e1;
      while nsp<=size(Spikes,1) & t>=Spikes(nsp,2)
         s1=s1+1;
         nsp=nsp+1;
      end
      k=k+1;
      TR(k)=t/1000;
      FR(k)=s2/frtau*1000;
      OR(k)=or*e2+(1-e2)*FR(k);
      or=OR(k);
   end
   if 1
      figure
      bar(TR,FR,ls1)
      title(['Kernel smoothed firing rate for cell ' num2str(fromId)])
      ylabel('Firing rate (Hz)')
   end
   if 1
      figure
      plot(TR,OR,ls1)
      xlabel('Time (sec)')
      ylabel('Ongoing rate (Hz)')
   end
end

if 1
   figure
   if 0 & fromId==toId
      hbins=-144:72:144;
   else
      hbins=-170:20:170;
   end
   histTheta=hist(Theta,hbins);
   bar(hbins,histTheta,ls1)
   xlabel('Theta phase (degree)')
   ylabel('Count')
   title(['Theta rhythm firing for ' idTitle]);
   if 0
      figure
      polar([hbins hbins(1)]/180*pi,[histTheta histTheta(1)],'.-')
      title('Theta rhythm firing histogram in polar coordinates')
   end
end

if 0
   figure
   bar(-10:5:10,hist(mod(Spikes(:,2)+25/2,25)-25/2,-10:5:10),ls1)
   xlabel('Gamma phase offset (ms)')
   ylabel('Count')
   title(['Gamma rhythm firing for ' idTitle]);
end

if 0 & fromId==toId
   figure
   plot(Spikes(2:end,2),Spikes(2:end,2)-Spikes(1:end-1,2),'.')
   xlabel('time (ms)')
   ylabel('ISI (ms)')
end

if 0 & fromId==toId
   figure
   plot(Spikes(2:end,2)/1000,1000./(Spikes(2:end,2)-Spikes(1:end-1,2)),'.')
   xlabel('time (s)')
   ylabel('ISI-based freq (Hz)')
end

if 0
   % sort by cell id and time - approx ISI for all cells
   [sval, sidx] = sort(Spikes(:,1)*1e7+Spikes(:,2));
   isi = Spikes(sidx(2:end),2)-Spikes(sidx(1:end-1),2);
   isi(find(Spikes(sidx(2:end),1)~=Spikes(sidx(1:end-1),1)))=0;
   figure
   hist(isi(find(isi>0 & isi<250)),40);
   title('Spike train ISI histogram')
   xlabel('ISI (ms)')
   ylabel('Count')
end

% ---------------------------------------------------------------

function plot_dendrite_spikes(tprange)
global DM SRCDIR NCOMP DATADIR TSTRT TEND PLOTGRAY

% Params
% tprange = range of theta phases used for processing ([low high])

disp('Plotting dendrite spikes')
tprange

load_dm % get dendrite morphology table

ls0='o';
ls1='-';
if PLOTGRAY
   ls0='ko';
   ls1='k-';
end

% Read dendrite spike records. Times are in msec.
% cols: 1-neuron id, 2-spike time, 3-dendrite comp, 4-time since cell spike
SP = dlmread([DATADIR 'test-baker-dendrite-spikes.txt'],',',1,0);

% Try to drop a partial final record, if any
% Then make sure there is something to process.
if size(SP,1)>1 & SP(end-1,2)>=SP(end,2)
   SP=SP(1:end-1,:);
end
if isempty(SP)
   disp('No dendrite spikes could be read')
   return
end

% Apply theta phase mask
Theta=mod(SP(:,2)+125/2,125)/125*360-180;
idx=find(Theta>=tprange(1) & Theta<=tprange(2));
SP=SP(idx,:);
Theta=Theta(idx);

% Get time step size plus start and end times
idx=find(SP(1:end-1,2)<SP(2:end,2));
h=min(SP(idx+1,2)-SP(idx,2));
Tstrt=max(TSTRT,min(SP(:,2)));
Tend=min(TEND,max(SP(:,2)));

% Categorize dendrite spikes vs BPAP spikes based on time difference
DS = SP(find(SP(:,4)>=10 & SP(:,2)>=TSTRT & SP(:,2)<=TEND),:);
BP = SP(find(SP(:,4)<10 & SP(:,2)>=TSTRT & SP(:,2)<=TEND),:);

nDS=size(DS,1);
nBP=size(BP,1);
nDend=max(SP(:,3));
if nDS==0 & nBP==0
   disp('No dendrite spikes in time range')
   return
end

% Put spikes in bins of 1 ms for counting by time and compartment
nBin=1;
if nBP>0
   nBin=1+floor(BP(end,2)-Tstrt);
end
BPM=sparse(1+floor(BP(:,2)-Tstrt),BP(:,3),ones(nBP,1),nBin,nDend);
BPCnt1=sum(BPM,1);
BPCnt2=sum(BPM,2);
total_BPAP=nBP

nBin=1;
if nDS>0
   nBin=1+floor(DS(end,2)-Tstrt);
end
DSM=sparse(1+floor(DS(:,2)-Tstrt),DS(:,3),ones(nDS,1),nBin,nDend);
DSCnt1=sum(DSM,1);
DSCnt2=sum(DSM,2);
total_dendrite_spikes=nDS

if 0
   figure
   plot(1:length(BPCnt1),BPCnt1,ls0)
   xlabel('Dendrite compartment number')
   ylabel('BPAP count')
   figure
   plot(1:length(DSCnt1),DSCnt1,ls0)
   xlabel('Dendrite compartment number')
   ylabel('Dendritic spike count')
end

if 1
   figure
   scatter(DM(:,9),100-BPCnt1/max(BPCnt1)*100,DM(:,5)*100,ls0)
   xlabel('Laminar coordinate (micron)')
   ylabel('Backpropagation failures (%)')
end

if 1
   figure
   scatter(DM(:,9),DSCnt1/(Tend-Tstrt)*1000,DM(:,5)*100,ls0)
   xlabel('Laminar coordinate (micron)')
   ylabel('Dendrite compartment spike freq (Hz)')
end

figure
k=0;
dr=0.4;
for r=dr:dr:1.6
   k=k+1;
   R(k)=r-dr/2;
   ridx=find(DM(:,5)>r-dr & DM(:,5)<=r);
   dsf=DSCnt1(ridx)/(Tend-Tstrt)*100;
   if isempty(dsf)
      DSFmean(k)=nan;
      DSFcnt(k)=nan;
   else
      DSFmean(k)=mean(dsf);
      DSFcnt(k)=length(dsf);
   end
end
bar(R,DSFmean,ls1)
ylabel('Mean dendrite spike freq (Hz)')
xlabel('Dendrite radius (micron)')

tbin=1e3;
nbin=floor(tbin/h);
k=0;
for i=1:nbin:length(DSCnt2)
   k=k+1;
   DSbin(k) = sum(DSCnt2(i:min(length(DSCnt2),i+nbin-1)));
   Tbin(k)=(i-1)*h+Tstrt;
end

figure
bar(Tbin/1000,DSbin/tbin*1000,ls1)
ylabel('Total dendrite spike frequency (spikes/sec)')
xlabel('Time (sec)')


% ---------------------------------------------------------------
% Support functions

function plot_color(X,Y,Z,r,n)
% plot color graph clipped to the radius provided

global CELLSZ MAZETYPE PLOTGRAY PLOTCONTOUR

if PLOTCONTOUR %plot via contour
   maxz=max(max(Z));
   contour_levels=maxz*[0.9 0.7 0.5 0.2];
   if PLOTGRAY
      contour_levels=floor(contour_levels);
   end
   contour_levels
   [cs,h]=contourf(X,Y,Z,contour_levels,'k-');
   if PLOTGRAY
      clabel(cs,h,'color','w');
   end
   xlabel('X (cm)')
   ylabel('Y (cm)')
   hold on
   axis([-30 30 -30 30])
   axis square
   if PLOTGRAY
      colormap(gray)
   end
   caxis([0 maxz])
   colorbar
   return
end

if MAZETYPE==1
   % draw square maze
   pcolor(X,Y,Z)
   axis square
   colorbar
   xlabel('X (cm)')
   ylabel('Y (cm)')
   return
end

% Draw circle maze by dividing cells
% up into small squares so that circle outline
% can be made evident even though pcolor is used.
csz=CELLSZ;
y=Y(1)+ csz/n*(0:n*length(Y));
x=X(1)+ csz/n*(0:n*length(X));
z=nan*zeros(length(y),length(x));

for i1=1:length(Y)
   for i2=1:length(X)        
      k1=1+n*(i1-1);
      k2=1+n*(i2-1);
      for j1=k1:k1+n-1
         for j2=k2:k2+n-1
            if (x(j2)+csz/(2*n))^2 + (y(j1)+csz/(2*n))^2 <= (r+csz/(2*n))^2
               z(j1,j2) = Z(i1,i2);
            end
         end
      end
   end
end

if 1
   pcolor(x,y,z)
   shading flat
   
   if PLOTGRAY 
      colormap(gray)
      brighten(0.5)
   end
   
   xlim([X(1)-csz/2 X(end)+csz/2]);
   ylim([Y(1)-csz/2 Y(end)+csz/2]);
   axis square
   colorbar
   xlabel('X (cm)')
   ylabel('Y (cm)')
end

% ---------------------------------------------------------------

function [sbits,skinfo,sparsity]=infoMetrics(NS,NC)

% This computes information metrics associated with
% place cell firing. This version does not have corrections
% for low firing-rate place cells but does conform to
% the usual methods of measuring spatial information.

% Inputs
% NS = array with number of spikes in a spatial bin
% NC = array with occupancy counts in a spatial bin
%
% Outputs
% sbits = a spatial information metric similar to the Skaggs
%   metric except that the resulting information content of a
%   spike can never exceed the information needed to specify a
%   visited spatial location, something that can theoretically 
%   happen for the Skaggs mutual information metric when spatial 
%   occupancies are not uniformly distributed.
%
% skinfo = Skaggs mutual information metric computed without smoothing.
% sparsity = sparsity metric described by Skaggs et al.

% See Skaggs WE, McNaughton BL, Wilson MA, Barnes C (1996) 
% Theta phase precession in hippocampal populations and the compression 
% of temporal sequences. Hippocampus 6:149-172 for a discussion of
% information metrics.

ns=reshape(NS,1,prod(size(NS))); % get number of spikes as array
nc=reshape(NC,1,prod(size(NC))); % get occupancy count as array

% Take out cases where 0 denominators arise. Note that NC<1
% should be treated the same as zero since low values (1e-8)
% are used to avoid zero denominators elsewhere.

nzi=find(nc>=1);
nco=nc(nzi); % occupancy counts within visited bins
ncv=length(nco); % number of cells visited

osi=find(nc>=1 & ns>0); % note: check on nc should be redundant
ncos=nc(osi); % occupancy counts in visited bins with spikes
nsos=ns(osi); % spike counts in visited bins with spikes

totCnt=sum(nco); % total occupancy counts
totSpk = sum(ns); % total spikes
if totSpk==0 | ncv==0
   sbits=0;
   skinfo=0;
   sparsity=0;
   return
end

po=nco/totCnt; % probability of occupancy in all bins visited
pos=ncos/totCnt; % probabilty of occupancy given a spike occurred
pfs=nsos/totSpk; % probability of a given spike being in a bin
frs=(nsos./ncos); % firing rate in a bin for bins where a spike occurred
mfr=(totSpk/totCnt); % mean firing rate in spikes per time-tick units
lambda=frs/mfr; % ratio of bin firing rate to mean rate

sbits=(sum(pfs.*log(pfs))+log(ncv))/log(2); % entropy delta for a spike
skinfo=sum(pos.*lambda.*log(lambda))/log(2); % Skaggs mutual info metric
sparsity=sum(pos.*lambda)^2 / sum(pos.*(lambda.^2)); % sparsity (see Skaggs)

% -------------------------------------------------------------------

function [m,r,v]=cmean(x,c)
% circular mean and variance
% x = input values
% c = circular range of x (usually 2*pi or 360)
% m = circular mean of x
% r = mean radius of x
% v = circular variance of x

if find(size(x)>1)~=1
   if prod(size(x))==0
      warning('Sample empty')
      m=nan;
      r=nan;
      v=nan;
      return
   end
   error('sample must be one-dimensional')
end

xs=sum(sin(2*pi/c*x));
xc=sum(cos(2*pi/c*x));
if abs(xc)>eps
   m=c*atan2(xs,xc)/(2*pi);
else
   m=c/4*sign(xs);
end
r=sqrt(xs^2+xc^2)/length(x);
v=1-r;

