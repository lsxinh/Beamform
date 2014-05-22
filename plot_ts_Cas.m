%% beamforming output for whole year
% fl=dir(['/data/geophys/scratch/ep8g10/pbm_MSci/20*']);
% sang=zeros(138,length(fl)*33); %138 is number of frequencies.  33 is number of timewindows each day
% spwr=sang;ssl=sang;
% for i=1:length(fl)
%  load (['/data/geophys/scratch/ep8g10/pbm_MSci/',fl(i).name]);
%  %size of beam is 181x138x41x33: direction:freq:slowness:window
% for ic= [5:length(I)-5];  
%    for iday=1:size(beam,4);
%        %average/smooth over frequency range and cut down to one time window
%         tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday),2)));
%        % tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday+[-7:1:7]))));
%         tre=10*log10(squeeze((tre)));
%         tre(isnan(tre))=0;
%        %tre is now a second order tensor direction:slowness
% %        tre=tre-max(max(tre));
%         [ii,j]=find(tre==max(max(tre)),1); %finds max of tre (and associated direction and slowness)
%         sang(ic,(i-1)*33+iday)=theta(ii); %angle at max beampower, for this frequency and timestep
%         spwr(ic,(i-1)*33+iday)=max(max(tre)); %max beampower
%         ssl(ic,(i-1)*33+iday)=SL(j);  %slowness at max beampower, for this frequency and timestep
%         %isub=isub+1; %islow(isub)=SL(j(1));
%    end
%    
% end
% display(fl(i).name);
% end
% 
% eval(['save /data/geophys/scratch/jn6g09/Beamforming/ts_Cas ' 'sang spwr ssl']);

%% contour plot
clear all;close all;
load (['/data/geophys/scratch/jn6g09/Beamforming/ts_Cas.mat']); %output from section 1
load (['/data/geophys/scratch/jn6g09/Beamforming/WY_Cas_var.mat']) % beamforming variables- slowness, frequency

fl=dir(['/data/geophys/scratch/ep8g10/pbm_MSci/20*']);

%choose dates to plot
start=[2013 5 1 0 0 0];
stop=[2013 5 31 0 0 0];

%choose frequencies to plot
freqs=frq(I)';
startf=0.0313;
stopf=0.15;

%choose which plots to make (=1 to plot)
plotint=1;
plotsl=0;
plotaz=0;

%create a time vector
clear year day time_pre yearstrt
for i=1:length(fl)
yeari=str2double(fl(i).name(1:4));year(i,:)=yeari;
dayi=str2double(fl(i).name(6:8));day(i,:)=dayi;
end
for j=1:length(fl)
yearstrt(j)=datenum(year(j),1,1,0,0,0);
time_pre(j)=yearstrt(j)+day(j)-1;
end 
nt=size(beam,4)
for i=1:length(fl)
    time((i-1)*nt+1)=time_pre(i);
    for j=1:nt-1
        time((i-1)*nt+1+j)=time_pre(i)+(1/nt)*j;
    end
end
time=time';
dv1=datevec(time);


%cut data down into chosen date
startnum=datenum(start);
stopnum=datenum(stop);
starti=find(round(10.*time)==round(10.*startnum),1,'first');
stopi=find(round(10.*time)==round(10.*stopnum),1,'last');
time=time(starti:stopi);
dv=dv1(starti:stopi);
spwr=spwr(:,starti:stopi);
ssl=ssl(:,starti:stopi);
sang=sang(:,starti:stopi);

%cut data into chosen frequnecies
[c startfi] = min(abs(freqs-startf));
[c stopfi]= min(abs(freqs-stopf));
freqs=freqs(startfi:stopfi);
spwr=spwr(startfi:stopfi,:);
ssl=ssl(startfi:stopfi,:);
sang=sang(startfi:stopfi,:);

%remove zero data at beginning and end frequencies caused by smoothing
%over a frequency range
val=find(any(spwr,2)==1);
freqs=freqs(val);
spwr=spwr(val,:);
ssl=ssl(val,:);
sang=sang(val,:);

%tidy up
clear c dayi frq I SL day dv1 i j nt startfi stopfi startnum stopnum
clear starti stopfi stopi stopnum theta time_pre val yeari yearstrt year


%power time series
if plotint==1
figure(1)
pcolor(time,freqs,spwr); shading flat;
hold on;
y=zeros(length(time))+min(freqs);%plotting dots at each time step somehow
%forces it to plot it correctly
plot(time,y,'o','MarkerSize',1)
ylabel('frequency (Hz) '); 
xlabel('date');
caxis([9 16])
colorbar;
ylabel(colorbar,'power (dB)');
datetick('x','dd/mm/yy','keeplimits','keepticks')
title('time-series of maximum intensity');
set(gcf,'Color',[1,1,1]);
end

%slowness time series
if plotsl==1
figure(2)
pcolor(time,freqs,ssl); shading flat;
hold on;
plot(time,y,'o','MarkerSize',1)
ylabel('frequency (Hz)');
xlabel('date');
colorbar
ylabel(colorbar,'slowness (s/km)');
datetick('x','dd/mm/yy','keeplimits','keepticks')
title('time-series of slowness at maximum intensity');
set(gcf,'Color',[1,1,1]);
end

%azimuth time series
if plotaz==1
figure(3)
pcolor(time,freqs,sang); shading flat;
hold on;
plot(time,y,'o','MarkerSize',1)
ylabel('frequency (Hz)');
xlabel('date');
colorbar
caxis([0 360])
ylabel(colorbar,'azimuth (\circ)');
datetick('x','dd/mm/yy','keeplimits','keepticks')
title('time-series of azimuth at maximum intensity');
set(gcf,'Color',[1,1,1]);
end
