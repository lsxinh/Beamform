%% radar plots of 1 day

%% Actual beamforming data
% load day 82 (from start of dataset)
%load (['/data/geophys/scratch/ep8g10/pbm_MSci/2012_302.mat']);
% day 183
load (['/data/geophys/scratch/ep8g10/pbm_MSci/2012_302.mat']);

% %% data reshaped for radar plot (slowness against azimuth)
% isub=0;
 for ic= [5:5:length(I)-5];  
    for iday=3:30;
         tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday+[-2:1:2]),2)));
%        % tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday+[-7:1:7]))));
         tre=10*log10(squeeze(mean(tre,3)));
% %        tre=tre-max(max(tre));
%         [i,j]=find(tre==max(max(tre)),1);  
%         isub=isub+1; %islow(isub)=SL(j(1));
% 
% %% Plot
% %make a radar plot of the beamformer output, slowness concentric circles,
% %azimuthal angle
% 
time=((iday)-1)*24/33
% %beamRt=10*log((tre./min(min(tre))));for actual data
% % setup grid
% %vector of slowness in s/km
% SL=0:0.0098:0.4;
% % azimuth in deg
% theta=0:2:360;
% 
 [XX,RAD]=meshgrid(SL,theta*pi/180);
 [X1,Y1]=pol2cart(RAD,XX);
%  
% figure(1)
%  
%  %set background colour
 set(gcf,'Color',[1,1,1]);
% 
 h=polar([0 2*pi], [0 max(SL)]); axis('ij'), delete(h), hold on, view([-90 90]); 
 pcolor(X1,Y1,real(tre));shading flat;colorbar;
 per=1./frq(I(ic)) 
 title(['T= ',num2str(per),'s','Time=',num2str(time),' hr']); %for plane wave
 pause(.1)
    end
 end
   
%% for set frequencies and times

return
figure(12)
 %set background colour
set(gcf,'Color',[1,1,1]);

%ic =70 ; % for T = 6s => f=  0.1667Hz
ic =11 ; % f =0.05Hz, T=20s
%ic = 36;% f = 0.1Hz, T=10s
%ic = 31;
%iday =15 ; % 
%iday=27;
iday = 17;
time_b = ((iday -2)./33).*24; %start time 
time_e = ((iday +2)./33).*24; %end time

    tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday+[-2:1:2]),2)));
    tre=10*log10(squeeze(mean(tre,3)));
    [i,j]=find(tre==max(max(tre)),1);  
  
h=polar([0 2*pi], [0 max(SL)]); axis('ij'), delete(h), hold on, view([-90 90]); 
pcolor(X1,Y1,real(tre));shading flat;
colorbar; ylabel(colorbar,'power/ dB');
per=1./frq(I(ic)) 
title(['T= ',num2str(round(per)),'s, day 183, time average between ',num2str((round(time_b.*10))./10), ' and ', num2str((round(time_e.*10))./10), 'hr']); 
        

        
