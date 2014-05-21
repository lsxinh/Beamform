%% beamforming output for whole year (different method)
% fl=dir(['../pbm_MSci/20*']);
% sang=zeros(138,length(fl)*33);
% spwr=sang;ssl=sang;
% for i=1:length(fl)
%  load (['../pbm_MSci/',fl(i).name]);
%  
% for ic= [5:length(I)-5];  
%    for iday=1:size(beam,4);
%         tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday),2)));
%        % tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday+[-7:1:7]))));
%         tre=10*log10(squeeze((tre)));
%         tre(isnan(tre))=0;
% %        tre=tre-max(max(tre));
%         [ii,j]=find(tre==max(max(tre)),1); 
%         sang(ic,(i-1)*33+iday)=theta(ii);
%         spwr(ic,(i-1)*33+iday)=max(max(tre));
%         ssl(ic,(i-1)*33+iday)=SL(j);
%         isub=isub+1; %islow(isub)=SL(j(1));
%    end
%    
% end
% display(fl(i).name);
% end
% 
% eval(['save /data/geophys/scratch/ep8g10/CODES/dud ' 'sang spwr ssl']);

%% contour plot
 load (['/data/geophys/scratch/ep8g10/CODES/dud.mat']); %output from section 1
load (['/data/geophys/scratch/ep8g10/pbm_MSci/WY_Cas_var.mat']) % beamforming variables- slowness, frequency
% 
% %power time series
% % figure(1)
timev = [0:298/9834:297.9697]; %time vector
% pcolor(timev,frq(I(5:133)),spwr(5:133,:)); shading flat;
% ylabel('frequency/ Hz '); xlabel('time /days');
% colorbar;ylabel(colorbar,'power/ dB');
% caxis([2 16]);
% title('time-series of intensity for f = 0.03Hz to 0.3Hz (method (2))');
% %set background colour
% set(gcf,'Color',[1,1,1]);
% 
% % % %
% figure(2)
% hold on
% %subplot(3,2,1)
% pcolor(timev(1981:3631),frq(I(5:60)),spwr(5:60,1981:3631)); shading flat;
% ylabel('frequency/ Hz '); xlabel('time /days');colorbar;
% ylabel(colorbar,'power/ dB');
% caxis([8 16]);
% title(['time-series of intensity for days ',num2str((round(timev(1981)))),'- ',num2str((round(timev(3631)))) ]);
% %set background colour
% set(gcf,'Color',[1,1,1]);
% 
% figure(3)
% mind=5776;
% maxd=6271;
% %subplot(3,2,2)
% pcolor(timev(mind:maxd),frq(I(5:60)),spwr(5:60,mind:maxd)); shading flat;
% ylabel('f/ Hz '); xlabel('time /days');%colorbar;ylabel(colorbar,'power/ dB');
% caxis([8 16]);
% title(['intensity for days ',num2str((round(timev(mind)))),'- ',num2str((round(timev(maxd)))) ]);
% %set background colour
% set(gcf,'Color',[1,1,1]);

figure(34)
mind= 6023 %82.5d  %2641 = 80d; %2674  = 81d
maxd=6089 %6073 84d  %2772 84d %2740 =83d;
%subplot(3,2,2)
pcolor(timev(mind:maxd),frq(I(5:60)),spwr(5:60,mind:maxd)); shading flat;
ylabel('f/ Hz '); xlabel('time /days');%colorbar;
caxis([8 16]);
title(['intensity for days ',num2str((round(timev(mind)))),'- ',num2str((round(timev(maxd)))) ]);
%set background colour
set(gcf,'Color',[1,1,1]);
% 
% % %slowness time series
% figure(4)
% pcolor(timev,frq(I(5:133)),ssl(5:133,:)); shading flat;
% ylabel('frequency/ Hz'); xlabel('time /days');colorbar;
% title('time-series of slowness for f = 0.03Hz to 0.3Hz (method(2))');
% ylabel(colorbar,'slowness/ s(km)^{-1}');
% %set background colour
%  set(gcf,'Color',[1,1,1]);
% %  
% %figure(5)
% subplot(3,2,3)
% pcolor(timev(1981:3631),frq(I(5:60)),ssl(5:60,1981:3631)); shading flat;
% ylabel('frequency/ Hz '); xlabel('time /days');colorbar;
% %caxis([8 16]);
% ylabel(colorbar,'slowness/ s(km)^{-1}');
% title(['slowness for days ',num2str((round(timev(1981)))),'- ',num2str((round(timev(3631)))) ]);
% %set background colour
% set(gcf,'Color',[1,1,1]);
% 
% %figure(6)
% subplot(3,2,4)
% mind=5776;
% maxd=6271;
% pcolor(timev(mind:maxd),frq(I(5:60)),ssl(5:60,mind:maxd)); shading flat;
% ylabel('f/ Hz '); xlabel('time /days')%colorbar;
% %ylabel(colorbar,'slowness/ s(km)^{-1}');
%  title(['slowness for days ',num2str((round(timev(mind)))),'- ',num2str((round(timev(maxd)))) ]);
% %set background colour
% set(gcf,'Color',[1,1,1]);
% 
% % azimuths time series
% figure(4)
% pcolor(timev,frq(I(5:133)),sang(5:133,:)); shading flat;
% ylabel('frequency/ Hz'); xlabel('time /days');colorbar;
% title('time-series of azimuth for f = 0.03Hz to 0.3Hz (method(2))');
% ylabel(colorbar,'azimuth/ \circ'); caxis([0 360]);
% set(colorbar,'YTick',[0:45:360]);
% %set background colour
%  set(gcf,'Color',[1,1,1]);


% figure(5)
% %subplot(3,2,5)
% pcolor(timev(1981:3631),frq(I(5:60)),sang(5:60,1981:3631)); shading flat;
% ylabel('frequency/ Hz '); xlabel('time /days');colorbar;
% %caxis([8 16]);
% ylabel(colorbar,'azimuth/ \circ'); caxis([0 360]);
% title(['azmiuth for days ',num2str((round(timev(1981)))),'- ',num2str((round(timev(3631)))) ]);
% %set background colour
% set(gcf,'Color',[1,1,1]);
% 
% figure(6)
% %subplot(3,2,6)
% mind=5776;
% maxd=6271;
% pcolor(timev(mind:maxd),frq(I(5:60)),sang(5:60,mind:maxd)); shading flat;
% ylabel('frequency/ Hz'); xlabel('time /days');colorbar;
% ylabel(colorbar,'azimuth/ \circ'); 
% caxis([0 360]);
% title(['azimuth for days ',num2str((round(timev(mind)))),'- ',num2str((round(timev(maxd)))) ]);
% %set background colour
% set(gcf,'Color',[1,1,1]);

  
% %frequency spectrum
% figure(7); xlabel('frequency /Hz'); ylabel('power?');
% figure(7); plot(spwr(:,1),'k*');title('day 1 (part) frequency spectrum?');