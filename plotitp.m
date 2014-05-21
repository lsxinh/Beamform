clear all
for day=1 %change date here
%      load(['../processed1/backup',num2str(day),'.mat']);  
     load(['/Volumes/JFNexternaldrive/Seismic/CASCADIAresponse/Jan/TestDay/Event_2012_309/processedbm_cutstations/LHZ.2012.309'],'-mat');
%     load beamallyear beamsave
%     beam=beamsave/366;
%     beam=sum(beam,4);
figure(3)
y=SL;  x=theta*pi/180;
Y1=ones(length(x),1)*y; X1= x*ones(1,length(y));
[X2,Y2]=pol2cart(X1,Y1);

isub=0
clear iang islow inoise
icnt=0;
for ii=[5:3:50];
icnt=icnt+1;
iss=find(frq(I)>=1/ii);
issc(icnt)=iss(1);
end
for ic= [61 98 133 166 201 253 335 405]
        
    %88 49 32 22 16 11 8 5];%6:20:500 %80 %101; 45:45 %Nfrq %70 %:70 %10%:5:30 %:size(beam,2)
    %ic=20  ic=70 ic=60 ic=107

    for iday=2%:2 %2:Ntime %270:330 %300 %10:1:Ntime-5 %1:Ntime% 40 %1:1:210 %Ntime-5 %1:Nday

      %    tre=double(squeeze(mean(beam(:,ic+[-2:1:2],:,iday+[-5:5]),2)));
     %   tre=double(squeeze(mean(beam(:,ic+[-5:1:5],:,iday+[-1:1]),2)));
        tre=double(squeeze(mean(beam(:,ic+[-6:1:6],:,iday),2)));
        %        tre=double(squeeze(mean(beam(:,60:77,:,iday+[-1:1]),2)));
        %tre=squeeze(mean(beam(:,ic+[-20:1:20],:,iday),2));
        tre=10*log10(squeeze(mean(tre,3)));
       % for ik=1:14
        %    tre(:,ik)=tre(:,ik);%-5; %0.25*(8-ik);
        %end
        tre=tre-max(max(tre));
        [JJ,J]=find(tre==0); JJ=JJ(1);J=J(1);
         isub=isub+1; iang(isub)=JJ(1); islow(isub)=SL(J(1));inoise(isub)=mean(mean(tre));
 %       subplot(4,4,isub)
         hold off
        h=polarpeter([0 2*pi], [0 0.4]); axis('ij'), delete(h), hold on, view([-90 90])
        pcolor(X2,Y2,tre),caxis([-4 0]);colormap(jet),shading('flat')
        colorbar
        %grid
        % polarpeter([x x x x],ones(size(x))*[0.05 0.1 0.15 ],'w')
        polarpeter([x x x x x x x x],ones(size(x))*[.05:.05:.4],'w')
        polarpeter(pi/4*[1 1]'*[1:8] ,[0.0 0.4]'*ones(1,8),'w')
   
        c82 = cos(82*pi/180); s82 = sin(82*pi/180); rinc=0.00015;
        for ii=[0.04   0.1]
            text((ii+rinc/20)*c82,(ii+rinc/20)*s82, ...
                ['  ' num2str(ii)],'verticalalignment','bottom','handlevisibility','off')
        end
        title(['Per: ' num2str(1./frq(I(ic))) 's bin' num2str(ic) '      ''Timebin', num2str(iday) ' ' num2str(1+iday*2^9*20/60/60/24,3) ' A' num2str(theta(JJ)) ' S ' num2str(J),'day=',num2str(day)],'fontsize',14)
        pause(3)
%        drawnow,
        %if  (rem(isub,5)==0), drawnow, pause(1), end;
        %     title(['Freq: ' num2str(frq(ic),3) ' bin' num2str(ic)],'fontsize',14)
        %%
    end
end
end %day
