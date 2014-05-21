clear all
global EV
global U
global epsilon0

ev=dir(['/Volumes/JFNexternaldrive/Seismic/CASCADIAresponse/Jan/TestDay/Event*']); %list of Event directories
ddrr='/Volumes/JFNexternaldrive/Seismic/CASCADIAresponse/Jan/TestDay/'; %directory containing Event directories

for kk=1:1 %:length(ev); %loop over Events
    
    clearvars -except EV U epsilon0 ev ddrr kk

stnfile=[ddrr,ev(kk).name,'/processed1/stations.LHZ']

load(stnfile, '-mat')



%load -mat /Volumes/JFNexternaldrive/Seismic/CASCADIAresponse/Jan/Testday/Event_2012_266/processed1/stations.LHZ  %change this to pick up the date you want
%figure(1);clf;hold on; plot([infom(:).slon], [infom(:).slat],'r*');grid on
% FOR CAR LOcation is not in file
%if (strmatch(infom(154).staname(1:3),'CAR'))
%infom(154).slat=35.30819;
%infom(154).slon=	-119.84583;
%end
%%

% full grid
%Japan
LonLref=-180; LonUref= 180;  LatLref= 0; LatUref= 90;  
%LonLref=-125; LonUref= -112;  LatLref= 23; LatUref= 46;  
%
plot([LonLref LonLref LonUref LonUref LonLref],[LatLref LatUref LatUref ...
      LatLref LatLref]);   

stacoord=[ [infom.slon]; [infom.slat]];
%%Find the stations which belong to this grid 
ISTA=find( (stacoord(1,:)>=LonLref ) &  (stacoord(1,:)<=LonUref )& ...
   (stacoord(2,:)>=LatLref ) &  (stacoord(2,:)<=LatUref )  );
size(ISTA)
plot([infom(ISTA).slon], [infom(ISTA).slat],'g*')

meanlat=mean([infom(ISTA).slat]); meanlon=mean([infom(ISTA).slon]);

% readstaresp
 
ic=0;
for ista=ISTA
  ic=ic+1;
 [ran aa b]=dist_wh([meanlat [infom(ista).slat]],[meanlon [infom(ista).slon]]);
 if (aa<0),    aa=360+aa;   end  
 xsta(ic)=ran*sin(aa*pi/180); ysta(ic)=ran*cos(aa*pi/180);
end
coord=[xsta; ysta];
%plot(xsta,ysta,'*');figure; plot([infom(ISTA).slon], [infom(ISTA).slat],'g*')
%%


iyr=str2num(ev(kk).name(7:10))
for imonth=str2num(ev(kk).name(12:14)) %Change julian days here!!!!!!!!!!!!!!  0+[1:4] %:-1:6
%%359 
%% Change the julian date above
  iday1=1;
%Japan
%kcomp={'LHZ' 'BHN' 'BHE'};
%kcomp={'LHZ' 'BHN' 'BHE'};
inpath =   [ddrr,ev(kk).name,'/processed1/'] ; %change the directory to the sac2matfreq beam
%output directory
    %  sta1=stations{ista};
    for ista=ISTA
      sta1=infom(ista).staname(1:5);
      filename= [inpath sta1,'.',num2str(iyr),'.',  num2str(imonth) '.LHZ'];
      if  exist(filename,'file')
	eval(['load -mat ' filename]); %seis1(:,:,ic)=squeeze(fseis(1:end,iday1,:)); 
%	      eval(['load ' filename]);
%	      seis1(:,:,ic)=squeeze(fseis(1:end,:)); 
        break
    end
  end
 Nfreq=size(frq(Imin:Imax),2);
 Nsta=length(ISTA);
 Nframes= max(Nsample(imonth,:));

  seis1=zeros(Nfreq, Nframes,length(ISTA),'single');
    ic=0
    for ista=ISTA
      ic=ic+1;
      %  sta1=stations{ista};
      sta1=infom(ista).staname(1:5);
      j1=1;
      %    filename= [inpath sta1 kcomp{j1} num2str(234+iday) '.mat'];
       filename= [inpath sta1,'.',num2str(iyr),'.',  num2str(imonth) '.LHZ'];
      if  exist(filename,'file')
	eval(['load  -mat ' filename]); seis1(:,1:size(fseis,3),ic)=squeeze(fseis(1:end,iday1,:)); 
%	      eval(['load ' filename]); seis1(:,:,ic)=squeeze(fseis(1:end,:)); 
 size(fseis,3)
      else
  	display([filename ' did not exist' ])
      end
    end
 %   correct for instrument response
    % I=Imin:Imax;
%    for ista=ISTA
%      if sum(ista==ISTA1)
%         icc=1;
%       elseif sum(ista==ISTA2)%
%	 icc=2;
%       elseif sum(ista==ISTA3)%
%	 icc=3;
%       else%
%	 icc=0  %then in blow up!
%       end
%       for ii=1:Nfreq%
%	 seis1(ii,:,ista)=seis1(ii,:,ista)/stationresp(I(ii),icc);
%       end
%     end
%%
    I=Imin:Imax;
    theta = (0:2:360)';
    Nfrq= floor(length(I));
    Nsta=length(coord);
    Ntheta=length(theta);
    SL=0.0:0.01:0.4;
  %  SL=0.15:0.01:0.4;
    thetarad=theta*pi/180;
    projection=[sin(thetarad) cos(thetarad)]*coord;
    timestep=5
    Ntime=floor(Nframes/timestep);
   % beam=zeros(Ntheta,Nfrq,length(SL),Ntime,'single');
 %    beamab=zeros(Ntheta,Nfrq,length(SL),Ntime,'single');
  % beam=zeros(Ntheta,Nfrq,length(SL));
    icc=0;

 %Plane wave
per=6;
f=1/per;
omegaf=f*2*pi;
slowt=0.33;
ang=70;
angrad=ang*pi/180;
projection2=[sin(angrad) cos(angrad)]*coord;

    
    %keyboard
 for sl=SL   %1000:50:4000
     icc=icc+1
     disp(cputime)
     %for ifreq=1:Nfrq
         % this should be more correct 30 nov 8 am      freq=frq(I(ifreq)-1);
         %freq=frq(I(ifreq));
         %omega=freq*2*pi;
         rep=exp((j*omegaf*sl/1000)*projection)/sqrt(Nsta);
         DD=exp((j*omegaf*slowt/1000)*projection2)/sqrt(Nsta);
         beam(:,icc)= sum(abs(rep*DD').^2,2);
         %for itime=1:Ntime
            % itc=(itime-1)*timestep+1;
%change in here to do rotation 

             %vhelp= squeeze(seis1(ifreq,itc:itc+timestep,:));
             
             %beam(:,ifreq,icc,itime)= sum(abs(rep*vhelp').^2,2);
       %      C=vhelp'*vhelp;
             %           beam(:,ifreq,icc,itime)= sum(abs(rep*vhelp').^2,2);

             %clear EV
%             beamab(:,ifreq,icc,itime)= wncfast(C,rep(:,:)',2,0.2).';
         %end
   
     %end  %ifreq-v7.3
 end     %SL
%%   
    eval(['save  -v7.3 /Volumes/JFNexternaldrive/Seismic/CASCADIAresponse/Jan/TestDay/', ev(kk).name,'/processedbm_planewavetest/LHZ.',num2str(iyr),'.', num2str(imonth) ' beam  frq I  imonth timestep Ntime Nfrq SL theta Ntheta'])
%%
end

end
figure(1)
beamRt=10*log(beam./max(max(beam)));
 [XX,RAD]=meshgrid(SL,thetarad);
 [X1,Y1]=pol2cart(RAD,XX);
figure(per*100)
 h=polar([0 2*pi], [0 max(SL)]); axis('ij'), delete(h), hold on, view([-90 90]); 
  pcolor(X1,Y1,real(beamRt));shading flat;colorbar;caxis([-20 0])
