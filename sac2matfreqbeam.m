clear all

ev=dir(['/Volumes/JFNexternaldrive/Seismic/CASCADIADATA/Event*']); %list of Event directories
ddrr='/Volumes/JFNexternaldrive/Seismic/CASCADIADATA/'; %directory containing Event directories


%Nsample=zeros(366,36,18);
%stations={}

for kk=1:length(ev); %loop over Events
    stations={}
    clearvars -except ev ddrr kk 
    outpath =   ['/Volumes/JFNexternaldrive/Seismic/CASCADIADATA/',ev(kk).name,'/processed1/']; % this is the directory where the output is saved

    names=dir([ddrr,ev(kk).name,'/rsp/*LHZ.SAC.rsp']); %this is the date that you want to work with
    nfiles = length(names)  %number of stations

    Ndays=1;
    Ts=1;  %SACdata(1).times.delta
    Fe=1/Ts;

    %Create filter
    freq_int=[0.002 0.4];  
    [BB,AA]=butter(4,[freq_int]/Fe*2);


    ista=0
    ip=0

    for p=1:nfiles %loop over stations
    
        fprintf(1,'\r file= %i ',p)
    
        %read in sac file
        dataname =[ ddrr,ev(kk).name,'/rsp/', names(p).name] ;
        sac=readsac(dataname);
    
        %get info (station name, location, sampling rate)
        info.staname=sac.staname;
        info.kcomp=sac.kcomp;
        info.slat=sac.slat;
        info.slon=sac.slon;
        fs=round(1/sac.tau);
        info.fs=fs;
    
        %save info to "stations" and "infom" variables
        if (p==1 | isempty(strmatch(sac.staname(1:end),stations)))
            ista=ista+1;
           
            stations{ista}=sac.staname(1:end);
            infom(ista)=info;
            istacur=ista;
    
        else
            istacur=(strmatch(sac.staname(1:end),stations));
       
        end
    
        %get date info (year and day of year)
        imonth=sac.jr;
        iyr=sac.an;
    
        fprintf(1,'\r file %i istacur= %i imonth %i day-length %f ',p, istacur,imonth, sac.npts/3600/24)
    
    
        %get seismogram
        seis=sac.trace;
      
        %filter seismogram
        seis= filtfilt(BB,AA,detrend(seis));
    
        %calculate standard deviation
        sigmatrace=std(seis);

        %remove all values above one standard deviation
        seisz= seis;
        TEMP_FILTER='y';
        if TEMP_FILTER=='y'
            Threshold_STD=1;
            Threshold_Temp= Threshold_STD*sigmatrace; 
            II=find(abs(seisz(:))>Threshold_Temp);
            seisz(II)=Threshold_Temp*sign(seisz(II));
                
            % Removes mean
            meanfac=mean(seisz); seisz=seisz-ones(size(seisz))*meanfac;
        end;  % end temporal filtering
    
        seis=seisz;
     
        
         %spectrogram(seis,2^11,[],2^11,[]);
        
        
        %split the data into smaller chunks    
        p=11; %data split into 2048 second timeseries   
        frq=Fe/2*([1:2^(p-1)]-1)/2^(p-1); %set up frequency axis
        %frqtest=Fe/2*linspace(0,1,2^11/2+1);
        Ismall=2^p;  %length of each sub sample
        display(['Ismall' num2str(Ismall)])
        ipick=1:Ismall;
        Nsub=floor(size(seis,1)/2^p); %number of sub samples
        Nsample(imonth,istacur,iyr-1995)=Nsub;
        seissmall=zeros(2^p,Ndays,Nsub);
            
        for i=1:Ndays
            for jj=1:Nsub
                seissmall(:,i,jj)=seis((jj-1)*Ismall+ipick,i);
            end
        end
    
        %fft each sub sample
        fseis=fft(seissmall, 2^p); %plot(abs(fseis(:,1,1)))
    
        %Keep only phase of signal (remove amplitude)
        fseis(find(abs(fseis)==0))=1;%remove zeros
       % I=find(frq>0.0 & frq<3);
        I=find(frq>freq_int(1) & frq<freq_int(2));
        fseis=fseis(I,:,:)./abs(fseis(I,:,:));
            
        %define frequency min and max      
        frq2=frq(I); %positive frequency axis
        Imin=min(I);
        Imax=max(I);
     
        %save spectra
        outfile= strcat(sac.staname(1:end), '.', num2str(iyr), '.',  num2str(imonth), '.LHZ');
        save([outpath outfile], 'info', 'fseis', 'Imin', 'Imax', 'p', 'frq', 'frq2');
       
        clear sac
    end %for block (loop over stations)

    %save station info
    save([outpath 'stations.LHZ'],'infom', 'Nsample');

end %for block (loop over events)


