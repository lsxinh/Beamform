fid=fopen('trinet_resp','r')

%%
dum=fgetl(fid);
dum=fgetl(fid);
dum=fgetl(fid);
ista=0;
 while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if (strcmp(tline(1:2),'CI') & strcmp(tline(10:12),'LHZ'))
       if (str2num(tline(81:84))>=2006)
%          disp(tline)
          ista=ista+1;
   %       keyboard
          staresp(ista).sta=tline(4:6);
          staresp(ista).f0  =str2num(tline(21:29));
          staresp(ista).damp=str2num(tline(33:38));
       end
    end
 end
 
       fclose(fid);
%%
f0=[staresp(:).f0];
% 7+1= 8 stations: 
[val,I1]=find(f0==0.0027778);
[val,I4]=find(f0==0.002778);
% 114+7=121
[val,I2]=find(f0==0.0083333);
 [val,I5]=find(f0==0.0083);
 % 27
  [val,I3]=find(f0==0.033);
 J=sort([I1 I2 I3 I4 I5]);
 J1=sort([I2 I5]); % STS-2
 J2=sort([I1 I4]); %
 J3=sort(I3);   %
 %%
 nsta1=length(J1)
 nsta2=length(J2)
 nsta3=length(J3)
 ISTA1=[]; ISTA2=[]; ISTA3=[];
 for ista_ct=1:length(ISTA)
     ista=ISTA(ista_ct)
     for ista1=1:nsta1
     if strcmp(infom(ista).staname(1:3),staresp(J1(ista1)).sta)
         ISTA1=[ISTA1 ista];
       % keyboard
        break
      end
    end
    
     for ista1=1:nsta2
   if strcmp(infom(ista).staname(1:3),staresp(J2(ista1)).sta)
         ISTA2=[ISTA2 ista];
       % keyboard
        break
       end
    end
    
     for ista1=1:nsta3
     if strcmp(infom(ista).staname(1:3),staresp(J3(ista1)).sta)
         ISTA3=[ISTA3 ista];
       % keyboard
        break
     end
     end
 end

 ISTA=unique(sort([ISTA1 ISTA2 ISTA3]));
 nsta=length(ISTA);

fid =fopen('./stationresp/list1','r');
  dum=fgetl(fid); dum=fgetl(fid); dum=fgetl(fid);
cc1=fscanf(fid,'%f',[5 inf]) ;
size(cc1)
fid =fopen('./stationresp/list2','r');
 dum=fgetl(fid); dum=fgetl(fid); dum=fgetl(fid);
cc2=fscanf(fid,'%f',[5 inf]) ;
size(cc2)
fid =fopen('./stationresp/list3','r');
 dum=fgetl(fid); dum=fgetl(fid); dum=fgetl(fid);
cc3=fscanf(fid,'%f',[5 inf]);
size(cc3)

plot(cc1(1,2:end),cc1(5,2:end),cc2(1,2:end),cc2(5,2:end),cc3(1,2:end),cc3(5,2:end))
aaa=(cc1(2,:)+i*cc1(3,:));
stationresp(:,1)=aaa./abs(aaa);
aaa=(cc3(2,:)+i*cc1(3,:));
stationresp(:,2)=aaa./abs(aaa);
aaa=(cc2(2,:)+i*cc1(3,:));
stationresp(:,3)=aaa./abs(aaa);
