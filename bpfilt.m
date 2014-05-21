function [y]=bpfilt(x,dt,lf,hf);
%  BPFILT Bandpass filter time series. 
%         BPFILT(X,DT,LF,HF) takes a time series sampled at DT 
%	  and filters it with a 2nd order, 2 pass butterworth 
%         filter between frequencies LF and HF. If X is a matrix
%         BPFILT filters the individual rows of X.
%
%   nyq=0.5/dt;
%   wn=[lf/nyq,hf/nyq];
%   [b,a]=butter(2,wn);
%   [nx,mx]=(size(x)); 
%   y=zeros(size(x));
%   for ix=1:nx
%     y(ix,:)=filtfilt(b,a,x(ix,:));
%   end

x=a(:,2);
dt=0.02;
lf=0;
hf=0.5;

  nyq=0.5/dt;
  wn=[lf/nyq,hf/nyq];
  %cut=wn[1,2];
  [c,b]=butter(2,0.02);%0.02 should be wn
  [nx,mx]=(size(x)); 
  y=zeros(size(x)); 
  for ix=1:nx %have changed nx to mx
    y=filtfilt(c,b,x);
  end