function [vnvsig,vnvevi,zf] = voice_decision(s,fs,winlen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:  [vnv,vnvevidence,filtsig] = vnvseg(wav,fs, ntrend);
%
% Description: This function computes the voicing strength of the given speech
%	signal 'wav' with sampling rate 'fs'. 'nmean' specifies the window
%	size in 'ms' for removing local mean.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('winlen'))
    
 winlen = 10;
 
end;

s=s(:);
ds=diff(s);
ds(end+1)=ds(end);
%figure;plot(ds);
zf=zff_method(ds,fs);

%gci=find(diff((zf>0))==1); % +ve zero crossings
%es=abs(zf(gci+1)-zf(gci-1));


% s=s/sqrt(sum(s.^2)/length(s));
% ds=ds/sqrt(sum(ds.^2)/length(ds));
zf=zf/sqrt(sum(zf.^2)/length(zf));
zfe=mean_smooth(zf.^2,winlen*fs/1000);

%zfbys=zfe./se;
%zfbysevi=-tansig(zfbys-10);
zf_evi=1-exp(-10*zfe);

%rse10=res2sig(ds,fs,10);
%rse10=mean_smooth(rse10,20*fs/1000);
%rse10evi=-tansig(10*(rse10-.5));

vnvf_evi=zf_evi;
%vnvfevi=zfbys;
%vnvevi = [zfevi(:) zfbysevi(:) rse10evi(:)];
vnvevi = 10;	%my. just for argument sending
vnvsig=vnvf_evi > 0.4;
vnvsig=medfilt1(double(vnvsig),double(20*fs/1000));
%vnvsig=medfilt1(vnvsig,10*fs/1000);	%my
vnvsig=vnvsig > 0.0;

return;
