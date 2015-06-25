function [normzfSig, zfSig, gcisig, gci, es, esSig, fs, voiced_segments, winLength]=gci_location(wav,fs,vnv)
%	function [winLength]=epochExtract(wavFile)
% 	function [gci]=epochExtract(wav, fs)
% 	zfSig is the zero frequency signal derived from speech signal.
% 	gci is the glottal closure instant timing in seconds.

% 	Read the speech signal...
%	[wav fs]=wavread(wavFile);
	wav=resample(wav,8000,fs);
	fs=8000;
	disp('Sampling Frequency:');
	lwav=length(wav);

	vnv(1)=0;vnv(2)=0;		%exceptions
	vnv(end)=0;vnv(end-1)=0;	%exceptions

%	figure;plot([1:lwav]/fs,wav);hold on;plot([1:lwav]/fs,vnv);
	vnv_diff = diff(vnv);
	voiced_segments_start = find(vnv_diff==1);
	voiced_segments_end = find(vnv_diff==-1);

%	winLength = zeros();
% 	Find the average pitch period (in ms) for trend removal for every voiced segment...
	for i = 1:length(voiced_segments_start)

		voiced_segments(:,i)=[voiced_segments_start(i), voiced_segments_end(i)];

		wav_temp = wav(voiced_segments_start(i):voiced_segments_end(i));
		if( (voiced_segments_end(i)-voiced_segments_start(i)) > (5*fs/1000))
			winLength(i) = modifiedxcorrWinLen(wav_temp,fs);
		else
			winLength(i) = 2;
		end
	end
%	winLength

% Derive the zero-frequency filered signal...
	[zfSig normzfSig]=zeroFreqFilter(wav,fs,winLength,voiced_segments_start,voiced_segments_end,lwav);

	vad=getVad(zfSig,fs);

%%%%%%%%%%%% Detect the polarity of the signal and correct it %%%%%%%%%%%%%%%%
	
    [pzc pslope]=zerocros(zfSig,'p');
	[nzc nslope]=zerocros(zfSig,'n');
	
	l=min(length(pzc),length(nzc));

	dslope=abs(pslope(1:l))-abs(nslope(1:l));
	if(sum(dslope>0)>length(dslope)/2)
		gci=pzc;
		disp('Polarity of the signal: +ve');
	else
		gci=nzc;
		wav=-wav;
		zfSig=-zfSig;
        normzfSig = -normzfSig;
		disp('Polarity of the signal : -ve');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	gci=gci+3;
	es=-(normzfSig(gci-1)-normzfSig(gci+1))/2; %%% --- changing zfSig to normzfSig  - 08/10/12 - sathya %%%%%
    es = abs(es);
	%es = es/max(es);
	esSig = zeros(1,length(wav));
	esSig(gci) = es;


	gcisig=zeros(1,length(wav));
	gcisig(gci)=1;	

function [zf normzf]=zeroFreqFilter(wav,fs,winLength,voiced_segments_start,voiced_segments_end,lwav)

	% Difference the speech signal...
	dwav=diff(wav);
	dwav=dwav/max(abs(dwav));
	N=length(dwav);

	% Pass the differenced speech signal twice through zero-frequency resonator..	
% 	I commented all the plots

%	figure;plot(wav(1:voiced_segments_start(1)+40));
	if(voiced_segments_start(1) > 160)
		[normzf_temp, zf_temp,gci,es,wlen] = zfsig(wav(1:voiced_segments_start(1)+min(320,voiced_segments_start(1))),fs,5);
		zf(1:voiced_segments_start(1)) = zf_temp(1:voiced_segments_start(1));
		normzf(1:voiced_segments_start(1)) = normzf_temp(1:voiced_segments_start(1));
	end
 
%	zfSig(1:voiced_segments_start(1))=remTrend(zfSig(1:voiced_segments_start(1)),5);	
	for i = 1:length(voiced_segments_start)-1
		[normzf_temp, zf_temp,gci,es,wlen] = zfsig(wav(voiced_segments_start(i)-min(320,voiced_segments_start(i)-1) : voiced_segments_end(i)+min(320,lwav-voiced_segments_end(i))),fs,winLength(i));
		zf(voiced_segments_start(i):voiced_segments_end(i)) = zf_temp(min(321,voiced_segments_start(i)):end-min(320,lwav-voiced_segments_end(i)));
		normzf(voiced_segments_start(i):voiced_segments_end(i)) = normzf_temp(min(321,voiced_segments_start(i)):end-min(320,lwav-voiced_segments_end(i)));
	
		[normzf_temp, zf_temp,gci,es,wlen] = zfsig(wav(voiced_segments_end(i)-min(320,voiced_segments_end(i)-1):voiced_segments_start(i+1)+min(320,lwav-voiced_segments_start(i+1))),fs,6);
		zf(voiced_segments_end(i):voiced_segments_start(i+1)) = zf_temp(min(321,voiced_segments_end(i)):end-min(320,lwav-voiced_segments_start(i+1)));
		normzf(voiced_segments_end(i):voiced_segments_start(i+1)) = normzf_temp(min(321,voiced_segments_end(i)):end-min(320,lwav-voiced_segments_start(i+1)));

	end

	[normzf_temp, zf_temp,gci,es,wlen] = zfsig(wav(voiced_segments_start(end)-min(320,voiced_segments_start(end)) : voiced_segments_end(end)+min(320,N-voiced_segments_end(end))),fs,winLength(end));
	zf(voiced_segments_start(end):voiced_segments_end(end)) = zf_temp(321:end-min(320,N-voiced_segments_end(end)));
	normzf(voiced_segments_start(end):voiced_segments_end(end)) = normzf_temp(321:end-min(320,N-voiced_segments_end(end)));

	[normzf_temp, zf_temp,gci,es,wlen] = zfsig(wav(voiced_segments_end(end)-320:length(wav)),fs,5);
	zf(voiced_segments_end(end):length(wav)) = zf_temp(321:end);
	normzf(voiced_segments_end(end):length(wav)) = normzf_temp(321:end);

%%My function start
function [f]=peakpick(x)

	dx=diff(x);
	bdx=dx>0;
	dbdx=diff(bdx);
	f=find(dbdx==-1);
	
%%My function end
	
	
function [f,s]=zerocros(x,m)
	if nargin<2
    		m='b';
	end
	s=x>=0;
	k=s(2:end)-s(1:end-1);
	if any(m=='p')
  		f=find(k>0);
	elseif any(m=='n')
	    f=find(k<0);
	else
	    f=find(k~=0);
	end
	s=x(f+1)-x(f);

function [idx]=xcorrWinLen(wav,fs)

	frameSize=30*fs/1000;
	frameShift=20*fs/1000;

	en=conv(wav.^2,ones(frameSize,1));
	en=en(frameSize/2:end-frameSize/2);
	en=en/frameSize;
	en=sqrt(en);
	en=en>max(en)/5;

	b=buffer(wav,frameSize,frameShift,'nodelay');
	vad=sum(buffer(en,frameSize,frameShift,'nodelay'));

	FUN=@(x) xcorr((x-mean(x)).*hamming(length(x)),'coeff')./xcorr(hamming(length(x)),'coeff');
	out=blkproc(b,[frameSize,1],FUN);

	out=out(frameSize:end,:);
	
	minPitch=1;  %2 ms == 500 Hz.
       	maxPitch=16; %16 ms == 66.66 Hz.	

	[maxv maxi]=max(out(minPitch*fs/1000:maxPitch*fs/1000,:));

	%h=hist(maxi(vad>frameSize/2)+minPitch,(3:15)*8-4);
	x=(minPitch:0.5:maxPitch)*fs/1000+2;
	pLoc=maxi(vad>frameSize*0.8)+minPitch*fs/1000;
	y=hist(pLoc,x);
	y=y/length(pLoc);
	
	%bar(x,y,1,'EdgeColor',[1 1 1],'FaceColor',[0 0 0]);
	%set(gca,'xTick',(1:maxPitch)*fs/1000+0.5*fs/1000, 'xTickLabel',(1:maxPitch));
	%set(gca,'yTick',[0 0.1 0.2 0.3 0.4],'yTickLabel',[0 0.1 0.2 0.3 0.4]);
	%xlabel('Time (s)');
	%ylabel('Normalized frequency');
        %allText   = findall(h, 'type', 'text');
        %allAxes   = findall(h, 'type', 'axes');
        %allFont   = [allText; allAxes];
	%xlim([1 maxPitch+1]*fs/1000)
        %set(allFont,'FontSize',18);

	%advexpfig(h,'hist.eps','-deps2c','w',20,'h',20);

	%close(h);

	[val idx]=max(y);
	idx=round(idx/2)+minPitch+2;

function [idx]=modifiedxcorrWinLen(wav,fs)

	frameSize=30*fs/1000;
	frameShift=20*fs/1000;
	minPitch=1;  %2 ms == 500 Hz.
       	maxPitch=16; %16 ms == 66.66 Hz.

	if(length(wav) < 30*fs/1000)
		frameSize = 5*fs/1000;
		frameShift = 2*fs/1000;
		minPitch = 1;
		maxPitch = 5;
				
	end
	
	en=conv(wav.^2,ones(frameSize,1));
	en=en(frameSize/2:end-frameSize/2);
	en=en/frameSize;
	en=sqrt(en);
	en=en>max(en)/5;

	b=buffer(wav,frameSize,frameShift,'nodelay');
	vad=sum(buffer(en,frameSize,frameShift,'nodelay'));

	FUN=@(x) xcorr((x-mean(x)).*hamming(length(x)),'coeff')./xcorr(hamming(length(x)),'coeff');
	out=blkproc(b,[frameSize,1],FUN);

	out=out(frameSize:end,:);


	[maxv maxi]=max(out(minPitch*fs/1000:maxPitch*fs/1000,:));

	%h=hist(maxi(vad>frameSize/2)+minPitch,(3:15)*8-4);
	x=(minPitch:0.5:maxPitch)*fs/1000+2;
	pLoc=maxi(vad>frameSize*0.8)+minPitch*fs/1000;
	y=hist(pLoc,x);
	y=y/length(pLoc);


%	t=[1:length(y)];
%	subplot(211);plot(y);
%	y=y.*round(t/2);
%	subplot(212);plot(y);
%	pause;
	
	%bar(x,y,1,'EdgeColor',[1 1 1],'FaceColor',[0 0 0]);
	%set(gca,'xTick',(1:maxPitch)*fs/1000+0.5*fs/1000, 'xTickLabel',(1:maxPitch));
	%set(gca,'yTick',[0 0.1 0.2 0.3 0.4],'yTickLabel',[0 0.1 0.2 0.3 0.4]);
	%xlabel('Time (s)');
	%ylabel('Normalized frequency');
        %allText   = findall(h, 'type', 'text');
        %allAxes   = findall(h, 'type', 'axes');
        %allFont   = [allText; allAxes];
	%xlim([1 maxPitch+1]*fs/1000)
        %set(allFont,'FontSize',18);

	%advexpfig(h,'hist.eps','-deps2c','w',20,'h',20);

	%close(h);

	[val idx]=max(y);
	idx=round(idx/2)+minPitch;	

function [vad]=getVad(sig,fs)

	winLength=20*fs/1000;
	en=conv(abs(sig),ones(winLength,1));
	en=en(winLength/2:length(en)-winLength/2);
	en=en/max(en);

	%figure; plot(sig); hold on; plot(en/max(abs(en)),'r');

	vad=en>0.1;
