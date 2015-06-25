clc;
clear all;
close all;

data_dir = '/home/sriharsha/Desktop/scripts/scripts/';
fileid = '2.wav';

filename = sprintf('%s%s',data_dir,fileid);

[speech2,fs1] = wavread(filename);
speech2=speech2(:,1);
speech2 = speech2/abs(max(speech2)); %normalise the speech signal
speech1=resample(speech2,8000,fs1);  %resample the speech signal to 8000hz
lenspeech1=length(speech1);
fs=8000;

%s=resample(speech2,8000,fs1);
[vnv,y,z]=voice_decision(speech1,fs,10);
vnv = vnv';
% figure;
% plot(speech1,'k');
% hold on
% plot(0.9*vnv);
%[normzfSig, zfSig, gcis, gci, st, stSig, fs, voiced_regions, winLength]=gci_location(filename,vnv);
[normzfSig, zfSig, gcis, gci, st, stSig, fs, voiced_regions, winLength]=gci_location(speech2,fs1,vnv);		%for demo

% normst = (normzfSig(gci-1)-normzfSig(gci+1))/2;
%% normst = normst/max(normst);
% normstSig = zeros(1,lenspeech1);
% normstSig(gci) = normst;
%
% st = abs(normst);	%check

%%%%%%%%%%%%%%% Features %%%%%%%%%%%%%%%%%%%%

pitchperiods = diff(gci);  %% - Pitch periods

pitchperiods=[pitchperiods pitchperiods(end)];
ppdiff = diff(pitchperiods);   %%% -- Difference in pitch periods
ppdiff(end+1) = ppdiff(end);

parr=zeros(1,lenspeech1); %% -- Pitch period array
parr(gci)=pitchperiods;
parr=parr.*vnv(1:end);  %considering only voiced pitch periods

parr_ms=parr*1000/fs;   %pitch periods in milli seconds

starr=zeros(1,lenspeech1); %% -- SoE array
starr(gci)=st;
starr=starr.*vnv(1:end);

arr=zeros(1,lenspeech1); %% -- vnv array
arr(gci)=1;
arr=arr.*vnv(1:end);

ratioofStrandPit = st./pitchperiods;
ratioofStrandPitarr = zeros(1,lenspeech1);
ratioofStrandPitarr(gci) = ratioofStrandPit;
ratioofStrandPitarr = ratioofStrandPitarr.*vnv(1:end);

pitchwindow = 5;
pitchoverlap = 4;
pitchperiodsbuffer = buffer(pitchperiods,pitchwindow,pitchoverlap,'nodelay'); %% -- buffering of pitch periods
ppdiffbuffer = buffer(ppdiff,pitchwindow,pitchoverlap,'nodelay'); %% -- buffering of diff. in pitch periods

dur_buffer = sum(pitchperiodsbuffer(1:pitchoverlap,:)); %% summing the pitch periods in a buffer to find the duration of the segment
%pitchperiodsdiff = (max(pitchperiodsbuffer)-min(pitchperiodsbuffer))./max(pitchperiodsbuffer);  %%% finding the difference in pitch periods by taking diff. of max and min pitch period in a buffer and normalizing it
pitchperiodsdiff = (max(pitchperiodsbuffer)-min(pitchperiodsbuffer));%./max(pitchperiodsbuffer);  %%% finding the difference in pitch periods by taking diff. of max and min pitch period in a buffer and normalizing it
%pitchperiodsdiff = (max(pitchperiodsbuffer)-min(pitchperiodsbuffer))/(fs/1000); %%% finding the difference in pitch periods by taking diff. of max and min pitch period in a buffer and normalizing it
pitchperiodsdiff(end:end+pitchoverlap) = pitchperiodsdiff(end);
pitchperiodsdiffarr=zeros(1,lenspeech1);
pitchperiodsdiffarr(gci)=pitchperiodsdiff;
pitchperiodsdiffarr=pitchperiodsdiffarr.*vnv(1:end);   %%%% pitch period difference array

ppcumdiff = (sum(ppdiffbuffer))./(dur_buffer.*(max(pitchperiodsbuffer)));  %%% basically trying to average pitch period difference in each buffer -- could be thought of as smoothing
ppcumdiff(end:end+pitchoverlap) = ppcumdiff(end);
ppcumdiffarr=zeros(1,lenspeech1);
ppcumdiffarr(gci)=ppcumdiff;
ppcumdiffarr=ppcumdiffarr.*vnv(1:end); %%%% not being used currently


stdiff = diff(st);
stdiff(end+1) = stdiff(end);
stbuffer = buffer(st,pitchwindow,pitchoverlap,'nodelay');  %%% buffering of SoE
stdiffbuffer = buffer(stdiff,pitchwindow,pitchoverlap,'nodelay');  %%% buffering of diff. in SoE
%stdiff = (max(stbuffer)-min(stbuffer))./max(stbuffer);  %%% finding the diff. in SoE by taking the diff. of max and min SoE in each buffer and normalizing it
stdiff = (max(stbuffer)-min(stbuffer));%./max(stbuffer);  %%% finding the diff. in SoE by taking the diff. of max and min SoE in each buffer
stdiff(end:end+pitchoverlap) = stdiff(end);
stdiffarr=zeros(1,lenspeech1);
stdiffarr(gci)=stdiff;
stdiffarr=stdiffarr.*vnv(1:end);

stdiff1 = (max(stbuffer)-min(stbuffer))./(dur_buffer.*(max(stbuffer)));  %%% basically trying to average SoE difference in each buffer -- could be thought of as smoothing
stdiff1(end:end+pitchoverlap) = stdiff(end);
stdiffarr1=zeros(1,lenspeech1);
stdiffarr1(gci)=stdiff1;
stdiffarr1=stdiffarr1.*vnv(1:end);  %%%%

stcumdiff = (sum(stdiffbuffer))./(dur_buffer.*(max(stbuffer)));
stcumdiff(end:end+pitchoverlap) = stcumdiff(end);
stcumdiffarr=zeros(1,lenspeech1);
stcumdiffarr(gci)=stcumdiff;
stcumdiffarr=stcumdiffarr.*vnv(1:end); %%%%
%%%%%%%%%%%%%%% Features %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Algorithm Start %%%%%%%%%%%%%

%%%%% Setting thresholds for various features %%%%

pitchperiodsThreshold = (7*fs)/1000;	% less 40 percent
stThreshold = 0.1;		% more 20 percent  %%% changed from 0.1 to 0.09
ratioofStrandPitThreshold = 0.0016;	% more 10 percent  %%%%% changed from 0.002 to 0.0016 0n 25/6/15
%pitchperiodsdiffThreshold = 0.001*160;	% more 20 percent
%pitchperiodsdiffThreshold = 0.2;	% more 20 percent
pitchperiodsdiffThreshold = (2*fs)/1000;	% more 20 percent %%%
%pitchperiodsdiffThreshold = 0.001;	% more 70 percent
%stdiffThreshold = 0.004*160;	% more 20 percent
stdiffThreshold = 0.08;	% more 20 percent %%%  note: chnaged from 0.1 to 0.08 on 25/6/15


laugh_regions = [];
% laughter_regions = [];
laughter_regions_pp = [];
laughter_regions_st = [];
laughter_regions_ratio = [];
laughter_regions_ppdiff = [];
laughter_regions_stdiff = [];
laughter_regions_ppslope = [];
laughter_regions_stslope = [];
y_lms = [];
y_polyfit = [];
x_tot = [];
% laughter_regions_num = [];		%voiced segment number (in order) which is laughter
laugh_region_num = [];

%%%% Extracting all features for each voiced region and finding if that
%%%% voiced region is a laughter region using that feature
%%%% We are also keeping a count of how many features have approved of the
%%%% region as a laugh region.

[rows_voiced_regions columns_voiced_regions] = size(voiced_regions);


for i=1:columns_voiced_regions
    
    temp_pitchperiods =  pitchperiods(find( (gci>voiced_regions(1,i)) & (gci<voiced_regions(2,i) )));
    temp_st =  st(find( (gci>voiced_regions(1,i)) & (gci<voiced_regions(2,i) )));
    temp_ratio =  ratioofStrandPit(find( (gci>voiced_regions(1,i)) & (gci<voiced_regions(2,i) )));
    temp_ppdiff =  pitchperiodsdiff(find( (gci>voiced_regions(1,i)) & (gci<voiced_regions(2,i) )));
    temp_stdiff =  stdiff(find( (gci>voiced_regions(1,i)) & (gci<voiced_regions(2,i) )));
    
    laughcount=0;
    
    %%% If pitch period values exceed the pitch period threshold for more than
    %%% 40% of the region, then consider it as laugh
    
    if( length(find(temp_pitchperiods<pitchperiodsThreshold)) > 0.5*length(temp_pitchperiods) )	%pitch period constraint
        laughter_regions_pp = [laughter_regions_pp voiced_regions(:,i)];
        laughcount=laughcount+1;
    end
    
    %%% If SoE values exceed the SoE threshold for more than
    %%% 5% of the region, then consider it as laugh
    flag_st = 0;
    if( length(find(temp_st>stThreshold)) > 0.2*length(temp_st) )	%strength constraint  %% changed from 0.3 to 0.2  25/6/15
        laughter_regions_st = [laughter_regions_st voiced_regions(:,i)];
        laughcount=laughcount+1;
        flag_st = 1;   %% flag_st is made '1' if it satisfies the threshold criterion of SoE
    end
    %%% If maximum SoE value in a segment is higher than a
    %%% threshold, then give it a bias that it might be a laugh segment by
    %%% adding the count by 0.5
    if flag_st == 0
        if max(temp_st) > 0.16  %%%%% changed from 0.2 to 0.16
            laughcount = laughcount+0.5;
        end
    end
    
    %%% If ratio of SoE and pitch period values exceed the ratio threshold for more than
    %%% 10% of the region, then consider it as laugh
    
    if( length(find(temp_ratio>ratioofStrandPitThreshold)) > 0.25*length(temp_ratio) )	%ratio constraint %%% changed from 0.3 to 0.25 on 25/6/15
        laughter_regions_ratio = [laughter_regions_ratio voiced_regions(:,i)];
        laughcount=laughcount+1;
    end
    
    %%% If 'differences in pitch period' values exceed the 'differences in pitch period' threshold for more than
    %%% 40% of the region, then consider it as laugh
    
    if( length(find(temp_ppdiff>pitchperiodsdiffThreshold)) > 0.08*length(temp_ppdiff) )	%ppdiff constraint %%% chnaged from 0.4 to 0.2 25/6/15
        laughter_regions_ppdiff = [laughter_regions_ppdiff voiced_regions(:,i)];
        laughcount=laughcount+1;
        %		laughter_regions = [laughter_regions voiced_regions(:,i)];
    end
    
    %%% If 'differences in SoE' values exceed the 'differences in SoE' threshold for more than
    %%% 40% of the region, then consider it as laugh
    
    if( length(find(temp_stdiff>stdiffThreshold)) > 0.15*length(temp_stdiff) )	%strengthdiff constraint  %%% changed from 0.2 to 0.15 on 25/6/15
        laughter_regions_stdiff = [laughter_regions_stdiff voiced_regions(:,i)];
        laughcount=laughcount+0.5;
    end
    
    %% -- if the condition (threshold) is satisfied for more than 3 out of 5 features, then we consider it as a laugh region %%
    %disp('length of voiced regions');
    length_voiced_regions = voiced_regions(2,i)-voiced_regions(1,i);
    if laughcount>=3
        if length_voiced_regions > (30*fs/1000)   %%%% Removing spurious small regions lesser than 30 ms. in duration under the assumption that a call of such small duration rarely exists
            laugh_regions = [laugh_regions voiced_regions(:,i)];
            laugh_region_num = [laugh_region_num i];
        end
    end
    
end

%%% Laugh regions in my decision method  %%%%

laugh_regions_sig = zeros(1,lenspeech1);
if(length(laugh_regions) ~=0)
    for i = 1:length(laugh_regions(1,:))
        laugh_regions_sig(laugh_regions(1,i):laugh_regions(2,i)) = 1;
    end
end

%%%%% Finding the intersections of laughter regions from different features
%%%%% using intersection method

laughter_regions_temp1 = intersect(laughter_regions_pp',laughter_regions_ppdiff','rows')';
laughter_regions_temp2 = intersect(laughter_regions_temp1', laughter_regions_stdiff','rows')';
laughter_regions_temp3 = intersect(laughter_regions_temp2',laughter_regions_st','rows');
laughter_regions = intersect(laughter_regions_temp3,laughter_regions_ratio','rows')';
[temp, ia, ib] = intersect(laughter_regions',voiced_regions','rows');
laughter_regions_num = ib;

%%%%%%%%%%%%%%%%  regions detected based on different features  %%%%%%%%%%%

laughter_regions_sig = zeros(1,lenspeech1);
if(length(laughter_regions) ~=0)
    for i = 1:length(laughter_regions(1,:))
        laughter_regions_sig(laughter_regions(1,i):laughter_regions(2,i)) = 1;
    end
end

laughter_regions_pp_sig = zeros(1,lenspeech1);				%pp decision signal
if(length(laughter_regions_pp) ~=0)
    for i = 1:length(laughter_regions_pp(1,:))
        laughter_regions_pp_sig(laughter_regions_pp(1,i):laughter_regions_pp(2,i)) = 1;
    end
end

laughter_regions_st_sig = zeros(1,lenspeech1);				%st decision signal
if(length(laughter_regions_st) ~=0)
    for i = 1:length(laughter_regions_st(1,:))
        laughter_regions_st_sig(laughter_regions_st(1,i):laughter_regions_st(2,i)) = 1;
    end
end

laughter_regions_ratio_sig = zeros(1,lenspeech1);			%ratio decision signal
if(length(laughter_regions_ratio) ~=0)
    for i = 1:length(laughter_regions_ratio(1,:))
        laughter_regions_ratio_sig(laughter_regions_ratio(1,i):laughter_regions_ratio(2,i)) = 1;
    end
end

laughter_regions_ppdiff_sig = zeros(1,lenspeech1);			%ppdiff decision signal
if(length(laughter_regions_ppdiff) ~=0)
    for i = 1:length(laughter_regions_ppdiff(1,:))
        laughter_regions_ppdiff_sig(laughter_regions_ppdiff(1,i):laughter_regions_ppdiff(2,i)) = 1;
    end
end

laughter_regions_stdiff_sig = zeros(1,lenspeech1);			%stdiff decision signal
if(length(laughter_regions_stdiff) ~=0)
    for i = 1:length(laughter_regions_stdiff(1,:))
        laughter_regions_stdiff_sig(laughter_regions_stdiff(1,i):laughter_regions_stdiff(2,i)) = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Post-processing  %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%  Code snippet to remove very small segments  %%%%%%%%%%%%%%

laughter_regions_mod = laughter_regions_sig;
laughter_regions_mod(1) = 0;
laughter_regions_mod(end) = 0;

laughter_seg = zeros(lenspeech1,1);
laughter_seg(2:end) = diff(laughter_regions_mod);
laughter_seg_st = find(laughter_seg == 1);
laughter_seg_end = find(laughter_seg == -1);
no_seg = length(laughter_seg_st);

if no_seg >= 1
    for k = 1 : no_seg
        
        if laughter_seg_end(k) - laughter_seg_st(k) <= 200 %%% 200 samples = 25 msec
            
            laughter_regions_mod(laughter_seg_st(k) : laughter_seg_end(k)) = 0;
            
        end
        
    end
    
end

laughter_regions_sig = laughter_regions_mod;

%%%%%%%%%  Code snippet to include short missed alughter segments  %%%%%%%%

laughter_regions_mod1 = laughter_regions_sig;
laughter_regions_mod1(1) = 0;
laughter_regions_mod1(end) = 0;

laughter_seg = zeros(lenspeech1,1);
laughter_seg(2:end) = diff(laughter_regions_mod1);
laughter_seg_st = find(laughter_seg == 1);
laughter_seg_end = find(laughter_seg == -1);
no_seg = length(laughter_seg_st);

voiced_reg = voiced_regions';
voi_reg_st = voiced_reg(:,1);
voi_reg_end = voiced_reg(:,2);

laughter_pp_seg = zeros(lenspeech1,1);
laughter_regions_pp_sig(1) = 0;
laughter_regions_pp_sig(end) = 0;
laughter_pp_seg(2:end) = diff(laughter_regions_pp_sig);
laughter_pp_st = find(laughter_pp_seg == 1);
laughter_pp_end = find(laughter_pp_seg == -1);

laughter_st_seg = zeros(lenspeech1,1);
laughter_regions_st_sig(1) = 0;
laughter_regions_st_sig(end) = 0;
laughter_st_seg(2:end) = diff(laughter_regions_st_sig);
laughter_st_st = find(laughter_st_seg == 1);
laughter_st_end = find(laughter_st_seg == -1);

laughter_ratio_seg = zeros(lenspeech1,1);
laughter_regions_ratio_sig(1) = 0;
laughter_regions_ratio_sig(end) = 0;
laughter_ratio_seg(2:end) = diff(laughter_regions_ratio_sig);
laughter_ratio_st = find(laughter_ratio_seg == 1);
laughter_ratio_end = find(laughter_ratio_seg == -1);

laughter_ppdiff_seg = zeros(lenspeech1,1);
laughter_regions_ppdiff_sig(1) = 0;
laughter_regions_ppdiff_sig(end) = 0;
laughter_ppdiff_seg(2:end) = diff(laughter_regions_ppdiff_sig);
laughter_ppdiff_st = find(laughter_ppdiff_seg == 1);
laughter_ppdiff_end = find(laughter_ppdiff_seg == -1);

laughter_stdiff_seg = zeros(lenspeech1,1);
laughter_regions_stdiff_sig(1) = 0;
laughter_regions_stdiff_sig(end) = 0;
laughter_stdiff_seg(2:end) = diff(laughter_regions_stdiff_sig);
laughter_stdiff_st = find(laughter_stdiff_seg == 1);
laughter_stdiff_end = find(laughter_stdiff_seg == -1);

if no_seg >=2
    for k = 2 : no_seg
        
        if laughter_seg_st(k) - laughter_seg_end(k-1) <= 6800 %%% 4400 samples = 850 msec
            
            miss_reg1 = find(voi_reg_st <= laughter_seg_st(k) & voi_reg_st >= laughter_seg_end(k-1));
            lmr = length(miss_reg1);
            
            if lmr >= 1
                
                for mr = 1 : lmr
                    
                    len_req_pp = length(find(laughter_pp_st == voi_reg_st(miss_reg1(mr))));
                    len_req_st = length(find(laughter_st_st == voi_reg_st(miss_reg1(mr))));
                    len_req_ratio = length(find(laughter_ratio_st == voi_reg_st(miss_reg1(mr))));
                    len_req_ppdiff = length(find(laughter_ppdiff_st == voi_reg_st(miss_reg1(mr))));
                    len_req_stdiff = length(find(laughter_stdiff_st == voi_reg_st(miss_reg1(mr))));
                    
                    tot_len_feat = len_req_pp + len_req_st + len_req_ratio + len_req_ppdiff + len_req_stdiff;
                    
                    if (voi_reg_end(miss_reg1(mr)) - voi_reg_st(miss_reg1(mr)) <= 1200) && (tot_len_feat >= 3) %%% 1200 samples = 125 msec
                        
                        laughter_regions_mod1(voi_reg_st(miss_reg1(mr)) : voi_reg_end(miss_reg1(mr))) = 1;
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

%%%%%%%%%%%%%  Code snippet to join detected laughter segments  %%%%%%%%%%%

laughter_regions_sig = laughter_regions_mod1;
laughter_regions_mod2 = laughter_regions_sig;
laughter_regions_mod2(1) = 0;
laughter_regions_mod2(end) = 0;

laughter_seg = zeros(lenspeech1,1);
laughter_seg(2:end) = diff(laughter_regions_mod2);
laughter_seg_st = find(laughter_seg == 1);
laughter_seg_end = find(laughter_seg == -1);
no_seg = length(laughter_seg_st);

if no_seg >=2
    
    for k = 2 : no_seg
        
        if laughter_seg_st(k) - laughter_seg_end(k-1) <= 1600 %%% 1600 samples = 200 msec
            
            laughter_regions_mod2(laughter_seg_end(k-1) : laughter_seg_st(k)) = 1;
            
        end
        
    end
    
end

laughter_regions_sig = laughter_regions_mod2;

%%%%%%%%%%%  Code to remove short segments detected as laughter  %%%%%%%%%%

laughter_regions_mod3 = laughter_regions_sig;
laughter_regions_mod3(1) = 0;
laughter_regions_mod3(end) = 0;

laughter_seg = zeros(lenspeech1,1);
laughter_seg(2:end) = diff(laughter_regions_mod3);
laughter_seg_st = find(laughter_seg == 1);
laughter_seg_end = find(laughter_seg == -1);
no_seg = length(laughter_seg_st);

if no_seg >= 1
    for k = 1 : no_seg
        
        if laughter_seg_end(k) - laughter_seg_st(k) <= 640 %%% 640 samples = 80 msec
            
            laughter_regions_mod3(laughter_seg_st(k) : laughter_seg_end(k)) = 0;
            
        end
        
    end
    
end

laughter_regions_sig = laughter_regions_mod3;

%%%%%%%%%%%%%%%  Code to write laughter boundaries top a file  %%%%%%%%%%%%

laughter_bound = laughter_regions_sig;
laughter_bound(1) = 0;
laughter_bound(end) = 0;

laughter_seg = zeros(lenspeech1,1);
laughter_seg(2:end) = diff(laughter_bound);
laughter_seg_st = find(laughter_seg == 1);
laughter_seg_end = find(laughter_seg == -1);

laughter_st_time = laughter_seg_st./fs;
laughter_end_time = laughter_seg_end./fs;
no_seg = length(laughter_st_time);

non_laughter_end_time = zeros(no_seg + 1,1);
non_laughter_st_time = zeros(no_seg + 1,1);

if no_seg == 0
    
    non_laughter_end_time = lenspeech1/fs;
    non_laughter_st_time = 0;
    
else
    
    non_laughter_end_time(1:end - 1) = (laughter_seg_st)./fs;
    non_laughter_end_time(end) = lenspeech1/fs;
    non_laughter_st_time(2:end) = (laughter_seg_end)./fs;
    
end

fname1 = fileid(1:end-4);
fname = sprintf('%s%s',fname1,'.txt');
fid = fopen(fname,'w');

if no_seg >= 1
    
    for k = 1 : no_seg
        
        fprintf(fid,'%f\t%f\t%s\n',non_laughter_st_time(k),non_laughter_end_time(k),'non laughter');
        fprintf(fid,'%f\t%f\t%s\n',laughter_st_time(k),laughter_end_time(k),'laughter');
        
    end

    fprintf(fid,'%f\t%f\t%s',non_laughter_st_time(k+1),non_laughter_end_time(k+1),'non laughter');
    
else
    
    fprintf(fid,'%f\t%f\t%s',non_laughter_st_time(1),non_laughter_end_time(1),'non laughter');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Post-Processing End  %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Algorithm End  %%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Plotting for individual features and all the features combined %%%%%

%%% Assigning values for the sake of plotting %%%

laughter_regions_pp_sig = zeros(1,lenspeech1);				%pp decision signal
if(length(laughter_regions_pp) ~=0)
	for i = 1:length(laughter_regions_pp(1,:))
		laughter_regions_pp_sig(laughter_regions_pp(1,i):laughter_regions_pp(2,i)) = 10;
	end
end

laughter_regions_st_sig = zeros(1,lenspeech1);				%st decision signal
if(length(laughter_regions_st) ~=0)
	for i = 1:length(laughter_regions_st(1,:))
		laughter_regions_st_sig(laughter_regions_st(1,i):laughter_regions_st(2,i)) = max(starr);
	end
end

laughter_regions_ratio_sig = zeros(1,lenspeech1);			%ratio decision signal
if(length(laughter_regions_ratio) ~=0)
	for i = 1:length(laughter_regions_ratio(1,:))
		laughter_regions_ratio_sig(laughter_regions_ratio(1,i):laughter_regions_ratio(2,i)) = 0.015;
	end
end

laughter_regions_ppdiff_sig = zeros(1,lenspeech1);			%ppdiff decision signal
if(length(laughter_regions_ppdiff) ~=0)
	for i = 1:length(laughter_regions_ppdiff(1,:))
		laughter_regions_ppdiff_sig(laughter_regions_ppdiff(1,i):laughter_regions_ppdiff(2,i)) = max(pitchperiodsdiffarr)*0.9;
	end
end

laughter_regions_stdiff_sig = zeros(1,lenspeech1);			%stdiff decision signal
if(length(laughter_regions_stdiff) ~=0)
	for i = 1:length(laughter_regions_stdiff(1,:))
		laughter_regions_stdiff_sig(laughter_regions_stdiff(1,i):laughter_regions_stdiff(2,i)) = max(stdiffarr)*0.9;
	end
end

%x_alphas = (lenspeech1/fs)*(1.025);
x_alphas = 5.4*(1.025);

h=figure;
ax(1)=subplot(611);plot([1:lenspeech1]/fs,speech1/max(abs(speech1)),'k');%hold on;stem(find(arr>0)/fs,arr(find(arr>0)),'k');
ylim([-1 1]);grid;text(x_alphas,0,'(a)');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
hold on;plot((1:lenspeech1)/fs,0.6 * laugh_regions_sig,'r');
%hold on;plot([1:lenspeech1]/fs,actual_laugharr,'--k');
ax(2)=subplot(612);plot(find(parr~=0)/fs,parr(find(parr~=0))/(fs/1000),'k.');grid;ylim([0 15]);ylabel('T_{0} (ms)');text(x_alphas,7.5,'(b)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_pp_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
ax(3)=subplot(613);plot(find(starr~=0)/fs,starr(find(starr~=0)),'k.');grid;text(x_alphas,0.2,'(c)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_st_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
%ax(4)=subplot(614);plot(find(ratioofStrandPit~=0)/fs,ratioofStrandPit(find(ratioofStrandPit~=0)),'k.');grid;ylabel('');text(x_alphas,0.01,'(d)');
ax(4)=subplot(614);plot(find(ratioofStrandPitarr~=0)/fs,ratioofStrandPitarr(find(ratioofStrandPitarr~=0)),'k.');grid;text(x_alphas,0.01,'(d)');ylim([0 0.02]);
hold on;plot([1:lenspeech1]/fs,laughter_regions_ratio_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
ax(5)=subplot(615);plot(find(pitchperiodsdiffarr~=0)/fs,pitchperiodsdiffarr(find(pitchperiodsdiffarr~=0)),'k.');grid;text(x_alphas,0.0075,'(e)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_ppdiff_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
ax(6)=subplot(616);plot(find(stdiffarr~=0)/fs,stdiffarr(find(stdiffarr~=0)),'k.');grid;text(x_alphas,0.01,'(f)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_stdiff_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
linkaxes(ax,'x');xlabel('Time(s)');
xlim([1 lenspeech1]/fs);
%xlim([3.5 5.5]);

h=figure;
ax(1)=subplot(611);plot([1:lenspeech1]/fs,speech1/max(abs(speech1)),'k');%hold on;stem(find(arr>0)/fs,arr(find(arr>0)),'k');
ylim([-1 1]);grid;text(x_alphas,0,'(a)');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
hold on;plot((1:lenspeech1)/fs,0.8 * laughter_regions_sig,'r');
%hold on;plot([1:lenspeech1]/fs,actual_laugharr,'--k');
ax(2)=subplot(612);plot(find(parr~=0)/fs,parr(find(parr~=0))/(fs/1000),'k.');grid;ylim([0 15]);ylabel('T_{0} (ms)');text(x_alphas,7.5,'(b)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_pp_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
ax(3)=subplot(613);plot(find(starr~=0)/fs,starr(find(starr~=0)),'k.');grid;text(x_alphas,0.2,'(c)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_st_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
%ax(4)=subplot(614);plot(find(ratioofStrandPit~=0)/fs,ratioofStrandPit(find(ratioofStrandPit~=0)),'k.');grid;ylabel('');text(x_alphas,0.01,'(d)');
ax(4)=subplot(614);plot(find(ratioofStrandPitarr~=0)/fs,ratioofStrandPitarr(find(ratioofStrandPitarr~=0)),'k.');grid;text(x_alphas,0.01,'(d)');ylim([0 0.02]);
hold on;plot([1:lenspeech1]/fs,laughter_regions_ratio_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
ax(5)=subplot(615);plot(find(pitchperiodsdiffarr~=0)/fs,pitchperiodsdiffarr(find(pitchperiodsdiffarr~=0)),'k.');grid;text(x_alphas,0.0075,'(e)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_ppdiff_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
ax(6)=subplot(616);plot(find(stdiffarr~=0)/fs,stdiffarr(find(stdiffarr~=0)),'k.');grid;text(x_alphas,0.01,'(f)');
hold on;plot([1:lenspeech1]/fs,laughter_regions_stdiff_sig,'r');
%set(gca,'xTick',[3.5:0.5:5.5],'xTickLabel',[0:0.5:2.0]);
linkaxes(ax,'x');xlabel('Time(s)');
xlim([1 lenspeech1]/fs);