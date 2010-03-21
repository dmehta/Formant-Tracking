function [wave , SampleRate] = TIMITread(file_dir,file_name)

% TIMITread.m
% read TIMIT wave file
% [wave , SampleRate] = TIMITread(filename)
% Date : 1/1/2002
% Author : Guan-Ming Su
%...modified by M.S...

% point to right directory %
initial_dir=pwd;
cd(file_dir);

% get the extent file name
ind = find(file_name == '.');
ext = lower(file_name(ind+1:length(file_name))); 

% if the file ext is not TIMIT, return
if strcmp(ext,'wav')== 0; return; end
   
fp = fopen(file_name,'r');
if fp == -1 ; %error('Error on opening file'); 
    wave=[];
    SampleRate=[];
    return;
end;

fseek(fp,0,'bof');
nist=fscanf(fp,'%s',1);
if strcmp(nist(1:4),'NIST') ==0
   %error('Error reading voice file header, not TIMIT ..');
    wave=[];
    SampleRate=[];
   return;
end 
    
% number of bytes for file header
nbyte_header=fscanf(fp,'%d',1);

% skip database_id, database_version, utterance_id,channel_count,sample_count
for i=1:5; X=fscanf(fp,'%s %s %s',3);   end

% obtain sample rate
X=fscanf(fp,'%s %s',2);  SampleRate =fscanf(fp,'%d',1);

% skip sample_max and sample_min
for i=1:2;  X=fscanf(fp,'%s %s %s',3); end;

% obtain bytes per sample 	
X=fscanf(fp,'%s %s',2); bps=fscanf(fp,'%d',1);
if bps == 2;   ftype='short'; end  

% skip sample_byte_format
X=fscanf(fp,'%s %s %s',3);

% obtain bits per sample 	
X=fscanf(fp,'%s %s',2); bips=fscanf(fp,'%d',1);	

% skip the header		
fseek(fp,nbyte_header,'bof'); 

wave = fread(fp,inf,ftype);
wave=wave/2^(bips-1);
fclose(fp);

cd(initial_dir)