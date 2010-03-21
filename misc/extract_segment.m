function [segment,Fs]=extract_segment(dialect,speaker,phrase,ind_beg,ind_end)

% This function returns the speech segment of the TIMIT database specified
% by:   
%   - the dialect: 'dr1', 'dr2',...,'dr8'
%   - the speaker ID: 4 letters (1st onw is m or f) + 1 digit EXAMPLE:'fcjf0'
%   - the phrase text label: EXAMPLES 'sa1', 'sx127', 'sl1697'
%   - the samples to take, given by ind_beg and ind_end

%%%%% Locate TIMIT CD %%%%%
initial_dir=pwd;
[useless,timit_dir]=uigetfile('*.doc','Please select the "Read me" file in the Timit CD directory');
cd(initial_dir)

%%%%% Point to the right directory %%%%%
file_dir=strcat(timit_dir,'/timit/train/',dialect,'/',speaker);
file_name=strcat(phrase,'.wav');

%%%%% Read the wavefile and return desired segment %%%%% 
[segment,Fs]=TIMITread(file_dir,file_name);
segment=segment(ind_beg:ind_end);
