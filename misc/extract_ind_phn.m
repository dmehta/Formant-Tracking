function [ind_beg,ind_end]=extract_ind_phn(file_dir,phrase,phoneme);

% This function returns the indices corresponding to the first phoneme
% 'phoneme' in the phrase 'phrase.phn', in the directory 'file_dir'

%%%%% Open the .phn file %%%%%
initial_dir=pwd;
phn_file=strcat(file_dir,'/',phrase,'.phn');
fp=fopen(phn_file,'r');

%%%%% Select the line in the file corresponding to the desired phoneme %%%%
current_phoneme='notaphoneme';
line_number=0;
while strcmp(current_phoneme,phoneme)==0
    X=fscanf(fp,'%s %s',2); 
    current_phoneme=fscanf(fp,'%s',1);
    line_number=line_number+1;
    if(strcmp(X,'')) % Fixes endless loop
        break;
    end
end
fclose(fp);

%%%%% Retrieve the phonemes indexes %%%%%
fp=fopen(phn_file,'r');
for ii=1:line_number-1
    X=fscanf(fp,'%s %s %s',3);
end
ind_beg=fscanf(fp,'%s',1);
ind_end=fscanf(fp,'%s',1);
fclose(fp);
