function [out] = pulseq_scp_skyra(pulseq_user, external_path)

    temp_file1 = strrep(which('pulseq_scp_hydra'),'\','/');
    temp_file1 = fopen([temp_file1(1:end-19), '/hydra1.txt']);
    temp_file2 = strrep(which('pulseq_scp_hydra'),'\','/');
    temp_file2 = fopen([temp_file2(1:end-19), '/hydra2.txt']);
    temp_cmd  = ['pscp -pw "' char(str2num(fgetl(temp_file1))-ceil(exp(pi))) '" -sftp ' external_path ' ' pulseq_user fgetl(temp_file2) pulseq_user '/'];
    system(temp_cmd);
    fclose(temp_file1);
    fclose(temp_file2);
    clear temp_file1 temp_file2 temp_line1 temp_line2 temp_cmd;
    out = 0;

end