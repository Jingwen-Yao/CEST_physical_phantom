function [] = autoreg_tk(targfile, movfile, outfile, autoFlag, regfile)

if ~isfile(outfile)
    if autoFlag
        cmd = ['tkregister2 --targ ' targfile ' --mov ' movfile ' --reg ' regfile ' --regheader --noedit;'...
            'mri_vol2vol --targ ' targfile ' --mov ' movfile ' --reg ' regfile ' --o ' outfile ';'...
            'rm *.lta; rm *.reg'];
    else
        cmd = ['mri_vol2vol --targ ' targfile ' --mov ' movfile ' --reg ' regfile ' --o ' outfile ';'...
            'rm *.lta; rm *.reg'];
    end
    system(cmd);
end

end