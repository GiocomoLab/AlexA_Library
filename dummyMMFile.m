scanbox_config
istep=1;
ncol=796;
nlines=512;
if(sbconfig.mmap>0)     % make sure file exists...
    cprintf('*comment','[%02d] Setting up memory mapped files\n',istep); istep=istep+1;
    fnmm = which('scanbox');
    fnmm = strsplit(fnmm,'\');
    fnmm{end-1} = 'mmap';
    fnmm{end} = 'scanbox.mmap';
    sbconfig.fnmm = strjoin(fnmm,'\');    % name of memory mapped file
    
    if(~exist(sbconfig.fnmm,'file'))
        fidmm = fopen(sbconfig.fnmm,'w');
        fwrite(fidmm,zeros(1,16+1024*796*2*2,'int16'),'uint16'); % 16 words of header + max frame size * 2 channels
        fclose(fidmm);
    end
    
    % launch sbxplugin_server in separate matlab instance
    cprintf('*comment','[%02d] Launch plugin server\n',istep); istep=istep+1;
    system('matlab -nosplash -r sbxplugin_server &');
    
end


    try
        clear mmfile;
    catch
    end
    global mmfile;
    mmfile = memmapfile(sbconfig.fnmm,'Writable',true,'Format', ...
        { 'int16' [1 16] 'header' ; 'uint16' [nlines ncol] 'chA' ; 'uint16' [nlines ncol] 'chB'} , 'Repeat', 1);
    mmfile.Data.header(1) = -1;                 % semaphore or frame #
                                                % -1, not started
                                                % -2, stopped
                                                % 0...N frame number
    mmfile.Data.header(2) = int16(nlines);      % number of lines
    mmfile.Data.header(3) = int16(ncol);        % number of columns
    mmfile.Data.header(4) = 0;                  % TTL corresponding to stimulus
    mmfile.Data.header(5) = 0;   % volumetric scanning flag
    mmfile.Data.header(6) = int16(0);   % period of volumetric wave
    mmfile.Data.header(7) = int16(1); % code for plugin id #
    mmfile.Data.header(8) = int16(1); 
    mmfile.Data.header(9) = int16(1);
    mmfile.Data.header(10)= int16(0);
    mmfile.Data.header(11) = int16(0);
    mmfile.Data.header(12) = int16(0);



   % if grabbing or onfocus flag is on
   buffersCompleted=0;
   while true
       buffersCompleted=buffersCompleted+1;
            if(mmfile.Data.header(1)<0)                           % data was consumed?  If not, move on...  Server loses a frame
                mmfile.Data.header(4) = 1;               % ttl flag 
                mmfile.Data.chA = randi([0 100],nlines,ncol,'uint16');
               % mmfile.Data.chB = chB';
                mmfile.Data.header(1) = buffersCompleted;
                if mod(buffersCompleted,20)==3
             
                    mmfile.Data.header(10) = int16(randi([-100 100]));
                    mmfile.Data.header(11) = int16(randi([-100 100]));
                    mmfile.Data.header(12) = int16(randi([-100 100]));
                   
                end
            end
            pause(0.005);
   end
        