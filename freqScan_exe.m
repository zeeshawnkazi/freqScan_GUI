function freqScan_exe()
    % freqScan_exe() sets global parameters, cleans up stuff and executes the freqScan GUI.
  
    clc
    imaqreset
    delete(instrfind)
     
    if exist('hardware', 'var') == 1 % if necessary, kill most recent hardware object 
        hardware.kill();
        disp('killed old hardware object')
    end
    
    clear all
    close all
    
    % launch GUI
    freqScan();
end