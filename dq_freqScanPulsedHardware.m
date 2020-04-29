classdef dq_freqScanPulsedHardware < handle
    % freqScanHardware is an object that encapsulates hardware interactions
    % including:
    %   -WINDFREAK SynthHD dual channel RF generator
    %   -Image acquisition
    %    %%% NOT YET %%% -R&H RF generator
    % 
    % Only one Hardware object needs to be created for a GUI; simply kill()
    % and redo the init() call to re-program the hardware.
    
    
    
    properties
        rf_sweeper                      % rf_sweeper object (through DAQ)
        dq_rf_sweeper                   % Dual Channel rf_sweeper object (through USB)
        vid                             % camera object
        is_initialized                  % boolean flag; '1' if this is already initialized
        exposure_time_sec               % exposure time of camera image in seconds
        rf_frequency                    % rf_frequency
        rf_power                        % RF power in dBm
        ccd_size                        % image size
        capture_taken                   % boolean flag; '1' if a picture was taken with this initialized object
        images_per_freq                 % number of images taken at each frequency step
        binning                         % binning of camera pixels
    end
    
    properties (Constant)
        % Hamamatsu camera constants (from data sheet)
        DARK_OFFSET = 100;              % hamamatsu dark offset: 100 extra counts each time a picture is taken
        CONVERSION_FACTOR = 0.46;       % hamamatsu average photons/count
    end
    
    methods
        function obj = init(obj, binning, ccd_size, exposure_time, images_per_freq, rf_frequency, rf_power, num_freqs, num_sweeps, varargin)
            % initializes the Hamamatsu camera, DAQ system, and RF sweeper
            % binning_index: index corresponding to obj.binning setting
            % ccd_size_index: index corresponding to ccd size
            
            if ~isempty(varargin)
                ddq = varargin{1};
            else
                ddq = 0;
            end
            
            if obj.is_initialized
                error('hardware is already initialized')
            end
            
            %% SET OBJECT FIELDS
            % see "Properties" for description of fields
            obj.exposure_time_sec = exposure_time;
            obj.images_per_freq = images_per_freq;
            obj.rf_frequency = rf_frequency;
            obj.rf_power = rf_power;
            
            
            %% Set up obj.dq RF Sweeper
            if ddq
                obj.dq_rf_sweeper = visa('ni', 'ASRL7::INSTR');

                % Create the visa object if it does not exist
                if isempty(obj.dq_rf_sweeper)
                    % this is the key line. works without the rest unless program halts or something
                    obj.dq_rf_sweeper = visa('ni', 'ASRL7::INSTR');
                else
                    fclose(obj.dq_rf_sweeper);
                    obj.dq_rf_sweeper = obj.dq_rf_sweeper(1);
                end

                try
                    fopen(obj.dq_rf_sweeper);    % Connect to rf generator object, obj.rf_sweeper. 
                catch ME_OUTSIDE
                    error('turn on the RF sweeper!')
                end
                    set(obj.dq_rf_sweeper, 'Timeout', 2);
            else               
                obj.dq_rf_sweeper = visa('ni', 'ASRL4::INSTR');

                % Create the visa object if it does not exist
                if isempty(obj.dq_rf_sweeper)
                    % this is the key line. works without the rest unless program halts or something
                    obj.dq_rf_sweeper = visa('ni', 'ASRL4::INSTR');
                else
                    fclose(obj.dq_rf_sweeper);
                    obj.dq_rf_sweeper = obj.dq_rf_sweeper(1);
                end

                try
                    fopen(obj.dq_rf_sweeper);    % Connect to rf generator object, obj.rf_sweeper. 
                catch ME_OUTSIDE
                    error('turn on the RF sweeper!')
                end
                    set(obj.dq_rf_sweeper, 'Timeout', 2);
            end


            fprintf(obj.dq_rf_sweeper, 'C0');                                           % Channel A setup
            fprintf(obj.dq_rf_sweeper, ['f' num2str(obj.rf_frequency * 10^3)]);         % Sets RF frequency as rf_frequency
            fprintf(obj.dq_rf_sweeper, 'h0');                                           % Mute the RF signal
            fprintf(obj.dq_rf_sweeper, 'E1r1');                                         % RF off
            fprintf(obj.dq_rf_sweeper, 'c0');                                           % Sets power sweep mode off
            fprintf(obj.dq_rf_sweeper, ['W' num2str(obj.rf_power)]);                    % Sets the power level 
            fprintf(obj.dq_rf_sweeper, 'g0');                                           % Turn off sweep

                
            %% SET UP CCD AREA
            obj.binning = binning;
            binning_index = binning_pixels_to_index(binning);
            obj.ccd_size = ccd_size;

            if obj.binning * obj.ccd_size > 2048 % if obj.binning * ccd area is greater than total camera pixel number
                error('invalid image size requested. check obj.binning and ccd size.')
            end

            imagingM = {'MONO16_2048x2048_FastMode','MONO16_BIN2x2_1024x1024_FastMode','MONO16_BIN4x4_512x512_FastMode'};
            imagingMode = imagingM{binning_index};

            % region of interest (intensity sampling region)    
            ROIPosition = [2048/obj.binning/2 - obj.ccd_size/2 2048/obj.binning/2 - obj.ccd_size/2 obj.ccd_size obj.ccd_size];

            
            %% INITIALIZE VIDEO OBJECT
            try
                obj.vid                     = videoinput('hamamatsu', 1, imagingMode);
                 
            catch ME
                error('turn on Hamamatsu camera!')
            end
                
            src                         = getselectedsource(obj.vid);
            MAX_EXPOSURE_TIME           = 1;

            if obj.exposure_time_sec > MAX_EXPOSURE_TIME
                error('desired exposure time is too large for Hamamatsu camera')
            end
            
            src.ExposureTime            = obj.exposure_time_sec; % assign exposure time
    
            triggerconfig(obj.vid, 'hardware', 'RisingEdge', 'EdgeTrigger');
            set(obj.vid,'Timeout',50);
            obj.vid.FramesPerTrigger    = 1;
            obj.vid.TriggerRepeat       = (obj.images_per_freq * num_freqs * num_sweeps) - 1;
            obj.vid.ROIPosition         = ROIPosition; 
            obj.vid 
            start(obj.vid);
            %% set flag to initialized
            obj.is_initialized                     = 1;
        end

        
        
        function image = capture(obj)
            % With the initialized parameters, capture and return an image
            
            % check if this object is initialized
            if not(obj.is_initialized)
                error('Cannot capture image on Hardware object that is not initialized');
            end
            
            % Capture an image using the initialized Hardware object and it
            % corresponding capture parameters.

            %% TAKE PIX
            % Sets RF power ON
            fprintf(obj.dq_rf_sweeper, 'C0');
            fprintf(obj.dq_rf_sweeper, 'h1');  
                
            % START PULSE BLASTER
            if PB.start() < 0 % trigger laser, RF switch, and/or DAQ
                error('Pulseblaster could not start')
            end
            
            while get(obj.vid,'FramesAvailable') < obj.images_per_freq  %Wait until at least 1 frame is available
                  unavailable=1;     
            end
            PB.stop();

            fprintf(obj.dq_rf_sweeper, 'C0');
            fprintf(obj.dq_rf_sweeper, 'h0');    

            % now get images from camera
            try   
                image = squeeze(double(getdata(obj.vid, obj.vid.FramesAvailable)));
            catch me
                disp(me)
                obj.vid.FramesAcquired
                obj.vid.TriggerRepeat
                kill(obj)
                error('frames acquired may not be equal to frames requested')
            end
            
            
        end
        
           
        % CHANGE RF FREQUENCY
        function set_freq(obj, frequency)

            if not(obj.is_initialized)
                error('initialize the Hardware object first')
            end
            
            obj.rf_frequency = frequency;
            fprintf(obj.dq_rf_sweeper, 'C0'); 
            fprintf(obj.dq_rf_sweeper, ['f' num2str(obj.rf_frequency * 10^3)]);  
            
            
        end            
        
        % KILL
        function kill(obj)
            % called after init() and measurement_script() to kill current
            % hardware connections
            
            if obj.is_initialized == 0
                error('Hardware is not initialized')
            end
                
            % KILL HAMAMATSU CAMERA
            stop(obj.vid)
            delete(obj.vid)

            % KILL DAQ AND RF SWEEPER
            fprintf(obj.dq_rf_sweeper, 'C0');   % Channel A
            fprintf(obj.dq_rf_sweeper, 'h0');   % RF muted
            fprintf(obj.dq_rf_sweeper, 'E0r0'); % RF off
            fprintf(obj.dq_rf_sweeper, 'C1');   % Channel B
            fprintf(obj.dq_rf_sweeper, 'h0');   % RF muted
            fprintf(obj.dq_rf_sweeper, 'E0r0'); % RF off
            
            fclose(obj.dq_rf_sweeper);          % close 
            delete(obj.dq_rf_sweeper);          % delete
            
            obj.is_initialized = 0;
        end
    end
    
    methods(Static)
        
        function photons = counts2photons(camera_counts)
            % converts counts (from Hamamatsu camera) to photons
            % camera counts: number of counts over integration region
            % rf_off_counts: measured counts during laser excitation integration with RF off
                actual_counts = camera_counts - (Hardware.DARK_OFFSET);
                photons = actual_counts * Hardware.CONVERSION_FACTOR;
        end
    end
end