function [Error_Stats] = Error_Statistics(Data_Structure,seed_instances,snr_instances,Parameter_Average,Average_Variable)
    % --------- Inputs -------------
    % Seed Instances: Number of Seeds Explored
    % snr_instances: Number of SNR values explroed
    % Parameter Average: Avaerage of Seed or SNR Values
    % Average Variable: Modes or Modal Parameters (MPS)  
    
    % Number of Modes (Needs to be Automated)
    % v1.0 2023-12-21 Added the EMA-related terms, mps_gcvd_r_mat,
    % mps_gcvd_r_err_mat, and mps_gcvd_r_err_stat 

    n = 10;
    % Extracting Data from Structures
    switch Average_Variable
        case 'Modes'
            ephi_cvd_mat = Data_Structure.cvd;
            ephi_gcvd_mat = Data_Structure.gcvd;
            ephi_dmd_mat = Data_Structure.dmd;
            ephi_era_mat = Data_Structure.era;
            ephi_itd_mat = Data_Structure.itd;            
        case 'MPS'
            mps_cvd_mat = Data_Structure.cvd;
            mps_gcvd_mat = Data_Structure.gcvd;
            mps_gcvd_r_mat = Data_Structure.gcvd_r;
            mps_dmd_mat = Data_Structure.dmd;
            mps_era_mat = Data_Structure.era;
            mps_itd_mat = Data_Structure.itd; 
        case 'MPS_Error'
            mps_cvd_err_mat = Data_Structure.cvd;
            mps_gcvd_err_mat = Data_Structure.gcvd;
            mps_gcvd_r_err_mat = Data_Structure.gcvd_r;
            mps_dmd_err_mat = Data_Structure.dmd;
            mps_era_err_mat = Data_Structure.era;
            mps_itd_err_mat = Data_Structure.itd; 
        case 'MAC'
            mac_cvd_mat = Data_Structure.cvd;
            mac_gcvd_mat = Data_Structure.gcvd;
            mac_dmd_mat = Data_Structure.dmd;
            mac_era_mat = Data_Structure.era;
            mac_itd_mat = Data_Structure.itd; 
    end 
   
    switch Average_Variable
        case 'Modes'
            if isnumeric(Parameter_Average)
                % Compute the statistics (TODO: this may be included into a function)
                % CVD
                ephi_cvd_stat(:,1,:) = median(ephi_cvd_mat, Parameter_Average);
                ephi_cvd_stat(:,2,:) = std(ephi_cvd_mat,[], Parameter_Average);
                ephi_cvd_stat(:,3,:) = min(ephi_cvd_mat,[], Parameter_Average);
                % GCVD
                ephi_gcvd_stat(:,1,:) = median(ephi_gcvd_mat, Parameter_Average);
                ephi_gcvd_stat(:,2,:) = std(ephi_gcvd_mat,[], Parameter_Average);
                ephi_gcvd_stat(:,3,:) = min(ephi_gcvd_mat,[], Parameter_Average);
               
                % DMD
                ephi_dmd_stat(:,1,:) = median(ephi_dmd_mat, Parameter_Average);
                ephi_dmd_stat(:,2,:) = std(ephi_dmd_mat,[], Parameter_Average);
                ephi_dmd_stat(:,3,:) = min(ephi_dmd_mat,[], Parameter_Average);
                % ERA
                ephi_era_stat(:,1,:) = median(ephi_era_mat, Parameter_Average);
                ephi_era_stat(:,2,:) = std(ephi_era_mat,[], Parameter_Average);
                ephi_era_stat(:,3,:) = min(ephi_era_mat,[], Parameter_Average);
                % ITD
                ephi_itd_stat(:,1,:) = median(ephi_itd_mat, Parameter_Average);
                ephi_itd_stat(:,2,:) = std(ephi_itd_mat,[], Parameter_Average);
                ephi_itd_stat(:,3,:) = min(ephi_itd_mat,[], Parameter_Average);
            else
                switch Parameter_Average
                    case 'Seed'
                        % Allocate Error Matrices
                        ephi_cvd_stat = zeros(n,3,snr_instances);
                        ephi_gcvd_stat = zeros(n,3,snr_instances);
                        ephi_dmd_stat = zeros(n,3,snr_instances);
                        ephi_era_stat = zeros(n,3,snr_instances);
                        ephi_itd_stat = zeros(n,3,snr_instances);
                        % Compute the statistics (TODO: this may be included into a function)
                        % CVD
                        ephi_cvd_stat(:,1,:) = median(ephi_cvd_mat,2);
                        ephi_cvd_stat(:,2,:) = std(ephi_cvd_mat,[],2);
                        ephi_cvd_stat(:,3,:) = min(ephi_cvd_mat,[],2);
                        % GCVD
                        ephi_gcvd_stat(:,1,:) = median(ephi_gcvd_mat,2);
                        ephi_gcvd_stat(:,2,:) = std(ephi_gcvd_mat,[],2);
                        ephi_gcvd_stat(:,3,:) = min(ephi_gcvd_mat,[],2);
                        % DMD
                        ephi_dmd_stat(:,1,:) = median(ephi_dmd_mat,2);
                        ephi_dmd_stat(:,2,:) = std(ephi_dmd_mat,[],2);
                        ephi_dmd_stat(:,3,:) = min(ephi_dmd_mat,[],2);
                        % ERA
                        ephi_era_stat(:,1,:) = median(ephi_era_mat,2);
                        ephi_era_stat(:,2,:) = std(ephi_era_mat,[],2);
                        ephi_era_stat(:,3,:) = min(ephi_era_mat,[],2);
                        % ITD
                        ephi_itd_stat(:,1,:) = median(ephi_itd_mat,2);
                        ephi_itd_stat(:,2,:) = std(ephi_itd_mat,[],2);
                        ephi_itd_stat(:,3,:) = min(ephi_itd_mat,[],2);

                    case 'SNR'
                        % Allocate Error Matrices
                        ephi_cvd_stat = zeros(n,3,seed_instances);
                        ephi_gcvd_stat = zeros(n,3,seed_instances);
                        ephi_dmd_stat = zeros(n,3,seed_instances);
                        ephi_era_stat = zeros(n,3,seed_instances);
                        ephi_itd_stat = zeros(n,3,seed_instances);
                        % Compute the statistics (TODO: this may be included into a function)
                        M = size(ephi_cvd_mat,1);
                        N = size(ephi_cvd_mat,2);
                        % CVD
                        %                     ephi_cvd_stat(:,1,:) = reshape(median(ephi_cvd_mat,3),M,[],N);
                        ephi_cvd_stat(:,1,:) = median(ephi_cvd_mat,3);
                        ephi_cvd_stat(:,2,:) = std(ephi_cvd_mat,[],3);
                        ephi_cvd_stat(:,3,:) = min(ephi_cvd_mat,[],3);
                        % GCVD
                        ephi_gcvd_stat(:,1,:) = median(ephi_gcvd_mat,3);
                        ephi_gcvd_stat(:,2,:) = std(ephi_gcvd_mat,[],3);
                        ephi_gcvd_stat(:,3,:) = min(ephi_gcvd_mat,[],3);
                        % DMD
                        ephi_dmd_stat(:,1,:) = median(ephi_dmd_mat,3);
                        ephi_dmd_stat(:,2,:) = std(ephi_dmd_mat,[],3);
                        ephi_dmd_stat(:,3,:) = min(ephi_dmd_mat,[],3);
                        % ERA
                        ephi_era_stat(:,1,:) = median(ephi_era_mat,3);
                        ephi_era_stat(:,2,:) = std(ephi_era_mat,[],3);
                        ephi_era_stat(:,3,:) = min(ephi_era_mat,[],3);
                        % ITD
                        ephi_itd_stat(:,1,:) = median(ephi_itd_mat,3);
                        ephi_itd_stat(:,2,:) = std(ephi_itd_mat,[],3);
                        ephi_itd_stat(:,3,:) = min(ephi_itd_mat,[],3);

                end
            end
            % Structure Containing the Error Statistics
            Error_Stats.cvd = ephi_cvd_stat;
            Error_Stats.gcvd = ephi_gcvd_stat;
            Error_Stats.dmd = ephi_dmd_stat;
            Error_Stats.era = ephi_era_stat;
            Error_Stats.itd = ephi_itd_stat;

        case 'MPS'
            switch Parameter_Average
                case 'Seed'
                    % Initialize the MPS Error Statistics Matrices
                    mps_cvd_stat = zeros(3,n,3,snr_instances);
                    mps_gcvd_stat = zeros(3,n,3,snr_instances);
                    mps_gcvd_r_stat = zeros(3,n,3,snr_instances);
                    mps_dmd_stat = zeros(3,n,3,snr_instances);
                    mps_era_stat = zeros(3,n,3,snr_instances);
                    mps_itd_stat = zeros(3,n,3,snr_instances);
                    % Compute the statistics 
                    % CVD
                    mps_cvd_stat(:,:,1,:) = median(mps_cvd_mat,3);
                    mps_cvd_stat(:,:,2,:) = std(mps_cvd_mat,[],3);
                    mps_cvd_stat(:,:,3,:) = min(mps_cvd_mat,[],3);
                    % GCVD
                    mps_gcvd_stat(:,:,1,:) = median(mps_gcvd_mat,3);
                    mps_gcvd_stat(:,:,2,:) = std(mps_gcvd_mat,[],3);
                    mps_gcvd_stat(:,:,3,:) = min(mps_gcvd_mat,[],3);
                    % GCVD-EMA
                    mps_gcvd_r_stat(:,:,1,:) = median(mps_gcvd_r_mat,3);
                    mps_gcvd_r_stat(:,:,2,:) = std(mps_gcvd_r_mat,[],3);
                    mps_gcvd_r_stat(:,:,3,:) = min(mps_gcvd_r_mat,[],3);
                    % DMD
                    mps_dmd_stat(:,:,1,:) = median(mps_dmd_mat,3);
                    mps_dmd_stat(:,:,2,:) = std(mps_dmd_mat,[],3);
                    mps_dmd_stat(:,:,3,:) = min(mps_dmd_mat,[],3);
                    % ERA
                    mps_era_stat(:,:,1,:) = median(mps_era_mat,3);
                    mps_era_stat(:,:,2,:) = std(mps_era_mat,[],3);
                    mps_era_stat(:,:,3,:) = min(mps_era_mat,[],3);
                    % ITD
                    mps_itd_stat(:,:,1,:) = median(mps_itd_mat,3);
                    mps_itd_stat(:,:,2,:) = std(mps_itd_mat,[],3);
                    mps_itd_stat(:,:,3,:) = min(mps_itd_mat,[],3);
                case 'SNR'
                    % Initialize the MPS Error Statistics Matrices
                    mps_cvd_stat = zeros(3,n,3,seed_instances);
                    mps_gcvd_stat = zeros(3,n,3,seed_instances);
                    mps_gcvd_r_stat = zeros(3,n,3,seed_instances);
                    mps_dmd_stat = zeros(3,n,3,seed_instances);
                    mps_era_stat = zeros(3,n,3,seed_instances);
                    mps_itd_stat = zeros(3,n,3,seed_instances);
                    % Compute the statistics 
                    % CVD
                    mps_cvd_stat(:,:,1,:) = median(mps_cvd_mat,4);
                    mps_cvd_stat(:,:,2,:) = std(mps_cvd_mat,[],4);
                    mps_cvd_stat(:,:,3,:) = min(mps_cvd_mat,[],4);
                    % GCVD
                    mps_gcvd_stat(:,:,1,:) = median(mps_gcvd_mat,4);
                    mps_gcvd_stat(:,:,2,:) = std(mps_gcvd_mat,[],4);
                    mps_gcvd_stat(:,:,3,:) = min(mps_gcvd_mat,[],4);
                    % GCVD-EMA
                    mps_gcvd_r_stat(:,:,1,:) = median(mps_gcvd_r_mat,4);
                    mps_gcvd_r_stat(:,:,2,:) = std(mps_gcvd_r_mat,[],4);
                    mps_gcvd_r_stat(:,:,3,:) = min(mps_gcvd_r_mat,[],4);
                    % DMD
                    mps_dmd_stat(:,:,1,:) = median(mps_dmd_mat,4);
                    mps_dmd_stat(:,:,2,:) = std(mps_dmd_mat,[],4);
                    mps_dmd_stat(:,:,3,:) = min(mps_dmd_mat,[],4);
                    % ERA
                    mps_era_stat(:,:,1,:) = median(mps_era_mat,4);
                    mps_era_stat(:,:,2,:) = std(mps_era_mat,[],4);
                    mps_era_stat(:,:,3,:) = min(mps_era_mat,[],4);
                    % ITD
                    mps_itd_stat(:,:,1,:) = median(mps_itd_mat,4);
                    mps_itd_stat(:,:,2,:) = std(mps_itd_mat,[],4);
                    mps_itd_stat(:,:,3,:) = min(mps_itd_mat,[],4);
            end
            % Structure Containing the Error Statistics
            Error_Stats.cvd = mps_cvd_stat;
            Error_Stats.gcvd = mps_gcvd_stat;
            Error_Stats.gcvd_r = mps_gcvd_r_stat;
            Error_Stats.dmd = mps_dmd_stat;
            Error_Stats.era = mps_era_stat;
            Error_Stats.itd = mps_itd_stat;

        case 'MPS_Error'
            if isnumeric(Parameter_Average)
                mps_cvd_err_stat(:,:,1,:) = mean(mps_cvd_err_mat, Parameter_Average);
                mps_cvd_err_stat(:,:,2,:) = std(mps_cvd_err_mat,[], Parameter_Average);
                mps_cvd_err_stat(:,:,3,:) = min(mps_cvd_err_mat,[], Parameter_Average);

                mps_gcvd_err_stat(:,:,1,:) = median(mps_gcvd_err_mat, Parameter_Average);
                mps_gcvd_err_stat(:,:,2,:) = std(mps_gcvd_err_mat,[], Parameter_Average);
                mps_gcvd_err_stat(:,:,3,:) = min(mps_gcvd_err_mat,[], Parameter_Average);

                mps_gcvd_r_err_stat(:,:,1,:) = median(mps_gcvd_r_err_mat, Parameter_Average);
                mps_gcvd_r_err_stat(:,:,2,:) = std(mps_gcvd_r_err_mat,[], Parameter_Average);
                mps_gcvd_r_err_stat(:,:,3,:) = min(mps_gcvd_r_err_mat,[], Parameter_Average);

                mps_dmd_err_stat(:,:,1,:) = median(mps_dmd_err_mat, Parameter_Average);
                mps_dmd_err_stat(:,:,2,:) = std(mps_dmd_err_mat,[], Parameter_Average);
                mps_dmd_err_stat(:,:,3,:) = min(mps_dmd_err_mat,[], Parameter_Average);

                mps_era_err_stat(:,:,1,:) = median(mps_era_err_mat, Parameter_Average);
                mps_era_err_stat(:,:,2,:) = std(mps_era_err_mat,[], Parameter_Average);
                mps_era_err_stat(:,:,3,:) = min(mps_era_err_mat,[], Parameter_Average);

                mps_itd_err_stat(:,:,1,:) = median(mps_itd_err_mat, Parameter_Average);
                mps_itd_err_stat(:,:,2,:) = std(mps_itd_err_mat,[], Parameter_Average);
                mps_itd_err_stat(:,:,3,:) = min(mps_itd_err_mat,[], Parameter_Average);
            else
                switch Parameter_Average
                    case 'Seed'
                        mps_cvd_err_stat(:,:,1,:) = median(mps_cvd_err_mat,3);
                        mps_cvd_err_stat(:,:,2,:) = std(mps_cvd_err_mat,[],3);
                        mps_cvd_err_stat(:,:,3,:) = min(mps_cvd_err_mat,[],3);

                        mps_gcvd_err_stat(:,:,1,:) = median(mps_gcvd_err_mat,3);
                        mps_gcvd_err_stat(:,:,2,:) = std(mps_gcvd_err_mat,[],3);
                        mps_gcvd_err_stat(:,:,3,:) = min(mps_gcvd_err_mat,[],3);

                        mps_gcvd_r_err_stat(:,:,1,:) = median(mps_gcvd_r_err_mat, Parameter_Average);
                        mps_gcvd_r_err_stat(:,:,2,:) = std(mps_gcvd_r_err_mat,[], Parameter_Average);
                        mps_gcvd_r_err_stat(:,:,3,:) = min(mps_gcvd_r_err_mat,[], Parameter_Average);

                        mps_dmd_err_stat(:,:,1,:) = median(mps_dmd_err_mat,3);
                        mps_dmd_err_stat(:,:,2,:) = std(mps_dmd_err_mat,[],3);
                        mps_dmd_err_stat(:,:,3,:) = min(mps_dmd_err_mat,[],3);

                        mps_era_err_stat(:,:,1,:) = median(mps_era_err_mat,3);
                        mps_era_err_stat(:,:,2,:) = std(mps_era_err_mat,[],3);
                        mps_era_err_stat(:,:,3,:) = min(mps_era_err_mat,[],3);

                        mps_itd_err_stat(:,:,1,:) = median(mps_itd_err_mat,3);
                        mps_itd_err_stat(:,:,2,:) = std(mps_itd_err_mat,[],3);
                        mps_itd_err_stat(:,:,3,:) = min(mps_itd_err_mat,[],3);
                    case 'SNR'
                        mps_cvd_err_stat(:,:,1,:) = median(mps_cvd_err_mat,4);
                        mps_cvd_err_stat(:,:,2,:) = std(mps_cvd_err_mat,[],4);
                        mps_cvd_err_stat(:,:,3,:) = min(mps_cvd_err_mat,[],4);

                        mps_gcvd_err_stat(:,:,1,:) = median(mps_gcvd_err_mat,4);
                        mps_gcvd_err_stat(:,:,2,:) = std(mps_gcvd_err_mat,[],4);
                        mps_gcvd_err_stat(:,:,3,:) = min(mps_gcvd_err_mat,[],4);

                        mps_gcvd_r_err_stat(:,:,1,:) = median(mps_gcvd_r_err_mat, 4);
                        mps_gcvd_r_err_stat(:,:,2,:) = std(mps_gcvd_r_err_mat,[], 4);
                        mps_gcvd_r_err_stat(:,:,3,:) = min(mps_gcvd_r_err_mat,[], 4);

                        mps_dmd_err_stat(:,:,1,:) = median(mps_dmd_err_mat,4);
                        mps_dmd_err_stat(:,:,2,:) = std(mps_dmd_err_mat,[],4);
                        mps_dmd_err_stat(:,:,3,:) = min(mps_dmd_err_mat,[],4);

                        mps_era_err_stat(:,:,1,:) = median(mps_era_err_mat,4);
                        mps_era_err_stat(:,:,2,:) = std(mps_era_err_mat,[],4);
                        mps_era_err_stat(:,:,3,:) = min(mps_era_err_mat,[],4);

                        mps_itd_err_stat(:,:,1,:) = median(mps_itd_err_mat,4);
                        mps_itd_err_stat(:,:,2,:) = std(mps_itd_err_mat,[],4);
                        mps_itd_err_stat(:,:,3,:) = min(mps_itd_err_mat,[],4);
                end
            end
            % Structure Containing the Error Statistics
            Error_Stats.cvd = mps_cvd_err_stat;
            Error_Stats.gcvd = mps_gcvd_err_stat;
            Error_Stats.gcvd_r = mps_gcvd_r_err_stat;
            Error_Stats.dmd = mps_dmd_err_stat;
            Error_Stats.era = mps_era_err_stat;
            Error_Stats.itd = mps_itd_err_stat;
        case 'MAC'
            if isnumeric(Parameter_Average)
                mac_cvd_stat(:,:,1,:) = median(mac_cvd_mat,Parameter_Average);
                mac_cvd_stat(:,:,2,:) = std(mac_cvd_mat,[],Parameter_Average);
                mac_cvd_stat(:,:,3,:) = min(mac_cvd_mat,[],Parameter_Average);

                mac_gcvd_stat(:,:,1,:) = median(mac_gcvd_mat,Parameter_Average);
                mac_gcvd_stat(:,:,2,:) = std(mac_gcvd_mat,[],Parameter_Average);
                mac_gcvd_stat(:,:,3,:) = min(mac_gcvd_mat,[],Parameter_Average);

                mac_dmd_stat(:,:,1,:) = median(mac_dmd_mat,Parameter_Average);
                mac_dmd_stat(:,:,2,:) = std(mac_dmd_mat,[],Parameter_Average);
                mac_dmd_stat(:,:,3,:) = min(mac_dmd_mat,[],Parameter_Average);

                mac_era_stat(:,:,1,:) = median(mac_era_mat,Parameter_Average);
                mac_era_stat(:,:,2,:) = std(mac_era_mat,[],Parameter_Average);
                mac_era_stat(:,:,3,:) = min(mac_era_mat,[],Parameter_Average);

                mac_itd_stat(:,:,1,:) = median(mac_itd_mat,Parameter_Average);
                mac_itd_stat(:,:,2,:) = std(mac_itd_mat,[],Parameter_Average);
                mac_itd_stat(:,:,3,:) = min(mac_itd_mat,[],Parameter_Average);
            else
                switch Parameter_Average
                    case 'Seed'
                        mac_cvd_stat(:,:,1,:) = median(mac_cvd_mat,3);
                        mac_cvd_stat(:,:,2,:) = std(mac_cvd_mat,[],3);
                        mac_cvd_stat(:,:,3,:) = min(mac_cvd_mat,[],3);

                        mac_gcvd_stat(:,:,1,:) = median(mac_gcvd_mat,3);
                        mac_gcvd_stat(:,:,2,:) = std(mac_gcvd_mat,[],3);
                        mac_gcvd_stat(:,:,3,:) = min(mac_gcvd_mat,[],3);

                        mac_dmd_stat(:,:,1,:) = median(mac_dmd_mat,3);
                        mac_dmd_stat(:,:,2,:) = std(mac_dmd_mat,[],3);
                        mac_dmd_stat(:,:,3,:) = min(mac_dmd_mat,[],3);

                        mac_era_stat(:,:,1,:) = median(mac_era_mat,3);
                        mac_era_stat(:,:,2,:) = std(mac_era_mat,[],3);
                        mac_era_stat(:,:,3,:) = min(mac_era_mat,[],3);

                        mac_itd_stat(:,:,1,:) = median(mac_itd_mat,3);
                        mac_itd_stat(:,:,2,:) = std(mac_itd_mat,[],3);
                        mac_itd_stat(:,:,3,:) = min(mac_itd_mat,[],3);
                    case 'SNR'
                        mac_cvd_stat(:,:,1,:) = median(mac_cvd_mat,4);
                        mac_cvd_stat(:,:,2,:) = std(mac_cvd_mat,[],4);
                        mac_cvd_stat(:,:,3,:) = min(mac_cvd_mat,[],4);

                        mac_gcvd_stat(:,:,1,:) = median(mac_gcvd_mat,4);
                        mac_gcvd_stat(:,:,2,:) = std(mac_gcvd_mat,[],4);
                        mac_gcvd_stat(:,:,3,:) = min(mac_gcvd_mat,[],4);

                        mac_dmd_stat(:,:,1,:) = median(mac_dmd_mat,4);
                        mac_dmd_stat(:,:,2,:) = std(mac_dmd_mat,[],4);
                        mac_dmd_stat(:,:,3,:) = min(mac_dmd_mat,[],4);

                        mac_era_stat(:,:,1,:) = median(mac_era_mat,4);
                        mac_era_stat(:,:,2,:) = std(mac_era_mat,[],4);
                        mac_era_stat(:,:,3,:) = min(mac_era_mat,[],4);

                        mac_itd_stat(:,:,1,:) = median(mac_itd_mat,4);
                        mac_itd_stat(:,:,2,:) = std(mac_itd_mat,[],4);
                        mac_itd_stat(:,:,3,:) = min(mac_itd_mat,[],4);
                end
            end
            % Structure Containing the Error Statistics
            Error_Stats.cvd = mac_cvd_stat;
            Error_Stats.gcvd = mac_gcvd_stat;
            Error_Stats.dmd = mac_dmd_stat;
            Error_Stats.era = mac_era_stat;
            Error_Stats.itd = mac_itd_stat;
    end

end 