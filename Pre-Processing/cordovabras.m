
classdef cordovabras
    methods(Static)
        
        function rainfall = preprocess_rainfall(data)
            p_raw = readtable(data);
            dt_raw =  datetime(p_raw.DATE,'InputFormat','yyyyMMdd HH:mm');
            for i=1:size(dt_raw)
                if p_raw.HPCP(i) > 2000
                    p_raw.HPCP(i) = 0;
                end
            end
            precip_temp =  timetable(dt_raw, p_raw.HPCP);
            rainfall = retime(precip_temp,'hourly','fillwithconstant');
        end
        function event_based_precip = event_def(rainfall)
            % Incorporater your code for rainfall event definition
            % output of the function would be a timeseries that
            % rainfall intensity during events will be constant at each
            % hour
            % and the dry periods will have zero rainfall
            raw_date = rainfall.dt_raw;
            hr_rain  = rainfall.Var1;
            count = 0; % dry hour counter
            [date_N,~] = size(raw_date);
            rain_accu = 0; % mm
            rain_duration = 0; % hr
            duration_date = zeros(date_N,1); % each rainfall duration
            starting_idx = 1; % start of rain period
            avg_precip = zeros(date_N,1);
            for i = 1:date_N
                if (hr_rain(i) == 0)
                    count = count + 1;
                    if count == 48 % if no rain for 48 hrs
                        if (rain_accu == 0 & rain_duration == 0)
                            avg_precip(starting_idx:i) = 0;
                            duration_date(starting_idx:i) = 0;
                        else
                            avg_precip(starting_idx:i) = rain_accu / rain_duration;
                            duration_date(starting_idx:i) = i-starting_idx;
                            continue
                        end
                    end
                    continue
                elseif count > 48
                    count = 0; % set dry hour counter to 0
                end
                if (count == 0)
                    starting_idx = i;
                end
                rain_accu = rain_accu + hr_rain(i);
                rain_duration = rain_duration + 1;
                count = 0;
            end
            
            event_based_precip = table(raw_date,avg_precip,duration_date);
        end

        function infiltration = infilt(theta_i, theta_s, theta_wp, t_r, i)
            %Infiltration rate
            n = 0.35; %Porosity [%]
            k_s = 30/24; %Hydraulic conductivity at saturation [mm/d]
            psi_s = 190; %Saturated soil matric potential [mm]
            d = 5.5; %Diffusivity index
            m = 0.286; %Pore Size index
            c = 10; %pore connectivity index
            %pwp = 173 %Permanent Wilting Point [mm]
            %theta_s = 747; %Soil moisture at saturation [mm]
            omega = 0;
            s_0 = (theta_i)/(theta_s);% + theta_wp);

            function phi = phi(d, s_0)
                function sum_term = sum_term(n_i)
                    sum_term = 1/(d + 5/3 - n_i) * (gamma(d+1)/(gamma(n+1)*gamma(d-n))) * (s_0/(1 - s_0))^n_i;
                end
                sum_all = 0;
                for n_i= 1:d
                    sum_all = sum_all + sum_term(n_i);
                end
                phi = (1 - s_0)^d * (1/(d + 5/3) + sum_all);
            end
            
            s = 2 * (1 - (s_0)) * ((5 * n * k_s * psi_s * phi(d, s_0))/(3 * m * pi))^0.5;
            a = 0.5 * k_s * (1 + (s_0)^c) - omega;
            t_0 = s^2/(2 * (i-a)^2);
            %print(t, a, s)
            if (t_r <= t_0)
                infiltration = i * t_r;
            else
                infiltration = a * t_r + s * sqrt(t_r/2);
            end
        end

        function perc = perc(theta_i, theta_s, theta_wp)
            k_s = 30/24; % Hydraulic conductivity at saturation [mm/d]
            c = 10; % pore connectivity index
            omega = 0;
            % _theta_s = theta_s - theta_wp % Soil moisture at saturation [mm]
            if (theta_i > theta_wp)
                perc =  k_s * ((theta_i)/(theta_s))^c - omega;
            else
                perc = 0;
            end
        end

        function evap = evap(theta_i, theta_star, theta_wp, et_p, a)
            %input:
            % theta_i: initial soil moisture at time t
            % Potential ET at time t
            % Output:
            % et_a: Actual ET
            % Parameters:
%             k_c = 0.5;      % crop coefficient
%             e_0 = 5;         % potential evaporative flux
            b = 1.0;          % coefficient
            theta_i = theta_i - theta_wp;
%             yaron_c = 0.019;
%             a = yaron_c;% 0.019%k_c * et_p/(theta_star**b)
            if ((theta_i) >= theta_star)
                et_a = et_p;
                evap = et_a;
            elseif (theta_i <= 0)
                evap = 0;
            else
                et_a = a * ((theta_i))^b;
                evap =  et_a/24; % convert unit: mm/d -> mm/hr
            end
        end
        
        function sm_norm = normalize(soil_moisture, theta_w, theta_s)
            sm_norm = (soil_moisture) ./ (theta_s);
        end
    end
end


%%

