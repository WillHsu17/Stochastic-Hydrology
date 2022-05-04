clear

method = cordovabras;
rainfall = method.preprocess_rainfall("Colorado.csv");
event_based_precip = method.event_def(rainfall);
% growth period 1, growth period 2, growth period 3, mm
theta_s = [246; 574; 574];

% assume initial soil moisture = 0.1 - 1 * theta_s mm
theta_i = [];
for i = 0.1:0.1:1
    theta_i(end+1:end+3) = i .* theta_s;
end
theta_i = theta_i';

% expand the saturated soil moisture to 30 times
theta_s = repmat(theta_s,10,1);

% soil moisture with star which is used for delpletion
theta_star = 0.58 .* theta_s;

% PWP for each growth period
theta_w = [74; 173; 173];
theta_w = repmat(theta_w,10,1);

% potential ET
ET_p = [3.1; 6.3; 4.6] ./ 24; % convert unit: mm/d -> mm/hr
ET_p = repmat(ET_p,10,1);

% yaron coefficient, a
a = [0.022; 0.019; 0.014];
a = repmat(a,10,1);

% build a table for params at three different stage
T_params = table(theta_i, theta_s, theta_w, theta_star, ET_p, a);
output = cell(30,5);
% row denotes different initial soil moisture at different stage
% 0.1 thetas at stage 1
% 0.1 thetas at stage 2
% 0.1 thetas at stage 3
% 0.2 thetas at stage 1, etc
% column denotes calculated data 
% 1st col: date
% 2nd col: normalized soil moisture
% 3rd col: infiltration
% 4th col: percolation
% 5th col: evapotranspiration


%%

% loop throughout all date
for j = 1:30
    thetai = T_params.theta_i(j);
    thetas = T_params.theta_s(j);
    thetaw = T_params.theta_w(j);
    thetastar = T_params.theta_star(j);
    et_p = T_params.ET_p(j);
    a = T_params.a(j);
    
    
    % initiailzed the storage for soil moisture, infiltration, percolation, and
    % evapotranspiration
    [date_N,~] = size(rainfall);
    date = event_based_precip.raw_date;
    soil_moisture = zeros(date_N,1);
    soil_moisture(1) = thetai; % initial soil moisture
    infiltration = zeros(date_N,1);
    percolation = zeros(date_N,1);
    evapotranspiration = zeros(date_N,1);
    
    for i = 1:date_N-1
        d = event_based_precip.raw_date(i);
        thetai = soil_moisture(i);
        t_r = event_based_precip.duration_date(i);
        I = event_based_precip.avg_precip(i);
        % infiltration instantaneously
        if (I ~= 0) & (event_based_precip.avg_precip(i-1) == 0)
            infiltration(i) = method.infilt(thetai,thetas,thetaw,t_r,I);
        else
            infiltration(i) = 0;
        end    

        % percolation
        percolation(i) = method.perc(thetai,thetas,thetaw);

        % ET
        evapotranspiration(i) = method.evap(thetai,thetastar,thetaw,et_p,a);

        % soil mositure at next time step
        thetai = thetai + infiltration(i) - percolation(i) - evapotranspiration(i);
        if thetai >= thetas
            thetai = thetas;
        end
        soil_moisture(i+1) = thetai;
        
    end
    % normalized the soil mositure w.r.t PWP and thetas
      normalized_soil_moisture = method.normalize(soil_moisture,thetaw,thetas);
      
    % build a output table for each condition
    output(j,:) = {date, normalized_soil_moisture, infiltration, percolation, evapotranspiration};
end

