function [median_curve,lowprc_curve,upprc_curve,mean_curve,std_curve] = ...
    prctiles_curve(curves_vector,times_vector,center_bin_time,low_perc,up_perc)
%PRCTILES_CURVE computes median, lower/upper percentiles, mean and SD of a
%set of curves in time bins. The first and last bins are open-ended.
% INPUT:
% curves_vector: values of all the curves (stacked in a single vector)
% times_vector: times of the values in curves_vector
% center_bin_time: centers of the time bins
% low_perc, up_perc: lower and upper percentiles (e.g. 25, 75)
% OUTPUT: median_curve, lowprc_curve, upprc_curve, mean_curve, std_curve
% (one value per bin)

edges = (center_bin_time(1:end-1)+center_bin_time(2:end))/2;
m = length(center_bin_time);

median_curve = zeros(m,1);
lowprc_curve = zeros(m,1);
upprc_curve  = zeros(m,1);
mean_curve   = zeros(m,1);
std_curve    = zeros(m,1); 

% first bin
temp = curves_vector(times_vector<edges(1));
median_curve(1) = prctile(temp,50);
lowprc_curve(1) = prctile(temp,low_perc);
upprc_curve(1)  = prctile(temp,up_perc);
mean_curve(1)   = mean(temp,1,'omitnan');
std_curve(1)    = std(temp,1,'omitnan');

% middle bins
for i = 1:length(edges)-1
    temp = curves_vector(times_vector>=edges(i) & times_vector<edges(i+1));
    median_curve(i+1) = prctile(temp,50);
    lowprc_curve(i+1) = prctile(temp,low_perc);
    upprc_curve(i+1)  = prctile(temp,up_perc);
    mean_curve(i+1)   = mean(temp,1,'omitnan');
    std_curve(i+1)    = std(temp,1,'omitnan');
end

% last bin
temp = curves_vector(times_vector>=edges(i+1));
median_curve(i+2) = prctile(temp,50);
lowprc_curve(i+2) = prctile(temp,low_perc);
upprc_curve(i+2)  = prctile(temp,up_perc);
mean_curve(i+2)   = mean(temp,1,'omitnan');
std_curve(i+2)    = std(temp,1,'omitnan');

end
