%Peter Yong, Total protein of samples from volume quantification function
function[total_protein]=total_protein_quantification(sample_total_volume, load_volume, intensity, control_ratio)  %sample_total_volume is the total
%volume of the sample (mL), load_volume is the amount of sample loaded to
%the lane for the western (mL), intensity is the intensity measured for the
%band in the lane, and  control_ratio is the ratio of protein:intensity as
%determined by a control (mg/intensity)
    lane_protein=control_ratio*intensity; %calculating the total protein in the sample using a ratio between protein amount and intensity found using a control
    total_protein=(sample_total_volume*lane_protein)/load_volume; %calculating the total protein in a given sample group given the volume of the sample
    %, the volume of protein loaded for the samples, and the amount of
    %protein in the lane
