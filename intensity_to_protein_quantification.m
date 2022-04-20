%Peter Yong, Protein Quantification function, Final Project 12/19/19
function[sample_total_protein]=intensity_to_protein_quantification(FL_intensity, cleaved_intensity, control_ratio) %defining function to quantify total protein, control ratio is protein:intensity
    total_intensity=FL_intensity+cleaved_intensity; %adding up the intensities of the full length and cleaved protein bands
    sample_total_protein=control_ratio*total_intensity; %calculating the total protein in the sample using a ratio between protein amount and intensity found using a control
    
