%Peter Yong, GCD 5005 Final Project: Western Analysis
%Western Analysis of samples from different steps in a nickel affinity
%column protein purification, using intensity to quantify protein and
%determine protien yield efficiency
clear;
clc;
close all;

offset=10; %setting an offset value to make upper and lower lines to create a region box during image analysis
cleaved_offset=120; %setting a value that represents how far below on the image the cleaved protein at 15 kDa is from the full length
%protein at 48 kDa
FL_intensity_sum=0; %initially setting the variable that collects the intensity of the whole region selected for the full length protein band to 0
cleaved_intensity_sum=0; %initially setting the variable that collects the intensity of the whole region selected for the cleaved protein band to 0
control_ratio=5E-8; %mg of protein per intensity as determined by a loading control of known concentration and intensity measurement
load_vol=input('Please enter the volume loaded into each lane (uL): '); %enter the volume of the sample loaded on the gel for western blot
load_vol=load_vol*1E-3; %converting uL of load volume to mL
extract_vol=input('Please enter the extract volume (mL): '); %enter the total volume of the extract sample
FT_vol=input('Please enter the flow through volume (mL): '); %enter the total volume of the flow through (usually the same as extract volume as all
%of the extract is allowed to flow through)
wash_vol=input('Please enter the volume of the wash aliquots (mL): '); %enter the volume of each wash collected, in this case washes were collected in 1 mL portions
elute_vol=input('Please enter the volume of the elution aliquots (mL): '); %enter the volume of each elution  collected, in  this case elutions were collected in 1 mL portions

full_intensity_array=zeros(1,14); %preallocating a full length intensity array for each band
cleaved_intensity_array=zeros(1,14); %preallocating the cleaved intensity array for each band
total_protein_in_lanes=zeros(1,14); %preallocating an array that represents the amount of cleaved and full length protein combined in each lane
intensity_ratio_array=zeros(1,14); %preallocating array that represents the ratio between the intensities of the cleaved protein band and full length protein band
western_1=imread('Western_1.tif'); %reading in loading control image
imshow(western_1, [], 'InitialMagnification', 'fit'); %showing the loading control image
for i=1:14 %for loop going through 14 times, once for each blot band
    FL_intensity_sum=0; %resetting full length protein band intensity sum to 0 at the beginning of the loop for a new band
    cleaved_intensity_sum=0; %resetting cleaved protein band intensity sum to 0 at the beginning of the loop for a new band
    hold on; %hold
    [x_coord, y_coord, intensities]=impixel(); %impixel to select points around each side of the protein band and find pixel coordinates/intensities
    y_plus=mean(y_coord)+offset; %creating a y value adding the offset to make box region  that will encompass the band
    y_minus=mean(y_coord)-offset; %creating a y value subtracting the offset to make box region that will encompass the band
    plot([x_coord(1) x_coord(2)], [mean(y_coord) mean(y_coord)], 'w.-'); %plotting the mean y coord selected on the image for each selected x coord
    plot([x_coord(1) x_coord(2)], [y_plus y_plus], 'w.-'); %plotting the lower side of box on western blot image
    plot([x_coord(1) x_coord(2)], [y_minus y_minus], 'w.-'); %plotting the upper side of box on western blot image
    for n=x_coord(1,1):x_coord(2,1) %for loop going through all pixels within box in the x direction
       for j=y_minus:y_plus %for loop going through all pixels within box in the y direction
           [x_coord_2, y_coord_2, intensities_2]=impixel(western_1, n, j); %for the pixel, getting info
           FL_intensity_sum=FL_intensity_sum+intensities_2(1); %adding the red intensity for each pixel to a variable collecting the sum of all intensities in the region
       end
    end
    cleaved_y_coord=y_coord+cleaved_offset; %selecting coordinates at the same x positions around the band but changing y coordinates to be around the
    % cleaved 15 kDa protein band
    cleaved_y_plus=mean(cleaved_y_coord)+offset; %creating lower side of a box to encompass the fluorescene intensities for the band
    cleaved_y_minus=mean(cleaved_y_coord)-offset; %creating upper side of a box to encompass the fluorescene intensities for the band
    plot([x_coord(1) x_coord(2)], [mean(cleaved_y_coord) mean(cleaved_y_coord)], 'b.-'); %plotting the mean cleaved y coord at 15 kDa on the image
    plot([x_coord(1) x_coord(2)], [cleaved_y_plus cleaved_y_plus], 'b.-'); %plotting the lower side of box on western blot image
    plot([x_coord(1) x_coord(2)], [cleaved_y_minus cleaved_y_minus], 'b.-'); %plotting the upper side of box on western blot image
    for p=x_coord(1,1):x_coord(2,1) %for loop going through all pixels within box in the x direction
       for q=cleaved_y_minus:cleaved_y_plus %for loop going through all pixels within box in the y direction
           [x_coord_3, y_coord_3, intensities_3]=impixel(western_1, p, q); %for the pixel, getting info
           cleaved_intensity_sum=cleaved_intensity_sum+intensities_3(1); %adding the red intensity for each pixel
       end
    end

    full_intensity_array(i)=FL_intensity_sum; %adding the total intensity for a band to an array
    cleaved_intensity_array(i)=cleaved_intensity_sum; %adding the total intensity for the cleaved protein band for each lane
    total_protein_in_lanes(i)=intensity_to_protein_quantification(full_intensity_array(i), cleaved_intensity_array(i), control_ratio); %calculating the total protein in each lane and adding it to an array using a function, protein in mg
    intensity_ratio=1-(cleaved_intensity_sum/FL_intensity_sum);  %calculating the relative intensity of full length protein to full length + cleaved protein in each lane
    intensity_ratio_array(i)=intensity_ratio; %adding the intensity ratio calculated above to an array for each lane

end
extract_total_protein=total_protein_quantification(extract_vol, load_vol, full_intensity_array(1), control_ratio); %extract_total_protein in mg,
%calculated using the function total_protein_quantification that takes in the volume
%of the sample (mL), gel loading volume (uL), intensity of the sample,  and
%the control ratio between protein (mg) and intensity
FT_total_protein=total_protein_quantification(FT_vol, load_vol, full_intensity_array(2), control_ratio); %same as above but with flow through
W1_total_protein=total_protein_quantification(wash_vol, load_vol, full_intensity_array(3), control_ratio); %same as above
W2_total_protein=total_protein_quantification(wash_vol, load_vol, full_intensity_array(4), control_ratio); %but
W3_total_protein=total_protein_quantification(wash_vol, load_vol, full_intensity_array(5), control_ratio); %with
W8_total_protein=total_protein_quantification(wash_vol, load_vol, full_intensity_array(6), control_ratio); %the
W9_total_protein=total_protein_quantification(wash_vol, load_vol, full_intensity_array(7), control_ratio); %washes
E1_total_protein=total_protein_quantification(elute_vol, load_vol, full_intensity_array(8), control_ratio); %same
E2_total_protein=total_protein_quantification(elute_vol, load_vol, full_intensity_array(9), control_ratio); %as
E3_total_protein=total_protein_quantification(elute_vol, load_vol, full_intensity_array(10), control_ratio); %above
E4_total_protein=total_protein_quantification(elute_vol, load_vol, full_intensity_array(11), control_ratio); %but
E5_total_protein=total_protein_quantification(elute_vol, load_vol, full_intensity_array(12), control_ratio); %with
E6_total_protein=total_protein_quantification(elute_vol, load_vol, full_intensity_array(13), control_ratio); %the
E7_total_protein=total_protein_quantification(elute_vol, load_vol, full_intensity_array(14), control_ratio); %elutions

percent_bound=1-(FT_total_protein/extract_total_protein); %calculating the amount of protein that binded to the column when running the extract through by
%comparing the total amount of protein in the extract with the total amount
%of protein left in the flow through after going through the column
sum_elution_protein=E1_total_protein+E2_total_protein+E3_total_protein+E4_total_protein+E5_total_protein+E6_total_protein+E7_total_protein;
%calculating the sum of the protein amount in the elutions collected
percent_eluted=(sum_elution_protein/extract_total_protein); %calculating the amount of protein in the elutions relative to the amount of protein we started
%with in the extract

step_array=zeros(1,3);  %preallocating for array that represents step in purification
for f=1:3 %for loop  going through the 3 purification steps of interest
    step_array(f)=f; %setting values in array between 1 and 3
end

percentage_array=[]; %creating an empty array
percentage_array(1)=1; %setting the total protein value to 1
percentage_array(2)=percent_bound;  %setting the second spot in the array to be the amount of protein that bound to the column as the extract flowed through
percentage_array(3)=percent_eluted; %setting the third spot in the array to be the amount of protein that was collected in elutions in relation to the total protein in the extract
figure; %opening new figure window
plot(step_array, percentage_array, 'g*-'); %plotting the steps on x and the percentage yields on y with green astericks and lines connecting
axis([1 3 0 1]); %setting the axis limits
xlabel('Purification Step'); %adding an x axis labbel
ylabel('Relative Protein Retained'); %adding a y axis label
title('Protein Purification Yield'); %adding a graph title

total_steps_array=zeros(1,14); %making an array that counts for each sample in the purification
for r=1:14 %for loop going through for each band
    total_steps_array(r)=r; %increasing value in each place by one
end

figure; %opening a new figure window
plot(total_steps_array, intensity_ratio_array, 'b*-'); %plotting the intensity ratio between full length protein and all protein for each step
axis([1 14 0 1]); %setting axis limits
xlabel('Purification Step'); %adding an x axis label
ylabel('Relative Amount of Full Length Protein'); %adding a y axis label
title('Relative Amount of Full Length Protein to Cleaved Protein in Purification Steps'); %adding a title
