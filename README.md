# MS_2022
Code used in the completion of my master's thesis in bioengineering on the effects of obesity and type 2 diabetes on cardiac structure and arterial stiffness

Daniel Enge MS Thesis analyses

Cardiac dictionary:

•	Strain: Ratio of total deformation to the initial dimension of what is being deformed.  In this case, measured throughout entire cardiac cycle.
•	Longitudinal: Along the long axis of the heart, measuring from base to apex.
•	Circumferential: Circles the heart, how circumference of the heart (or width of each AHA section) changes.
•	Radial: Goes into the heart, from epicardium to endocardium.  How much the heart squeezes in towards itself.
•	Strain rate: Derivative of strain.  The speed at which the heart contracts (systolic) and relaxes (diastolic).
•	End-diastolic volume: Maximum blood volume of the LV during the cardiac cycle
•	End-systolic volume: Minimum blood volume of the LV during the cardiac cycle
•	Septal wall thickness: Thickness at the center of the septal wall separating the left and right ventricles
•	Free wall thickness: Thickness of the LV wall opposite of the septal wall.
Left ventricle – tissue tracking (Circle software)
Volume
1.	In Bi/Triplanar Module in Circle, trace LV in both 2 chamber views and 4 chamber views for one slice/phase.
2.	Draw manual endocardial and epicardial contours for first one, then propagate to other slices.
3.	Repeat for both phases, end-systole and end-diastole.
Wall thickness
1.	In Bi/Triplanar Module, use single line measurements on the center of the septal wall and the opposite side for a measure of free wall thickness.
Strain

1.	In Tissue Tracking Module in Circle, trace LV endo- and epicardial contours in both 2 chamber views and 4 chamber views for both end-diastole and end-systole.
2.	Place a reference point at the top of the septum, where the LV and RV meet.
3.	Perform strain analysis
4.	Save results to a scientific report and average global strain results together.
Arterial dictionary:

•	2D Phase contrast MRI: MRI encoded with both magnitude and phase in order to measure thru-plane flow.
•	Relative area change (RAC): (max area – min area)/max area *100.  Basic measure of stiffness.
•	Distensibility:  RAC/pulse pressure.	 A measure of stiffness.
•	Flow area pulse wave velocity (2D PC-MRI QA-PWV): Slope of flow-area curve during contraction (a.k.a. Wave Intensity, dI).  Represents the flow speed normalized to aortic area.
•	Wave reflection index: -AUC of negative flow/AUC of positive flow.  How much reflection is there overall in the cardiac cycle?
•	Forward compression wave (FCW): Max of positive flow. Indicates ventricular contractility, accelerates blood flow.  Early forward-traveling pushing wave.
•	Backward compression wave (BCW): Min of negative flow (most negative). Indicates initial reflection of wave, decelerates blood flow.  Late backward traveling pushing wave.
•	Forward expansion wave (FEW): Indicates late systolic flow deceleration.  Forward traveling suction wave.
•	Backward expansion wave (BEW): Indicates late systolic flow acceleration. Backward-traveling suction wave.
•	2D Transit Time PWV: The path length between ascending and descending aorta planes in 2D PC-MRI divided by the time it takes for blood to travel.
•	4D Flow PWV: A cross correlation of the flow waveforms computed at each plane along the centerline of the aorta to calculate an aggregate PWV.

3D segmentation protocol (4D flow MRI):

Redcap
1.	Log onto redcap with your designated redcap username and password
a.	Select the correct study: Retro Chart Review
b.	Under reports, select ‘Al Segmentation dashboard’ to review all active patients and summaries of where they sit currently on the dashboard 
c.	Select the record ID of the patient of interest to pull up the record home page of the patient.
i.	Status Icons exist to show the completion of the patient
ii.	Be aware of which scanner was used (3T or 1.5T)
iii.	Patients with a good, or excellent scan quality should be focused on
d.	Select MRI baseline data to search where to find the patient in the server’s database, so you can pull that patients data to your own private server.
i.	Organize the dogma and flow of data in a manner that is clear and will allow you to copy the post-processing data back to the server 
2.	Find the patient’s 4d-flow data on the lab server and move it to a local folder in your computer’s documents.
a.	Move only the 4d-flow data into a local server in this specific direction:
i.	RL
ii.	AP
iii.	FH

Matlab-Velomap
3.	Open Matlab:
a.	In the command window:
i.	>filesort_philips 
1.	select all three directions of the 4D_flow data that you copied over to convert to seamens files, via the already written code 
ii.	[yes]
iii.	>velomap_tool
1.	this allows you to pull up the newly converted 4d_flow into velomap to convert to a pc-MRA image to upload to slicer.
4.	In Velomap:
a.	Expand the window to full screen
b.	Load Data: select the ‘mag’ folder from the 4d_flow data that was synthesized in Matlab. Hit open.
c.	Select the correct user ID (00N), so the system will know who completed the pre-processing 
d.	Scan through slices in different views to ensure that everything looks normal and was transferred correctly 
i.	The 4D-flow image is best visualized at ‘time-frame’ 5, while only moving through the ‘slice number’ to visualize anatomy.
ii.	The main velocity direction window allows you to switch from Vj, Vk, and Vi which correlates too scan velocity direction. 
iii.	Find a slice number that visualizes the pulmonary trunk and Ascending Aorta to begin the pre-processing functions below:
e.	In the preprocessing functions:
i.	Select ‘Eddy Current Correction’
1.	Change the value from .1 to .02 or until a proper value looks correct to remove eddy currence from the ROI. 
2.	Select ‘preview’
3.	Crop the box surrounding the ROI to the entire scan and not just on the thorax
4.	Repeat this process for a handful of slice numbers to ensure the eddy current value is correct and the yellow color scale is not covering an anatomy of interest.
ii.	Select noise Masking, and then ‘noise filter’ into the command window on the right, and select preview to mask out all unwanted noise and tissues on the scan.
1.	Adjust noise filter threshold if needed (usually not needed)
iii.	Select Anti-aliasing
1.	Once added to the command window, select preview to apply the anti-aliasing package to the images.
2.	Scan through the slice numbers to ensure that there is no alias artifacts still, especially in the ascending and descending aorta, and the MPA.
a.	Be aware of the areas of aliasing that are uncorrectable and note them on redcap. 
b.	You will often have to change the ‘main-velocity direction’ to determine if there is aliasing in different flow velocity directions in the regions of interest.
iv.	Select PC-MRA 
1.	Under ‘Data Conversion to’ select:
a.	Ensight
b.	Avi Movie
c.	MR-Struct
d.	PC-MRA to DICOM
2.	In the middle window, select “PC-MRA 1: Mean mag weighted velocity” and “PC-MRA 5: Square root of mean mag weighted velocity”
a.	PC-MRA is the conversion file that you will need for the next steps in 3D-Slicer.
v.	Select apply to data (NOT LOAD DATA), and crop the image to your ROI (entirety of thorax and proximal neck) and wait for it to run fully without clicking out of the window. 
1.	This step takes the longest to run to make all of the required image outputs.
vi.	Once all of the data is converted to the needed files, close out of velomap.

3-D Slicer
5.	Open 3D-Slicer
a.	Load DICOM data
i.	From your local folder with the pre-process data you saved from velomap:
1.	Import mode_SqrtMeanSquares file into Slicer (this is the data that was pre-processed in velomaps) from the PCmra_dicom_subject_n_user_n folder
2.	Load that data into slicer. It should be visualized on your home page with the standard 4-window view (axial, sagittal, coronal, and 3d viewer)
b.	Adjust windowing of the scans to visualize the anatomy with a slight grey-score 
i.	Grey score and slight shadows help to find the true vessel walls of the anatomy of interest.
c.	Under modules: Open Vascular Modeling Toolkit-Level Set Segmentation 
i.	Set ‘Vesselness Volume’ to the 4D flow data set (N: 4DFLOW_AP)
ii.	Under ‘seeds’ create new MarkupFiducial as the sections of the anatomy you are segmenting. i.e DA, AscAa, Arch, BCA, LCC, LSC, MPA, RPA, and LPA
iii.	Create a new LabelMapVolume for each of these different ‘seeds’
1.	Ie. LabelmapVolume for DA, labelmapvolume1 for AscAa.
2.	The labelmapVolumes will be renamed to ‘LevelSetSegmentation’
iv.	Set two seeds on the anatomy apart from each other to grow the section in a 3D visual
1.	For example, for the DA seed, set one fiducial at the distal Aortic Arch and the 2nd fiducial at the distal descending thoracic aorta
2.	Adjust thresholding and inflation to get a good view
a.	Inflation should be at a negative number: -10-35 works the best
3.	Hit ‘preview’ and then ‘Start’
a.	Maximizing the amount of anatomy that is picked up here by the seeds will make the next steps easier.
v.	Repeat these steps for every new section of anatomy you are segmenting 
1.	Under the “…” box, hide the previous markup fiducials so they do not distract you on the next piece of anatomy.
2.	Always remember to create a new output labelmap for each new piece of anatomy or it will cancel out the previous region. You should have two markup fiducials for each region.
vi.	When segmenting the Aorta:
1.	When seeding the left common carotid artery (LCC) and left subclavian artery (LSC), segment only a few centimeters up from the arch of the aorta… about the length of the arch of the aorta up. 
vii.	When finished with these steps, turn off the markup point tool so you can click around the scan on next steps without adding a fiducial.
d.	Switch Modules to ‘Image Label Combine’
i.	Set ‘Input Label Map A’ to level set segmentation 
ii.	Set ‘Input Label Map B’ to your next label set segmentation 
iii.	Set ‘Output Label Map’ to output label map
iv.	Combine all segmentation to output label map by ‘applying’
1.	You should be able to scroll through the images and find all separate seed coloring together on the output label map.
v.	Add all different level set segmentations to the output label map before moving onto the next step.
e.	Switch modules to ‘Editor’
i.	Change Master Volume to: N: 4DFLOW_AP
ii.	Change Merge Volume to: Output Label Map
iii.	First ‘link’ all scan views together using the chain-link symbol and change the region label map segments as outlines rather than filled color to best visualize the anatomy. Both options are under the pin in the top left of each viewing frame.
iv.	Use this step to remove all unwanted segmentations in surrounding anatomy not of interest. Use the ‘paint effect tool’ to paint in the rest of the anatomy that you need by scrolling through each MRI slide until the whole anatomy is painted. 
1.	Change the paint effect label number to ‘5’ which is for vessels (red color).
2.	Change the paint brush size in the radius section to have the appropriate size brush to ‘paint’ in the anatomy.
v.	Use the erase tool, while clicked on the ‘paint’ tool to erase unwanted areas.
1.	Command’ E’ is a shortcut to switch back and forth between paint and erase
2.	You will need to change the size of the paint brush
vi.	Make sure all wanted anatomy is outlined in the correct places before moving onto the next step.

f.	Select “Model Maker’ under modules
i.	Input Volume: Output label Map
ii.	Models: Create new model Hierarchy 
iii.	Model Name: Aorta or Pulmonary
iv.	Adjust parameters to smooth until the segmentation looks clean
1.	Smooth at ‘30’ works great
v.	Apply
1.	Check the model to see if it looks correct and well done… if not, repeat previous steps
g.	Save Data
i.	Change model save to .stl
1.	Change the model name to Aorta or Pulmonary
ii.	Save to the right directory in your files
1.	ResultsN  slicer  Aorta or Pulmonary

4D Flow PWV

Need from above:
Aorta mask struct (convert from Aorta.stl with Runconversion.m code)
Mag struct
Vel struct

1.	Convert stl for aorta into a struct using Runconversion.m code
2.	Organize into folders for batch processing (only 3 above files with one folder per scan)
3.	Run batch4D_PWV code and results are saved to excel file in Batch_Analysis folder.

2D PC-MRI PWV and WIA (Wave Intensity Analysis)
1.	Trace ascending and descending aorta throughout whole cardiac cycle in Circle in the Flow module.
2.	Can use automated segmentation as a starting point to edit.
3.	Export results into scientific reports for each scan/subject.
4.	Organize excel files of scientific report for each scan from Circle into a folder.
5.	Using WIA code, select first half of data for ascending aorta or the second half of it for the descending aorta
6.	Run WIA code in this folder.
7.	For each scan, it will have a breakpoint to select 4-6 points to calculate PWV with. You will be able to see the 10 points to select from in a graph as well. If the upslope of the flow area curve is not shown like it should be for PWV calculation, you may need to check the flow and area data to figure out where to select. 
8.	You will have to stop the code and restart it each time you need to change the points selected. As you go through your folder, the results will be added onto an output variable.
2D Transit time PWV
1.	Obtain ascending and descending aorta waveforms (save in excel)
a.	You do this in 2D PWV analysis (code named WIA). Didn’t save waveforms though. They are from “flow” in start of loop in WIA code. Use raw data. This code smoothes and interpolates too.
2.	Run code (skip close all, clear all, clc step since you’ve already copied in your ascending and descending vectors). It is called PC2D_PWV_analysis_20211215.m.
3.	Code will prompt you to enter vectors for ascending and descending. Copy exactly from your excel spreadsheet (e.g. descending aorta flow should be negative)
4.	Output gives you a variety of different offsets depending on the range of data you use. There should be some point where it kind of stabilizes, but I ended up selecting 0.25-0.6 as the fit range that seemed to be most reliable.
5.	Copy over offset timing.
6.	Determine temporal resolution from 2D PC data (should be in image header, can use TR_for+PC2D_PWV.m to get TR from dicoms as you select them in MatLab)
a.	(example) TR = 20.2, TE = 2.8
7.	Multiply offset by temporal resolution
8.	Determine path length
a.	Use candy cane view and get intersection with 2D PC-MRI from Circle. Transcribe with your best estimation into the candy cane view in ImageJ. Use measuring tool for estimate. 
9.	Calculate PWV (distance/time).

AI segmentations of 4D Flow MRIs
Library of training data with:
•	mag_struct
•	vel_struct
•	aorta_mask_struct
•	pulmonary_mask_struct

2/3 of data will be used as training set with the other 1/3 used as a testing set

