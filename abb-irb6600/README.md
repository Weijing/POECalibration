the data is measured by Leica tracker

given m is the number of measured data, n is the number of joints.

The folder should consist of the following item:
	1. joints.txt: contains the joints angle data (m*n) 
	2. t_measured.txt: 4n*4 measured data (where is m*(4*4 SE(3)) )
	3. s_vector.txt: 6m*1 data defines the rotation axis for each joints
	4. initial_kinematics.txt: 4m*4 data (m * (4*4 SE(3))), the inital kinematics 
