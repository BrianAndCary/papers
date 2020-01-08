Important code used for analysis and experimentation in a 2020 paper from Brian A. Cary

common_dependencies:
	Directory with functions that are used or might be used by multiple functions/programs

mini_FI_GUI:
	contains "Mini_FI_GUI_v#" program which can be ran and used as a method for mini EPSC/IPSC detection and extraction.
Right now, data needs to be structured in a particular way in which the program expects for matlab files. Can also read in 
".ibw" files from IGOR, but this is far less tested. The program also has a rudimentary ability to detect and analysis spike
number elicted from current injection.

Online_Rec_Classifier:
	This is a program used for Automated and online recording and analysis of EEG/EMG data used for real time detection
of sleep/wake dense episodes. Dependent on NIDAQ board connected to external amplifier. It is dependent on camera settings to some extent. The timing of the program and its ability to 
record and analyze video in real time is crucial to it working. Would require some initial testing to make sure the python script
that records video works on a new computer; the video codec must match what the computer can handle.

StateCoder_GUI:
	This is a GUI that was used for post-hoc classification of behavioral state using both an automated threshold-based
algorithm and a random forest machine learning classifier with enough training data. This program again assumes that the raw
EEG/EMG data is structured in a particular way as  a ".mat" file, and assumes certain structure for video files. Creates additional
folder "tempPyVidcodeFiles" if not there already.