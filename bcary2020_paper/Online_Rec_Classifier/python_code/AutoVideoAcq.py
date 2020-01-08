"""

Written by Keith Hengen with Major edits by
Alejandro Pacheco and Brian Cary

"""

import sys
import cv2
#import cv2.cv as cv
import time
import numpy as np
import csv
import os

#sys.path.insert(0, 'C:\Users\EEG\Desktop\Data\Video_files\Python_Code')

from PYTHON_VIDEOTRACK_MAIN import video_track_function

#########################
#Instructions: give the python an argument which is the sample number
# arg2 = path to save dir
# arg3 = number of seconds to run video





#########################

print '\n\n____Running Video Tracking____\n\n'
#time.sleep(0.5)


#################

sample_num = sys.argv[1]
datdir = sys.argv[2]
nsecs = int(sys.argv[3])
###Change these for setup/runs
animal = 'video'

frate = 2. # this is the number of frames collected per second

# files will turn over/save every nmin minutes
#nmin        = secs./60 # how many minutes you store in each file
ncycles     = nsecs*frate

###################



direc = os.path.join(datdir,animal)
#print 'The video data for {} will be saved in directory: {}'.format(animal,direc)
#isOkDir = input('Is this OK? (1=yes,0=no):  ')

#direc = os.path.join(datdir,animal)  
#print '\nOkay. Your video data for {} will be saved in: {}'.format(animal,direc)
#print 'I will create that folder if it does not exist.'
#direc = os.path.join('/Users/keithhengen/Desktop/animal_video_frames', 'KH67_68')
#animal = 'KH67_68'

try:
    os.stat(direc)
    # if is already a folder, find the highest number file:
    #print 'I found previous video files for this animal. \nAdding to the end of the list'
    s = [ ]
    prev = os.listdir(direc)
    if not prev:
        s = 1
    else:    
        for x in prev:
            if animal in x:
                print(x)
                s.append(int(x[x.rfind('_')+1 : x.index('.')]))
        s = max(s)+1
       
except:
#    os.mkdir(direc)
    s = 1

# set up the camera and the display window to be used
#print 'Setting up the camera'
cv2.namedWindow('ratcam', 1)

cap = cv2.VideoCapture(0)
#set acquisition frame size
cap.set(cv2.CAP_PROP_FRAME_WIDTH, 360)
cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 360)

# make a function to update the video writer - this includes the filename
def vidsaver(direc,animal,cyclecount,dim):
    
    path    = os.path.join(direc,animal+'_'+str(sample_num)+'.avi')
    #fourcc  = cv2.VideoWriter_fourcc('h', '2', '6', '3')
    fourcc  = cv2.VideoWriter_fourcc(*'XVID') # may depend on camera/graphics setup
    
    video_writer = cv2.VideoWriter(path,fourcc, 4, dim)
    
    return video_writer, path
    
    
def resetter(cyclecount,acqtime):
        # write the frametimes to a csv file:
    csvfile = os.path.join(direc,animal+'_'+str(sample_num)+'.csv')
    with open(csvfile, "w") as output:
        csvwriter = csv.writer(output, lineterminator='\n')
        for val in acqtime:
            csvwriter.writerow([val])
            
    print 'CSV file of times written successfully'
        # reset values, add to counter  
    frames  = 0
    acqtime = [ ]
    cyclecount +=1
    return(cyclecount,frames)
    


# Deal with the Metadata
#piDir = 'C:\Users\EEG\Desktop\Data\Video_files\Metadata'
piDir = datdir
checkerDir = os.path.join(piDir,'RECORD')
try:
    os.stat(checkerDir)
except:
#    os.system(cmd_string)
    print('check dir error')
    
metafile = os.path.join(checkerDir,animal+'_vid_metadata.csv')


metadat = [animal,time.strftime("%d/%m/%Y"),time.strftime("%H:%M:%S"),
           direc,str(frate),str(nsecs)]

with open(metafile, "w") as output:
    csvwriter = csv.writer(output, lineterminator='\n')
    for val in metadat:
        csvwriter.writerow([val]) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

macminichecker = os.path.join(checkerDir,'macministat.bin')
try:
    os.stat(macminichecker)
    os.remove(macminichecker)
    open(macminichecker, 'w+')
except:
    open(macminichecker, 'w+')

    

acqtime     = [ ]
frames      = 0
cyclecount  = s
c = 0
first_run = 1


while True:
    
    tic1 = time.time()
    
    acqtime.append(tic1)
      
    ret, frame = cap.read()
    frame = np.asarray(frame[:,:])
    r = 360.0 / frame.shape[1]    
    dim = (360, int(frame.shape[0]*r))
    frame = cv2.resize(frame, dim, interpolation = cv2.INTER_AREA)
    
    # update the videowriter
    if frames == 0:
        vidwriter, vid_path = vidsaver(direc,animal,cyclecount,dim)
        
    
    tic2 = time.time()
    vidwriter.write(frame)
    toc2 = time.time()
    
    cv2.imshow('ratcam', frame)
    print "Frame", frames
    
    frames += 1
    
    toc1 = time.time()
    #print 'that took', toc1-tic1
    #print toc2-tic2, ' was writing.'
    
    paws = int( round( 1000* (1./frate - (toc1-tic1) )))
    if paws<=0 or paws>=1000:
        print 'Timing is crappy.'
        (cyclecount,frames) = resetter(cyclecount,acqtime)
    else:
        c = cv2.waitKey(paws)     # wait for keyboard input
        
    #print 'Exiting waitkey'
    sys.stdout.flush()
    
       
    if frames % 10 == 0:
        timearr = np.array(time.time())
        timearr.tofile(macminichecker)
    #        with open(macminichecker,'a') as f:
    #            f.write(repr(time.time()) + '\n')
        #print 'updated PI time'
       
    if c == 27: # esc key to stop video acquisition immediately
        cv2.destroyAllWindows()
        cv2.waitKey(0)
        # write the frametimes to a csv file:
        csvfile = os.path.join(direc,animal+'_'+str(cyclecount)+'.csv')
        with open(csvfile, "w") as output:
            csvwriter = csv.writer(output, lineterminator='\n')
            for val in acqtime:
                csvwriter.writerow([val]) 
        try:
            os.stat(macminichecker)
            os.remove(macminichecker)
        except:
            exit(0)
    
        cv2.destroyAllWindows()
        vidwriter.release()
        cap.release()
        break
        
    # reset etc. after ncycles
    if frames == ncycles:
        (cyclecount,frames) = resetter(cyclecount,acqtime)
        
        vidwriter.release()
        acqtime = []
        
        if first_run == 1:
            print('first run')
            global roi_window
            roi_window = video_track_function(vid_path,0,[])
            print(roi_window)
            first_run = 0
        else:
            print(roi_window)
            roi_window = video_track_function(vid_path,0,roi_window)
        
        break
            
    
    

    
    


        
        

  
  









