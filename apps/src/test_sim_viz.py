from paraview.simple import *
import os

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ViewSize = [720, 720]


renderView1.CenterAxesVisibility = 0
renderView1.OrientationAxesVisibility = 0

# Properties modified on renderView1
renderView1.Background = [0.5, 0.5, 0.5]

# create a new 'Box'
box1 = Box()

# Properties modified on box1
box1.ZLength = 0.0
box1.Center = [0.5, 0.5, 0.0]

# show data in view
box1Display = Show(box1, renderView1)
box1Display.DiffuseColor = [0.0, 0.0, 0.0]

# create a new 'PVD Reader'
resultspvd = PVDReader(FileName='/scratch/Cooper/chaste_test_output/Exe_VaryAdhesionOfSingleCell/sim/0_0_0/results_from_time_7.5/results.pvd')


# show data in view
resultspvdDisplay = Show(resultspvd, renderView1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 1.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.ResetCamera()
renderView1.CameraParallelScale = 0.6


time_step_info = resultspvd.TimestepValues

scene = GetAnimationScene()
# scene.UpdateAnimationUsingDataTimeSteps()
scene.NumberOfFrames = len(time_step_info)
scene.StartTime = time_step_info[0]
scene.EndTime = time_step_info[-1]

# save animation images/movie
output_directory = '/scratch/Cooper/chaste_test_output/Exe_VaryAdhesionOfSingleCell/sim/0_0_0/animation_from_time_7.5/'
WriteAnimation(output_directory + 'results.png', Magnification=1, Compression=False)

framerate = str(len(time_step_info) / 10) # make the video last 10 seconds
#For background on this line and the magic -crf value see https://trac.ffmpeg.org/wiki/Encode/H.264
os.system('avconv -r ' + framerate + ' -f image2 -i ' + output_directory + 'results.%04d.png -c:v h264 -crf 0 -y ' +output_directory + 'movie.mp4')
