from Tracking import *
import numpy as np
import math
import matplotlib.pyplot as plt



# 3-d tracking
init_point = np.array([100, 1 , 0]).astype(float)
end_point = np.array([100, 0.5, 0.5]).astype(float)
dX = init_point - end_point
m = Tracking(3, 64, 0, 250, 1e-5, 1 - 1e-5, init_point, dX)
points, All_boxes = m.iterative_tracking()

# plotting
m.plot_trajectory(0, 250, 1e-5, 1 - 1e-5, 64)
