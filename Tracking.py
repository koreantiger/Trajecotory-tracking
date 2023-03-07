import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from darts.engines import value_vector, index_vector
from darts.models.physics.dead_oil import DeadOil

"""
This class, we track newton trajectory in an arbitrary dimension of OBL parameter space
"""
class Tracking:
     def __init__(self, nc, n_point, min_p, max_p, min_z, max_z, X, dX):
         # constructor which initialize my parameters
         self.eps = 1e-5
         self.n_point = n_point    # OBL resolution
         self.nc = nc
         self.p_vec = np.linspace(min_p, max_p, self.n_point)
         self.z_vec = np.linspace(min_z, max_z, self.n_point)
         self.init_point = X
         self.init_pointt = X    # for plotting
         self.end_point = X - dX
         self.C = self.end_point - self.init_point
         self.point_out_boundary(self.init_point)
         self.point_out_boundary(self.end_point)
         self.idx_obl_node_init = self.idx_obl_node(self.init_point)
         self.idx_obl_node_end = self.idx_obl_node(self.end_point)
         self.points = []




     def iterative_tracking(self):
         # method where iteratively we increment from the init_point to reach to the last point
         self.points = []
         self.All_boxes = []
         ctr = 0
         while (np.array_equal(self.idx_obl_node_init, self.idx_obl_node_end) == False):
             self.points.append(self.init_point.tolist())
             self.sigma = self.compute_sigma()
             self.init_point = self.next_point()
             if len(self.points) > self.n_point*10:
                 # criteria to stop while loop in case of not converging
                  break
             self.All_boxes.append(self.idx_obl_node_init.tolist())
             x=  np.array(self.idx_obl_node_init, copy = True)
             self.box()
             ctr+=1
         return self.points, self.All_boxes


     def point_out_boundary(self, X):
         # method to block the points from getting out of OBL limits
         # X is the point
         for i in range(self.nc):
                 if (i == 0):
                         if X[i] < self.p_vec[0]:
                                 X[i] = self.p_vec[0] + self.eps
                         elif X[i] > self.p_vec[-1]:
                                 X[i] = self.p_vec[-1] - self.eps
                 else:
                         if X[i] < self.z_vec[0]:
                                 X[i] = self.z_vec[0] + self.eps
                         elif X[i] > self.z_vec[-1]:
                                 xx =self.z_vec[-1] - self.eps
                                 X[i] = xx



     def box(self):
         # Method to detect the obl node where our point[s] is located in the parameter space
         # Note that we extract only one node of the cube or hypercube which is the upper node since using math.ceil rounding to the upper integer

         for i in range(self.nc):
             if (i == 0):
                 self.idx_obl_node_init[i] = math.ceil((self.init_point[i]- self.p_vec[0])/(self.p_vec[1]-self.p_vec[0]))
                 if(self.idx_obl_node_init[i] < 0):
                      self.idx_obl_node_init[i] = 0
                 elif(self.idx_obl_node_init[i] > len(self.p_vec) - 1):
                      self.idx_obl_node_init[i] = len(self.p_vec) - 1

             else:
                self.idx_obl_node_init[i] = int(math.ceil((self.init_point[i]- self.z_vec[0])/(self.z_vec[1]-self.z_vec[0])))
                if (self.idx_obl_node_init[i] < 0):
                    self.idx_obl_node_init[i] = 0
                elif (self.idx_obl_node_init[i] > len(self.z_vec) - 1):
                    self.idx_obl_node_init[i] = len(self.z_vec) - 1



     def compute_sigma(self):
         # compute sigma (the distance from each interface of obl cube[hypercube])
         # find the minimum of that since we move to that direction

         sigma_ar = np.zeros(self.nc)  # initiliaze the array of sigma which is the distance from each interface

         for i in range(self.nc):
             if (i == 0):
             # for pressure interface since limit of the array is different
                 if self.C[i] == 0:
                     # no movement in that direction,( make sigma a big number )
                     sigma_ar[i] = 100000
                 else:
                     sigma_ar[i] = max((self.p_vec[self.idx_obl_node_init[i]]- self.init_point[i])/self.C[i] , (self.p_vec[self.idx_obl_node_init[i]-1] - self.init_point[i]) /self.C[i])
             else:
                 if self.C[i] == 0:
                     # no movement in that direction ( make sigma a big number )
                     sigma_ar[i] = 100000
                 else:
                     sigma_ar[i] = max((self.z_vec[self.idx_obl_node_init[i]] - self.init_point[i]) / self.C[i], (self.z_vec[self.idx_obl_node_init[i] - 1] - self.init_point[i]) / self.C[i])




         sigmaa = np.min(sigma_ar)  # Find the minimum of sigma since we move in that direction

         return sigmaa

     def compute_sigma_vectorize(self):
         # compute sigma (the distance from each interface of obl cube[hypercube])
         # find the minimum of that since we move to that direction

         sigma_ar = np.zeros(self.nc)  # initiliaze the array of sigma which is the distance from each interface
         idx_zero = np.where(self.C == 0)[0]
         if (idx_zero.size != 0):
             self.C[idx_zero]  = 1e5
         sigma_ar[0] = max((self.p_vec[self.idx_obl_node_init[0]]- self.init_point[0])/self.C[0] , (self.p_vec[self.idx_obl_node_init[0]-1] - self.init_point[0]) /self.C[0])
         sigma_ar[1:] = np.maximum((self.z_vec[self.idx_obl_node_init[1:]] - self.init_point[1:]) / self.C[1:],
                    (self.z_vec[self.idx_obl_node_init[1:] - 1] - self.init_point[1:]) / self.C[1:])
         sigmaa = np.min(sigma_ar)  # Find the minimum of sigma since we move in that direction
         return sigmaa


     def idx_obl_node(self, point):
         # Method to extract the obl node where our point is located in the parameter space
         # Note that we extract only one node of the cube or hypercube which is the upper node since using math.ceil
         idx_array = np.zeros(self.nc).astype(int)
         for i in range(self.nc):
             if (i == 0):
                 idx_array[i] = math.ceil((point[i] - self.p_vec[0]) / (self.p_vec[1] - self.p_vec[0]))
                 if (idx_array[i] < 0):
                     idx_array[i] = 0
                 elif (idx_array[i] > len(self.p_vec) - 1):
                     idx_array[i] = len(self.p_vec) - 1

             else:
                 idx_array[i] = math.ceil((point[i] - self.z_vec[0]) / (self.z_vec[1] - self.z_vec[0]))
                 if (idx_array[i] <=  0):
                     idx_array[i] = 1 # always rounded to the top one
                 elif (idx_array[i] > len(self.z_vec) - 1):
                     idx_array[i] = len(self.z_vec) -1
         return idx_array

     def idx_obl_node_vectorize(self, point):
        # Method to extract the obl node  where our point is located in the parameter space
        # Note that we extract only one node of the cube or hypercube which is the upper node since using math.ceil

        idx_array = np.zeros(self.nc).astype(int)
        xx= np.where(idx_array < 0)[0]

        idx_array[0] = math.ceil((point[0] - self.p_vec[0]) / (self.p_vec[1] - self.p_vec[0]))
        idx_array[1:] = np.ceil(point[1:] - self.z_vec[0]) / (self.z_vec[1] - self.z_vec[0])
        if (idx_array[i] < 0):
            idx_array[i] = 0
        elif (idx_array[i] > len(self.p_vec) - 1):
            idx_array[i] = len(self.p_vec) - 1


        for i in range(self.nc):
            if (i == 0):
                idx_array[i] = math.ceil((point[i] - self.p_vec[0]) / (self.p_vec[1] - self.p_vec[0]))
                if (idx_array[i] < 0):
                    idx_array[i] = 0
                elif (idx_array[i] > len(self.p_vec) - 1):
                    idx_array[i] = len(self.p_vec) - 1

            else:
                idx_array[i] = math.ceil((point[i] - self.z_vec[0]) / (self.z_vec[1] - self.z_vec[0]))
                if (idx_array[i] < 0):
                    idx_array[i] = 0
                elif (idx_array[i] > len(self.z_vec) - 1):
                    idx_array[i] = len(self.z_vec) - 1
        return idx_array

     def next_point(self):
        # Increment our point [sigma * C] to reach next point in new cube
        nextpoint = self.init_point + (self.sigma + self.eps) * self.C

        return nextpoint


     def plot_trajectory(self, p_min, p_max, z_min, z_max, n_points):
         # plotting trajectory in each newton iteration
         p_vec = np.linspace(p_min, p_max, n_points)
         z_vec = np.linspace(z_min, z_max, n_points)
         fig = plt.figure()
         if (len(self.init_pointt) < 3):
             # Testing in 2 Dimension
             ax1 = fig.add_subplot(111)
             ax1.plot([self.init_pointt[1], self.end_point[1]], [self.init_pointt[0], self.end_point[0]])
             ax1.scatter(self.init_pointt[1], self.init_pointt[0], facecolor='red')
             ax1.scatter(self.end_point[1], self.end_point[0], facecolor='red')
             majors = [z_vec[0], z_vec[-1], n_points]
             ax1.set_xlim(0, 1)
             ax1.set_ylim(p_vec[0], p_vec[-1])
             ax1.xaxis.set_minor_locator(MultipleLocator(max(np.diff(z_vec))))
             ax1.yaxis.set_minor_locator(MultipleLocator(max(np.diff(p_vec))))
             ax1.grid(which='minor')
             for i in range(len(self.points)):
                 x, y = self.points[i]
                 ax1.scatter(y, x, s=10, facecolor='green')
                 plt.pause(0.1)
             """
                for i in range(len(self.All_boxes)):
                    p, z = self.All_boxes[i]
                    ax1.scatter(z_vec[z],p_vec[p], s=10, facecolor='blue', marker="x")
             """
         elif len(self.init_pointt) == 3:
             # Testing in 3 Dimension
             ax1 = fig.add_subplot(111, projection='3d')
             ax1.plot([self.init_pointt[2], self.end_point[2]], [self.init_pointt[1], self.end_point[1]], [self.init_pointt[0], self.end_point[0]])
             ax1.scatter(self.init_pointt[2], self.init_pointt[1], self.init_pointt[0], facecolor='red')
             ax1.scatter(self.end_point[2], self.end_point[1], self.end_point[0], facecolor='red')
             ax1.set_xlim(0, 1)
             ax1.set_ylim(0, 1)
             ax1.set_zlim(p_vec[0], p_vec[-1])
             ax1.xaxis.set_minor_locator(MultipleLocator(max(np.diff(z_vec))))
             ax1.yaxis.set_minor_locator(MultipleLocator(max(np.diff(z_vec))))
             ax1.zaxis.set_minor_locator(MultipleLocator(max(np.diff(p_vec))))
             ax1.grid(which = 'minor')
             for i in range(len(self.points)):
                 x, y, z = self.points[i]
                 ax1.scatter(z, y, x, s=10, facecolor='green')
                 plt.pause(0.1)
         plt.show()
         return 0


