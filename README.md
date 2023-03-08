# Trajecotory-tracking
Newton’s trajectory passes several cells or hyper cells in a higher dimension of parameter space.
In order to detect those cells or hyper cells, we implement the iterative tracking algorithm that trace a random trajectory in an arbitrary dimension of parameter space.
The algorithm is incrementing to the next cell iteratively by finding the minimum distance between the initial point with
respect to interfaces and moving on in a gradient direction to find the next point until reaching the last hypercubes. This algorithm has been used for netwon trajectory tracking in OBL space to detect trust region zones. Write me in case you need the paper related to this code.

