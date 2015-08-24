import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


if __name__ == '__main__':
    # Model parameters
    K = 121    # Number of times steps
    dt = 0.01  # Delta t
    c = 0.1    # Coefficient of drag
    g = -9.81  # Gravity
    # For constructing matrix B
    bx = 0
    by = 0
    bz = 0
    bvx = 0.25
    bvy = 0.25
    bvz = 0.1
    # For constructing matrix C
    rx = 1.0
    ry = 5.0
    rz = 5.0

    # Create 3D axes for plotting
    ax = Axes3D(plt.figure())

    #####################
    # Part 1:
    #
    # Load true trajectory and plot it
    # Normally, this data wouldn't be available in the real world
    #####################

    s_true = np.loadtxt('P4_trajectory.txt', delimiter=',')
    x_coords = s_true[:,0]
    y_coords = s_true[:,1]
    z_coords = s_true[:,2]

    ax.plot(x_coords, y_coords, z_coords,
            '--b', label='True trajectory')

    #####################
    # Part 2:
    #
    # Read the observation array and plot it (Part 2)
    #####################

    s_obs = np.loadtxt('P4_measurements.txt', delimiter=',')
    stretch_mat = np.matrix([[1 / rx, 0, 0],
                             [0, 1 / ry, 0],
                             [0, 0, 1 / rz]])
    s_obs_unstretched = np.array((stretch_mat * (s_obs).transpose()))
    x_coords = s_obs_unstretched[0,:]
    y_coords = s_obs_unstretched[1,:]
    z_coords = s_obs_unstretched[2,:]

    ax.plot(x_coords, y_coords, z_coords,
            '.g', label='Observed trajectory')

    #####################
    # Part 3:
    # Use the initial conditions and propagation matrix for prediction
    # A = ?
    # a = ?
    # s = ?
    #####################

    # Initial conditions for s0
    A = np.matrix([[1, 0, 0, dt, 0, 0],
                   [0, 1, 0, 0, dt, 0],
                   [0, 0, 1, 0, 0, dt],
                   [0, 0, 0, 1 - c * dt, 0, 0],
                   [0, 0, 0, 0, 1 - c * dt, 0],
                   [0, 0, 0, 0, 0, 1 - c * dt]])
    a = np.matrix([0, 0, 0, 0, 0, g * dt]).transpose()
    s_simple = np.zeros([6, K])
    s_simple[:,0] = s_true[0,:]

    # Compute the rest of sk using Eq (1)
    for i in range(1, K):
        s_simple[:,i] = (A * np.matrix(s_simple[:,i - 1]).transpose() + a).transpose()

    s_simple = np.array(s_simple)
    x_coords = s_simple[0,:]
    y_coords = s_simple[1,:]
    z_coords = s_simple[2,:]

    ax.plot(x_coords, y_coords, z_coords,
            '-k', label='Blind trajectory')

    #####################
    # Part 4:
    # Use the Kalman filter for prediction
    # B = ?
    # C = ?
    #
    # Initial conditions for s0 and Sigma0
    # Compute the rest of sk using Eqs (2), (3), (4), and (5)
    #####################

    # ax.plot(x_coords, y_coords, z_coords,
    #         '-r', label='Filtered trajectory')

    # Show the plot
    ax.legend()
    plt.show()
