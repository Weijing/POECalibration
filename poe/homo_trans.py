
import numpy as np
import math


def rot3_x(a):
    """
        generate 3 * 3 rotation matrix from angle a,
        angle in radius
    """
    t = np.array([[1.0, 0.0, 0.0], [0.0, math.cos(a), math.sin(a)],
                  [0.0, -math.sin(a), math.cos(a)]])
    return t


def rot3_z(a):
    t = np.array([[math.cos(a), math.sin(a), 0.0], 
                    [-math.sin(a), math.cos(a), 0.0],
                  [0.0, 0.0, 1.0]])
    return t


def rot3_y(a):
    t = np.array([[math.cos(a), 0.0, -math.sin(a), ], [0.0, 1.0, 0.0],
                  [math.sin(a), 0.0,  math.cos(a)]])
    return t


def rot4_x(a):
    r = rot3_x(a)
    b = np.zeros((3, 1))
    c = np.zeros((1, 3))
    d = np.ones((1, 1))
    t = np.vstack([np.hstack([r, b]),
                   np.hstack([c, d])])
    return t


def rot4_y(a):
    r = rot3_y(a)
    b = np.zeros((3, 1))
    c = np.zeros((1, 3))
    d = np.ones((1, 1))
    t = np.vstack([np.hstack([r, b]),
                   np.hstack([c, d])])
    return t


def rot4_z(a):
    r = rot3_z(a)
    b = np.zeros((3, 1))
    c = np.zeros((1, 3))
    d = np.ones((1, 1))
    t = np.vstack([np.hstack([r, b]),
                   np.hstack([c, d])])
    return t


def trans4(pos):
    """
        return a 4*4 transformation matrix, pos must be a 1*3 np array
    """
    t = np.eye(4)
    t[0:3, 3] = pos.transpose()
    return t


def rp_homo4(r, p):
    """
        convert 3*1 rotation and 1*3 postion to 4*4 homogeneous 
        transformation matrix
    """
    c = np.zeros((1, 3))
    d = np.ones((1, 1))
    # print r, p.transpose(), p
    p = p[:, None]   # convert 1D numpy array to 2D
    # print p
    t = np.vstack([np.hstack([r, p]),
                   np.hstack([c, d])])
    return t


def homo4_pos(t):
    """
        return 1*3 position from 4*4 homogeneous transformation matrix
    """
    p = t[:3, 3]
    return p.flatten()


def homo4_r3(t):
    """
        return 3*3 rotaion matrix from 4*4 homogeneous transformation matrix
    """
    return t[:3,:3]

def homo4_euler(t):
    """
        convert 4*4 homogeneous transformation matrix to 1*6 position
        and euler angle rotaion
        res[:3] is pos, res[3:] is euler angle
    """
    pass

# ====================================================================
if __name__ == "__main__":
    print(__file__)
    # print rot3_x(math.pi/4)
    # print rot3_y(math.pi/4)
    # print rot3_z(math.pi/4)
    # print rot4_x(math.pi/4)
    # print trans4(np.array([1, 2, 3]))
    # print homo4_pos(np.eye(4))
