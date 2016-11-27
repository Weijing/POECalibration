import numpy as np
import math
# import warnings

# import modules from the same package
from homo_trans import *
from poe_model import *


# ==========================================

def adjoint_rep(t):
    """
        the adjoint representation of matrix, refer to lie algebra
    """
    p = homo4_pos(t)
    r = homo4_r3(t)
    p_skew = skew_sym(p)
    # print p,r,p_skew
    res = np.vstack([np.hstack([r, np.dot(p_skew, r)]),
                     np.hstack([np.zeros((3, 3)), r])])
    return res


def log_v(t):
    """
        maps the 4*4 homogenous transformation matrix to 6*1 vector form
    """
    p = homo4_pos(t)
    r = homo4_r3(t)
    r_tr = 0.5 * (np.trace(r) - 1)

    SMALL_NUM = 0.00001

    if math.fabs(r_tr) >= 1.0:
        phi = 0.0
    else:
        phi = math.acos(r_tr)

    if math.fabs(phi) < SMALL_NUM:
        w_hat = 0.5 * (r - r.transpose())
    else:
        w_hat = phi * (r - r.transpose()) / (2 * math.sin(phi))

    w = np.array([w_hat[2, 1], w_hat[0, 2], w_hat[1, 0]])

    w_norm = norm(w)

    if w_norm < SMALL_NUM:
        inv_ap = p
    else:

        inv_a = np.eye(3) - 0.5 * w_hat + ((2 * math.sin(w_norm)) - w_norm *
                            (1 + math.cos(w_norm))) / (2 * math.sin(w_norm) * 
                            w_norm * w_norm) * np.dot(w_hat, w_hat)
        inv_ap = np.dot(inv_a, p)

    return np.hstack([inv_ap, w])

#------------------------------------------------
def robot_ad_mat(num_j, t, s, q):
    """
        find the adjoint matrx of robot given the fixed frame 
        and joint positions
        --input:
                num_j:  number of joints
                t:      num_j * 4 * 4 matrix, initial robot configuration
                s:      num_j * 6 matrix, rotation/prasmatic joints axis
                q:      1 * num_j joints
        --output:
                tmp_a: 6*num_j adjoint matrix
    """
    tmp_a = adjoint_rep(t[0])
    for i in range(1, num_j):
        tmp_t = np.dot(t[0], twist_exp(s[0], q[0]))
        tmp = adjoint_rep(tmp_t)
        for j in range(1, i):
            tmp_t2 = np.dot(t[j], twist_exp(s[j], q[j]))
            tmp = np.dot(tmp, adjoint_rep(tmp_t2))
        tmp = np.dot(tmp, adjoint_rep(t[i]))
        tmp_a = np.hstack([tmp_a, tmp])

    return tmp_a


def robot_poe_cali(t_m, t0, s, q, beta=0.5):
    """
        the poe calibration framework
        --input:
                t_m: the measured pose of end effector, m*4*4 matrix
                t0: the robot inital configuration, num_j*4*4 matrix
                s:  the rotaion/translation axis,   num_j*6 matrix
                q:  joints angles, m*6/7 matrix
        --output:
                t_cali: calibrated robotic configuration
    """
    # number of measured data
    num_m = q.shape[0]

    # number of joints
    num_j = q.shape[1] #- 1

    # 7 instead of 6 because of the tool on the end-effector
    if q.shape[1] != 7:
        print("warrning: joint size is not 7, may cause error")

    ad = np.zeros((num_m*6, 6*num_j))

    dy = np.zeros(num_m*6)

    for i in range(num_m):

        t_t2b = poe_trans(t0, s, q[i])
        t_t2b_inv = np.linalg.inv(t_t2b)
        # print t_m[i]
        # print t_t2b_inv
        dy[6*i:6*(i+1)] = log_v(np.dot(t_m[i], t_t2b_inv))
        ad[6*i:6*(i+1)] = robot_ad_mat(num_j, t0, s, q[i])

    print "calibration: ad size: ", ad.shape
    
    # rcond is critical
    x = np.dot(np.linalg.pinv(ad, rcond=1e-12), dy)

    t_cali = np.copy(t0)

    for i in range(num_j):
        t_cali[i] = np.dot(t_cali[i], exp_se3(beta * x[6*i:6*(i+1)]) )

    return t_cali
# -----------------------------------------------------
def robot_poe_cali2(y_m, t0, s, q, beta=0.5):
    """
        the poe calibration framework
        --input:
                y_m: the measured pose of end effector, 
                t0: the robot inital configuration, num_j*4*4 matrix
                s:  the rotaion/translation axis,   num_j*6 matrix
                q:  joints angles, m*6/7 matrix
        --output:
                t_cali: calibrated robotic configuration
    """
    # number of measured data
    num_m = q.shape[0]

    # number of joints
    num_j = q.shape[1] #- 1

    if q.shape[1] != 7:
        print("warrning: joint size is not 7, may cause error")

    ad = np.zeros((num_m*6, 6*num_j))

    dy = np.zeros(num_m*6)

    for i in range(num_m):

        t_t2b = poe_trans(t0, s, q[i])
        t_t2b_inv = np.linalg.inv(t_t2b)
        # print t_m[i]
        # print t_t2b_inv
        t_m = exp_se3(y_m[i])
        dy[6*i:6*(i+1)] = log_v(np.dot(t_m, t_t2b_inv))
        ad[6*i:6*(i+1)] = robot_ad_mat(num_j, t0, s, q[i])

    print "calibration: ad size: ", ad.shape
    x = np.dot(np.linalg.pinv(ad, rcond=1e-14), dy)
    t_cali = np.copy(t0)

    for i in range(num_j):
        t_cali[i] = np.dot(t_cali[i], exp_se3(beta * x[6*i:6*(i+1)]) )

    return t_cali
# -----------------------------------------------------
def robot_tb_cali(t_m, t0, s, q, beta=0.5):
    """
        the poe calibration framework
        --input:
                t_m: the measured pose of end effector, m*4*4 matrix
                t0: the robot inital configuration, num_j*4*4 matrix
                s:  the rotaion/translation axis,   num_j*6 matrix
                q:  joints angles, m*6/7 matrix
        --output:
                t_cali: calibrated robotic configuration
    """
    # number of measured data
    num_m = q.shape[0]

    # number of joints
    num_j = q.shape[1] #- 1

    if q.shape[1] != 7:
        print("warrning: joint size is not 7, may cause error")

    ad = np.zeros((num_m*6, 6*num_j))

    dy = np.zeros(num_m*6)

    for i in range(num_m):

        t_t2b = poe_trans(t0, s, q[i])
        t_t2b_inv = np.linalg.inv(t_t2b)
        # print t_m[i]
        # print t_t2b_inv
        dy[6*i:6*(i+1)] = log_v(np.dot(t_m[i], t_t2b_inv))
        ad[6*i:6*(i+1)] = robot_ad_mat(num_j, t0, s, q[i])

    print "calibration: ad size: ", ad.shape
    x = np.dot(np.linalg.pinv(ad, rcond=1e-12), dy)
    t_cali = np.copy(t0)

    
    t_cali[0] = np.dot(t_cali[0], exp_se3(beta * x[0:6]) )
    # t_cali[6] = np.dot(t_cali[6], exp_se3(beta * x[36:]) )
    t_cali[6] = np.dot(t_cali[6], exp_se3(beta * x[36:]) )

    return t_cali
# -----------------------------------------------------
def robot_tb_cali2(y_m, t0, s, q, beta=0.5):
    """
        the poe calibration framework
        --input:
                y_m: the measured pose of end effector, m*4*4 matrix
                t0: the robot inital configuration, num_j*4*4 matrix
                s:  the rotaion/translation axis,   num_j*6 matrix
                q:  joints angles, m*6/7 matrix
        --output:
                t_cali: calibrated robotic configuration
    """
    # number of measured data
    num_m = q.shape[0]

    # number of joints
    num_j = q.shape[1] #- 1

    if q.shape[1] != 7:
        print("warrning: joint size is not 7, may cause error")

    ad = np.zeros((num_m*6, 6*num_j))

    dy = np.zeros(num_m*6)

    for i in range(num_m):

        t_t2b = poe_trans(t0, s, q[i])
        t_t2b_inv = np.linalg.inv(t_t2b)
        # print t_m[i]
        # print t_t2b_inv
        t_m = exp_se3(y_m[i])
        dy[6*i:6*(i+1)] = log_v(np.dot(t_m, t_t2b_inv))
        ad[6*i:6*(i+1)] = robot_ad_mat(num_j, t0, s, q[i])

    print "calibration: ad size: ", ad.shape
    x = np.dot(np.linalg.pinv(ad, rcond=1e-12), dy)
    t_cali = np.copy(t0)

    
    t_cali[0] = np.dot(t_cali[0], exp_se3(beta * x[0:6]) )
    # t_cali[6] = np.dot(t_cali[6], exp_se3(beta * x[36:]) )
    t_cali[6] = np.dot(t_cali[6], exp_se3(beta * x[6:12]) )

    return t_cali

#-------------------------------------------------------
def robot_poe_cali_load(t_m, t0, s, q, m_in, m, k, p, i_s, j_e, it=1, beta=0.5):
    """
        the poe calibration framework, NOT WORKING
        --input:
                t_m: the measured pose of end effector, m*4*4 matrix
                t0: the robot inital configuration, num_j*4*4 matrix
                s:  the rotaion/translation axis,   num_j*6 matrix
                q:  joints angles, m*6/7 matrix
        --output:
                t_cali: calibrated robotic configuration
    """
    
    # number of measured data
    num_m = q.shape[0]

    # number of joints
    num_j = q.shape[1] #- 1

    if q.shape[1] != 7:
        print("warrning: joint size is not 7, may cause error")

    ad = np.zeros((num_m*6, 6*num_j))

    dy = np.zeros(num_m*6)

    for i in range(num_m):
        m[-1] = m_in[i]
        # t_t2b = poe_trans(t0, s, q[i])
        t_t2b = poe_trans_load_ij(t0, s, q[i], m, k, p, i_s, j_e, it=1)
        t_t2b_inv = np.linalg.inv(t_t2b)
        # print t_m[i]
        # print t_t2b_inv
        dy[6*i:6*(i+1)] = log_v(np.dot(t_m[i], t_t2b_inv))
        ad[6*i:6*(i+1)] = robot_ad_mat(num_j, t0, s, q[i])

    print "calibration: ad size: ", ad.shape
    x = np.dot(np.linalg.pinv(ad, rcond=1e-14), dy)
    t_m = np.copy(t0)

    for i in range(num_j):
        t_m[i] = np.dot(t_m[i], exp_se3(beta * x[6*i:6*(i+1)]) )

    return t_m

# ------------------------------------------------------ 
def poe_dy(t0,s,q,t_m):
    """
        compute the difference of end effector of measured results 
        and kinematci results
    """
    if q.shape[1] != 7:
        print("warrning: joint size is not 7, may cause error")

    num_m = q.shape[0]
    # dy = np.zeros(num_m*6)
    dy = np.zeros((num_m, 6))
    for i in range(num_m):
        # print t0,q[i]
        t_t2b = poe_trans(t0, s, q[i])
        # print t_t2b
        t_t2b_inv = np.linalg.inv(t_t2b)
        # dy[6*i:6*(i+1)] = log_v(np.dot(t_m[i], t_t2b_inv))
        dy[i] = log_v(np.dot(t_m[i], t_t2b_inv))
    return dy

def poe_dy2(t0,s,q,y_m):
    """
        compute the difference of end effector of measured results 
        and kinematci results
    """
    if q.shape[1] != 7:
        print("warrning: joint size is not 7, may cause error")
    num_m = q.shape[0]
    dy = np.zeros((num_m, 6))
    for i in range(num_m):
        t_t2b = poe_trans(t0, s, q[i])
        t_t2b_inv = np.linalg.inv(t_t2b)
        t_m = exp_se3(y_m[i])
        # print t_m, t_t2b
        dy[i] = log_v(np.dot(t_m, t_t2b_inv))
        # print dy[i]
    return dy

# ------------------------------------------------------   

class RobotPoe:
    """
        ## TODO: put calibrations in the class
    """
    pass

# ===============================================
if __name__ == "__main__":
    print "poe calibration module"
    print adjoint_rep(np.eye(4))
    print log_v(np.eye(4))
    log_v_test = np.loadtxt("test.txt")
    print log_v_test
    print log_v(log_v_test)
