from __future__ import division
import numpy as np
import math


# import modules from the same package
from homo_trans import *

G = np.array([0,0,-9.81])               # gravity constant
SMALL_NUM = 0.0000001   # small constant    
#==================================================


def norm(x):
    return np.linalg.norm(x)
#------------------------------------------------


def skew_sym(w):
    """
        input: 
                1*3 numpy array contains rotation angles
        output:
                3*3 numpy array skew symmetric matrix
    """
    res = np.array(
        [[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]])
    return res


def skew_exp(w, a):
    """
        input: 
            w: 1*3 numpy array contains rotation angles
            a: angle 
        return 3*3 skew exponential matrix
    """
    # print w,a 
    w2 = skew_sym(w)
    res = np.eye(3) + math.sin(a) * w2 + np.dot(w2, w2) * (1 - math.cos(a))
    return res


def twist_exp(rp, a):
    """
        the twist exponential, returns a 4*4 homogenous transformation matrix
    """
    w = rp[3:6]
    v = rp[0:3]
    r = skew_exp(w, a)
    p = a * v
    # print r,p
    return rp_homo4(r, p)

def exp_se3(x):
    """
        exponential of special eulcilean group 3
        input:
                1*6 np array
        output
                4*4 transformation matrix
    """
    SMALL_NUM = 0.0000001
    i = np.eye(3)
    v = x[:3]
    w = x[3:]

    w_norm = norm(w)
    w_hat = skew_sym(w)

    if math.fabs(w_norm) > SMALL_NUM:
        a = i + (1-math.cos(w_norm))/(w_norm**2) * w_hat \
                + (w_norm - math.sin(w_norm)) \
                /(w_norm ** 3) * np.dot(w_hat, w_hat)
        r = i + (math.sin(w_norm)/w_norm) * w_hat \
                + (1-math.cos(w_norm))/(w_norm ** 2) * np.dot(w_hat, w_hat)
    else:
        a = i
        r = i
    p = np.dot(a, v)
    p = p[:, None]

    res = np.vstack([np.hstack([r, p]),
                     np.hstack([np.zeros((3)), np.ones(1)])])
    return res    

# ----------------------------------------------------

# def adjoint_rep(m):
#     """
#         the adjoint representation of matrix, refer to lie algebra
#          defined in poe_calibration module, to delete
#     """
#     pass

# ----------------------------------------------------
def poe_trans(t0, s, q):
    """
        compute the overall transformation matrix of poe
        INPUT:
            s: the configuration of each joints
            q: the joints angles        
    """
    t = np.eye(4)
    # num_q = size(q)
    num_q = q.shape[0]
    if num_q != 7:
        print("warnning, robot model may not be right")

    for i in range(num_q):
        t = np.dot(np.dot(t, t0[i]), twist_exp(s[i], q[i]))
    return t
# --------------------------------------------------


def robot_kinematics(t0, s, q):
    """
        robot kinematics, given the joint rotations, calculate
        the 4*4 tool-base transformation matrix 
        INPUT:
            s: the configuration of each joints
            q: the joints angles
    """
    return poe_trans(t0,s,q)

# --------------------------------------------------
def poe_trans_load_ij(t0, s, q, m, k, p, i_s, j_e, it=1):
    """
        load transformation, with starting and ending joints, 
        ignore this for normal calibration
        INPUT:
            t0,s,q:     the same as non load problem
            i_s, j_e:   int, starting and ending index of joints
            m:          mass
            k:          inertia
        OUTPUT:
            res:        the resulting 4*4 transformation matrix
    """
    t_ij = np.eye(4)
    if q.shape[0] != 7:
        print("warrning: joint size is not 7, may cause error")
    if it <= 0:
        for i in range(i_s, j_e+1):
            t_ij = np.dot(np.dot(t_ij, t0[i]), twist_exp(s[i], q[i]))
    else:
        it = it -1
        for i in range(i_s,j_e+1):
            p_cg_in = m[i] * p[i]
            sum_m = m[i]

            for j in range(i+1,j_e+1):
                p_cg_in = p_cg_in + m[j] * np.dot(poe_trans_load_ij(t0, s, q, m, 
                            k, p, i_s+1, j_e, it=it) , p[j])
                sum_m = sum_m + m[j]

            if sum_m < SMALL_NUM:
                p_cg_in = 1.0 * np.array([0,0,0,1])
            else:
                p_cg_in = p_cg_in / sum_m

            p_cg_in = p_cg_in[:3]

            rt = poe_trans_load_ij(t0, s, q, m, k, p, 0, i, it=it)
            # print rt
            G_i = np.dot(homo4_r3(rt), G)
            # print "Gi:", G_i
            q_tilde = k[i] * sum_m * np.inner(s[i,3:], np.cross(p_cg_in, G_i))
            # print "qtilde:", q_tilde
            t_temp =  np.dot(np.dot(t_ij, t0[i]), twist_exp(s[i], q[i]))
            t_ij = np.dot(t_temp, twist_exp(s[i], q_tilde))
        # print "test"
    return t_ij



def fk_poe_load(t0, s, q, m, k, p):
    """
        compute the overall transformation matrix of poe, with
        consideration of mass of each link
    """
    j_e = 5
    temp = poe_trans_load_ij(t0, s, q, m, 
                        k, p, 0, j_e)
    t_e = np.dot(temp, t0[j_e+1])
    return t_e

def robot_kine_load():
    pass

# ================================================
if __name__ == "__main__":
    print(__file__)
    # print rot3_x(math.pi/2)
    # print skew_sym(np.array([1, 2, 3]))
    # print norm(np.array([1, 2, 3]))

    # # exp_se3 passed sanaty test of all 1s and 0s
    # print exp_se3(np.ones(6))  

         
