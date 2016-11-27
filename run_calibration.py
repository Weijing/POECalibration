#!/usr/bin/python
"""
    this is used to calibration robotic arm with Local POE method

    usage: this is meant to be a standalone module, run the following command:
            >>> python run_calibration.py -h 
            for help

    Author:        Wei Jing
    Data:           Nov 2016
"""
from __future__ import division
import numpy as np
from time import gmtime, strftime
import argparse
import math
import sys, os
import scipy.io as sio

sys.path.append("./poe")

from homo_trans import *
from load_config import *
from poe_model import *
from poe_calibration import *
# ============================================================================
def error_check(q, t_m, fk, s):
    """ error 
    
    [error in se(3)]
    
    Arguments:
        q {[type]} -- [description]
        t_m {[type]} -- [description]
        fk {[type]} -- [description]
        s {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    error = poe_dy(fk, s, q, t_m)
    return error

def calibration_poe(q, t_m, fk_init, s, iter=10, step=0.3):
    fk_cali = np.copy(fk_init)
    fk_cali = robot_poe_cali(t_m, fk_init, s, q[:,:], beta=0.1)
    for i in range(15):
        print i
        fk_cali = robot_poe_cali(t_m, fk_cali, s, q[:,:])
    return fk_cali

def generate_data():
    parser = argparse.ArgumentParser(description='Calibration using Local POE')

    parser.add_argument('-m', dest="file_mat", default="data.mat",
                        help="input file, matlab .mat")

    # parser.add_argument('-n', dest="N",  type=int, default=300,
    #                     help="number of data to generate")
    # parser.add_argument('-r', dest="r",  type=float, 
    #                 default=math.pi/2, help="range of data")

    parser.add_argument('-f', dest="data_folder",  type=str, 
                    default="abb-irb6600", help="choose the folder ")

    print ("parameters: ")
    print( parser.parse_args())
    in_arg = parser.parse_args()

    crt_path = os.path.dirname(os.path.realpath(__file__))

    data_path = os.path.join(crt_path, in_arg.data_folder)
    print(data_path)
    try:
        file_s = os.path.join(data_path, "s_vector.txt")

        file_fk_init = os.path.join(data_path, "fk_init.txt")
        file_q = os.path.join(data_path, "joints.txt")
        file_t_m = os.path.join(data_path, "t_measured.txt")

        s = load_robot_config_s(file_s)
        fk_init = load_robot_config_t0(file_fk_init)
        q = np.loadtxt(file_q)
        t_m0 = np.loadtxt(file_t_m)

    except Exception as e:
        raise 
    num = q.shape[0]
    t_m = t_m0.reshape(num,4,4)
    fk_cali = calibration_poe(q, t_m, fk_init, s)

    np.set_printoptions(formatter={'float_kind':'{:.5f}'.format})
    print fk_cali

    print error_check(q, t_m, fk_cali, s)

# =============================================================================
if __name__ == "__main__":

    generate_data()