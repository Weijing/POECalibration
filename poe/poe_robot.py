from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import sys, os

from sklearn.gaussian_process import GaussianProcess
from sklearn.cross_validation import KFold

from homo_trans import *
from load_config import *
from poe_model import *
from poe_calibration import *
from poe_ik import poe_ik_gp2, poe_ik, poe_ik_load

class RobotPOE(object):
    """[summary]
    
    [description]
    
    Variables:
        print("test") {[type]} -- [description]
    """

    def __init__(self, file_ta=None, file_load=None):
        self.NUM_JOINTS = 7
        fn_pwd = os.path.dirname(__file__)
        # print fn_pwd
        fn_root = os.path.join(os.path.dirname(__file__), os.pardir)
        # print fn_root

        file_robot_conf_s = "../robot_config/sVector2.txt"
        file_robot_conf_t = "../robot_config/nominalMatrices2.txt"
        if file_ta is None:
            file_ta = "../robot_config/TPerturbed4.txt"

        self.s = load_robot_config_s(os.path.join(fn_pwd, file_robot_conf_s))
        self.t0 = load_robot_config_t0(os.path.join(fn_pwd,file_robot_conf_t))
        self.t_a = load_robot_config_t0(os.path.join(fn_pwd,file_ta))

        if file_load is None:
            file_load = "../robot_config/loading_model.mat"

        self.t_cali = None
        self.gp = None
        self.gp_ik = None

        data = sio.loadmat(os.path.join(fn_pwd,file_load))
        self.k = data["K"][0]
        self.p = data["P"]
        self.mass = data['mass'].T
        self.m = 1.0*(data["M"][0])

        # immutable attributes
        self.m.setflags(write=False)
        self.k.setflags(write=False)
        self.p.setflags(write=False)
        self.s.setflags(write=False)
        self.t0.setflags(write=False)



    def fk(self, q, t=None):
        """Forward kinematics for 1 joint set input
            :param q: joint angles, 1d array, size=self.NUM_JOINTS
        """
        if q.shape[0] != self.NUM_JOINTS:
            raise ValueError
        if t is None:
            t = self.t0
        return poe_trans(t, self.s, q)

    def fk_load(self, q, m_load=None):
        """[summary]
        
        [description]
        
        Arguments:
            q {[np array]} -- [joint angles]
        
        Keyword Arguments:
            m_load {[float]} -- [description] (default: {None})
        
        Returns:
            [type] -- [description]
        """

        if m_load is None:
            m_temp = self.m
        else:
            m_temp = np.zeros(5)
            m_temp[:4] = self.m[:4]
            m_temp[-1] = m_load
        
        return fk_poe_load(self.t_a, self.s, q, m_temp, self.k, self.p)

    def ik_unload(self, p_e, q0, t=None):

        if t is None:
            t = self.t0 

        return poe_ik(q0, p_e, t, self.s)


    def ik_gp(self, p_e, q0, gp, t=None):
        if t is None:
            t = self.t0
        return poe_ik_gp2(q0, p_e, t, self.s, gp)

    def ik_load(self,p_e, q0, t=None, m_load=None):   
        if t is None:
            t = self.t_a
        if m_load is None:
            m_temp = self.m
        else:
            m_temp = np.zeros(5)
            m_temp[:4] = self.m[:4]
            m_temp[-1] = m_load
        return poe_ik_load(q0, p_e, t, self.s, m_temp, self.k, self.p)

    def calibrate(self, q, dy, beta=0.5, iter=10):

        # don't recalibrate
        if self.t_cali is not None:
            return

        self.t_cali = np.copy(self.t0)
        
        if dy.shape[1] == 6:            
            for _ in range(iter):
                self.t_cali = robot_poe_cali2(dy, self.t_cali, 
                                        self.s, q, beta=0.2)
            return
        # in case of 2d t_m
        if dy.shape[1] == 4:
            print ("warning: the measured input to calibration is 2D")            
            for _ in range(iter):
                self.t_cali = robot_poe_cali(dy, self.t_cali, 
                                        self.s, q, beta=0.2)
            return
    

    def gp_fk_train(self):
        pass

    def gp_ik_train(self):



# ====================================================
if __name__ == '__main__':
    print("test")
