
import numpy as np 
import sys

def load_robot_config_s(filename):
    """
        load robotic configuration: rotation axis
    """
    data = np.loadtxt(filename)
    n_row = data.shape[0]
    if n_row % 6 != 0:
        print("Wrong format of s configuration")
        sys.exit()
    num_j = int(n_row/6)
    s = np.zeros((num_j,6))
    for i in range(num_j):
        s[i] = data[i*6:(i+1)*6]
    return s

def load_robot_config_t0(filename):
    """
        load robotic configuration: POE T{0}
    """
    data = np.loadtxt(filename)
    n_row = data.shape[0]
    if n_row % 4 != 0:
        print("Wrong format of s configuration")
        sys.exit()
    
    num_j = int(n_row/4)
    t = np.zeros((num_j,4,4))
    for i in range(num_j):
        t[i] = data[i*4:(i+1)*4,:]
    return t

def read_tr_data(filename):
    pass

def read_tt_data(filename):
    pass

def save_res(filename):
    pass

class CaliLog(object):

    def __init__(self, filename):
        try:
            self.f = open(filename, "w")
        except:
            print("Open log file failed")

    def write_log(self, s):
        self.f.write(s)
        self.f.write('\n')

    def close_log(self):
        self.f.close()



# ==============================================
if __name__ == "__main__":
    log1 = CaliLog("test.txt")
    print(log1.f)
    
    a = np.array([1, 2, 3])
    log1.write_log(str(a))
    log1.close_log()

