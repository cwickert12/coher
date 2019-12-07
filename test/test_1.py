import unittest
#import sys
#sys.path.append('../')
#from acon2 import *
import numpy as np
 
class Test_1(unittest.TestCase):
 
    def test1(self):


        x = [5,7]
        y = [5.00000000001,7]


        np.testing.assert_almost_equal( x,y,5)
        np.testing.assert_almost_equal( x[0],y[0],5)

 


    
if __name__ == '__main__':
    unittest.main()
