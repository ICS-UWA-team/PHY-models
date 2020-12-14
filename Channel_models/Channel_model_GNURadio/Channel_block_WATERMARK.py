"""
Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
"""
# ========================================================================================
# >> IMPORTS
# ========================================================================================
import numpy as np
from collections import deque
from gnuradio import gr
import scipy.io
import scipy.signal
import math

# ========================================================================================
# >> my_test_block
# ========================================================================================
class blk(gr.sync_block):  # other base classes are basic_block, decim_block, interp_block
    """Embedded Python Block example - a simple multiply const"""

    def __init__(self, sampling_frequency=6e5, file_path='/home/froztgal/BCH1_001.mat'):  # only default arguments here
        """arguments to this function show up as parameters in GRC"""
        gr.sync_block.__init__(
            self,
            name='Channel_block_WATERMARK',   # will show up in GRC
            in_sig=[np.complex64],
            out_sig=[np.complex64]
        )
        # if an attribute with the same name as a parameter is found,
        # a callback is registered (properties work, too).
        self.fc = sampling_frequency
        self.read_file(file_path)
        self.hsamp = self.h.shape[0]
        self.upsamp = math.ceil(self.fc / self.fs_tau)
        self.buffer = deque([0] * self.upsamp)
        self.counters_limits = [(self.hsamp * self.upsamp), self.h.shape[1]]
        print(self.counters_limits)
        self.coeff = [x / self.counters_limits[0] for x in range(1, int(self.counters_limits[0]))]
        self.counter1 = 0
        self.counter2 = 1
        
    def shift(self, new_element):
        self.buffer.popleft()
        self.buffer.append(new_element)
        
    def read_file(self, h_file):
        mat = scipy.io.loadmat(h_file)
        self.fs_tau = float(mat['fs_tau'][0][0])
        self.h = mat['h'].T

    def work(self, input_items, output_items):
        # Beginning of signal processing
        
        signal_out = np.zeros(len(input_items[0]), dtype=complex)
        for i in range(len(input_items[0])):
                      
            # Send new input sample in buffer
            self.shift(input_items[0][i])

            # Preparing and resampling of h
            col1 = self.h[:, self.counter2 - 1]
            htmp1 = (1 - self.coeff[self.counter1]) * col1
            col2 = self.h[:, self.counter2]
            htmp2 = self.coeff[self.counter1] * col2
            htmp = htmp1 + htmp2
            tmp = scipy.signal.resample(htmp, self.upsamp)
            

            # Convolution
            signal_out[i] = scipy.signal.fftconvolve(self.buffer, tmp, 'valid')

            # Checking counters
            # Counter for interpolation h (sample time = 1/Fs)
            if self.counter1 == self.counters_limits[0] - 2:
                self.counter1 = 0
                # Counter for reading h (sample time = dt)
                if self.counter2 == self.counters_limits[1]:
                    self.counter2 = 0
                else:
                    self.counter2 += 1
            else:
                self.counter1 += 1
                           
        # Ending of signal processing
        output_items[0][:] = signal_out
        return len(output_items[0])
       
       
       
