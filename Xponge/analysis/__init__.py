"""
This **module** gives functions and classes for analysis
"""
import numpy as np


class MdoutReader():
    """
This **class** is used to read the mdout file generated by SPONGE
    """
    def __init__(self, filename):
        with open(filename) as f:
            self.content = f.readline().split()
            self.data = np.loadtxt(f)
        self.content_index = {self.content[i]: i for i in range(len(self.content))}

    def __getattribute__(self, attr):
        if attr not in ("content", "content_index", "data") and attr in self.content_index.keys():
            return self.data[:, self.content_index[attr]]
        return super().__getattribute__(attr)
