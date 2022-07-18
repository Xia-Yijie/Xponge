import numpy as np


class mdout():
    def __getattribute__(self, attr):
        if attr not in ("content", "content_index", "data") and attr in self.content_index.keys():
            return self.data[:, self.content_index[attr]]
        else:
            return super().__getattribute__(attr)

    def __init__(self, filename):
        with open(filename) as f:
            self.content = f.readline().split()
            self.data = np.loadtxt(f)
        self.content_index = {self.content[i]: i for i in range(len(self.content))}
