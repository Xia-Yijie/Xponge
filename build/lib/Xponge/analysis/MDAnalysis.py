try:
    import MDAnalysis as mda
except:
    raise Exception("'MDAnalysis' package needed. Maybe you need 'pip install MDANalysis'")

from MDAnalysis.coordinates import  base
from MDAnalysis.lib import util
import os.path 
import numpy as np
class SPONGE_Trajectory_Reader(base.ReaderBase):
    def __init__(self, dat_file_name, box, n_atoms, **kwargs):
        super().__init__(dat_file_name, **kwargs)
        if type(box) == type(""):
            self.boxname = box
            self.box = None
            self._get_box_offset()
        else:
            self.boxname = None
            self.box = box
        self._n_atoms = n_atoms
        self._n_frames = os.path.getsize(dat_file_name)//12//self.n_atoms
        self.trajfile = None
        self.boxfile = None
        self.ts = base.Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()        

    def _read_frame(self, frame):
        if self.trajfile is None:
            self.open_trajectory()
        if self.boxfile is not None:
            self.boxfile.seek(self._offsets[frame])
        self.trajfile.seek(self.n_atoms * 12 * frame)
        self.ts.frame = frame - 1
        return self._read_next_timestep()

    def _read_next_timestep(self):
        ts = self.ts
        if self.trajfile is None:
            self.open_trajectory()
        t = self.trajfile.read(12 * self.n_atoms)
        if not t:
            raise EOFError
        ts._pos = np.frombuffer(t, dtype=np.float32).reshape(self.n_atoms,3)
        
        if self.box is not None:
            ts.dimensions = self.box
        else:
            ts.dimensions = list(map(float, self.boxfile.readline().split()))
        ts.frame += 1
        return ts

    @property
    def n_frames(self):
        return self._n_frames

    @property
    def n_atoms(self):
        return self._n_atoms
    
    def close(self):
        if self.trajfile is not None:
            self.trajfile.close()
            self.trajfile = None
        if self.boxfile is not None:
            self.boxfile.close()
            self.boxfile = None

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        self.trajfile = util.anyopen(self.filename, "rb")
        if self.box is None:
            self.boxfile = util.anyopen(self.boxname)
        ts = self.ts
        ts.frame = -1
        return self.trajfile, self.boxfile

    def _get_box_offset(self):
        self._offsets = [0]
        with util.openany(self.boxname) as f:
            line = f.readline()
            while line:
                self._offsets.append(f.tell())
                line = f.readline()
        self._offsets.pop()

def Universe(topname, trajname, box, **kwargs):
    import warnings
    warnings.filterwarnings("ignore", message='No coordinate reader found for')
    u = mda.Universe(topname)
    warnings.filterwarnings("default", message='No coordinate reader found for')
    u._trajectory = SPONGE_Trajectory_Reader(trajname, box, len(u.atoms), **kwargs)
    return u
   

