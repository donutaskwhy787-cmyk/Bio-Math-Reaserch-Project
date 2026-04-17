import numpy as np
#idh 1 = mutant
class Participant:
  def __init__(self):
    self.values = dict()
  def initalize(self, headers, values):
    hds = headers.split(" ")
    vals = values.split(" ")
    self.values = dict()

    for i in range(len(hds)):
      self.values[hds[i]] = vals[i]

  def deepcopy(self):
    c = dict()
    p = Participant()
    for i in self.values:
      c[i] = self.values[i]
    p.values = c
    return p

    
