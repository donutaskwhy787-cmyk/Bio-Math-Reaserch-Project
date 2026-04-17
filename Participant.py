import numpy as np
#idh 1 = mutant
class Participant:
  def __init__(self, headers, values):
    hds = headers.split(" ")
    vals = values.split(" ")
    self.values = dict()

    for i in range(len(hds)):
      self.values[hds[i]] = vals[i]


    
