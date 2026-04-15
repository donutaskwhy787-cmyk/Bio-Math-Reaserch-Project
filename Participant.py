import numpy as np

class Participant:
  def __init__(self, headers, values):
    hds = headers.split(" ")
    vals = values.split(" ")
    self.values = dict()

    for i in range(len(hds)):
      self.values[hds[i]] = vals[i]


    
