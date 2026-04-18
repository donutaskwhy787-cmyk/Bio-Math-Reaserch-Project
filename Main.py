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

  def test(self, a, b, c, d, e):
        self.values = dict()
        self.values["IDH"] = a
        self.values["EGFR"] = c
        self.values["MGMT"] = b
        self.values["GRADE"] = e
        self.values["MONTHS"] = d


class Consistancy:
    def rDR(self, idh):
        return 1 - idh
    def rDC(self, idh, mgmt):
        if(idh == 0 or mgmt == 0):
            return 1
        else:
            return 0
    def rRC(self, egfr):
        return egfr
    def severity(self, grade, months):
        if(grade == 1):
            return 1
        if(months < 24):
            return 1
        return 0
    def eDR(self, idh, egfr):
        return abs(self.rDR(idh) - egfr)
    def eDC(self, idh, mgmt, grade, months):
        return abs(self.rDC(idh, mgmt) - self.severity(grade, months))
    def eRC(self, egfr, grade, months):
        return abs(self.rRC(egfr) - self.severity(grade, months))
    
    def consistancyVector(self, p):
        vals = p.vals
        idh = vals["IDH"] 
        egfr = vals["EGFR"]
        mgmt = vals["MGMT"]
        grade = vals["GRADE"]
        surv = vals["MONTHS"]
        eDR = self.eDR(idh, egfr)
        eDC = self.eDC(idh, mgmt, grade, surv)
        eRC = self.eRC(egfr, grade, surv)

        return str(eDR) + " " + str(eDC) + " " + str(eRC)
    
    def consistancy(self, p):
      str = self.consistancyVector(p)
      sum = 0
      for i in str:
         if(i == '1'):
            sum += 1
      return sum/3


class EncodedCases:
  #Use raw molecular call, ignore WHO text label if conflicting
  #Exclude sample if missing (required for subtype classification)
  def IDHMutationStatus(self, par):
    if(par.values.get("IDH") is None):
      return None
    p = par.deepcopy()
    
    if(par.values["IDH"] == "mutant"):
      p.values["IDH"] = 1
    else:
      p.values["IDH"] = 0
    return p

  #Binarize from beta value: threshold ≥ 0.30 = methylated
  #If MGMT data is missing, fill it in with whatever value appears most in the cohort, 
  #and mark those samples so we know they were filled in and not real measurements
  def MGMTMethylation(self, par):
    p = par.deepcopy()
    if(par.values.get("MGMT") is None):
      print("remind michael that he still has to do this part, and is just waiting to see what the data looks like")
    elif(par.values["MGMT"] >= 0.3):
      p.values["MGMT"] = 1
    else:
      p.values["MGMT"] = 0
    return p

  #Log2(FPKM+1) normalize; amplified = CNV≥2 OR FPKM ≥ 90th pct
  #Exclude sample from transcriptomic node if RNA-seq missing
  def EGFRExpression(self, par):
    print("")

  #Convert days to months if needed (÷30.44), floor at 0
  #Exclude if OS completely absent; note censored vs. deceased separately
  def overallSurvival(self, par):
    p = par.deepcopy()
    if(p.values.get("MONTHS")is None):
        return None
    if(p.values["MONTHS"] < 24):
       p.values["MONTHS"] = 0
    else:
       p.values["MONTHS"] = 1
    return p
    

  #Standardize: 1=Deceased, 0=Living; resolve string variants
  #Exclude if missing (needed for survival interpretation)
  def vitalStatus(self, par):
    p = par.deepcopy()
    if(p.values.get("STATUS") is None):
       return None
    print(p.values["STATUS"] == 0)
    if(p.values["MONTHS"] == 0 and p.values["STATUS"] == "0"):
       return None
    return p
  
  #Map numeric: 2→2, 3→3, 4→4; cross check IDH for label consistency
  #Use IDH status to infer if WHO label absent or ambiguous
  def WHOGrade(self, par):
    print("")

class OTOutputs:
   def cost(self, i, j):
      print("AHHH")

    

p = Participant()
p.test("mutant", None, 0, 0, 0)
p.values["STATUS"] = "0"
ec = EncodedCases()

print(ec.vitalStatus(p))
