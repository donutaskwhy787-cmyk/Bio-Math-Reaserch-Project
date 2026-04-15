import numpy as numpy

class EncodedCases:
  #Use raw molecular call, ignore WHO text label if conflicting
  #Exclude sample if missing (required for subtype classification)
  def IDHMutationStatus(self, par):
    print("")

  #Binarize from beta value: threshold ≥ 0.30 = methylated
  #If MGMT data is missing, fill it in with whatever value appears most in the cohort, 
  #and mark those samples so we know they were filled in and not real measurements
  def MGMTMethylation(self, par):
    print("")   

  #Log2(FPKM+1) normalize; amplified = CNV≥2 OR FPKM ≥ 90th pct
  #Exclude sample from transcriptomic node if RNA-seq missing
  def EGFRExpression(self, par):
    print("")

  #Convert days to months if needed (÷30.44), floor at 0
  #Exclude if OS completely absent; note censored vs. deceased separately
  def overallSurvival(self, par):
    print("")

  #Standardize: 1=Deceased, 0=Living; resolve string variants
  #Exclude if missing (needed for survival interpretation)
  def vitalStatus(self, par):
    print("")
  
  #Map numeric: 2→2, 3→3, 4→4; cross check IDH for label consistency
  #Use IDH status to infer if WHO label absent or ambiguous
  def WHOGrade(self, par):
    print("")

