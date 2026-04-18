class EncodedCases:
  #Use raw molecular call, ignore WHO text label if conflicting
  #Exclude sample if missing (required for subtype classification)
  def IDHMutationStatus(self, par):
    if(par.values["IDH"] == None):
      return None
    p = par.deepcopy()
    if(par.values["IDH"] == "mutant"):
      p.values["IDH"] = 1
    else:
      p.values["IDH"] = 0

  #Binarize from beta value: threshold ≥ 0.30 = methylated
  #If MGMT data is missing, fill it in with whatever value appears most in the cohort, 
  #and mark those samples so we know they were filled in and not real measurements
  def MGMTMethylation(self, par):
    p = par.deepcopy()
    if(par.values["MGMPT"] == None):
      print("remind michael that he still has to do this part, and is just waiting to see what the data looks like")
    elif(par.values["MGMT"] >= 0):
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
    print("")

  #Standardize: 1=Deceased, 0=Living; resolve string variants
  #Exclude if missing (needed for survival interpretation)
  def vitalStatus(self, par):
    print("")
  
  #Map numeric: 2→2, 3→3, 4→4; cross check IDH for label consistency
  #Use IDH status to infer if WHO label absent or ambiguous
  def WHOGrade(self, par):
    print("")
