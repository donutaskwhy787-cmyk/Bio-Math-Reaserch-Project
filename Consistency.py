
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
    def severity(self, grade, survival):
        if(grade == 1):
            return 1
        if(survival < 24):
            return 1
        return 0
    def eDR(self, idh, egfr):
        return abs(self.rDR(idh) - egfr)
    def eDC(self, idh, mgmt, grade, survival):
        return abs(self.rDC(idh, mgmt) - self.severity(grade, survival))
    def eRC(self, egfr, grade, survival):
        return abs(self.rRC(egfr) - self.severity(grade, survival))
    
    def consistancy(self, p):
        vals = p.vals
        idh = vals["IDH"] 
        egfr = vals["EGFR"]
        mgmt = vals["MGMT"]
        grade = vals["GRADE"]
        surv = vals["SURVIVAL"]
        eDR = self.eDR(idh, egfr)
        eDC = self.eDC(idh, mgmt, grade, surv)
        eRC = self.eRC(egfr, grade, surv)

        print(eDR)
        print(eDC)
        print(eRC)
        return (eDR + eDC + eRC)/3.0
    
class tester:
    def __init__(self, a, b, c, d, e):
        self.vals = dict()
        self.vals["IDH"] = a
        self.vals["EGFR"] = c
        self.vals["MGMT"] = b
        self.vals["GRADE"] = e
        self.vals["SURVIVAL"] = d

