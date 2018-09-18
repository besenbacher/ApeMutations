
    
class family():
    def __init__(self, child, father, mother, child_gender, grandchildren):
        self.child = child
        self.gender = child_gender
        self.father = father
        self.mother = mother
        self.grandchildren = grandchildren
    def __iter__(self):
        yield self.child
        yield self.father
        yield self.mother
        for x in self.grandchildren:
            yield x
    def trio(self):
        yield self.father
        yield self.mother
        yield self.child

def read_families(fname):
    f = open(fname)
    L = []
    for line in f:
        child, father, mother, gender, grandchildren = line.split()
        if grandchildren == "NA":
            grandchildren = []
        else:
            grandchildren = grandchildren.split(",")
        L.append(family(child, father, mother, gender, grandchildren))
    f.close()
    return L

def replace_empty(L):
    if L==None or len(L) == 0:
        return "PASS"
    else: 
        return "_".join(L)

