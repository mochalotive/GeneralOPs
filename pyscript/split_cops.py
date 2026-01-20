import sys
from typing import Union, Dict, List

#typedefs
Number = Union[float, int] #for the semi-likely event this needs complex numbers later
Value = Union[Number, str] #mostly for grid

class CopFile:
    def __init__(self):
        self.tags : List[str] = [] #trying this for debugging :^)
        self.data : Dict[str, List[Value] ] = {} #eg !BOND 1,2,3
    def __eq__(self, other):
        if (self.tags != other.tags):
            return False
        return True
        #barebones check if they're at all compatible
        #doesn't really check if all the members are compatible
    def __add__(self, other):
          #not changing switch width, cutoff, or grid!
        for i, tag in enumerate(self.tags):
            #finitely many exceptions
            if tag  == "SWITCH_WIDTH":
                continue
            if tag == "CUTOFF":
                continue
            if tag == "GRID":
                continue
            if tag == "NCELLS":
                self.data[tag][0]+= other.data[tag][0]
                continue
            self.data[tag].extend(other.data[tag]) 
        #now fix grid, this is lazy, not the real grid
        self.data["GRID"] = [self.data["NCELLS"][0], 1, 1, 0, 0, 0]
        return self         
        
    def to_string(self):
        myStr: str  = ""
        for i, tag in enumerate(self.tags):
            myStr+=(f"!{tag}: \n")
            if tag == "GRID":
                for j, num in enumerate(self.data[tag]):
                    myStr+=(f"{num} ")
                myStr+=(f"\n")
                continue
            for j, num in enumerate(self.data[tag]):
                if (i == len(self.tags)-1 and j == len(self.data[tag])-1):
                    myStr+=(f"{num} ")  
                    continue
                myStr+=(f"{num} \n")  
        return myStr
    @classmethod
    def open(cls, filename) -> str:
        obj = cls() #instance of my class
        current_section = None
        with open(filename) as f:
        
            for line in f:
                line = line.strip()
                if line.startswith("!"):
                    section = line.rstrip(":").lstrip("!")
                    current_section = section
                    obj.tags.append(section)
                    obj.data[current_section] = []
                    continue #move to next line
                temp = line.split()
                for chunk in temp:
                    try:
                        if "." in chunk or "e" in chunk.lower():
                            addme = float(chunk)
                        else:
                            addme = int(chunk)
                    except ValueError:
                        addme = chunk #string fallback
                    obj.data[current_section].append(addme)
        return obj #return pointer to my class instance
    #setters and getters
    @property
    def switch_width(self):
        return float ( self.data["SWITCH_WIDTH"][0] )
    @switch_width.setter
    def switch_width(self, w: float):
        try:
            assert w >= 0, "switch width must be geq zero."
            self.data["SWITCH_WIDTH"] = [float(w)]
        except TypeError:
            print("expecting a positive, real number for switch_width")
        
    @property
    def cutoff(self):
        return float( self.data["CUTOFF"][0])
    @cutoff.setter
    def cutoff(self, w: float):
        try:
            assert w > 0, "cutoff must be greater than zero."
            self.data["CUTOFF"] = [float(w)]
        except TypeError:
            print("expecting a positive, real number for cutoff")

    @property 
    def grid(self):
        temp = self.data["GRID"][0]
        temp = temp.split()
        return [int(line) for line in temp]
    @property
    def ncells(self):
       return int( self.data["NCELLS"][0]) 
    @ncells.setter
    def ncells(self, n: int):
        try:
            assert n >= 0, "ncells must be geq than zero."
            self.data["NCELLS"] = [int(n)]
        except TypeError:
            print("expecting natural number for ncells")
       
ops = ["DISTANCE_OPS", "BOND_ORIENTATION_OPS" , "RELATIVE_ORIENTATION_OPS" , "INTERNAL_OPS", "LOCAL_DENSITIES" , "TOTAL_OPS"] 
assert len(sys.argv) > 1,  "Need to supply one COP file, and how many to split it into!"
copFile = CopFile.open(sys.argv[1])
cops: List[CopFile] = []
numCells : List[int] = []
for l in sys.argv[2:len(sys.argv)]:
    try:
        numCells.append(int(l))
    except ValueError:
        print("expecting an int for lengths!")
totalCell = copFile.ncells
index = 0
for l in numCells:
    cops.append(CopFile())
    cops[-1].tags = copFile.tags
    cops[-1].switch_width = copFile.switch_width
    cops[-1].data["GRID"] = copFile.data["GRID"]
    cops[-1].data["GRID"][0] = l #l as in lima
    cops[-1].data["GRID"][1] = 1
    cops[-1].data["GRID"][2] = 1
    cops[-1].cutoff = copFile.cutoff
    cops[-1].ncells = int(l) #should probably have an input file for the other things
    for op in ops:
        cops[-1].data[op] = copFile.data[op][index:l+index]
    index = l 
ind = 0
for i in cops:
    with open(f"cop{ind}", "w") as temp:
        temp.write(i.to_string())
        temp.close()
    ind+=1
