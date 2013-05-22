import sys

from ase import atoms

amSymToInt = {
  "s" : 0,
  "p" : 1,
  "d" : 2,
  "f" : 3,
  "g" : 4
}

aufbau = {
   1: "s",                   
   2: "s",  3: "p",          
   4: "s",  5: "p",  7: "d", 
   6: "s",  8: "p", 10: "d", 13: "f",
   9: "s", 11: "p", 14: "d", 17: "f", 21: "g",
  12: "s", 15: "p", 18: "d", 22: "f",
  16: "s", 19: "p", 23: "d",
  20: "s", 24: "p",
  25: "s"
}

def getMaxNumberEl(sym):
  """ returns the maximal allowed number of electrons
      in a "sym" shell
  """
  am = amSymToInt[sym]
  return (2*am + 1) * 2


def readTMToDict(pathToFile):
  """ Reads all basissets in the given file 
      to a dict and returns it.
  """
  basisset = {}
  
  print "Convert all basissets from", pathToFile

  inFileHandle = open(pathToFile, "r")
  for line in inFileHandle:
  
    items = line.split()

    # skip empty lines
    if items == []: continue
    
    # search atom names otherwise skip that
    if items[0].title() not in atoms.chemical_symbols: continue

    symbol = items[0].title()
    basisset[symbol] = []
    
      
    for item in items:
      if "ecp" in item.strip(): 
        print "EEEEECCCCCCCCCPPPPPPPPPP"


    print "reading basisset for", symbol   
    # there should be a '*' next.. skip comments
    for line in inFileHandle:
      # skip comments
      if line.strip()[0] == "#": continue
      # found the right smybol
      if line.strip()[0] == "*": break

    # write the actual basis function
    for line in inFileHandle:
      items = line.split()

      if items == []: 
        print "WARNING: Empty line in atomic block!"
        continue
      if items[0].strip()[0] == "*": 
        print "finished basisset for ", symbol
        break
    
      functype = items[1].strip()
      function = { "type"    : functype, 
                   "numbers" : []      ,
                   "numEl"   : 0       ,}

      numLinesToRead = int(items[0])
      for line in inFileHandle:
        items = line.split()
        function["numbers"].append( (float(items[0]), float(items[1])) )
        numLinesToRead -= 1
        if numLinesToRead == 0: break
      
      basisset[symbol].append(function)

  inFileHandle.close()
  return basisset

def writeAsCrystal(basisset):
  """ Writes the data given in the basiset
      dict to crystal files.
  """

  for key in basisset.keys():
    filename = "%02d_%s.crystal" % (atoms.atomic_numbers[key], key)
    print "Writing ", filename

    fileH = open(filename, "w")
    # writing the header to the file
    funcs = basisset[key]
    fileH.write("%d %d\n" % (atoms.atomic_numbers[key], len(funcs)))

    # distribute the electrons with the aufbau principle
    numUndistributedEl = atoms.atomic_numbers[key]
    counter = 1
    while (numUndistributedEl > 0):
      shellSym     = aufbau[counter]
      maxElInShell = getMaxNumberEl(shellSym) 
      # how many electrons will be distributed in this step
      if (maxElInShell < numUndistributedEl):
        numEl = maxElInShell
      else : 
        numEl = numUndistributedEl
      numUndistributedEl -= numEl
      # search for the apropriate shell
      for shell in funcs:
        if shell["type"] == shellSym.strip() and shell["numEl"] == 0:
          shell["numEl"] = numEl
          break
      counter += 1


    for func in funcs:
      fileH.write("0 %d %d %.1f 1.0\n" % 
                   (amSymToInt[func["type"].lower()],
                    len(func["numbers"]),
                    func["numEl"]
                   )
                 )
      for prim in func["numbers"]:
        fileH.write("  %f   %f\n" % prim)

    fileH.close()

if __name__ == "__main__":
  basisset = readTMToDict(sys.argv[1])
  writeAsCrystal(basisset)

