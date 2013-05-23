import sys

from ase import atoms

amSymToInt = {
  "s" : 0,
  "p" : 1,
  "d" : 2,
  "f" : 3,
  "g" : 4
}

amSymToIntCrystal = {
  "s" : 0,
  "sp": 1,
  "p" : 2,
  "d" : 3,
  "f" : 4,
  "g" : 5
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

amOrder = ["s", "p", "d", "f", "g"]

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

  def readBasis(inFileHandle, basisset):
    """ to read a basis
    """
    # read the actual basis function
    for line in inFileHandle:
      items = line.strip().split()

      if items[0][0] == "#": continue # skip comments
      if items == []: 
        print "WARNING: Empty line in atomic block!"
        continue
      if items[0].strip()[0] == "*": 
        print "finished basisset for ", symbol
        break
    
      functype = items[1].strip()
      function = { "type"    : functype, 
                   "exp"     : []      ,
                   "coeff"   : []      ,
                   "numEl"   : 0       ,}

      numLinesToRead = int(items[0])
      for line in inFileHandle:
        items = line.split()
        function["exp"].append( float(items[0]) )
        function["coeff"].append( float(items[1]) )
        numLinesToRead -= 1
        if numLinesToRead == 0: break
      
      basisset[symbol]["basis"].append(function)

  def readECP(inFileHandle, basisset, symbol):
    for line in inFileHandle:
      items = line.strip()
      if items[0] == "#": continue # skip comments
      
      # if we found this we can break
      if items[0:5] == "ncore":        
        items = items.replace("=", " ").split()
        basisset[symbol]["ecp_ncore"] = items[1]
        basisset[symbol]["ecp_lmax"]  = items[3]
        break 

      items = line.split()
      if items == []: 
        print "WARNING: Empty line in ecp atomic block!"
        continue


    # found the shell part.. read that mofo    
    ecp = None
    for line in inFileHandle:
      items = line.strip().split()
      
      if items == []: 
        print "WARNING: Empty line in ecp atomic block!"
        continue
      if items[0].strip()[0] == "*": 
        print "finished ecp for ", symbol
        break
      if items[0] == "#": continue  # skip comments
      if items[0][0] == "*" : break # this is the end

      # found new shell thingy .. write and  build new ecp 
      if len(items) == 1:
        if ecp is not None : basisset[symbol]["ecp"].append(ecp)
        functype = items[0].strip()
        ecp = { "type"    : functype, 
                "exp"     : []      ,
                "coeff"   : []      ,
                "r^n"     : []      }
        continue

      ecp["exp"].append(items[0])
      ecp["coeff"].append(items[2])
      ecp["r^n"].append(items[1])

    # the forgotten one
    if ecp is not None: basisset[symbol]["ecp"].append(ecp)
    
    print symbol, basisset[symbol]["ecp"]

  #############################################################################
  ##                                                                         ##
  ##                            actual function                              ##
  ##                                                                         ##
  #############################################################################

  basisset = {}
  
  print "Read all basissets from", pathToFile

  inFileHandle = open(pathToFile, "r")

  for line in inFileHandle:
  
    items = line.split()

    # skip empty lines
    if items == []: continue
    
    # search atom names otherwise skip that
    if items[0].title() not in atoms.chemical_symbols: continue

    symbol = items[0].title()
    
    # if not existent creat the key
    if symbol not in basisset.keys():
      basisset[symbol] = {"basis"     : [],
                          "ecp"       : [],
                          "ecp_ncore" : 0 ,
                          "ecp_lmax"  : 0 }
    
    # check if we got an ecp case
    isECP = False 
    for item in items:
      if "ecp" in item.strip(): 
        isECP = True


    print "reading for", symbol   
    # there should be a '*' next.. skip comments
    for line in inFileHandle:
      # skip comments
      if line.strip() == "" : continue
      if line.strip()[0] == "#": continue
      # found the right smybol
      if line.strip()[0] == "*": break
    
    if not isECP:
      readBasis(inFileHandle, basisset)
    else:
      readECP(inFileHandle, basisset, symbol)

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
    funcs = basisset[key]["basis"]
    ecps  = basisset[key]["ecp"]
    number = atoms.atomic_numbers[key]
    if ecps != []: number += 200
    fileH.write("%d %d\n" % (number, len(funcs)))

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
    
    if ecps != []:
      fileH.write("INPUT\n")

      # get ecp information 
      nECP = [0]*6
      for ecp in ecps:
        if ecp["type"][0] == "s": nECP[1] += len(ecp["exp"])
        if ecp["type"][0] == "p": nECP[2] += len(ecp["exp"])
        if ecp["type"][0] == "d": nECP[3] += len(ecp["exp"])
        if ecp["type"][0] == "f": nECP[4] += len(ecp["exp"])
        if ecp["type"][0] == "g": nECP[5] += len(ecp["exp"])

      fileH.write("%d. %d %d %d %d %d %d\n" % 
                  (number-200-int(basisset[key]["ecp_ncore"]),
                  nECP[0],nECP[1],nECP[2],nECP[3],nECP[4],nECP[5]
                  )
                 )

      for sym in amOrder:
        for ecp in ecps:
          if ecp["type"][0] == sym:
            for coeff, exp, rn in zip(ecp["coeff"],ecp["exp"],ecp["r^n"]):
              fileH.write("  %4.6f   % 4.6f   %d\n" % 
                          (float(coeff),float(exp),int(rn))
                         )


    for func in funcs:
      fileH.write("0 %d %d %.1f 1.0\n" % 
                   (amSymToIntCrystal[func["type"].lower()],
                    len(func["exp"]),
                    func["numEl"]
                   )
                 )
      for exp, coeff in zip(func["exp"],func["coeff"]):
        fileH.write("  %f   %f\n" % (exp, coeff))

    fileH.close()

if __name__ == "__main__":
  basisset = readTMToDict(sys.argv[1])
  writeAsCrystal(basisset)

