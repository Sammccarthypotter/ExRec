import os
import sys
from pathlib import Path

#help function 
def HelpCommandLine():
  print("Version 1.0 October 6, 2023")
  print("To run PhyCombine.py, please perform the following steps:\n")
  print("1) place your input PHYLIP files into a folder named \"phylip_files.\"")
  print("2) place \"phylip_files\" into the folder containing the Phycombine.py application.")
  print("3) if your input files are in strict sequential or strict interleaved formats, then type the following commands on the command line followed by enter:\n")
  print("python3 Phycombine.py \n")
  print("If your input files are in relaxed sequential or relaxed interleaved formats, then type the following commands on the command line followed by enter:\n")
  print("python3 Phycombine.py r \n")
  print("Note: because of the various PHYLIP file formats in use, we highly recommend that you read \nPages 5-9 in the Manual.pdf file, which illustrates many of the different \nPHYLIP format types and the corresponding command line syntax you need to use to successfully \nrun Phycombine.py. \n")
  exit()

def find_number(lis):
  str_list = list(lis)
  output = []
  for letter in str_list:
    if letter.isdigit():
      output.append(letter)
  if len(output) > 0:
    return int(''.join(output))
  else: 
    return len(lis)
 
#checks if file is openable
def inputPath(inputFilePath):
   try:
       return open(inputFilePath)
   except Exception as e:
       print('Cannot open', inputFilePath, file=sys.stderr)
       exit()
 
#prints the final file in the nexus format 
def printOutFile(finalFile,outputFilePath):
   ncha = 0
   with open(outputFilePath, 'w') as f:
       for locus in range(len(finalFile.charset)):
           ncha += len(finalFile.Fullfilecode2[0][locus])
       sys.stdout = f  # Change the standard output to the file we created.
       print("#NEXUS")
       print("begin data;")
       print("\tdimensions ntax=" + str(finalFile.numberOfspecies) + " nchar=" +
             str(ncha) + ";")
       print("\tformat datatype=dna missing=? gap=- interleave;")
       print("matrix")
       for locus in range(len(finalFile.charset)):
           for i in range(0, len(finalFile.Fullfilecode2[0][locus]), 70):
               for species in range(finalFile.numberOfspecies):
                   row = finalFile.Fullfilecode2[species][locus][i:i + 70].copy()
                   for x in range(len(row)):
                     if row[x] == 'N':
                       row[x] = '?'
                            
                   name = list(finalFile.speciesName[species][0:finalFile.correction])
                   for ele in range(correction - len(name)):
                     name.append(" ")
                   print(''.join(name),''.join(row))
               print()
       print(";")
       print("end;")
       if len(finalFile.charset) > 1:
           print("#nexus")
           print("begin sets;")
           outputlength = 0
           backSlash = "\\"
           for locus in range(len(finalFile.charset)):
               print("charset " + "'" + fullLociNames[locus].split(backSlash)[-1].split('.')[0] + ".phy" + "'" + " = " + str(outputlength+1) + "-" + str(len(finalFile.Fullfilecode2[0][locus])+outputlength) + ";")
               outputlength += len(finalFile.Fullfilecode2[0][locus])
           print("end;")
 
class FinalFile:
 def __init__ (self,charset,Fullfilecode2,numberOfspecies,speciesName,correction):
       self.charset = charset
       self.Fullfilecode2 = Fullfilecode2
       self.numberOfspecies = numberOfspecies
       self.speciesName = speciesName
       self.correction = correction

class PhlipFile:
 def IdentifyFormat(self):
  return "Relaxed "
 def Speciesfinder(self,speciesList):
    global fullspecieslist
    try:
      fullspecieslist
    except:
      fullspecieslist =[]
    for species in speciesList:
      if species not in fullspecieslist:
        fullspecieslist.append(species)
 def nucleotideaccepted(self,Scan):
   if (Scan == "-" or Scan == "C" or Scan == "G" or Scan == "T" or Scan == "A" or Scan == "?" or Scan == "N"):
     return True
   else:
     return False

 def Search_Loci(self,TestInputFile):
  charset = [[0,0]]
  speciesNumber=[]
  for line in TestInputFile:
    Scan = line.split()
    if len(Scan) == 2:
      charset.append([charset[-1][-1]+1,charset[-1][1]+int(Scan[1])])
      speciesNumber.append(Scan[0])
      break
  
  if len(charset) > 1:
    charset.pop(0)
  return [speciesNumber,charset]

 def readFileListUpper(self,fileList):
   #Reads relaxed sequential with name on upper line 
   speciesNameLine = True
   code=0
   numChar = 0
   charsetRead = False
   for line in (fileList):
    if charsetRead:
      if speciesNameLine and len(line) > 0:
        self.speciesName[code]=(list(line.strip(" ")))
        speciesNameLine = False
      else:
        if len(line) == BaseNumChar:
          self.Fullfilecode[code].extend(line)
          speciesNameLine = True
          code += 1
    else:
      lineList = line.split()
      if len(lineList) == 2:
        BaseNumChar=int(lineList[-1])
        charsetRead = True

 def readFileListStrict(self,fileList):
   #Reads strict sequential and strict interleaced with name on same line  
   seenCharsetLine = False
   code=0
   speciesNameCorrection = 10
   for line in (fileList):
    fileTextStream = list(line)
    if seenCharsetLine == True and len(fileTextStream) > 0:
      for block in line[speciesNameCorrection:].split():
        self.Fullfilecode[code].extend(list(block))
      if  speciesNameCorrection > 1:
        self.speciesName[code]=((fileTextStream[0:10]))
      code +=1
    if seenCharsetLine == False and len(line.split()) == 2:
      seenCharsetLine = True
    if code >= self.numberOfspecies:
      code =0
      speciesNameCorrection = 0

 def readFileListWraped(self,fileList):
   #Reads wraped with name on upper line 
   speciesNameLine = True
   code=0
   numChar = 0
   for line in (fileList):
    if speciesNameLine:
      lineList = line.split()
      if len(lineList)<2:
        self.speciesName[code]=(list(line.strip(" ")))
        speciesNameLine = False
      else:
        BaseNumChar=int(lineList[-1])
    else:
      numChar += len(line)
      self.Fullfilecode[code].extend(line)
      if numChar == BaseNumChar:
        numChar = 0
        speciesNameLine = True
        code += 1

 def readFileListRelaxed(self,fileList):
  #Reads relaxed interleaved files and relaxed sequential 
   seenCharsetLine = False
   code=0
   speciesListNotFull = True
   for line in (fileList):
    fileTextStream = line.split()
    if seenCharsetLine == True and len(fileTextStream) > 0:
      for block in fileTextStream[speciesListNotFull:]:
        self.Fullfilecode[code].extend(list(block))
      if speciesListNotFull:
        self.speciesName[code]=((fileTextStream[0]))
      code +=1
    if seenCharsetLine == False and len(fileTextStream) == 2:
      seenCharsetLine = True
    if code >= self.numberOfspecies:
      code =0
      speciesListNotFull = False

 def readFileListRelaxedWraped(self,fileList):
   seenCharsetLine = False
   code=0
   speciesNameLine = True
   numChar = 0
   for line in (fileList):
    fileTextStream = line.split()
    if seenCharsetLine == True and len(line) > 0:
      if speciesNameLine:
        self.speciesName[code]=fileTextStream[0]
        self.Fullfilecode[code].extend(list(fileTextStream[-1]))
        numChar += len(list(fileTextStream[-1]))
        speciesNameLine = False
      else:
        numChar += len(line)
        self.Fullfilecode[code].extend(list(line))
        if numChar == BaseNumChar:
          numChar = 0
          speciesNameLine = True
          code += 1
    if seenCharsetLine == False and len(fileTextStream ) == 2:
      seenCharsetLine = True
      lineList = line.split()
      BaseNumChar=int(lineList[-1])

 def readFileListStrictWraped(self,fileList):
   seenCharsetLine = False
   code=0
   speciesNameCorrection = 10
   speciesNameLine = True
   numChar = 0
   for line in (fileList):
    if seenCharsetLine == True and len(line) > 0:
      if speciesNameLine:
        self.speciesName[code]=line[0:speciesNameCorrection]
        self.Fullfilecode[code].extend(list(line[speciesNameCorrection:]))
        numChar += len(list(line[speciesNameCorrection:]))
        speciesNameLine = False
      else:
        numChar += len(line)
        self.Fullfilecode[code].extend(list(line))
        if numChar == BaseNumChar:
          numChar = 0
          speciesNameLine = True
          code += 1
    if seenCharsetLine == False and len(line.split()) == 2:
      seenCharsetLine = True
      lineList = line.split()
      BaseNumChar=int(lineList[-1])

 def __init__ (self, inputFilePath):
 
   fileList = [line.rstrip('\n') for line in inputPath(inputFilePath)]
 
   TestInputFile = inputPath(inputFilePath)
 
   TestInputFile2 = inputPath(inputFilePath)
 
   numberOfLines = 0
 
   self.correction = 0
   self.correction2 = 0
   self.line_To_Start = 0
 
   #checks if the file is a Philip file
   FirstLociData = self.Search_Loci(TestInputFile)
   self.numberOfspecies = int(FirstLociData[0][0])
   self.nchar=FirstLociData[-1][-1][-1]

   numberOfLines = self.line_To_Start
   self.Fullfilecode = []
   self.speciesName = []
   self.Fullfilecode2 = [] 
   #initializing lists
   for spec3 in range(self.numberOfspecies):
       self.Fullfilecode.append([])
       self.speciesName.append([])
       self.Fullfilecode2.append([])

   format = "Strict"
   
   flags = []

   for string in sys.argv:
    flags.append(string.upper())

   if len(sys.argv) >= 2:
    if "R" in flags:
      format = "Relaxed"
    if "U" in flags:
        format = "Upper"
    if "W" in flags:
      format = "Strict-Wraped"
      if "R" in flags:
        format = "Relaxed-Wraped"
      if "U" in flags:
        format = "Wraped-Upper"

    print(flags)
   
   print(format)
   if  format == "Upper":
    self.readFileListUpper(fileList)
   elif format == "Wraped-Upper":
    self.readFileListWraped(fileList)
   elif format == "Relaxed":
    self.readFileListRelaxed(fileList)
   elif format == "Strict":
    self.readFileListStrict(fileList)
   elif format == "Strict-Wraped":
    self.readFileListStrictWraped(fileList)
   elif format == "Relaxed-Wraped":
    self.readFileListRelaxedWraped(fileList)

   Locis = [[]]
   locusNumber = 0
   LociName = []
   charset=FirstLociData[1]
 
   #BREAKS stream of nucleotide into separate Loci
   for locus in charset:
       if locusNumber == 0:
        LociName.append(str(inputFilePath).split('/')[-1])
       else:
        LociName.append(str(inputFilePath).split('/')[-1]+"_"+str(locusNumber))
       for species in range(self.numberOfspecies):
           Locis[locusNumber].append(self.Fullfilecode[species][int(locus[0]) - 1:int(locus[-1])])
       locusNumber += 1
       Locis.append([])

   for locus in range(len(charset)):
       for spec6 in range(self.numberOfspecies):
           self.Fullfilecode2[spec6].append(Locis[locus][spec6])

   self.charset = charset
   self.LociName = LociName 
   self.Speciesfinder(self.speciesName)

if __name__ == "__main__":
    #comand line arguments
  if len(sys.argv) == 2:
    helpCheck = sys.argv[1]
    if helpCheck.lower() == "help":
       HelpCommandLine()

  original_stdout = sys.stdout  # Save a reference to the original standard output
 
  this = Path(__file__).absolute().parent

  inputFolder = 'phylip_files'

  if not os.path.isdir(this / inputFolder):
    os.makedirs(inputFolder)
  
  FileList = os.listdir(this / inputFolder)

  #this removes hidden files 
  for file in FileList:
      if file[:1] == ".":
          FileList.remove(file)

  FileList.sort()

  FileList.sort(key=find_number)
  
  print(FileList)

  inputFilePath = this / inputFolder
  filelist=[]
  for file in FileList:
    filelist.append(PhlipFile(inputFilePath / file))
  try:
    Nexus1 = filelist[0]
  except: 
    raise IndexError("No phylip_files")
  
  fullCharset =[]
  fullLociNames = []
 
  for file in filelist:
    fullLociNames = fullLociNames + file.LociName
    if len(fullCharset) == 0:
      fullCharset = fullCharset + file.charset
    else:
      newSet = []
      for loci in file.charset:
          newLoci = [int(loci[0])+int(fullCharset[-1][-1]),int(loci[-1]) + int(fullCharset[-1][-1])]
          newSet = newSet + [newLoci]
      fullCharset = fullCharset + newSet

  Fullfilecode= []
  for file in filelist:
    if len(Fullfilecode) == 0:
      Fullfilecode = []
      for species in range(len(fullspecieslist)):
        Fullfilecode.append([])
    speciesIndex = 0
    for species in fullspecieslist:
      if species in file.speciesName:
        finalIndex = file.speciesName.index(species)
        NewData = file.Fullfilecode2[finalIndex]
      else:
        questionlist = []
        for l in range(len((file.Fullfilecode2[0][0]))):
          questionlist.append("?")
        NewData = [questionlist]
      Fullfilecode[speciesIndex].extend(NewData)
      speciesIndex += 1


  numberOfspecies = len(fullspecieslist)
  speciesName = fullspecieslist
  max_len = -1
  correction = 0
  for element in fullspecieslist:
      if len(element) > max_len:
          max_len = len(element)
          correction = len(element)
  
  Fullfilecode2 = Fullfilecode
  charset = fullCharset
  
  finalFile = FinalFile(charset,Fullfilecode2,numberOfspecies,fullspecieslist,correction)
  
  outputFilePath = "combineloci.nex"
  with open(outputFilePath, 'w') as f:
    sys.stdout = f
    print()
  sys.stdout = original_stdout
  
  printOutFile(finalFile,outputFilePath)
  sys.stdout = original_stdout
