import os
import sys
from pathlib import Path

def HelpCommandLine():
  print("Version 1.0 October 6, 2023")
  print("To run Physplit.py, please perform the following steps:")
  print("1) place the \"Trunc_combineloci.nex\" file (which was output by FGT.py) into the folder containing PhySplit.py.")
  print("2) type the following commands on the command line followed by enter: \n")
  print("python3 PhySplit.py")
  exit()

original_stdout = sys.stdout  # Save a reference to the original standard output
 
this = Path(__file__).absolute().parent

inputFolder = 'nexus_files'

FileList = ["Trunc_combineloci.nex"]

#this removes hidden files 
for file in FileList:
    if file[:1] == ".":
        FileList.remove(file)
 
#print(FileList)
 
#checks if file is openable
def inputPath(inputFilePath):
   try:
       return open(inputFilePath)
   except Exception as e:
       print('Cannot open', inputFilePath, file=sys.stderr)
       exit()
 
#this looks for the first valid line of the DNA code
def printOutFile(finalFile,outputFilePath,locus):
  numberOfspecies = finalFile.numberOfspecies
  Fullfilecode = finalFile.Fullfilecode2.copy()
  speciesName = finalFile.speciesName.copy()
  qustionlist = []
  for l in range(len((Fullfilecode[0][locus]))):
    qustionlist.append("?")
  for species in range(finalFile.numberOfspecies-1,0,-1):
    if qustionlist == Fullfilecode[species][locus][0:]:
      Fullfilecode.pop(species)
      speciesName.pop(species)
      numberOfspecies = numberOfspecies -1
  count = 0
  with open(outputFilePath, 'w') as f:
    sys.stdout = f
    print()
    print(" ",numberOfspecies ,"     ",len(Fullfilecode[0][locus]))
    print()
    for species in range(numberOfspecies):
        count += 1
        row = Fullfilecode[species][locus][0:]
        name = ''.join(finalFile.speciesName[species][0:finalFile.correction]).strip()
        name += "-^"
        name += fullLociNames[locus]
        print(''.join(name))
        print(''.join(row))
    print()


 
class FinalFile:
 def __init__ (self,charset,Fullfilecode2,numberOfspecies,speciesName,correction):
       self.charset = charset
       self.Fullfilecode2 = Fullfilecode2
       self.numberOfspecies = numberOfspecies
       self.speciesName = speciesName
       self.correction = correction
 
 
class NexusFile:
 def Speciesfinder(self,inputFilePath):
    global fullspecieslist
    try:
      fullspecieslist
    except:
      fullspecieslist =[]

    correction = self.correction
    fileList = [line.strip('\n') for line in inputPath(inputFilePath)]
    lengthlist = fileList.index('end;')
    for line in range(self.line_To_Start,lengthlist):
      Nspecies = list((fileList[line][0:correction]))
      if len(Nspecies) == (correction):
        Nspecies = list((fileList[line][0:correction]).strip())
        Nspeciesup = list((fileList[line-1][0:correction]).strip())
        if line+1 < lengthlist:
         Nspeciesdown = list((fileList[line+1][0:correction]).strip())
        if Nspecies not in fullspecieslist:
          if Nspeciesup in fullspecieslist:
           fullspecieslist.insert(fullspecieslist.index(Nspeciesup)+1,Nspecies)
          elif Nspeciesdown in fullspecieslist:
           fullspecieslist.insert(fullspecieslist.index(Nspeciesdown),Nspecies)
          else:
           fullspecieslist.insert(0,Nspecies)
    
 def nucleotideaccepted(self,Scan):
   if (Scan == "-" or Scan == "C" or Scan == "G" or Scan == "T" or Scan == "A" or Scan == "?" or Scan == "N"):
     return True
   else:
     return False

 def Search_Starting_Line(self,TestInputFile):
   self.line_To_Start =0
   self.correction
   for line2 in TestInputFile:
       Scan = list(line2)
       for letter in range(len(Scan)):
           if(len(Scan) - letter > 5) and (self.nucleotideaccepted(Scan[letter]) == True):
             if (self.nucleotideaccepted(Scan[letter + 1])) and (self.nucleotideaccepted(Scan[letter + 2])) and (self.nucleotideaccepted(Scan[letter + 4]) == True):
                   self.correction = letter
                   return self.correction
       self.line_To_Start += 1
      
 def __init__ (self, inputFilePath):
 
   fileList = [line.rstrip('\n') for line in inputPath(inputFilePath)]
 
   TestInputFile = inputPath(inputFilePath)
 
   TestInputFile2 = inputPath(inputFilePath)

   numberOfLines = 0
 
   self.correction = 0
   self.correction2 = 0
   self.line_To_Start = 0
 
   charset = []
   LociName = []
 
   #checks if the file is a Nexus file
   if fileList[0] == "#NEXUS":
       line = 0
       while (fileList[line]).upper() != "BEGIN DATA;" and line < 10:
           line += 1
       comandLine = fileList[line + 1].split()
       comandLineSplit = comandLine[1].split("=")
       numberOfspecies = int(comandLineSplit[1])
       comandLineSplit = comandLine[2].split("=")
       nchar = int(comandLineSplit[1].strip(';'))
   if fileList[0] != "#NEXUS":
       self.Search_Starting_Line()
   else:
       while (fileList[line]).upper().strip() != "MATRIX":
           line += 1
       self.line_To_Start = line + 1
   if fileList[0] == "#NEXUS":
       line = -1
       while (fileList[line]) != ";":
           line -= 1
           if fileList[line].upper().strip() == "BEGIN SETS;":
               for x in range(len(fileList) + line, len(fileList) - 1):
                   lineToCheck = fileList[x].strip(";").split()
                   if len(lineToCheck) > 0 and lineToCheck[0] == "charset":
                       LociName.append(lineToCheck[1].strip("'"))
                       if len(lineToCheck) == 4:
                         charset.append(lineToCheck[-1].split("-"))
                       elif len(lineToCheck) == 6:
                         charset.append([lineToCheck[-3],lineToCheck[-1]])
                  
       if len(charset) < 1:
           charset.append([1, nchar])
   if len(LociName) == 0:
    LociName = [inputFilePath]
   for i in range(len(LociName)):
    if not LociName[i].split(".")[-1] == ".phy":
      LociName[i] = LociName[i].split(".")[0] + ".phy"



   numberOfLines = self.line_To_Start
   self.Search_Starting_Line(TestInputFile)
   SortingList = []
   Fullfilecode = []
   spicesName = []
   Fullfilecode2 = [] 
   #initializing lists
   for spec3 in range(numberOfspecies):
       Fullfilecode.append([])
       spicesName.append([])
       Fullfilecode2.append([])
   firstemptyline = True
   #counts the number of species
   if fileList[0] != "#NEXUS":
       for line in TestInputFile2:
           numberOfLines += 1
           if numberOfLines >= 0 and len(line) <= 3 and firstemptyline:
               numberOfspecies = numberOfLines
               firstemptyline = False
               break
 
   firstemptyline = True
 
     #Reads the input file breaks it down in to list of nucleotides
   while fileList[numberOfLines] != ";":
       line = fileList[numberOfLines]
       numberOfLines += 1
       if numberOfLines >= 0 and len(line) > 4:
           SortingList.append(line)
           if len(SortingList) > numberOfspecies:
               SortingList.clear()
               SortingList.append(line)
           if len(SortingList) == numberOfspecies:
               for code in range(numberOfspecies):
                   dnaOnly = list(
                       SortingList[code])[self.correction:len(list(SortingList[code]))]
                   Fullfilecode[code].extend(dnaOnly)
               if firstemptyline:
                   for code5 in range(numberOfspecies):
                       dnaOnly = list(SortingList[code5])[0:self.correction]
                       spicesName[code5].extend(dnaOnly)
                   firstemptyline = False
 
   Locis = [[]]
   locusNumber = 0

   #BREAKS stream of nucleotide into separate Loci
   for locus in charset:
       for species in range(numberOfspecies):
           Locis[locusNumber].append(Fullfilecode[species][int(locus[0]) - 1:int(locus[-1])])
       locusNumber += 1
       Locis.append([])

 
   for locus in range(len(charset)):
       for spec6 in range(numberOfspecies):
           Fullfilecode2[spec6].append(Locis[locus][spec6])
 
   self.Fullfilecode= Fullfilecode
   self.charset = charset
   self.nchar = nchar
   self.Fullfilecode2 = Fullfilecode2
   self.numberOfspecies = numberOfspecies
   self.spicesName = spicesName
   self.LociName = LociName 

   self.Speciesfinder(inputFilePath)
 
if __name__ == "__main__":
  #comand line arguments
  if len(sys.argv) == 2:
    helpCheck = sys.argv[1]
    if helpCheck.lower() == "help":
       HelpCommandLine()


inputFilePath = "Trunc_combineloci.nex"

if len(sys.argv) == 2:
  inputFilePath = sys.argv[1]

filelist=[]
for file in FileList:
 filelist.append(NexusFile(inputFilePath))
try:
  Nexus1 = filelist[0]
except: 
  raise IndexError("No nexus_files in folder nexus_files")
  
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
 back = 0
 for species in range(len(fullspecieslist)):
      if species-back >= len(file.spicesName):
        back += 1
      if (list((''.join(file.spicesName[species-back])).strip()) == fullspecieslist[species]):
        Fullfilecode[species] = Fullfilecode[species]+(file.Fullfilecode2[species-back])
      else:
        back += 1
        qustionlist = []
        for l in range(len((file.Fullfilecode2[species-back][0]))):
          qustionlist.append("?")
        Fullfilecode[species] = Fullfilecode[species] + [qustionlist]


numberOfspecies = len(fullspecieslist)
spicesName = fullspecieslist

max_len = -1
for ele in fullspecieslist:
    if len(ele) > max_len:
        max_len = len(ele)
        correction = len(ele)
 
Fullfilecode2 = Fullfilecode
charset = fullCharset
 
finalFile = FinalFile(charset,Fullfilecode2,numberOfspecies,fullspecieslist,correction)
 
outputFolder = "phylip_split_files"
 
if not os.path.isdir(this / outputFolder):
   os.makedirs(outputFolder)

for loci in range(len(finalFile.charset)):
  outputFilePath = this / outputFolder / fullLociNames[loci]
  print(len(finalFile.Fullfilecode2))
  with open(outputFilePath, 'w') as f:
    sys.stdout = f
    printOutFile(finalFile,outputFilePath,loci)
  sys.stdout = original_stdout
