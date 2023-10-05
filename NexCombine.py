import os
import sys
from pathlib import Path

def HelpCommandLine():
  print()
  print("To run NexCombine.py, please perform the following steps:\n")
  print("1) place your input NEXUS files into a folder named \"nexus_files.\"")
  print("2) place \"nexus_files\" into the folder containing the NexCombine.py application.")
  print("3) type the following commands on the command line followed by enter:\n")
  print(">python3 NexCombine.py")
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

#this looks for the first valid line of the DNA code
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
                            
                   name = finalFile.speciesName[species][0:finalFile.correction]
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
               print("charset " + "'" + fullLociNames[locus].split(backSlash)[-1].split('.')[0]+ "'" + " = " + str(outputlength+1) + "-" + str(len(finalFile.Fullfilecode2[0][locus])+outputlength) + ";")
               outputlength += len(finalFile.Fullfilecode2[0][locus])
           print("end;")
 
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
    fileList = [line.strip('\t').strip('\n').strip() for line in inputPath(inputFilePath)]
    if 'end;' in fileList:
      lengthlist = fileList.index('end;')
    elif 'End;' in fileList:
      lengthlist = fileList.index('End;')
    elif 'END;' in fileList:
      lengthlist = fileList.index('END;')
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
 
   inputFile = inputPath(inputFilePath)
 
   TestInputFile = inputPath(inputFilePath)
 
   TestInputFile2 = inputPath(inputFilePath)
 
   TestInputFile3 = inputPath(inputFilePath)
 
   numberOfLines = 0
 
   self.correction = 0
   self.correction2 = 0
   self.line_To_Start = 0
 
   charset = []
 
   LocusNumber = []
   LocusName  = []
   startinglength = []
   S = []
   infiniteSites = []
   Rm = []
   locationsRM = []
   longestBlock = []
   finalLength = []
 
   #checks if the file is a Nexus file
   if fileList[0] == "#NEXUS":
       line = 0
       while (fileList[line].strip().strip('\t')).upper() != "BEGIN DATA;" and line < 10:
           line += 1
       commandLine = fileList[line + 1].split()
       commandLineSplit = commandLine[1].split("=")
       numberOfspecies = int(commandLineSplit[1])
       commandLineSplit = commandLine[2].split("=")
       nchar = int(commandLineSplit[1].strip(';'))
   if fileList[0] != "#NEXUS":
       self.Search_Starting_Line()
   else:
       while (fileList[line]).upper().strip() != "MATRIX":
           line += 1
       self.line_To_Start = line + 1
   if fileList[0] == "#NEXUS":
       line = -1
       while (fileList[line].strip().strip('\t')) != ";":
           line -= 1
           if fileList[line].upper().strip().strip('\t') == "BEGIN SETS;":
               for x in range(len(fileList) + line, len(fileList) - 1):
                   lineToCheck = fileList[x].strip(";").split()
                   if len(lineToCheck) > 0 and lineToCheck[0] == "charset":
                       if len(lineToCheck) == 4:
                         charset.append(lineToCheck[-1].split("-"))
                       elif len(lineToCheck) == 6:
                         charset.append([lineToCheck[-3],lineToCheck[-1]])
                  
       if len(charset) < 1:
           charset.append([1, nchar])
 
   print(charset)
   numberOfLines = self.line_To_Start
   self.Search_Starting_Line(TestInputFile)
   SortingList = []
   Fullfilecode = []
   speciesName = []
   Fullfilecode2 = [] 
   #initializing lists
   for spec3 in range(numberOfspecies):
       Fullfilecode.append([])
       speciesName.append([])
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
   while fileList[numberOfLines].strip().strip('\t') != ";":
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
                       speciesName[code5].extend(dnaOnly)
                   firstemptyline = False
 
   Locis = [[]]
   locusNumber = 0
   LociName = []
   
   #BREAKS stream of nucleotide into separate Loci
   for locus in charset:
       if locusNumber == 0:
        LociName.append(str(inputFilePath).split('/')[-1])
       else:
        LociName.append(str(inputFilePath).split('/')[-1]+"_"+str(locusNumber))
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
   self.speciesName = speciesName
   self.LociName = LociName 

   self.Speciesfinder(inputFilePath)

if __name__ == "__main__":
  #comand line arguments
  if len(sys.argv) == 2:
    helpCheck = sys.argv[1]
    if helpCheck.lower() == "help":
       HelpCommandLine()

  original_stdout = sys.stdout  # Save a reference to the original standard output
 
  this = Path(__file__).absolute().parent

  inputFolder = 'nexus_files'

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
    filelist.append(NexusFile(inputFilePath / file))
  try:
    Nexus1 = filelist[0]
  except: 
    raise IndexError("No nexus_files")
 
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
          if species-back >= len(file.speciesName):
            back += 1
          if (list((''.join(file.speciesName[species-back])).strip()) == fullspecieslist[species]):
            Fullfilecode[species] = Fullfilecode[species]+(file.Fullfilecode2[species-back])
          else:
            back += 1
            questionlist = []
            for l in range(len((file.Fullfilecode2[species-back][0]))):
              questionlist.append("?")
            Fullfilecode[species] = Fullfilecode[species] + [questionlist]

numberOfspecies = len(fullspecieslist)
speciesName = fullspecieslist

max_len = -1
for ele in fullspecieslist:
    if len(ele) > max_len:
        max_len = len(ele)
        correction = len(ele)
 
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