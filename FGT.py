import sys
from time import time
from random import choice
t0= time()

original_stdout = sys.stdout  # Save a reference to the original standard output

def HelpCommandLine():
    ContinueHelp = True
    print("you have asked for help")
    print("I will do my best but reading the online manual might be better")
    while ContinueHelp:
        print("\n","I can help with \"Flags\"(F), \"default-settings\"(D), and, \"online-support\"(O)")
        print("you can exit this any time with \"X\"")
        UserInput = input("What can I help with: ")
        if UserInput == "X":
            ContinueHelp = False
        elif UserInput == "F" or UserInput == "Flags":
            print("To use flags the file name of the input file must be placed in the command line")
            print("python3 FGT.py combineloci.nex")
            print("Flags can be typed after the file name order does not matter")
            print("python3 FGT.py combineloci.nex MS F- R")
            print()
            print("MS")
            print("F-")
            print("R")
            print("Nucleotide symbols to be ignored: \"N\" , \"?\" ")
        elif UserInput == "D" or UserInput == "default-settings":
            print("If no flags or file names are placed in the command line the program will use a combineloci.nex")
        else:
            print("I am sorry response has not been written for that")


#checks if file is openable
def inputPath(inputFilePath):
    try:
        return open(inputFilePath)
    except Exception as e:
        print('Cannot open', inputFilePath, file=sys.stderr)
        exit()
 
#confirms the nucleotide is not a type missing data 
def nucleotideKnown(nucleotide):
    return nucleotide != "-" and nucleotide != "?" and nucleotide != "*" 

#checks if nucleotide symbol is in the list of nucleotides to be ignored
def nucleotideIgnored(nucleotide):
    return nucleotide not in ignore 

def LookForMissingSequence(SortingList,numberOfSpecies):
    for species in range(numberOfSpecies):
        if "?" in SortingList[species]:
            if "A" not in SortingList[species] and "C" not in SortingList[species] and "G" not in SortingList[species] and "T" not in SortingList[species]:
                for nucleotide2 in range(len(SortingList[species])):
                    SortingList[species][nucleotide2] = "*"


def removeStar(SortingList,numberOfSpecies):
    for species in range(numberOfSpecies):
        for nucleotide2 in range(len(SortingList[species])):
            if "*" in SortingList[species]:
                SortingList[species][nucleotide2] = "?"

#creates a dictonary that shows the number of nucleotide polymorphism at each site
def InitializeRedundantNucleotideDictionary(SortingList,numberOfSpecies,nucleotide):
    NucleotideDictionary = {}
    for species in range(numberOfSpecies):
        NucleotideDictionary.update({SortingList[species][nucleotide]:species})
    return NucleotideDictionary
                    
#looks for nucleotide Polymorphism in list
def nucleotidePolymorphism(SortingList,numberOfSpecies ,multiMutation,nucleotideNumber,mutations,undefinedNucleotides):
   for nucleotide in range(len(SortingList[0])):
       NucleotideDictionary = InitializeRedundantNucleotideDictionary(SortingList,numberOfSpecies,nucleotide)
       listNucleotide3 = list(NucleotideDictionary.keys())
       for item in ignore:
            if item in listNucleotide3:
                NucleotideDictionary.clear()
                if len(undefinedNucleotides) > 0:
                    if undefinedNucleotides[-1] != (nucleotide+1):
                         undefinedNucleotides.append(nucleotide+1)
                else: 
                    undefinedNucleotides.append(nucleotide+1)
       listNucleotide3 = list(NucleotideDictionary.keys())
       if "-" in listNucleotide3:
            NucleotideDictionary.pop('-')
       if "?" in listNucleotide3:
            NucleotideDictionary.pop("?")
       if "*" in listNucleotide3:
            NucleotideDictionary.pop("*") 
       listNucleotide3 = list(NucleotideDictionary.keys())
       if len(listNucleotide3) == 2:
            nucleotideNumber.append(nucleotide+1)
            for spec in range(numberOfSpecies):
                mutations[spec].append(list(SortingList[spec])[nucleotide])
       elif len(listNucleotide3) > 2:
        multiMutation.append(nucleotide+1)
        
def nucleotideaccepted(Scan):
   if (Scan == "-" or Scan == "C" or Scan == "G" or Scan == "T" or Scan == "A" or Scan == "?" or Scan == "N"):
     return True
   else:
     return False
#this looks for the first valid line of the DNA code
def Search_Starting_Line(TestInputFile):
   global line_To_Start
   global correction
   for line2 in TestInputFile:
       Scan = list(line2)
       for letter in range(len(Scan)):
           if(len(Scan) - letter > 5) and (nucleotideaccepted(Scan[letter]) == True):
             if (nucleotideaccepted(Scan[letter + 1])) and (nucleotideaccepted(Scan[letter + 2])) and (nucleotideaccepted(Scan[letter + 4]) == True):
               Scan = list(line2)
               correction = letter
               return correction
               line_To_Start += 1
 
#this creates a dictionary of the mutations
def MutationCompression(mutations,FGTtrunk,nucleotideNumber):
   FGTcount = -1
   for section in range(len(mutations[0])):
       for section2 in range(section + 1, len(mutations[0])):
           FGTtrunk.append({})
           FGTcount += 1
           for spec2 in range(len(mutations)):
               if nucleotideKnown(mutations[spec2][section]) and nucleotideKnown(mutations[spec2][section2]):
                   Key = str(mutations[spec2][section]) + " " + str(mutations[spec2][section2])
                   #this is the place where the mutations were not included in the sections
                   FGTtrunk[FGTcount][Key] = [ nucleotideNumber[section], nucleotideNumber[section2]]
 
#this identifies recombinant events
def FourGameteTest (FGTtrunk,FGTcombine,FTGrecombinant):
   for dict in FGTtrunk:
       #if there are three or more gametes there is a recombinant event
       if len(dict) > 3:
           Key2 = list(dict.keys())
           FTGrecombinant.append(dict.get(Key2[0]))
       elif len(dict) > 0:
           Key1 = list(dict.keys())
           FGTcombine.append(dict.get(Key1[0]))

def CheckIfRegionIsInternal(BaseRegion,TestRegion):
    if BaseRegion[0] <= TestRegion[0] and TestRegion[-1]  <= BaseRegion[-1]:
        return TestRegion
    else: return BaseRegion

def CheckOverlap(templateRegion,testRegion):
    #takes in list of two numbers that label the border of the template region and checks for over lap with the list of two numbers in the test region
    if templateRegion[-1] > testRegion[0] and templateRegion[0] < testRegion[-1]:
        return True
    else: 
        return False 

def SelectShortestNonOverlapingRegions(ListFGTrecombinantRegions):
    FinalList= []
    print(ListFGTrecombinantRegions)
    BaseRegion = ListFGTrecombinantRegions[0]
    for TestRegion in ListFGTrecombinantRegions:
        if CheckOverlap(BaseRegion,TestRegion):
            BaseRegion =  CheckIfRegionIsInternal(BaseRegion,TestRegion)
        else: 
            FinalList.append(BaseRegion)
            BaseRegion = TestRegion
    if BaseRegion not in FinalList:
        FinalList.append(BaseRegion)
    return FinalList

def lengthRegion(Region):
    #takes in a list of two numbers and returns the diffrence
    return Region[-1]-Region[0]

def listToStr(list):
 strList = ' '.join(str(v) for v in list)
 return strList
#prints out spreadsheet of the data
def printData(LocusNumber,LocusName,startinglength,undefinedSites,S,infiniteSites,Rm,locationsRM,locationsRandom,longestBlock,finalLength):
 fileName = inputFilePath.split("/")[-1]
 outputFilePath = "Results Summary Table.txt"
 
 with open(outputFilePath, 'w') as f:
       sys.stdout = f
       tableFormat = "longest block"+"\t"
       if "R" in ignore:
         tableFormat = "Locations without recombinantion events"+"\t"+"Selected block"+"\t"
       print("Locus number"+"\t"+"Locus name"+"\t"+"starting length (bp)"+"\t"+"Length excluding gaps\missing data"+"\t"+"S"+"\t"+"Infinite sites no"+"\t"+"Rm"+"\t"+"Locations of recombinantion events"+"\t"+tableFormat +"final length (bp)")
       for item in range(len(LocusNumber)):
        if "R" in ignore:
         print(str(LocusNumber[item])+"\t"+LocusName[item]+"\t"+str(startinglength[item])+"\t"+str(startinglength[item]-len(undefinedSites[item]))+"\t"+str(S[item])+"\t"+infiniteSites[item]+"\t"+str(Rm[item])+"\t"+locationsRM[item]+"\t"+str(locationsRandom[item])+"\t"+longestBlock[item]+"\t"+str(finalLength[item]))
        else:
         print(str(LocusNumber[item])+"\t"+LocusName[item]+"\t"+str(startinglength[item])+"\t"+str(startinglength[item]-len(undefinedSites[item]))+"\t"+str(S[item])+"\t"+infiniteSites[item]+"\t"+str(Rm[item])+"\t"+locationsRM[item]+"\t"+longestBlock[item]+"\t"+str(finalLength[item]))
 
       sys.stdout = original_stdout

def writeRawResults(outputFilePath,LocusName,locus,charset,outPutString,rawResultsString):
    rawResultsString.append(str(LocusName[locus]) + "\n" + "Loci length:" + str(1+int(charset[locus][1])-int(charset[locus][0])) + "\n" + "\n" + outPutString + "\n" + "\n")

def printRawResults(outputFilePath,rawResultsString):
    with open(outputFilePath, 'a') as f:
        sys.stdout = f  # Change the standard output to the file we created.
        for line in rawResultsString:
            print(line)
    sys.stdout = original_stdout

def SetUp(inputFilePath):
 
   fileList = [line.rstrip('\n') for line in inputPath(inputFilePath)]
 
   TestInputFile = inputPath(inputFilePath)
 
   TestInputFile2 = inputPath(inputFilePath)
 
   numberOfLines = 0
 
   global correction
   correction = 0
   global line_To_Start
   line_To_Start = 0
 
   charset = []
 
   LocusNumber = []
   LocusName  = []
   startinglength = []
   S = []
   infiniteSites = []
   undefinedSites =[]
   Rm = []
   locationsRM = []
   locationsRandom = []
   longestBlock = []
   finalLength = []

   rawResultsString = []
 
   #checks if the file is a Nexus file
   if fileList[0] == "#NEXUS":
       line = 0
       while (fileList[line]).upper() != "BEGIN DATA;" and line < 10:
           line += 1
       commandline = fileList[line + 1].split()
       commandlineSplit = commandline[1].split("=")
       numberOfSpecies  = int(commandlineSplit[1])
       commandlineSplit = commandline[2].split("=")
       nchar = int(commandlineSplit[1].strip(';'))
   if fileList[0] != "#NEXUS":
       Search_Starting_Line()
   else:
       while (fileList[line]).upper().strip() != "MATRIX":
           line += 1
       line_To_Start = line + 1
   if fileList[0] == "#NEXUS":
       line = -1
       while (fileList[line]) != ";":
           line -= 1
           if fileList[line].upper().strip() == "BEGIN SETS;":
               for i in range(len(fileList) + line, len(fileList) - 1):
                   lineToCheck = fileList[i].strip(";").split()
                   if len(lineToCheck) > 0 and lineToCheck[0] == "charset":
                       if len(lineToCheck) == 4:
                         charset.append(lineToCheck[-1].split("-"))
                       elif len(lineToCheck) == 6:
                         charset.append([lineToCheck[-3],lineToCheck[-1]])
                       LocusName.append(lineToCheck[1])
                  
       if len(charset) < 1:
           charset.append([1, nchar])
           LocusName.append("Loci_1")
 
 
   numberOfLines = line_To_Start
   Search_Starting_Line(TestInputFile)
   SortingList = []
   Fullfilecode = []
   speciesName = []
   global Fullfilecode2
   Fullfilecode2 = [] 
   #initializing lists
   for spec3 in range(numberOfSpecies ):
       Fullfilecode.append([])
       speciesName.append([])
       Fullfilecode2.append([])
   firstemptyline = True
   #counts the number of species
   if fileList[0] != "#NEXUS":
       for line in TestInputFile2:
           numberOfLines += 1
           if numberOfLines >= 0 and len(line) <= 3 and firstemptyline:
               numberOfSpecies  = numberOfLines
               firstemptyline = False
               break
 
   firstemptyline = True
 
     #Reads the input file breaks it down in to list of nucleotides
   while fileList[numberOfLines] != ";":
       line = fileList[numberOfLines]
       numberOfLines += 1
       if numberOfLines >= 0 and len(line) > 4:
           SortingList.append(line)
           if len(SortingList) > numberOfSpecies :
               SortingList.clear()
               SortingList.append(line)
           if len(SortingList) == numberOfSpecies :
               for code in range(numberOfSpecies ):
                   dnaOnly = list(
                       SortingList[code])[correction:len(list(SortingList[code]))]
                   Fullfilecode[code].extend(dnaOnly)
               if firstemptyline:
                   for code5 in range(numberOfSpecies ):
                       dnaOnly = list(SortingList[code5])[0:correction]
                       speciesName[code5].extend(dnaOnly)
                   firstemptyline = False
 
   Locis = [[]]
   locusNumber = 0
 
   #BREAKS stream of nucleotide into separate Loci
   for locus in charset:
       for species in range(numberOfSpecies ):
           Locis[locusNumber].append(Fullfilecode[species][int(locus[0]) - 1:int(locus[-1])])
       locusNumber += 1
       Locis.append([])
 
   outputFilePath = "Raw Results Summary.txt"
 
   #truncates the Loci
   for locus in range(len(charset)):
       LocusNumber.append(locus+1)
       print(LocusName[locus])
       nucleotideNumber = []
       mutations = []
       multiMutation = []
       undefinedNucleotides = []
       for spec3 in range(numberOfSpecies ):
         mutations.append([])
 
       if ms == True:
            LookForMissingSequence(Locis[locus],numberOfSpecies)

       nucleotidePolymorphism(Locis[locus],numberOfSpecies ,multiMutation,nucleotideNumber,mutations,undefinedNucleotides)
       
       removeStar(Locis[locus],numberOfSpecies)

       if len(multiMutation) > 0:
         ListMultiMutation = [str(mutation) for mutation in multiMutation]
         stringMultiMutation = ", ".join(ListMultiMutation)
       else:
         stringMultiMutation = "Zero"
 
       infiniteSites.append(stringMultiMutation)
       undefinedSites.append(undefinedNucleotides)
       S.append(str(len(nucleotideNumber)+len(multiMutation)))
 
       outPutString = "sites that violate infinite sites assumption: " + stringMultiMutation + "\n" + "single nucleotide Polymorphism: " + str(len(nucleotideNumber)) + "\n" + listToStr(nucleotideNumber) + "\n"
 
       FGTtrunk = []
 
       MutationCompression(mutations,FGTtrunk,nucleotideNumber)

       FGTcombine = []
       FGTrecombinants = []
 
       FourGameteTest (FGTtrunk,FGTcombine,FGTrecombinants)

       outPutString += "Number of pairs of sites with four gametic types: " + str(len(FGTrecombinants)) + "\n" + listToStr(FGTrecombinants)
  
       empty = 0
       LongFGTrecombinants = []
       
       for Region in FGTrecombinants:
            if lengthRegion(Region) >= lengthRegionFilter and Region != [0,0]:
                LongFGTrecombinants.append(Region)
       FGTrecombinants = LongFGTrecombinants
 
       #checks if there are recombinant
       if len(FGTrecombinants) < 1:
           FGTrecombinants.append([0, 0])
           empty = 1

       FinalRM =SelectShortestNonOverlapingRegions(FGTrecombinants)

       outPutString += "\n" "Rm " + str(len(FinalRM) - empty) + "\n"
       outPutString += str(FinalRM)

       if len(FinalRM) > 1 and FinalRM[0] == [0,0]:
        FinalRM.pop(0)
 
       XloutPutString = ""

       Rm.append(str(len(FinalRM) - empty ))

       XloutPutString += str(FinalRM)
  
       locationsRM.append(XloutPutString)
       Output = [[1,FinalRM[0][0] - 1]]

       for recombinant in range(len(FinalRM)):
           if recombinant == len(FinalRM) - 1:
               Output.append([FinalRM[recombinant][-1] + 1,len(Locis[locus][0])])
           elif FinalRM[recombinant][-1] != FinalRM[recombinant + 1][0] and FinalRM[recombinant][-1] + 1 != FinalRM[recombinant][0]:
               Output.append([FinalRM[recombinant][-1] + 1,FinalRM[recombinant + 1][0] - 1])

       final = [0, 0]

       #picks random section of DNA without a recombinant region
       if "R" in ignore:
            if Output[0][-1] < 0:
                Output = Output[1:]
            final = choice(Output)
       else:
        for lists in Output:
           if (lists[0] - lists[-1]) < (final[0] - final[-1]):
               final = lists
        if final[0]== 0:
           final[0]=1

       if final[0]!=1:
         final[0]= final[0]-1

       if final[-1] != 1+int(charset[locus][1])-int(charset[locus][0]):
         final[-1]= final[-1]+1

       locationsRandom.append(Output)
       outPutString += "\n" + "final" + str(final) + "\n" + "final truncation length: " + str(final[1]-final[0]+1)
       longestBlock.append(str(final))
       finalLength.append(str(str(final[1]-final[0]+1)))
       startinglength.append(1+int(charset[locus][1])-int(charset[locus][0]))
       writeRawResults(outputFilePath,LocusName,locus,charset,outPutString,rawResultsString)
  
       for spec6 in range(numberOfSpecies ):
               Fullfilecode2[spec6].append(Locis[locus][spec6][final[0]-1:final[-1]])
 
   fileName = inputFilePath.split("/")[-1]
  
   outputFilePath = "Trunc_" + fileName[0:-4] + ".nex"
 
   ncha = 0
 
   with open(outputFilePath, 'w') as f:
       for locus in range(len(charset)):
           ncha += len(Fullfilecode2[0][locus]) 
       sys.stdout = f  # Change the standard output to the file we created.
       print("#NEXUS")
       print("begin data;")
       print("\tdimensions ntax=" + str(numberOfSpecies ) + " nchar=" + str(ncha) + ";")
       print("\tformat datatype=dna missing=? gap=- interleave;")
       print("matrix")
       for locus in range(len(charset)):
           for i in range(0, len(Fullfilecode2[0][locus]), 70):
               for species in range(numberOfSpecies ):
                   row = Fullfilecode2[species][locus][i:i + 70]
                   name = speciesName[species][0:correction]
                   print(''.join(name), " ", ''.join(row))
               print()
       print(";")
       print("end;")
       if len(charset) > 1:
           print("#nexus")
           print("begin sets;")
           outputlength = 0
           for locus in range(len(charset)):
               print("charset " + LocusName[locus] + " = " + str(outputlength+1) + "-" + str(len(Fullfilecode2[0][locus])+outputlength) + ";")
               outputlength += len(Fullfilecode2[0][locus])
           print("end;")
 
       sys.stdout = original_stdout  # Reset the standard output to its original value
 
   outputFilePath = "Trunc_" + fileName[0:-4] + ".phy"
 
   count = 0
 
   with open(outputFilePath, 'w') as f:
       sys.stdout = f
       for locus in range(len(charset)):
           numbercount = 1
           print()
           print(" ",numberOfSpecies ,"     ",len(Fullfilecode2[species][locus]))
           print()
           for species in range(numberOfSpecies ):
               count += 1
               row = Fullfilecode2[species][locus][0:]
               name = LocusName[locus].strip("'")
               name += "^"
               name += ''.join(speciesName[species]).strip()
               numbercount +=1
               print(''.join(name))
               print(''.join(row))
           print()
 
       sys.stdout = original_stdout
 
   printData(LocusNumber,LocusName,startinglength,undefinedSites,S,infiniteSites,Rm,locationsRM,locationsRandom,longestBlock,finalLength)
   outputFilePath = "Raw Results Summary.txt"

if __name__ == "__main__":
    ignore = []
    ms = False
    lengthRegionFilter = 0

    #comand line arguments
    if len(sys.argv) == 2:
        inputFilePath = sys.argv[1]
        if inputFilePath.lower() == "help":
            HelpCommandLine()
    elif len(sys.argv) > 2:
        if "help" in sys.argv or "HELP" in sys.argv or "Help" in sys.argv:
            HelpCommandLine()
        inputFilePath = sys.argv[1]
        flags=sys.argv[2:]
        ignore = list(map(lambda x: x.upper(),flags))
        if "MS" in ignore:
            ms = True
        if "F-" in ignore:
            try:
                lengthRegionFilter = int(ignore[ignore.index("F-")+1])
            except:
                print("Could not identify Filter length")
    else:
        inputFilePath = "combineloci.nex"

    SetUp(inputFilePath)

    t1 = time() 
    # get the execution time
    elapsed_time = t1 - t0
    timeInMinutes = int(elapsed_time//60)
    timeInSeconds = elapsed_time % 60 
    print('Run time:', str(timeInMinutes), 'minutes' , str(timeInSeconds)[0:3], 'seconds')