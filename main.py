AminoAcidsWeights = {}  # dictionary contains the amino acids and their weights #key:AA value:weight of AA
AminoAcidsWeightsForLinearSpectrum = {}
SpectrumWeights = []
with open("weight.txt") as file:
    for line in file:
        line = line.rstrip()
        RowSpectrumWeight = line.split(" ")
        AminoAcidsWeights[int(RowSpectrumWeight[1])] = RowSpectrumWeight[0]
        # RowSpectrumWeight[0]=key=AA    #RowSpectrumWeight[1]=value=weight of AA
        AminoAcidsWeightsForLinearSpectrum[RowSpectrumWeight[0]] = int(RowSpectrumWeight[1])

with open("spectrum-weights.txt") as RetrieveSpectrumWeights:
    for line in RetrieveSpectrumWeights:
        SpectrumWeights.append(int(line.rstrip()))


def GetInitialList(spectrum):
    # InitialListWeights = list(filter(lambda weight: weight <= 186 and weight != 0, spectrum))
    # InitialListAmino = [AminoAcidsWeights[it] for it in InitialListWeights]
    InitialList = []
    i = 0
    # Find the initial List by comparing the weights to the spectrum
    while i < len(spectrum):
        for key, value in AminoAcidsWeights.items():
            if key == spectrum[i]:
                InitialList.append(value)
        i += 1
    # return the initial list that contains the characters to form the subPeptides
    return InitialList


def LinearSpectrum(peptides):
    LenPeptide = len(peptides)
    SumWeights = 0
    weights = []
    peptidesList = []
    # Get the weight of the subPeptide than is less than 3 characters by the sum of the two characters of the subPeptide
    if LenPeptide < 3:
        for peptide in peptides:
            SumWeights = SumWeights + AminoAcidsWeightsForLinearSpectrum[peptide]
        weights.append(SumWeights)
    else:
        # Get the weight of every possible combinations of the subPeptide if its lengths is more than 3 characters
        # Make all the possible Combinations and store them in the list
        for i in range(LenPeptide):
            for j in range(LenPeptide):
                if peptides[i: i + j + 1] not in peptidesList:
                    peptidesList.append(peptides[i: i + j + 1])
        # find the weights of all the combinations by the sum of weights of the characters in each combination
        for peptide in peptidesList:
            SumWeights = 0
            for subPeptide in peptide:
                SumWeights = SumWeights + AminoAcidsWeightsForLinearSpectrum[subPeptide]
            weights.append(SumWeights)
    # return the list of weights of all the combination of the subPeptide
    listOfWeights = sorted(weights)
    return listOfWeights


def IsConsistent(subPeptide, spectrum):
    found = []
    weights = LinearSpectrum(subPeptide)
    TheoreticalSpectrum = spectrum.copy()
    # Compare the weights returned from the linear spectrum with that in the Theoretical Spectrum
    for weight in weights:
        # if it is found in the spectrum, 1 is added to the found list
        if weight in TheoreticalSpectrum:
            found.append(1)
            TheoreticalSpectrum.remove(weight)
        # if it is not found in the spectrum, 0 is added to the found list
        else:
            found.append(0)
    # if there is any 0 in the found list, then this subPeptide is not consistent
    # if the found list contains 1, then this subPeptide is consistent
    IsConsistentPeptide = 'false'
    for f in found:
        if f == 1:
            IsConsistentPeptide = 'true'
        elif f == 0:
            IsConsistentPeptide = 'false'
            break
    # return if this subPeptide is consistent or not
    return IsConsistentPeptide


def MainFunction(spectrum):
    InitList = GetInitialList(spectrum)
    TempList = InitList.copy()
    items = []
    AllLinearRepresentation = []
    flag = 0

    # loop in the initial list, and make a combination with the TempList that contains a copy of the initial list
    for f in range(len(InitList) - 1):
        for x in range(0, len(InitList)):
            for z in range(len(TempList)):
                # make the subPeptide by combining the one in the initial list with that in the TempList
                item = TempList[z] + InitList[x]
                if item not in items:
                    items.append(item)

                # this loop finds if there is 2 repeated characters beside each other
                for i in range(len(item) - 1):
                    if item[i] == item[i + 1]:
                        flag = 1
                    else:
                        flag = 0
                # if there is any repeated characters, then the item is removed from the list of items
                if flag == 1:
                    items.remove(item)

        TempList.clear()
        # this loop finds if the item is consistent or not bya calling the is_consistent function
        for it in items:
            is_consistent = IsConsistent(it, spectrum=spectrum)
            # if it is consistent, then it is added to the TempList to be used in the next iteration
            if is_consistent == 'true':
                TempList.append(it)
                if len(it) == len(InitList):
                    AllLinearRepresentation.append(it)
        items.clear()
        if len(TempList) == 0:
            break
    # return the list that contains all the linear representations of the cyclic sequence of the
    # protein (Peptide Sequence)
    return ' '.join(AllLinearRepresentation)


if __name__ == '__main__':
    # calling the main function and print the list of that contains all the linear representations
    LinearRepresentation = MainFunction(spectrum=SpectrumWeights)
    print(LinearRepresentation)
