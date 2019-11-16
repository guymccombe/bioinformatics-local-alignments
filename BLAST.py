import numpy as nm


def heuralign(alphabet, scoringMatrix, seq1, seq2, seedLength=3, threshold=0):
    alphabet += '-'
    highScoringSubstrings = []
    for i in range(len(seq2) - seedLength + 1):
        potentialSeed = seq2[i:i+seedLength]
        if getScoreOfMatchingStrings(alphabet, scoringMatrix, potentialSeed, potentialSeed) > threshold:
            highScoringSubstrings += [potentialSeed]

    allStringsInAlphabet = generateAllStrings(alphabet[:-1], seedLength)

    matches = []
    for str1 in highScoringSubstrings:
        for str2 in allStringsInAlphabet:
            if getScoreOfMatchingStrings(alphabet, scoringMatrix, str1, str2) > threshold:
                if str2 in seq1 and str2 not in matches:
                    matches += [str2]

    # TODO Match extension to generate alignment.


def getScoreOfMatchingStrings(alphabet, scoringMatrix, str1, str2):
    score = 0
    for i in range(len(str1)):
        score += getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, str1[i], str2[i])
    return score


def getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, character1, character2):
    rowIndex = alphabet.index(character2)
    colIndex = alphabet.index(character1)

    row = scoringMatrix[rowIndex]
    out = row[colIndex]

    return out


def generateAllStrings(alphabet, length, prefix=""):
    if length == 0:
        return [prefix]

    strings = []
    for char in alphabet:
        strings += generateAllStrings(alphabet, length - 1, prefix + char)
    return strings


a = heuralign("ABC", [[1, -1, -2, -1], [-1, 2, -4, -1], [-2, -4, 3, -2], [-1, -1, -2, 0]], "AABBAACA", "CBACCCBA")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
