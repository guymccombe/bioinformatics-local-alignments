import numpy as np


def dynproglin(alphabet, scoringMatrix, sequence1, sequence2):
    if len(sequence1) < len(sequence2):
        out = dynproglin(alphabet, scoringMatrix, sequence2, sequence1)
        return out[0], out[2], out[1]

    alphabet += '-'
    _, end, _ = twoColumnSearch(alphabet, scoringMatrix, sequence1, sequence2)

    reversed1, reversed2 = trimAndFlipSequences(sequence1, sequence2, end)
    score, start, _ = twoColumnSearch(
        alphabet, scoringMatrix, reversed1, reversed2)

    trimmedSequence1, trimmedSequence2 = trimAndFlipSequences(
        reversed1, reversed2, start)
    start = len(reversed1) - start[0] - 1, len(reversed2) - start[1] - 1

    indices = hirschberg(alphabet, scoringMatrix,
                         trimmedSequence1, trimmedSequence2, start[0], start[1])
    return score, indices[0], indices[1]


def twoColumnSearch(alphabet, scoringMatrix, sequence1, sequence2, isLocal=True):
    bestScore = 0
    coordinateOfBest = (0, 0)
    column = [0 for _ in range(len(sequence2) + 1)]

    for i in range(len(sequence1)):
        lastCol = column
        column = [0]

        for j in range(1, len(sequence2) + 1):
            upLeft = lastCol[j-1] + getScoreOfMatchingCharactersFromScoringMatrix(
                alphabet, scoringMatrix, sequence1[i], sequence2[j-1])
            up = column[j-1] + getScoreOfMatchingCharactersFromScoringMatrix(
                alphabet, scoringMatrix, '-', sequence2[j-1])
            left = lastCol[j] + getScoreOfMatchingCharactersFromScoringMatrix(
                alphabet, scoringMatrix, sequence1[i], '-')
            if (isLocal):
                score = max(upLeft, up, left, 0)
            else:
                score = max(upLeft, up, left)
            column.append(score)
            if score > bestScore:
                bestScore = score
                coordinateOfBest = (i, j-1)
    return bestScore, coordinateOfBest, column


def getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, character1, character2):
    rowIndex = alphabet.index(character2)
    colIndex = alphabet.index(character1)

    row = scoringMatrix[rowIndex]
    out = row[colIndex]

    return out


def trimAndFlipSequences(seq1, seq2, trimPoint):
    s1 = ''.join(reversed(seq1[:trimPoint[0]+1]))
    s2 = ''.join(reversed(seq2[:trimPoint[1]+1]))
    return s1, s2


def hirschberg(alphabet, scoringMatrix, seq1, seq2, offsetX=0, offsetY=0, isLeftTree=True):
    if len(seq1) == 0 or len(seq2) == 0:
        return [], []

    if len(seq1) == 1:
        if seq1 in seq2:
            return [offsetX], [offsetY + seq2.index(seq1)]
        else:
            return [], []

    if len(seq2) == 1:
        if seq2 in seq1:
            return [offsetX + seq1.index(seq2)], [offsetY]
        else:
            return [], []

    else:
        midX = len(seq1) // 2

        _, _, scoreL = twoColumnSearch(
            alphabet, scoringMatrix, seq1[:midX], seq2, isLocal=False)
        _, _, scoreR = twoColumnSearch(
            alphabet, scoringMatrix, seq1[midX:][::-1], seq2[::-1], isLocal=False)

        sumOfScores = [x+y for x, y in zip(scoreL, scoreR[::-1])]
        if isLeftTree:
            reverseSum = sumOfScores[::-1]
            midY = len(reverseSum) - np.argmax(reverseSum) - 1
        else:
            midY = np.argmax(sumOfScores)

        left = hirschberg(alphabet, scoringMatrix,
                          seq1[:midX], seq2[:midY], offsetX, offsetY)
        right = hirschberg(alphabet, scoringMatrix,
                           seq1[midX:], seq2[midY:], offsetX + midX, offsetY + midY, isLeftTree=False)
        return left[0] + right[0], left[1] + right[1]


alphabet = "ABCD"
matrix = [
    [1, -5, -5, -5, -1],
    [-5, 1, -5, -5, -1],
    [-5, -5, 5, -5, -4],
    [-5, -5, -5, 6, -4],
    [-1, -1, -4, -4, -9]
]

# Score:    39
# Indices:  [5, 6, 7, 8, 9, 10, 11, 12, 18, 19] [0, 1, 5, 6, 11, 12, 16, 17, 18, 19]
a = dynproglin(alphabet, matrix, "AAAAACCDDCCDDAAAAACC",
               "CCAAADDAAAACCAAADDCCAAAA")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])

# Score:    17
# Indices:  [2, 6, 11, 14, 17] [0, 1, 2, 3, 4]
a = dynproglin(alphabet, matrix, "AACAAADAAAACAADAADAAA", "CDCDDD")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])

# Score:    81
# Indices:  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58] [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 61, 62, 63, 64, 65, 66, 67, 68, 69]
a = dynproglin(alphabet, matrix, "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
               "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
