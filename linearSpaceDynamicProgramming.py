import numpy as np

def dynproglin(alphabet, scoringMatrix, sequence1, sequence2):
    alphabet += '_'
    _, end, _ = twoColumnSearch(alphabet, scoringMatrix, sequence1, sequence2)

    reversed1, reversed2 = trimAndFlipSequences(sequence1, sequence2, end)
    score, start, _ = twoColumnSearch(alphabet, scoringMatrix, reversed1, reversed2)

    trimmedSequence1, trimmedSequence2 = trimAndFlipSequences(reversed1, reversed2, start)
    start = len(reversed1) - start[0] - 1, len(reversed2) - start[1] - 1

    indices = hirschberg(alphabet, scoringMatrix, trimmedSequence1, trimmedSequence2, start[0], start[1])
    return score, indices[0], indices[1]

# TODO Investigate why this is returning incorrect scores for longer sequences.
def twoColumnSearch(alphabet, scoringMatrix, sequence1, sequence2):
    column = [0 for i in range(len(sequence2) + 1)]
    bestScore = 0
    coordinateOfBest = (0,0)

    for i in range(len(sequence1)):
        lastCol = column
        column = [0]
        for j in range(1, len(sequence2) + 1):
            score = max([
                lastCol[j-1] + getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, sequence1[i], sequence2[j-1]),
                column[j-1]  + getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, sequence1[i], '_'),
                lastCol[j]   + getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, '_', sequence2[j-1]),
                0
            ])
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

def hirschberg(alphabet, scoringMatrix, seq1, seq2, offsetX = 0, offsetY = 0):
    print(seq1, seq2)

    if len(seq1) == 0 or len(seq2) == 0:
        return [], []

    if len(seq1) == 1:
        if (seq1[0] in seq2):
            return [offsetX], [offsetY + seq2.index(seq1[0])]
        else:
            return [], []

    if len(seq2) == 1:
        if (seq2[0] in seq1):
            return [offsetX + seq1.index(seq2[0])], [offsetY]
        else:
            return [], []
    
    else:
        midX = len(seq1) // 2

        _, _, scoreL = twoColumnSearch(alphabet, scoringMatrix, seq1[:midX], seq2)
        _, _, scoreR = twoColumnSearch(alphabet, scoringMatrix, seq1[midX:][::-1], seq2[::-1])

        sumOfScores = [x+y for x, y in zip(scoreL, scoreR[::-1])]
        midY = np.argmax(sumOfScores)

        left = hirschberg(alphabet, scoringMatrix, seq1[:midX], seq2[:midY], offsetX, offsetY)
        right = hirschberg(alphabet, scoringMatrix, seq1[midX:], seq2[midY:], offsetX + midX, offsetY + midY)
        return left[0] + right[0], left[1] + right[1]

a = dynproglin ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

## Expected output:
# Score:    5
# Indices:  [3,5,6][1,2,3]


alphabet = "ABCD"
matrix = [
    [ 1,-5,-5,-5,-1],
    [-5, 1,-5,-5,-1],
    [-5,-5, 5,-5,-4],
    [-5,-5,-5, 6,-4],
    [-1,-1,-4,-4,-9]
]

# Score:    39
# Indices:  [5, 6, 7, 8, 9, 10, 11, 12, 18, 19] [0, 1, 5, 6, 11, 12, 16, 17, 18, 19]
a = dynproglin(alphabet, matrix, "AAAAACCDDCCDDAAAAACC", "CCAAADDAAAACCAAADDCCAAAA")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

# Score:    17
# Indices:  [2, 6, 11, 14, 17] [0, 1, 2, 3, 4]
a = dynproglin(alphabet, matrix, "AACAAADAAAACAADAADAAA", "CDCDDD")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

# Score:    81
# Indices:  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58] [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 61, 62, 63, 64, 65, 66, 67, 68, 69]
a = dynproglin(alphabet, matrix, "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD", "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])
