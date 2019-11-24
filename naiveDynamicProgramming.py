def dynprog(alphabet, scoringMatrix, sequence1, sequence2):
    alphabet += '-'
    alignmentScoreMatrix, alignmentBacktrackMatrix = initialiseMatrices(
        alphabet, scoringMatrix, sequence1, sequence2)
    score = 0
    coordinateOfBest = (0, 0)

    for i in range(1, len(alignmentScoreMatrix)):
        for j in range(1, len(alignmentScoreMatrix[0])):
            alignmentScoreMatrix[i][j], alignmentBacktrackMatrix[i][j] = max([
                (alignmentScoreMatrix[i-1][j-1] + getScoreOfMatchingCharactersFromScoringMatrix(alphabet,
                                                                                                scoringMatrix, sequence1[i-1], sequence2[j-1]), 'D'),
                (alignmentScoreMatrix[i-1][j] + getScoreOfMatchingCharactersFromScoringMatrix(
                    alphabet, scoringMatrix, sequence1[i-1], '-'), 'U'),
                (alignmentScoreMatrix[i][j-1] + getScoreOfMatchingCharactersFromScoringMatrix(
                    alphabet, scoringMatrix, '-', sequence2[j-1]), 'L'),
                (0, 'G')
            ])

            if alignmentScoreMatrix[i][j] > score:
                score = alignmentScoreMatrix[i][j]
                coordinateOfBest = (i, j)

    sequence1Indices, sequence2Indices = backtrackAndReturnIndicesOfAlignment(
        alignmentBacktrackMatrix, coordinateOfBest)
    return score, sequence1Indices, sequence2Indices


def initialiseMatrices(alphabet, scoringMatrix, sequence1, sequence2):
    alignmentScore = [
        [-1 for _ in range(len(sequence2) + 1)] for _ in range(len(sequence1) + 1)]
    backtrackMatrix = [
        ['x' for _ in range(len(sequence2) + 1)] for _ in range(len(sequence1) + 1)]

    alignmentScore[0][0], backtrackMatrix[0][0] = 0, 'G'

    for i in range(1, len(alignmentScore)):
        alignmentScore[i][0], backtrackMatrix[i][0] = max(
            (alignmentScore[i-1][0] + getScoreOfMatchingCharactersFromScoringMatrix(
                alphabet, scoringMatrix, '-', sequence1[i-1]), 'U'),
            (0, 'G')
        )

    for i in range(1, len(alignmentScore[0])):
        alignmentScore[0][i], backtrackMatrix[0][i] = max(
            (alignmentScore[0][i-1] + getScoreOfMatchingCharactersFromScoringMatrix(
                alphabet, scoringMatrix, '-', sequence2[i-1]), 'L'),
            (0, 'G')
        )

    return alignmentScore, backtrackMatrix


def getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, character1, character2):
    rowIndex = alphabet.index(character2)
    colIndex = alphabet.index(character1)

    row = scoringMatrix[rowIndex]
    out = row[colIndex]

    return out


def backtrackAndReturnIndicesOfAlignment(backtrackMatrix, startCoordinate):
    i, j = startCoordinate
    indices1, indices2 = [], []
    currentState = backtrackMatrix[i][j]

    while (currentState != 'G'):
        if currentState == 'D':
            i -= 1
            j -= 1
            indices1.insert(0, i)
            indices2.insert(0, j)
        elif currentState == 'U':
            i -= 1
        else:
            j -= 1
        currentState = backtrackMatrix[i][j]

    return indices1, indices2


# Expected output:
# Score:    5
# Indices:  [3,5,6][1,2,3]
a = dynprog("ABC", [[1, -1, -2, -1], [-1, 2, -4, -1],
                    [-2, -4, 3, -2], [-1, -1, -2, 0]], "AABBAACA", "CBACCCBA")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])


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
a = dynprog(alphabet, matrix, "AAAAACCDDCCDDAAAAACC",
            "CCAAADDAAAACCAAADDCCAAAA")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])

# Score:    17
# Indices:  [2, 6, 11, 14, 17] [0, 1, 2, 3, 4]
a = dynprog(alphabet, matrix, "AACAAADAAAACAADAADAAA", "CDCDDD")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])

# Score:    81
# Indices:  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58] [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 61, 62, 63, 64, 65, 66, 67, 68, 69]
a = dynprog(alphabet, matrix, "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
            "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
