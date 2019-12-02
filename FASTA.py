import numpy as np


def heuralign(alphabet, scoringMatrix, seq1, seq2, seedLength=3, bandWidth=32):
    seeds = generateSeeds(seedLength, seq1, seq2)
    if seeds.size == 0:
        if seedLength > 1:
            return heuralign(alphabet, scoringMatrix, seq1, seq2, seedLength - 1)
        else:
            return (0, [], [])

    alphabet += '-'

    expandedSeeds = expandSeeds(seeds, alphabet, scoringMatrix, (seq1, seq2))
    expandedSeeds.view("int8, uint8, uint8, uint8").sort(order=['f0'], axis=0)
    highestScoringSeeds = np.unique(expandedSeeds)[-10:][::-1]

    return bandedSmithWaterman(highestScoringSeeds, bandWidth, alphabet, scoringMatrix, (seq1, seq2))


def generateSeeds(seedLength, seq1, seq2):
    seeds = []

    for i in range(len(seq1) - seedLength + 1):
        substring = seq1[i:i+seedLength]
        for j in range(i, len(seq2) - seedLength + 1):
            if substring == seq2[j:j+seedLength]:
                seeds += [(0, i, j, seedLength)]

    return np.array(seeds, dtype="int8, uint8, uint8, uint8")


def expandSeeds(seeds, alphabet, scoringMatrix, sequences, greedyExpansionThreshold=0):
    seq1, seq2 = sequences
    for i in range(len(seeds)):
        scoreOfSeed = 0
        seed = seeds[i]
        for j in range(seed[1], seed[1] + seed[3]):
            scoreOfSeed += getScoreOfMatchingCharactersFromScoringMatrix(
                alphabet, scoringMatrix, seq1[j], seq1[j])
        seed[0] = scoreOfSeed

        # Expand forwards
        seed = expandInDirection(
            seed, +1, alphabet, scoringMatrix, sequences, greedyExpansionThreshold)

        # Expand backwards
        seed = expandInDirection(
            seed, -1, alphabet, scoringMatrix, sequences, greedyExpansionThreshold)

        seeds[i] = seed

    return seeds


def expandInDirection(seed, direction, alphabet, scoringMatrix, sequences, threshold):
    seq1, seq2 = sequences
    if direction > 0:
        i, j = seed[1] + seed[3], seed[2] + seed[3]
    else:
        i, j = seed[1] - 1, seed[2] - 1

    score = seed[0]
    maxExpansion = (seed[0], -1)
    while i > -1 and i < len(seq1) and j > -1 and j < len(seq2) and score > threshold:
        score += getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, seq1[i], seq2[j])

        if score > maxExpansion[0]:
            maxExpansion = (score, i)

        i += direction
        j += direction

    if maxExpansion[0] > seed[0]:
        if direction > 0:
            seed[3] = maxExpansion[1] - seed[1] + 1
            seed[0] = maxExpansion[0]
        else:
            seed[3] = seed[1] + seed[3] - maxExpansion[1]
            seed[2] -= seed[1] - maxExpansion[1]
            seed[1] = maxExpansion[1]
            seed[0] = maxExpansion[0]

    return seed


def getScoreOfMatchingCharactersFromScoringMatrix(alphabet, scoringMatrix, character1, character2):
    rowIndex = alphabet.index(character2)
    colIndex = alphabet.index(character1)

    row = scoringMatrix[rowIndex]
    out = row[colIndex]

    return out


def bandedSmithWaterman(seeds, width, alphabet, scoringMatrix, sequences):
    seq1, seq2 = sequences[0], sequences[1]

    diagonal, sumOfScores = 0, 0
    for seed in seeds:
        # int() converts from unsigned
        diagonal += seed[0] * (seed[1] - int(seed[2]))
        sumOfScores += seed[0]
    diagonal = diagonal // sumOfScores

    alignmentScoreMatrix = np.zeros((len(seq2) + 1, width), dtype="int8")
    indexBacktrackMatrix = np.chararray((len(seq2) + 1, width))

    currentBestScore = -1
    coordOfBestScore = (-1, -1)
    for j in range(-1, len(seq2)):
        for i in range(width):
            alignmentScoreMatrix[j+1, i], indexBacktrackMatrix[j+1, i] = align(
                j, i, diagonal, width, alignmentScoreMatrix, alphabet, scoringMatrix, sequences)

            if alignmentScoreMatrix[j+1, i] > currentBestScore:
                currentBestScore = alignmentScoreMatrix[j+1, i]
                coordOfBestScore = (j+1, i)

    indices1, indices2 = backtrackFrom(indexBacktrackMatrix,
                                       coordOfBestScore, diagonal, width)
    return currentBestScore, indices1, indices2


def align(j, i, diagonal, width, scores, alphabet, scoringMatrix, sequences):
    seq1, seq2 = sequences
    seq1Index = j + diagonal + (i - width // 2)

    # ie index out of range
    if seq1Index < -1 or seq1Index >= len(seq1):
        return -1, 'X'

    # ie at position (0,0) in a smith-waterman matrix
    if j == -1 and seq1Index == -1:
        score = getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, '-', '-')
        return max(0, score), 'G'

    parameters = i, j, seq1Index, scores, width, alphabet, scoringMatrix, seq1, seq2

    # ie the first row of a smith-waterman matrix
    if j == -1:
        return max(left(*parameters), (0, 'G'))

    # ie the first column of a smith-waterman matrix
    if seq1Index == -1:
        return max(up(*parameters), (0, 'G'))

    # general case
    return max(left(*parameters), up(*parameters), diag(*parameters), (0, 'G'))


def left(i, j, seq1Index, scores, width, alphabet, scoringMatrix, seq1, seq2):
    if i == 0:  # ie no value to the left of this one
        return -1, 'X'

    return (getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, '-', seq1[seq1Index]) + scores[j+1, i-1], 'L')


def up(i, j, seq1Index, scores, width, alphabet, scoringMatrix, seq1, seq2):
    if i == width-1:  # ie no value above this one
        return -1, 'X'

    return (getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, seq2[j], '-') + scores[j, i+1], 'U')


def diag(i, j, seq1Index, scores, width, alphabet, scoringMatrix, seq1, seq2):
    return (getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, seq2[j], seq1[seq1Index]) + scores[j, i], 'D')


def backtrackFrom(matrix, startPoint, diagonal, width):
    j, i = startPoint
    seq1index = j + diagonal + (i - width // 2)
    seq2index = j
    indices1, indices2 = [], []
    state = matrix[startPoint]

    while state != b'G':
        if state == b'D':
            j -= 1
            seq1index -= 1
            seq2index -= 1
            indices1.insert(0, seq1index)
            indices2.insert(0, seq2index)
        elif state == b'U':
            j -= 1
            i += 1
            seq2index -= 1
        else:
            i -= 1
            seq1index -= 1

        state = matrix[(j, i)]

    return indices1, indices2


alphabet = "ABCD"
matrix = [
    [1, -5, -5, -5, -1],
    [-5, 1, -5, -5, -1],
    [-5, -5, 5, -5, -4],
    [-5, -5, -5, 6, -4],
    [-1, -1, -4, -4, -9]


]

# Expected output:
# Score:    5
# Indices:  [3,5,6][1,2,3]
a = heuralign("ABC", [[1, -1, -2, -1], [-1, 2, -4, -1],
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
a = heuralign(alphabet, matrix, "AAAAACCDDCCDDAAAAACC",
              "CCAAADDAAAACCAAADDCCAAAA")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])

# Score:    17
# Indices:  [2, 6, 11, 14, 17] [0, 1, 2, 3, 4]
a = heuralign(alphabet, matrix, "AACAAADAAAACAADAADAAA", "CDCDDD")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])

# Score:    81
# Indices:  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58] [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 61, 62, 63, 64, 65, 66, 67, 68, 69]
a = heuralign(alphabet, matrix, "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
              "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])
