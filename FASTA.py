import numpy as np


def heuralign(alphabet, scoringMatrix, seq1, seq2, seedLength=3, bandWidth=32, greedyExpansionThreshold=0):
    seeds = generateSeeds(seedLength, seq1, seq2)
    if seeds.size == 0:
        if seedLength > 1:
            return heuralign(alphabet, scoringMatrix, seq1, seq2, seedLength - 1)
        else:
            return (0, [], [])
    alphabet += '-'

    expandedSeeds = expandSeeds(
        seeds, alphabet, scoringMatrix, (seq1, seq2), greedyExpansionThreshold)
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


def expandSeeds(seeds, alphabet, scoringMatrix, sequences, greedyExpansionThreshold):
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

    if diagonal < 0:
        i = -1 - width
        j = -diagonal - width // 2
    else:
        i = diagonal - width // 2
        j = -1

    smithWatermanDict = {}
    currentBestScore = -1
    coordOfBestScore = (-1, -1)
    while i < len(seq1) and j < len(seq2):
        if j > -2:
            for k in range(width):
                smithWatermanDict[(i+k, j)] = align(
                    i+k, j, smithWatermanDict, alphabet, scoringMatrix, sequences)

                if smithWatermanDict[(i+k, j)]["score"] > currentBestScore:
                    currentBestScore = smithWatermanDict[(i+k, j)]["score"]
                    coordOfBestScore = (i+k, j)
        i += 1
        j += 1

    indices1, indices2 = backtrackFrom(smithWatermanDict, coordOfBestScore)
    return currentBestScore, indices1, indices2


def align(i, j, scores, alphabet, scoringMatrix, sequences):
    seq1, seq2 = sequences

    # ie index out of range
    if i < -1 or i >= len(seq1):
        return {"score": -1, "direction": 'G'}

    # ie at position (0,0) in a smith-waterman matrix
    if j == -1 and i == -1:
        return {"score": 0, "direction": 'G'}

    parameters = i, j, scores, alphabet, scoringMatrix, seq1, seq2

    # general case
    maxAsTuple = max(left(*parameters), up(*parameters),
                     diag(*parameters), (0, 'G'))
    return {"score": maxAsTuple[0], "direction": maxAsTuple[1]}


def left(i, j, scores, alphabet, scoringMatrix, seq1, seq2):
    return (getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, '-', seq1[i]) + scores.get((i-1, j), {"score": float("-inf")})["score"], 'L')


def up(i, j, scores, alphabet, scoringMatrix, seq1, seq2):
    return (getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, seq2[j], '-') + scores.get((i, j-1), {"score": float("-inf")})["score"], 'U')


def diag(i, j, scores, alphabet, scoringMatrix, seq1, seq2):
    return (getScoreOfMatchingCharactersFromScoringMatrix(
            alphabet, scoringMatrix, seq2[j], seq1[i]) + scores.get((i-1, j-1), {"score": float("-inf")})["score"], 'D')


def backtrackFrom(dictionary, startPoint):
    i, j = startPoint
    indices1, indices2 = [], []
    state = dictionary[startPoint]["direction"]

    while state != 'G':
        if state == 'D':
            indices1.insert(0, i)
            indices2.insert(0, j)
            i -= 1
            j -= 1
        elif state == 'U':
            j -= 1
        else:
            i -= 1
        state = dictionary[(i, j)]["direction"]

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
