def dynproglin(alphabet, scoringMatrix, sequence1, sequence2):
    alphabet += '_'
    score, end, _ = twoColumnSearch(alphabet, scoringMatrix, sequence1, sequence2)

    reversed1, reversed2 = trimAndFlipSequences(sequence1, sequence2, end)
    _, start, _ = twoColumnSearch(alphabet, scoringMatrix, reversed1, reversed2)

    trimmedSequence1, trimmedSequence2 = trimAndFlipSequences(reversed1, reversed2, start)
    indices = hirschberg(alphabet, scoringMatrix, trimmedSequence1, trimmedSequence2)
    return score, indices[0], indices[1]

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

def hirschberg(alphabet, scoringMatrix, seq1, seq2):
    print(seq1, seq2)
    if len(seq1) < 2 or len(seq2) < 2:
        return [], []

    else:
        midX = len(seq1) // 2

        _, _, scoreL = twoColumnSearch(alphabet, scoringMatrix, seq1[:midX+1], seq2)
        _, _, scoreR = twoColumnSearch(alphabet, scoringMatrix, seq1[midX+1:], seq2)

        sumOfScores = [x+y for x, y in zip(scoreL, scoreR[::-1])]
        midY = sumOfScores.index(max(sumOfScores))

        left = hirschberg(alphabet, scoringMatrix, seq1[:midX+1], seq2[:midY+1])
        right = hirschberg(alphabet, scoringMatrix, seq1[midX+1:], seq2[midY+1:])
        return #todo

a = dynproglin ("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
print("Score:   ", a[0])
print("Indices: ", a[1],a[2])

## Expected output:
# Score:    5
# Indices:  [3,5,6][1,2,3]
