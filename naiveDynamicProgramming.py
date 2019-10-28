def findAlignments(sequence1, sequence2):
    alignmentScoreMatrix, alignmentBacktrackMatrix = initialiseMatrices(sequence1, sequence2)
    score = 0
    coordinateOfBest = (0,0)

    for i in range (1, len(alignmentScoreMatrix)):
        for j in range (1, len(alignmentScoreMatrix[0])):
            alignmentScoreMatrix[i][j], alignmentBacktrackMatrix[i][j] = max([
                (alignmentScoreMatrix[i-1][j-1] + getScoreOfMatchingCharactersFromScoringMatrix(sequence1[i-1], sequence2[j-1]), 'D'),
                (alignmentScoreMatrix[i-1][ j ] + getScoreOfMatchingCharactersFromScoringMatrix(sequence1[i-1], '-'), 'U'),
                (alignmentScoreMatrix[ i ][j-1] + getScoreOfMatchingCharactersFromScoringMatrix('-', sequence2[j-1]), 'L'),
                (0, 'G')
            ])

            if alignmentScoreMatrix[i][j] > score:
                score = alignmentScoreMatrix[i][j]
                coordinateOfBest = (i,j)

    sequence1Indices, sequence2Indices = backtrackAndReturnIndicesOfAlignment(alignmentBacktrackMatrix, coordinateOfBest)

    print(alignmentScoreMatrix)
    print(alignmentBacktrackMatrix)
    print(score, sequence1Indices, sequence2Indices)

def initialiseMatrices(sequence1, sequence2):
    scoringMatrix = [[[] for i in range(len(sequence2) + 1)] for i in range(len(sequence1) + 1)]
    backtrackMatrix = [[[] for i in range(len(sequence2) + 1)] for i in range(len(sequence1) + 1)]

    scoringMatrix[0][0], backtrackMatrix[0][0] = 0, 'G'

    for i in range(1, len (scoringMatrix)):
        scoringMatrix[i][0], backtrackMatrix[i][0] = max(
            (scoringMatrix[i-1][0] + getScoreOfMatchingCharactersFromScoringMatrix('-', sequence1[i-1]), 'U'),
            (0, 'G')
        )

    for i in range(1, len (scoringMatrix[0])):
        scoringMatrix[0][i], backtrackMatrix[0][i] = max (
            (scoringMatrix[0][i-1] + getScoreOfMatchingCharactersFromScoringMatrix('-', sequence2[i-1]), 'L'),
            (0, 'G')
        )

    return scoringMatrix, backtrackMatrix

def getScoreOfMatchingCharactersFromScoringMatrix(character1, character2):
    global alphabet, scoringMatrix

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
        indices1.insert(0, i-1)
        indices2.insert(0, j-1)

        if currentState == 'D':
            i -= 1
            j -= 1
        elif currentState == 'U':
            i -= 1
        else:
            j -= 1
        currentState = backtrackMatrix[i][j]

    return indices1, indices2









alphabet = ['A', 'C', 'G', 'T']
scoringMatrix = [[ 0,-2,-2,-2,-2],
                 [-2, 1,-1,-1,-1],
                 [-2,-1, 1,-1,-1],
                 [-2,-1,-1, 1,-1],
                 [-2,-1,-1,-1, 1]]
alphabet.insert(0, '-')

findAlignments('TAATA','TACTAA')