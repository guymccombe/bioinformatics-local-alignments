def findAlignments(sequence1, sequence2):
    if sequence1 > sequence2:
        findAlignments(sequence2, sequence1)
    else:
        column = [0 for i in range(len(sequence2) + 1)]
        bestScore = 0
        coordinateOfBest = (0,0)

        for i in range(len(sequence1)):
            print (column)
            lastCol = column
            column = [0]
            for j in range(1, len(sequence2) + 1):
                score = max([
                    lastCol[j-1] + getScoreOfMatchingCharactersFromScoringMatrix(sequence1[i], sequence2[j-1]),
                    column[j-1]  + getScoreOfMatchingCharactersFromScoringMatrix(sequence1[i], '-'),
                    lastCol[j]   + getScoreOfMatchingCharactersFromScoringMatrix('-', sequence2[j-1]),
                    0
                ])
                column.append(score)
                if score > bestScore:
                    bestScore = score
                    coordinateOfBest = (i, j-1)

        print(bestScore, coordinateOfBest)
        #todo backtrack


def getScoreOfMatchingCharactersFromScoringMatrix(character1, character2):
    global alphabet, scoringMatrix

    rowIndex = alphabet.index(character2)
    colIndex = alphabet.index(character1)
    
    row = scoringMatrix[rowIndex]
    out = row[colIndex]

    return out

alphabet = ['A', 'C', 'G', 'T']
scoringMatrix = [[ 0,-2,-2,-2,-2],
                 [-2, 1,-1,-1,-1],
                 [-2,-1, 1,-1,-1],
                 [-2,-1,-1, 1,-1],
                 [-2,-1,-1,-1, 1]]
alphabet.insert(0, '-')

findAlignments('TAATA','TACTAA')