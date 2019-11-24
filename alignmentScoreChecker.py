def alignmentScoreChecker(sequences, alignment):
    score = 0
    print(sequences, alignment)
    for i, j in zip(alignment[0], alignment[1]):
        score += getScoreOfMatchingCharactersFromScoringMatrix(
            sequences[0][i], sequences[1][j])

    minBound, maxBound = min(alignment[0]), max(alignment[0])
    for i in [x for x in range(minBound, maxBound) if x not in alignment[0]]:
        score += getScoreOfMatchingCharactersFromScoringMatrix(
            sequences[0][i], '-')

    minBound, maxBound = min(alignment[1]), max(alignment[1])
    for i in [x for x in range(minBound, maxBound) if x not in alignment[1]]:
        score += getScoreOfMatchingCharactersFromScoringMatrix(
            '-', sequences[1][i])

    return score


def getScoreOfMatchingCharactersFromScoringMatrix(character1, character2):
    global alphabet, scoringMatrix
    rowIndex = alphabet.index(character2)
    colIndex = alphabet.index(character1)

    row = scoringMatrix[rowIndex]
    out = row[colIndex]

    return out


alphabet = "ABC-"
scoringMatrix = [[1, -1, -2, -1], [-1, 2, -4, -1],
                 [-2, -4, 3, -2], [-1, -1, -2, 0]]

print(alignmentScoreChecker(("AABBAACA", "CBACCCBA"), ([3, 5, 6], [1, 2, 3])))

alphabet = "ABCD-"
scoringMatrix = [
    [1, -5, -5, -5, -1],
    [-5, 1, -5, -5, -1],
    [-5, -5, 5, -5, -4],
    [-5, -5, -5, 6, -4],
    [-1, -1, -4, -4, -9]
]

print(alignmentScoreChecker(("AAAAACCDDCCDDAAAAACC", "CCAAADDAAAACCAAADDCCAAAA"),
                            ([5, 6, 7, 8, 9, 10, 11, 12, 18, 19], [0, 1, 5, 6, 11, 12, 16, 17, 18, 19])))

print(alignmentScoreChecker(("AACAAADAAAACAADAADAAA", "CDCDDD"),
                            ([2, 6, 11, 14, 17], [0, 1, 2, 3, 4])))

print(alignmentScoreChecker(("DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD",
                             "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"), ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 61, 62, 63, 64, 65, 66, 67, 68, 69])))
