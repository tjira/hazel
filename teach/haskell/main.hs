-- define some basic types
data Matrix = Matrix [[Double]] deriving(Eq, Show)
data Vector = Vector [Double] deriving(Eq, Show)

--- basic matrix operations
(Matrix a1) <<+>> (Matrix a2) = Matrix $ zipWith (zipWith (+)) a1 a2
(Matrix a1) <<->> (Matrix a2) = Matrix $ zipWith (zipWith (-)) a1 a2
(Matrix a1) <<*>> (Matrix a2) = Matrix $ zipWith (zipWith (*)) a1 a2

--- advanced matrix operations

--- basic vector operations
(Vector a1) <+> (Vector a2) = Vector $ zipWith (+) a1 a2
(Vector a1) <-> (Vector a2) = Vector $ zipWith (-) a1 a2
(Vector a1) <*> (Vector a2) = Vector $ zipWith (*) a1 a2

--- vector dot product
(Vector a1) `dot` (Vector a2) = sum $ zipWith (*) a1 a2

main = do
    let s = Vector [1, 1]
    let t = Vector [2, 1]
    let a = Matrix [[1, 2], [3, 4]]
    let b = Matrix [[2, 1], [3, 5]]
    -- print $ s `dot` t
    print $ a
