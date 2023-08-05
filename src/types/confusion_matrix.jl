struct ConfusionMatrix
    tp::Any
    tn::Any
    fp::Any
    fn::Any
    function ConfusionMatrix(a, b, c, d)
        s = a + b + c + d
        return new(a / s, b / s, c / s, d / s)
    end
end

ConfusionMatrix(A::Matrix{Float64}) = ConfusionMatrix(A[1, 1], A[2, 2], A[2, 1], A[1, 2])
ConfusionMatrix(A::Matrix) = ConfusionMatrix(A[1, 1], A[2, 2], A[2, 1], A[1, 2])
