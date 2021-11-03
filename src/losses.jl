function calculate_loss_robust(x; α=2, c=1)
    if α ≈ 2
        loss = 0.5 * (x / 2)^2
    elseif α ≈ 0
        loss = log((x / c)^2 + 1)
#     elseif α > 1e9
#         loss = 1 - exp(-0.5 * (x / c)^2)
    else
        A = abs(α - 2) / α
        B = ((x / c)^2) / abs(α - 2) + 1
        loss = A * (B^(α / 2) - 1)
    end
end


function RMSE(x, y, w = ones(length(x)))
    @assert length(x) == length(y)
    se = w .* (x .- y) .^2
    return √(sum(se)/sum(w))
end
