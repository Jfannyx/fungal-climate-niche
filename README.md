# Niche modeling predicts that soil fungi occupy a precarious climate in boreal forests
-thresholds_permuted <- findThresholds(neon_dob_permuted, lognets_permuted, pred_vars = terms)
These codes don't work for me with error being 'unused argument (pred_vars = terms)'
when I changed "pred_vars" to "terms", it works, but the following function getNicheMetrics2D doesn't work.
