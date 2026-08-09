# External-map comparison

`compare_with_kletenik_anosognosia.py` compares the public historical combined map with NeuroVault image 795501. It always reports the observed spatial Pearson correlation after aligning the external map to the historical MNI grid.

With authorized participant-level inputs, the script can additionally generate a revised two-sided Freedman–Lane null distribution. Each permutation recomputes the two partial-correlation maps, applies the historical nominal 101/181 Fisher-z weighting, and correlates the permuted combined map with the fixed external map.

This revised null is explicit and reproducible, but it should not be presented as the unidentified historical procedure that produced the manuscript-reported `p = 0.93`. The manuscript value should be updated to the documented result after execution, or the comparison should be described without that historical p-value.
