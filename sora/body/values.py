not_found_message = 'specified object was not found'
many_objects_message = 'specified query matched more than one object'

spectral_features = {
    "id1": "Very steep red slope shortward of 0.75 μm; moderately deep absorption feature longward of 0.75 μm.",
    "id2": ("Linear, generally featureless spectra. Differences in UV absorption features and presence/absence"
            " of narrow absorption feature near 0.7 μm."),
    "id3": ("Linear, generally featureless spectra. Differences in UV absorption features and presence/absence"
            " of narrow absorption feature near 0.7 μm."),
    "id4": "Relatively featureless spectrum with very steep red slope.",
    "id5": ("Generally featureless spectrum with reddish slope; differences in subtle absorption features"
            " and/or spectral curvature and/or peak relative reflectance."),
    "id6": "Reddish slope shortward of 0.7 μm; deep, rounded absorption feature longward of 0.75 μm.",
    "id7": "Moderate reddish slope downward of 0.7 μm; deep absorption longward of 0.75 μm.",
    "id8": ("Moderately steep reddish slope downward of 0.7 μm; moderate to steep absorption longward of 0.75 μm;"
            " peak of reflectance at 0.73 μm. Bus subgroups intermediate between S and A, K, L, Q, R classes."),
    "id9": "Moderately reddish shortward of 0.75 μm; flat afterward.",
    "id10": "Reddish shortward of 0.7 μm; extremely deep absorption longward of 0.75 μm.",
    "id11": ("Moderately steep red slope shortward of 0.75 μm; smoothly angled maximum and flat to blueish"
             " longward of 0.75 μm, with little or no curvature."),
    "id12": "Very steep red slope shortward of 0.75 μm; flat longward of 0.75 μm; differences in peak level.",
    "id13": "Peculiar trend, known so far for very few asteroids."
}

smass = {
    "A": spectral_features["id1"],
    "B": spectral_features["id2"],
    "C": spectral_features["id3"],
    "Cb": spectral_features["id3"],
    "Ch": spectral_features["id3"],
    "Cg": spectral_features["id3"],
    "Chg": spectral_features["id3"],
    "D": spectral_features["id4"],
    "X": spectral_features["id5"],
    "Xc": spectral_features["id5"],
    "Xe": spectral_features["id5"],
    "Xk": spectral_features["id5"],
    "Q": spectral_features["id6"],
    "R": spectral_features["id7"],
    "S": spectral_features["id8"],
    "Sa": spectral_features["id8"],
    "Sk": spectral_features["id8"],
    "Sl": spectral_features["id8"],
    "Sq": spectral_features["id8"],
    "Sr": spectral_features["id8"],
    "T": spectral_features["id9"],
    "V": spectral_features["id10"],
    "K": spectral_features["id11"],
    "L": spectral_features["id12"],
    "Ld": spectral_features["id12"],
    "O": spectral_features["id13"]
}

tholen = {
    "A": spectral_features["id1"],
    "B": spectral_features["id2"],
    "F": spectral_features["id2"],
    "C": spectral_features["id3"],
    "G": spectral_features["id3"],
    "D": spectral_features["id4"],
    "E": spectral_features["id5"],
    "M": spectral_features["id5"],
    "P": spectral_features["id5"],
    "Q": spectral_features["id6"],
    "R": spectral_features["id7"],
    "S": spectral_features["id8"],
    "T": spectral_features["id9"],
    "V": spectral_features["id10"]
}
