def b10film_ion_escape_probability( x ):
  # X is the distance between the point where the neutron is absorbed in a B4C
  # thin film, and the counting gas volume. Returns the probability of secondary
  # ions reaching the counting gas and depositing at least 120keV of energy.
  #
  # Parameterised based on Geant4 simulations by T. Kittelmann, October 2025.
  #
  # Details (results are not believed to be overly sensitive to these):
  #
  #   Geant4 11.3.2, NCrystal materials, QGSP_BIC_HP_EMZ list, Monochromatic
  #   neutrons hitting spherical geometry. B4C enrichment level 95%. 70-30
  #   Ar-CO2 counting gas.
  #
  if x < 1.01:
    return 0.984 - 0.6135644 * x
  if x < 2.5:
    if x < 1.3:
      if x < 1.1:
        return 0.7873778 - 0.4188889 * x
      return 0.6159 - 0.263 * x
    return 0.4934833 - 0.1688333 * x
  if x > 4.2199578:
    return 0.0
  if x < 2.7:
    return 0.40265 - 0.1325 * x
  if x < 2.9:
    return 0.30005 - 0.0945 * x
  if x < 3.1:
    return 0.2 - 0.06 * x
  if x < 3.435:
    return 0.1237493 - 0.03540299 * x
  if x < 3.85:
    return 0.01488675 - 0.003710843 * x
  return 0.006843243 - 0.001621622 * x
