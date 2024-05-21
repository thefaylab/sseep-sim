species <- list(
  summerflounder = list(
      ages = c(0:7),
      years = c(1:15),
      Rec_age0 = c(28416000,33088000,44582000,60598000,48689000),
      F = 0.422,
      M = 0.25,
      Linf = 83.6,
      k = 0.14
    ),
  scup = list(
    ages = c(0:7),
    years = c(1:15),
    Rec_age0 = c(107000,142000,75000,61000,112000),
    F = 0.1,
    M = 0.2,
    Linf = 46.6,
    k = 0.15
  ),
  blackseabass = list(
    ages = c(0:8),
    years = c(1:15),
    Rec_age0 = c(39629000,93799000,51186000,14872000,46198000),
    F = 0.4,
    M = 0.4,
    Linf = 58.9,
    k = 0.22
  ),
  atlanticmackerel = list(
    ages = c(0:10),
    years = c(1:15),
    Rec_age0 = c(83000,38000,91000,163000,455000),
    F = 0.9,
    M = 0.2,
    Linf = 39.18,
    k = 0.387
  ),
  atlanticherring = list(
    ages = c(0:10),
    years = c(1:15),
    Rec_age0 = c(1370270,1608170,776348,174758,392286),
    F = 0.51,
    M = 0.35,
    Linf = 29.92,
    k = 0.46
  ),
  butterfish = list(
    ages = c(0:4),
    years = c(1:15),
    Rec_age0 = c(8301000,10299000,4489000,7006000,9813000),
    F = 0.2,
    M = 1.3,
    Linf = 18.3,
    k = 0.8
  ))

saveRDS(species,  file = here("data", "species_data_052024.rds"))

