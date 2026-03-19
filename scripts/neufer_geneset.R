# ============================================================
# Neufer Mitochondrial Redox Framework — Curated Gene Set
# Based on: Vandiver & Neufer (2025) "Mitochondria: connecting
# oxygen to life" from On Oxygen (Elsevier)
# ============================================================

neufer_geneset <- list(

  # ETS Core
  ets_CI = c("Ndufs1","Ndufs2","Ndufs3","Ndufs7","Ndufs8",
             "Ndufv1","Ndufv2","Ndufa9","Ndufa10","Ndufb8"),
  ets_CIII = c("Uqcrc1","Uqcrc2","Uqcrfs1","Uqcrb","Cyc1"),
  ets_CIV = c("Cox4i1","Cox5a","Cox5b","Cox6a1","Cox6b1","Cox7a2"),
  cyt_c = c("Cycs"),

  # Q-Pool Feeder Enzymes (reductive stress entry points)
  q_pool = c("Sdha","Sdhb","Sdhc","Sdhd",
             "Etfdh","Etfa","Etfb",
             "Gpd2","Dhodh","Chdh",
             "Prodh","Prodh2","Sqor"),

  # ATP Synthase
  atp_synthase = c("Atp5f1a","Atp5f1b","Atp5f1c","Atp5f1d","Atp5f1e",
                   "Atp5pb","Atp5mc1","Atp5mc2","Atp5mc3","Atp5po"),

  # NNT — central to Neufer thesis
  nnt = c("Nnt"),

  # Glutathione buffering circuit
  gsh = c("Gpx1","Gpx4","Gsr","Gclc","Gclm","Gss","Glrx2"),

  # Thioredoxin buffering circuit
  trx = c("Txn2","Txnrd2","Prdx3","Prdx5"),

  # Superoxide handling
  sod = c("Sod2","Sod1"),

  # Non-NNT NADPH sources (for comparison)
  nadph_other = c("Idh2","Me2","Aldh1l2"),

  # Uncoupling proteins
  ucp = c("Ucp1","Ucp2","Ucp3"),

  # Shuttle systems
  shuttles = c("Mdh1","Mdh2","Got1","Got2","Ldha","Ldhb","Gpd1","Gpd2"),

  # PDH complex
  pdh = c("Pdha1","Pdhb","Dlat","Dld","Pdhx"),

  # Beta-oxidation (reductive stress generators)
  beta_ox = c("Acadm","Acadl","Acadvl","Hadha","Hadhb","Cpt1a","Cpt1b","Cpt2")
)

all_neufer_genes <- unique(unlist(neufer_geneset))
cat("Neufer framework gene set:", length(all_neufer_genes), "unique genes\n")
