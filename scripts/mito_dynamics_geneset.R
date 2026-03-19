# ============================================================
# Mito Dynamics Gene Set
# Three rings + redox bridge connecting to neufer_geneset.R
# ============================================================

mito_dynamics_geneset <- list(

  # RING 1: Lactate-mitochondria axis
  mct_transport        = c("Slc16a1", "Slc16a3", "Slc16a2", "Slc16a7"),
  ldh                  = c("Ldha", "Ldhb"),
  mito_pyruvate_carrier = c("Mpc1", "Mpc2"),
  pdh                  = c("Pdha1", "Pdhb", "Dlat", "Dld", "Pdhx"),
  pdk_pdp              = c("Pdk1", "Pdk2", "Pdk3", "Pdk4", "Pdp1", "Pdp2"),
  lactate_signaling    = c("Hcar1"),

  # RING 2: Mitochondrial dynamics
  fission   = c("Dnm1l", "Fis1", "Mff", "Mief1", "Mief2"),
  fusion    = c("Mfn1", "Mfn2", "Opa1"),
  mitophagy = c("Pink1", "Prkn", "Bnip3", "Bnip3l", "Fundc1",
                "Ulk1", "Ulk2", "Atg5", "Atg7", "Sqstm1"),
  biogenesis = c("Ppargc1a", "Ppargc1b", "Tfam", "Nrf1",
                 "Gabpa", "Esrra", "Sirt1", "Sirt3"),

  # RING 3: Mitochondrial signaling / retrograde communication
  mito_upr           = c("Lonp1", "Clpp", "Hspd1", "Hspe1", "Hspa9"),
  integrated_stress  = c("Atf4", "Atf5", "Ddit3", "Eif2ak4"),
  ampk_axis          = c("Prkaa1", "Prkaa2", "Stk11", "Camkk2"),
  calcium            = c("Mcu", "Micu1", "Micu2", "Vdac1", "Vdac2", "Vdac3"),
  mito_derived_signals = c("Fgf21", "Gdf15"),

  # BRIDGE: Links back to redox / Neufer analysis
  ros_signaling = c("Nfe2l2", "Keap1", "Hmox1", "Nqo1", "Txnip")
)

all_mito_dynamics_genes <- unique(unlist(mito_dynamics_geneset))
cat("Mito dynamics gene set:", length(all_mito_dynamics_genes), "unique genes\n")
