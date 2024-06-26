# In this code:
# Prepping cover data from the initial IDE dataset
# querying the AusTraits database
# calculating DSI for sites

# load libraries
library(Taxonstand)
library(tidyverse)
library(stringr)

install.packages("remotes")
remotes::install_github("traitecoevo/austraits", dependencies = TRUE, upgrade = "ask")
library(austraits) 

#### Prepping final species list ####

# read in raw IDE cover data
data=read.csv("./Raw.Data/IDE_cover_2023-01-02.csv")

# get taxon list
taxon.df=as.data.frame(table(data$Taxon))
# 1904 species, but some are sp. so will need to be excluded
# 35977 rows

# excluded species with sp.

data.2=subset(data, !(data$Taxon %in% c("ADESMIA SP.(eea.br)","AGOSERIS SP.(OREAC.US)","AGROSTIS sp.","AGROSTIS SP.",
                                "AGROSTIS SP.(SCRUZH.US)","ALLIUM SP.", "ALLIUM SP.(BFL.US)","ALOPECURUS  SP.(LCSOUTH.CL)",
                                "ALOPECURUS SP.(LCNORTH.CL)","ALOPECURUS SP.(PURDUE.US)","AMBROSIA SP.(PURDUE.US)",
                                "ANTHOXANTHUM  SP.(LCSOUTH.CL)","ARISTIDA SP.","ASTER SP.(BFL.US)","ASTER SP.1(MILPARINKA.AU)",
                                "ASTER SP.7(QUILPIE.AU)","ASTRAGALUS  SP.(LCSOUTH.CL)","ASTRAGALUS SP.(LCNORTH.CL)",
                                "ATRIPLEX SP","ATRIPLEX SP.(MILPARINKA.AU)","AUSTROSTIPA SP.(JILPANGER.AU)","AVENA SP.",
                                "AXONOPUS SP.(BIVENSARM.US)","BACCHARIS SP.(chacra.ar)","BIDENS Bidens","BOTHRIOCHLOA SP.(NAPOSTA.AR)",
                                "BRACHYCOME SP.","BROMUS SP.","BRYOPHYTE","BRYOPHYTE SP.(CHACRA.AR)","BRYOPHYTE SP.(KERNB.CA)",
                                "BRYOPHYTE SP.(KERNNU.CA)","BRYOPHYTE SP.(LYGRAINT.NO)","BRYOPHYTE SP.(lygraold.no)",
                                "BRYOPHYTE SP.1(LYGRAINT.NO)","BRYOPHYTE SP.1(lygraold.no)","BRYOPHYTE SP.2(LYGRAINT.NO)",
                                "BRYOPHYTE SP.2(lygraold.no)","BRYOPHYTE SP.3(lygraold.no)","BULBOSTYLIS SP.(UKULINGADRT.ZA)",
                                "CAREX SP.","CAREX SP.(BIDDULPH.CA)","CAREX SP.(BIVENSARM.US)","CAREX SP.(KERNB.CA)",
                                "CAREX SP.(KERNNU.CA)","CAREX SP.(KINSELLA.CA)","CAREX SP.(LYGRAINT.NO)","CAREX SP.(lygraold.no)",
                                "CAREX SP.(MATADOR.CA)","CAREX SP.(MATTHEIS.CA)","CAREX SP.(OREAC.US)","CAREX SP.1(LYGRAINT.NO)",
                                "CASTILLEJA SP.(OKLAH.US)","CENTAUREA sp.","CENTAUREA SP.","CENTAUREA SP.(ayora.es)",
                                "CERASTIUM  SP.(LCSOUTH.CL)","CERASTIUM SP.(LCNORTH.CL)","CLADONIA SP.(LYGRAINT.NO)",
                                "CLADONIA SP.1(LYGRAINT.NO)","CONVOLVULUS SP.(JILPANGER.AU)","CRASSULA SP.(JILPANGER.AU)",
                                "CREPIS SP.(OREAA.US)","CREPIS SP.(OREAC.US)","CRYPTANTHA SP.(JRNCHI.US)","CYPERUS SP.",
                                "CYPERUS SP.(BIVENSARM.US)","DANTHONIA SP.(eea.br)","DELPHINIUM SP.(OREAC.US)","DESCURAINIA SP.(SPVDRT.AR)",
                                "DIANTHUS SP.(qdtsouth.cl)","DICHANTHELIUM SP.(SLP.US)","DICHANTHELIUM SP.2(SLP.US)","DICRANUM SP.(LYGRAINT.NO)",
                                "DICRANUM SP.(lygraold.no)","DIOSCOREA  S2P.(LCSOUTH.CL)","DIOSCOREA  SP1.(LCSOUTH.CL)","DIOSCOREA SP1.(LCNORTH.CL)",
                                "DIOSCOREA SP1.(QDTNORTH.CL)","DIOSCOREA SP1.(qdtsouth.cl)","DIOSCOREA SP2.(qdtsouth.cl)","ELEOCHARIS SP.(OKLAH.US)",
                                "ERAGROSTIS SP.(OKLAH.US)","ERIGERON SP.(KERNB.CA)","ERIGERON SP.(KERNNU.CA)","ERIGERON SP.(MATADOR.CA)",
                                "ERIOGONUM SP.(NNSS.US)","EUONYMUS SP.","FORB sp.","FORB SP.(KERNB.CA)","FORB SP.(KERNNU.CA)","GAILLARDIA SP.(bfl.us)",
                                "GAMOCHAETA SP.(CHACRA.AR)","GAMOCHAETA SP.(NAPOSTA.AR)","GAMOCHAETA SP.(SPVDRT.AR)","GERANIUM  SP.(LCSOUTH.CL)",
                                "GERANIUM SP.(JILPANGER.AU)","GERANIUM SP.(LCNORTH.CL)","GERANIUM SP.(MILPARINKA.AU)","GILIA SP.(NNSS.US)",
                                "GILIA SP.(OREAA.US)","GILIA SP.(OREAC.US)","GUTIERREZIA SP.(BFL.US)","HELIANTHEMUM  SP.(LCSOUTH.CL)",
                                "HELIANTHEMUM SP.(LCNORTH.CL)","HELIANTHUS SP.(BROOKDALE.CA)","HIBBERTIA SP.(JILPANGER.AU)","HIERACIUM SP.(ELVADRT.EE)",
                                "HYMENOPAPPUS SP. (cdpt_drt.us)","HYPOCHAERIS SP.","HYPOCHAERIS SP.(SCRUZH.US)","HYPOCHAERIS SP.(SCRUZM.US)",
                                "IPOMOEA SP.(BIVENSARM.US)","JUNCUS SP.","JUNCUS SP.(BIVENSARM.US)","JUNCUS SP.(LYGRAINT.NO)","JUNCUS SP.(SLP.US)",
                                "JUNELLIA SP.(CHACRA.AR)","LATHYRUS SP.(BROOKDALE.CA)","LEONTODON SP.","LEPIDIUM SP.(NNSS.US)","LESQUERELLA SP.",
                                "LEUCOCORYNE SP1.(QDTNORTH.CL)","LEUCOCORYNE SP2.(QDTNORTH.CL)","LICHEN ","LICHEN","LICHEN SP.(LYGRAINT.NO)","LICHEN SP.(lygraold.no)",
                                "LICHEN SP.1(lygraold.no)","LINUM SP.(BFL.US)","LITHOPHRAGMA SP.(OREAA.US)","LITHOPHRAGMA SP.(OREAC.US)","LOLIUM SP.(CHACRA.AR)",
                                "LOMATIUM SP.(OREAA.US)","LOMATIUM SP.(OREAC.US)","LOTUS","LOTUS ","LOTUS SP.","LUZULA SP.","LUZULA SP.(FALLS.AU)","MADIA SP.(OREAC.US)",
                                "MALVA SP.","MALVA SP.(SCRUZL.US)","MELILOTUS SP.(ESW.CA)","MINURIA SP.(credoj.au)","MINURIA SP.(CREDOM.AU)","MONTIA SP.(OREAA.US)",
                                "MONTIA SP.(OREAC.US)","OLEARIA SP.","OPHIOGLOSSUM SP.(JILPANGER.AU)","OXALIS  SP.(LCSOUTH.CL)","OXALIS SP.(LCNORTH.CL)",
                                "OXALIS SP.(QDTNORTH.CL)","OXALIS SP.(qdtsouth.cl)","PAEPALANTHUS SP.GUARIBAS.BR","PANICUM SP.(SLP.US)","PASPALUM SP.",
                                "PASPALUM SP.(SLP.US)","PECTOCARYA SP.","PLANTAGO SP.","PLANTAGO SP.(OKLAH.US)","PLANTAGO SP.(SLP.US)",
                                "POA SP.(LYGRAINT.NO)","POLYGALA SP.","POLYGALA SP.GUARIBAS.BR","POLYGONUM SP.(OREAA.US)","POLYGONUM SP.(OREAC.US)",
                                "POLYPOGON SP.","POTENTILLA SP.(KERNNU.CA)","PRUNELLA SP.","PRUNUS SP.","PTILOTUS SP.","PTILOTUS SP.(MILPARINKA.AU)",
                                "RANUNCULUS sp.","RANUNCULUS SP.","ROSA SP.","RUBUS SP.(BIVENSARM.US)","RUBUS SP.(OKLAH.US)","RUMEX SP.",
                                "RYTIDOSPERMA SP.(JILPANGER.AU)","SCLEROLAENA SP.","SCLEROLAENA SP.1(MILPARINKA.AU)","SCLEROLAENA SP.2(MILPARINKA.AU)",
                                "SCUTELLARIA SP.(SLP.US)","SENECIO SP.(OREAC.US)","SIDA SP.","SIDA SP.(CREDOJ.AU)","SIDA SP.(CREDOM.AU)",
                                "SILENE SP.(KERNNU.CA)","SMILAX SP.(BIVENSARM.US)","SOLIDAGO SP.(NAPOSTA.AR)","SOLIDAGO SP.(OKLAH.US)",
                                "SONCHUS  SP1.(LCSOUTH.CL)","SONCHUS  SP2.(LCSOUTH.CL)","SPHAERALCEA SP.","SPOROBOLUS SP.","Sporobulus sp.",
                                "STELLARIA  SP.(LCSOUTH.CL)","STELLARIA SP.(LCNORTH.CL)","TALINUM SP.(CERRILLOS.AR)","TARAXACUM SP.",
                                "TARAXACUM SP.(ELVADRT.EE)","THELYMITRA SP.(JILPANGER.AU)","TOXICOSCORDION SP.","TRIFOLIUM SP.","TRIFOLIUM SP.(BROOKDALE.CA)",
                                "TRIFOLIUM sp.(jilpanger.au)","TRIFOLIUM SP.(PURDUE.US)","TRIFOLIUM SP.(SCRUZH.US)","ULMUS SP.(OKLAH.US)","UNKNOWN ",
                                "UNKNOWN  SP.(LCSOUTH.CL)","UNKNOWN  SP2.(LCSOUTH.CL)","UNKNOWN A(COWIDRT.CA)","UNKNOWN AMARYLLIDACEAE SP.(HARD.US)",
                                "UNKNOWN ASTERACEAE ","UNKNOWN ASTERACEAE  SP4.(LCSOUTH.CL)","UNKNOWN ASTERACEAE SP.(CREDOM.AU)","UNKNOWN ASTERACEAE SP.(HARD.US)",
                                "UNKNOWN ASTERACEAE SP.(OREAC.US)","UNKNOWN ASTERACEAE SP.2(MILPARINKA.AU)","UNKNOWN ASTERACEAE SP.2(OREAC.US)","UNKNOWN ASTERACEAE SP1.(qdtsouth.cl)",
                                "UNKNOWN ASTERACEAE SP2.(qdtsouth.cl)","UNKNOWN ASTERACEAE SP3.(LCNORTH.CL)","UNKNOWN ASTERACEAE SP3.(QDTNORTH.CL)",
                                "UNKNOWN ASTERACEAE SP3.(qdtsouth.cl)","UNKNOWN ASTERACEAE SP4.(LCNORTH.CL)","UNKNOWN BRASSICACEAE SP.(HARD.US)",
                                "UNKNOWN D(COWIDRT.CA)","UNKNOWN FABACEAE SP.(HARD.US)","UNKNOWN FORB(BIDDULPH.CA)","UNKNOWN G(COWIDRT.CA)","UNKNOWN GRASS",
                                "UNKNOWN GRASS ","UNKNOWN GRASS SP.","UNKNOWN GRASS(COWIDRT.CA)","UNKNOWN H(COWIDRT.CA)","UNKNOWN MINT(COWIDRT.CA)",
                                "UNKNOWN ONAGRACEAE SP.(QDTNORTH.CL)","UNKNOWN POACEAE SP.(HARD.US)","UNKNOWN POACEAE SP.(JILPANGER.AU)","UNKNOWN POACEAE SP.(MATADOR.CA)",
                                "UNKNOWN POACEAE SP.(MILPARINKA.AU)","UNKNOWN POACEAE SP.(QUILPIE.AU)","UNKNOWN POACEAE SP.(SPVDRT.AR)",
                                "UNKNOWN POACEAE SP.1(QUILPIE.AU)","UNKNOWN POACEAE SP1.(qdtsouth.cl)","UNKNOWN POACEAE SP2.(QDTNORTH.CL)",
                                "UNKNOWN POACEAE SP3.(qdtsouth.cl)","UNKNOWN POLEMONIACEAE SP.(OREAA.US)","UNKNOWN POLEMONIACEAE SP.(OREAC.US)",
                                "UNKNOWN SCLEROLAENA ","UNKNOWN SP.","UNKNOWN SP.(CERRILLOS.AR)","UNKNOWN SP.(GUARIBAS.BR)","UNKNOWN SP.(HARD.US)",
                                "UNKNOWN sp.(jilpanger.au)","UNKNOWN SP.(JRNCHI.US)","UNKNOWN SP.(KONZADRT.US)","UNKNOWN SP.(LCNORTH.CL)",
                                "UNKNOWN SP.(LYGRAINT.NO)","UNKNOWN SP.(MATADOR.CA)","UNKNOWN SP.(NNSS.US)","UNKNOWN SP.(OKLAH.US)","SELAGINELLA DENSA",
                                "UNKNOWN SP.(OREAA.US)","UNKNOWN SP.(OREAC.US)","UNKNOWN SP.(SCRUZH.US)","UNKNOWN SP.(SLP.US)","UNKNOWN SP.1(OKLAH.US)",
                                "UNKNOWN SP.11(OKLAH.US)","UNKNOWN SP.2(MILPARINKA.AU)","UNKNOWN SP.2(OKLAH.US)","UNKNOWN SP.4(OKLAH.US)","UNKNOWN SP.4(SPVDRT.AR)",
                                "UNKNOWN SP.5(SPVDRT.AR)","UNKNOWN SP.6(OKLAH.US)","UNKNOWN SP.6(SPVDRT.AR)","UNKNOWN SP.7(OKLAH.US)","UNKNOWN SP.9(OKLAH.US)",
                                "UNKNOWN SP.9(SPVDRT.AR)","UNKNOWN SP.D20(OKLAH.US)","UNKNOWN SP.D22(OKLAH.US)","UNKNOWN SP.D25(OKLAH.US)","UNKNOWN SP.D29(OKLAH.US)",
                                "UNKNOWN SP3.(LCNORTH.CL)","UNKNOWN VIOLACEAE  SP.(LCSOUTH.CL)","UNKNOWN WEED(COWIDRT.CA)","UTRICULARIA SP.(GUARIBAS.BR)",
                                "VERNONIA SP.(PURDUE.US)","VERONICA SP.(OREAA.US)","VICIA SP.","VIOLA sp.","VIOLA SP.","VIOLA SP.(ayora.es)","VIOLA SP.(OREAC.US)",
                                "WAHLENBERGIA SP.","XYRIS SP.","XYRIS SP.GUARIBAS.BR","ACALYPHA SP.(SLP.US)","ATRIPLEX SP.","LASERPITIUM SP.")))

# correcting some errors
which(data.2$Taxon == "CROTON potsii")
data.2[19480,13] = "CROTON POTTSII"
which(data.2$Taxon == "GOODENIA CYLCOPTERA")
data.2[c(7936,7950,7953,11914,11934,12671,12690,12708,12715,29342,29345,29358,29376),13] = "GOODENIA CYCLOPTERA"
which(data.2$Taxon == "LYTHRUM HYSSOPIFOLIUM")
data.2[c(23264,23320,23616,23627),13] = "LYTHRUM HYSSOPIFOLIA"

# get taxon list
taxon.2=as.data.frame(table(data.2$Taxon))
# 1600 species
# 33970 rows

# write.csv(taxon.2, "IDE.species.check.csv")

#### verify species names ####
# https://cran.r-project.org/web/packages/Taxonstand/Taxonstand.pdf
taxon.check=Taxonstand::TPL(taxon.2$Var1)

# write.csv(taxon.check, file="taxon.check.csv")

#### calculate cover change ####
# comparing control-control, drought-drought
# BACI design of (drought.after-drought.before)-(control.after-control.before)

before = data.2 %>%
  filter(n_treat_years == 0)
# plots with before data

before.1 = before %>%
  group_by(site_code, Taxon, year, block, plot,trt) %>%
  summarize(max_cover = max_cover)

before.2 = pivot_wider(before.1, names_from = trt, values_from = max_cover)

before.3 = before.2 %>%
  group_by(site_code, Taxon) %>%
  reframe(mean.before.control = mean(Control, na.rm = TRUE),
          mean.before.drought = mean(Drought, na.rm = TRUE)) %>%
  mutate(mean.before.control.0 = replace_na(mean.before.control,0),
         mean.before.drought.0 = replace_na(mean.before.drought,0))
# 1416 data points for 973 taxa, from 73 sites

after = data.2 %>%
  filter(n_treat_years == 1)
# plots with after data for year 1

after.1 = after %>%
  group_by(site_code, Taxon, year, block, plot, trt) %>%
  summarize(mean_cover = mean(max_cover))
# had to take mean here b/c some site have multiple measurements for same taxa in same trt in same block and plot

after.2 = pivot_wider(after.1, names_from = trt, values_from = mean_cover)

after.3 = after.2 %>%
  group_by(site_code, Taxon) %>%
  reframe(mean.after.control = mean(Control, na.rm = TRUE),
          mean.after.drought = mean(Drought, na.rm = TRUE)) %>%
  mutate(mean.after.control.0 = replace_na(mean.after.control,0),
         mean.after.drought.0 = replace_na(mean.after.drought,0))
# 1575 data points for 1045 taxa, from 81 sites

# merge before and after together

all.data = full_join(after.3,before.3)
# 1988 observations, 1274 species, 86 sites

all.data.2 = all.data[,c(1,2,5,6,9,10)]
all.data.2[is.na(all.data.2)] = 0
# some species in only after but not before and some in only before but not after

all.data.2 = all.data.2 %>%
  mutate(drought.after.before = mean.after.drought.0 - mean.before.drought.0,
         control.after.before = mean.after.control.0 - mean.before.control.0,
         cover.change = drought.after.before - control.after.before)

#write.csv(all.data.2, file = "./Formatted.Data/BACI.data.final.csv")

#### AusTraits ####
# getting trait data from AusTraits
austraits <- load_austraits(version = "4.1.0", path = "austraits")
austraits.2=readRDS("austraits-4.1.0.rds")

taxon.df=as.data.frame(unique(austraits.2$traits$taxon_name))

cover.response=read.csv("./Formatted.Data/BACI.data.final.csv")
cover.species=as.data.frame(unique(cover.response$Taxon))

cover.species.list=as.vector(cover.species$`unique(cover.response$Taxon)`)
cover.species.list.2=str_to_sentence(cover.species.list)

austraits.subset <- extract_taxa(austraits.2, taxon_name = cover.species.list.2)
austraits.traits <- extract_trait(austraits.subset, c("leaf_N_per_dry_mass","plant_height",
                                                      "root_diameter","root_N_per_dry_mass",
                                                      "root_specific_root_length",
                                                      "leaf_mass_per_area"))
austraits.subset.traits=austraits.traits$traits

# write.csv(austraits.subset.traits, file="AusTraits.subset.traits.csv")

#### Calculate drought severity index for sites ####

# has the amount of precipitation that fell during year 1 of the study at each site
yr1.ppt=read.csv("./Raw.Data/cover_ppt_2023-05-10.csv")
# ppt.1 column is the amount of ppt that fell in the plot in the 365 days prior to biomass collection
# MAP column is the mean annual precipitation at the site

drt.trt=read.csv("./Raw.Data/Site_Elev-Disturb.csv")
drt.trt.2 = drt.trt[,c(2,25)]
# drought_trt column has the targeted % water removal for the site

# formula: 
# 1. precip.drt = amount ppt received * % targeted removal (as decimal)
# 2. (precip.drt - MAP)/MAP

# merge drt.trt with yr1.ppt

all.drt.data = merge(yr1.ppt, drt.trt.2, by = c("site_code"))

all.drt.data = all.drt.data %>%
  filter(trt == "Drought") # just want drought plots since they had targeted removal

# turn drought_trt into decimal
all.drt.data$drought_trt_dec = as.numeric(all.drt.data$drought_trt)/100

all.drt.data$precip.drt = all.drt.data$ppt.1*all.drt.data$drought_trt_dec
all.drt.data$drt.sev.index = (all.drt.data$precip.dr - all.drt.data$map)/all.drt.data$map

# slim dataset to plots with 1 year of drought

all.drt.data.2 = subset(all.drt.data, all.drt.data$n_treat_years == 1)

# get the mean drt.sev.index for each site, this is just mean of drought plots
site.drt.sev.index = all.drt.data.2 %>%
  group_by(site_code) %>%
  reframe(mean.drt.sev.index = mean(drt.sev.index, na.rm = TRUE))

# write.csv(site.drt.sev.index, "./Formatted.Data/site.drt.dev.index.csv")

#### get functional group and life_form information for new BACI species ####

new.species = read.csv("./New.dfs/BACI.species.list.csv", row.names = 1)
colnames(new.species)="Taxon"
data=read.csv("./Raw.Data/IDE_cover_2023-01-02.csv") %>%
  select(Taxon,local_lifeform,functional_group)

species.merge = left_join(new.species,data, by = "Taxon")
species.merge.2 = unique(species.merge)

write.csv(species.merge.2, file = "./New.dfs/new.species.info.csv")
