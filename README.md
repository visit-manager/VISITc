# VISITc

Developer: Akihiko Ito (Professor at The University of Tokyo)

README prepared by: Akihiko Ito and Takuma Sakauchi (Graduate Student at The University of Tokyo)

## Description

Tags: Biogeochemical, land-surface-models, site, continent-n-global, 1D

VISIT (Vegetation Integrative Simulator for Trace gases) is a process-based terrestrial ecosystem model designed to simulate the coupled carbon, nitrogen, and water cycles. It integrates plant and soil processes within a single framework to represent biogeochemical interactions between land ecosystems and the atmosphere.

The model estimates major greenhouse-gas fluxes from soils, including CO2, CH4 (wetland emissions and upland uptake), and N2O, across scales ranging from site studies to continental and global applications. In large-scale simulations, each grid cell can be subdivided into natural upland, natural wetland, and cropland components, which are calculated separately.

Soil carbon dynamics are represented with box-flow schemes whose complexity ranges from a 2-box structure in global applications to a 9-box structure in point-scale studies. Soil nitrogen dynamics are likewise simulated with a box-flow approach that includes both organic and inorganic nitrogen pools. The soil water budget is represented with a simple two-layer hydrological scheme, with water-holding capacity determined by soil texture.

VISIT was among the earliest ecosystem models to account explicitly for the impact of soil erosion on the carbon budget. The model has been applied to a wide range of climate-change studies, including diagnosis of present-day greenhouse-gas budgets, scenario development, climate projections with carbon-cycle feedbacks, impact assessment, and evaluation of mitigation and land-management options. It has been evaluated against field observations such as chamber measurements and tower fluxes, and large-scale simulations have also been examined using atmospheric data. The model is written in standard C and can be compiled on any platform with a C compiler.

## Architecture

![Architecture](docs/overall-architecture.png)

## Scientific Articles

- Inatomi, M., Ito, A., Ishijima, K., and Murayama, S.: Greenhouse gas budget of a cool temperate deciduous broadleaved forest in Japan estimated using a process-based model, Ecosystems, 13, 472-483, doi:10.1007/s10021-010-9332-7, 2010.
- Ito, A.: Simulated impacts of climate and land-cover change on soil erosion and implication for the carbon cycle, 1901 to 2100, Geophys. Res. Lett., 34, doi:10.1029/2007GL029342, 2007.
- Ito, A., and Inatomi, M.: Use of a process-based model for assessing the methane budgets of global terrestrial ecosystems and evaluation of uncertainty, Biogeosciences, 9, 759-773, doi:10.5194/bg-9-759-2012, 2012.
- Ito, A., and Oikawa, T.: A simulation model of the carbon cycle in land ecosystems (Sim-CYCLE): A description based on dry-matter production theory and plot-scale validation, Ecological Modelling, 151, 147-179, 2002.
- Ito, A., Nishina, K., and Noda, H. M.: Evaluation of global warming impacts on the carbon budget of terrestrial ecosystems in monsoon Asia: a multi-model analysis, Ecological Research, 31, 459-474, doi:10.1007/s11284-016-1354-y, 2016.
- Ito, A., Nishina, K., Ishijima, K., Hashimoto, S., and Inatomi, M.: Emissions of nitrous oxide (N2O) from soil surfaces and their historical changes in East Asia: a model-based assessment, Progress in Earth and Planetary Science, 5, doi:10.1186/s40645-40018-40215-40644, 2018.

## Technical Information

- Operating system(s): MacOS, UNIX, Windows (C compiler required)
- Licence: freely available
- Output(s): plant carbon stock, productivity, soil carbon stock, soil respiration, ecosystem productivity, water budget, greenhouse gases (CO2, CH4, and N2O) budget, biomass burning emission, BVOC emission, soil erosion, DOC discharge
- Export format(s): binary or text output files

## Questions

For questions, please contact akihikoito@g.ecc.u-tokyo.ac.jp or use GitHub Issues or the organization Discussions feature.