# VISITc

Developer: Akihiko Ito (Professor at The University of Tokyo)

README prepared by: Akihiko Ito and Takuma Sakauchi (Graduate Student at The University of Tokyo)

## Description

Tags: Biogeochemical, land-surface-models, site, continent-n-global, 1D

VISIT (Vegetation Integrative SImulator for Trace gases) is a process-based terrestrial ecosystem model simulating carbon, nitrogen, and water cycles. The model comprises both plant and soil components in an ecosystem in an integrated manner, allowing us to simulate land-atmosphere biogeochemical interactions. The model evaluates soil emissions of major greenhouse gases (CO2, CH4 [wetland emission and upland absorption], and N2O) at point to global scales, taking account of climate and land-use conditions. In broad-scale simulations, each grid mesh is sub-divided into natural upland, natural wetland, and cropland, and calculated separately. The soil carbon dynamics is simulated by box-flow schemes, which differ in complexity from a 2-box scheme at global scale to 9-box scheme at point scale. The soil nitrogen dynamics is also simulated by a box-flow scheme, which is composed of organic and inorganic nitrogen pools. The soil water budget is simulated with a simple 2-layer hydrological scheme, considering water-holding capacity determined by soil texture. The model is one of the earliest models that account for the impact of soil erosion on carbon budget (Ito, 2007). The model has been used for various kinds of studies related to climate change, spanning from diagnosis of current greenhouse gas budget, scenario development, climate projection including carbon-cycle feedback, climatic impact assessment, and assessment of mitigation and management options. The model has been validated with field data such as chamber and tower-flux measurements; at broad scales, outputs of the model have been examined by using atmospheric data. The model is written in plain C language and can be compiled on any platform with a C compiler.

## Overall Architecture

![Overall architecture](docs/overall-architecture.png)

## Questions

For questions, please contact akihikoito@g.ecc.u-tokyo.ac.jp or use Issues and the organization Discussions feature.