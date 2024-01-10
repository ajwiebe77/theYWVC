# theYWVC
The Yukon Well Vulnerability Calculator RShiny App
10 Mar 2023

Licence: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)

Description

This project contains code related to the Yukon Well Vulnerability Calculator, an app available within the RShiny framework (https://ajwiebe77.shinyapps.io/the-yukon-well-vulnerability-calculator/; Wiebe and McKenzie, 2022 - http://dx.doi.org/10.22541/essoar.167267811.10671930/v1). The code is written in the R programming language (R Core Team, 2022).

The app employs the Nagheli et al. (2020) method to visualize a groundwater flow system (flow lines and lines of constant hydraulic head or equipotential) for aquifers bounded by surface water features (constant head boundaries) based on a small number of aquifer and well parameters. The method represents an aquifer as a simple geometric shape (e.g., rectangle or triangle). The goal of this project was to create a web application to allow users to employ the complicated mathematical solutions by Nagheli et al. (2020) without needing to code the method themselves. The app may be useful to illustrate what a groundwater flow system and estimate potential capture zone boundaries based on simplifying assumptions. It may also be used as an educational tool about groundwater flow to wells. The method draws capture zone boundaries for the wells represented, and the app allows users to download the results as GIS (geographical information systems) shapefiles (to be used in making maps, for instance, using QGIS - https://www.qgis.org/en/site/). Given the fact that generating a visualization with this app takes only a few minutes if parameters are available or can be estimated, it could serve as an initial guess at what a flow system looks like. Constructing a fully-distributed groundwater flow model to perform the same task could take months!

The Nagheli et al. (2020) method was originally coded in MATLAB, which offers some advantages such as identifying flow line segments related to capture zone boundaries. The decison to code the method in the R programming language was related to the ability to deploy the method as a web app with an interactive map in the RShiny framework.

The app is still under development and is not yet as robust as it could be. It also does not include all of the geometrical aquifer configurations that are discussed by Nagheli et al. (2020).


Usage

Copying all files to a folder and then opening and running the app.R file within RStudio software (https://posit.co/download/rstudio-desktop/) will allow testing and modification of the code on a user's computer.





References

Nagheli, S., Samani, N., Barry, D.A., 2020. Capture zone models of multi-well system in aquifers bounded with regular and irregular inflow boundaries. J. Hydrol. X 7, 100053. https://doi.org/10.1016/j.hydroa.2020.100053.

R Core Team, 2022. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/ (Accessed 17 Mar 2022).

Wiebe, A.J.*, and McKenzie, J.M., 2022. An Open-Source Web Tool for Visualizing Estimates of Well Capture Zones Near Surface Water Features. Poster presentation at: AGU Fall Meeting 2022, Chicago, IL, USA, 12-16 Dec 2022. http://dx.doi.org/10.22541/essoar.167267811.10671930/v1.
