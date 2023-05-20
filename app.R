# app.R
#
# Author: Andrew J. Wiebe, 20 May 2023
#
# The Yukon Well Vulnerability Calculator, version 1.3
#
# Code history: https://github.com/ajwiebe77/theYWVC
#
# Version 1.3 - Notes:
#     holzbecher_v10 works for case 1 option 1 and case 1 option 2 verification but not well for case 2 option 1 (could be because of the injection well)
#     Added a button on the setup tab so that the user can view the current aquifer geometry based on the settings
#     * Known Bug: need to remove the jagged flow lines
#     * Known Bug: changing cases in the verification results does not clear the graph if the flowlines etc have been drawn
#     The selection of flowline segments for the capture zones is now more robust (but not perfect).
# Version 1.2 - Notes:
#     Added ability to export shapefiles (in a zipped file) via a "Download" link on the Visualization tab
#     Synced the map centre and zoom level of the map in the Visualize tab with the map in the Setup tab
# Version 1.1 - Notes:
#     Regarding exporting shapefile points: the "Points to path" tool in, for example, QGIS (https://www.qgis.org/en/site/) can be used to generate lines from the sets of points using the "group" and "order" attributes.
#     Results graph on Visualize tab now matches the scale of the map in that tab
#
# Code ideas from the websites listed below (and in the code):
# https://www.r-bloggers.com/2018/11/accessing-openstreetmap-data-with-r/
# https://rstudio.github.io/leaflet/shiny.html
# https://rstudio.github.io/leaflet/basemaps.html


library(shiny)
library(leaflet)
library(stats) # for nagheli_case1test1test.R
library(nleqslv) # for nagheli_case1test1test.R
library(sp)
library(raster) # for exporting shapefiles (points)
library(foreign) # https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile
library(sf) ## transforming coordinates, loading shapefiles
library(data.table)

# FUNCTIONS -------------------------------------------

# List function files
# https://www.r-bloggers.com/2013/01/storing-a-function-in-a-separate-file-in-r/
## https://stackoverflow.com/questions/40716589/shiny-app-with-own-function
source("nagheli_case1fcn.R", local = TRUE)
source("findStg_case1.R", local = TRUE)
source("dispot_str1_v1deriv.R", local = TRUE)
source("dispot_str1_v1.R", local = TRUE)
source("nagheli_case2fcn.R", local = TRUE)
source("dispot_str2_v4.R", local = TRUE)
source("dispot_str2_v3deriv_imagX.R", local = TRUE)
source("holzbecher_v10.R", local = TRUE)
source("mergecheck.R", local = TRUE)
source("mergecheck2.R", local = TRUE)
source("checkSplits.R", local = TRUE)
source("trimAtBoundaries.R", local = TRUE)
source("long2UTM.R", local = TRUE)

# UI --------------------------------------------------

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "darkly"), # theme code from https://mastering-shiny.org/action-layout.html
  titlePanel("The Yukon Well Vulnerability Calculator"),
  tabsetPanel(
    tabPanel("Setup",
             sidebarLayout(
               sidebarPanel(
                 
                 sliderInput(inputId = "alpha_setup",
                             label = "Aquifer Angle (degrees)",
                             min = 10,
                             max = 180,
                             value = 149),
                 # radioButtons("radioComm_setup", h3("Community"),
                 #              choices = list("Carmacks" = 1, "Champagne Aishihik First Nation" = 2,
                 #                             "Dawson" = 3, "Haines Junction" = 4, "Kluane Lake First Nation" = 5,
                 #                             "Liard First Nation" = 6, "Mayo" = 7, "Selkirk First Nation" = 8,
                 #                             "Watson Lake" = 9, "Whitehorse" = 10, "White River First Nation" = 11),selected = 10),
                 radioButtons("radioComm_setup", h3("Community"),
                              choices = list("Whitehorse" = 1),selected = 1),
                 radioButtons("radioCase_setup", h3("Case"),
                              # choices = list("Case 1" = 1, "Case 2" = 2), selected=2),
                              choices = list("Case 1" = 1, "Case 2" = 2), selected=2)
                 # radioButtons("radioOpt_setup", h3("Option"),
                 #              choices = list("Option 1" = 1, "Option 2" = 2, "Option 3" = 3, "Option 4" = 4), selected=4)
               ),
               
               mainPanel(
                 leafletOutput("setupMap"), 
                 # add text box to display selected origin point on map
                 # verbatimTextOutput("origin"),
                 # Try adding text boxes specifying how to draw the graph overlain on the map
                 # ideas from: https://shiny.rstudio.com/articles/layout-guide.html
                 fluidRow(
                   column(3,
                          numericInput("originx_setup",
                                       h6("X (Graph Origin)"),
                                       value = 497811.8)
                   ),
                   column(3,
                          numericInput("originy_setup",
                                       h6("Y (Graph Origin)"),
                                       value = 6728905.4)
                   ),
                   column(3,
                          numericInput("rotate1",
                                       #h6("Rotation angle (degrees ccw from East)"),
                                       h6("Rotation angle (deg)"),
                                       value = -53.5)
                   ),
                   column(3,
                          numericInput("distance1",
                                       h6("Boundary length (m)"),
                                       value = 1595.5)
                   )
                 ),
                 fluidRow(
                   column(3,
                          numericInput("hyd_cond1",
                                       h6("Hydraulic conductivity (m/s)"),
                                       value = 0.001)
                   ),
                   column(3,
                          numericInput("hyd_grad1",
                                       h6("Hydraulic gradient (m/m)"),
                                       value = 0.0008)
                   ),
                   column(3,
                          numericInput("beta1",
                                       h6("Direction of gradient (deg)"),
                                       value = 180)
                   ),
                   column(3,
                          numericInput("aqf_thick1",
                                       h6("Aquifer thickness (m)"),
                                       value = 23)
                   )
                 ),
                 fluidRow(
                   column(3,
                          numericInput("Wellx1",
                                        h6("Well X (m)"),
                                        value = 497922)
                   ),
                   column(3,
                         numericInput("Welly1",
                                        h6("Well Y (m)"),
                                         value = 6730096)
                   ),
                   column(3,
                          numericInput("WellQ1",
                                       h6("Well Q (m^3/s)"),
                                       value = 0.115)
                   ),
                   column(3,
                          actionButton("setup_updateaqf", "Update Aquifer")
                   )
                 ),
                 fluidRow(
                   column(3,
                          numericInput("Wellx2",
                                       h6("Well X (m)"),
                                       value = 498494)
                   ),
                   column(3,
                          numericInput("Welly2",
                                       h6("Well Y (m)"),
                                       value = 6728987)
                   ),
                   column(3,
                          numericInput("WellQ2",
                                       h6("Well Q (m^3/s)"),
                                       value = 0.040)
                   )
                 ),
                 fluidRow(
                   column(3,
                          numericInput("Wellx3",
                                       h6("Well X (m)"),
                                       value = 0)
                   ),
                   column(3,
                          numericInput("Welly3",
                                       h6("Well Y (m)"),
                                       value = 0)
                   ),
                   column(3,
                          numericInput("WellQ3",
                                       h6("Well Q (m^3/s)"),
                                       value = 0)
                   )
                 ),
                 fluidRow(
                   column(3,
                          numericInput("Wellx4",
                                       h6("Well X (m)"),
                                       value = 0)
                   ),
                   column(3,
                          numericInput("Welly4",
                                       h6("Well Y (m)"),
                                       value = 0)
                   ),
                   column(3,
                          numericInput("WellQ4",
                                       h6("Well Q (m^3/s)"),
                                       value = 0)
                   )
                 ),
                 fluidRow(
                   column(3,
                          numericInput("Wellx5",
                                       h6("Well X (m)"),
                                       value = 0)
                   ),
                   column(3,
                          numericInput("Welly5",
                                       h6("Well Y (m)"),
                                       value = 0)
                   ),
                   column(3,
                          numericInput("WellQ5",
                                       h6("Well Q (m^3/s)"),
                                       value = 0)
                   )
                 )
               )
             )
    ),
    tabPanel("Visualize",
      sidebarLayout(
          sidebarPanel(
            
            verbatimTextOutput("vis_parameters_list")
          ),
         mainPanel(
           leafletOutput("outputMap"), 
           fluidRow(
             radioButtons("radioSHP_vis", h4("Export shapefile"),                             # use this line and the next for exporting
                                             choices = list("Yes" = 1, "No" = 2),selected=2),
             
            ),
           fluidRow(
             column(3,
                    actionButton("Calculate", "Calculate")
             ),
             column(3,
                    downloadLink("downloadSHPlink", label = "Download Shapefiles") # https://shiny.rstudio.com/reference/shiny/latest/downloadbutton
             )
           )
         )
      )
    ),
    tabPanel("Instructions",
             ### show image of Nagheli et al. (2020) geometric shapes
             h1("Overview"),
             p("The method used here to draw capture zones is based on a 2D analysis of a simplified flow system that contains surface water features (constant hydraulic head boundaries) and wells. The lines of constant hydraulic head (equipotential contours, analogous to the elevation contours on a topography map) and the streamlines that show water flow toward pumping wells or away from injection wells are drawn using equations from a paper by Nagheli et al. (2020). The user must enter the locations and pumping rates of the pumping or injection wells (+ for pumping, - for injection, perhaps counter-intuitively), the hydraulic conductivity of the aquifer (assumed to be homogeneous), and the aquifer thickness."),
             p("The geometric shape representing the aquifer must also be selected. This is currently restricted to simple triangles or rectangles. The geometry cases are illustrated below. The user may view the shape specified by the current parameters on the map in the Setup tab by clicking the Update Aquifer button. Blue lines are surface water (constant water level) boundaries, and black lines are the other outer edges of the aquifer, which is assumed to continue past the black lines, indefinitely."),
             p("The analytical equations used to calculate the equipotential and flow lines were implemented from Nagheli et al. (2020). The sequential contouring method of Holzbecher(2018) was used in an attempt to remove the branch cuts, i.e., spurious mathematical artifacts in the flowline contours."),
             h1("Instructions"),
             p("From the setup tab, enter the following information:"),
             p("The radio buttons on the left side can be used to zoom in to a particular community in the Yukon Territory. Zoom and pan operations in the map allow the user to select other geographical locations."),
             p("The slider on the left side at the top is used to set the angle of the triangle (alpha) for a case 2 aquifer (see below for aquifer types)."),
             p("Graph origin point on the map (this is where the (0,0) point on the dimensionless graph would be)."),
             p("A rotation angle for the geometric shape may also be specified to align the triangle or rectangle aquifer shape with the surface water features set as boundaries. The rotation angle must be specified as the angle between the x-axis in the triangle or rectangle and the angle of the x-axis in the rotated shape. The angle is specified in degrees counterclockwise from the horizontal; negative angles will rotate the shape clockwise."),
             p("The boundary length is the length of the x-axis at map scale."),
             p("Wells should not be placed too close together (either in terms of their y-coordinates or their x-coordinates). Flowline contours may not be drawn smoothly if this occurs."),
             p("After the input parameters are entered, the user may navigate to the 'Visualization tab' and click on the 'Calculate' button to initiate the calculations that will draw the flowfield (the user may first wish to choose the option of downloading the resulting shapefiles - see below). The app will convert the UTM coordinates (in metres) of the parameters to dimensionless units for the calculations and then convert them back to the actual map units when a shapefile is exported."),
             p("If the 'Export shapefile' radio button 'yes' is selected on the Visualization tab prior to clicking the Calculate button, then the user may download a zipped folder containing the shapefiles by clicking on the 'Download Shapefiles' link beside the Calculate button after the results are displayed on the map. The user should specify the extension .zip when downloading the folder. A text file listing the parameters used during the calculations is included in the downloaded folder."),
             h1("Parameters"),
             p("alpha - the angle between the two constant head boundaries (surface water features)"),
             p("Origin X - the UTM x-coordinate (easting) corresponding to the graph origin (0,0)"),
             p("Origin Y - the UTM y-coordinate (northing) corresponding to the graph origin (0,0)"),
             p("Rotation angle - the rotation angle by which to rotate the x-axis of the aquifer shape, in degrees counterclockwise from east"),
             p("boundary length - the x-axis distance on the map."),
             p("K (Hydraulic conductivity) - average hydraulic conductivity for the aquifer"),
             p("dh/dl (Hydraulic gradient) - the magnitude of the regional hydraulic gradient in the aquifer"),
             p("beta (Direction of gradient) - the angle toward which the regional gradient plunges, in degrees counterclockwise from east"),
             p("b (Aquifer thickness) - the average vertical thickness of the aquifer"),
             p("Well X - the UTM x-coordinate (easting) of the well (m)"),
             p("Well Y - the UTM y-coordinate (northing) of the well (m)"),
             p("Well Q - the average well pumping rate (m^3/s, that is, cubic metres per second)"),
             h1("Parameter Guidelines"),
             p("Reasonable hydraulic gradient values:"),
             p(tagList("Reasonable ranges of hydraulic conductivity values (", a("Freeze and Cherry, 1979", href = "https://fc79.gw-project.org/english/chapter-2/#2.3"), "):")),
             p("Gravel\t\t1E-3 to 1 m/s"),
             p("Clean Sand\t\t1E-5 to 1E-2 m/s"),
             p("Silty Sand\t\t1E-7 to 1E-3 m/s"),
             p("Reasonable pumping rates: BC MOE (2000) gives an example reasonable household pumping rate as 2270 L/day (0.000026 m^3/s); this would depend on the number of people in the household. One of many public supply wells for a town or city in a productive sand and gravel aquifer could have a pumping rate of around 0.080 m^3/s, for example (Matrix and SSPA, 2014)."),
             h1("Aquifer Geometry Options"),
             fluidRow(
               column(3,
                      img(src = "geometric aquifer shapes - case1.jpg", height = 100, width = 150)
                      
               ),
               column(9,
                      p(h2("Case 1")),
                      p("Nagheli et al. (2020) Case 1 - Rectangular Aquifer with constant head along one boundary and semi-infinite boundaries along the other sides.")
               )
             ),
             fluidRow(
               column(3,
                      img(src = "geometric aquifer shapes - case2a.jpg", height = 100, width = 150)
               ),
               column(9,
                      p(h2("Case 2a")),
                      p("Nagheli et al. (2020) Case 2a - Triangular Aquifer with constant head boundaries along two sides and a semi-infinite boundary along the third side.")
               )
             )
             
             
             
    ),
    tabPanel("Data for Yukon Territory",
             mainPanel(
                h1("Yukon Government Source Water Protection Reports and Information"),
                p(),
                h2("Whitehorse"),
                p(),
                p(tagList("", a("https://open.yukon.ca/data/datasets/water-well-capture-zones", href = "https://open.yukon.ca/data/datasets/water-well-capture-zones"), ".")),
                h2("Additional Reports"),
                p(tagList("", a("https://ftp-env-public.gov.yk.ca/WR/YWWR/", href = "https://ftp-env-public.gov.yk.ca/WR/YWWR/"), "."))
             )
    ),
             
    tabPanel("References",
             mainPanel(
               p(),
               h1("References"),
               p(tagList("British Columbia Ministry of Environment (BC MOE), 2000. Well protection toolkit.", a("https://www.env.gov.bc.ca/wsd/plan_protect_sustain/groundwater/wells/well_protection/wellprotect.html", href = "https://www.env.gov.bc.ca/wsd/plan_protect_sustain/groundwater/wells/well_protection/wellprotect.html")," (last access 22 July 2021).")),
               p(tagList("Freeze, A., Cherry, J., 1979. Groundwater. Prentice-Hall Inc., Englewood Cliffs, NJ, USA. Available from The Groundwater Project: ", a("https://fc79.gw-project.org/", href = "https://fc79.gw-project.org/"), " (last access 8 Feb 2023); CC BY-NC-ND 4.0. ")),
               p(tagList("Health Canada, 2020. Guidelines for Canadian Drinking Water Quality-Summary Table. Water and Air Quality Bureau, Healthy Environments and Consumer Safety Branch, Health Canada: Ottawa, ON.", a("https://www.canada.ca/en/health-canada/services/environmental-workplace-health/reports-publications/water-quality/guidelines-canadian-drinking-water-quality-summary-table.html", href = "https://www.canada.ca/en/health-canada/services/environmental-workplace-health/reports-publications/water-quality/guidelines-canadian-drinking-water-quality-summary-table.html"), ".")),
               p(tagList("Holzbecher, E., 2018. Streamline visualization of potential flow with branch cuts, with applications to groundwater. J. Flow Vis. Image Process. 25(2):119-144, ",a("https://doi.org/10.1615/JFlowVisImageProc.2018025918", href="https://doi.org/10.1615/JFlowVisImageProc.2018025918"),".")),
               ### https://stackoverflow.com/questions/42047422/create-url-hyperlink-in-r-shiny
               p(tagList("Matrix Solutions Inc. (Matrix), S.S. Papadopulos and Associates Inc. (SSPA), 2014. Region of Waterloo Tier Three Water Budget and Local Area Risk Assessment. Final Report, Sep. 2014. Prepared for: Region of Waterloo. ", a("https://www.sourcewater.ca/en/source-protection-areas/region-of-waterloo-tier-3.aspx", href="https://www.sourcewater.ca/en/source-protection-areas/region-of-waterloo-tier-3.aspx"), ".")),
               p(tagList("Nagheli, S., Samani, N., Barry, D.A., 2020. Capture zone models of a multi-well system in aquifers bounded with regular and irregular inflow boundaries. J. Hydrol. X 7, 100053,",a("https://doi.org/10.1016/j.hydroa.2020.100053", href="https://doi.org/10.1016/j.hydroa.2020.100053"), "."))
             )
    ),
    tabPanel("Verification",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("verifradioCase", h3("Case"),
                              # choices = list("Case 1 (square aquifer)" = 1, "Case 2 (triangular aquifer)" = 2), selected=1),
                              choices = list("Case 1 (square aquifer)" = 1, "Case 2 (triangular aquifer)" = 2), selected=1),
                 radioButtons("verifradioOpt", h3("Option"),
                              choices = list("Option 1" = 1, "Option 2" = 2, "Option 3" = 3, "Option 4" = 4), selected=1)
               ),
               
               mainPanel(
                 h2("Verification examples from Nagheli et al. (2020)"),
                 plotOutput(outputId = "verifplot"),
                 fluidRow(
                   column(3,
                          actionButton("verifCalculate", "Calculate")
                   )
                 ),
                 fluidRow(
                   column(12,
                          verbatimTextOutput("verif_aqf_details")
                   )
                 ),
                 fluidRow(
                   column(12,
                          verbatimTextOutput("verif_well_loc")
                   )
                 )
               )
             )
    ),
    tabPanel("About/Disclaimer",
             mainPanel(
               h1("About"),
               p("The Yukon Well Vulnerability Calculator - version 1.3 - 16 Mar 2022."),
               p("The goal of the Yukon Well Vulnerability Calculator app is to facilitate the estimation of risk to drinking water wells located near surface water features. The current version allows users to visualize capture zone boundaries for simple aquifer shapes (rectangles and triangles) bounded by surface water features on one or two sides."),# when potential contaminant sources are present within the capture zone. The first step allows the user to estimate/visualize the well capture zone."),# The second step requests locations of potential contaminant sources. The third step uses a 2D analysis to project what contaminant breakthrough curves might look like at the well under advective and simple dispersive transport. The focus in this case is on the contaminant benzene due to the prevalence of releases of petroleum hydrocarbons in the Yukon Territory (Yukon Government, 2022) and the strict standard for this substance in the Canadian Drinking Water Quality Guidelines (Health Canada, 2020)."),
               p("Well capture zones are often visualized after calculating water flow and/or solute transport with a 3D numerical model. Such models can take months to construct and hours to months to run. Analytical solutions such as the ones used in this app can run in minutes. For reference, other simplified capture zone estimation methods that involve solving equations are discussed by BC MOE (2000). The solutions of Nagheli et al. (2020) incorporate boundaries posed by surface water features."),
               p("Instructions for the use of the app are listed under the 'Instructions' tab."),
               p("Possible data sources for public well fields in communities in the Yukon Territory are listed under the 'Data for Yukon Territory' tab."),
               p("Literature references are listed under the 'References' tab."),
               p("The verification examples from Nagheli et al. (2020) are outlined in the 'Verification' tab."),
               p("The App is still under active development. If something does not work, check back in a week or two!"),
               h1("Code"),
               p(tagList("An archive of previous versions of the code is available at ", a("https://github.com/ajwiebe77/theYWVC", href = "https://github.com/ajwiebe77/theYWVC"),".")),
               h1("Disclaimer"),
               p("Disclaimer: The Yukon Well Vulnerability Calculator app ('The software') is intended as a visualization tool for showing what flow systems and capture zones might theoretically look like under highly simplified conditions. This software is provided 'AS IS' with no guarantee that the information generated is accurate or reasonable for the site or aquifer that is assessed herein. User discretion is advised. Neither the authors nor the funding agencies nor any affiliated organizations or institutions accept any liability related to the use or misuse of The Software. The user of the Yukon Well Vulnerability Calculator assumes all responsibilities for how this app or the code are used."),
               h1("Acknowledgements"),
               p("Funding for the development of this app was provided by the Canada First Research Excellence Fund (CFREF) Global Water Futures (GWF) program. Funding for the initial development was awarded to J.M. McKenzie, McGill University, for the project entitled, 'New Tools for Northern Groundwater Vulnerability Assessment.' Subsequent funding was through the GWF Core Modelling Team, under the supervision of Dr. D.L. Rudolph, University of Waterloo."),
               p(" "),
               fluidRow(
                 column(8,
                        img(src = "CFREF_logo.jpg", height = 94, width = 300),
                 ),
                 column(3,
                        img(src = "GWF_logo.jpg", height = 94, width = 94)
                 )
               ),
               p(" "),
               p("The idea to present this work as an app on the RShiny platform was provided by T. Hora."),
               p("Thanks to Dr. S. Nagheli, who provided insights into the process for calculating the capture zones that is outlined in the Nagheli et al. (2020) reference."),
               p("Support from the online communities discussing R coding is gratefully acknowledged. There are notes in the code where methods used were based on ideas from online discussions or blog posts."),
               h1("Author Contributions"),
               p("Conceptualization, methodology, and software coding: Andrew J. Wiebe, Postdoctoral Researcher, McGill University (2021-2022); Research Associate, University of Waterloo"),
               p("Supervision and funding acquisition: Dr. Jeffrey M. McKenzie, Professor, McGill University"),
               fluidRow(
                 column(3,
                        img(src = "McGill logo.png", height = 94, width = 400)
                 )
               ),
               p("Supervision and funding acquisition: Dr. David L. Rudolph, Professor, University of Waterloo"),
               fluidRow(
                 column(3,
                        img(src = "UWaterloo logo.png", height = 107, width = 400)
                 )
               ),
               h1("Citation"),
               p("The recommended citation for this app is:"),
               p("Wiebe, A.J., McKenzie, J.M., and Rudolph, D.L. 2022. The Yukon Well Vulnerability Calculator. Version: <version # above>. https://ajwiebe77.shinyapps.io/the-yukon-well-vulnerability-calculator/.")
            )
    )
    
  )
)

# SERVER ----------------------------------------------

server <- function(input, output, session) {
  
  SFlong <- 0
  
  ### SETUP -------------------------------------------------------------------
  
  output$vis_parameters_list <- renderText({
    if(input$radioCase_setup == 1){
      paste("Parameters\n------------\nCase =", input$radioCase_setup, "\nx-axis boundary length =", input$distance1 , "m\nRotation angle =", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
    }else if(input$radioCase_setup == 2){
      paste("Parameters\n------------\nCase =", input$radioCase_setup, "\nx-axis boundary length =", input$distance1 , "m\nalpha =",input$alpha_setup, "degrees\nRotation angle =", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
    }
  })
  
  output$setupMap <- renderLeaflet({
    
    if(input$radioComm_setup == 1){
      # Whitehorse:
      m_setup <- leaflet() %>% setView(lng = -135.056839, lat = 60.71, zoom = 12) # lat = 60.721188
      m_setup %>% addTiles()
    }
    
    
  })
  
  ### OUTPUT MAP ---------------------------------------------------------------
  output$outputMap <- renderLeaflet({
    
    ### generalize: load the same view as on the setup map (Setup tab)
    ### https://stackoverflow.com/questions/34985889/how-to-get-the-zoom-level-from-the-leaflet-map-in-r-shiny
    m_setup <- leaflet() %>% setView(lng = input$setupMap_center$lng, lat = input$setupMap_center$lat, zoom = input$setupMap_zoom)
    m_setup %>% addTiles()
    
    
  })
  
  output$verifplot <- renderPlot({
    # session$resetBrush("verifplot") # 14 Feb 2023
    plot(c(0,0), c(0,0), xlim=c(0,1),ylim=c(0,1))
    rect(0,0,1,1, density = NA, col="white", border="white") # test 16 Mar 2023
    # verifplot <- NULL # 14 Feb 2023
    # reset(verifplot) # 16 Mar 2023
    print("verifplot_A")
    
    K <- 5 / (60*60*24)
    b <- 20
    r <- 1
    
    if (input$verifradioCase == 1 && input$verifradioOpt == 1){
      xwD <- rbind(c(0.125, 0.9), c(0.375, 0.2), c(0.625, 0.5), c(0.75, 0.65), c(0.875, 0.3))
      QwDa <- cbind(c(0.01,0.02, 0.03, 0.01, 0.02))
      Qwa <- QwDa * K * b * r
      
      wells <- cbind(xwD[,1], xwD[,2], Qwa)
      extwells <- subset(wells, wells[,3]>0)
      
      plot(c(0,1), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))
      lines(c(0,0), c(0,1), lty = 2)
      lines(c(0,1), c(1,1), lty = 2)
      lines(c(1,1), c(1,0), lty = 2)
      points(extwells[,1], extwells[,2],pch=18, col="black")
      
    }else if(input$verifradioCase == 1 && input$verifradioOpt >= 2){ # in case options 3 or 4 are selected for case 1
      xwD <- rbind(c(0.125, 0.9), c(0.375, 0.2), c(0.625, 0.5), c(0.75, 0.65), c(0.875, 0.3))
      QwDb <- cbind(c(0.01, 0.02, -0.03, 0.01, 0.02))
      Qwb <- QwDb * K * b * r
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- subset(wells, wells[,3]>0)
      injwells <- subset(wells, wells[,3]<0)
      
      plot(c(0,1), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))
      lines(c(0,0), c(0,1), lty = 2)
      lines(c(0,1), c(1,1), lty = 2)
      lines(c(1,1), c(1,0), lty = 2)
      points(extwells[,1], extwells[,2],pch=18, col="black")
      points(injwells[,1], injwells[,2],pch=5)
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 1){
      xwD <- c(0.5, 0.15)
      QwDb <- c(0.01)
      Qwb <- QwDb * K * b * r
      alpha <- 0.19 * pi
      dh_dl <- 0.001
      beta <- 0.09 * pi
      
      wells <- matrix(c(xwD[1], xwD[2], Qwb), 1, 3) # wells needs to be recognized as a matrix (even if only one row)
      extwells <- wells
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
      
      plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))
      lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
      lines(outerarc[,1], outerarc[,2], lty = 2)
      points(extwells[1], extwells[2],pch=18, col="black")
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 2){
      xwD <- rbind(c(0.3, 0.2), c(0.5, 0.7))
      QwDb <- cbind(c(0.01, 0.02))
      Qwb <- QwDb * K * b * r
      alpha <- 0.42 * pi
      dh_dl <- 0.0
      beta <- 0.0
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- wells
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
      
      plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))
      lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
      lines(outerarc[,1], outerarc[,2], lty = 2)
      points(extwells[,1], extwells[,2],pch=18, col="black")
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 3){
      xwD <- rbind(c(0.3, 0.5), c(0.5, 0.2), c(0.7, 0.55))
      QwDb <- cbind(c(0.01, 0.02, 0.03))
      Qwb <- QwDb * K * b * r
      alpha <- pi / 2
      dh_dl <- 0.001
      beta <- pi / 4
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- wells
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
      
      plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))
      lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
      lines(outerarc[,1], outerarc[,2], lty = 2)
      points(extwells[,1], extwells[,2],pch=18, col="black")
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 4){
      xwD <- rbind(c(-0.2, 0.7), c(0, 0.3), c(0.2, 0.5), c(0.4, 0.4), c(0.7, 0.1))
      QwDb <- cbind(c(0.01, 0.02, 0.03, 0.01, 0.02))
      Qwb <- QwDb * K * b * r
      alpha <- 0.7 * pi
      dh_dl <- 0.0
      beta <- 0
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- wells
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
      
      plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(cos(alpha),1), ylim=c(0,1))
      lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
      lines(outerarc[,1], outerarc[,2], lty = 2)
      points(extwells[,1], extwells[,2],pch=18, col="black") 
    }
  })
  
  output$verif_aqf_details <- renderText({
    
    if (input$verifradioCase == 1 && input$verifradioOpt == 1){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: pi/2 (North)")
    }else if(input$verifradioCase == 1 && input$verifradioOpt == 2){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: pi/2 (North)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 1){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: 0.09pi (East 16.2 degrees North)\nTriangle angle alpha: 0.19pi (34.2 degrees)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 2){
      paste0("Dimensionless regional gradient magnitude: 0.0\nRegional gradient direction: 0.0 (East)\nTriangle angle alpha: 0.42pi (75.6 degrees)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 3){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: pi/4 (East 45 degrees North)\nTriangle angle alpha: pi/2 (90 degrees)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 4){
      paste0("Dimensionless regional gradient magnitude: 0.0\nRegional gradient direction: 0.0 (East)\nTriangle angle alpha: 0.7pi (126 degrees)")
    }
  })
  
  output$verif_well_loc <- renderText({
    # well_loc <- rbind(c("well number", 1,2,3,4,5), c("xwD", 0.125, 0.375, 0.625, 0.75, 0.875))
    if (input$verifradioCase == 1 && input$verifradioOpt == 1){
      paste0("Well number\txwD\tywD\tQwD\n1\t\t0.125\t0.9\t0.01\n2\t\t0.375\t0.2\t0.02\n3\t\t0.625\t0.5\t0.03\n4\t\t0.75\t0.65\t0.01\n5\t\t0.875\t0.3\t0.02")
    }else if(input$verifradioCase == 1 && input$verifradioOpt == 2){
      paste0("Well number\txwD\tywD\tQwD\n1\t\t0.125\t0.9\t0.01\n2\t\t0.375\t0.2\t0.02\n3\t\t0.625\t0.5\t-0.03\n4\t\t0.75\t0.65\t0.01\n5\t\t0.875\t0.3\t0.02")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 1){
      paste0("Well number\txwD\tywD\tQwD\n1\t\t0.5\t0.15\t0.01")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 2){
      paste0("Well number\txwD\tywD\tQwD\n1\t\t0.3\t0.2\t0.01\n2\t\t0.5\t0.7\t0.02")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 3){
      paste0("Well number\txwD\tywD\tQwD\n1\t\t0.3\t0.5\t0.01\n2\t\t0.5\t0.2\t0.02\n3\t\t0.7\t0.55\t0.03")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 4){
      paste0("Well number\txwD\tywD\tQwD\n1\t\t-0.2\t0.7\t0.01\n2\t\t0.0\t0.3\t0.02\n3\t\t0.2\t0.5\t0.03\n4\t\t0.4\t0.4\t0.01\n5\t\t0.7\t0.1\t0.02")
    }
  })
  
  ### VERIFICATION CALCULATIONS ------------------------------------------------
  observeEvent(input$verifCalculate, {
    
    print("verifCalc")
    
    yFLIP <- FALSE
    mirrorIMPERMEABLE <- FALSE
    
    K <- 5 / (60*60*24)
    b <- 20
    r <- 1
    
    if (input$verifradioCase == 1 && input$verifradioOpt == 1){
      option <- 1
      xwD <- rbind(c(0.125, 0.9), c(0.375, 0.2), c(0.625, 0.5), c(0.75, 0.65), c(0.875, 0.3))
      QwDa <- cbind(c(0.01,0.02, 0.03, 0.01, 0.02))
      Qwa <- QwDa * K * b * r
      dh_dl <- 0.001
      beta <- pi/2
      wells <- cbind(xwD[,1], xwD[,2], Qwa)
      extwells <- subset(wells, wells[,3]>0)
      injwells <- subset(wells, wells[,3]<0)
      
    }else if(input$verifradioCase == 1 && input$verifradioOpt >= 2){ # in case options 3 or 4 are selected for case 1
      option <- 2
      xwD <- rbind(c(0.125, 0.9), c(0.375, 0.2), c(0.625, 0.5), c(0.75, 0.65), c(0.875, 0.3))
      QwDb <- cbind(c(0.01, 0.02, -0.03, 0.01, 0.02))
      Qwb <- QwDb * K * b * r
      dh_dl <- 0.001
      beta <- pi/2
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- subset(wells, wells[,3]>0)
      injwells <- subset(wells, wells[,3]<0)
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 1){
      option <- 1
      xwD <- c(0.5, 0.15)
      QwDb <- c(0.01)
      Qwb <- QwDb * K * b * r
      alpha <- 0.19 * pi
      dh_dl <- 0.001
      beta <- 0.09 * pi
      
      # wells <- cbind(xwD[1], xwD[2], Qwb) # this is not recognized as a matrix, so nrow() crashes in holzbecher_v10() function
      wells <- matrix(c(xwD[1], xwD[2], Qwb), 1, 3)
      extwells <- wells
      injwells <- subset(wells, wells[,3]<0)
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 2){
      option <- 2
      xwD <- rbind(c(0.3, 0.2), c(0.5, 0.7))
      QwDb <- cbind(c(0.01, 0.02))
      Qwb <- QwDb * K * b * r
      alpha <- 0.42 * pi
      dh_dl <- 0.0
      beta <- 0.0
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- wells
      injwells <- subset(wells, wells[,3]<0)
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 3){
      option <- 3
      xwD <- rbind(c(0.3, 0.5), c(0.5, 0.2), c(0.7, 0.55))
      QwDb <- cbind(c(0.01, 0.02, 0.03))
      Qwb <- QwDb * K * b * r
      alpha <- pi / 2
      dh_dl <- 0.001
      beta <- pi / 4
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- wells
      injwells <- subset(wells, wells[,3]<0)
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 4){
      option <- 4
      xwD <- rbind(c(-0.2, 0.7), c(0, 0.3), c(0.2, 0.5), c(0.4, 0.4), c(0.7, 0.1))
      QwDb <- cbind(c(0.01, 0.02, 0.03, 0.01, 0.02))
      Qwb <- QwDb * K * b * r
      alpha <- 0.7 * pi
      dh_dl <- 0.0
      beta <- 0
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- wells
      injwells <- subset(wells, wells[,3]<0)
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
    }
    
    q0 <- dh_dl * K * b # regional uniform flow rate per unit width [L2/T]
    
    numpts <- 1000 # all cases - 17 May 2023
    
    if(input$verifradioCase == 1){
      consts <- cbind(c(r, beta, K, b, q0))
    }else if(input$verifradioCase == 2){
      consts <- cbind(c(r, alpha, beta, K, b, q0))
    }
    
    # call a function to get the stagnation points, discharge potential, and streamlines
    if(input$verifradioCase == 1){
      data1 <- nagheli_case1fcn(numpts,wells,consts)
    }else if(input$verifradioCase == 2){
      data1 <- nagheli_case2fcn(numpts,option,wells,consts)
    }
    
    xy_stg <- data1[[1]]
    zD_xy <- data1[[2]]
    phi <- data1[[3]]
    psi <- data1[[4]]
    
    # print(xy_stg)
    # print(min(psi[is.finite(psi)])) # https://stackoverflow.com/questions/35974779/finding-the-max-of-a-r-dataframe-column-ignoring-inf-and-na
    # print(max(psi[is.finite(psi)]))
    
    if(input$verifradioCase == 2 && option == 4){
      zD_xy2 <- matrix(0.0, (2*numpts + 1)*(2*numpts + 1),2)
    }else{
      zD_xy2 <- matrix(0.0, (numpts + 1)*(numpts + 1),2)
    }

    print("entering Holzbecher function")
    
    str_holz <- holzbecher_v10(psi, zD_xy, wells, consts, xy_stg) # 25 Apr 2023
    
    print("left Holzbecher function")
    psi_list <- str_holz[[1]]
    if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
      capzones <- str_holz[[2]]
    }

    print(length(psi_list))

    ### select equipotential lines

    numContoursDis <- 20
    phi_list <- contourLines(zD_xy[,1], zD_xy[,2], phi, nlevels=numContoursDis)

    ### trim the points to the aquifer boundaries
    psi_list <- trimAtBoundaries(psi_list, input$verifradioCase, option, alpha, TRUE) # labelled x and y vectors
    if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
      capzones <- trimAtBoundaries(capzones, input$verifradioCase, option, alpha, TRUE) # no labels on the x and y vectors
    }
    phi_list <- trimAtBoundaries(phi_list, input$verifradioCase, option, alpha, TRUE) # labelled x and y vectors
    
    
    ### PLOT VERIFICATION EXAMPLE -----------------------------------------------
    
    output$verifplot <- renderPlot({
      print("verifplot_B")
      plot(c(0,0), c(0,0), xlim=c(0,1),ylim=c(0,1))
      rect(0,0,1,1, density=5, col="red", border="white") # test 16 Mar 2023
      
      if (input$verifradioCase == 1 && input$verifradioOpt == 1){
        plot(c(0,1), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))
        
        for(i in 1:length(phi_list)){ # plot the streamlines
          lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
        }
        
        for(i in 1:length(psi_list)){ # plot the streamlines
          lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
        }
        
        if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
          for(i in 1:length(capzones)){ # plot the capture zone boundaries
            lines(capzones[[i]]$x, capzones[[i]]$y, col="black", lwd=2)
          }
        }
        
        points(extwells[,1],extwells[,2])
        if(nrow(injwells) > 1){
          points(injwells[,1], injwells[,2],pch=5)
        }else if(nrow(injwells) == 1){
          points(injwells[1], injwells[2],pch=5)
        }
        points(xy_stg[,1], xy_stg[,2], pch=18, col="green")
        
        lines(c(0,0), c(0,1), lty = 2)
        lines(c(0,1), c(1,1), lty = 2)
        lines(c(1,1), c(1,0), lty = 2)
        points(extwells[,1], extwells[,2],pch=18, col="black")
        
      }else if(input$verifradioCase == 1 && input$verifradioOpt >= 2){ # in case options 3 or 4 are selected for case 1
        
        plot(c(0,1), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))
        
        for(i in 1:length(phi_list)){ # plot the streamlines
          lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
        }
        
        for(i in 1:length(psi_list)){ # plot the streamlines
          lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
        }
        
        if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
          for(i in 1:length(capzones)){ # plot the capture zone boundaries
            lines(capzones[[i]]$x, capzones[[i]]$y, col="black", lwd=2)
          }
        }
        
        points(extwells[,1],extwells[,2])
        if(nrow(injwells) > 1){
          points(injwells[,1], injwells[,2],pch=5)
        }else if(nrow(injwells) == 1){
          points(injwells[1], injwells[2],pch=5)
        }
        points(xy_stg[,1], xy_stg[,2], pch=18, col="green")

        lines(c(0,0), c(0,1), lty = 2)
        lines(c(0,1), c(1,1), lty = 2)
        lines(c(1,1), c(1,0), lty = 2)
        points(extwells[,1], extwells[,2],pch=18, col="black")

      }else if(input$verifradioCase == 2 && input$verifradioOpt == 1){

        plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))

        for(i in 1:length(phi_list)){ # plot the streamlines
          lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
        }

        for(i in 1:length(psi_list)){ # plot the streamlines
          lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
        }

        if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
          for(i in 1:length(capzones)){ # plot the capture zone boundaries
            lines(capzones[[i]]$x, capzones[[i]]$y, col="black", lwd=2)
          }
        }

        if(nrow(extwells) == 1){
          points(extwells[1],extwells[2],pch=18, col="black")
        }else if(nrow(extwells) > 1){
          points(extwells[,1],extwells[,2],pch=18, col="black")
        }
        points(xy_stg[,1], xy_stg[,2], pch=18, col="green")
        
        alpha <- 0.19 * pi
        outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
        
        lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
        lines(outerarc[,1], outerarc[,2], lty = 2)

      }else if(input$verifradioCase == 2 && input$verifradioOpt == 2){

        plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))

        for(i in 1:length(phi_list)){ # plot the streamlines
          lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
        }

        for(i in 1:length(psi_list)){ # plot the streamlines
          lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
        }

        if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
          for(i in 1:length(capzones)){ # plot the capture zone boundaries
            lines(capzones[[i]]$x, capzones[[i]]$y, col="black", lwd=2)
          }
        }

        if(nrow(extwells) == 1){
          points(extwells[1],extwells[2],pch=18, col="black")
        }else if(nrow(extwells) > 1){
          points(extwells[,1],extwells[,2],pch=18, col="black")
        }
        points(xy_stg[,1], xy_stg[,2], pch=18, col="green")
        
        alpha <- 0.42 * pi
        outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
        
        lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
        lines(outerarc[,1], outerarc[,2], lty = 2)

      }else if(input$verifradioCase == 2 && input$verifradioOpt == 3){

        plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(0,1), ylim=c(0,1))

        for(i in 1:length(phi_list)){ # plot the streamlines
          lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
        }

        for(i in 1:length(psi_list)){ # plot the streamlines
          lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
        }

        if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
          for(i in 1:length(capzones)){ # plot the capture zone boundaries
            lines(capzones[[i]]$x, capzones[[i]]$y, col="black", lwd=2)
          }
        }

        if(nrow(extwells) == 1){
          points(extwells[1],extwells[2],pch=18, col="black")
        }else if(nrow(extwells) > 1){
          points(extwells[,1],extwells[,2],pch=18, col="black")
        }
        points(xy_stg[,1], xy_stg[,2], pch=18, col="green")
        
        alpha <- pi / 2
        outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
        
        lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
        lines(outerarc[,1], outerarc[,2], lty = 2)

      }else if(input$verifradioCase == 2 && input$verifradioOpt == 4){
        
        alpha <- 0.7 * pi
        plot(c(1,0), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis", xlim=c(cos(alpha),1), ylim=c(0,1))

        for(i in 1:length(phi_list)){ # plot the streamlines
          lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
        }

        for(i in 1:length(psi_list)){ # plot the streamlines
          lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
        }

        if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
          for(i in 1:length(capzones)){ # plot the capture zone boundaries
            lines(capzones[[i]]$x, capzones[[i]]$y, col="black", lwd=2)
          }
        }

        if(nrow(extwells) == 1){
          points(extwells[1],extwells[2],pch=18, col="black")
        }else if(nrow(extwells) > 1){
          points(extwells[,1],extwells[,2],pch=18, col="black")
        }
        points(xy_stg[,1], xy_stg[,2], pch=18, col="green")
        
        outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
        
        lines(c(0,cos(alpha)), c(0,sin(alpha)), type = "l", cex = 2)
        lines(outerarc[,1], outerarc[,2], lty = 2)
      }
    })
  })
  
  observeEvent(input$setup_updateaqf, {
    ### CODE TO ADD AQUIFER OUTLINE AND WELL POINTS TO THE SETUP MAP ----------------
    yFLIP <- FALSE
    # if(input$radioGeom_setup == 2){ # Next version
    #   yFLIP <- TRUE
    # }
    
    mirrorIMPERMEABLE <- FALSE
    
    r <- input$distance1
    numwells <- 0
    
    if(abs(input$WellQ1) > 0){
      numwells <- numwells + 1
      wells_temp <- cbind(input$Wellx1, input$Welly1, input$WellQ1)
    }
    
    if(abs(input$WellQ2) > 0){
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx2, input$Welly2, input$WellQ2))
    }
    
    if(abs(input$WellQ3) > 0){
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx3, input$Welly3, input$WellQ3))
    }
    
    if(abs(input$WellQ4) > 0){
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx4, input$Welly4, input$WellQ4))
    }
    
    if(abs(input$WellQ5) > 0){ # 16 Mar 2023
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx5, input$Welly5, input$WellQ5))
    }
    
    if(numwells == 1){
      wells_temp <- matrix(wells_temp, nrow=1)
    }
    
    wells_temp[,1] <- (wells_temp[,1] - input$originx_setup) / input$distance1 # note: also need to account for rotation!
    wells_temp[,2] <- (wells_temp[,2] - input$originy_setup) / input$distance1
    
    rot_matrix <- rbind(c(cos(- input$rotate1 * pi/180), - sin(- input$rotate1 * pi / 180)), c(sin(- input$rotate1 * pi / 180), cos(- input$rotate1 * pi / 180)))
    wells_temp2 <- matrix(0,numwells,3)
    
    for(i in 1:numwells){
      wells_temp2[i,1] <- rot_matrix[1,1] * wells_temp[i,1] + rot_matrix[1,2] * wells_temp[i,2]
      wells_temp2[i,2] <- rot_matrix[2,1] * wells_temp[i,1] + rot_matrix[2,2] * wells_temp[i,2]
    }
    
    wells_temp2[,3] <- wells_temp[,3]
    
    ### Determine case 2 option (based on Nagheli et al. examples)
    option <- 1
    if(input$radioCase_setup == 2 && input$alpha_setup < 90 && numwells == 1){
      option <- 1
    }else if(input$radioCase_setup == 2 && input$alpha_setup < 90 && numwells == 2){
      option <- 2
    }else if(input$radioCase_setup == 2 && input$alpha_setup == 90){
      option <- 3
    }else if(input$radioCase_setup == 2 && input$alpha_setup > 90 && input$alpha_setup < 180){
      option <- 4
    }else if(input$radioCase_setup == 2 && input$alpha_setup > 180 && input$alpha_setup <= 270){
      option <- 5 # not yet implemented
    }
    
    alpha <- input$alpha_setup * pi / 180
    beta <- input$beta1 * pi / 180
    # user may need to specify the ultimate angle, then the code will adjust (unrotate) for the dimensionless calcuations
    
    rot_angle <- input$rotate1 * pi/180
    rot_matrix <- rbind(c(cos(rot_angle), - sin(rot_angle)), c(sin(rot_angle), cos(rot_angle)))
    
    if(input$radioCase_setup == 1){
      ptA <- c(0, 0)
      ptB <- c(1, 0)
      boundary <- rbind(ptA, ptB) # adjusted points 6 Feb 2023
    }else if(input$radioCase_setup == 2){
      if(yFLIP == FALSE){
        ptA <- c(cos(alpha), sin(alpha))
        ptB <- c(0, 0)
        ptC <- c(1, 0)
        boundary <- rbind(ptA, ptB, ptC)
      }else{
        ### assume that the impermeable boundary is vertical for now
        ptA <- c(-1, tan(alpha))
        ptB <- c(0, 0)
        ptC <- c(-1, 0)
        boundary <- rbind(ptA, ptB, ptC)
      }
    }
    
    if(input$radioCase_setup == 2){
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
    }else if(input$radioCase_setup == 1){
      ptA <- c(0, 0)
      ptB <- c(0, 1)
      ptC <- c(1, 1)
      ptD <- c(1, 0)
      outerarc <- rbind(ptA, ptB, ptC, ptD)
    }
    
    extwells <- subset(wells_temp2, wells_temp2[,3]>0)
    injwells <- subset(wells_temp2, wells_temp2[,3]<0)
    
    ### rotate all points ---------------------
    boundary2 <- boundary
    
    for (i in 1:nrow(boundary)){
      boundary2[i,1] <- rot_matrix[1,1] * boundary[i,1] + rot_matrix[1,2] * boundary[i,2]
      boundary2[i,2] <- rot_matrix[2,1] * boundary[i,1] + rot_matrix[2,2] * boundary[i,2]
    }
    
    if(mirrorIMPERMEABLE == FALSE && (input$radioCase_setup == 1 || input$radioCase_setup == 2)){
      outerarc2 <- outerarc
      
      for(i in 1:nrow(outerarc)){
        outerarc2[i,1] <- rot_matrix[1,1] * outerarc[i,1] + rot_matrix[1,2] * outerarc[i,2]
        outerarc2[i,2] <- rot_matrix[2,1] * outerarc[i,1] + rot_matrix[2,2] * outerarc[i,2]
      }
    }
    
    extwells2 <- extwells
    
    for(i in 1:nrow(extwells)){
      extwells2[i,1] <- rot_matrix[1,1] * extwells[i,1] + rot_matrix[1,2] * extwells[i,2]
      extwells2[i,2] <- rot_matrix[2,1] * extwells[i,1] + rot_matrix[2,2] * extwells[i,2]
    }
    
    if(nrow(injwells) > 0){
      injwells2 <- injwells
      
      for(i in 1:nrow(injwells)){
        injwells2[i,1] <- rot_matrix[1,1] * injwells[i,1] + rot_matrix[1,2] * injwells[i,2]
        injwells2[i,2] <- rot_matrix[2,1] * injwells[i,1] + rot_matrix[2,2] * injwells[i,2]
      }
    }
    
    # ### scale the coordinates and translate about the UTM coordinates corresponding to the origin
    
    boundary2[,1] <- boundary2[,1] * r + input$originx_setup
    boundary2[,2] <- boundary2[,2] * r + input$originy_setup
    
    if(mirrorIMPERMEABLE == FALSE && (input$radioCase_setup == 1 || input$radioCase_setup == 2)){
      outerarc2[,1] <- outerarc2[,1] * r + input$originx_setup
      outerarc2[,2] <- outerarc2[,2] * r + input$originy_setup
    }
    
    extwells2[,1] <- extwells2[,1] * r + input$originx_setup
    extwells2[,2] <- extwells2[,2] * r + input$originy_setup
    
    if(nrow(injwells) > 0){
      injwells2[,1] <- injwells2[,1] * r + input$originx_setup
      injwells2[,2] <- injwells2[,2] * r + input$originy_setup
    }
    
    
    # z1 <- long2UTM(input$setupMap_center$lng) # NOT WORKING HERE
    z1 <- 8 # Need to adjust based on the map somehow # currently hardcoded for UTM Zone 8
    proj <- CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
    
    ### Add well points to map --------------
    
    df_extwells <- data.frame(x=extwells2[,1], y = extwells2[,2])
    if(nrow(injwells) > 0){
      df_injwells <- data.frame(x=injwells2[,1], y = injwells2[,2])
    }
    
    #Conversion of data frame to sf object
    sf_extwells <- st_as_sf(x = df_extwells,
                            coords = c("x", "y"),
                            crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
    )
    if(nrow(injwells) > 0){
      sf_injwells <- st_as_sf(x = df_injwells,
                              coords = c("x", "y"),
                              crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      )
    }
    
    #Projection transformation
    extwells_t = st_transform(sf_extwells, crs = "+proj=longlat +datum=WGS84")
    if(nrow(injwells) > 0){
      injwells_t = st_transform(sf_injwells, crs = "+proj=longlat +datum=WGS84")
    }
    ### Add aquifer boundaries to map ------------
    
    #### How to plot polylines: https://stackoverflow.com/questions/50867215/leaflet-add-multiple-polylines
    
    df_bdy_temp <- data.frame(x=boundary2[,1], y = boundary2[,2])
    
    #Conversion of data frame to sf object
    sf_bdy_temp <- st_as_sf(x = df_bdy_temp,
                            coords = c("x", "y"),
                            crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
    )
    
    #Projection transformation
    sf_bdy_temp_t = st_transform(sf_bdy_temp, crs = "+proj=longlat +datum=WGS84")
    sf_bdy_temp_tm <- as.matrix(sf_bdy_temp_t) # prefer to deal with matrices; access elements via sfc_m[[i]][1] for point i, x-coordinate, etc.
    
    bdy_segments <- matrix(0,(nrow(boundary2) - 1), 4)
    
    for(i in 1:(nrow(boundary2) - 1)){
      bdy_segments[i,1:2] <- c(sf_bdy_temp_tm[[i]][1], sf_bdy_temp_tm[[i]][2])
      bdy_segments[i,3:4] <- c(sf_bdy_temp_tm[[(i+1)]][1], sf_bdy_temp_tm[[(i+1)]][2])
    }
    
    df_bdy <- data.frame(x1 = bdy_segments[,1], y1 = bdy_segments[,2], x2 = bdy_segments[,3], y2 = bdy_segments[,4])
    colnames(df_bdy) <- c('x1', 'y1', 'x2', 'y2')
    
    setDT(df_bdy)
    
    ## create an 'id' / index value. This assumes each row of your data is a separate line.
    df_bdy[, idx := .I]
    
    ## create an `sfc` column (where each row is an `sfg` object)
    sf <- df_bdy[
      , {
        geometry <- sf::st_linestring(x = matrix(c(x1, y1, x2, y2), ncol = 2, byrow = T))
        geometry <- sf::st_sfc(geometry)
        geometry <- sf::st_sf(geometry = geometry)
      }
      , by = idx
    ]
    
    ## convert to sf
    sf_aqf <- sf::st_as_sf(sf)
    
    ### plot outer arc or outer line on map -----------------
    
    # if(input$radioCase_setup == 2 && mirrorIMPERMEABLE == FALSE){
    if(mirrorIMPERMEABLE == FALSE && (input$radioCase_setup == 1 || input$radioCase_setup == 2)){
      df_outer_temp <- data.frame(x = outerarc2[,1], y = outerarc2[,2])
    }else if(input$radioCase_setup == 2 && mirrorIMPERMEABLE == TRUE){
      df_outer_temp <- data.frame(x=c(boundary2[1,1],boundary2[3,1]), y = c(boundary2[1,2],boundary2[3,2]))
    }
    
    #Conversion of data frame to sf object
    sf_outer_temp <- st_as_sf(x = df_outer_temp,
                              coords = c("x", "y"),
                              crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
    )
    
    #Projection transformation
    sf_outer_temp_t = st_transform(sf_outer_temp, crs = "+proj=longlat +datum=WGS84")
    sf_outer_temp_tm <- as.matrix(sf_outer_temp_t) # prefer to deal with matrices; access elements via sfc_m[[i]][1] for point i, x-coordinate, etc.
    
    # outer_segments <- matrix(0,(nrow(boundary2) - 1), 4)
    if(input$radioCase_setup == 2 && mirrorIMPERMEABLE == FALSE){  # 6 Feb 2023
      outer_segments <- matrix(0,(nrow(outerarc2) - 1), 4)
    }else if(input$radioCase_setup == 2 && mirrorIMPERMEABLE == TRUE){
      outer_segments <- matrix(0,(nrow(boundary2) - 1), 4)
    }else if(input$radioCase_setup == 1 && mirrorIMPERMEABLE == FALSE){  # 6 Feb 2023
      outer_segments <- matrix(0,nrow(outerarc2), 4)
    }
    
    for(i in 1:(nrow(sf_outer_temp_tm) - 1)){ # 3 Feb 2023 - outer boundary was only plotting for a short section
      outer_segments[i,1:2] <- c(sf_outer_temp_tm[[i]][1], sf_outer_temp_tm[[i]][2])
      outer_segments[i,3:4] <- c(sf_outer_temp_tm[[(i+1)]][1], sf_outer_temp_tm[[(i+1)]][2])
    }
    
    df_outer <- data.frame(x1 = outer_segments[,1], y1 = outer_segments[,2], x2 = outer_segments[,3], y2 = outer_segments[,4])
    colnames(df_outer) <- c('x1', 'y1', 'x2', 'y2')
    
    setDT(df_outer)
    
    ## create an 'id' / index value. This assumes each row of your data is a separate line.
    df_outer[, idx := .I]
    
    ## create an `sfc` column (where each row is an `sfg` object)
    sf <- df_outer[
      , {
        geometry <- sf::st_linestring(x = matrix(c(x1, y1, x2, y2), ncol = 2, byrow = T))
        geometry <- sf::st_sfc(geometry)
        geometry <- sf::st_sf(geometry = geometry)
      }
      , by = idx
    ]
    
    ## convert to sf
    sf_outer <- sf::st_as_sf(sf)
    
    ### Update the leaflet map --------------------------------------------
    
    ### generalize: load the same view as on the setup map (Setup tab)
    ### https://stackoverflow.com/questions/34985889/how-to-get-the-zoom-level-from-the-leaflet-map-in-r-shiny
    output$setupMap <- renderLeaflet({
      m <- leaflet() %>% setView(lng = input$setupMap_center$lng, lat = input$setupMap_center$lat, zoom = input$setupMap_zoom)
      # m %>% addTiles() %>% addCircles(data=extwells_t, col = 'black', stroke = FALSE, fillOpacity = 0.8, radius = 6)
      if(nrow(injwells) > 0){
        m %>% addTiles()  %>% addPolylines(data = sf_aqf, col = 'blue') %>% addPolylines(data = sf_outer, col = 'black') %>% addCircles(data=injwells_t, col = 'black', fillOpacity = 1.0, radius = 6) %>% addCircles(data=extwells_t, col = 'black', stroke = FALSE, fillOpacity = 0.8, radius = 6)
      }else{
        m %>% addTiles()  %>% addPolylines(data = sf_aqf, col = 'blue') %>% addPolylines(data = sf_outer, col = 'black') %>% addCircles(data=extwells_t, col = 'black', stroke = FALSE, fillOpacity = 0.8, radius = 6)
      }
    })
    
    
  })
  
  ### CALCULATIONS -------------------------------------------------------------------
  observeEvent(input$Calculate, {
    
    print("calculate")
    
    yFLIP <- FALSE
    # if(input$radioGeom_setup == 2){
    #   yFLIP <- TRUE
    # }
    
    mirrorIMPERMEABLE <- FALSE
    
    # outer <- c(0,0)
    boundary <- c(0,0)
    
    #### test Whitehorse
    
    numwells <- 0
    
    if(abs(input$WellQ1) > 0){
      numwells <- numwells + 1
      wells_temp <- cbind(input$Wellx1, input$Welly1, input$WellQ1)
    }
    
    if(abs(input$WellQ2) > 0){
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx2, input$Welly2, input$WellQ2))
    }
    
    if(abs(input$WellQ3) > 0){
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx3, input$Welly3, input$WellQ3))
    }
    
    if(abs(input$WellQ4) > 0){
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx4, input$Welly4, input$WellQ4))
    }
    
    if(abs(input$WellQ5) > 0){ # 16 Mar 2023
      numwells <- numwells + 1
      wells_temp <- rbind(wells_temp, c(input$Wellx5, input$Welly5, input$WellQ5))
    }
    
    wells_temp <- as.matrix(wells_temp)
    print(wells_temp)
    
    wells_temp[,1] <- (wells_temp[,1] - input$originx_setup) / input$distance1 # note: also need to account for rotation!
    wells_temp[,2] <- (wells_temp[,2] - input$originy_setup) / input$distance1
    rot_matrix <- rbind(c(cos(- input$rotate1 * pi/180), - sin(- input$rotate1 * pi / 180)), c(sin(- input$rotate1 * pi / 180), cos(- input$rotate1 * pi / 180)))
    wells_temp2 <- matrix(0,numwells,3)
    for(i in 1:numwells){
      wells_temp2[i,1] <- rot_matrix[1,1] * wells_temp[i,1] + rot_matrix[1,2] * wells_temp[i,2]
      wells_temp2[i,2] <- rot_matrix[2,1] * wells_temp[i,1] + rot_matrix[2,2] * wells_temp[i,2]
    }
    
    wells_temp2[,3] <- wells_temp[,3]
    
    ### Determine case 2 option (based on Nagheli et al. examples)
    option <- 1
    if(input$radioCase_setup == 2 && input$alpha_setup < 90 && numwells == 1){
      option <- 1
    }else if(input$radioCase_setup == 2 && input$alpha_setup < 90 && numwells == 2){
      option <- 2
    }else if(input$radioCase_setup == 2 && input$alpha_setup == 90){
      option <- 3
    }else if(input$radioCase_setup == 2 && input$alpha_setup > 90 && input$alpha_setup < 180){
      option <- 4
    }else if(input$radioCase_setup == 2 && input$alpha_setup > 180 && input$alpha_setup <= 270){
      option <- 5 # not yet implemented
    }
      
    if(input$radioCase_setup == 2 && option == 4){
      numpts <- 500
    }else{
      numpts <- 1000
    }
    
    r <- input$distance1 # 1675.19
    K <- input$hyd_cond1 # 0.000243
    b <- input$aqf_thick1 # fixed 4 Jan 2023
    
    q0 <- input$hyd_grad1 * K * b # regional uniform flow rate per unit width [L2/T]
    
    alpha <- input$alpha_setup * pi / 180
    beta <- input$beta1 * pi / 180 # pi - (2.75 - pi/2)
    beta <- beta - input$rotate1 * pi/180
    
    if(input$radioCase_setup == 1){
      consts <- cbind(c(r, beta, K, b, q0))
    }else if(input$radioCase_setup == 2){
      consts <- cbind(c(r, alpha, beta, K, b, q0))
    }
    
    # call a function to get the stagnation points, discharge potential, and streamlines
    if(input$radioCase_setup == 1){
      data1 <- nagheli_case1fcn(numpts,wells_temp2,consts)
    }else if(input$radioCase_setup == 2){
      data1 <- nagheli_case2fcn(numpts,option,wells_temp2,consts)
    }
    
    xy_stg <- data1[[1]]
    zD_xy <- data1[[2]]
    phi <- data1[[3]]
    psi <- data1[[4]]
    
    print("zD_xy")
    str(zD_xy)
    
    rot_angle <- input$rotate1 * pi/180
    rot_matrix <- rbind(c(cos(rot_angle), - sin(rot_angle)), c(sin(rot_angle), cos(rot_angle)))
    
    if(input$radioCase_setup == 2 && option == 4){
      zD_xy2 <- matrix(0.0, (2*numpts + 1)*(2*numpts + 1),2)
    }else{
      zD_xy2 <- matrix(0.0, (numpts + 1)*(numpts + 1),2)
    }
    
    counter <- 1
    
    for(i in 1:nrow(zD_xy)){
      for(ii in 1:nrow(zD_xy)){ # there are only two columns, so use the number of rows again to make a square matrix
        zD_xy2[counter,1] <- rot_matrix[1,1] * zD_xy[ii,1] + rot_matrix[1,2] * zD_xy[i,2]
        zD_xy2[counter,2] <- rot_matrix[2,1] * zD_xy[ii,1] + rot_matrix[2,2] * zD_xy[i,2]
        counter <- counter + 1
      }
    }

    # ### scale the coordinates and translate about the UTM coordinates corresponding to the origin
    zD_xy2[,1] <- zD_xy2[,1] * r + input$originx_setup
    zD_xy2[,2] <- zD_xy2[,2] * r + input$originy_setup
    
    if(input$radioCase_setup == 1){
      ptA <- c(0, 0)
      ptB <- c(1, 0)
      boundary <- rbind(ptA, ptB) # adjusted points 6 Feb 2023
    }else if(input$radioCase_setup == 2){
      ptA <- c(cos(alpha), sin(alpha))
      ptB <- c(0, 0)
      ptC <- c(1, 0)
      boundary <- rbind(ptA, ptB, ptC)
    }
    
    if(input$radioCase_setup == 1){ # 6 Feb 2023
      ptA <- c(0, 0)
      ptB <- c(0, 1)
      ptC <- c(1, 1)
      ptD <- c(1, 0)
      outerarc <- rbind(ptA, ptB, ptC, ptD)
    }else if(input$radioCase_setup == 2){
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
    }
    
    xlimits <- c(-1,1)
    ylimits <- c(0,2)
    
    xlimits2 <- c(-2,2)
    ylimits2 <- c(-2,2)
    
    # ### Plot the data for the selected Case -----------------------
    
    extwells <- subset(wells_temp2, wells_temp2[,3]>0)
    injwells <- subset(wells_temp2, wells_temp2[,3]<0)
    
    ### remove the branch cuts
    print("entering Holzbecher function")
    str_holz <- holzbecher_v10(psi, zD_xy, wells_temp2, consts, xy_stg) # test - 25 Apr 2023
    print("left Holzbecher function")
    psi_list <- str_holz[[1]]
    if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
      capzones <- str_holz[[2]]
    }
    
    ### select equipotential lines
    
    numContoursDis <- 20
    phi_list <- contourLines(zD_xy[,1], zD_xy[,2], phi, nlevels=numContoursDis)
    
    ### trim the points to the aquifer boundaries
    psi_list <- trimAtBoundaries(psi_list, input$radioCase_setup, option, alpha, TRUE) # labelled x and y vectors
    if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
      capzones <- trimAtBoundaries(capzones, input$radioCase_setup, option, alpha, TRUE) # no labels on the x and y vectors
    }
    phi_list <- trimAtBoundaries(phi_list, input$radioCase_setup, option, alpha, TRUE) # labelled x and y vectors
    
    print(length(capzones))
    
    ### plot dimensionless data in RStudio window
    
    # numContoursStr <- 20
    # numContoursDis <- 20
    
    plot(boundary[,1], boundary[,2], lwd = 3, type = "l", col = "blue", xlim=xlimits, ylim=ylimits, xlab="x", ylab="y")
    
    for(i in 1:length(phi_list)){ # plot the streamlines
      lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
    }
    
    for(i in 1:length(psi_list)){ # plot the streamlines
      lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
    }
    
    if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
      for(i in 1:length(capzones)){ # plot the capture zone boundaries
        # lines(capzones[[i]][,1], capzones[[i]][,2], col="black", lwd=2)
        lines(capzones[[i]]$x, capzones[[i]]$y, col="black", lwd=2) # 19 Apr 2023
      }
    }
    
    lines(boundary[,1], boundary[,2], lwd = 3, col = "blue")
    lines(outerarc[,1],outerarc[,2], lwd = 3, lty = 2, col = "black")
    points(extwells[,1],extwells[,2])
    if(nrow(injwells) > 1){
      points(injwells[,1], injwells[,2],pch=5)
    }else if(nrow(injwells) == 1){
      points(injwells[1], injwells[2],pch=5)
    }
    points(xy_stg[,1], xy_stg[,2], pch=18, col="green")
    
    
    # ### rotate all points ---------------------
    boundary2 <- boundary
    
    for (i in 1:nrow(boundary)){
      boundary2[i,1] <- rot_matrix[1,1] * boundary[i,1] + rot_matrix[1,2] * boundary[i,2]
      boundary2[i,2] <- rot_matrix[2,1] * boundary[i,1] + rot_matrix[2,2] * boundary[i,2]
    }
    
    if(input$radioCase_setup == 1 || input$radioCase_setup == 2){ # 6 Feb 2023
      outerarc2 <- outerarc
      
      for(i in 1:nrow(outerarc)){
        outerarc2[i,1] <- rot_matrix[1,1] * outerarc[i,1] + rot_matrix[1,2] * outerarc[i,2]
        outerarc2[i,2] <- rot_matrix[2,1] * outerarc[i,1] + rot_matrix[2,2] * outerarc[i,2]
      }
    }
    
    extwells2 <- extwells
    
    for(i in 1:nrow(extwells)){
      extwells2[i,1] <- rot_matrix[1,1] * extwells[i,1] + rot_matrix[1,2] * extwells[i,2]
      extwells2[i,2] <- rot_matrix[2,1] * extwells[i,1] + rot_matrix[2,2] * extwells[i,2]
    }
    
    if(nrow(injwells) > 0){
      injwells2 <- injwells
      
      for(i in 1:nrow(injwells)){
        injwells2[i,1] <- rot_matrix[1,1] * injwells[i,1] + rot_matrix[1,2] * injwells[i,2]
        injwells2[i,2] <- rot_matrix[2,1] * injwells[i,1] + rot_matrix[2,2] * injwells[i,2]
      }
    }
    
    xy_stg2 <- xy_stg
    
    for(i in 1:nrow(xy_stg2)){
      xy_stg2[i,1] <- rot_matrix[1,1] * xy_stg[i,1] + rot_matrix[1,2] * xy_stg[i,2]
      xy_stg2[i,2] <- rot_matrix[2,1] * xy_stg[i,1] + rot_matrix[2,2] * xy_stg[i,2]
    }
    
    psi_list2 <- psi_list
    
    for(i in 1:length(psi_list2)){
      psi_list2[[i]]$x <- rot_matrix[1,1] * psi_list[[i]]$x + rot_matrix[1,2] * psi_list[[i]]$y
      psi_list2[[i]]$y <- rot_matrix[2,1] * psi_list[[i]]$x + rot_matrix[2,2] * psi_list[[i]]$y
    }
    
    if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
      capzones2 <- capzones
      
      for(i in 1:length(capzones2)){
        capzones2[[i]]$x <- rot_matrix[1,1] * capzones[[i]]$x + rot_matrix[1,2] * capzones[[i]]$y
        capzones2[[i]]$y <- rot_matrix[2,1] * capzones[[i]]$x + rot_matrix[2,2] * capzones[[i]]$y
      }
    }
    
    phi_list2 <- phi_list
    
    for(i in 1:length(phi_list2)){
      phi_list2[[i]]$x <- rot_matrix[1,1] * phi_list[[i]]$x + rot_matrix[1,2] * phi_list[[i]]$y
      phi_list2[[i]]$y <- rot_matrix[2,1] * phi_list[[i]]$x + rot_matrix[2,2] * phi_list[[i]]$y
    }
    
    ### plot rotated dimensionless data in RStudio window
    
    # plot(boundary2[,1], boundary2[,2], lwd = 3, type = "l", col = "black", xlim=xlimits2, ylim=ylimits2)
    # lines(outerarc2[,1],outerarc2[,2], lwd = 3, lty = 2, col = "black")
    # points(extwells2[,1],extwells2[,2])
    # if(nrow(injwells) > 0){
    #   points(injwells2[,1], injwells2[,2],pch=5)
    # }
    # points(xy_stg2[,1], xy_stg2[,2], pch=18, col="green")
    
    # ### scale the coordinates and translate about the UTM coordinates corresponding to the origin
    
    ### determine x and y bounds for graph to match the map zoom level
    z1 <- long2UTM(input$outputMap_center$lng)
    proj <- CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
    bnds <- input$outputMap_bounds
    print(bnds)
    
    # xycen <- as.data.frame(cbind(input$outputMap_center$lng, input$outputMap_center$lat))
    xycen <- st_point(c(input$outputMap_center$lng, input$outputMap_center$lat))
    xycen_sfc <- st_sfc(xycen, crs=4326)
    trnsf_xycen <- st_transform(xycen_sfc, crs=proj)
    
    xy1 <- st_point(c(bnds$west, bnds$south))
    xy1_sfc <- st_sfc(xy1, crs=4326)
    trnsf_xy1 <- st_transform(xy1_sfc, crs=proj)
    
    xy2 <- st_point(c(bnds$east, bnds$north))
    xy2_sfc <- st_sfc(xy2, crs=4326)
    trnsf_xy2 <- st_transform(xy2_sfc, crs=proj)
    
    xdist <- st_coordinates(trnsf_xy2)[1] - st_coordinates(trnsf_xy1)[1]
    ydist <- st_coordinates(trnsf_xy2)[2] - st_coordinates(trnsf_xy1)[2]
    
    xlimits2 <- c(st_coordinates(trnsf_xycen)[1] - floor(xdist/2), st_coordinates(trnsf_xycen)[1] + floor(xdist/2))
    ylimits2 <- c(st_coordinates(trnsf_xycen)[2] - floor(ydist/2), st_coordinates(trnsf_xycen)[2] + floor(ydist/2))
    
    ### 
    
    boundary2[,1] <- boundary2[,1] * r + input$originx_setup
    boundary2[,2] <- boundary2[,2] * r + input$originy_setup
    
    if(input$radioCase_setup == 1 || input$radioCase_setup == 2){
      outerarc2[,1] <- outerarc2[,1] * r + input$originx_setup
      outerarc2[,2] <- outerarc2[,2] * r + input$originy_setup
    }
    
    extwells2[,1] <- extwells2[,1] * r + input$originx_setup
    extwells2[,2] <- extwells2[,2] * r + input$originy_setup
    
    if(nrow(injwells) > 0){
      injwells2[,1] <- injwells2[,1] * r + input$originx_setup
      injwells2[,2] <- injwells2[,2] * r + input$originy_setup
    }
    
    xy_stg2[,1] <- xy_stg2[,1] * r + input$originx_setup
    xy_stg2[,2] <- xy_stg2[,2] * r + input$originy_setup
    
    for(i in 1:length(psi_list2)){
      psi_list2[[i]]$x <- psi_list2[[i]]$x * r + input$originx_setup
      psi_list2[[i]]$y <- psi_list2[[i]]$y * r + input$originy_setup
    }
    
    if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
      for(i in 1:length(capzones2)){
        capzones2[[i]]$x <- capzones2[[i]]$x * r + input$originx_setup
        capzones2[[i]]$y <- capzones2[[i]]$y * r + input$originy_setup
      }
    }
    
    for(i in 1:length(phi_list2)){
      phi_list2[[i]]$x <- phi_list2[[i]]$x * r + input$originx_setup
      phi_list2[[i]]$y <- phi_list2[[i]]$y * r + input$originy_setup
    }
    
    if(input$radioSHP_vis == 1){ ### export shapefiles
      # ### Save data as shapefiles -----------------------------------------
      z1 <- long2UTM(input$outputMap_center$lng)
      proj <- CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      
      ### Need lists of (x,y,phi) and (x,y,psi) points, not just the phi and psi matrices
      
      ### export equipotential lines
      counter <- 1
      xyphi <- matrix(data=0.0,numpts*numpts,4)
      
      for(i in 1:length(phi_list2)){
        for(ii in 1:length(phi_list2[[i]]$x)){
          if(is.na(phi_list2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
            xyphi[counter,1] <- phi_list2[[i]]$x[ii]
            xyphi[counter,2] <- phi_list2[[i]]$y[ii]
            xyphi[counter,3] <- i # group number
            xyphi[counter,4] <- ii # order number
            counter <- counter + 1
          }
        }
      }
      
      xyphi <- xyphi[c(1:counter-1), c(1,2,3,4),drop=FALSE] # remove zero rows
      
      equipotpts <- SpatialPoints(xyphi, proj4string = proj)
      # shapefile(equipotpts, filename = "equipotential.shp", overwrite = TRUE) # use the code below in the download handler
      # ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
      # dbfdata <- read.dbf("equipotential.dbf", as.is = TRUE)
      # # ## add new attribute data (just the numbers 1 to the number of objects)
      # dbfdata$group <- xyphi[,3] ### creates a new field called "group"
      # dbfdata$order <- xyphi[,4] ### creates a new field called "order"
      # write.dbf(dbfdata, "equipotential.dbf") ## overwrite the file with this new copy
      
      ### export flowlines
      counter <- 1
      xypsi <- matrix(data=0.0,numpts*numpts,4)
      
      for(i in 1:length(psi_list2)){
        for(ii in 1:length(psi_list2[[i]]$x)){
          if(is.na(psi_list2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
            xypsi[counter,1] <- psi_list2[[i]]$x[ii]
            xypsi[counter,2] <- psi_list2[[i]]$y[ii]
            xypsi[counter,3] <- i # group number
            xypsi[counter,4] <- ii # order number
            counter <- counter + 1
          }
        }
      }
      
      xypsi <- xypsi[c(1:counter-1), c(1,2,3,4),drop=FALSE] # remove zero rows
      
      flowlinepts <- SpatialPoints(xypsi, proj4string = proj)
      # shapefile(flowlinepts, filename = "flowlines.shp", overwrite = TRUE) # use the code below in the download handler
      ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
      # dbfdata <- read.dbf("flowlines.dbf", as.is = TRUE)
      # # ## add new attribute data (just the numbers 1 to the number of objects)
      # dbfdata$group <- xypsi[,3] ### creates a new field called "group"
      # dbfdata$order <- xypsi[,4] ### creates a new field called "order"
      # write.dbf(dbfdata, "flowlines.dbf") ## overwrite the file with this new copy
      
      ### export capture zone boundaries
      if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
        counter <- 1
        xycap <- matrix(data=0.0,numpts*numpts,4)
        
        for(i in 1:length(capzones2)){
          for(ii in 1:length(capzones2[[i]]$x)){
            if(is.na(capzones2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function

              xycap[counter,1] <- capzones2[[i]]$x[ii]
              xycap[counter,2] <- capzones2[[i]]$y[ii]
              xycap[counter,3] <- i # group number
              xycap[counter,4] <- ii # order number
              counter <- counter + 1
            }
          }
        }
        
        xycap <- xycap[c(1:counter-1), c(1,2,3,4),drop=FALSE] # remove zero rows
        
        cappts <- SpatialPoints(xycap, proj4string = proj)
        # shapefile(cappts, filename = "capzones.shp", overwrite = TRUE) # use the code below in the download handler
        ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
        # dbfdata <- read.dbf("capzones.dbf", as.is = TRUE)
        # # ## add new attribute data (just the numbers 1 to the number of objects)
        # dbfdata$group <- xycap[,3] ### creates a new field called "group"
        # dbfdata$order <- xycap[,4] ### creates a new field called "order"
        # write.dbf(dbfdata, "capzones.dbf") ## overwrite the file with this new copy
      }
      
      ### export stagnation points
      
      stgpts <- SpatialPoints(xy_stg2, proj4string = proj)
      
      ### export aquifer boundary points ### added 4 Jan 2023
      
      if(input$radioCase_setup == 1){
        aqfpts <- SpatialPoints(outerarc2, proj4string = proj)
      }
      
      if(input$radioCase_setup == 2){
        aqf <- matrix(0.0,(nrow(boundary2) + nrow(outerarc2) - 1),2) #rbind(boundary2, outerarc2[2:(length(outerarc2) - 1),])
        counter <- 1
        
        for(i in 1:nrow(boundary2)){
          aqf[counter,1] <- boundary2[i,1]
          aqf[counter,2] <- boundary2[i,2]
          counter <- counter + 1
        }
        
        for(i in 2:nrow(outerarc2)){ # include the last point so that the shape will close after the user uses "points to path" (it does not seem to matter if the unclosed path is used to generate a polygon)
          aqf[counter,1] <- outerarc2[i,1]
          aqf[counter,2] <- outerarc2[i,2]
          counter <- counter + 1
        }
        
        aqfpts <- SpatialPoints(aqf, proj4string = proj)
      }
      
      ### Note: the "Points to path" tool in QGIS can be used to generate lines from the sets of points using the "group" and "order" attributes
      
      ### Inspired by https://stackoverflow.com/questions/70386165/how-to-download-multiple-files-from-r-shiny-app
      ### and https://shiny.rstudio.com/reference/shiny/latest/downloadbutton
      
      output$downloadSHPlink <- downloadHandler(
        filename = function(){
          paste("gis_data_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file){
          
          temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
          dir.create(temp_directory)
          
          shapefile(equipotpts, filename = paste0(temp_directory, "/equipotential.shp"), overwrite = TRUE)
          # 4 Jan 2023
          ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
          dbfdata <- read.dbf(paste0(temp_directory, "/equipotential.dbf"), as.is = TRUE)
          ## add new attribute data (just the numbers 1 to the number of objects)
          dbfdata$group <- xyphi[,3] ### creates a new field called "group"
          dbfdata$order <- xyphi[,4] ### creates a new field called "order"
          write.dbf(dbfdata, paste0(temp_directory, "/equipotential.dbf")) ## overwrite the file with this new copy
          
          shapefile(flowlinepts, filename = paste0(temp_directory, "/flowlines.shp"), overwrite = TRUE)
          # 4 Jan 2023
          ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
          dbfdata <- read.dbf(paste0(temp_directory, "/flowlines.dbf"), as.is = TRUE)
          # ## add new attribute data (just the numbers 1 to the number of objects)
          dbfdata$group <- xypsi[,3] ### creates a new field called "group"
          dbfdata$order <- xypsi[,4] ### creates a new field called "order"
          write.dbf(dbfdata, paste0(temp_directory, "/flowlines.dbf")) ## overwrite the file with this new copy
          
          if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
            shapefile(cappts, filename = paste0(temp_directory, "/capzones.shp"), overwrite = TRUE)
            ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
            dbfdata <- read.dbf(paste0(temp_directory, "/capzones.dbf"), as.is = TRUE)
            # ## add new attribute data (just the numbers 1 to the number of objects)
            dbfdata$group <- xycap[,3] ### creates a new field called "group"
            dbfdata$order <- xycap[,4] ### creates a new field called "order"
            write.dbf(dbfdata, paste0(temp_directory, "/capzones.dbf")) ## overwrite the file with this new copy
          }
          
          shapefile(stgpts, filename = paste0(temp_directory, "/stagnationpts.shp"), overwrite = TRUE)
          
          shapefile(aqfpts, filename = paste0(temp_directory, "/aquiferpts.shp"), overwrite = TRUE)
          
          # https://statisticsglobe.com/write-lines-of-text-to-txt-file-in-r
          if(input$radioCase_setup == 1){ # add this 4 Jan 2023
            # text_lines <- paste("Parameters\n------------\nCase =", input$radioCase_setup, "\nx-axis boundary length =", input$distance1 , "m\nRotation angle =", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
            text_lines <- paste("The Yukon Well Vulnerability Calculator - v1.3\nParameters\n------------\nCase =", input$radioCase_setup, "\nx-axis boundary length =", input$distance1 , "m\nRotation angle =", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
          }else if(input$radioCase_setup == 2){
            # text_lines <- paste("Parameters\n------------\nCase =", input$radioCase_setup, "\nx-axis boundary length =", input$distance1 , "m\nalpha =",input$alpha_setup, "degrees\nRotation angle =", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
            text_lines <- paste("The Yukon Well Vulnerability Calculator - v1.3\nParameters\n------------\nCase =", input$radioCase_setup, "\nx-axis boundary length =", input$distance1 , "m\nalpha =",input$alpha_setup, "degrees\nRotation angle =", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
          }
          
          writeLines(text_lines, paste0(temp_directory, "/parameters.txt"))
          
          zip::zip(
            zipfile = file,
            files = dir(temp_directory),
            root = temp_directory
          )
        },
        contentType = "application/zip"
      )
      
    }
    
    
    z1 <- long2UTM(input$outputMap_center$lng)
    proj <- CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
    
    ### Prepare data for plotting on leaflet map
    
    # move the code out of the renderLeaflet command below to here ****
    
    
    ### LEAFLET OUTPUT MAP ---------------------------------------------------------------
    # Note: assumptions about the boundary make some of the following code specific to Case 2 (Triangular aquifer)
    print("update output map with shapefiles")
    
    output$outputMap <- renderLeaflet({
      ### update the output map on the visualization tab
      
      
      ### https://stackoverflow.com/questions/33045388/projecting-my-shapefile-data-on-leaflet-map-using-r
      # https://gis.stackexchange.com/questions/349317/working-in-projected-coordinates-rather-than-lat-long-inleaflet
      
      ### Add stagnation points to map --------------
      
      df_stg <- data.frame(x=xy_stg2[,1], y = xy_stg2[,2]) # works for one point
      
      #Conversion of data frame to sf object
      sf_stg <- st_as_sf(x = df_stg,
                         coords = c("x", "y"),
                         crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      )
      
      #Projection transformation
      stg_t = st_transform(sf_stg, crs = "+proj=longlat +datum=WGS84")
      
      
      ### Add well points to map --------------
      
      df_extwells <- data.frame(x=extwells2[,1], y = extwells2[,2])
      if(nrow(injwells) > 0){
        df_injwells <- data.frame(x=injwells2[,1], y = injwells2[,2])
      }
      
      #Conversion of data frame to sf object
      sf_extwells <- st_as_sf(x = df_extwells,
                              coords = c("x", "y"),
                              crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      )
      if(nrow(injwells) > 0){
        sf_injwells <- st_as_sf(x = df_injwells,
                                coords = c("x", "y"),
                                crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
        )
      }
      
      #Projection transformation
      extwells_t = st_transform(sf_extwells, crs = "+proj=longlat +datum=WGS84")
      if(nrow(injwells) > 0){
        injwells_t = st_transform(sf_injwells, crs = "+proj=longlat +datum=WGS84")
      }
      ### Add aquifer boundaries to map ------------
      
      #### How to plot polylines: https://stackoverflow.com/questions/50867215/leaflet-add-multiple-polylines
      
      df_bdy_temp <- data.frame(x=boundary2[,1], y = boundary2[,2])
      
      #Conversion of data frame to sf object
      sf_bdy_temp <- st_as_sf(x = df_bdy_temp,
                              coords = c("x", "y"),
                              crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      )
      
      #Projection transformation
      sf_bdy_temp_t = st_transform(sf_bdy_temp, crs = "+proj=longlat +datum=WGS84")
      sf_bdy_temp_tm <- as.matrix(sf_bdy_temp_t) # prefer to deal with matrices; access elements via sfc_m[[i]][1] for point i, x-coordinate, etc.
      
      bdy_segments <- matrix(0,(nrow(boundary2) - 1), 4)
      
      for(i in 1:(nrow(boundary2) - 1)){
        bdy_segments[i,1:2] <- c(sf_bdy_temp_tm[[i]][1], sf_bdy_temp_tm[[i]][2])
        bdy_segments[i,3:4] <- c(sf_bdy_temp_tm[[(i+1)]][1], sf_bdy_temp_tm[[(i+1)]][2])
      }
      
      df_bdy <- data.frame(x1 = bdy_segments[,1], y1 = bdy_segments[,2], x2 = bdy_segments[,3], y2 = bdy_segments[,4])
      colnames(df_bdy) <- c('x1', 'y1', 'x2', 'y2')
      
      setDT(df_bdy)
      
      ## create an 'id' / index value. This assumes each row of your data is a separate line.
      df_bdy[, idx := .I]
      
      ## create an `sfc` column (where each row is an `sfg` object)
      sf <- df_bdy[
        , {
          geometry <- sf::st_linestring(x = matrix(c(x1, y1, x2, y2), ncol = 2, byrow = T))
          geometry <- sf::st_sfc(geometry)
          geometry <- sf::st_sf(geometry = geometry)
        }
        , by = idx
      ]
      
      ## convert to sf
      sf_aqf <- sf::st_as_sf(sf)
      
      if(mirrorIMPERMEABLE == FALSE && (input$radioCase_setup == 1 || input$radioCase_setup == 2)){
        df_outer_temp <- data.frame(x = outerarc2[,1], y = outerarc2[,2])
      }else if(input$radioCase_setup == 2 && mirrorIMPERMEABLE == TRUE){
        df_outer_temp <- data.frame(x=c(boundary2[1,1],boundary2[3,1]), y = c(boundary2[1,2],boundary2[3,2]))
      }
      
      #Conversion of data frame to sf object
      sf_outer_temp <- st_as_sf(x = df_outer_temp,
                                coords = c("x", "y"),
                                crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      )
      
      #Projection transformation
      sf_outer_temp_t = st_transform(sf_outer_temp, crs = "+proj=longlat +datum=WGS84")
      sf_outer_temp_tm <- as.matrix(sf_outer_temp_t) # prefer to deal with matrices; access elements via sfc_m[[i]][1] for point i, x-coordinate, etc.
      
      # outer_segments <- matrix(0,(nrow(boundary2) - 1), 4)
      if(input$radioCase_setup == 2 && mirrorIMPERMEABLE == FALSE){  # 6 Feb 2023
        outer_segments <- matrix(0,(nrow(outerarc2) - 1), 4)
      }else if(input$radioCase_setup == 2 && mirrorIMPERMEABLE == TRUE){
        outer_segments <- matrix(0,(nrow(boundary2) - 1), 4)
      }else if(input$radioCase_setup == 1 && mirrorIMPERMEABLE == FALSE){  # 6 Feb 2023
        outer_segments <- matrix(0,nrow(outerarc2), 4)
      }
      
      for(i in 1:(nrow(sf_outer_temp_tm) - 1)){ # 3 Feb 2023 - outer boundary was only plotting for a short section
        outer_segments[i,1:2] <- c(sf_outer_temp_tm[[i]][1], sf_outer_temp_tm[[i]][2])
        outer_segments[i,3:4] <- c(sf_outer_temp_tm[[(i+1)]][1], sf_outer_temp_tm[[(i+1)]][2])
      }
      
      df_outer <- data.frame(x1 = outer_segments[,1], y1 = outer_segments[,2], x2 = outer_segments[,3], y2 = outer_segments[,4])
      colnames(df_outer) <- c('x1', 'y1', 'x2', 'y2')
      
      setDT(df_outer)
      
      ## create an 'id' / index value. This assumes each row of your data is a separate line.
      df_outer[, idx := .I]
      
      ## create an `sfc` column (where each row is an `sfg` object)
      sf <- df_outer[
        , {
          geometry <- sf::st_linestring(x = matrix(c(x1, y1, x2, y2), ncol = 2, byrow = T))
          geometry <- sf::st_sfc(geometry)
          geometry <- sf::st_sf(geometry = geometry)
        }
        , by = idx
      ]
      
      ## convert to sf
      sf_outer <- sf::st_as_sf(sf)
      
      
      ### plot equipotential lines -------------------------
      
      ### First, calculate how many pairs of points there are
      counter <- 0
      for(i in 1:length(phi_list2)){
        counter <- counter + (length(phi_list2[[i]]$x) - 1) # there are (length(phi_list2[[i]]$x) - 1) line segments between the length(phi_list2[[i]]$x) number of points
      }
      
      phi_segments <- matrix(0,(counter * 2), 2) # Need to compile as one single list of x and y pairs for conversion of spatial projection: pair the endpoints of the segments as follows: row 1 and row counter + 1, row 2 and row counter + 2, ..., row counter and row counter + counter; reconstruct the pairs later
      
      pcount <- 1
      
      ### create a list of pairs of points between which there are segments # could make this faster possibly by using blocks of x and y values for the sublists
      for(i in 1:length(phi_list2)){
        for(ii in 1:(length(phi_list2[[i]]$x) - 1)){
          if(is.na(phi_list2[[i]]$x[ii]) == FALSE && is.na(phi_list2[[i]]$x[(ii + 1)]) == FALSE){
            phi_segments[pcount,1:2] <- c(phi_list2[[i]]$x[ii], phi_list2[[i]]$y[ii])
            phi_segments[(pcount + counter),1:2] <- c(phi_list2[[i]]$x[(ii + 1)], phi_list2[[i]]$y[(ii + 1)])
            pcount <- pcount + 1
          }
        }
      }
      
      pcount <- pcount - 1
      
      phi_segments[(pcount + 1):(pcount * 2),1:2] <- phi_segments[(counter+1):(pcount + counter), 1:2] # overwrite missing values where NaNs would be
      phi_segments <- phi_segments[1:(pcount * 2),1:2] # remove empty part of revised matrix
      
      df_phi_temp <- data.frame(x = phi_segments[1:(pcount * 2),1], y = phi_segments[1:(pcount * 2),2])
      
      #Conversion of data frame to sf object
      sf_phi_temp <- st_as_sf(x = df_phi_temp,
                              coords = c("x", "y"),
                              crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      )
      
      #Projection transformation
      sf_phi_temp_t = st_transform(sf_phi_temp, crs = "+proj=longlat +datum=WGS84")
      sf_phi_temp_tm <- as.matrix(sf_phi_temp_t) # prefer to deal with matrices; access elements via sfc_m[[i]][1] for point i, x-coordinate, etc.
      
      phi_segments2 <- matrix(0, pcount, 4)
      
      for(i in 1:pcount){
        phi_segments2[i,1:2] <- c(sf_phi_temp_tm[[i]][1], sf_phi_temp_tm[[i]][2])
        phi_segments2[i,3:4] <- c(sf_phi_temp_tm[[(i + pcount)]][1], sf_phi_temp_tm[[(i + pcount)]][2])
      }
      
      df_phi <- data.frame(x1 = phi_segments2[,1], y1 = phi_segments2[,2], x2 = phi_segments2[,3], y2 = phi_segments2[,4]) # reconstruct pairs of segment endpoints
      colnames(df_phi) <- c('x1', 'y1', 'x2', 'y2')
      
      setDT(df_phi)
      
      ## create an 'id' / index value. This assumes each row of your data is a separate line.
      df_phi[, idx := .I]
      
      ## create an `sfc` column (where each row is an `sfg` object)
      sf <- df_phi[
        , {
          geometry <- sf::st_linestring(x = matrix(c(x1, y1, x2, y2), ncol = 2, byrow = T))
          geometry <- sf::st_sfc(geometry)
          geometry <- sf::st_sf(geometry = geometry)
        }
        , by = idx
      ]
      
      ## convert to sf
      sf_phi <- sf::st_as_sf(sf)
      
      ### plot streamlines ------------------------------------------------------
      
      ### First, calculate how many pairs of points there are
      counter <- 0
      for(i in 1:length(psi_list2)){
        counter <- counter + (length(psi_list2[[i]]$x) - 1) # there are (length(psi_list2[[i]]$x) - 1) line segments between the length(psi_list2[[i]]$x) number of points
      }
      
      psi_segments <- matrix(0,(counter * 2), 2) # Need to compile as one single list of x and y pairs for conversion of spatial projection: pair the endpoints of the segments as follows: row 1 and row counter + 1, row 2 and row counter + 2, ..., row counter and row counter + counter; reconstruct the pairs later
      
      pcount <- 1
      
      ### create a list of pairs of points between which there are segments # could make this faster possibly by using blocks of x and y values for the sublists
      for(i in 1:length(psi_list2)){
        for(ii in 1:(length(psi_list2[[i]]$x) - 1)){
          if(is.na(psi_list2[[i]]$x[ii]) == FALSE && is.na(psi_list2[[i]]$x[(ii + 1)]) == FALSE){
            psi_segments[pcount,1:2] <- c(psi_list2[[i]]$x[ii], psi_list2[[i]]$y[ii])
            psi_segments[(pcount + counter),1:2] <- c(psi_list2[[i]]$x[(ii + 1)], psi_list2[[i]]$y[(ii + 1)])
            pcount <- pcount + 1
          }
        }
      }
      
      pcount <- pcount - 1
      
      psi_segments[(pcount + 1):(pcount * 2),1:2] <- psi_segments[(counter+1):(pcount + counter), 1:2] # overwrite missing values where NaNs would be
      psi_segments <- psi_segments[1:(pcount * 2),1:2] # remove empty part of revised matrix
      
      df_psi_temp <- data.frame(x = psi_segments[1:(pcount * 2),1], y = psi_segments[1:(pcount * 2),2])
      
      #Conversion of data frame to sf object
      sf_psi_temp <- st_as_sf(x = df_psi_temp,
                              coords = c("x", "y"),
                              crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
      )
      
      #Projection transformation
      sf_psi_temp_t = st_transform(sf_psi_temp, crs = "+proj=longlat +datum=WGS84")
      sf_psi_temp_tm <- as.matrix(sf_psi_temp_t) # prefer to deal with matrices; access elements via sfc_m[[i]][1] for point i, x-coordinate, etc.
      
      psi_segments2 <- matrix(0, pcount, 4)
      
      for(i in 1:pcount){
        psi_segments2[i,1:2] <- c(sf_psi_temp_tm[[i]][1], sf_psi_temp_tm[[i]][2])
        psi_segments2[i,3:4] <- c(sf_psi_temp_tm[[(i + pcount)]][1], sf_psi_temp_tm[[(i + pcount)]][2])
      }
      
      df_psi <- data.frame(x1 = psi_segments2[,1], y1 = psi_segments2[,2], x2 = psi_segments2[,3], y2 = psi_segments2[,4]) # reconstruct pairs of segment endpoints
      colnames(df_psi) <- c('x1', 'y1', 'x2', 'y2')
      
      setDT(df_psi)
      
      ## create an 'id' / index value. This assumes each row of your data is a separate line.
      df_psi[, idx := .I]
      
      ## create an `sfc` column (where each row is an `sfg` object)
      sf <- df_psi[
        , {
          geometry <- sf::st_linestring(x = matrix(c(x1, y1, x2, y2), ncol = 2, byrow = T))
          geometry <- sf::st_sfc(geometry)
          geometry <- sf::st_sf(geometry = geometry)
        }
        , by = idx
      ]
      
      ## convert to sf
      sf_psi <- sf::st_as_sf(sf)
      
      if(length(str_holz[[2]]) > 0){
        
        ### First, calculate how many pairs of points there are
        counter <- 0
        for(i in 1:length(capzones2)){
          counter <- counter + (length(capzones2[[i]]$x) - 1) # there are (length(capzones2[[i]][,1]) - 1) line segments between the length(capzones2[[i]][,1]) number of points
        }
        
        capzo_segments <- matrix(0,(counter * 2), 2) # Need to compile as one single list of x and y pairs for conversion of spatial projection: pair the endpoints of the segments as follows: row 1 and row counter + 1, row 2 and row counter + 2, ..., row counter and row counter + counter; reconstruct the pairs later
        
        pcount <- 1
        
        ### create a list of pairs of points between which there are segments # could make this faster possibly by using blocks of x and y values for the sublists
        for(i in 1:length(capzones2)){
          for(ii in 1:(length(capzones2[[i]]$x) - 1)){
            if(is.na(capzones2[[i]]$x[ii]) == FALSE && is.na(capzones2[[i]]$x[(ii + 1)]) == FALSE){
              capzo_segments[pcount,1:2] <- c(capzones2[[i]]$x[ii], capzones2[[i]]$y[ii])
              capzo_segments[(pcount + counter),1:2] <- c(capzones2[[i]]$x[(ii + 1)], capzones2[[i]]$y[(ii + 1)])
              pcount <- pcount + 1
            }
          }
        }
        
        pcount <- pcount - 1
        
        capzo_segments[(pcount + 1):(pcount * 2),1:2] <- capzo_segments[(counter+1):(pcount + counter), 1:2] # overwrite missing values where NaNs would be
        capzo_segments <- capzo_segments[1:(pcount * 2),1:2] # remove empty part of revised matrix
        
        df_capzo_temp <- data.frame(x = capzo_segments[1:(pcount * 2),1], y = capzo_segments[1:(pcount * 2),2])
        
        #Conversion of data frame to sf object
        sf_capzo_temp <- st_as_sf(x = df_capzo_temp,
                                  coords = c("x", "y"),
                                  crs = CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
        )
        
        #Projection transformation
        sf_capzo_temp_t = st_transform(sf_capzo_temp, crs = "+proj=longlat +datum=WGS84")
        sf_capzo_temp_tm <- as.matrix(sf_capzo_temp_t) # prefer to deal with matrices; access elements via sfc_m[[i]][1] for point i, x-coordinate, etc.
        
        capzo_segments2 <- matrix(0, pcount, 4)
        
        for(i in 1:pcount){
          capzo_segments2[i,1:2] <- c(sf_capzo_temp_tm[[i]][1], sf_capzo_temp_tm[[i]][2])
          capzo_segments2[i,3:4] <- c(sf_capzo_temp_tm[[(i + pcount)]][1], sf_capzo_temp_tm[[(i + pcount)]][2])
        }
        
        df_capzo <- data.frame(x1 = capzo_segments2[,1], y1 = capzo_segments2[,2], x2 = capzo_segments2[,3], y2 = capzo_segments2[,4]) # reconstruct pairs of segment endpoints
        colnames(df_capzo) <- c('x1', 'y1', 'x2', 'y2')
        
        setDT(df_capzo)
        
        ## create an 'id' / index value. This assumes each row of your data is a separate line.
        df_capzo[, idx := .I]
        
        ## create an `sfc` column (where each row is an `sfg` object)
        sf <- df_capzo[
          , {
            geometry <- sf::st_linestring(x = matrix(c(x1, y1, x2, y2), ncol = 2, byrow = T))
            geometry <- sf::st_sfc(geometry)
            geometry <- sf::st_sf(geometry = geometry)
          }
          , by = idx
        ]
        
        ## convert to sf
        sf_capzo <- sf::st_as_sf(sf)
      }
      
      ### Plot the results on the leaflet map --------------------------------------------
      
      ### generalize: load the same view as on the setup map (Setup tab)
      ### https://stackoverflow.com/questions/34985889/how-to-get-the-zoom-level-from-the-leaflet-map-in-r-shiny
      m <- leaflet() %>% setView(lng = input$setupMap_center$lng, lat = input$setupMap_center$lat, zoom = input$setupMap_zoom)
      
      if(nrow(injwells) > 0 &&  length(str_holz[[2]]) > 0){
        m %>% addTiles()  %>% addPolylines(data = sf_aqf, col = 'black') %>% addPolylines(data = sf_phi, col = 'blue', weight = 1) %>% addPolylines(data = sf_psi, col = 'red', weight = 1) %>% addPolylines(data = sf_capzo, col = 'black', weight = 2) %>% addCircles(data = stg_t, col = 'green', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addCircles(data=injwells_t, col = 'black', fillOpacity = 1.0, radius = 6) %>% addCircles(data=extwells_t, col = 'black', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addPolylines(data = sf_outer, dashArray = "10 10", col = 'black')
      }else if(nrow(injwells) > 0 && length(str_holz[[2]]) == 0){
        m %>% addTiles()  %>% addPolylines(data = sf_aqf, col = 'black') %>% addPolylines(data = sf_phi, col = 'blue', weight = 1) %>% addPolylines(data = sf_psi, col = 'red', weight = 1) %>% addCircles(data = stg_t, col = 'green', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addCircles(data=injwells_t, col = 'black', fillOpacity = 1.0, radius = 6) %>% addCircles(data=extwells_t, col = 'black', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addPolylines(data = sf_outer, dashArray = "10 10", col = 'black')
      }else if(length(str_holz[[2]]) > 0){
        m %>% addTiles()  %>% addPolylines(data = sf_aqf, col = 'black') %>% addPolylines(data = sf_phi, col = 'blue', weight = 1) %>% addPolylines(data = sf_psi, col = 'red', weight = 1) %>% addPolylines(data = sf_capzo, col = 'black', weight = 2) %>% addCircles(data = stg_t, col = 'green', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addCircles(data=extwells_t, col = 'black', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addPolylines(data = sf_outer, dashArray = "10 10", col = 'black')
      }else if(length(str_holz[[2]]) ==  0){
        m %>% addTiles()  %>% addPolylines(data = sf_aqf, col = 'black') %>% addPolylines(data = sf_phi, col = 'blue', weight = 1) %>% addPolylines(data = sf_psi, col = 'red', weight = 1) %>% addCircles(data = stg_t, col = 'green', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addCircles(data=extwells_t, col = 'black', stroke = FALSE, fillOpacity = 0.8, radius = 6) %>% addPolylines(data = sf_outer, dashArray = "10 10", col = 'black')
      }
      
      
      
    })
    
    #   arrowLen <- 0.2
    #   arrows(0, 0, arrowLen * cos(beta[option]), arrowLen * sin(beta[option]))
    
    
  })
  
}

shinyApp(ui, server)