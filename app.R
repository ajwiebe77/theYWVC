# app.R
#
# Author: Andrew J. Wiebe, 19 Dec 2022
#
# The Yukon Well Vulnerability Calculator, version 1.1
#
# Version 1.1 - Notes:
#     Regarding exporting shapefile points: the "Points to path" tool in, for example, QGIS (https://www.qgis.org/en/site/) can be used to generate lines from the sets of points using the "group" and "order" attributes.
#     Shapefile export now supported
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
source("holzbecher_v3.R", local = TRUE)
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
                              choices = list("Case 1" = 1, "Case 2" = 2), selected=2),
                 radioButtons("radioOpt_setup", h3("Option"),
                              choices = list("Option 1" = 1, "Option 2" = 2, "Option 3" = 3, "Option 4" = 4), selected=4)
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
                 # fluidRow(
                 #   column(12,
                 #          verbatimTextOutput("info")
                 #   )
                 # ),
                 # fluidRow(
                 #   column(12,
                 #          verbatimTextOutput("info2")
                 #   )
                 # )
                 # plotOutput(outputId = "plot1")
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
                   # ),
                   # column(3,
                   #        actionButton("plotpts", "Plot Points")
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
                                       h6("Well Q (m)"),
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
                                       h6("Well Q (m)"),
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
                                       h6("Well Q (m)"),
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
                                       h6("Well Q (m)"),
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
            
      #       sliderInput(inputId = "alpha_vis",
      #                   label = "Aquifer Angle (degrees)",
      #                   min = 10,
      #                   max = 180,
      #                   value = 149),
      #       radioButtons("radioComm_vis", h3("Community"),
      #                    choices = list("Carmacks" = 1, "Champagne Aishihik First Nation" = 2,
      #                                   "Dawson" = 3, "Haines Junction" = 4, "Kluane Lake First Nation" = 5,
      #                                   "Liard First Nation" = 6, "Mayo" = 7, "Selkirk First Nation" = 8,
      #                                   "Watson Lake" = 9, "Whitehorse" = 10, "White River First Nation" = 11),selected = 10),
      #       radioButtons("radioCase_vis", h3("Case"),
      #                    choices = list("Case 1" = 1, "Case 2" = 2), selected=2),
      #       radioButtons("radioOpt_vis", h3("Option"),
      #                    choices = list("Option 1" = 1, "Option 2" = 2, "Option 3" = 3, "Option 4" = 4), selected=4)
          ),
         mainPanel(
           leafletOutput("outputMap"), 
           # add text box to display selected origin point on map
            # verbatimTextOutput("origin"),
           # Try adding text boxes specifying how to draw the graph overlain on the map
           # ideas from: https://shiny.rstudio.com/articles/layout-guide.html
           p(),
           p(),
           fluidRow(
             # column(3,
             #        numericInput("originx",
             #            h6("X (Graph Origin)"),
             #            value = 0)
             #          ),
             #  column(3,
             #      numericInput("originy",
             #                   h6("Y (Graph Origin)"),
             #                   value = 0)
             #      ),
             column(3,
                   actionButton("Calculate", "Calculate")
             )
           ),
           p(),
           p(),
           fluidRow(
             #  column(3,
             #        numericInput("rotate",
             #              #h6("Rotation angle (degrees ccw from East)"),
             #              h6("Rotation angle (deg)"),
             #              value = 0)
             #        ),
             # column(3,
             #        numericInput("distance",
             #            h6("Boundary length (m)"),
             #            value = 400)
             #        ),
             # column(3,
                    # actionButton("exportSHP", "Export Shapefile")
             # )
             # radioButtons("radioSHP_vis", h3("Export shapefile"),
             #                                 choices = list("Yes" = 1, "No" = 2),selected=2)
            ),
           # fluidRow(
           #   column(12,
           #            verbatimTextOutput("info")
           #          )
           # ),
           # fluidRow(
           #   column(12,
           #          verbatimTextOutput("info2")
           #   )
           # )
          plotOutput(outputId = "visplot1")
         )
      )
    ),
    tabPanel("Instructions",
             ### show image of Nagheli et al. (2020) geometric shapes
             h1("Overview"),
             p("The method used here to draw capture zones is based on a 2D analysis of a simplified flow system that contains surface water features (constant hydraulic head boundaries) and wells. The lines of constant hydraulic head (equipotential contours, analogous to the elevation contours on a topography map) and the streamlines that show water flow toward pumping wells or away from injection wells are drawn using equations from a paper by Nagheli et al. (2020). The user must enter the locations and pumping rates of the pumping or injection wells (+ for pumping, - for injection, perhaps counter-intuitively), the hydraulic conductivity of the aquifer (assumed to be homogeneous), and the aquifer thickness."),
             p("The geometric shape representing the aquifer must also be selected. This is currently restricted to simple triangles or rectangles. The geometry cases are illustrated below."),
             p("The analytical equations used to calculate the equipotential and flow lines were implemented from Nagheli et al. (2020). The sequential contouring method of Holzbecher(2018) was used in an attempt to remove the branch cuts, i.e., spurious mathematical artifacts in the flowline contours."),
             h1("Instructions"),
             p("From the setup tab, enter the following information:"),
             p("The radio buttons on the left side can be used to zoom in to a particular community in the Yukon Territory. Zoom and pan operations in the map allow the user to select other geographical locations."),
             p("The slider on the left side at the top is used to set the angle of the triangle (alpha) for a case 2 aquifer (see below for aquifer types)."),
             p("Graph origin point on the map (this is where the (0,0) point on the dimensionless graph would be)."),
             p("A rotation angle for the geometric shape may also be specified to align the triangle or rectangle aquifer shape with the surface water features set as boundaries. The rotation angle must be specified as the angle between the x-axis in the triangle or rectangle and the angle of the x-axis in the rotated shape. The angle is specified in degrees counterclockwise from the horizontal; negative angles will rotate the shape clockwise."),
             p("The boundary length is the length of the x-axis at map scale."),
             p("Wells should not be placed too close together (either in terms of their y-coordinates or their x-coordinates). Flowline contours may not be drawn smoothly if this occurs."),
             p("After the input parameters are entered, the user may navigate to the 'Visualization tab' and click on the 'Calculate' button to initiate the calculations that will draw the flowfield. The app will convert the UTM coordinates (in metres) of the parameters to dimensionless units for the calculations and then convert them back to the actual map units when a shape file is exported."),
             h1("Parameters"),
             p("alpha - the angle between the two constant head boundaries (surface water features)"),
             p("K - average hydraulic conductivity for the aquifer"),
             p("dh/dl - the magnitude of the regional hydraulic gradient in the aquifer"),
             p("beta - the angle toward which the regional gradient plunges, in degrees counterclockwise from east"),
             p("aquifer thickness - the average vertical thickness of the aquifer"),
             h1("Parameter Guidelines"),
             p("Reasonable hydraulic gradient values:"),
             p("Reasonable hydraulic conductivity values:"),
             p("Reasonable pumping rates: BC MOE (2000) gives an example reasonable household pumping rate as 2270 L/day (0.026 L/s); this would depend on the number of people in the household. One of many public supply wells for a town or city in a productive sand and gravel aquifer could have a pumping rate of around 80 L/s, for example (Matrix and SSPA, 2014)."),
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
                # h2("Carmacks"),
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
               # p(tagList("Health Canada, 2020. Guidelines for Canadian Drinking Water Quality-Summary Table. Water and Air Quality Bureau, Healthy Environments and Consumer Safety Branch, Health Canada: Ottawa, ON.", a("https://www.canada.ca/en/health-canada/services/environmental-workplace-health/reports-publications/water-quality/guidelines-canadian-drinking-water-quality-summary-table.html", href = "https://www.canada.ca/en/health-canada/services/environmental-workplace-health/reports-publications/water-quality/guidelines-canadian-drinking-water-quality-summary-table.html"), ".")),
               p(tagList("Holzbecher, E., 2018. Streamline visualization of potential flow with branch cuts, with applications to groundwater. J. Flow Vis. Image Process. 25(2):119-144, ",a("https://doi.org/10.1615/JFlowVisImageProc.2018025918", href="https://doi.org/10.1615/JFlowVisImageProc.2018025918"),".")),
               ### https://stackoverflow.com/questions/42047422/create-url-hyperlink-in-r-shiny
               p(tagList("Matrix Solutions Inc. (Matrix), S.S. Papadopulos and Associates Inc. (SSPA), 2014. Region of Waterloo Tier Three Water Budget and Local Area Risk Assessment. Final Report, Sep. 2014. Prepared for: Region of Waterloo. ", a("https://www.sourcewater.ca/en/source-protection-areas/region-of-waterloo-tier-3.aspx", href="https://www.sourcewater.ca/en/source-protection-areas/region-of-waterloo-tier-3.aspx"), ".")),
               p(tagList("Nagheli, S., Samani, N., Barry, D.A., 2020. Capture zone models of a multi-well system in aquifers bounded with regular and irregular inflow boundaries. J. Hydrol. X 7, 100053,",a("https://doi.org/10.1016/j.hydroa.2020.100053", href="https://doi.org/10.1016/j.hydroa.2020.100053"), "."))
               # p(tagList("Yukon Government, 2022. Contaminated sites,", a("https://map-data.service.yukon.ca/GeoYukon/Environmental_Monitoring/Contaminated_Sites/", href = "https://map-data.service.yukon.ca/GeoYukon/Environmental_Monitoring/Contaminated_Sites/"), " (last access: 12 August 2022."))
             )
    ),
    tabPanel("Verification",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("verifradioCase", h3("Case"),
                              choices = list("Case 1 (square aquifer)" = 1, "Case 2 (triangular aquifer)" = 2), selected=1),
                 radioButtons("verifradioOpt", h3("Option"),
                              choices = list("Option 1" = 1, "Option 2" = 2, "Option 3" = 3, "Option 4" = 4), selected=1)
               ),
               
               mainPanel(
                 h1("Verification examples from Nagheli et al. (2020)"),
                 plotOutput(outputId = "verifplot"),
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
               p("The Yukon Well Vulnerability Calculator - version 1.1 - 19 Dec 2022."),
               p("The goal of the Yukon Well Vulnerability Calculator app is to facilitate the estimation of risk to drinking water wells located near surface water features. The current version allows users to visualize capture zone boundaries for simple aquifer shapes (rectangles and triangles) bounded by surface water features on one or two sides."),# when potential contaminant sources are present within the capture zone. The first step allows the user to estimate/visualize the well capture zone."),# The second step requests locations of potential contaminant sources. The third step uses a 2D analysis to project what contaminant breakthrough curves might look like at the well under advective and simple dispersive transport. The focus in this case is on the contaminant benzene due to the prevalence of releases of petroleum hydrocarbons in the Yukon Territory (Yukon Government, 2022) and the strict standard for this substance in the Canadian Drinking Water Quality Guidelines (Health Canada, 2020)."),
               p("Well capture zones are often visualized after calculating water flow and/or solute transport with a 3D numerical model. Such models can take months to construct and hours to months to run. Analytical solutions such as the ones used in this app can run in minutes. For reference, other simplified capture zone estimation methods that involve solving equations are discussed by BC MOE (2000). The solutions of Nagheli et al. (2020) incorporate boundaries posed by surface water features."),
               p("Instructions for the use of the app are listed under the 'Instructions' tab."),
               p("Possible data sources for public well fields in communities in the Yukon Territory are listed under the 'Data for Yukon Territory' tab."),
               p("Literature references are listed under the 'References' tab."),
               p("The verification examples from Nagheli et al. (2020) are outlined in the 'Verification' tab."),
               p("The App is still under active development. If something does not work, check back in a week or two!"),
               h1("Disclaimer"),
               # p("Disclaimer: The Yukon Well Vulnerability Calculator app ('The software') is intended as a visualization tool for showing what flow systems and capture zones might theoretically look like under highly simplified conditions, and for estimating a simple risk rating for a well based on specified potential contaminant locations and aquifer/well properties. This software is provided 'AS IS' with no guarantee that the information generated is accurate or reasonable for the site or aquifer that is assessed herein. User discretion is advised. Neither the authors nor the funding agencies nor any affiliated organizations or institutions accept any liability related to the use or misuse of The Software."),
               p("Disclaimer: The Yukon Well Vulnerability Calculator app ('The software') is intended as a visualization tool for showing what flow systems and capture zones might theoretically look like under highly simplified conditions. This software is provided 'AS IS' with no guarantee that the information generated is accurate or reasonable for the site or aquifer that is assessed herein. User discretion is advised. Neither the authors nor the funding agencies nor any affiliated organizations or institutions accept any liability related to the use or misuse of The Software."),
               h1("Acknowledgements"),
               p("Funding for the development of this app was provided by the Canada First Research Excellence Fund (CFREF) Global Water Futures (GWF) program. The specific project was awarded to J. McKenzie, McGill University, and entitled, 'New Tools for Northern Groundwater Vulnerability Assessment.'"),
               p(" "),
               fluidRow(
                 column(8,
                        img(src = "CFREF_logo.jpg", height = 94, width = 300),
                        # img(src = "GWF_logo.jpg", height = 94, width = 94)
                 ),
                 column(3,
                        img(src = "GWF_logo.jpg", height = 94, width = 94)
                 )
                 # ),
                 # column(10,
                 #        img(src = "GWF_logo.jpg", height = 94, width = 94)
                 # )
               ),
               p(" "),
               p("The idea to present this work as an app on the RShiny platform was provided by T. Hora."),
               p("Thanks to Dr. S. Nagheli, who provided insights into the process for calculating the capture zones that is outlined in the Nagheli et al. (2020) reference."),
               p("Support from the online communities discussing R coding is gratefully acknowledged. There are notes in the code where methods used were based on ideas from online discussions or blog posts."),
               h1("Author Contributions"),
               p("Conceptualization, methodology, and software coding: Andrew J. Wiebe, Postdoctoral Researcher, McGill University"),
               p("Supervision and funding acquisition: Dr. Jeffrey M. McKenzie, Professor, McGill University"),
               fluidRow(
                 column(3,
                        img(src = "McGill logo.png", height = 94, width = 400)
                 )
               ),
               h1("Citation"),
               p("The recommended citation for this app is:"),
               p("Wiebe, A.J., and McKenzie, J.M. 2022. The Yukon Well Vulnerability Calculator. CC-BY-... https://ajwiebe77.shinyapps.io/ywvc_v1/.")
            )
    )
    
  )
)

# SERVER ----------------------------------------------

server <- function(input, output, session) {
  
  # setupMap <- createLeafletMap(session, "setupMap")
  
  SFlong <- 0
  
  # cenLongit <- 0
  # cenLat <- 0
  
  ### SETUP -------------------------------------------------------------------
  
  output$vis_parameters_list <- renderText({
    if(input$radioCase_setup == 1){
      paste("Parameters\n------------\nCase =", input$radioCase_setup ,"\nRotation angle", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
    }else if(input$radioCase_setup == 2){
      paste("Parameters\n------------\nCase =", input$radioCase_setup ,"\nalpha =",input$alpha_setup, "degrees\nRotation angle", input$rotate1, "degrees\nK =", input$hyd_cond1, "m/s\ndh/dl =",input$hyd_grad1, "\nbeta =", input$beta1,"degrees\nAquifer thickness =", input$aqf_thick1, "m\n\nWell_x (m)  Well_y (m)  Q (m^3/s)\n", input$Wellx1, input$Welly1, input$WellQ1, "\n", input$Wellx2, input$Welly2, input$WellQ2, "\n", input$Wellx3, input$Welly3, input$WellQ3, "\n", input$Wellx4, input$Welly4, input$WellQ4, "\n", input$Wellx5, input$Welly5, input$WellQ5, "\n")
    }
  })
  
  output$setupMap <- renderLeaflet({
    
    if(input$radioComm_setup == 1){
      # Whitehorse:
      m_setup <- leaflet() %>% setView(lng = -135.056839, lat = 60.721188, zoom = 12)
      m_setup %>% addTiles()
    }
    # if(input$radioComm_setup == 1){
    #   # Carmacks:
    #   m_setup <- leaflet() %>% setView(lng = -136.2872, lat = 62.0883, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 2){
    #   # Champagne Aishihik First Nation:
    #   m_setup <- leaflet() %>% setView(lng = -135.75, lat = 60.85, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 3){
    #   # Dawson:
    #   m_setup <- leaflet() %>% setView(lng = -139.432037, lat = 64.060066, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 4){
    #   # Haines Junction:
    #   m_setup <- leaflet() %>% setView(lng = -137.50, lat = 60.76, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 5){
    #   # Kluane Lake First Nation, Burwash Landing:
    #   m_setup <- leaflet() %>% setView(lng = -139, lat = 61.36, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 6){
    #   # Liard First Nation, 2 Mile Community:
    #   m_setup <- leaflet() %>% setView(lng = -128.75, lat = 60.09, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 7){
    #   # Mayo
    #   m_setup <- leaflet() %>% setView(lng = -135.917, lat = 63.6, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 8){
    #   # Selkirk First Nation, Pelly Crossing
    #   m_setup <- leaflet() %>% setView(lng = -136.55, lat = 62.83, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 9){
    #   # Watson Lake
    #   m_setup <- leaflet() %>% setView(lng = -128.710913, lat = 60.062806, zoom = 12)
    #   m_setup %>% addTiles()
    # }else if(input$radioComm_setup == 10){
    #   # Whitehorse:
    #   m_setup <- leaflet() %>% setView(lng = -135.056839, lat = 60.721188, zoom = 12)
    #   m_setup %>% addTiles()
    # }else{
    #   # White River First Nation, Beaver Creek
    #   m_setup <- leaflet() %>% setView(lng = -140.875, lat = 62.38, zoom = 12)
    #   m_setup %>% addTiles()
    # }
    
  })
  
  ### OUTPUT MAP ---------------------------------------------------------------
  output$outputMap <- renderLeaflet({
    
    #lng <- 0
    
    #leaflet() %>%
      #addProviderTiles(providers$Stamen.TonerLite,
      #                 options = providerTileOptions(noWrap = TRUE)
      #) %>%
      #addMarkers(data = points())
    
    if(input$radioComm_setup == 1){
      # Whitehorse:
      m_setup <- leaflet() %>% setView(lng = -135.056839, lat = 60.721188, zoom = 12)
      m_setup %>% addTiles()
    }
    
    # if(input$radioComm_setup == 1){
    #   # Carmacks:
    #   m <- leaflet() %>% setView(lng = -136.2872, lat = 62.0883, zoom = 12)
    #   m %>% addTiles()
    #   #output$value <- renderPrint(input$radioComm)
    # }else if(input$radioComm_setup == 2){
    #   # Champagne Aishihik First Nation:
    #   m <- leaflet() %>% setView(lng = -135.75, lat = 60.85, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 3){
    #   # Dawson:
    #   m <- leaflet() %>% setView(lng = -139.432037, lat = 64.060066, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 4){
    #   # Haines Junction:
    #   m <- leaflet() %>% setView(lng = -137.50, lat = 60.76, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 5){
    #   # Kluane Lake First Nation, Burwash Landing:
    #   m <- leaflet() %>% setView(lng = -139, lat = 61.36, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 6){
    #   # Liard First Nation, 2 Mile Community:
    #   m <- leaflet() %>% setView(lng = -128.75, lat = 60.09, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 7){
    #   # Mayo
    #   m <- leaflet() %>% setView(lng = -135.917, lat = 63.6, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 8){
    #   # Selkirk First Nation, Pelly Crossing
    #   m <- leaflet() %>% setView(lng = -136.55, lat = 62.83, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 9){
    #   # Watson Lake
    #   m <- leaflet() %>% setView(lng = -128.710913, lat = 60.062806, zoom = 12)
    #   m %>% addTiles()
    # }else if(input$radioComm_setup == 10){
    #   # Whitehorse:
    #   m <- leaflet() %>% setView(lng = -135.056839, lat = 60.721188, zoom = 12)
    #   m %>% addTiles()
    # }else{
    #   # White River First Nation, Beaver Creek
    #   m <- leaflet() %>% setView(lng = -140.875, lat = 62.38, zoom = 12)
    #   m %>% addTiles()
    # }
    
    # Add the origin as a point?
    
    # the following does not work
    
    #z1 <- long2UTM(input$outputMap_center$lng)
    #utm.prj <- "+proj=utm +zone=z1 ellps=WGS84"
    #latlong.prj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    
    #utmbounds <- spTransform(input$outputMap_bounds, CRS(utm.prj))
    # utmbounds <- st_transform(input$outputMap_bounds, crs=CRS(utm.prj)) # would need to fix this
    
    #bx1 <- input$outputMap_bounds_lng1
    #bx2 <- input$outputMap_bounds_lng2
    #by1 <-input$outputMap_bounds_lng2
    #by2 <- input$outputMap_bounds_lng2
    
    #if and(bx1 < originx, originx < bx2, by1 < originy, originy < by2){
      #x1 <- 
      #y1 <- 
      #m %>% addCircles(x1, y1)
    #}
  })
  

  
  output$verifplot <- renderPlot({
    
    if (input$verifradioCase == 1 && input$verifradioOpt == 1){
      xwD <- rbind(c(0.125, 0.9), c(0.375, 0.2), c(0.625, 0.5), c(0.75, 0.65), c(0.875, 0.3))
      QwDa <- cbind(c(0.01,0.02, 0.03, 0.01, 0.02))
      
      K <- 5 / (60*60*24)
      b <- 20
      r <- 1
      Qwa <- QwDa * K * b * r
      
      wells <- cbind(xwD[,1], xwD[,2], Qwa)
      extwells <- subset(wells, wells[,3]>0)
      
      plot(c(0,1), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis")
      lines(c(0,0), c(0,1), lty = 2)
      lines(c(0,1), c(1,1), lty = 2)
      lines(c(1,1), c(1,0), lty = 2)
      points(extwells[,1], extwells[,2],pch=18, col="black")
      
    }else if(input$verifradioCase == 1 && input$verifradioOpt == 2){
      xwD <- rbind(c(0.125, 0.9), c(0.375, 0.2), c(0.625, 0.5), c(0.75, 0.65), c(0.875, 0.3))
      QwDb <- cbind(c(0.01, 0.02, -0.03, 0.01, 0.02))
      
      K <- 5 / (60*60*24)
      b <- 20
      r <- 1
      Qwb <- QwDb * K * b * r
      
      wells <- cbind(xwD[,1], xwD[,2], Qwb)
      extwells <- subset(wells, wells[,3]>0)
      injwells <- subset(wells, wells[,3]<0)
      
      plot(c(0,1), c(0,0), type = "l", cex = 2, xlab="Dimensionless x-axis", ylab="Dimensionless y-axis")
      lines(c(0,0), c(0,1), lty = 2)
      lines(c(0,1), c(1,1), lty = 2)
      lines(c(1,1), c(1,0), lty = 2)
      points(extwells[,1], extwells[,2],pch=18, col="black")
      points(injwells[,1], injwells[,2],pch=5)
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 1){
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 2){
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 3){
      
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 4){
      
    }
  })
  
  output$verif_aqf_details <- renderText({
    
    if (input$verifradioCase == 1 && input$verifradioOpt == 1){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: pi/2 (North)")
    }else if(input$verifradioCase == 1 && input$verifradioOpt == 2){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: pi/2 (North)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 1){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: 0.09pi (East 16.2 degrees North)\nTriangle angle alpha: 0.19pi (East 34.2 degrees North)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 2){
      paste0("Dimensionless regional gradient magnitude: 0.0\nRegional gradient direction: 0.0 (East)\nTriangle angle alpha: 0.42pi (East 75.6 degrees North)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 3){
      paste0("Dimensionless regional gradient magnitude: 0.001\nRegional gradient direction: pi/4 (East 45 degrees North)\nTriangle angle alpha: pi/2 (North)")
    }else if(input$verifradioCase == 2 && input$verifradioOpt == 4){
      paste0("Dimensionless regional gradient magnitude: 0.0\nRegional gradient direction: 0.0 (East)\nTriangle angle alpha: 0.7pi (North 36 degrees West")
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
  
  # add origin point clicked on to text box - did not work to add click code to the leaflet map
  #output$origin <- renderText({
  #  paste0("x=", input$map_click$x, "\ny=", input$map_click$y)
  #})
  
  # observeEvent(input$updateCentre, {
  # xy <- eventReactive(input$updateCentre, {
  #   cenLongit <- input$setupMap_center$lng
  #   cenLat <- input$setupMap_center$lat
  # })
    
  output$info <- renderText({
      # Convert the lat long to UTM zone
      # https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
  
      # https://rstudio.github.io/leaflet/shiny.html
  
      #paste0(input$outputMap_center)
      # WORKS!
  
      # cenLong <- input$setupMap_center$lng
      # paste0(long2UTM(cenLong)) ## This prints the UTM Zone based on the function at the top of this file
      # WORKS!
      
      #bnds <- input$outputMap_bounds
      #paste0(bnds)
      # WORKS!
      
      # if(input$updateCentre == 0){
      #   return()
      # }
      
      input$updateCentre
    
      isolate(cenLongit <- input$setupMap_center$lng)
      isolate(cenLat <- input$setupMap_center$lat)
      zone = long2UTM(cenLongit)
      # paste0(cenLongit)
      # xy <- data.frame(ID = 1:length(cenLat), X = cenLongit, Y = cenLat)
      
      if(input$updateCentre == 0){
        return()
      }else{
        # xy <- as.data.frame(cbind(cenLongit, cenLat))
        # # str(xy)
        # colnames(xy) <- c('X', 'Y')
        # coordinates(xy) <- c("X", "Y")
        xy <- st_point(c(cenLongit, cenLat))
        # proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
        # trnsf_xy <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," +ellps=WGS84 +datum=WGS84", sep='')))
        xy_sfc <- st_sfc(xy, crs=4326)
        trnsf_xy <- st_transform(xy_sfc, crs=CRS(paste("+proj=utm +zone=",zone," +ellps=WGS84 +datum=WGS84", sep='')))
        # paste0(trnsf_xy$X, ", ", trnsf_xy$Y)
        paste0(st_coordinates(trnsf_xy)[1], ", ", st_coordinates(trnsf_xy)[2])
      }
      
      
      # if(input$updateCentre)
      # observeEvent(input$updateCentre){
      #   paste0(cenLongit)
      # }
      # })
      # cenLongit <- input$setupMap_center$lng
      # paste0(cenLongit)
      # zone = long2UTM(cenLongit)
      
      # cenLat <- input$setupMap_center$lat
      # paste0(cenLongit, ", ", cenLat)
      
      # print(cenLongit)
      # print(cenLat)
      # print(zone)
      # print(length(cenLat))
      
      ### The following code was modified from:
      ### https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
      # data.df<-as.data.frame(cbind(intp$x, intp$y, intp$z)) # example
      # colnames(data.df)<-c('x','y','z')  # example
      # xy <- data.frame(ID = 1:length(cenLat), X = cenLongit, Y = cenLat)
      # xy <- as.data.frame(cbind(cenLongit, cenLat))
      # colnames(xy) <- c('X', 'Y')
      # coordinates(xy) <- c("X", "Y")
      # proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
      # trnsf_xy <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
      # trnsf_xy <- st_transform(xy, crs=CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep=''))) # replace with this?
      # paste0(trnsf_xy$X, ", ", trnsf_xy$Y)
      
    # })
  })

  
  
  # output$info2 <- renderText({
  #   #click <- input$outputMap_click
  #   #if(is.null(click)){
  #   #  return()
  #   #}
  # 
  #   #paste0(click$lng)
  # })
  # WORKS
  output$info2 <- renderText({
    click <- input$setupMap_click
    if(is.null(click)){
     return()
    }

    # paste0(click$lng)
    
    ### try to get utm coord into text box after click
    long2 <- click$lng
    lat2 <- click$lat
    zone = long2UTM(long2)
    
    # xy <- as.data.frame(cbind(long2, lat2))
    # # str(xy)
    # colnames(xy) <- c('X', 'Y')
    # coordinates(xy) <- c("X", "Y")
    xy <- st_point(c(long2, lat2))
    # proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    xy_sfc <- st_sfc(xy, crs=4326)
    # trnsf_xy <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," +ellps=WGS84 +datum=WGS84", sep='')))
    trnsf_xy <- st_transform(xy_sfc, crs=CRS(paste("+proj=utm +zone=",zone," +ellps=WGS84 +datum=WGS84", sep='')))
    # paste0(trnsf_xy$X, ", ", trnsf_xy$Y)
    
    # output$setupMap <- renderLeaflet({
    #   # m_setup <- leaflet() %>% addMarkers(lng = click$lng, lat = click$lat) ## does not work
    #   m_setup %>% addTiles() %>% addMarkers(lng = click$lng, lat = click$lat)
    # })
    
    # output$setupMap$showPopup(click$lat, click$lng, "hello")
    
  })
  
  # observe({ ## does not work
  #   click <- input$setupMap_click
  #   if(is.null(click)){
  #     return()
  #   }
  #   
  #   (output$setupMap)$showPopup(click$lat, click$lng, "hello")
  #   
  # })
  
  
  observeEvent(input$plotpts, {
    
  })
  
  ### CALCULATIONS -------------------------------------------------------------------
  observeEvent(input$Calculate, {
    
    # print("calculate")
    
    # outer <- c(0,0)
    boundary <- c(0,0)
    
    #### test Whitehorse
    
    
    numwells <- 0
    
    if(abs(input$WellQ1) > 0){
      numwells <- numwells + 1
      wells_temp <- c(input$Wellx1, input$Welly1, input$WellQ1)
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
    
    # print(wells_temp)
    
    wells_temp[,1] <- (wells_temp[,1] - input$originx_setup) / input$distance1 # note: also need to account for rotation!
    wells_temp[,2] <- (wells_temp[,2] - input$originy_setup) / input$distance1
    
    # print(wells_temp)
    
    rot_matrix <- rbind(c(cos(- input$rotate1 * pi/180), - sin(- input$rotate1 * pi / 180)), c(sin(- input$rotate1 * pi / 180), cos(- input$rotate1 * pi / 180)))
    wells_temp2 <- matrix(0,numwells,3)
    for(i in 1:numwells){
      wells_temp2[i,1] <- rot_matrix[1,1] * wells_temp[i,1] + rot_matrix[1,2] * wells_temp[i,2]
      wells_temp2[i,2] <- rot_matrix[2,1] * wells_temp[i,1] + rot_matrix[2,2] * wells_temp[i,2]
    }
    
    wells_temp2[,3] <- wells_temp[,3]
    
    ### Determine case 2 option (based on Nagheli et al. examples)
    # option <- input$radioOpt_setup
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
    
    r_wedge <- input$distance1 # 1675.19
    
    K <- input$hyd_cond1 # 0.000243
    b <- 23
    
    q0 <- input$hyd_grad1 * K * b # regional uniform flow rate per unit width [L2/T]
    
    alpha <- input$alpha_setup * pi / 180
    beta <- input$beta1 * pi / 180 # pi - (2.75 - pi/2)
    beta <- beta - input$rotate1 * pi/180
    
    if(input$radioCase_setup == 1){
      consts <- cbind(c(r_wedge, beta, K, b, q0))
    }else if(input$radioCase_setup == 2){
      consts <- cbind(c(r_wedge, alpha, beta, K, b, q0))
    }
    # print(wells_temp2)
    # print(consts)
    
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
    
    # print(xy_stg)
    
    rot_angle <- input$rotate1 * pi/180
    rot_matrix <- rbind(c(cos(rot_angle), - sin(rot_angle)), c(sin(rot_angle), cos(rot_angle)))
    
    if(input$radioCase_setup == 2 && option == 4){
      zD_xy2 <- matrix(0.0, (2*numpts + 1)*(2*numpts + 1),2)
    }else{
      # zD_xy2 <- matrix(0.0, numpts*numpts,2)
      zD_xy2 <- matrix(0.0, (numpts + 1)*(numpts + 1),2)
    }
    
    # print(str(zD_xy))
    # print(nrow(zD_xy2))
    
    counter <- 1
    
    for(i in 1:nrow(zD_xy)){
      for(ii in 1:nrow(zD_xy)){ # there are only two columns, so use the number of rows again to make a square matrix
        zD_xy2[counter,1] <- rot_matrix[1,1] * zD_xy[ii,1] + rot_matrix[1,2] * zD_xy[i,2]
        zD_xy2[counter,2] <- rot_matrix[2,1] * zD_xy[ii,1] + rot_matrix[2,2] * zD_xy[i,2]
        counter <- counter + 1
      }
    }
    
    # print(counter)

    # ### scale the coordinates and translate about the UTM coordinates corresponding to the origin
    zD_xy2[,1] <- zD_xy2[,1] * r_wedge + input$originx_setup
    zD_xy2[,2] <- zD_xy2[,2] * r_wedge + input$originy_setup
    
    
    
    # if(input$radioCase_setup == 1 && nrow(wells_temp2) > 1){
    if(input$radioCase_setup == 1){
      ptA <- c(0, 1)
      ptB <- c(0, 0)
      ptC <- c(1, 0)
      ptD <- c(1, 1)
      # boundary <- rbind(ptA, ptB, ptC, ptA)
      boundary <- rbind(ptA, ptB, ptC, ptD)
    }else if(input$radioCase_setup == 2){
      ptA <- c(cos(alpha), sin(alpha))
      ptB <- c(0, 0)
      ptC <- c(1, 0)
      boundary <- rbind(ptA, ptB, ptC)
    }
    
    if(input$radioCase_setup == 2){
      outerarc <- cbind(cos(seq(0,alpha,alpha/10)), sin(seq(0,alpha,alpha/10)))
    }
    
    xlimits <- c(-1,1)
    ylimits <- c(0,2)
    
    xlimits2 <- c(-2,2)
    ylimits2 <- c(-2,2)
    
    # print(ylimits2)

    # ### Plot the data for the selected Case -----------------------
    
    extwells <- subset(wells_temp2, wells_temp2[,3]>0)
    injwells <- subset(wells_temp2, wells_temp2[,3]<0)
    
    ### remove the branch cuts
    
    str_holz <- holzbecher_v3(psi, zD_xy, wells_temp2, consts, xy_stg)
    psi_list <- str_holz[[1]]
    capzones <- str_holz[[2]]
    
    # print(length(psi_list))
    
    ### select equipotential lines
    
    numContoursDis <- 20
    phi_list <- contourLines(zD_xy[,1], zD_xy[,2], phi, nlevels=numContoursDis)
    
    ### trim the points to the aquifer boundaries
    psi_list <- trimAtBoundaries(psi_list, input$radioCase_setup, option, alpha, TRUE) # labelled x and y vectors
    capzones <- trimAtBoundaries(capzones, input$radioCase_setup, option, alpha, FALSE) # no labels on the x and y vectors
    phi_list <- trimAtBoundaries(phi_list, input$radioCase_setup, option, alpha, TRUE) # labelled x and y vectors
    
    # print(length(capzones))
    
    ### plot dimensionless data in RStudio window
    
    # numContoursStr <- 20
    # numContoursDis <- 20
    
    plot(boundary[,1], boundary[,2], lwd = 3, type = "l", col = "black", xlim=xlimits, ylim=ylimits, xlab="x", ylab="y")
    
    # contour(zD_xy[,1], zD_xy[,2], phi, nlevels=numContoursDis, col="blue", xlim=xlimits, ylim=ylimits)
    
    for(i in 1:length(phi_list)){ # plot the streamlines
      lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue")
    }
    
    for(i in 1:length(psi_list)){ # plot the streamlines
      lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red")
    }
    
    for(i in 1:length(capzones)){ # plot the capture zone boundaries
      lines(capzones[[i]][,1], capzones[[i]][,2], col="black", lwd=2)
    }
    
    lines(boundary[,1], boundary[,2], lwd = 3, col = "black")
    lines(outerarc[,1],outerarc[,2], lwd = 3, lty = 2, col = "black")
    points(extwells[,1],extwells[,2])
    if(nrow(injwells) > 0){
      points(injwells[,1], injwells[,2],pch=5)
    }
    points(xy_stg[,1], xy_stg[,2], pch=18, col="green")
    
    
    # ### rotate all points ---------------------
    boundary2 <- boundary
    
    for (i in 1:nrow(boundary)){
      boundary2[i,1] <- rot_matrix[1,1] * boundary[i,1] + rot_matrix[1,2] * boundary[i,2]
      boundary2[i,2] <- rot_matrix[2,1] * boundary[i,1] + rot_matrix[2,2] * boundary[i,2]
    }
    
    if(input$radioCase_setup == 2){
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
    
    capzones2 <- capzones
    
    for(i in 1:length(capzones2)){
      capzones2[[i]][,1] <- rot_matrix[1,1] * capzones[[i]][,1] + rot_matrix[1,2] * capzones[[i]][,2]
      capzones2[[i]][,2] <- rot_matrix[2,1] * capzones[[i]][,1] + rot_matrix[2,2] * capzones[[i]][,2]
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
    # print(bnds)
    
    # xycen <- as.data.frame(cbind(input$outputMap_center$lng, input$outputMap_center$lat))
    xycen <- st_point(c(input$outputMap_center$lng, input$outputMap_center$lat))
    # colnames(xycen) <- c('X', 'Y')
    # coordinates(xycen) <- c("X", "Y")
    # proj4string(xycen) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    xycen_sfc <- st_sfc(xycen, crs=4326)
    # trnsf_xycen <- spTransform(xycen, proj)
    trnsf_xycen <- st_transform(xycen_sfc, crs=proj)
    
    # print(st_coordinates(trnsf_xycen)[1])
    # print(st_coordinates(trnsf_xycen)[2])
    
    # xy1 <- as.data.frame(cbind(bnds$west, bnds$south))
    # colnames(xy1) <- c('X', 'Y')
    # coordinates(xy1) <- c("X", "Y")
    xy1 <- st_point(c(bnds$west, bnds$south))
    # proj4string(xy1) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    # trnsf_xy1 <- spTransform(xy1, proj)
    # print(paste0(trnsf_xy1$X, ", ", trnsf_xy1$Y))
    xy1_sfc <- st_sfc(xy1, crs=4326)
    trnsf_xy1 <- st_transform(xy1_sfc, crs=proj)
    
    # xy2 <- as.data.frame(cbind(bnds$east, bnds$north))
    # colnames(xy2) <- c('X', 'Y')
    # coordinates(xy2) <- c("X", "Y")
    xy2 <- st_point(c(bnds$east, bnds$north))
    # proj4string(xy2) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    # trnsf_xy2 <- spTransform(xy2, proj)
    # print(paste0(trnsf_xy2$X, ", ", trnsf_xy2$Y))
    xy2_sfc <- st_sfc(xy2, crs=4326)
    trnsf_xy2 <- st_transform(xy2_sfc, crs=proj)
    
    # xdist <- trnsf_xy2$X - trnsf_xy1$X
    # ydist <- trnsf_xy2$Y - trnsf_xy1$Y
    xdist <- st_coordinates(trnsf_xy2)[1] - st_coordinates(trnsf_xy1)[1]
    ydist <- st_coordinates(trnsf_xy2)[2] - st_coordinates(trnsf_xy1)[2]
    
    # print(xdist)
    # print(ydist)
    
    # xlimits2 <- xlimits2 * r_wedge + input$originx_setup
    # ylimits2 <- ylimits2 * r_wedge + input$originy_setup
    # xlimits2 <- c(trnsf_xycen$X - floor(xdist/2), trnsf_xycen$X + floor(xdist/2))
    # ylimits2 <- c(trnsf_xycen$Y - floor(ydist/2), trnsf_xycen$Y + floor(ydist/2))
    xlimits2 <- c(st_coordinates(trnsf_xycen)[1] - floor(xdist/2), st_coordinates(trnsf_xycen)[1] + floor(xdist/2))
    ylimits2 <- c(st_coordinates(trnsf_xycen)[2] - floor(ydist/2), st_coordinates(trnsf_xycen)[2] + floor(ydist/2))
    
    ### 
    
    boundary2[,1] <- boundary2[,1] * r_wedge + input$originx_setup
    boundary2[,2] <- boundary2[,2] * r_wedge + input$originy_setup
    
    if(input$radioCase_setup == 2){
      outerarc2[,1] <- outerarc2[,1] * r_wedge + input$originx_setup
      outerarc2[,2] <- outerarc2[,2] * r_wedge + input$originy_setup
    }
    
    extwells2[,1] <- extwells2[,1] * r_wedge + input$originx_setup
    extwells2[,2] <- extwells2[,2] * r_wedge + input$originy_setup
    
    if(nrow(injwells) > 0){
      injwells2[,1] <- injwells2[,1] * r_wedge + input$originx_setup
      injwells2[,2] <- injwells2[,2] * r_wedge + input$originy_setup
    }
    
    xy_stg2[,1] <- xy_stg2[,1] * r_wedge + input$originx_setup
    xy_stg2[,2] <- xy_stg2[,2] * r_wedge + input$originy_setup
    
    for(i in 1:length(psi_list2)){
      psi_list2[[i]]$x <- psi_list2[[i]]$x * r_wedge + input$originx_setup
      psi_list2[[i]]$y <- psi_list2[[i]]$y * r_wedge + input$originy_setup
    }
    
    for(i in 1:length(capzones2)){
      capzones2[[i]][,1] <- capzones2[[i]][,1] * r_wedge + input$originx_setup
      capzones2[[i]][,2] <- capzones2[[i]][,2] * r_wedge + input$originy_setup
    }
    
    for(i in 1:length(phi_list2)){
      phi_list2[[i]]$x <- phi_list2[[i]]$x * r_wedge + input$originx_setup
      phi_list2[[i]]$y <- phi_list2[[i]]$y * r_wedge + input$originy_setup
    }
    
    output$visplot1 <- renderPlot({
      plot(boundary2[,1], boundary2[,2], lwd = 3, type = "l", col = "black",xlim=xlimits2, ylim=ylimits2,xlab="x",ylab="y")
      
      for(i in 1:length(phi_list2)){ # plot the streamlines
        lines(phi_list2[[i]]$x, phi_list2[[i]]$y, col="blue")
      }
      
      for(i in 1:length(psi_list2)){ # plot the streamlines
        lines(psi_list2[[i]]$x, psi_list2[[i]]$y, col="red")
      }
      
      for(i in 1:length(capzones2)){ # plot the capture zone boundaries
        lines(capzones2[[i]][,1], capzones2[[i]][,2], col="black", lwd=2)
      }
      
      lines(boundary2[,1], boundary2[,2], lwd = 3, col = "black")
      lines(outerarc2[,1],outerarc2[,2], lwd = 3, lty = 2, col = "black")
      points(extwells2[,1],extwells2[,2])
      if(nrow(injwells) > 0){
        points(injwells2[,1], injwells2[,2],pch=5)
      }
      points(xy_stg2[,1], xy_stg2[,2], pch=18, col="green")
      
    })
    
    # if(input$radioSHP_vis == 1){ ### export shapefiles
    #   # ### Save data as shapefiles -----------------------------------------
    #   z1 <- long2UTM(input$outputMap_center$lng)
    #   #utm.prj <- "+proj=utm +zone=z1 ellps=WGS84"
    #   # proj <- CRS("+proj=utm +zone=8 +datum=WGS84") # Carmacks is in UTM Zone 8N
    #   proj <- CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))
    #   
    #   ### Need lists of (x,y,phi) and (x,y,psi) points, not just the phi and psi matrices
    #   
    #   ### export equipotential lines
    #   counter <- 1
    #   xyphi <- matrix(data=0.0,numpts*numpts,4)
    #   
    #   for(i in 1:length(phi_list2)){
    #     for(ii in 1:length(phi_list2[[i]]$x)){
    #       if(is.na(phi_list2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
    #         xyphi[counter,1] <- phi_list2[[i]]$x[ii]
    #         xyphi[counter,2] <- phi_list2[[i]]$y[ii]
    #         xyphi[counter,3] <- i # group number
    #         xyphi[counter,4] <- ii # order number
    #         counter <- counter + 1
    #       }
    #     }
    #   }
    #   
    #   xyphi <- xyphi[c(1:counter-1), c(1,2,3,4),drop=FALSE] # remove zero rows
    #   
    #   # print(counter)
    #   # write.csv(xyphi, file = "xyphi.txt")
    #   
    #   equipotpts <- SpatialPoints(xyphi, proj4string = proj)
    #   shapefile(equipotpts, filename = "equipotential.shp", overwrite = TRUE)
    #   # ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
    #   dbfdata <- read.dbf("equipotential.dbf", as.is = TRUE)
    #   # ## add new attribute data (just the numbers 1 to the number of objects)
    #   dbfdata$group <- xyphi[,3] ### creates a new field called "group"
    #   dbfdata$order <- xyphi[,4] ### creates a new field called "order"
    #   write.dbf(dbfdata, "equipotential.dbf") ## overwrite the file with this new copy
    #   
    #   ### export flowlines
    #   counter <- 1
    #   xypsi <- matrix(data=0.0,numpts*numpts,4)
    #   
    #   for(i in 1:length(psi_list2)){
    #     for(ii in 1:length(psi_list2[[i]]$x)){
    #       if(is.na(psi_list2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
    #         xypsi[counter,1] <- psi_list2[[i]]$x[ii]
    #         xypsi[counter,2] <- psi_list2[[i]]$y[ii]
    #         xypsi[counter,3] <- i # group number
    #         xypsi[counter,4] <- ii # order number
    #         counter <- counter + 1
    #       }
    #     }
    #   }
    #   
    #   xypsi <- xypsi[c(1:counter-1), c(1,2,3,4),drop=FALSE] # remove zero rows
    #   
    #   flowlinepts <- SpatialPoints(xypsi, proj4string = proj)
    #   shapefile(flowlinepts, filename = "flowlines.shp", overwrite = TRUE)
    #   ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
    #   dbfdata <- read.dbf("flowlines.dbf", as.is = TRUE)
    #   # ## add new attribute data (just the numbers 1 to the number of objects)
    #   dbfdata$group <- xypsi[,3] ### creates a new field called "group"
    #   dbfdata$order <- xypsi[,4] ### creates a new field called "order"
    #   write.dbf(dbfdata, "flowlines.dbf") ## overwrite the file with this new copy
    #   
    #   ### export capture zone boundaries
    #   
    #   counter <- 1
    #   xycap <- matrix(data=0.0,numpts*numpts,4)
    #   
    #   for(i in 1:length(capzones2)){
    #     for(ii in 1:length(capzones2[[i]][,1])){
    #       if(is.na(capzones2[[i]][ii,1]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
    #         xycap[counter,1] <- capzones2[[i]][ii,1]
    #         xycap[counter,2] <- capzones2[[i]][ii,2]
    #         xycap[counter,3] <- i # group number
    #         xycap[counter,4] <- ii # order number
    #         counter <- counter + 1
    #       }
    #     }
    #   }
    #   
    #   xycap <- xycap[c(1:counter-1), c(1,2,3,4),drop=FALSE] # remove zero rows
    #   
    #   cappts <- SpatialPoints(xycap, proj4string = proj)
    #   shapefile(cappts, filename = "capzones.shp", overwrite = TRUE)
    #   ### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
    #   dbfdata <- read.dbf("capzones.dbf", as.is = TRUE)
    #   # ## add new attribute data (just the numbers 1 to the number of objects)
    #   dbfdata$group <- xycap[,3] ### creates a new field called "group"
    #   dbfdata$order <- xycap[,4] ### creates a new field called "order"
    #   write.dbf(dbfdata, "capzones.dbf") ## overwrite the file with this new copy
    #   
    #   ### export stagnation points
    #   
    #   stgpts <- SpatialPoints(xy_stg2, proj4string = proj)
    #   shapefile(stgpts, filename = "stagnationpts.shp", overwrite = TRUE)
    #   
    #   ### Note: the "Points to path" tool in QGIS can be used to generate lines from the sets of points using the "group" and "order" attributes
    #   
    #    
    # }
    
    
    #   arrowLen <- 0.2
    #   arrows(0, 0, arrowLen * cos(beta[option]), arrowLen * sin(beta[option]))
    
    
  })
  
}

shinyApp(ui, server)