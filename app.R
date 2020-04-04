### Code written by Pranav Minasandra from 30 Mar to 4 April
### in collaboration with Vishwesha Guttal.
### Let us know if you use this code.

################# Basic Stuff
library(shiny)
library(deSolve)

TOT_POP <- 6.56e07
PROP_YOUNG <- 0.9
BASAL_TRANSMISSION_RATE <- 0.3

recompute_init_vals <- function(){
init_vals <<- c(
  S_young = PROP_YOUNG*TOT_POP,#S_young    1
  E_young = 5,#E_young    2
  I_a_young =3,#I_a_young  3
  I_s_young =1,#I_s_young  4
  I_h_young =0,#I_h_young  5
  I_c_young =0,#I_c_young  6
  R_young =0,#R_young    7
  D_young =0,#D_young    8
  S_old = (1-PROP_YOUNG)*TOT_POP,#S_old    1
  E_old = 0,#E_old    2
  I_a_old =0,#I_a_old  3
  I_s_old =0,#I_s_old  4
  I_h_old =0,#I_h_old  5
  I_c_old =0,#I_c_old  6
  R_old =0,#R_old    7
  D_old =0,#D_old    8
  CumInf=5,
  CumInf_S=4
)}

recompute_init_vals()

####
#Empirical Values
Tau_I <- 5
Tau_a <- 7
Tau_s <- 7
Tau_h <- 4
Tau_c <- 10
P_a <- 0.25
P_s <- 0.55
P_h <- 0.15
P_c <- 0.05
f_c_old <- 0.35
f_c_young <- 0.01

params <- c(
  AsymSusRate =BASAL_TRANSMISSION_RATE,
	SymSusRate =2*BASAL_TRANSMISSION_RATE,
	HospitalSusRate =0.001,
	SevereSusRate =0.001,
	YoungFracSym =1 - P_a,
	YoungExposedInfectedRate =1/Tau_I,#0.2
	YoungAsymRecoveryRate =1/Tau_a,
	YoungSymRecoveryRate =P_s/Tau_s,#0.0785714286
	YoungHospitalRecoveryRate =1/Tau_h *P_h/(P_h + P_c),
	YoungSevereRecoveryRate =1/Tau_c * (1-f_c_young)/P_c,
	YoungSymHospitalRate =(1-P_s)/Tau_s ,
	YoungHospitalSevereRate =1/Tau_h * P_c/(P_h + P_c),
	YoungSevereMortRate =1/Tau_c * (f_c_young)/P_c,
	OldFracSym =1 - P_a,
	OldExposedInfectedRate =2*Tau_I,#0.4
	OldAsymRecoveryRate =0.5*Tau_a,
	OldSymRecoveryRate =0.5*P_s/Tau_s,#0.5*0.0785714286
	OldHospitalRecoveryRate =0.5*1/Tau_h *P_h/(P_h + P_c),
	OldSevereRecoveryRate =0.5*1/Tau_c * (1-f_c_old)/P_c,
	OldSymHospitalRate =2*(1-P_s)/Tau_s,
	OldHospitalSevereRate =2*1/Tau_h * P_c/(P_h + P_c),
	OldSevereMortRate =1/Tau_c * (f_c_old)/P_c
)

### Empirically observed data (defunct largely)
kar_data <- c(7,8,11,14,15,15,20,26,
          33,41,51,55,64,76,83,88,101,110)
mah_data <- c(33,39,41,45,48,52,64,74,97,107,
              122,130,153,186,203,220,302,335)
tn_data <- c(1,1,1,1,3,3,6,7,9,15,
             23,29,38,42,50,67,124,234)
in_data <- c(97,107,118,137,151,173,223,283,
             360,434,519, 606, 694, 834, 918, 1024,1251,1397,1834,2069)

kar_ori <- as.Date('150320', format='%d%m%y')
mah_ori <- as.Date('150320', format='%d%m%y')
in_ori <- as.Date('150320', format='%d%m%y')
tn_ori <- as.Date('150320', format='%d%m%y')
kl_ori <- tn_ori

### Derivative functions
tot_pop <- TOT_POP #Remnant of old code. Keep intact.
none_der_statevector <- function(t, vars, params){
  with(as.list(c(vars,params)),{
  return(list(c(
	#1. Derivative of S_young
	-(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_young,

	#2. Derivative of E_young
	(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_young
	- YoungExposedInfectedRate*E_young,

	#3. Derivative of I_a_young
	(1 - YoungFracSym) * YoungExposedInfectedRate * E_young
	- YoungAsymRecoveryRate * I_a_young,

	#4. Derivative of I_s_young
	YoungFracSym * YoungExposedInfectedRate * E_young
	- YoungSymRecoveryRate * I_s_young
	- YoungSymHospitalRate * I_s_young,

	#5. Derivative of I_h_young
	YoungSymHospitalRate * I_s_young
	- YoungHospitalRecoveryRate * I_h_young
	- YoungHospitalSevereRate * I_h_young,

	#6. Derivative of I_c_young
	YoungHospitalSevereRate * I_h_young
	- YoungSevereRecoveryRate * I_c_young
	- YoungSevereMortRate * I_c_young,

	#7. Derivative of R_young
	YoungAsymRecoveryRate * I_a_young
	+ YoungSymRecoveryRate * I_s_young
	+ YoungHospitalRecoveryRate * I_h_young
	+ YoungSevereRecoveryRate * I_c_young,

	#8. Derivative of D_young
	YoungSevereMortRate * I_c_young,

	#9. Derivative of S_old
	-(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_old,

	#10. Derivative of E_young
	(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_old
	- OldExposedInfectedRate*E_old,

	#11. Derivative of I_a_old
	(1 - OldFracSym) * OldExposedInfectedRate * E_old
	- OldAsymRecoveryRate * I_a_old,

	#12. Derivative of I_s_old
	OldFracSym * OldExposedInfectedRate * E_old
	- OldSymRecoveryRate * I_s_old
	- OldSymHospitalRate * I_s_old,

	#13. Derivative of I_h_old
	OldSymHospitalRate * I_s_old
	- OldHospitalRecoveryRate * I_h_old
	- OldHospitalSevereRate * I_h_old,

	#14. Derivative of I_c_old
	OldHospitalSevereRate * I_h_old
	- OldSevereRecoveryRate * I_c_old
	- OldSevereMortRate * I_c_old,

	#15. Derivative of R_old
	OldAsymRecoveryRate * I_a_old
	+ OldSymRecoveryRate * I_s_old
	+ OldHospitalRecoveryRate * I_h_old
	+ OldSevereRecoveryRate * I_c_old,

	#16. Derivative of D_old
	OldSevereMortRate * I_c_old,

	#17. Derivative of Cumulative Infections	
	YoungExposedInfectedRate*E_young
	+ OldExposedInfectedRate*E_old,

	#18. Derivative of Cumulative Symptomatic Infections
	YoungExposedInfectedRate*E_young*YoungFracSym
	+ OldExposedInfectedRate*E_old*OldFracSym
      )))
  })
}

symquar_der_statevector <- function(t, vars, params){
  with(as.list(c(vars,params)),{
  return(list(c(
	#1. Derivative of S_young
	-(AsymSusRate*(I_a_young + I_a_old*(1-lockdown_effectiveness))
		+ SymSusRate*(I_s_young + I_s_old*(1-lockdown_effectiveness))
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_young,

	#2. Derivative of E_young
	(AsymSusRate*(I_a_young + I_a_old*(1-lockdown_effectiveness))
		+ SymSusRate*(I_s_young + I_s_old*(1-lockdown_effectiveness))
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_young
	- YoungExposedInfectedRate*E_young,

	#3. Derivative of I_a_young
	(1 - YoungFracSym) * YoungExposedInfectedRate * E_young
	- YoungAsymRecoveryRate * I_a_young,

	#4. Derivative of I_s_young
	YoungFracSym * YoungExposedInfectedRate * E_young
	- YoungSymRecoveryRate * I_s_young
	- YoungSymHospitalRate * I_s_young,

	#5. Derivative of I_h_young
	YoungSymHospitalRate * I_s_young
	- YoungHospitalRecoveryRate * I_h_young
	- YoungHospitalSevereRate * I_h_young,

	#6. Derivative of I_c_young
	YoungHospitalSevereRate * I_h_young
	- YoungSevereRecoveryRate * I_c_young
	- YoungSevereMortRate * I_c_young,

	#7. Derivative of R_young
	YoungAsymRecoveryRate * I_a_young
	+ YoungSymRecoveryRate * I_s_young
	+ YoungHospitalRecoveryRate * I_h_young
	+ YoungSevereRecoveryRate * I_c_young,

	#8. Derivative of D_young
	YoungSevereMortRate * I_c_young,

	#9. Derivative of S_old
	-(AsymSusRate*(I_a_young + I_a_old)*(1-lockdown_effectiveness)
		+ SymSusRate*(I_s_young + I_s_old)*(1-lockdown_effectiveness)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_old,

	#10. Derivative of E_young
	(AsymSusRate*(I_a_young + I_a_old)*(1-lockdown_effectiveness)
		+ SymSusRate*(I_s_young + I_s_old)*(1-lockdown_effectiveness)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_old
	- OldExposedInfectedRate*E_old,

	#11. Derivative of I_a_old
	(1 - OldFracSym) * OldExposedInfectedRate * E_old
	- OldAsymRecoveryRate * I_a_old,

	#12. Derivative of I_s_old
	OldFracSym * OldExposedInfectedRate * E_old
	- OldSymRecoveryRate * I_s_old
	- OldSymHospitalRate * I_s_old,

	#13. Derivative of I_h_old
	OldSymHospitalRate * I_s_old
	- OldHospitalRecoveryRate * I_h_old
	- OldHospitalSevereRate * I_h_old,

	#14. Derivative of I_c_old
	OldHospitalSevereRate * I_h_old
	- OldSevereRecoveryRate * I_c_old
	- OldSevereMortRate * I_c_old,

	#15. Derivative of R_old
	OldAsymRecoveryRate * I_a_old
	+ OldSymRecoveryRate * I_s_old
	+ OldHospitalRecoveryRate * I_h_old
	+ OldSevereRecoveryRate * I_c_old,

	#16. Derivative of D_old
	OldSevereMortRate * I_c_old,

	#17. Derivative of Cumulative Infections	
	YoungExposedInfectedRate*E_young
	+ OldExposedInfectedRate*E_old,

	#18. Derivative of Cumulative Symptomatic Infections
	YoungExposedInfectedRate*E_young*YoungFracSym
	+ OldExposedInfectedRate*E_old*OldFracSym

      )))
  })
}

lockdown_begin_date <- as.Date('240320', format='%d%m%y')
lockdown_duration <- 21
lockdown_effectiveness <- 0.75
lockdown <- function(t, begin, duration){
  t <- ori+t
	if(t<(begin+duration) & t>begin){
		return(1-lockdown_effectiveness)
	}
	else{return(1)}
}

lockdown_der_statevector <- function(t, vars, params){
  with(as.list(c(vars,params)),{
  return(list(c(
	#1. Derivative of S_young
	-(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_young *lockdown(t, lockdown_begin_date, lockdown_duration),

	#2. Derivative of E_young
	(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_young *lockdown(t, lockdown_begin_date, lockdown_duration)
	- YoungExposedInfectedRate*E_young,

	#3. Derivative of I_a_young
	(1 - YoungFracSym) * YoungExposedInfectedRate * E_young
	- YoungAsymRecoveryRate * I_a_young,

	#4. Derivative of I_s_young
	YoungFracSym * YoungExposedInfectedRate * E_young
	- YoungSymRecoveryRate * I_s_young
	- YoungSymHospitalRate * I_s_young,

	#5. Derivative of I_h_young
	YoungSymHospitalRate * I_s_young
	- YoungHospitalRecoveryRate * I_h_young
	- YoungHospitalSevereRate * I_h_young,

	#6. Derivative of I_c_young
	YoungHospitalSevereRate * I_h_young
	- YoungSevereRecoveryRate * I_c_young
	- YoungSevereMortRate * I_c_young,

	#7. Derivative of R_young
	YoungAsymRecoveryRate * I_a_young
	+ YoungSymRecoveryRate * I_s_young
	+ YoungHospitalRecoveryRate * I_h_young
	+ YoungSevereRecoveryRate * I_c_young,

	#8. Derivative of D_young
	YoungSevereMortRate * I_c_young,

	#9. Derivative of S_old
	-(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_old*lockdown(t, lockdown_begin_date, lockdown_duration),

	#10. Derivative of E_young
	(AsymSusRate*(I_a_young + I_a_old)
		+ SymSusRate*(I_s_young + I_s_old)
		+ HospitalSusRate*(I_h_young + I_h_old)
		+ SevereSusRate*(I_c_young + I_c_old))/tot_pop * S_old*lockdown(t, lockdown_begin_date, lockdown_duration)
	- OldExposedInfectedRate*E_old,

	#11. Derivative of I_a_old
	(1 - OldFracSym) * OldExposedInfectedRate * E_old
	- OldAsymRecoveryRate * I_a_old,

	#12. Derivative of I_s_old
	OldFracSym * OldExposedInfectedRate * E_old
	- OldSymRecoveryRate * I_s_old
	- OldSymHospitalRate * I_s_old,

	#13. Derivative of I_h_old
	OldSymHospitalRate * I_s_old
	- OldHospitalRecoveryRate * I_h_old
	- OldHospitalSevereRate * I_h_old,

	#14. Derivative of I_c_old
	OldHospitalSevereRate * I_h_old
	- OldSevereRecoveryRate * I_c_old
	- OldSevereMortRate * I_c_old,

	#15. Derivative of R_old
	OldAsymRecoveryRate * I_a_old
	+ OldSymRecoveryRate * I_s_old
	+ OldHospitalRecoveryRate * I_h_old
	+ OldSevereRecoveryRate * I_c_old,

	#16. Derivative of D_old
	OldSevereMortRate * I_c_old,

	#17. Derivative of Cumulative Infections	
	YoungExposedInfectedRate*E_young
	+ OldExposedInfectedRate*E_old,

	#18. Derivative of Cumulative Symptomatic Infections
	YoungExposedInfectedRate*E_young*YoungFracSym
	+ OldExposedInfectedRate*E_old*OldFracSym

      )))
  })
}

### UI
ui <- fluidPage(
  theme = shinythemes::shinytheme("darkly"),
  titlePanel("IISc COVID-19 Age-Structured SEIR Model"),
  
  sidebarLayout(
    
    sidebarPanel(
      theme = shinythemes::shinytheme("darkly"),
      #br(),
      radioButtons(
        "quar",
        "Intervention measures",
        c("Full-scale lockdown"="lockdown", 
          "Senior citizens home-quarantined"="symquar", 
          "None"="none")
      ),
      sliderInput("dur", "Days plotted",
                  min=Sys.Date()-5,max=as.Date('150520', format='%d%m%y'),step=5,value=as.Date('200420',format='%d%m%y')),
      sliderInput("lockdowndur", "Days of lockdown",
                  min=0,max=50,value = c(lockdown_duration)),
      sliderInput("lockdowneff", "Intervention Effectiveness", min=0.0, max=1.0, value = lockdown_effectiveness),
      helpText("By what proportion are transmission rates reduced by staying home only?"),
      br(),
      selectInput(
        "state",
        "State",
        c("All India"="in",
          "Karnataka"="kar",
          "Maharashtra"="mah",
          "Tamil Nadu"="tn"
          #"Kerala"="kl",
        )
      ),br(),
      
      imageOutput("dispImage")
      ),
    
      mainPanel(
        h3(textOutput("captionr0", container = span)),br(),
        tabsetPanel(
          type="tabs",
		#tabPanel("Susceptible", plotOutput("plot_s")),
		#tabPanel("Exposed", plotOutput("plot_e")),
		tabPanel("Symptomatic Infections", plotOutput("plot_cis")),
          #tabPanel("Asymptomatic Infections", plotOutput("plot_ia")),
          #tabPanel("Symptomatic Infections", plotOutput("plot_is")),
          tabPanel("Hospitalisations", plotOutput("plot_ih")),
          #tabPanel("ICU Patients", plotOutput("plot_ic")),
          #tabPanel("Recoveries", plotOutput("plot_r")),
          #tabPanel("Deaths", plotOutput("plot_d")),
		#tabPanel("Documented Cases", plotOutput("plot_em")),
		tabPanel("Caveats",tags$iframe(style="height:600px; width:100%", src="Caveats.pdf")),
		tabPanel("FAQs & Report", tags$iframe(style="height:600px; width:100%", src="FAQ.pdf")),
		tabPanel("About", 
		         h2("The authors"),
		         p("This work was completed by",
		         a("Vishwesha Guttal", href = "http://teelabiisc.wordpress.com/"),
		         "and",
		         a("Pranav Minasandra", href="https://pminasandra.weebly.com/"),
		         "from the Indian Institute of Science, Bangalore, India."
		         ),
		         br(),
		         h2("Short summary"),
		         p("In this small age-stuctured SEIR model, we consider the states of 
		           Karnataka, Tamil Nadu, and Maharashtra, as well as India itself, 
		           and try to arrive at quantitative interpretations of
		           the epidemiology of COVID-19."),
		         br(),
		         p(strong("This is only an illustrative model, and should not be used for 
		           forecasting.")),br(),
		         h2("Inspirations"),
		         p("There are several different models that served as inspiration for this one,
		           including",
		           a("Alison Hill's model", href="https://alhill.shinyapps.io/COVID19seir/"),
		           "and", a("Joshua Weitz's model.", href="https://github.com/jsweitz/covid-19-ga-summer-2020")),
		         br(),
		         h2("Code availability"),
		         p("All code used in this app is available on",
		           a("github", href="https://github.com/tee-lab/AS-SEIR-COVID19-India/"))
		         )
          #tabPanel("Cumulative Infections", plotOutput("plot_ci")),
          
        )
      
    )
  )
)
#Quarantining <- FALSE
#sym_quar <- function(t, quarantining){
#  if(quarantining==TRUE){
#    if(t>20){
#      return(0)
#    }
#  }
#  else{return(1)}
#}


### Server
server <- function(input, output){

d <- reactive({
	f <- switch(
		input$quar,
		none = none_der_statevector,
		symquar = symquar_der_statevector,
		lockdown = lockdown_der_statevector
	)
	if(input$state=="kar"){
	  TOT_POP <<- 6.56e07
	  tot_pop <<- TOT_POP
	  PROP_YOUNG <<- 0.9
	  ori <<- kar_ori
	  e_data <<- kar_data
	  recompute_init_vals()
	  init_vals[2:4] <- init_vals[2:4]*e_data[1]
	}
	
	else if(input$state=="mah"){
	  TOT_POP <<- 1.19e08
	  tot_pop <<- TOT_POP
	  PROP_YOUNG <<- 0.89
	  ori <<- mah_ori
	  e_data <<- mah_data
	  recompute_init_vals()
	  init_vals[2:4] <- init_vals[2:4]*e_data[1]
	}
	
	else if(input$state=="in"){
	  TOT_POP <<- 1.387e09
	  tot_pop <<- TOT_POP
	  PROP_YOUNG <<- 0.911
	  ori <<- in_ori
	  e_data <<- in_data
	  recompute_init_vals()
	  init_vals[2:4] <- init_vals[2:4]*e_data[1]
	}
	
	else if(input$state=="tn"){
	  TOT_POP <<- 7.52e07
	  PROP_YOUNG <<- 0.88
	  ori <<- tn_ori
	  e_data <<- tn_data
	  recompute_init_vals()
	  init_vals[2:4] <- init_vals[2:4]*e_data[1]
	}
	
	params <- c(
	  AsymSusRate =BASAL_TRANSMISSION_RATE,
	  SymSusRate =2*BASAL_TRANSMISSION_RATE,
	  HospitalSusRate =0.001,
	  SevereSusRate =0.001,
	  YoungFracSym =1 - P_a,
	  YoungExposedInfectedRate =1/Tau_I,#0.2
	  YoungAsymRecoveryRate =1/Tau_a,
	  YoungSymRecoveryRate =P_s/Tau_s,#0.0785714286
	  YoungHospitalRecoveryRate =1/Tau_h *P_h/(P_h + P_c),
	  YoungSevereRecoveryRate =1/Tau_c * (1-f_c_young)/P_c,
	  YoungSymHospitalRate =(1-P_s)/Tau_s ,
	  YoungHospitalSevereRate =1/Tau_h * P_c/(P_h + P_c),
	  YoungSevereMortRate =1/Tau_c * (f_c_young)/P_c,
	  OldFracSym =1 - P_a,
	  OldExposedInfectedRate =2*Tau_I,#0.4
	  OldAsymRecoveryRate =0.5*Tau_a,
	  OldSymRecoveryRate =0.5*P_s/Tau_s,#0.5*0.0785714286
	  OldHospitalRecoveryRate =0.5*1/Tau_h *P_h/(P_h + P_c),
	  OldSevereRecoveryRate =0.5*1/Tau_c * (1-f_c_old)/P_c,
	  OldSymHospitalRate =2*(1-P_s)/Tau_s,
	  OldHospitalSevereRate =2*1/Tau_h * P_c/(P_h + P_c),
	  OldSevereMortRate =1/Tau_c * (f_c_old)/P_c
	)
	tot_pop <<- TOT_POP
	t <- seq(0, as.double(input$dur - ori), by=1)
	lockdown_effectiveness <<- input$lockdowneff
	lockdown_duration <<- input$lockdowndur
	data.frame(ode(y=init_vals, func=f, times=t, parms=params))
  
})

output$dispImage <- renderImage({
 filename <- paste("./www/",
                   input$state, ".gif", sep="") 
 
 list(src=filename, alt="Map of State",width='200')
}, deleteFile = F)


output$plot_ih <- renderPlot({
  out <<- d()
  dates <- ori + out$time
  plot(out$I_h_young + out$I_h_old ~ dates, 
       col="blue", type='l', 
       xlab="", ylab="Population size", 
       main=" Hospitalised Individuals",
       xaxt="n",
       lwd=1.5)
  if(input$quar=="lockdown"){
    abline(v=lockdown_begin_date)
    abline(v=lockdown_begin_date + input$lockdowndur)
    
  }
  if(input$quar != "none"){
    t <- seq(0, as.double(input$dur - ori), by=1)
    recompute_init_vals()
    init_vals[2:4] <- init_vals[2:4]*e_data[1]
    out2 <- data.frame(ode(y=init_vals, func=none_der_statevector, times=t, parms=params))
    points(out2$I_h_young + out2$I_h_old ~ dates, type='l', col='blue', lty='dotted')
  }
  axis(1, dates, format(dates, '%d/%m'), cex.axis=0.6)
})


output$plot_cis <- renderPlot({
  out <<- d()
  dates <- ori + out$time
  plot(out$CumInf_S ~ dates, 
       col="red", type='l', 
       xlab="", ylab="Population size", 
       main=" Cumulative Symptomatic Individuals",
       xaxt="n",
       lwd=1.5)
  axis(1, dates, format(dates, '%d/%m'), cex.axis=0.6)
  if(input$quar=="lockdown"){
    abline(v=lockdown_begin_date)
    abline(v=lockdown_begin_date + input$lockdowndur)
  }
  if(input$quar != "none"){
    t <- seq(0, as.double(input$dur - ori), by=1)
    recompute_init_vals()
    init_vals[2:4] <- init_vals[2:4]*e_data[1]
    out2 <- data.frame(ode(y=init_vals, func=none_der_statevector, times=t, parms=params))
    points(out2$CumInf_S ~ dates, type='l', col='red', lty='dotted')
  }
})

output$captionr0 <- renderText({
  if(input$quar=="none"){
    paste("R_0 =", round(P_a*BASAL_TRANSMISSION_RATE*Tau_a + P_s*BASAL_TRANSMISSION_RATE*2*Tau_s,3))
  }
  else if(input$quar=="symquar"){
    paste("R_0 =", round(P_a*BASAL_TRANSMISSION_RATE*(Tau_a*(PROP_YOUNG*(input$lockdowneff) + 1 - input$lockdowneff)) + P_s*BASAL_TRANSMISSION_RATE*2*Tau_s*(PROP_YOUNG*(input$lockdowneff) + 1 - input$lockdowneff),3))
  }
  else if(input$quar=="lockdown"){
    paste("R_0 =", (1-input$lockdowneff)*round(P_a*BASAL_TRANSMISSION_RATE*Tau_a + P_s*BASAL_TRANSMISSION_RATE*2*Tau_s,3))
  }
})


#output$faqviewer <- renderText({
#  return(paste('<iframe style="height:600px; width:100%" src="', './', '"></iframe>', sep = ""))
#})
}

shinyApp(ui, server)
