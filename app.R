### Code written by Pranav Minasandra from 30 Mar to 4 April
### at IISc, in the middle of the coronavirus pandemic and lockdown,
### in collaboration with Vishwesha Guttal.

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

kar_data <- c(7,8,11,14,15,15,20,26,
          33,41,51,55,64,76,83,88,101,110)
mah_data <- c(33,39,41,45,48,52,64,74,97,107,
              122,130,153,186,203,220,302,335)
tn_data <- c(1,1,1,1,3,3,6,7,9,15,
             23,29,38,42,50,67,124,234)
in_data <- c(97,107,118,137,151,173,223,283,
             360,434,519, 606, 694, 834, 918, 1024,1251,1397,1834,2069)

tot_pop <- TOT_POP
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
ui <- fluidPage(
  theme = shinythemes::shinytheme("darkly"),
  titlePanel("COVID-19 Age-Structured SEIR Model"),
  
  sidebarLayout(
    
    sidebarPanel(
      theme = shinythemes::shinytheme("darkly"),
      #br(),
      radioButtons(
        "quar",
        "Intervention measures",
        c("Quarantine vulnerable"="symquar",
          "Lockdown"="lockdown", 
          "None"="none")
      ),
      sliderInput("dur", "Days plotted",
                  min=as.Date('150320', format='%d%m%y'),max=as.Date('150520', format='%d%m%y'),step=5,value=as.Date('200420',format='%d%m%y')),
      sliderInput("lockdowndur", "Days of lockdown (NA for quaratine vulnerable)",
                  min=0,max=50,value = c(lockdown_duration)),
      sliderInput("lockdowneff", "Intervention Effectiveness", min=0.0, max=1.0, value = lockdown_effectiveness),
      helpText("(By what proportion are transmission rates reduced by the intervention?)"),
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
        #h3(textOutput("captionr0", container = span)),br(),
        tabsetPanel(
          type="tabs",
		#tabPanel("Susceptible", plotOutput("plot_s")),
		#tabPanel("Exposed", plotOutput("plot_e")),
		tabPanel("Introduction", 
		         h2("Caveats"),
		         p("We begin with a caveat, which is really important during the times of an ongoing epidemic. 
		         
		         The model of covid-19 that we show is simple, ignores many features of real disease 
		         dynamics and demographic features of Indian population. 
		         
             Hence, it is to be used only for understanding concepts but NOT for 
		         making quantitative forecasts about the pandemic." ),
		         br(),
		         h2("A mathematical model"),
		         p("The purpose of this R-Shiny-App is to get a 'hands-on' feeling of how the 
		         epidemic grows. We can also try to understand how intervention measures (e.g. quarantine, lockdown) help reduce the growth of disease."),
		         br(),
		         p("To do this, we use a simple mathematical model - which are simplistic representations of 
		           real world scenario to make predictions and get useful insights."),
		         br(),
		         p("We know that novel coronavirus affects younger (usually less than 60 years) and older population
		           differently, with the older being much vulberable. So, we assume that population consists of two types of individuals: fit and vulnerable."),
		         br(),
		         p("We assume that all individuals are susceptible (S) to the disease. Over time, they may become exposed (E) to virus 
		           because they came in contact with a sick person. After a while these exposed individuals become infected (I) with virus."),
		         br(),
		         p("We now know that covid-19 can infect different people in different ways. Some never show symptoms and recover (R). 
		           Some recover with mild illness. However, some of them become quite sick and may need hospitalisation and ICU. Some of these may even die (D)."),
		         br(),
		         p("We incorporate all these features of the disease via mathematical equations. Such models are SEIR models. 
		           We then solve them using computers (see Report for technical details)."),
		         br(),
		         h2("Using the App"),
		         p("Using this app, you can see what happens to epidemic when the virus is left to itself, with no interventions. 
		           And how the situation can improve with interventions."),
		         br(),
		         p("To do this, go to Symptomatic infections, where you see a plot of total infected individuals (who also show symptoms) over time. 
		           You can chose three different 'Interevention measures' on the left panel."),
		         br(),
		         p("Quarantine the vulnerable means that vulnerable (elderly and those with other health issues) minimise
		           contact with rest of the population."),
		         br(),
		         p("Lockdown means that the entire population shows physical distancing, and minimizes interactions."),
		         br(),
		         p("Intervention effectiveness is a heuristic for how effectively you think people will follow the rules 
		           or advisory of the government."),
		         br(),
		         h2("Flattening the curve"),
		         p("By now you may hvae heard of flatten the curve. It just means that slow down the spread of disease."),
		         br(),
		         p("Chose 'Lockdown' as the intervention measure. You will see two curves in both Infections and Hospitalisation tabs.
		           You will see that the steepness of these curves is smaller when intervention is present"),
		          br(), 
		          p("Play with sliding bars on the left panel to change how effective the interventions are. Observe how it affects the curves")
		         ),
		tabPanel("Total Infections", plotOutput("plot_cis")),
          #tabPanel("Asymptomatic Infections", plotOutput("plot_ia")),
          #tabPanel("Infections", plotOutput("plot_is")),
    tabPanel("Hospitalisations", plotOutput("plot_ih")),
          #tabPanel("ICU Patients", plotOutput("plot_ic")),
          #tabPanel("Recoveries", plotOutput("plot_r")),
          #tabPanel("Deaths", plotOutput("plot_d")),
		#tabPanel("Documented Cases", plotOutput("plot_em")),
		#tabPanel("Caveats",tags$iframe(style="height:600px; width:100%", src="Caveats.pdf")),
		tabPanel("Technical Report", tags$iframe(style="height:600px; width:100%", src="age-seir-covid-report.pdf")),
		tabPanel("FAQ",
		         h2("FAQ"), 
		         h3("Disclaimer"),
		         p("This model is simple, ignores many features of real disease dynamics and demographic features of Indian
population. Hence, it is to be used only for understanding concepts but NOT for making quantitative
forecasts about the pandemic."),
		         h3("What is the use of this model?"),
		         p("This model can be used to understand how the epidemic spreads 
		         without and with interventions. Specically, we incorporate two types of 
		         interventions  one where the vulnerable population is isolated and the
              other is the complete lockdown. 
              Therefore, we can investigate how different types of physical distancing
              and their effectiveness affects disease spread."),
		         br(),
            p("Tab 2 (Total Infections), one can view the total number of individuals who have been."),
		         br(),
		         p("The model illustrates the role of incubation time, which provides the disease a sort of inertia."),
		         br(),
		         p("On Tab 3 (Hospitalisations), we show the number of hospitalised individuals - based on our model
analysis. There is a gap between when a person gets infected and when she/he is brought to a hospital;
therefore, we see that this number reacts slowly to a lockdown. This eect could be seen in reality as
well."),
		         h3("Can I use these results of your model for policy in specific states of India?"),
		         p("No. This model is for illustrative purposes only, and should not be used for quantitative forecasting."),
		         h3("Can I play with the code for this model?"),
		         p("Yes, please. The code for this model is available at 
		         our GitHub repository (https://github.com/tee-lab/seir-covid-india) and is available via CC-BY-SA-4.0 licence (i.e, use, modify but credit our work and share using the same license.  
		         If you have any questions regarding the modelling approaches used, please send a well-drafted email to pranavm@iisc.ac.in"),
		         h3("Why don't you show actually observed data anywhere?"),
		         p("Empirical data is subject to the dynamics of societal awareness of the disease, testing rates, stochasticity,
and other numbers. Data from tests conducted do not represent all cases of infections, due to 
somewhat restrictive testing protocol followed by the ICMR. In contrast, our model depicts all cases -- and hence will be much larger
shown by real data. Therefore, showing data together with model results in our plots could have been misleading. 
Moveover, there are many online sources where real data
are beautifully shown. We therefore refrain from showing this data."),
		         ),
		tabPanel("About", 
		         h2("Codes"),
		         p("All codes are available at our",
		           a("GitHub repository", href="https://github.com/tee-lab/seir-covid-india"), 
		           "and is available via", a("CC-BY-SA-4.0 licence.",href="https://creativecommons.org/licenses/by-sa/4.0/" )),
		         br(),
		         h2("The authors"),
		         p("This app was created by",
		         a("Pranav Minasandra", href="https://pminasandra.weebly.com/"),
		           "a BS-MS student from the Centre for Ecological Sciences (CES)
		            Indian Institute of Science (IISc), Bengaluru, India, working with", 
		            a("Vishwesha Guttal (CES, IISc)", href="https://teelab.iisc.ac.in"),
		         "and with critical inputs from Prateek Sharma (Department of Physics, IISc) 
		         and Shivakumar Jolad (FLAME University, Pune)"),
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
		           but in particular",
		           a("Alison Hill's model", href="https://alhill.shinyapps.io/COVID19seir/"),
		           "and", a("Joshua Weitz's model.", href="https://github.com/jsweitz/covid-19-ga-summer-2020"))
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

kar_ori <- as.Date('150320', format='%d%m%y')
mah_ori <- as.Date('150320', format='%d%m%y')
in_ori <- as.Date('150320', format='%d%m%y')
tn_ori <- as.Date('150320', format='%d%m%y')
kl_ori <- tn_ori
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
       xlab="", ylab="Hospitalised Individuals", 
       main="Hospitalisations",
       xaxt="n",
       lwd=3.0)
  if(input$quar=="lockdown"){
    abline(v=lockdown_begin_date)
    abline(v=lockdown_begin_date + input$lockdowndur)
  }
  if(input$quar != "none"){
    t <- seq(0, as.double(input$dur - ori), by=1)
    recompute_init_vals()
    init_vals[2:4] <- init_vals[2:4]*e_data[1]
    out2 <- data.frame(ode(y=init_vals, func=none_der_statevector, times=t, parms=params))
    points(out2$I_h_young + out2$I_h_old ~ dates, type='l', col='blue', lty='dotted',lwd=3)
  }
  axis(1, dates, format(dates, '%d/%m'), cex.axis=0.6)
  legend("topleft",
    col="blue",
         lty=c("dotted","solid"),
         legend=c("No intervention", "With intervention"),
         lwd=3.0)
})

output$plot_cis <- renderPlot({
  out <<- d()
  dates <- ori + out$time
  plot(out$CumInf_S ~ dates, 
       col="red", type='l', 
       xlab="", ylab="Population size", 
       main=" Cumulative Symptomatic Individuals",
       xaxt="n",
       lwd=3.0)
  axis(1, dates, format(dates, '%d/%m'), cex.axis=0.6)
  legend("topleft",
        col="red",
         lty=c( "dotted","solid"),
         legend=c("No intervention","With intervention"),
         lwd=3.0)
  if(input$quar=="lockdown"){
    abline(v=lockdown_begin_date)
    abline(v=lockdown_begin_date + input$lockdowndur)
  }
  if(input$quar != "none"){
    t <- seq(0, as.double(input$dur - ori), by=1)
    recompute_init_vals()
    init_vals[2:4] <- init_vals[2:4]*e_data[1]
    out2 <- data.frame(ode(y=init_vals, func=none_der_statevector, times=t, parms=params))
    points(out2$CumInf_S ~ dates, type='l', col='red', lty='dotted',lwd=3)
  }
})

# output$captionr0 <- renderText({
#   if(input$quar=="none"){
#     paste("R_0 =", round(P_a*BASAL_TRANSMISSION_RATE*Tau_a + P_s*BASAL_TRANSMISSION_RATE*2*Tau_s,3))
#   }
#   else if(input$quar=="symquar"){
#     paste("R_0 =", round(P_a*BASAL_TRANSMISSION_RATE*(Tau_a*(PROP_YOUNG*(input$lockdowneff) + 1 - input$lockdowneff)) + P_s*BASAL_TRANSMISSION_RATE*2*Tau_s*(PROP_YOUNG*(input$lockdowneff) + 1 - input$lockdowneff),3))
#   }
#   else if(input$quar=="lockdown"){
#     paste("R_0 =", (1-input$lockdowneff)*round(P_a*BASAL_TRANSMISSION_RATE*Tau_a + P_s*BASAL_TRANSMISSION_RATE*2*Tau_s,3))
#   }
# })

}

shinyApp(ui, server)
