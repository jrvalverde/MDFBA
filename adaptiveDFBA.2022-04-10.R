################################################
#' adaptiveDFBA
#'
#' \code{adaptiveDFBA} performs an adaptive dynamic flux balance analysis
#'
#'         Calculate concentrations of metabolites of exchange reactions at defined
#'         time points  given the initial concentrations and optional dynamically
#'         changing  constraints and conditions.
#' 
#' To accomplish its task this function calls \code{\link{optimizeProb}} function to 
#'         get the fluxes then updates the concentrations and the reaction
#'         boundaries ..etc. as in Varma et al. (1994). It also applies a number of
#'         corrections that increase reliability, resilience and control over the
#'         simulation, as well as the possiility to produce additional informative
#'         output.
#' 
#' It is used like sybilDynFBA with the possibility of adding various extra optional
#'         parameters: one to specify the actual biomass reaction, one to control
#'         reaction  limits during the simulation, one to control concentration of
#'         metabolites and one to control output of additional information.
#' 
#' This function does away with a number of significant restrictions present in the
#'         original implementation, the OpenCOBRA implementation in MatLab, and the
#'         sybilDynFBA implementation in R (which is itself a direct translation of
#'         the OpenCOBRA reference implementation).
#' 
#' The code has been rearranged to avoid duplications and simplify the algorithm. FBA
#'         calculation is now done using a model that has been updated at each step
#'         to meet optional adaptive constraints. The model is reverted back to the
#'         original at the end of the calculation.
#' 
#' You no longer need to provide substrates and initial concentrations in the same
#' 	order as they are listed in the model.
#' 
#' Biomass is now computed using the Biomass reaction that you specify, not the
#' 	objective function (as it was before). This allows you to specify an
#'         objective function other than growth (e.g. secretion) and still 
#'         calculate biomass, growth and concentrations properly.
#' 
#' Input concentrations that have been explicitly set to zero are now respected and
#'         not set to infinity, plus it is ensured that only exchange reactions are
#'         affected.
#' 
#' Uptake rate actualization at each step now ensures that invalid values are
#'         always avoided
#' 
#' You can modify any reaction limit at any time, thus making it possible to 
#'         steer the  simulation at your will: you may change one or both limits on
#'         one or more reactions. This is useful to control systems where
#'         constraints need to be coerced (such as introducing a lag-phase at the
#'         start of growth) or where multiple goals must evolve independently (e.g.
#'         secretion is active at specific growth phases only), or when you have
#'         multiple objectives but want to evolve them independently of each other
#'         (e.g. you want to maximize growth and secretion, but the cell does not
#'         secrete always at the same rate with respect to growth). 
#' 
#' There is provision for modifying concentrations as well by either
#'         supplying a table with incremental changes to be applied at each time
#'         step, or a function that will get the current concentrations and return
#'         the new ones. This is useful to model situations where some nutrients
#'         are added at specific times (using the table) or feed-back fed-batch
#'         models where a sensor detects nutrient changes to trigger a modification
#'         of nutrient concentrations (using the sensor-feeder function).
#' 
#' If no initial concentration is given for a substrate that has an open
#'     uptake in the model (i.e. model.lb. < 0) the concentration is assumed to be
#'     high enough to not be limiting. If the uptake rate for a nutrient is
#'     calculated to exceed the maximum uptake rate for that nutrient specified in
#'     the model and the max \emph{uptake} rate specified is > 0, the maximum
#'     uptake  rate specified in the model is used instead of the calculated uptake
#'     rate.
#' 
#' \emph{dynamicConstraints} is a table providing the values of the reaction
#'     rate for each dynamically adapated reaction at each of the nSteps time
#'     points. During the simulation, the reaction rate will be coerced at each
#'     time point to the value provided in this table. Thus, in order to use this
#'     table, one must provide a value for each time point. This will typically not
#'     be available for all but a few experientally measured time points, hence, 
#'     values at intermediate time points must be somehow interpolated prior to
#'     calling this function using existing \code{R} facilities. Users are advised
#'     to carry out experiments at key points, visualize the results, and decide on
#'     an interpolation method that is BIOLOGICALLY MEANINGFUL to compute the table
#'     values.
#' 
#' \emph{BE CAREFUL WITH THE MODEL YOU CHOOSE FOR GENERATING YOUR VALUES!}
#'     
#' As an example, you could interpolate intermediate growth points using
#'     a sigmoid, but model other reactions (such as secretion) using a method
#'     that makes less assumptions (such as linear interpolation) and is
#'     scientifically safer.
#' 
#' You may modify only one rate or both. The table must consist of columns
#'     with one cell per time point. Each column pertains to one reaction. If the
#'     column is named after a reaction followed by \code{[.low.]} or \code{[.upp.]}
#'     (e.g. "\code{Ex_gly(e)[low]}" or "\code{Ex_gly(e).upp.}") then the values in
#'     that column will be used to modify only the corresponding (lower or upper)
#'     bound of that reaction. Note that you may have both, one column ending in
#'     .low. and another in .upp. for the same reaction to modify both limits. As a
#'     short-hand (but be careful and ensure it makes sense), you can set both
#'     limits to the same value by using a column named after the reaction without
#'     extension (e.g. a column called "\code{Ex_gly(e)}" would set both bounds of
#'     reaction  "\code{EX_gly(e)}" to the same value as specified in that column.
#'     If a column or bound appears more than once, the last instance
#'     (right-most) will be applied last and prevail.
#' 
#' There is also an (as yet untested) possibility of having a finer control:
#'     you can provide a function instead of a table. If you do this function must
#'     take as arguments the current model, concentrations, fluxes and step number
#'     and return a new model to substitute the original one. Inside the function
#'     you will have all the information needed to take decisions, such as, e.g. if
#'     we detect that time < 10h, induce a lag phase of growth; if we detect that
#'     some metabolite is low, repress coupled uptake reactions or activate
#'     alternative uptake mechanisms... We also get the last fluxes, so that inside
#'     the function we may also detect (if desired) activation or repression of a
#'     key reaction/route and trigger the corresponding changes. This is 
#'     intentionally left undocumented until we have mode experience in using
#'     it.
#' 
#' \emph{nutrientChanges}: We consider two possibilities for modifying concentrations
#'     (experimentwise): one is that you add nutrients at predefined time points,
#'     which is implemented as a table, and the other is that you have a sensor to
#'     detect when nutrients are below a given value and then add the nutrients,
#'     which is implemented as a function:
#' 
#' nutrientChanges may be a table providing the delta changes in nutrient
#'     concentration at each time point. The changes at every time point should be
#'     provided sorted in time order. Only final time points may be missing, not
#'     intermediate time points. Deltas may positive (add nutrient), negative
#'     (remove nutrient) or zero (leave unchanged).
#' 
#' nutrientChanges may be a function, designed to address the second case by
#'     simulating a fed-batch environent:  this function must take as only input
#'     argument \emph{all} the current  concentrations (i.e. the function will act
#'     as a "sensor" that is able to read concentrations at each time step) and
#'     must return \emph{all} the new concentrations (i.e. it also acts as a
#'     "feeder" able to modify concentrations at each time step; and note that it
#'     may "suck" nutrients if you want, much as you can with the table of delta
#'     changes at fixed times). This is left undocumented intentionally until we
#'     have more expertise using it.
#' 
#' @title adaptive dynamic flux balance analysis
#
#PARAMETERS:
#===========
#
#' @param model                 An object of class \code{\link{modelorg}} containing the
#'      			metabolic model. It should have been read from a TSV or XML
#'				file using sybil tools, and modified appropriately to define
#'				the simulation constrains.
#' @param substrateRxns         List of exchange reaction names for substrates
#'                      	initially in the media that may change (e.g. not
#'                      	h2o or co2)x
#' @param initConcentrations    Initial concentrations of substrates (in the same
#'                      	structure as substrateRxns)
#' @param initBiomass           Initial biomass (must be non zero)
#' @param timeStep              Time step size. Define the time intervals to evaluate the 
#'				problem at.
#' @param nSteps                Maximum number of time steps (time intervals) to evaluate, 
#'		      		the procedure may stop before completing this number when 
#'      			the FBA fails to find a viable solution (indicating that the
#'      			available substrates are not enough to support a solution / viability).
#'
#
#OPTIONAL PARAMETERS
#===================
#
#' @param retOptSol		Boolean. Indicates if optsol class will be returned or a simple list.
#'			Default: \code{FALSE}
#
#' @param fld                   Boolean. Indicates if all fluxes at all steps will be returned.
#'			Default: \code{FALSE}
#
#' @param exclUptakeRxns        List of uptake reactions whose substrate
#'                      concentrations do not change 
#'			Default: \code{c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)')}
# 
#' @param biomassRxn		Name of the biomass reaction so we can use it
#'			to compute growth, concentrations, etc...
#
#' @param dynamicConstraints    Values of dynamically changing reaction rates
#'			Default: \code{NULL}
#
#' @param nutrientChanges	A facility to simulate changing concentrations
#'			of nutrients: it may be either an array of delta values
#'			to add at each time point, or a function that will
#'			take the current concentrations and return the new ones (the old concentrations
#'      		will be substituted by the new ones, so ensure you return all).
#'			Default: \code{NULL}
#
#' @param method	Which method should be used to solve the system at each
#'			step. Currently it can be 'FBA', 'lpFBA', 'directFBA',
#'                      'lpMTF', 'directMTF' or 'MTF'. Default:
#'			\code{"FBA"}
#
#' @param verboseMode 	An integer value indicating the amount of output to stdout:\cr
#'				0 (quiet)\cr
#'				1 (input and status messages, plus a progress idicator)\cr
#'				2 (like 1 plus additional information and a grahical progress indicator)\cr
#'				3 (information on the computations being made plus a
#'      			  graphical progress bar.)\cr
#'				4 (additional information on the computations being made plus a
#'      			  graphical progress bar.)\cr
#'				  Default: \code{0}
#
#' @param \dots		Further arguments passed to \code{\linkS4class{sysBiolAlg}}.  Argument
#'      		\code{solverParm} is a good candidate.
#
#RETURN VALUES:
#=============
#' @return	returns \code{\link{optsol_dynamicFBA}} if retOptSol was set to \code{TRUE}
#'    or a list with the results containing\cr
#'    
#'    nprob = the number of steps actually computed\cr
#'    
#'    ncols = the number of columns (reactions) in the problem solved\cr
#'    
#'    nrows = the number of rows (metabolites) in the problem solved\cr
#'    
#'    all_fluxes = all the fluxes computed at each step\cr
#'    
#'    concentrationMatrix = the concentrations computed at each step\cr
#'    
#'    excRxnNames = the names of excluded reactions\cr
#'    
#'    timeVec = a vector with the times at which the problem was evaluated\cr
#'    
#'    biomassVec = a vector with the biomass computed at each time\cr
#
#% Add one or more standard keywords, see file 'KEYWORDS' in the
#% R documentation directory.
#' @keywords optimize 
# keywords Adaptive DFBA, DFBA, FBA, MTF, Metabolism, Metabolic analysis
#
#' @usage   adaptiveDFBA(model, substrateRxns, initConcentrations, initBiomass, 
#'          timeStep, nSteps, 
#'          exclUptakeRxns,
#'          retOptSol = TRUE,
#'          fld = FALSE, verboseMode = 2, 
#'          verboseMode = 2,
#'          biomassRxn = "",
#'          dynamicConstraints,
#'          nutrientChanges, ...)
#
#' @references     Valverde, J. R. et al. (2018) Manuscript in preparation.
#'
#'    Varma, A. and Palsson, B.O. 1994. Stoichiometric flux balance models
#'    quantitatively predict growth and metabolic by-product secretion in wild-type
#'    Escherichia coli W3110.  \emph{Appl Environ Microbiol} 60: 3724-3731.\cr
#'
#'    Valverde, Jose R., Sonia Gullon, and Rafael P. Mellado.  (2018) Modelling the
#'    metabolism of protein secretion through the Tat route in Streptomyces lividans.
#'    \emph{BMC microbiology} 18.1: 59.\cr
#'
#'    Quantitative prediction of cellular metabolism with constraint-based models: the
#'    COBRA Toolbox. \emph{Nat Protoc} \bold{2}, 727--738.\cr 
#'
#' @seealso 	\code{\link{modelorg}}, \code{\link{optsol_dynamicFBA}}, 
#'		\code{\link{optimizeProb}},  \code{\link{sysBiolAlg}}, 
#'		\code{\link{SYBIL_SETTINGS}}
#' 
#' @author	(C) Jose R Valverde, EMBnet/CNB, CSIC, 2018
#'
#	Released under EU-GPL
#
#' @examples 
#'	## The examples here require and optimization package, such as glpkAPI,
#'	## installed. If that package is not available, you have to set
#'	## the argument 'solver' (the default is: solver = "glpk") to the
#'	## linear programming solver you have installed (e.g. "clpAPI",
#'	## "lpSolveAPI", "cplexAPI").
#'
#'	## load the example data set
#'	data(Ec_core)
#'	  lowbnd(Ec_core)[react_id(Ec_core)=='EX_glc(e)']=-10;
#'	  lowbnd(Ec_core)[react_id(Ec_core)=='EX_o2(e)']=-18;
#'	## run adaptiveDFBA(), Ec_df will be an object of class \code{\link{optsol_dynamicFBA}}
#'	Ec_df <- adaptiveDFBA(Ec_core,substrateRxns={'EX_glc(e)'},initConcentrations=10,
#'		initBiomass=.035,timeStep=.25,nSteps=17,verboseMode=3)
#'
#'	## plot biomass and reactions
#'	plot(Ec_df,plotRxns=c('EX_glc(e)','EX_ac(e)'));
#'
#'      # Run \code{adaptiveDFBA} using the specified dynamic constraints 
#'	data(SlividansWT)
#'	verbose <- 2
#'      af_sol <- adaptiveDFBA(model, substrateRxns=substrateRxns, 
#'   		    initConcentrations=initConcentrations,
#'   		    initBiomass = initBiomass,
#'   		    timeStep = timeStep,
#'   		    nSteps = nSteps,
#'   		    exclUptakeRxns=exclUptakeRxns,
#'   		    retOptSol=TRUE,
#'   		    fld=TRUE,
#'   		    biomassRxn=biomassRxn,
#'   		    dynamicConstraints=dynamicConstraints,
#'   		    verboseMode=verbose);
#'       plot(af_sol, plotRxns=plotRxns);
#'
#'	\dontrun{
#'	## open S.lividans
#'	data(Slividans)
#'	# the distributed model has all exchange reactions without any limits
#'	model <- Slividans;
#'	biomass <- 'Biomass_SLI'	
#'	nSteps=50
#'	tStep=1
#'
#'	## load a tab-separated file with a list of dynamic constraints
#'	## e.g. given 'exp.tab'
#'	## Time	Biomass_SLI	SEC_AML
#'	## 0	0.01	0
#'	## 10	0.011	0
#'	## 16	0.14	1.63e-6
#'	## 24	0.135	1.63e-5
#'	## 48   0.001	2.23e-7
#'	dat <- read.table('exp.tab', sep='\t', header=TRUE)
#'	reactions <- names(dat)
#'	reactions <- reactions[ -which(reactions %in% 'Time') ]
#'
#'	## create dynamicConstraints data frame with target times
#'	dc = data.frame( Time = seq(1, nSteps * tStep, by=tStep))
#'	## 'Time' will be ignored, but helps readability and checking
#'	
#'	x <- dat['Time'][,1]	## Prepare a vector with time values
#'	for (r in reactions) {
#'	    y <- dat[r][,1]		## Prepare a vector with reaction values
#'	    ## make a linear interpolation
#'	    interdata <- data.frame(
#'		    approx(x, y, xout=seq(1, nSteps * tStep, by=tStep))
#'	            )
#'	    # add the corresponding column to the output data frame
#'	    dc[r] <- interdata['y']
#'	}
#'	  
#'	## Perform the simulation
#'	#   Note that we provide an insufficient medium for growth
#'	adfba <- adaptiveDFBA(model, 
#'	          substrateRxns=c('EX_mnl(e)'),	# we should better use NMMP
#'	          initConcentrations=54,
#'	          initBiomass=.01,
#'	          timeStep=tStep,
#'                nSteps=nSteps,
#'	          verboseMode=5,
#'	          biomassRxn=biomass,
#'	          dynamicConstraints=dc,
#'	          retOptSol=TRUE, fld=TRUE)
#'	
#'	print(plot(adfba, plotRxns=c('EX_mnl(e)')))
#'	}
#'
#
#' @export
#
# Function adaptiveDFBA() implements additional logic to provide
# controls of a dynamic FBA simulation. It is based on function
# dynamicFBA() from the sybilDynFBA package.
#
# The function dynamicFBA() is inspired by the function
# dynamicFBA() contained in the COBRA Toolbox.
# The algorithm is the same.
#
#	" Oui, je vois l'heure; il est l'Éternité! "
#			L'Horloge. Charles Baudelaire.
#
adaptiveDFBA <- function (models,
			substrateRxns,
                        initConcentrations,
                        initBiomass,
                        timeStep,
                        nSteps,
                        exclUptakeRxns,
			retOptSol = TRUE,
                  	fld = FALSE,
                        verboseMode = 2, 
			biomassRxn = "",
                        dynamicConstraints = NULL,
                        nutrientChanges = NULL,
                        method = "FBA",
                        ...){
                        
    

    # aids for debugging/verbose messages
    myName <- 'adaptiveDFBA IRENE 01'
    cat(myName)
    info <- paste(myName, 'INFO: ')
    warn <- paste(myName, 'WARNING: ')
    err  <- paste(myName, 'ERROR: ')

    # MAGIC MAGIC MAGIC MAGIC
    # we use here a 'magic' number, it should be checked periodically
    #	to verify it remains the same.
    # OK return status of optsol
    optsol.OK <- 5
    # and here we set the default method when only "FBA" is specified
    defaultFBA <- 'lpFBA'

    ##no tengo muy claro que esto esté muy bien hecho pero lo que he intentado es: si en 'models' se ha
    ##añadido un solo modelo, compruebo que es un objeto que pertenece a la clase modelorg???????,
    ## genero una lista de un solo elemento para que no haya que duplicar el resto del script y
    ## me devuelve el nombre del modelo.
    ##Si se han metido varios modelos en una lista, se comprueba que cada uno pertenezca a la clase 
    ##modelorg y que nos devuelva el nombre del modelo
    ##Si no es ninguna de las dos cosas, es que alguno no es un modelo y rompe el bucle con un
    ##mensaje de alerta
    
    if (class(models)=="modelorg"){
    	models <- list(models)
        return(models)
    }else if (is.list(models)){
    	for (model in models) {
        	if(class(model)=="modelorg"){ 
        		return(model)
                }else(
                	cat(model(name), "\s is not an object of the class 'modelorg'")
        }
    }else {
    	cat(models(name), "\s is not an object of the class 'modelorg'")
    	break
    }
    
    
    
    
    
    
    for(m in models) {
        ##--------------------------------------------------------------------------##
        # check prerequisites 
        if (!is(m, "modelorg")) {
          stop("adaptiveDFBA needs an object of class modelorg!")
        }

        ##--------------------------------------------------------------------------##
        #
        if (biomassRxn != '') {
	    # find biomass reaction index
	    # if the reaction is not found, we'll have an error
	    biomassIdx(m) <- which(react_id(m) %in% biomassRxn(m))
	    if (length(biomassIdx(m)) == 0) {
                # the reaction was not found
	        cat(warn, 'The specified biomass reaction', biomassRxn, '\n')
	        cat(warn, 'has not been found in the model: Ignoring it.\n')
	        biomassRxn(m) <- ''
	    }
        }

        # if the biomass reaction has not been specified, we'll first look for
        # a reaction called 'biomass' and, if none is found, we'll use the
        # first objective function.
        #
        # NOTE that there may be more than one objective function: we'll
        # settle for the first, but this may not be the correct one!
        # e.g. if obj <- c("SECRETION", "BIOMASS"), then growth will be
        # tied to secretion, not to biomass!!!
        if (biomassRxn == '' || missing(biomassRxn)) {
            # Try to find the real Biomass function before accepting to use only
            # the objective function (which might even be multiple objectives)
            #
            # THIS IS VERY, VERY RISKY
            #	It might be that the biomass reaction is called something else
            #	and the one named *biomass* is not the real Biomass reaction
            #	or that there is more than one reaction called '*biomass*'
            #   and the proper one is not the first. Hopefully we will get it
            #	right
            #
            biomassIdx(m) <- grep('biomass(m)', react_id(m), ignore.case=TRUE)[1]
            # if found use it
            if (length(biomassIdx(m) != 0) {
                biomassRxn(m) <- react_id(m)[biomassIdx]
                cat(warn, "Setting Biomass reaction to the first reaction containing\n")
                cat('	BIOMASS in its name:', biomassRxn(m), '(', biomassIdx(m), ')\n')
            }    
            # if not, revert to objective function
	    else {
	        cat(warn, "Setting Biomass reaction to default:\n")
                cat(warn, "	FIRST objective function!!!\n")
	        biomassIdx(m) <- which(sybil::obj_coef(m) != 0)[1]
	        biomassRxn(m) <- sybil::react_id(m)[biomassIdx][1]
	    }
            cat(warn, 'IF THIS IS WRONG, USE PARAMETER "biomassRxn" TO SET IT\n\n')
        }

        cat(info, 'Biomass reaction:', biomassRxn(m), "(",biomassIdx(m) ')\n\n')

        # If no initial concentration is given for a substrate that has an open
        # uptake in the model (i.e. model.lb < 0) the concentration is assumed to
        # be high enough to not be limiting. If the possible uptake rate calculated
        # for a nutrient exceeds the maximum uptake rate specified in the model 
        # for that nutrient and the max uptake rate specified is > 0, the maximum 
        # uptake limit specified in the model is used instead of the calculated 
        # maximum possible rate.

        if (verboseMode > 0) {
            cat(info, 'Concentrations provided\n')
            #print(substrateRxns)
            #print(initConcentrations)
            for (i in 1:length(substrateRxns))
                cat(paste(i, substrateRxns[i], initConcentrations[i], '\n', sep="	"))
	    cat('\n')
        }

        ##--------------------------------------------------------------------------##

        # Uptake reactions whose substrate concentrations do not change
        if (missing(exclUptakeRxns)){
            exclUptakeRxns = c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)');
        }
        if (verboseMode > 0){
            cat(info, 
                'Excluded uptake reactions\n',
       	        '	(substrate concentration does not change):\n')
            for (i in exclUptakeRxns) {
                cat(i, '	', sep='')
            }
            cat('\n\n')
        }

        ##--------------------------------------------------------------------------##
        # Find exchange reactions (all)

        excReact <- findExchReact(model);
        # represent exchange reactions as a boolean (index) vector
        excReactInd <- (react_id(model) %in% react_id(excReact));
        # Exclude reactions with concentrations that will not be changed 
        exclUptakeRxnsInd <- is.element(react_id(m) ,exclUptakeRxns);
        # update exchange reaction index
        excReactInd <- excReactInd & !exclUptakeRxnsInd;
        # get reaction names
        excRxnNames <- react_id(m)[excReactInd];

        ##--------------------------------------------------------------------------##
        # Find substrate reactions
        substrateRxnsInd=(react_id(m) %in% substrateRxns)
        if (verboseMode > 1) {
            cat(info, "Substrate reactions provided\n")
            for (i in substrateRxns) {
                cat(i, '	', sep='')
            }
            cat('\n\n')
        }
        # Figure out if substrate reactions are correct:
        #	all substrate reactions should be exchange reactions.
        #	A substrate that is not in an exchange reaction cannot be used.
        missingSub = substrateRxnsInd & !excReactInd;
        if (sum(missingSub) != 0){
	    cat(warn, 'All substrate reactions should be exchange reactions!\n')
            print(sum(missingSub));
            print(react_id(m)[missingSub]);
            print('Invalid substrate uptake reaction!');
        }

        ##--------------------------------------------------------------------------##
        # Initialize concentrations
        #
        # we'll use one concentration fir each possible reaction
        # and we'll start from a concentration of zero by default
        concentrations=rep(0,length(react_id(m)))
        #
        # assign each concentration its reaction name. This is useless in this
        # function, but will come handy when passing concentrations to auxiliary
        # functions at each step!
        names(concentrations) <- react_id(m)
        #
        # XXX j XXX
        # THERE WAS AN INCONSISTENCY HERE: WHEN SUBSTRATE RXNS ARE NOT IN THE SAME
        # ORDER AS THEY ARE IN THE MODEL, THIS ASSIGNMENT WILL ASSIGN THE
        # CONCENTRATIONS ERRONEOUSLY SINCE THE ORDER IN initConcentrations IS 
        # NOT SORTED BUT THE INDEX OF SUBSTRATE REACTIONS IS!!!
        #concentrations[substrateRxnsInd] = initConcentrations;
        #   using this instead is correct:
        # for each substrateRxn provided
    }

    #
    # XXX j XXX
    # THERE WAS AN INCONSISTENCY HERE: WHEN SUBSTRATE RXNS ARE NOT IN THE SAME
    # ORDER AS THEY ARE IN THE MODEL, THIS ASSIGNMENT WILL ASSIGN THE
    # CONCENTRATIONS ERRONEOUSLY SINCE THE ORDER IN initConcentrations IS 
    # NOT SORTED BUT THE INDEX OF SUBSTRATE REACTIONS IS!!!
    #concentrations[substrateRxnsInd] = initConcentrations;
    #   using this instead is correct:
    # for each substrateRxn provided
    ##irene: esto es para las concentraciones iniciales nada mas
    #cojo el primer elemento de la lista pero supongo que habra que añadir algún tipo de control para
    ## controlar que las concentraciones iniciales sean las mismas en todos los modelos
    for (i in 1:length(substrateRxns)) {
        # note that there MUST be a concentration for each substrate
        # and that concentrations and substrates MUST have been provided 
        # in the same order!!!
        concentrations[react_id(models[1]) == substrateRxns[i]] <- initConcentrations[i]
    }

    # Deal with uptake (lowbnd < 0) reactions for which there are no initial 
    # concentrations.
    # If the initial concentration for an uptake is zero, then set it to 1000
    #	This implies we cannot set an initial concentration to zero!
    #   Corollaries: for every uptake reaction there MUST be an initial
    #   concentration or, if not, 1000 will be used instead,
    #	and that an initial concentration can be provided to a secretion 
    #	reaction (and will not be touched here)
    #(concentrations == 0 & originalUptake > 0);
    #
    #	a second problem is that bidirectional reactions will be given
    # 	a concentration as well.
    #noInitConcentration <- (concentrations==0)&(lowbnd(model)<0)

    # it is most likely better to use the following
    #	concentrations is defined for all reactions in the model
    #	select those that
    #		are zero
    #		had a lower bound < 0
    #		are exchange reactions
    #		and were NOT explicitly provided
    ##irene: aquí lo mismo que arriba, cojo los del primer modelo porque deberían ser iguales
    noInitConcentration <- (concentrations==0) &
    			   (lowbnd(models[1])<0) &
                           (excReactInd) &
                           (!(substrateRxnsInd))
                          
    # set these to a concentration of 1000
    concentrations[noInitConcentration] <- 1000;

    if (verboseMode > 2) {
	cat(info, '\n\nConcentrations of substrates\n')
	#for (i in substrateRxns) {
        for (i in 1:length(concentrations)) {
          if (concentrations[i] != 0) {
            cat(react_id(models[1])[i], concentrations[i], '\n', sep="	")
          }
        }
        cat('\n')
    }
    if (verboseMode > 3) {
	cat(info, 'Concentrations of all exchange reactions\n')
	#print(substrateRxns)
	#print(initConcentrations)
	for (i in 1:length(excReactInd)) {
            if (excReactInd[i]) 
                cat(paste(i, react_id(models[1])[i], concentrations[i], sep="	" ), '\n')
        }
	cat('\n')
    }

    ##--------------------------------------------------------------------------##
    # set biomass to initial biomass
    for (m in models){
    	biomass(m) = initBiomass(m);

    	if (verboseMode > 0) {
        	cat(info, 'Initial Biomass', biomass(m), '(', biomassRxn(m), ')\n')
        }
    

        ##--------------------------------------------------------------------------##
        # Initialize bounds
        #
        # save original lower bounds for all reactions (MAKING uptake positive
        # and secretion negative !!!)
        originalUptake(m) <- -lowbnd(m)

        # these are to restore the model to its original bounds before returning
        originalLow(m) <- lowbnd(m)
        originalUpp(m) <- uppbnd(m)
    
        #	compute maximal uptake bound as the concentration available per
        #	biomass unit per unit of time lapse
        #	I.e.: this is the maximum amount of each metabolite available,
        #	hence cells cannot consume more than this even if the maximum
        #	uptake limit in the model permits higher potential consumption
        #
        #	May it be that this could cause a problem if the computed uptake
        #	where smaller than the upperbound?
        #	e.g. if EX_X is between (-)1.5 and (-)0.5 and the computed availability
        #   is < 0.5, then the new bounds would be (-)<0.5 and (-0.5) with
        #	lowbnd > upbnd
        #
        #	As long as the uptake limits are of the form <0 .. >=0 this shouldn't
        #   be a problem. If it were, the simulation should fail and then we
        #	can fix this below.
        # we can take at most biomass*timeStep per hour
        uptakeBound(m) =  concentrations(m)/(biomass(m)*timeStep);

        # Make sure bounds are NOT higher than those specified in the model
        #
        # identify uptake reactions whose potential UPTAKE bound is larger than
        # what the model allows
        aboveOriginal(m) = (uptakeBound(m) > originalUptake(m)) & (originalUptake(m) > 0);
        # then, for these, we will use the smaller original bound in the model
        # instead of the maximum permitted by the medium concentrations
        uptakeBound(m)[aboveOriginal] = originalUptake(m)[aboveOriginal];
        # uptake us defined by a negative sign
        lowbnd(m)[excReactInd]  = -uptakeBound(m)[excReactInd];    

        ##--------------------------------------------------------------------------##
        # Initialize simulation
        concentrationMatrix <- cbind(concentrations(m)[excReactInd]);
        biomassVec(m) <- biomass(m);
        timeVec <- 0;

        # start with fluxes at time zero = 0 and status = OK
        ##todo esto debería estar relacionado pq al final es flujo tambien entre especies
        all_fluxes <- c(rep(0, length(react_id(m))));	# start with flux=0 at t=0
        all_stat <- optsol.OK;		# start with an OK status

        #all_fluxes=NULL;
        #all_stat=NULL;


        ##------------------------------------- Prepare Problem object --------------------##
        # get OptObj instance for LP problem solving (FBA)
        #	This should not be needed, we'll do it before each FBA calculation 
        #	below so we always use the latest modified version of the model
        #   but since we may need to return an lpmod data object at the end...
        lpmod <- sybil::sysBiolAlg(m, algorithm = "fba", ...)

        ##------------------------------------- Prepare Progress bar ----------------------##
        # Inititialize progress bar
        if (verboseMode < 2) { pb <- txtProgressBar(0, nSteps) }
        if (verboseMode > 1) {
            # requires package tcltk
            library(tcltk)
            pb <- tkProgressBar(title = "DFBA progress", min = 0,
                                max = nSteps, width = 300)
        }
        
        old.model <- m
       
     }

    ##-----------------------------------------------------------------------------##
    # Do the simulation
    for (stepNo in 1:nSteps){

        if (verboseMode < 2) { setTxtProgressBar(pb, stepNo) }
	# for verboseMode == 2 we will output data during calculations
        if (verboseMode > 1) { 
            setTkProgressBar(pb, stepNo, 
            	label=paste( round(stepNo/nSteps*100, 0),
                "% done"))
	}
        if (verboseMode > 1) cat('\n\n', info, ' Step ', stepNo, '\n', sep='')

	# before doing the calculation we will apply any needed modifications
        # to the model/problem/environment to solve
        
        # First, modify the models to reflect time-specific constraints
        for (m in models) {
                      m <- dynamicConstraints(m, concentrations, sol$fluxes, stepNo

             if ((! missing(dynamicConstraints)) && 
                 (! is.null(dynamicConstraints)) ) {
	         if (verboseMode > 3) {
                     cat(info, 'Applying dynamic constraints\n')
                 }
                 if (is.function(dynamicConstraints)) {
                     # This should allow modifying the model depending on current
                     # conditions, e.g. if we detect that time < 10h induce a lag
                     # phase of growth, if we detect that some metabolite is low,
                     # repress coupled uptake reactions or activate alternative
                     # uptake mechanisms...
                     # We also pass the last fluxes, so that inside the function
                     # we may also detect (if desired) activation or repression of
                     # a key reaction/route and trigger the corresponding changes.
	             if (verboseMode > 3) {
                         cat(info, 'Calling dynamicConstraints()\n')
                     }

                 } else {
	             ## Coerce exchange rates for this step to the values provided
	             #

                     # Remove the time column (may be in any position)
                     # We effectively ignore the time column!
                     #	No wonder, as there should not be any "Time" reaction
                     #	It is only an aid for the caller
                     # In fact the following line should be unneeded since we
                     # test for valid reaction names below, and indeed, perhaps
                     # we should remove it in case there is ever a "Time" reaction.
	             #	Oh, and this could be out of the loop!
                     dynamicallyConstrainedRxns <- names(dynamicConstraints)[ -which(names(dynamicConstraints) %in% "Time")]

	             if (verboseMode > 3) {
                         cat('\nDynamical constraints are\n')
	                 print(dynamicallyConstrainedRxns)
	                 print(dynamicConstraints)
                         cat('\n')
	             }

                     for (dcr in dynamicallyConstrainedRxns) {
                         # Currently we demand that all time-step values are provided.
                         # but we will check anyway, just in case
                         #     we should likely use instead
                         # if (dynamicConstraints[dynamicConstraints["Time",1]==stepNo, dcr]) {
                         #     but we can't be sure there is a column named "Time"
                         if (stepNo <= dim(dynamicConstraints[dcr])[1]) {
                             #get reaction name and value
                             rxnName <- sub("\\[(low|upp)\\]", "", dcr)	# remove [low] or [upp]
	                     rxnValue <- dynamicConstraints[stepNo, dcr]
                             if (is.na(rxnValue)) next ;

                             # check that the reaction is a valid reaction
                             if (! (rxnName %in% react_id(m))) {
                    	         cat('\n', warn, 'ignoring unknown reaction ', rxnName, '\n')
			         next;
                             }
		             if (verboseMode > 3) {
	                         cat(info, 'Setting ', dcr,
                                 ' (', rxnName, ') ', 
                                 'at step ', stepNo, 
                                 ' to ', dynamicConstraints[dcr][1,stepNo],
                                 ' = ', rxnValue,
                                 '\n', sep='')
		             }
		             # get old limits for verboseMode reporting below
                             oldlow(m) <- lowbnd(m)[react_id(m) == rxnName]
		             oldupp(m) <- uppbnd(m)[react_id(m) == rxnName]

		             # update corresponding limits according to original reaction name
                             if (grepl("\\[low\\]", dcr)) {
                                 # reaction.low.
                                 lowbnd(m)[react_id(m) == rxnName] <- rxnValue
                                 # a negative lower bound implies uptake, so we change the sign
			         # note: this will be done afterwards all at once
                                 #originalUptake[react_id(model) == rxnName] <- -rxnValue
		             } else if (grepl("\\[upp\\]", dcr)) {
                                 # reaction.upp.
                                 uppbnd(m)[react_id(m) == rxnName] <- rxnValue
                             } else {
                                 # reaction
                                 lowbnd(m)[react_id(m) == rxnName] <- rxnValue
                                 # a negative lower bound implies uptake, so we change the sign
                                 # note: this will be done afterwards all at once
                                 #originalUptake[react_id(model) == rxnName] <- -rxnValue
                    	         uppbnd(m)[react_id(m) == rxnName] <- rxnValue
		             }
	                 }
                         else {
                             #if no new values are available, we retain the old value
                             cat(warn, 'no constraints for ', dcr, ' at ', stepNo, '\n')
                             cat(warn, 'leaving constraints untouched\n')
		         }

                         if (verboseMode > 3) {
                             cat(info,
                                 'Step ', stepNo, 
                                 ' Rxn ',rxnName, 
                                 ' olo ', oldlow,
                                 ' oup ', oldupp,
                    	         ' exp ', rxnValue,
                                 ' nlo ', lowbnd(model)[react_id(m) == rxnName],
                                 ' nup ', uppbnd(model)[react_id(m) == rxnName],
                                 ' dim ', dim(dynamicConstraints[dcr])[1],
                                 '\n'
                                 )
                         }
                         # finally, update problem as well (using the model as master)
        	         # -- not needed because we remake lpmod before the FBA 
                         # calculation below
                         #changeColsBnds(problem(lpmod),which(excReactInd),
                         #	lb=lowbnd(model)[excReactInd],
                         #   ub=uppbnd(model)[excReactInd]);
                         # and if we do not use an LPMOD object, it's unneeded
	             } # for (dcr in dynamicallyConstrainedRxns)

                 } # if (is.function(dynamicConstraints))

             } # if there are dynamicConstraints
             
        }

        # the model limits have changed
        #	so we need to upfate originalUptake to reflect the new limits
##        originalUptake[excReactInd] <- -lowbnd(model)[excReactInd]
        # this is where we do the update referred to earlier in the loop
        modified <- lowbnd(m) != lowbnd(old.model)
        originalUptake[modified & excReactInd] = -lowbnd(m)[modified & excReactInd]
        

	# Second, modify concentrations with the incremental
        # changes of media applied at this step, using either
        # an static table of deltas, or a sensor-feeder function
        if ( ! missing(nutrientChanges) && ! is.null(nutrientChanges) )  {
	    if (verboseMode > 0) {
                cat(info, 'Dynamically modifying concentrations\n')
            }
            if (is.function(nutrientChanges)) {
                # the change may be a function of current concentration or
                # of time
                concentrations <- nutrientChanges(concentrations, tStep, stepNo)
            } else {
                # nutrientChanges MUST be a table of concentration deltas at
                # each and every time point. At least we will use it
                # as such. If it is not, either we will get an error or
                # nothing will be modified (unless terribly unlucky).
        	updNutrNames <- names(nutrientChanges)[ -which(names(nutrientChanges) %in% "Time")]
                for (nutrient in updNutrNames) {
		    # modify its concentration as requested
                    if (stepNo <= dim(nutrientChanges[nutrient])[1]) {
                        # NOTE that we do not check if the value is positive
                        # or not, i.e. we may both add and remove nutrients
                        # at any step (depending on the sign of the delta)
                        for (m in models) {
		            concentrations[react_id(m) == nutrient] <- 
                                concentrations[react_id(m) == nutrient] +
                                nutrientChanges[stepNo, nutrient]
                        }
                                
                    }
                }
	    }
        }


        # Update bounds for uptake reactions
        #uptakeBound[excReactInd] =  concentrations[excReactInd]/(biomass*timeStep);
	# or alternatively for all (we'll only use exchange reactions anyway)
        #
        # NOTE: consider applying a MM formula upon request
        for (m in models) {
            uptakeBound(m) =  concentrations/(biomass(m)*timeStep);
            cat('The uptake bound for the model ', m, 'in the time step = ', timeStep, ' is = ', uptakeBound(m))

            # Figure out the newly computed bounds were above the original bounds
            ### j ### in the original formula, a computed uptake was considered 
            ### above the original only if the original uptake was > 0
            ###
            ### That is, for secretion reactions, it was ignored, which means that
            ### if the concentration is > 0, the potential uptake is > 0 and, since
            ### in a secretion reaction the original uptake is < 0, it is not 
            ### corrected here and it will not be considered above the original, 
            ### and thus will not be removed from the aboveOriginal vector,
            ### which means that it will not be reset to the model value,
            ### but left as an uptake and included in the updated lower limits,
            ### resulting in secretion being converted to updatke.
            ###
            ### In other words, whenever an uptake is possible, it will be set
            ### to uptake, even if we desire to allow secretion, not
            ### permitting secretion as long as there is any metabolite 
            ### available outside.
            ###
            ### That may make initially sense from a maximal instantaneous efficiency
            ### point of view: use all resources you can to your inmmediate
            ### advantage.
            ###
            ### What we see, instead, at least in S. lividans, is that cells have a
            ### genetically programmed built-in "inefficiency" to excrete metabolites
            ### to the medium whenever they are in a surplus of nutrients, producing
            ### a "reserve" of nutrients that are not used but set aside.
            ###
            ### This works to the advantage of cells as it results in a "reserve"
            ### of nutrients that are "saved" for later, when the main nutrient is
            ### exhausted. 
            ###
            ### Then, when the main nutrients are exhausted, or close to, the cells
            ### can start using up this "reserve" of previously produced, less-desirable
            ### metabolites/excreted nutrients, thus extending their life and increasing
            ### the chances of survival until a new nutrient source is potentially
            ### available.
            ###
            ### I.e. producing this reserve is suboptimal from an instantaneous point
            ### of viewm as it requires an energy consumption, but this energy is in 
            ### excess during growth in excess nutrients, so despite being an 
            ### "unnecessary" expense, it can be done without otherwise compromising 
            ### growth or consuming nutrients faster (which would lead to faster 
            ### exhaustion and extinction). However, on the long term, the strategy
            ### is better as it allows creation of a reserve or "savings account" to
            ### extend survival chances.
            ###
            ### I suspect there is an economics lesson hidden here regarding resource
            ### overuse and starvation vs. sustainable economy. I also suspect this 
            ### is not specific of S. lividans alone and may hid a lesson related to 
            ### obesity, fruits, seed, etc... , but we better leave this aside here.
            ###
            ### In any case, here is a good reason to coerce some metabolic rates at 
            ### specific time points.

            # aboveOriginal = (uptakeBound > originalUptake) & (originalUptake > 0);
            ### This fixes it
            aboveOriginal = (uptakeBound > originalUptake);

            # Revert to original bounds if the rate was too high
            uptakeBound(m)[aboveOriginal] = originalUptake(m)[aboveOriginal];

            # This induces numerical instabilities in some cases (when the original
            # concentration was set to 1000 and the metabolite is excreted).
            ###uptakeBound[uptakeBound > 1000] = 1000;
            # This should not be needed
            ###uptakeBound=ifelse(abs(uptakeBound) < 1e-9,0,uptakeBound);

            ## Change lower bounds according to the result of last step
	    #     a) in model
            lowbnd(m)[excReactInd]  = -uptakeBound[m][excReactInd];  
            #     b) in LP problem
            #     this is no longer used since we update the whole lpmod
            #     from the complete model (to reflect all changes made to the
            #     model), and we have just corrected uptakes on the model. This
            #     way we update all the lpmod at once before actual calculation.
            #lowb_tmp <- getColsLowBnds(problem(lpmod), which(excReactInd));
            #uppb_tmp <- getColsUppBnds(problem(lpmod), which(excReactInd));
            #changeColsBnds(problem(lpmod),which(excReactInd),
            #               lb=-uptakeBound[excReactInd],ub=uppb_tmp);
            #
	    #     update LP-problem with new model for (LP-based) calculations;
            #     without this it wouldn't honor the new Biomass limits nor any
            #     other internal changes made to the model
	    lpmod <- sybil::sysBiolAlg(m, algorithm = "fba", ...)
               
            if (verboseMode > 3) {
                lplow <- getColsLowBnds(problem(lpmod),
                                        which(c(rep(TRUE, length(react_id(m))))))
                lpupp <- getColsUppBnds(problem(lpmod), 
                                        which(c(rep(TRUE, length(react_id(m))))))
                cat('\n\n', info, 'step', stepNo, 'uptake rate corrections:\n');
                cat('\n    Name	PreviousBound	UptakeBound	LowBnd	LpModLo	UppBnd	LpModUp	Above\n')
                for (i in 1:length(excReactInd)) {
                    if (excReactInd[i]) {
                        cat(paste(react_id(m)[i], 
                            originalUptake[i], 
                            -uptakeBound[i], 
                            lowbnd(m)[i],
                            lplow[i],
                            uppbnd(m)[i],
                            lpupp[i],
                            aboveOriginal[i],
                            '\n', sep="	"))
                        if (lowbnd(m)[i] != lplow[i]) cat('	BAD LOW	')
                        if (uppbnd(m)[i] != lpupp[i]) cat('	BAD UPP	')
                        cat('\n')
                    }
                }
            }
            
        
        
            # we are ready for the calculations now:
            # set default for "FBA"
            if (method == "FBA") { method <- defaultFBA }

	    if (method == "directFBA") {
                ## use direct FBA on the model instead of an lp opt object
                ### j ###
                #
	        # optimizeProb for class sysBiolAlg uses a different logic
                # from optimizeProb for class modelorg and both give different
                # flux distributions at some time points, hence concentrations 
                # and system may evolve (slightly) differently.
                #
                # Oddly, the two give coincident flux distributions in most
                # steps; even after diverging for a few steps, they revert to
                # identical fluxes (implying that the concentration changes tend to
                # compensate and the system oscillates around (or tends to) some
                # equilibrium). So the overall curve doesn't change appreciably,
                # but I will not feel comfortable until I have thoroughly
                # traced the internal logics and can see the source of these
                # differences. Be careful when choosing either method.
                #
                # For now, after many simulations, the overall curve always 
                # looks the same even if some individual calculation is slightly
                # different, so it "seems" harmless.
                # 
                # Once it has been investigated, we will run speed tests and
                # choose the best/fastest one as the default.
                #
	        fba.sol <- optimizeProb(m, algorithm="fba", retOptSol=FALSE)
                if (verboseMode > 3) {
                    # Check against lpmod-based FBA calculation.
                    # In effect this means we are doing TWO FBA calculations
                    # and checking one against the other.
            	    #lpmod <- sybil::sysBiolAlg(model, algorithm = "fba", ...) # already done above
	            sol <- sybil::optimizeProb(lpmod);
                    n = 0
		    for (i in 1:length(sol$fluxes)) {
    		        if (sol$fluxes[i] != fba.sol$fluxes[i]){
                            if (verboseMode > 3) 
        		        cat(err, 'DIFF', i, react_id(m)[i], sol$fluxes[i], fba.sol$fluxes[i], '\n', sep='	')
                            n = n+1
     		        }
	            }
                    if (n > 0)
                        cat(warn, 'DIFF', n, 'differences found between MODEL- and LPMOD-based FBA fluxes\n')
                }
	        sol <- fba.sol
                if (verboseMode > 2) {
                    cat('\n')
                    cat(info, 'FBA step =', stepNo, 'ok =', sol$ok, 'stat =', sol$stat, '\n')
	        }
	        ## check solution status (did we find a solution?)
                if ( (sol$ok != 0) || (sol$stat != optsol.OK) ){   ## checkSolStat
                    cat('\n')
                    cat(warn, 'FBA: No feasible solution - nutrients exhausted\n');
                    break;
	        }
            }
	    else if (method == "lpFBA") {
                sol = sybil::optimizeProb(lpmod);
	        if (verboseMode > 2) {
                    cat('\n')
                    cat(info, 'FBA step =', stepNo, 'ok =', sol$ok, 'stat =', sol$stat, '\n')
	        }
                ## checkSolStat
                if ( length(checkSolStat(sol$stat,solver(problem(lpmod))))!=0 ){
                    cat(warn, 'FBA: No feasible solution - nutrients exhausted\n');
                    if (verboseMode > 3) {
            	        # get indexes of problems for which there was no solution
                        notSolved <- checkSolStat(sol$stat,solver(problem(lpmod)))
                        for (i in notSolved) print(i)
                        cat('\n')
                    }
                    break;
                }
            }         
            else if ((method == "MTF") ||
                     (method == "directMTF") ||
                     (method == "lpMTF")) {
                # use MTF instead of FBA
                # choose one of these three
                if (method == "lpMTF") {
                    # this is the two-step lpmod-based approach:
                    # first solve FBA to get a value for the objective function(s)
                    sol = sybil::optimizeProb(lpmod);
                    if ( length(checkSolStat(sol$stat,solver(problem(lpmod))))!=0 ){
                        cat(warn, 'lpFBA: No feasible solution - nutrients exhausted\n');
                        break
                    }
		    # then solve MTF with the computed objective
	            mtf.sol <- optimizeProb(m, algorithm="mtf",  mtfobj=sol$obj)

	        } else if (method == "directMTF") {
                    # first solve FBA to get a value for the objective function(s)
                    sol <- optimizeProb(m, algorithm="fba", retOptSol=FALSE)
                    if ( (sol$ok != 0) || (sol$stat != optsol.OK) ){   ## checkSolStat
                        cat('\n')
                        cat(warn, 'directFBA: No feasible solution - nutrients exhausted\n');
                        break;
	            }
		    # then solve MTF with the computed objective
                    mtf.sol <- optimizeProb(m, algorithm="mtf",  mtfobj=sol$obj)

                } else { # if (method == "MTF")
                    # this one does an FBA internally and so does not need
                    # a previous FBA calculation
                    mtf.sol <- optimizeProb(m, algorithm="mtf") 
                    # should give a message saying FBA is computed
                    # but it doesn't, although it seems to work 
                    # correctly nevertheless, I need to trace the 
                    # internal logic of Sybil.
	        }
                # NOTE that we could also call MTF it with retOptSol=FALSE
                # as we did in FBA. This is only a matter of taste.

                # check and update the solution
	        mtf.ok <- exit_code(checkOptSol(mtf.sol, onlywarn=FALSE))
                #mtf.ok <- checkOptSol(mtf.sol, onlywarn=TRUE)
                mtf.stat <- status_code(checkOptSol(mtf.sol))
                mtf.obj <- mod_obj(mtf.sol)
	        if (verboseMode > 2) {
                    cat('\n')
                    cat(info, 'MTF step =', stepNo, 'ok =', mtf.ok, 'stat =', mtf.stat, '\n')
	        }
                if ( (mtf.ok != 0) || (mtf.stat != optsol.OK) ){   ## checkSolStat
                    cat('\n')
                    cat(warn, 'MTF: No feasible solution - nutrients exhausted\n');
	            break;
	        }
                mtf.fluxdist <- getFluxDist(mtf.sol)
	        mtf.fluxes <- c(rep(0.0, length(mtf.fluxdist)))
                #mtf.flux.exchRxn <- getFluxDist(mtf.sol, ex)
                for (i in 1:length(mtf.fluxdist)) mtf.fluxes <- mtf.fluxdist
                # convert to list
                sol <- list(
                    ok     = mtf.ok,
                    stat   = mtf.stat,
                    obj    = mtf.obj,
                    fluxes = mtf.fluxes
                    #, exchRxnFluxes = mtf.flux.exchRxn
                )
	        if (verboseMode > 2) {
                    cat('\n')
                    cat(info, 'MTF step =', stepNo, 'ok =', sol$ok, 'stat =', sol$stat, '\n')
	        }
            } else {
                cat(err, 'Unknown method, should be one of \n\t"FBA", "directFBA", "lpFBA", "MTF", "directMTF or "lpMTF"\n')
            }
        
	    ### j ###
            # UNDOCUMENTED FEATURE
            #	WORKS IN PROGRESS!!!	WORKS IN PROGRESS!!!
            #	Compute FVA at each time point
            #	THIS WILL SLOW DOWN CALCULATIONS HORRIBLY
            #	BUT WILL ALLOW US TO COLLECT ADDITIONAL INFORMATION FOR
            #	OUR STATISTICAL ANALYSES
            #	WE STILL NEED TO THINK OUT THE BEST WAY TO PRODUCE THE
            #	OUTPUT SO THAT IT IS AMENABLE FOR STATISTICAL ANALYSIS
            #	WITHOUT BREAKING COMPATIBILITY WITH sybil::dynamicFBA
            #
            #	Probably we should store them in a data frame with two
            #	matrices (low/high) and then save the data frame at the
            #	end of the calculation. Choice of time as rows or columns
            #	will depend on the statistics we finally use, but this
            #	still is a work in progress. Or maybe we should use only
            #	two long columns and add the time step as a factor, or
            #	use a list...
            #
            #	For now, we'll do with creating a series of output files.
            #
            if (verboseMode > 5) {
                # Do FVA and save the output
                ranges <- fluxVar(m, percentage=80, verboseMode=verbose)
                print("\n\nBeing playful, aren't you? ;-)\n\n")

	        tabFVAvalues <- paste(m@mod_name, 'ADFBA-FVA-', stepNo, '.tab', sep='')
                cat(paste('\nSaving optimized output from FVA to "', tabFVAvalues,'"\n', sep=""))
                cat(paste('\nTotal length: ', length((ranges)), 
    	                  '\nNo. of reactions: ', length(ranges)/2, '\n\n'))

                nreact=length(ranges)/2
                for (i in 1:nreact) {
                    if (verbose > 6) {
                        print(paste('FVA,', stepNo, ',', i, ',', react_id(m)[i], ',', 
    			            lp_obj(ranges)[i], ',', lp_obj(ranges)[i+nreact]))
                    }
                    out <- paste(i, '	', react_id(m)[i], '	', 
    		                lp_obj(ranges)[i], '	', lp_obj(ranges)[i+nreact])
                    cat(out, file=tabFVAvalues, sep="\n", append=TRUE)
                }
                print(paste('--- END FVA', stepNo, '---\n'))
            }

            # if the objective function is not Biomass, this is plainly wrong!
            mu_obj =  sol$obj;  ##objvalue sol.f
	    # XXX j XXX any of these two alternatives should fix the issue
	    mu_bmidx(m) <- sol$fluxes(m)[biomassIdx];
            mu_bmrxn(m) <- sol$fluxes[react_id(m) == biomassRxn]
            mu(m) <- mu_bmidx(m)
            if (verboseMode > 2) {
	        cat(info, 'SOL$obj:', sol$obj, '\n')
                cat(info, 'mu(obj)', mu_obj, 
            	    'mu(biomassIdx(m)/',biomassIdx(m),')', mu_bmidx(m), 
                    'mu(biomassRxn(m)/',biomassRxn(m),')', mu_bmrxn(m), '\n')
                cat(info, 'Using MU =', mu(m), '\n')
                cat(info, 'Max. Biomass rate', uppbnd(m)[biomassIdx], '\n')
                if (mu > (uppbnd(m)[biomassIdx]))
                    cat(err, 'MU', mu(m) '> biomass(m)[upp]', uppbnd(m)[biomassIdx], '\n')
            }
            all_stat = c(all_stat,sol$stat)
	    # compute new biomass
	    if (verboseMode > 2) {
                computed_biomass(m) = biomass(m) * exp(sol$obj * timeStep);
                cat(info, "Computed biomass =", computed_biomass(m), "\n")
	    }

            biomass(m) = biomass(m)*exp(mu(m)*timeStep);
            #biomass = biomass*(1+mu*timeStep);
            biomassVec(m) = c(biomassVec(m),biomass(m));
	    # add the fluxes to the table of flux changes over time
            if(fld){
                if (stepNo == 1) {
	                all_fluxes = sol$fluxes;
                }else{
	                all_fluxes = cbind(all_fluxes,sol$fluxes);
                }
            }

	    # get uptake fluxes after FBA/MTF
            uptakeFlux = sol$fluxes[excReactInd];
            # Update concentrations
            concentrations[excReactInd]= sum(concentrations(m)[excReactInd]) - sum(uptakeFlux/mu(m)*biomass(m)*(1-exp(mu(m)*timeStep));
            concentrations[concentrations <= 0] = 0;
            # save concentrations
            concentrationMatrix = cbind(concentrationMatrix,concentrations[excReactInd]);
            if (verboseMode > 2) {
                cat('\n\n')
                cat(info, 'Concentration of substrates after', stepNo, 'steps\n')
                for (i in substrateRxns)
                    cat(i, concentrations[react_id(m) == i], '\n', sep=" ")
                #for (i in react_id(excReact))
                #    print(concentrations[excReactInd])
                cat('\n')
            }

            timeVec = c(timeVec,stepNo*timeStep);

	    if (verboseMode > 1) cat('Objective at step', stepNo, '(', stepNo*timeStep ,'h ):', sol$obj, '\n');
            if (verboseMode > 1) cat('Biomass of', m, 'at t =', stepNo, '(', stepNo*timeStep, "h ):", biomass(m), '\n');
            
        }

    }# end loop
    # close progress bar
    close(pb)

    ##--------------------------------------------------------------------------##
    # Simulation completed... prepare RETURN output
    # browser();
    #
    # restore model to original state
    for (m in models) {
        lowbnd(m) <- originalLow
        uppbnd(m) <- originalUpp

        # add names to the concentrationMatrix
        row.names(concentrationMatrix) <- rbind(react_id(m)[excReactInd]);
        row.names(all_fluxes) <- rbind(react_id(m));
        colnames(concentrationMatrix) <- timeVec;
        colnames(all_fluxes) <- timeVec[-1];
    }

    
    ##irene: esto va despuesw de los bucles y todo pq ya es el output, habría que generar un documento con 
    #todo esto en vez de solo imprimirlo en la pantalla
    
    
    ## Preparing OUTPUT
    #concentrationMatrix,excRxnNames,timeVec,biomassVec
    if (isTRUE(retOptSol)) {
        # return an optSol object
        if(is.null(all_fluxes)) all_fluxes=as.matrix(NA);
	return (optsol_dynamicFBA(solver = solver(problem(lpmod)),
				  method = method(problem(lpmod)),
				  nprob  = stepNo,
				  ncols  = react_num,
				  nrows  = met_num,
				  fld    = fld,
				  all_fluxes = all_fluxes,
				  concmat= concentrationMatrix,
				  exRxn  = excRxnNames,
				  tmVec  = timeVec,  
				  bmVec  = biomassVec
				  )
	       )
    }else{ # return a list instead of an optSol object
        return(optsol <- list(nprob  = stepNo,
		              ncols  = react_num,
		              nrows  = met_num,
		              all_fluxes = all_fluxes,
		              all_stat= all_stat,
		              concentrationMatrix=concentrationMatrix,
		              excRxnNames=excRxnNames,
		              timeVec= timeVec,  
		              biomassVec=biomassVec
		              )
               )
    }
}
