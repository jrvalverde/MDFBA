#########################################################################################
#	NOTE: This will be the script to make a function to prepare te data before 
#	running the for loop of the function
#
#
#
#
#
data_preparation <- function(model,
			substrateRxns,
                        exclUptakeRxns,
                        verboseMode = 2, 
			biomassRxn = "",
                        nutrientChanges = NULL,
                        ...){

    ##--------------------------------------------------------------------------##
    # check prerequisites 
    if (!is(model, "modelorg")) {
      stop("adaptiveDFBA needs an object of class modelorg!")
    }

    ##--------------------------------------------------------------------------##
    #
    if (biomassRxn != '') {
	# find biomass reaction index
	# if the reaction is not found, we'll have an error
	biomassIdx <- which(react_id(model) %in% biomassRxn)
	if (length(biomassIdx) == 0) {
            # the reaction was not found
	    cat(warn, 'The specified biomass reaction', biomassRxn, '\n')
	    cat(warn, 'has not been found in the model: Ignoring it.\n')
	    biomassRxn <- ''
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
        biomassIdx <- grep('biomass', react_id(model), ignore.case=TRUE)[1]
        # if found use it
        if (length(biomassIdx) != 0) {
            biomassRxn <- react_id(model)[biomassIdx]
            cat(warn, "Setting Biomass reaction to the first reaction containing\n")
            cat('	BIOMASS in its name:', biomassRxn, '(', biomassIdx, ')\n')
        }    
        # if not, revert to objective function
	else {
	    cat(warn, "Setting Biomass reaction to default:\n")
            cat(warn, "	FIRST objective function!!!\n")
	    biomassIdx <- which(sybil::obj_coef(model) != 0)[1]
	    biomassRxn <- sybil::react_id(model)[biomassIdx][1]
	}
        cat(warn, 'IF THIS IS WRONG, USE PARAMETER "biomassRxn" TO SET IT\n\n')
    }

    cat(info, 'Biomass reaction:', biomassRxn, "(",biomassIdx, ')\n\n')

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
    exclUptakeRxnsInd <- is.element(react_id(model) ,exclUptakeRxns);
    # update exchange reaction index
    excReactInd <- excReactInd & !exclUptakeRxnsInd;
    # get reaction names
    excRxnNames <- react_id(model)[excReactInd];

    ##--------------------------------------------------------------------------##
    # Find substrate reactions
    substrateRxnsInd=(react_id(model) %in% substrateRxns)
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
        print(react_id(model)[missingSub]);
        print('Invalid substrate uptake reaction!');
    }

    ##--------------------------------------------------------------------------##
    # Initialize concentrations
    #
    # we'll use one concentration fir each possible reaction
    # and we'll start from a concentration of zero by default
    concentrations=rep(0,length(react_id(model)))
    #
    # assign each concentration its reaction name. This is useless in this
    # function, but will come handy when passing concentrations to auxiliary
    # functions at each step!
    names(concentrations) <- react_id(model)
    #
    # XXX j XXX
    # THERE WAS AN INCONSISTENCY HERE: WHEN SUBSTRATE RXNS ARE NOT IN THE SAME
    # ORDER AS THEY ARE IN THE MODEL, THIS ASSIGNMENT WILL ASSIGN THE
    # CONCENTRATIONS ERRONEOUSLY SINCE THE ORDER IN initConcentrations IS 
    # NOT SORTED BUT THE INDEX OF SUBSTRATE REACTIONS IS!!!
    #concentrations[substrateRxnsInd] = initConcentrations;
    #   using this instead is correct:
    # for each substrateRxn provided
