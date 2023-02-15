################################################
# Function: dynamicFBA
#
# Performs a dynamic flux balance analysis
# 
# The function dynamicFBA() is inspired by the function
# dynamicFBA() contained in the COBRA Toolbox.
# The algorithm is the same.

dynamicFBA.I <- function(models,
                         substrateRxns,
                         initConcentrations,
                         initBiomass,
                         timeStep,
                         nSteps,
                         exclUptakeRxns, 
                         fld = FALSE,
                         verboseMode = 3, ...) {
#PARAMETERS:
#===========
# model                 Sybil model structure (class modelorg)
# substrateRxns         List of exchange reaction names for substrates
#                       initially in the media that may change (e.g. not
#                       h2o or co2)
# initConcentrations    Initial concentrations of substrates (in the same
#                       structure as substrateRxns)(must be a list including the biomass of 
#                       each organism))
# initBiomass           Initial biomass (must be non zero and a list including the biomass of 
#                       each organism)
# timeStep              Time step size
# nSteps                Maximum number of time steps
# fld                   indicates if all fluxes at all steps will be returned.
# retOptSol             indicates if optsol calls will be returned or simple list
#
#OPTIONAL PARAMETERS
#===================
# exclUptakeRxns        List of uptake reactions whose substrate
#                       concentrations do not change (Default =
#                       {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'})
# 
#RETURN VALUES:
#=============
# concentrationMatrix   Matrix of extracellular metabolite concentrations
# excRxnNames           Names of exchange reactions for the EC metabolites
# timeVec               Vector of time points
# biomassVec            Vector of biomass values
# all_fluxes            Matrix containing the fluxes of all reactions at different steps

#
# If no initial concentration is given for a substrate that has an open
# uptake in the model (i.e. model.lb < 0) the concentration is assumed to
# be high enough to not be limiting. If the uptake rate for a nutrient is
# calculated to exceed the maximum uptake rate for that nutrient specified
# in the model and the max uptake rate specified is > 0, the maximum uptake 
# rate specified in the model is used instead of the calculated uptake
# rate.



# aids for debugging/verbose messages
    myName <- '\nIRENE 01'
     
    
    info <- paste(myName, 'INFO: ')
    warn <- paste(myName, 'WARNING: ')
    err  <- paste(myName, 'ERROR: ')

##--------------------------------------------------------------------------##
 # check prerequisites 
    
    if (is(models, "modelorg")){
        models <- list(models);
        for (i in models) {
                cat('It is a single model with name ', i@mod_name, '\n')
                
        }
        
    } else if (is.list(models)) {
        for (i in 1:length(models)) {
                if(is(models[[i]], "modelorg")) { 
                        cat('The model number ', i, ' is  ', models[[i]]@mod_name, '\n')
                } else if (!is(models[[i]], "modelorg")){
                        stop("needs an object of class modelorg! The model", models[[i]]@mod_name, " is not") 
                }
        }
    }
##--------------------------------------------------------------------------##

initDynWorld <- function(models, substrateRxns,initConcentrations,initBiomass, 
                        exclUptakeRxns, verboseMode = 3,...) {
        
     
     # Uptake reactions whose substrate concentrations do not change
     if (missing(exclUptakeRxns)){
         exclUptakeRxns = c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)');
         if (verboseMode > 2){
            print('Default extra cellular uptake reactions will be used: ')
            print(exclUptakeRxns);
         }
     }
     
     
     excRxnNamesByModel=list()
     excReactIndByModel=list()
     substrateRxnsInd=list()
     #excReactInd= vector(mode="character")
     
	
     this_species=0
     
     for (m in models) {
          this_species <- this_species +1
          # Find exchange reactions for this model
          
          excReact = findExchReact(models[[this_species]]);
          excReactInd = (react_id(m) %in% react_id(excReact));#excReact$exchange
          
          #represent extra cellular reaction with boolean vector.
          exclUptakeRxnsInd = is.element(react_id(m) ,exclUptakeRxns);
          #Exclude reactions with concentrations that will not be changed 
          excReactInd = excReactInd & !exclUptakeRxnsInd;   #excInd & ~ismember(model.rxns,exclUptakeRxns);
          excReactIndByModel[[this_species]] <- excReactInd;
          
          #get reaction names
          excRxnNames = react_id(m)[excReactInd];                #excRxnNames = model.rxns(excInd)
          excRxnNamesByModel[[this_species]] <- excRxnNames;
          
          #excluded reactions that are not in the list of reactions initially provided
         
          substrateRxnsInd[[this_species]] <- (react_id(m) %in% substrateRxns);
          # Figure out if substrate reactions are correct: all substrate reactions should be exchange reactions.
          missingSub = substrateRxnsInd[[this_species]] & !excReactInd;
          if (sum(missingSub)!=0){
              print(sum(missingSub));
              print(react_id(m)[missingSub]);
              print('Invalid substrate uptake reaction!');
          }
          
     }
     
          ##    ***********************************************************     ##
          # Initialize concentrations
          #substrateMatchInd = intersect(excRxnNames,substrateRxns);
          #IRENE: he cambiado lo de la longitud de la lista para poder
          #sacarlo del bucle de los modelos y que no se genere una cada vez
          #PROVISIONAL PERO NO VA A FUNCIONAR PARA MAS ESPECIES
          #
          concentrations=rep(0,length(substrateRxnsInd))        #concentrations=rep(0,length(react_id(model)))
          
          #table(excRxnNames
          ##vector(length=length(excRxnNames),mode="numeric");
          #concentrations[1:length(concentrations)]=0;
          concentrations[substrateRxnsInd[[this_species]]] = initConcentrations;
	 
     biomass = list();
     uptakeBound= list();
     originalBound=list();
     aboveOriginal=list();
     #reinitialize the counter
     this_species=0
     for (m in models){
          this_species= this_species +1
          ##IRENE_PROVISIONAL, CUANDO ESTÉ TODO BIEN HECHO DESDE LOS ARGUMENTOS
          ##INICIALES, EL initBiomass pasará a ser una lista y esto cambiará a
          ##initBiomass[[this_species]]
          
          biomass[[this_species]]= initBiomass;
          
          # Deal with reactions for which there are no initial concentrations
          originalBound[[this_species]] = -lowbnd(m);# take all to be able to directly update
          originalBound[[this_species]][is.na(originalBound[[this_species]])] <- 0
          noInitConcentration = (concentrations==0)&(lowbnd(m)<0)#(concentrations == 0 & originalBound > 0);
          concentrations[noInitConcentration] = 1000;

          # Initialize bounds
          #IRENE: ¿aquí va a funcionar bien?¿va a generar una matriz para cada concentración
          #y para cada microorganismo?
          uptakeBound[[this_species]]=  concentrations/(biomass[[this_species]]*timeStep);
	  uptakeBound[[this_species]][is.na(uptakeBound[[this_species]])] <-0;
          # Make sure bounds are not higher than what are specified in the model
          aboveOriginal[[this_species]] = (uptakeBound[[this_species]] > originalBound[[this_species]]) & (originalBound[[this_species]] > 0);
          aboveOriginal[[this_species]][is.na(aboveOriginal[[this_species]])] <-0;
          
          uptakeBound[[this_species]][aboveOriginal[[this_species]]] = originalBound[[this_species]][aboveOriginal[[this_species]]];
          lowbnd(m)[excReactIndByModel[[this_species]]]  = -uptakeBound[[this_species]][excReactIndByModel[[this_species]]];  
          biomassVec=list()
     	  biomassVec[[this_species]] = biomass[[this_species]]
     }
     
     #IRENE: esto habrá que cambiarlo al meter más especies pq sino cogerá
     #el de la última
     concentrationMatrix = matrix(concentrations[excReactInd]); 
      
     timeVec = 0;
     
     world = list("biomassVec"=biomassVec, 
               "biomass"=biomass, 
               "excReactIndByModel"=excReactIndByModel, 
               "concentrations"=concentrations,
               "concentrationMatrix"=concentrationMatrix,
               "uptakeBound"=uptakeBound,
               "aboveOriginal"=aboveOriginal,
               "excRxnNames"=excRxnNames,
               "timeVec"=timeVec);
     
     return(world);
    
}

world <- initDynWorld(models=models, 
	     substrateRxns = substrateRxns,
             initConcentrations = initConcentrations,
             initBiomass = initBiomass, 
             exclUptakeRxns = exclUptakeRxns, 
             verboseMode = 3,...)
             
             

##------------------------------------- Prepare Problem object --------------------##
# get OptObj instances
#lpmod <- prepProbObj(model,
#                         nCols      = react_num(model),
#                         nRows      = met_num(model),
#             #            alg        = "FBA",
#                         solver     = solver,
#                         method     = method,
#                         lpdir      = lpdir
#                         #solverParm = solverParm
#             )

##-----------------------------------------------------------------------------##
updateModels <- function(models, 
                         exclUptakeRxns, 
                         world,
                         timeStep, 
                         nSteps, 
                         fld = FALSE,
                         verboseMode = 2,
                         ...) {
                    
                  
             
 
         if (verboseMode > 2) print('Step number    Biomass\n');
         # Inititialize progress bar ...');
         #if (verboseMode == 2)  progr <- .progressBar();  
         
         biomassVec = world$biomassVec
         biomass = world$biomass 
         excReactIndByModel = world$excReactIndByModel
         concentrations = world$concentrations
         concentrationMatrix = world$concentrationMatrix
         uptakeBound = world$uptakeBound
         aboveOriginal = world$aboveOriginal
         excRxnNames = world$excRxnNames
         timeVec=world$timeVec
   
##################################
	 all_fluxes=list();
         all_stat=list();
         all_lpmod=list();
         
         #reinitialize the counter
         this_species=0;
        
         for (m in models) {     
             this_species= this_species +1;
	     
             #if (verboseMode == 2)  progr <- .progressBar(stepNo, nSteps, progr);
             # FBA can only be computed for one model
             # We have several models
             # So, we need to compute FBA for each model
             #Run FBA
             lpmod <- sybil::sysBiolAlg(m, algorithm = "fba", ...)
             all_lpmod[[this_species]]= lpmod;
         }
         
             #We need to have all the solutions together for later steps
             all_sol=list();
             all_mu=list();
         for (stepNo in 1:nSteps){

             
             #reinitialize the counter in each loop
             this_species = 0
             
             for (m in models) {
             
                 this_species = this_species + 1
                 sol = sybil::optimizeProb(all_lpmod[[this_species]])
                 all_sol[[this_species]]= sol;
                 summary(sol);
                 mu =  all_sol[[this_species]]$obj;  ##objvalue sol.f
                 all_mu[[this_species]]= mu;
                 
                 ############## JR ###############
                 if ( length(checkSolStat(all_sol[[this_species]]$stat,solver(problem(all_lpmod[[this_species]]))))!=0 ){## checkSolStat
                    cat('No feasible solution - nutrients exhausted', '\n');
                    cat('sol stat is',all_sol[[this_species]]$stat, '\n')
                    break;
                 }
                 ############## JR ###############
                 
                 all_stat[[this_species]] = c(all_stat[[this_species]],all_sol[[this_species]]$stat);
                
                 uptakeFlux = all_sol[[this_species]]$fluxes[excReactIndByModel[[this_species]]];
                 
                 # Now we can calculate how much biomass we have for this species               
                 biomass[[this_species]] = biomass[[this_species]]*exp(all_mu[[this_species]]*timeStep);
                 
                 # Now, we want to remember the biomass at each step
                 biomassVec[[this_species]] = c(biomassVec[[this_species]],biomass[[this_species]]);
         	 
                 
                 if(fld){
                         if (stepNo == 1) {
                                 all_fluxes = all_sol[[this_species]]$fluxes;
                         }else{
                                 all_fluxes = cbind(all_fluxes,all_sol[[this_species]]$fluxes);
                         }
                 }


                 # and now we can update the total concentrations by removing the
                 # amount consumed by this species
                 
                
                # Update concentrations
                concentrations[excReactIndByModel[[this_species]]]= concentrations[excReactIndByModel[[this_species]]] - uptakeFlux/all_mu[[this_species]]*biomass[[this_species]]*(1-exp(all_mu[[this_species]]*timeStep));
                #concentrations = concentrations + uptakeFlux*biomass*timeStep;
                concentrations[concentrations <= 0] = 0;
                concentrationMatrix = cbind(concentrationMatrix,concentrations[excReactIndByModel[[this_species]]]);

                # Update bounds for uptake reactions
                uptakeBound[excReactIndByModel[[this_species]]] =  concentrations[excReactIndByModel[[this_species]]]/(biomass[[this_species]]*timeStep);
                # This is to avoid any numerical issues
                uptakeBound[uptakeBound > 1000] = 1000;
                # Figure out if the computed bounds were above the original bounds
                aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
                # Revert to original bounds if the rate was too high
                uptakeBound[aboveOriginal] = originalBound[aboveOriginal];# uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
                uptakeBound=ifelse(abs(uptakeBound) < 1e-9,0,uptakeBound);
              ## Change lower bounds according to the result of last step
                #lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];  
                uppb_tmp <- getColsUppBnds(problem(all_lpmod[[this_species]]), which(excReactIndByModel[[this_species]]));
                changeColsBnds(problem(all_lpmod[[this_species]]),which(excReactIndByModel[[this_species]]),lb=-uptakeBound[excReactIndByModel[[i]]],ub=uppb_tmp);

                if (verboseMode > 2) print(paste(stepNo,sep="    ",biomass));
                #waitbar(stepNo/nSteps,h);
                timeVec = c(timeVec,stepNo*timeStep);
            }
         }# end loop
         
         #row.names(concentrationMatrix) = excReactIndByModel        #row.names(concentrationMatrix)=react_id(m)[excReactInd];
         
         updateWorld = list("biomassVec"=biomassVec, 
               "excReactIndByModel"=excReactIndByModel, 
               "concentrationMatrix"=concentrationMatrix,
               "all_fluxes"=all_fluxes,
               "all_stat"=all_stat,
               "excRxnNames"=excRxnNames,
               "all_lpmod"= all_lpmod,
               "stepNo"=stepNo,
               "timeVec"= timeVec);
    }

updateWorld <- updateModels(models=models,
                        exclUptakeRxns = exclUptakeRxns, 
                        world = world,
                        timeStep = timeStep,
                        nSteps = nSteps,
                        fld = FALSE,
                        verboseMode = 2, 
                        ...)
             
                    
# browser();


makeOutput <- function(models,
                       updateWorld, 
                       retOptSol = TRUE, 
                       fld = FALSE, 
                       verboseMode = 2, 
                       ...) {
        ## Preparing OUTPUT
        #concentrationMatrix,excRxnNames,timeVec,biomassVec
        
        
        biomassVec=updateWorld$biomassVec;
        excReactIndByModel= updateWorld$excReactIndByModel;
        concentrationMatrix=updateWorld$concentrationMatrix;
        all_fluxes= updateWorld$all_fluxes;
        all_stat=updateWorld$all_stat;
        all_lpmod=updateWorld$all_lpmod;
        excRxnNames=updateWorld$excRxnNames;
        stepNo= updateWorld$stepNo;
        timeVec=updateWorld$timeVec;
        
        total_react = NULL
        total_met = NULL
        
        for (m in models) {
        
                total_react = c(total_react, (react_id(m)[!(react_id(m) %in% total_react)]))
                total_met = c(total_met, (met_id(m)[!(met_id(m) %in% total_met)]))
                
        }
        
        
        
        if (isTRUE(retOptSol)) {
                
            for (m in 1:length(all_lpmod)) {
            
                 if(is.null(all_fluxes[[m]])) all_fluxes=as.matrix(NA);
                     return (optsol_dynamicFBA.I(solver = solver(problem(all_lpmod[[m]])),
                                               method = method(problem(all_lpmod[[m]])),
                                               nprob  = stepNo,
                                               ncols  = length(total_react),
                                               nrows  = length(total_met),
                                               fld    = fld,
                                               all_fluxes = all_fluxes,
                                               concmat=concentrationMatrix,
                                               exRxn=excRxnNames,
                                               tmVec=timeVec,  
                                               bmVec=biomassVec)
                       )
            }
        }else{
                        
             return(optsol <- list(nprob  = stepNo,
                                  ncols  = length(total_react),
                                  nrows  = length(total_met),
                                  all_fluxes = all_fluxes,
                                  all_stat=all_stat,
                                  concentrationMatrix=concentrationMatrix,
                                  excRxnNames=excRxnNames,
                                  timeVec=timeVec,  
                                  biomassVec=biomassVec
                                          ))
        
        }
}

 
 makeOutput(models= models,
            updateWorld = updateWorld, 
            retOptSol = TRUE, 
            fld = FALSE, 
            verboseMode = 2, ...)
 }
