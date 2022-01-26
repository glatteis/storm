#include "EqualParameterReducer.h"
#include <cstdint>
#include <queue>
#include <string>
#include "adapters/RationalFunctionAdapter.h"
#include "logic/UntilFormula.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/StandardRewardModel.h"
#include "storage/BitVector.h"
#include "storm-pars/transformer/SparseParametricDtmcSimplifier.h"
#include "utility/constants.h"
#include "utility/macros.h"

namespace storm {
namespace derivative {

models::sparse::Dtmc<RationalFunction> EqualParameterReducer::minimizeEqualParameters(models::sparse::Dtmc<RationalFunction> dtmc, modelchecker::CheckTask<logic::Formula, RationalNumber> const& checkTask) {
    storage::SparseMatrix<RationalFunction> transitionMatrix = dtmc.getTransitionMatrix();
    
    // Repeat this algorithm until we can't mimimize anything anymore
    bool somethingChanged = true;
    auto allParameters = storm::models::sparse::getAllParameters(dtmc);
    
    models::sparse::StateLabeling newLabels(dtmc.getStateLabeling());

    // Check the reward model - do not touch states with rewards
    boost::optional<std::vector<RationalFunction>> stateRewardVector;
    boost::optional<std::string> stateRewardName;
    if (checkTask.getFormula().isRewardOperatorFormula()) {
        if (checkTask.isRewardModelSet()) {
            dtmc.reduceToStateBasedRewards();
            stateRewardVector = dtmc.getRewardModel(checkTask.getRewardModel()).getStateRewardVector();
            stateRewardName = checkTask.getRewardModel();
        } else {
            dtmc.reduceToStateBasedRewards();
            stateRewardVector = dtmc.getRewardModel("").getStateRewardVector();
            stateRewardName = dtmc.getUniqueRewardModelName();
        }
    }
    
    while (somethingChanged) {
        somethingChanged = false;
        // Search backwards from parameters to where they join their paths
        // Tally up the probabilities
        // Multiple backwards transitions not supported yet
        //
        // Maps are from state that is the last one before the parameter transition occurs
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::vector<uint_fast64_t>>> alreadyVisitedStates;
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::vector<RationalNumber>>> talliedUpProbabilities;
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, bool>> doneSearching;
        // States visited from first with p to second
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, uint_fast64_t>> increasingSuccessors;
        // States visited from first with 1-p to second
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, uint_fast64_t>> decreasingSuccessors;
        
        // We only join states that have all the same labels.
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::set<std::string>>> joiningLabels;
        
        // As I can't create my own, save the ps from the matrix to re-use them
        std::map<storm::RationalFunctionVariable, RationalFunction> pFunctions;

        uint_fast64_t numberOfSearchingTransitions = 0;
        for (auto const& parameter : allParameters) {
            alreadyVisitedStates[parameter] = std::map<uint_fast64_t, std::vector<uint_fast64_t>>();
            talliedUpProbabilities[parameter] = std::map<uint_fast64_t, std::vector<RationalNumber>>();
            increasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
            decreasingSuccessors[parameter] = std::map<uint_fast64_t, uint_fast64_t>();
            joiningLabels[parameter] = std::map<uint_fast64_t, std::set<std::string>>();
        }
        // Search for all occurences of parameters and fill up the maps
        for (uint_fast64_t row = 0; row < transitionMatrix.getRowCount(); row++) {
            for (auto const& entry : transitionMatrix.getRow(row)) {
                if (entry.getValue().isConstant()) {
                    continue;
                }
                STORM_LOG_ASSERT(entry.getValue().gatherVariables().size() == 1, "Flip minimization only supports transitions with a single parameter.");
                auto parameter = *entry.getValue().gatherVariables().begin();
                auto cache = entry.getValue().nominatorAsPolynomial().pCache();
                RationalFunction parameterAsFunction = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(parameter), cache));
                if (entry.getValue() == parameterAsFunction) {
                    pFunctions[parameter] = entry.getValue();
                    numberOfSearchingTransitions++;
                    auto parameter = entry.getValue().nominatorAsPolynomial().getSingleVariable();
                    STORM_LOG_ASSERT(alreadyVisitedStates[parameter].count(row) == 0, "Flip minimization only supports simple pMCs.");
                    alreadyVisitedStates[parameter][row] = std::vector<uint_fast64_t>();
                    alreadyVisitedStates[parameter][row].push_back(row);
                    talliedUpProbabilities[parameter][row] = std::vector<RationalNumber>();
                    talliedUpProbabilities[parameter][row].push_back(utility::one<RationalNumber>());
                    increasingSuccessors[parameter][row] = entry.getColumn();
                    doneSearching[parameter][row] = false;
                    joiningLabels[parameter][row] = newLabels.getLabelsOfState(row);
                } else if (utility::one<RationalFunction>() - entry.getValue() == parameterAsFunction) {
                    decreasingSuccessors[parameter][row] = entry.getColumn();
                } else {
                    STORM_LOG_ASSERT(false, "Flip minimization only supports simple pMCs.");
                }
            }
        }
        // Step by step, do a backwards search
        uint_fast64_t transitionsDoneSearching = 0;
        auto backwardsTransitions = transitionMatrix.transpose(true);
        // Final results of the search: Sets of sets of states that we can transform together
        // The maps are <parameter> -> <state where the paths join> -> <set of states that leads there, with this probability>
        std::map<storm::RationalFunctionVariable, std::map<uint_fast64_t, std::map<uint_fast64_t, RationalNumber>>> joinedFirstStates;
        for (auto const& parameter : allParameters) {
            joinedFirstStates[parameter] = std::map<uint_fast64_t, std::map<uint_fast64_t, RationalNumber>>();
        }
        while (transitionsDoneSearching < numberOfSearchingTransitions) {
            for (auto const& pair : alreadyVisitedStates) {
                auto parameter = pair.first;
                auto visitedStateMap = pair.second;
                // Step one step back for parameter and firstState
                for (auto const& pair2 : visitedStateMap) {
                    auto firstState = pair2.first;
                    auto visitedStates = pair2.second;
                    if (doneSearching[parameter][firstState]) {
                        continue;
                    }
                    auto row = backwardsTransitions.getRow(visitedStates.back());
                    
                    storage::MatrixEntry<uint_fast64_t, RationalFunction> entry;
                    bool entryFound = false;
                    for (auto const& loopEntry : row) {
                        if (loopEntry.getValue().isZero()) {
                            continue;
                        }
                        // We only support one backlink per state for now, so if we already found an entry, abort
                        if (entryFound) {
                            entryFound = false;
                            break;
                        }
                        entry = loopEntry;
                        entryFound = true;
                    }
                    if (!entryFound) {
                        transitionsDoneSearching++;
                        doneSearching[parameter][firstState] = true;
                        continue;
                    }

                    // Only search until another parametric transition occurs, then we can't continue anympore
                    if (!entry.getValue().isConstant()) {
                        transitionsDoneSearching++;
                        doneSearching[parameter][firstState] = true;
                        continue;
                    }
                    // If the last label is different, we stop here (but still join if possible)
                    if (joiningLabels[parameter][firstState] != dtmc.getLabelsOfState(entry.getColumn())) {
                        transitionsDoneSearching++;
                        doneSearching[parameter][firstState] = true;
                    }
                    // If the state reward vector is not zero, we stop here (but still join if possible)
                    if (stateRewardVector && !(*stateRewardVector)[visitedStates.back()].isZero()) {
                        transitionsDoneSearching++;
                        doneSearching[parameter][firstState] = true;
                    }

                    auto constantProbability = utility::convertNumber<RationalNumber>(entry.getValue());
                    auto currentState = entry.getColumn();
                    
                    // The next probability is the previous one times the considered tranistion
                    talliedUpProbabilities[parameter][firstState].push_back(talliedUpProbabilities[parameter][firstState].back() * constantProbability);
                    alreadyVisitedStates[parameter][firstState].push_back(currentState);

                    // Check if any search of the same parameter has already gotten here
                    for (auto const& pair3 : alreadyVisitedStates[parameter]) {
                        auto otherFirstState = pair3.first;
                        // Don't join a state with itself
                        if (otherFirstState == firstState) {
                            continue;
                        }
                        // If the labels don't match, this is not a match
                        if (joiningLabels[parameter][otherFirstState] != joiningLabels[parameter][firstState]) {
                            continue;
                        }
                        auto otherVisitedStates = pair3.second;
                        // TODO maybe use a set if this is too slow
                        bool done = false;
                        for (uint_fast64_t i = 0; i < otherVisitedStates.size(); i++) {
                            auto state = otherVisitedStates[i];
                            if (state == currentState) {
                                // Wow!!! Another search of the same parameter has gotten here.
                                transitionsDoneSearching++;
                                doneSearching[parameter][firstState] = true;
                                transitionsDoneSearching++;
                                doneSearching[parameter][otherFirstState] = true;
                                
                                if (!joinedFirstStates[parameter].count(state)) {
                                    joinedFirstStates[parameter][state] = std::map<uint_fast64_t, RationalNumber>();
                                }
                                RationalNumber p1 = talliedUpProbabilities[parameter][firstState].back();
                                joinedFirstStates[parameter][state][firstState] = p1;
                                RationalNumber p2 = talliedUpProbabilities[parameter][otherFirstState][i];
                                joinedFirstStates[parameter][state][otherFirstState] = p2;

                                // We are all done. We can continue the search for the next, unrelated thing.
                                done = true;
                                break;
                            }
                        }
                        if (done) {
                            break;
                        }
                    }
                }
            }
        }
        
        // Check if we found any parameters that we can join into one transition, and if so, do it
        for (auto const& parameter : allParameters) {
            if (joinedFirstStates[parameter].size() > 0) {
                somethingChanged = true;
                break;
            }
        }
        if (!somethingChanged) {
            break;
        }

        storage::SparseMatrix<RationalFunction> wipMatrix(transitionMatrix); 
        for (auto const& parameter : allParameters) {
            // We add three new states to the DTMC for this procedure for every state where parameters join
            if (joinedFirstStates[parameter].size() == 0) {
                continue;
            }
            std::cout << "joining " << joinedFirstStates[parameter].size() << " cases for " << parameter << std::endl;
            uint_fast64_t runningCounter = wipMatrix.getRowCount();

            std::map<uint_fast64_t, std::vector<uint_fast64_t>> newStates;
            for (auto const& pair : joinedFirstStates[parameter]) {
                // The state that will come before the parametric flip
                auto firstNewState = runningCounter;
                // The state after p 
                auto secondNewState = runningCounter + 1;
                // The state after 1-p
                auto thirdNewState = runningCounter + 2;
                newStates[pair.first] = {firstNewState, secondNewState, thirdNewState};
                runningCounter += 3;
            }
            
            // Extend labeling to more states (Found no better way to do it)
            models::sparse::StateLabeling nextNewLabels(runningCounter);
            for (auto const& label : newLabels.getLabels()) {
                nextNewLabels.addLabel(label);
            }
            for (uint_fast64_t state = 0; state < wipMatrix.getRowCount(); state++) {
                for (auto const& label : newLabels.getLabelsOfState(state)) {
                    nextNewLabels.addLabelToState(label, state);
                }
            }
            newLabels = nextNewLabels;

            storage::SparseMatrixBuilder<RationalFunction> builder;
            for (uint_fast64_t row = 0; row < wipMatrix.getRowCount(); row++) {
                // If from row, we are starting with a transition we will flip...
                if (joinedFirstStates[parameter].count(row)) {
                    // We will join until these transitions
                    auto statesBeforeParams = joinedFirstStates[parameter][row];
                    RationalNumber summedProbability = 0;
                    for (auto const& entry2 : statesBeforeParams) {
                        summedProbability += entry2.second;
                    }
                    // Add the rest of the transitions before, because we need to preserve order
                    for (auto const& entry : wipMatrix.getRow(row)) {
                        // If this state was (the second last) visited in any of the joinedFirstStates, do not add it
                        // Else add it
                        bool addState = true;
                        for (auto const& pair : statesBeforeParams) {
                            auto id = pair.first;
                            auto visitedStates = alreadyVisitedStates[parameter][id];
                            if (std::find(visitedStates.begin(), visitedStates.end(), entry.getColumn()) != visitedStates.end()) {
                                // This state was visited in id's path and we have thus joined it above, so don't add it
                                addState = false;
                            }
                        }
                        if (addState) {
                            builder.addNextValue(row, entry.getColumn(), entry.getValue());
                        }
                    }
                    builder.addNextValue(row, newStates[row][0], utility::convertNumber<RationalFunction>(summedProbability));
                } else { // Else, just add the stuff back into the matrix
                    for (auto const& entry : wipMatrix.getRow(row)) {
                        builder.addNextValue(row, entry.getColumn(), entry.getValue());
                    }
                }
            }
            
            // RationalFunction parameterAsFunction = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(parameter), cache));
            RationalFunction parameterAsFunction = pFunctions[parameter];

            for (auto const& pair : joinedFirstStates[parameter]) {
                auto row = pair.first;
                RationalNumber summedProbability = 0;
                for (auto const& entry2 : joinedFirstStates[parameter][row]) {
                    summedProbability += entry2.second;
                }
                // Add parametric transitions 
                builder.addNextValue(newStates[row][0], newStates[row][1], parameterAsFunction);
                builder.addNextValue(newStates[row][0], newStates[row][2], utility::one<RationalFunction>() - parameterAsFunction);
                
                // Add transitions from p and 1-p states to actual successors
                std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>, std::greater<std::pair<uint_fast64_t, RationalFunction>>> insertInOrder1;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto increasingSuccessor = increasingSuccessors[parameter][id];
                    insertInOrder1.emplace(std::make_pair(increasingSuccessor, utility::convertNumber<RationalFunction>(probability)));
                }
                while (!insertInOrder1.empty()) {
                    auto pair = insertInOrder1.top();
                    insertInOrder1.pop();
                    builder.addNextValue(newStates[row][1], pair.first, pair.second);
                }
                
                std::priority_queue<std::pair<uint_fast64_t, RationalFunction>, std::vector<std::pair<uint_fast64_t, RationalFunction>>, std::greater<std::pair<uint_fast64_t, RationalFunction>>> insertInOrder2;
                for (auto const& pair2 : pair.second) {
                    auto id = pair2.first;
                    RationalNumber probability = pair2.second / summedProbability;
                    auto decreasingSuccessor = decreasingSuccessors[parameter][id];
                    insertInOrder2.emplace(std::make_pair(decreasingSuccessor, utility::convertNumber<RationalFunction>(probability)));
                }
                while (!insertInOrder2.empty()) {
                    auto pair = insertInOrder2.top();
                    insertInOrder2.pop();
                    builder.addNextValue(newStates[row][2], pair.first, pair.second);
                }
                
                for (auto const& label : joiningLabels[parameter][row]) {
                    if (!newLabels.getLabels().count(label)) {
                        newLabels.addLabel(label);
                    }
                    for (uint_fast64_t i = 0; i < 3; i++) {
                        newLabels.addLabelToState(label, newStates[row][i]);
                    }
                }
                if (stateRewardVector) {
                    for (uint_fast64_t i = 0; i < 3; i++) {
                        stateRewardVector->push_back(storm::utility::zero<RationalFunction>());
                    }
                }
            }
            wipMatrix = builder.build();
            
            storage::BitVector statesToKeep(wipMatrix.getRowCount());
            
            // Remove the states that are not needed anymore
            for (auto const& pair : joinedFirstStates[parameter]) {
                auto joiningState = pair.first;
                for (auto const& pair2 : pair.second) {
                    // The identifier - the state the search started from
                    auto firstState = pair2.first;
                    auto visitedStates = alreadyVisitedStates[parameter][firstState];

                    // Delete everything up to the point where the joining state, which is not deleted, is met
                    for (auto const& visitedState : visitedStates) {
                        if (visitedState == joiningState) {
                            break;
                        }
                        std::cout << "Deleting " << visitedState << std::endl;
                        statesToKeep.set(visitedState, true);
                    }

                }
            }

            statesToKeep.complement();
            wipMatrix = wipMatrix.getSubmatrix(false, statesToKeep, statesToKeep, true);
            newLabels = newLabels.getSubLabeling(statesToKeep);
            if (stateRewardVector) {
                std::vector<RationalFunction> newStateRewards(0);
                for (uint_fast64_t i = 0; i < stateRewardVector->size(); i++) {
                    if (statesToKeep[i]) {
                        newStateRewards.push_back((*stateRewardVector)[i]);
                    }
                }
                stateRewardVector = newStateRewards;
            }
        }
        transitionMatrix = wipMatrix;
    }
    models::sparse::Dtmc<RationalFunction> newDTMC(transitionMatrix, newLabels);
    
    if (stateRewardVector) {
        models::sparse::StandardRewardModel<RationalFunction> newRewardModel(*stateRewardVector);
        newDTMC.addRewardModel(*stateRewardName, newRewardModel);
    }

    // Use model simplifier to get rid of unnecessary states
    // storm::transformer::SparseParametricDtmcSimplifier<storm::models::sparse::Dtmc<RationalFunction>> simplifier(newDTMC);
    // STORM_LOG_ASSERT(simplifier.simplify(checkTask.getFormula()), "Could not simplify model.");
    // newDTMC = *simplifier.getSimplifiedModel()->template as<models::sparse::Dtmc<RationalFunction>>();
    return newDTMC;
}

class EqualParameterReducer;
}
}