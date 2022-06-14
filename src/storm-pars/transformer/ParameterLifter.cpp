#include "storm-pars/transformer/ParameterLifter.h"
#include <carl/core/DivisionResult.h>
#include <carl/core/Monomial.h>
#include <carl/core/MultivariatePolynomial.h>
#include <carl/core/PolynomialFactorizationPair.h>
#include <carl/core/Term.h>
#include <carl/core/Variable.h>
#include <carl/core/polynomialfunctions/Factorization.h>
#include <carl/numbers/adaption_native/operations.h>
#include <storm_wrapper.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "storage/geometry/NativePolytope.h"
#include "storage/geometry/Polytope.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/utility/vector.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "utility/constants.h"
#include "utility/logging.h"
#include "utility/macros.h"
#include "storm/storage/geometry/Halfspace.h"

namespace storm {
    namespace transformer {
        
        

        template<typename ParametricType, typename ConstantType>
        ParameterLifter<ParametricType, ConstantType>::ParameterLifter(storm::storage::SparseMatrix<ParametricType> const& pMatrix, std::vector<ParametricType> const& pVector, storm::storage::BitVector const& selectedRows, storm::storage::BitVector const& selectedColumns, bool generateRowLabels, bool useMonotonicityInFuture) {
            // get a mapping from old column indices to new ones
            oldToNewColumnIndexMapping = std::vector<uint_fast64_t>(selectedColumns.size(), selectedColumns.size());
            uint_fast64_t newIndex = 0;
            for (auto const& oldColumn : selectedColumns) {
                oldToNewColumnIndexMapping[oldColumn] = newIndex++;
            }

            // create vector, such that the occuringVariables for all states can be stored
            occurringVariablesAtState = std::vector<std::set<VariableType>>(pMatrix.getColumnCount());
            
            // Stores which entries of the original matrix/vector are non-constant. Entries for non-selected rows/columns are omitted
            auto nonConstMatrixEntries = storm::storage::BitVector(pMatrix.getEntryCount(), false); //this vector has to be resized later
            auto nonConstVectorEntries = storm::storage::BitVector(selectedRows.getNumberOfSetBits(), false);
            // Counters for selected entries in the pMatrix and the pVector
            uint_fast64_t pMatrixEntryCount = 0;
            uint_fast64_t pVectorEntryCount = 0;
            
            // The matrix builder for the new matrix. The correct number of rows and entries is not known yet.
            storm::storage::SparseMatrixBuilder<ConstantType> builder(0, selectedColumns.getNumberOfSetBits(), 0, true, true, selectedRows.getNumberOfSetBits());
            rowGroupToStateNumber = std::vector<uint_fast64_t>();
            uint_fast64_t newRowIndex = 0;
            uint_fast64_t countNonParam = 0;
            
            std::cout << pMatrix << std::endl;
            
            for (auto const& entry : pVector) {
                std::cout << entry << " ";
            }
            std::cout << std::endl;

            for (auto const& rowIndex : selectedRows) {
                builder.newRowGroup(newRowIndex);
                rowGroupToStateNumber.push_back(rowIndex);
                
                // Gather the occurring variables within this row and set which entries are non-constant
                std::set<VariableType> occurringVariables;
                bool constant = true;
                for (auto const& entry : pMatrix.getRow(rowIndex)) {
                    if (selectedColumns.get(entry.getColumn())) {
                        if (!storm::utility::isConstant(entry.getValue())) {
                            storm::utility::parametric::gatherOccurringVariables(entry.getValue(), occurringVariables);
                            nonConstMatrixEntries.set(pMatrixEntryCount, true);
                            constant = false;
                        }
                        ++pMatrixEntryCount;
                    } else {
                        if (!storm::utility::isConstant(entry.getValue())) {
                            storm::utility::parametric::gatherOccurringVariables(entry.getValue(), occurringVariables);
                        }
                    }
                }

                if (constant) {
                    countNonParam++;
                }
                ParametricType const& pVectorEntry = pVector[rowIndex];
                std::set<VariableType> vectorEntryVariables;
                if (!storm::utility::isConstant(pVectorEntry)) {
                    storm::utility::parametric::gatherOccurringVariables(pVectorEntry, vectorEntryVariables);
                    if (generateRowLabels) {
                        // If row labels are to be generated, we do not allow unspecified valuations.
                        // Therefore, we also 'lift' parameters that only occurr on a vector.
                        occurringVariables.insert(vectorEntryVariables.begin(), vectorEntryVariables.end());
                    }
                    nonConstVectorEntries.set(pVectorEntryCount, true);
                }
                ++pVectorEntryCount;

                const bool linearTransitions = onlyLinearTransitions(pMatrix.getRow(rowIndex));
                
                auto countPlaceHolders = 0;
                // Compute the (abstract) valuation for each row
                // Only linear transitions? Then we can do the standard PLA thing - make a rectangle! (First branch)
                if (linearTransitions) {
                    std::vector<RectangleAbstractValuation> rowValuations;
                    

                    // The vertices of the abstract region is the rectangle
                    for (auto const& vertex : getRectangleVertices(occurringVariables)) {
                        rowValuations.push_back(vertex);
                    }
                    

                    for (auto const& val : rowValuations) {
                        if (generateRowLabels && !occurringVariables.empty()) {
                            // FIXME Only storing rectangle valuations as row labels for now
                            rowLabels.push_back(val);
                          }

                        
                        // Insert matrix entries for each valuation. For non-constant entries, a dummy value is inserted and the function and the valuation are collected.
                        // The placeholder for the collected function/valuation are stored in the matrixAssignment. The matrixAssignment is completed after the matrix is finished

                        // Rectangle valuation -> You get the placeholder directly from the local rectangle valuation
                        for (auto const& entry: pMatrix.getRow(rowIndex)) {
                            if(selectedColumns.get(entry.getColumn())) {
                                if(storm::utility::isConstant(entry.getValue())) {
                                    builder.addNextValue(newRowIndex, oldToNewColumnIndexMapping[entry.getColumn()], storm::utility::convertNumber<ConstantType>(entry.getValue()));
                                } else {
                                    builder.addNextValue(newRowIndex, oldToNewColumnIndexMapping[entry.getColumn()], storm::utility::one<ConstantType>());
                                    ConstantType& placeholder = functionValuationCollector.addRectValuation(entry.getValue(), val);
                                    matrixAssignment.push_back(std::pair<typename storm::storage::SparseMatrix<ConstantType>::iterator, ConstantType&>(typename storm::storage::SparseMatrix<ConstantType>::iterator(), placeholder));
                                    countPlaceHolders++;
                                }
                            }
                        }
                        //Insert the vector entry for this row
                        if (storm::utility::isConstant(pVectorEntry)) {
                            vector.push_back(storm::utility::convertNumber<ConstantType>(pVectorEntry));
                        } else {
                            vector.push_back(storm::utility::one<ConstantType>());
                            RectangleAbstractValuation vectorVal(val);
                            for(auto const& vectorVar : vectorEntryVariables) {

                                if (occurringVariables.find(vectorVar) == occurringVariables.end()) {
                                    assert(!generateRowLabels);
                                    vectorVal.addParameterUnspecified(vectorVar);
                                }
                            }
                            ConstantType& placeholder = functionValuationCollector.addRectValuation(pVectorEntry, vectorVal);
                            vectorAssignment.push_back(std::pair<typename std::vector<ConstantType>::iterator, ConstantType&>(typename std::vector<ConstantType>::iterator(), placeholder));
                        }
                        ++newRowIndex;
                    }
                } else {
                    // Big-step valuation -> get all vertices from the valuation and insert them in some order
                    STORM_LOG_ERROR_COND(occurringVariables.size() == 1, "Transitions need to be univariate or linear.");
                    std::vector<RationalFunction> functions;
                    // std::cout << pMatrix << std::endl;
                    for (auto const& entry : pMatrix.getRow(rowIndex)) {
                        if (!entry.getValue().isZero()) {
                            functions.push_back(entry.getValue());
                        }
                    }
                    // if (!storm::utility::isConstant(pVectorEntry)) {
                    //     functions.push_back(pVectorEntry);
                    // }
                    
                    auto val = BigStepAbstractValuation(*occurringVariables.begin(), functions);
                    
                    // Create 2^(n+1) vertices
                    std::size_t const numOfVertices = std::pow(2, val.getVectorIndices().size() + 1);
                    
                    // A BigStep valuation has at most numOfVertices rows with functions.size() columns
                    // This means we have numOfVertices * functions.size() ConstantType references 
                    auto placeholders = functionValuationCollector.addBigStepValuation(val, numOfVertices * functions.size());
                    auto placeholdersIterator = placeholders.begin();

                    for (uint_fast64_t i = 0; i < numOfVertices; i++) {
                        for (auto const& entry: pMatrix.getRow(rowIndex)) {
                            if(selectedColumns.get(entry.getColumn())) {
                                if(storm::utility::isConstant(entry.getValue())) {
                                    STORM_LOG_ERROR("Constant values in big-step valuations currently not supported.");
                                } else {
                                    ConstantType& placeholder = *placeholdersIterator;
                                    builder.addNextValue(newRowIndex, oldToNewColumnIndexMapping[entry.getColumn()], storm::utility::one<ConstantType>());
                                    matrixAssignment.push_back(std::pair<typename storm::storage::SparseMatrix<ConstantType>::iterator, ConstantType&>(typename storm::storage::SparseMatrix<ConstantType>::iterator(), placeholder));
                                    countPlaceHolders++;
                                    placeholdersIterator++;
                                    
                                    if (entry.getValue() == pVectorEntry) {
                                        vector.push_back(storm::utility::one<ConstantType>());
                                        vectorAssignment.push_back(std::pair<typename std::vector<ConstantType>::iterator, ConstantType&>(typename std::vector<ConstantType>::iterator(), placeholder));
                                        countPlaceHolders++;
                                    }
                                }
                            }
                        }
                        // Insert the vector entry for this row
                        if (storm::utility::isConstant(pVectorEntry)) {
                            vector.push_back(storm::utility::convertNumber<ConstantType>(pVectorEntry));
                        } else {
                            // ConstantType& placeholder = *placeholdersIterator;
                            // vector.push_back(storm::utility::one<ConstantType>());
                            // vectorAssignment.push_back(std::pair<typename std::vector<ConstantType>::iterator, ConstantType&>(typename std::vector<ConstantType>::iterator(), placeholder));
                            // countPlaceHolders++;
                        }
                        ++newRowIndex;
                    }

                }

                    
                if (useMonotonicityInFuture) {
                    // Save the occuringVariables of a state, needed if we want to use monotonicity
                    for (auto& var : occurringVariables) {
                        occuringStatesAtVariable[var].insert(rowIndex);
                    }
                    occurringVariablesAtState[rowIndex] = std::move(occurringVariables);
                }
            }

            // Matrix and vector are now filled with constant results from constant functions and place holders for non-constant functions.
            matrix = builder.build(newRowIndex);
            // std::cout << matrix << std::endl;
            vector.shrink_to_fit();
            matrixAssignment.shrink_to_fit();
            vectorAssignment.shrink_to_fit();
            nonConstMatrixEntries.resize(pMatrixEntryCount);

            // Now insert the correct iterators for the matrix and vector assignment
            auto matrixAssignmentIt = matrixAssignment.begin();
            uint_fast64_t startEntryOfRow = 0;
            for (uint_fast64_t group = 0; group < matrix.getRowGroupCount(); ++group) {
                uint_fast64_t startEntryOfNextRow = startEntryOfRow + matrix.getRow(group, 0).getNumberOfEntries();
                for (uint_fast64_t matrixRow = matrix.getRowGroupIndices()[group]; matrixRow < matrix.getRowGroupIndices()[group + 1]; ++matrixRow) {
                    auto matrixEntryIt = matrix.getRow(matrixRow).begin();
                    for(uint_fast64_t nonConstEntryIndex = nonConstMatrixEntries.getNextSetIndex(startEntryOfRow); nonConstEntryIndex < startEntryOfNextRow; nonConstEntryIndex = nonConstMatrixEntries.getNextSetIndex(nonConstEntryIndex + 1)) {
                        matrixAssignmentIt->first = matrixEntryIt + (nonConstEntryIndex - startEntryOfRow);
                        ++matrixAssignmentIt;
                    }
                }
                startEntryOfRow = startEntryOfNextRow;
            }
            STORM_LOG_ASSERT(matrixAssignmentIt == matrixAssignment.end(), "Unexpected number of entries in the matrix assignment.");

            auto vectorAssignmentIt = vectorAssignment.begin();
            for(auto const& nonConstVectorEntry : nonConstVectorEntries) {
                for (uint_fast64_t vectorIndex = matrix.getRowGroupIndices()[nonConstVectorEntry]; vectorIndex != matrix.getRowGroupIndices()[nonConstVectorEntry + 1]; ++vectorIndex) {
                    vectorAssignmentIt->first = vector.begin() + vectorIndex;
                    ++vectorAssignmentIt;
                }
            }
            STORM_LOG_ASSERT(vectorAssignmentIt == vectorAssignment.end(), "Unexpected number of entries in the vector assignment.");
        }
    
        template<typename ParametricType, typename ConstantType>
        void ParameterLifter<ParametricType, ConstantType>::specifyRegion(storm::storage::ParameterRegion<ParametricType> const& region, storm::solver::OptimizationDirection const& dirForParameters) {
            // write the evaluation result of each function,evaluation pair into the placeholders
            functionValuationCollector.evaluateCollectedFunctions(region, dirForParameters);

            //apply the matrix and vector assignments to write the contents of the placeholder into the matrix/vector
            for (auto &assignment : matrixAssignment) {
                STORM_LOG_WARN_COND(!storm::utility::isZero(assignment.second), "Parameter lifting on region " << region.toString() << " affects the underlying graph structure (the region is not strictly well defined). The result for this region might be incorrect.");
                assignment.first->setValue(assignment.second);
            }

            for (auto &assignment : vectorAssignment) {
                *assignment.first = assignment.second;
            
            }
            std::cout << matrix << std::endl;
            for (auto const& entry : vector) {
                std::cout << entry << " ";
            }
            std::cout << std::endl;
            STORM_LOG_ASSERT(matrix.isProbabilistic(), "Matrix not probabilistic!");
        }

        template<typename ParametricType, typename ConstantType>
        uint_fast64_t ParameterLifter<ParametricType, ConstantType>::getRowGroupIndex(uint_fast64_t originalState) const {
            return matrix.getRowGroupIndices()[oldToNewColumnIndexMapping[originalState]];
        }

        template<typename ParametricType, typename ConstantType>
        uint_fast64_t ParameterLifter<ParametricType, ConstantType>::getOriginalStateNumber(uint_fast64_t newState) const {
            return rowGroupToStateNumber[newState];
        }

        template<typename ParametricType, typename ConstantType>
        uint_fast64_t ParameterLifter<ParametricType, ConstantType>::getRowGroupSize(uint_fast64_t originalState) const {
            return matrix.getRowGroupSize(oldToNewColumnIndexMapping[originalState]);
        }
        template<typename ParametricType, typename ConstantType>
        uint_fast64_t ParameterLifter<ParametricType, ConstantType>::getRowGroupCount() const {
            return matrix.getRowGroupCount();
        }
    
        template<typename ParametricType, typename ConstantType>
        storm::storage::SparseMatrix<ConstantType> const& ParameterLifter<ParametricType, ConstantType>::getMatrix() const {
            return matrix;

        }
    
        template<typename ParametricType, typename ConstantType>
        std::vector<ConstantType> const& ParameterLifter<ParametricType, ConstantType>::getVector() const {
            return vector;
        }
        
        template<typename ParametricType, typename ConstantType>
        std::vector<typename ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation> const& ParameterLifter<ParametricType, ConstantType>::getRowLabels() const {
            return rowLabels;
   
        }
        template<typename ParametricType, typename ConstantType>
        bool ParameterLifter<ParametricType, ConstantType>::onlyLinearTransitions(typename storage::SparseMatrix<ParametricType>::const_rows row) const {
            for (auto const& entry : row) {
                std::set<VariableType> occurringVariables;
                storm::utility::parametric::gatherOccurringVariables(entry.getValue(), occurringVariables);
                for (auto const& parameter : occurringVariables) {
                    if (!entry.getValue().derivative(parameter).isConstant()) {
                        return false;
                    }
                }
            }
            return true;
        }
    
        template<typename ParametricType, typename ConstantType>
        std::vector<typename ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation> ParameterLifter<ParametricType, ConstantType>::getRectangleVertices(std::set<VariableType> const& variables) const {
            std::size_t const numOfVertices = std::pow(2, variables.size());
            std::vector<RectangleAbstractValuation> result(numOfVertices);
            
            for (uint_fast64_t vertexId = 0; vertexId < numOfVertices; ++vertexId) {
                //interprete vertexId as a bit sequence
                //the consideredVariables.size() least significant bits of vertex will always represent the next vertex
                //(00...0 = lower boundaries for all variables, 11...1 = upper boundaries for all variables)
                uint_fast64_t variableIndex = 0;
                for (auto const& variable : variables) {
                    if ((vertexId >> variableIndex) % 2 == 0) {
                        result[vertexId].addParameterLower(variable);
                    } else {
                        result[vertexId].addParameterUpper(variable);
                    }
                    ++variableIndex;
                }
            }
            return result;
        }

        template<typename ParametricType, typename ConstantType>
        const std::vector<std::set<typename ParameterLifter<ParametricType, ConstantType>::VariableType>> & ParameterLifter<ParametricType, ConstantType>::getOccurringVariablesAtState() const {
            return occurringVariablesAtState;
        }

        template<typename ParametricType, typename ConstantType>
        std::map<typename ParameterLifter<ParametricType, ConstantType>::VariableType, std::set<uint_fast64_t>> ParameterLifter<ParametricType, ConstantType>::getOccuringStatesAtVariable() const {
            return occuringStatesAtVariable;
        }

        template<typename ParametricType, typename ConstantType>
        bool ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::operator==(RectangleAbstractValuation const& other) const {
            return this->lowerPars == other.lowerPars && this->upperPars == other.upperPars && this->unspecifiedPars == other.unspecifiedPars;
        }
        
        template<typename ParametricType, typename ConstantType>
        void ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::addParameterLower(VariableType const& var) {
            lowerPars.insert(var);
        }
    
        template<typename ParametricType, typename ConstantType>
        void ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::addParameterUpper(VariableType const& var) {
            upperPars.insert(var);
        }
    
        template<typename ParametricType, typename ConstantType>
        void ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::addParameterUnspecified(VariableType const& var) {
            unspecifiedPars.insert(var);
        }

        template<typename ParametricType, typename ConstantType>
        std::size_t ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::getHashValue() const {
            std::size_t seed = 0;
            for (auto const& p : lowerPars) {
                carl::hash_add(seed, p);
            }
            for (auto const& p : upperPars) {
                carl::hash_add(seed, p);
            }
            for (auto const& p : unspecifiedPars) {
                carl::hash_add(seed, p);
            }
            return seed;
        }
        
        template<typename ParametricType, typename ConstantType>
        typename ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::getSubValuation(std::set<VariableType> const& pars) const {
            RectangleAbstractValuation result;
            for (auto const& p : pars) {
                if (std::find(lowerPars.begin(), lowerPars.end(), p) != lowerPars.end()) {
                    result.addParameterLower(p);
                } else if (std::find(upperPars.begin(), upperPars.end(), p) != upperPars.end()) {
                    result.addParameterUpper(p);
                } else if (std::find(unspecifiedPars.begin(), unspecifiedPars.end(), p) != unspecifiedPars.end()) {
                    result.addParameterUnspecified(p);
                } else {
                    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Tried to obtain a subvaluation for parameters that are not specified by this valuation");
                }
            }
            return result;
        }
            
        template<typename ParametricType, typename ConstantType>
        std::set<typename ParameterLifter<ParametricType, ConstantType>::VariableType> const& ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::getLowerParameters() const {
            return lowerPars;
        }
        
        template<typename ParametricType, typename ConstantType>
        std::set<typename ParameterLifter<ParametricType, ConstantType>::VariableType> const& ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::getUpperParameters() const {
            return upperPars;
        }
        
        template<typename ParametricType, typename ConstantType>
        std::set<typename ParameterLifter<ParametricType, ConstantType>::VariableType> const& ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::getUnspecifiedParameters() const {
            return unspecifiedPars;
        }
        
        template<typename ParametricType, typename ConstantType>
        std::vector<storm::utility::parametric::Valuation<ParametricType>> ParameterLifter<ParametricType, ConstantType>::RectangleAbstractValuation::getConcreteValuations(storm::storage::ParameterRegion<ParametricType> const& region) const {
            auto result = region.getVerticesOfRegion(unspecifiedPars);
            for(auto& valuation : result) {
                for (auto const& lowerPar : lowerPars) {
                    valuation.insert(std::pair<VariableType, CoefficientType>(lowerPar, region.getLowerBoundary(lowerPar)));
                }
                for (auto const& upperPar : upperPars) {
                    valuation.insert(std::pair<VariableType, CoefficientType>(upperPar, region.getUpperBoundary(upperPar)));
                }
            }
            return result;
        }

        template<typename ParametricType, typename ConstantType>
        ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::BigStepAbstractValuation(VariableType const parameter, std::vector<storm::RationalFunction> const transitions): parameter(parameter), transitions(transitions) {
            for (uint_fast64_t i = 0; i < transitions.size(); i++) {
                RationalFunction transition = transitions[i];
                std::set<VariableType> occurringVariables;
                storm::utility::parametric::gatherOccurringVariables(transition, occurringVariables);

                STORM_LOG_ERROR_COND(occurringVariables.size() == 1, "Only univariate big-step transitions are supported");
                STORM_LOG_ERROR_COND(transition.denominator().isConstant(), "(Currently) only c_1 * p^a * (1-p)^b + c_2 big-step transitions supported");
                // 
            
                transition.simplify();
                auto constantPart = transition.constantPart();
                std::cout << transition << std::endl;
                std::cout << constantPart << std::endl;

                CoefficientType denominator = transition.denominator().constantPart();
                ConstantType constantDenom = utility::convertNumber<ConstantType>(denominator);
                
                std::cout << constantDenom << std::endl;

                auto nominator = RawPolynomial(transition.nominator());
                
                auto p = *occurringVariables.begin();
                

                
                // auto result1 = nominator.divideBy(parameter);
                // auto result2 = nominator.divideBy(oneMinusParameter);
                
                // std::cout << "p: " << result1.quotient << ", " << result1.remainder << std::endl;
                // std::cout << "1-p: " << result2.quotient << ", " << result2.remainder << std::endl;

                RawPolynomial currentNominator = nominator;
                
                // std::cout << result.quotient << ", " << result.remainder << std::endl;


                // ConstantType offset = utility::zero<ConstantType>();
                // if (nominator.trailingTerm().isConstant()) {
                //     auto coeff = nominator.trailingTerm().coeff();
                //     offset = utility::convertNumber<ConstantType>(coeff / denominator);
                //     transition -= coeff;
                // }

                // transition.simplify();
                
                
                // auto factorization = carl::factorization(nominator);
                
                // ConstantType offset = utility::zero<ConstantType>();
                // // if (factorization.size() == 1 && factorization.begin()->first != parameter && factorization.begin()->first != oneMinusParameter) {
                // if (factorization.size() == 1 && nominator.nrTerms() > 1) {
                //     auto offsetAsCln = nominator.constantPart();
                //     nominator -= offsetAsCln;
                //     factorization = carl::factorization(nominator);
                //     offset = utility::convertNumber<ConstantType>(offsetAsCln / denominator);
                // }
                
                // auto parameter = RawPolynomial(p);
                // auto oneMinusParameter = RawPolynomial(1) - parameter;
                // uint_fast64_t a = 0;
                // uint_fast64_t b = 0;
                // ConstantType constant = utility::one<ConstantType>();
                // for (auto const& element : factorization) {
                //     if (element.first == parameter) {
                //         a = element.second;       
                //     } else if (element.first == oneMinusParameter || element.first == -oneMinusParameter) {
                //         b = element.second;
                //     } else if (element.first.isConstant()) {
                //         constant = utility::abs(utility::convertNumber<ConstantType>(element.first.constantPart() / denominator));
                //     }
                
                // }
                // 
                auto result = tryDecomposing(nominator, true);
                
                STORM_LOG_ASSERT(result, "Polynomial had no decomposition.");
                
                uint_fast64_t a = result->first.first;
                uint_fast64_t b = result->first.second;
                ConstantType constant = result->second.first / constantDenom;
                ConstantType offset = result->second.second / constantDenom;
                
                // Polynomial is constant * p^a * (1-p)^b + offset
                //
                auto aAndBPair = std::make_pair(a, b);
                this->asAndBs.push_back(aAndBPair);
                this->constantsAndOffsets.push_back(std::make_pair(utility::convertNumber<ConstantType>(constant), utility::convertNumber<ConstantType>(offset)));
                
                // The maximum of the polynomial part lies at a / (a + b), so compute that
                // It is corrected for constant and offset later, not now
                this->vectorIndices[aAndBPair].push_back(asAndBs.size() - 1);
                if (!this->maxima.count(aAndBPair)) {
                    CoefficientType maximumCoeff;
                    if (a != 0 || b != 0) {
                        maximumCoeff = utility::convertNumber<CoefficientType>(a) / utility::convertNumber<CoefficientType>(a + b);
                    } else {
                        maximumCoeff = utility::zero<CoefficientType>();
                    }
                    ConstantType maximum = utility::convertNumber<ConstantType>(maximumCoeff);
                    

                    std::map<VariableType, CoefficientType> substitution;
                    substitution.emplace(p, maximumCoeff);

                    this->maxima.emplace(aAndBPair, std::make_pair(maximum, utility::convertNumber<ConstantType>(transition.evaluate(substitution))));
                }
            }
        }

        template<typename ParametricType, typename ConstantType>
        boost::optional<std::pair<std::pair<uint_fast64_t, uint_fast64_t>, std::pair<ConstantType, ConstantType>>> ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::tryDecomposing(RawPolynomial polynomial, bool firstIteration) {
            auto parameterPol = RawPolynomial(parameter);
            auto oneMinusParameter = RawPolynomial(1) - parameterPol;
            if (polynomial.isConstant()) {
                return std::make_pair(std::make_pair((uint_fast64_t) 0, (uint_fast64_t) 0), std::make_pair(utility::convertNumber<ConstantType>(polynomial.constantPart()), utility::zero<ConstantType>()));
            }
            auto byOneMinusP = polynomial.divideBy(oneMinusParameter);
            if (byOneMinusP.remainder.isZero()) {
                auto recursiveResult = tryDecomposing(byOneMinusP.quotient, false);
                if (recursiveResult) {
                    return std::make_pair(std::make_pair(recursiveResult->first.first, recursiveResult->first.second + 1), recursiveResult->second);
                }
            }
            auto byP = polynomial.divideBy(parameterPol);
            if (byP.remainder.isZero()) {
                auto recursiveResult = tryDecomposing(byP.quotient, false);
                if (recursiveResult) {
                    return std::make_pair(std::make_pair(recursiveResult->first.first + 1, recursiveResult->first.second), recursiveResult->second);
                }
            }
            if (!firstIteration) {
                return boost::none;
            }
            if (byOneMinusP.remainder.isConstant()) {
                auto rem1 = utility::convertNumber<ConstantType>(byOneMinusP.remainder.constantPart());
                auto recursiveResult = tryDecomposing(byOneMinusP.quotient, false);
                if (recursiveResult) {
                    STORM_LOG_ASSERT(recursiveResult->second.second == 0, "");
                    return std::make_pair(std::make_pair(recursiveResult->first.first, recursiveResult->first.second + 1), 
                        std::pair<ConstantType, ConstantType>(recursiveResult->second.first, rem1));
                }
            }
            if (byP.remainder.isConstant()) {
                auto rem2 = utility::convertNumber<ConstantType>(byP.remainder.constantPart());
                auto recursiveResult = tryDecomposing(byP.quotient, false);
                if (recursiveResult) {
                    STORM_LOG_ASSERT(recursiveResult->second.second == 0, "");
                    return std::make_pair(std::make_pair(recursiveResult->first.first + 1, recursiveResult->first.second), 
                        std::pair<ConstantType, ConstantType>(recursiveResult->second.first, rem2));
                }
            }
            return boost::none;
        }
        
        template<typename ParametricType, typename ConstantType>
        bool ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::operator==(BigStepAbstractValuation const& other) const {
            return this->parameter == other.parameter && this->transitions == other.transitions;
        }

        template<typename ParametricType, typename ConstantType>
        typename ParameterLifter<ParametricType, ConstantType>::VariableType const& ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getParameter() const {
            return parameter;
        }

        template<typename ParametricType, typename ConstantType>
        uint_fast64_t ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getNumTransitions() const {
            return this->transitions.size();
        }

        template<typename ParametricType, typename ConstantType>
        std::vector<RationalFunction> const& ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getTransitions() const {
            return this->transitions;
        }

        template<typename ParametricType, typename ConstantType>
        std::vector<std::pair<ConstantType, ConstantType>> const& ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getConstantsAndOffsets() const {
            return this->constantsAndOffsets;
        }

        template<typename ParametricType, typename ConstantType>
        std::vector<std::pair<uint_fast64_t, uint_fast64_t>> const& ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getAsAndBs() const {
            return this->asAndBs;
        }

        template<typename ParametricType, typename ConstantType>
        std::map<std::pair<uint_fast64_t, uint_fast64_t>, std::pair<ConstantType, ConstantType>> const& ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getMaxima() const {
            return this->maxima;
        }

        template<typename ParametricType, typename ConstantType>
       std::map<std::pair<uint_fast64_t, uint_fast64_t>, std::vector<uint_fast64_t>> const& ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getVectorIndices() const {
            return this->vectorIndices;
        }
        
        
        template<typename ParametricType, typename ConstantType>
        std::size_t ParameterLifter<ParametricType, ConstantType>::BigStepAbstractValuation::getHashValue() const {
            std::size_t seed = 0;
            carl::hash_add(seed, parameter);
            for (auto const& t : transitions) {
                carl::hash_add(seed, t);
            }
            return seed;
        }
        
        template<typename ParametricType, typename ConstantType>
        ConstantType& ParameterLifter<ParametricType, ConstantType>::FunctionValuationCollector::addRectValuation(ParametricType const& function, RectangleAbstractValuation const& valuation) {
            ParametricType simplifiedFunction = function;
            storm::utility::simplify(simplifiedFunction);
            std::set<VariableType> variablesInFunction;
            storm::utility::parametric::gatherOccurringVariables(simplifiedFunction, variablesInFunction);
            
            auto simplifiedValuation = valuation.getSubValuation(variablesInFunction);
            // If the given FunctionValuation already exists in the map, this means there is already a ConstantType reference
            // that needs to be filled with the value from the valuation. If this reference prevents insertion, it is returned
            // in insertionRes.first->second and will then be inserted into the matrix assignment. This causes the function to
            // only be calculated once, the result written to the ConstantType reference, and the ConstantType reference is then
            // contained in both places where this valuation occurs.
            auto insertionRes = collectedRectangleValuations.insert(std::pair<FunctionValuation, ConstantType>(FunctionValuation(std::move(simplifiedFunction), std::move(simplifiedValuation)), storm::utility::one<ConstantType>()));
            return insertionRes.first->second;
        }

        template<typename ParametricType, typename ConstantType>
        std::vector<std::reference_wrapper<ConstantType>> ParameterLifter<ParametricType, ConstantType>::FunctionValuationCollector::addBigStepValuation(BigStepAbstractValuation const& valuation, uint_fast64_t numOfPlaceholders) {
            std::vector<ConstantType> vertices(numOfPlaceholders, ConstantType(1));
            std::vector<std::reference_wrapper<ConstantType>> vertexReferences(vertices.begin(), vertices.end());
            auto insertionRes = collectedBigStepValuations.insert(std::make_pair(valuation, std::make_pair(std::move(vertices), std::move(vertexReferences))));
            return insertionRes.first->second.second;
        }
    
        template<typename ParametricType, typename ConstantType>
        void ParameterLifter<ParametricType, ConstantType>::FunctionValuationCollector::evaluateCollectedFunctions(storm::storage::ParameterRegion<ParametricType> const& region, storm::solver::OptimizationDirection const& dirForUnspecifiedParameters) {
            // Evaluate rectangle transitions
            for (auto &collectedFunctionValuationPlaceholder : collectedRectangleValuations) {
                ParametricType const &function = collectedFunctionValuationPlaceholder.first.first;
                RectangleAbstractValuation const &abstrValuation = collectedFunctionValuationPlaceholder.first.second;
                ConstantType &placeholder = collectedFunctionValuationPlaceholder.second;

                auto concreteValuations = abstrValuation.getConcreteValuations(region);
                auto concreteValuationIt = concreteValuations.begin();
                placeholder = storm::utility::convertNumber<ConstantType>(storm::utility::parametric::evaluate(function, *concreteValuationIt));
                for (++concreteValuationIt; concreteValuationIt != concreteValuations.end(); ++concreteValuationIt) {
                    ConstantType currentResult = storm::utility::convertNumber<ConstantType>(storm::utility::parametric::evaluate(function, *concreteValuationIt));
                    if (storm::solver::minimize(dirForUnspecifiedParameters)) {
                        placeholder = std::min(placeholder, currentResult);
                    } else {
                        placeholder = std::max(placeholder, currentResult);
                    }
                }
            }
            for (auto &collectedBigStepPlaceholders : collectedBigStepValuations) {
                BigStepAbstractValuation bigStepTransition = collectedBigStepPlaceholders.first;
                std::vector<std::reference_wrapper<ConstantType>> placeholders = collectedBigStepPlaceholders.second.second;
                
                // Compute lower bound and upper bound of function for every transition
                std::vector<ConstantType> lowerBounds(bigStepTransition.getVectorIndices().size());
                std::vector<ConstantType> upperBounds(bigStepTransition.getVectorIndices().size());
                std::vector<ConstantType> sumsOfConstants(bigStepTransition.getVectorIndices().size());
                
                auto maxima = bigStepTransition.getMaxima();
                auto transitions = bigStepTransition.getTransitions();
                auto constants = bigStepTransition.getConstantsAndOffsets();
                auto p = bigStepTransition.getParameter();
                
                // The functions we consider are always p^a * (1-p)^b.
                // These always have a maximum at a / (a + b), are monotone increasing before and monotone decreasing after.
                // So we are using a case distinction to compute the lower and upper bounds.
                
                // Populate lowerBounds and upperBounds at position i
                uint_fast64_t i = 0;
                std::vector<std::pair<uint_fast64_t, uint_fast64_t>> lowerUpperIndicesToPair;
                for (auto const& pair : bigStepTransition.getVectorIndices()) {
                    lowerUpperIndicesToPair.push_back(pair.first);
                    auto f = transitions[pair.second[0]];
                    
                    ConstantType sumOfConstants = 0;
                    for (uint_fast64_t index : pair.second) {
                        sumOfConstants += constants[index].first;
                    }
                    sumsOfConstants[i] = sumOfConstants;

                    CoefficientType lowerPCoeff = region.getLowerBoundary(p);
                    CoefficientType upperPCoeff = region.getUpperBoundary(p);
                    ConstantType lowerP = utility::convertNumber<ConstantType>(region.getLowerBoundary(p));
                    ConstantType upperP = utility::convertNumber<ConstantType>(region.getUpperBoundary(p));
                    
                    // Compute function values at left and right ends 
                    std::map<VariableType, CoefficientType> substitution;
                    substitution[p] = lowerPCoeff;
                    auto left = utility::convertNumber<ConstantType>(f.evaluate(substitution));
                    substitution[p] = upperPCoeff;
                    auto right = utility::convertNumber<ConstantType>(f.evaluate(substitution));

                    if (lowerP <= maxima[pair.first].first && upperP >= maxima[pair.first].first) {
                        // Case 1: The valuation is around the maximum of the function.
                        // Lower bound is the minimum of left and right
                        lowerBounds[i] = utility::min(left, right);
                        // Upper bound is easy
                        upperBounds[i] = maxima[pair.first].second;
                    } else if (lowerP > maxima[pair.first].first) {
                        // Case 2: The valuation is on the right of the maximum => monotone decreasing
                        lowerBounds[i] = right;
                        upperBounds[i] = left;
                    } else if (upperP < maxima[pair.first].first) {
                        // Case 3: The valuation is on the left of the maximum => monotone increasing
                        lowerBounds[i] = left;
                        upperBounds[i] = right;
                    }
                    i++;
                }
                
                // lowerBounds[i] *= constants[i].first;
                // upperBounds[i] *= constants[i].first;

                // lowerBounds[i] += constants[i].second;
                // upperBounds[i] += constants[i].second;
                
                // std::cout << "aa" << std::endl;
                // Build the polytope
                // 
                uint_fast64_t numHalfspaces = bigStepTransition.getVectorIndices().size();

                std::vector<storage::geometry::Halfspace<ConstantType>> halfspaces;
                // We assign new variables to the transitions: x_1, x_2, x_3, ..., x_n-1
                // The last transition has the value 1 - x_1 - x_2 - ... - x_n-1
                // Populate constraints for variables
                for (uint_fast64_t i = 0; i < numHalfspaces - 1; i++) {
                    std::vector<ConstantType> normalVectorLower(numHalfspaces - 1);
                    std::vector<ConstantType> normalVectorUpper(numHalfspaces - 1);
                    
                    normalVectorLower[i] = -utility::one<ConstantType>();
                    normalVectorUpper[i] = utility::one<ConstantType>();
                    
                    storage::geometry::Halfspace<ConstantType> lowerHalfspace(normalVectorLower, -lowerBounds[i] * sumsOfConstants[i]);
                    storage::geometry::Halfspace<ConstantType> upperHalfspace(normalVectorUpper, upperBounds[i] * sumsOfConstants[i]);
                    
                    halfspaces.push_back(lowerHalfspace);
                    halfspaces.push_back(upperHalfspace);
                }
                
                std::vector<ConstantType> normalVectorLowerLast(numHalfspaces - 1);
                std::vector<ConstantType> normalVectorUpperLast(numHalfspaces - 1);
                for (uint_fast64_t i = 0; i < numHalfspaces - 1; i++) {
                    normalVectorLowerLast[i] = utility::one<ConstantType>();
                    normalVectorUpperLast[i] = -utility::one<ConstantType>();
                }
                storage::geometry::Halfspace<ConstantType> lowerHalfspaceLast(normalVectorLowerLast, -lowerBounds[numHalfspaces - 1] + utility::one<ConstantType>());
                storage::geometry::Halfspace<ConstantType> upperHalfspaceLast(normalVectorUpperLast, upperBounds[numHalfspaces - 1] - utility::one<ConstantType>());

                halfspaces.push_back(lowerHalfspaceLast);
                halfspaces.push_back(upperHalfspaceLast);
                
                auto polytope = storage::geometry::NativePolytope<ConstantType>(halfspaces);

                auto vertices = polytope.getVertices();
                
                // for (auto const& vertex : polytope.getVertices()) {
                //     for (auto const& element : vertex) {
                //         std::cout << element << " ";
                //     }
                //     std::cout << std::endl;
                // }
                
                // std::cout << std::endl;

                // Compute how many unique rows there are in total
                uint_fast64_t numUniqueRows = 0;
                for (auto const& aAndB : lowerUpperIndicesToPair) {
                    numUniqueRows += bigStepTransition.getVectorIndices().at(aAndB).size();
                }
                
                // Modulo through the vertices so that we populate every row
                for (uint_fast64_t vertexIndex = 0; vertexIndex < vertices.size(); vertexIndex++) {
                    auto vertex = vertices[vertexIndex];

                    ConstantType lastValue = 1;
                    for (auto const& value : vertex) {
                        lastValue -= value;
                    }
                    vertex.push_back(lastValue);

                    uint_fast64_t numberOfRows = placeholders.size() / bigStepTransition.getNumTransitions();
                    for (uint_fast64_t row = vertexIndex; row < numberOfRows; row += numUniqueRows) {
                        for (uint_fast64_t j = 0; j < vertex.size(); j++) {
                            auto aAndB = lowerUpperIndicesToPair[j];
                            auto transitionIndices = bigStepTransition.getVectorIndices().at(aAndB);
                            for (auto const& i : transitionIndices) {
                                ConstantType& reference = placeholders[row * bigStepTransition.getNumTransitions() + i];
                                reference = (constants[i].first / sumsOfConstants[j]) * vertex[j] + constants[i].second;
                            }
                        }
                    }
                }
            }
        }
        
        template class ParameterLifter<storm::RationalFunction, double>;
        template class ParameterLifter<storm::RationalFunction, storm::RationalNumber>;
    }
}
